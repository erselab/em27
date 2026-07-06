"""Single-spectrum EM27/SUN retrieval against PROFFAST (M3 → M4).

Encapsulates the M3 retrieval so it can be driven both from the notebook and
from the M4 batch loop (``scripts/run_m4.py``).  Physics lives in ``gert``;
this only wires EM27 readers / instrument / scene into a ``GERTRetrieval``.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np

from gert.forward_model import ForwardModel
from gert.rt_solver import TransmissionSolver
from gert.retrieval import GERTRetrieval, StateVector

from .readers import read_spectrum, nearest_map, read_ils_list
from .scene import map_to_atmosphere
from .instrument import (build_em27_instrument, ils_physical, ils_gaussian,
                         EM27_OPD_CM, EM27_WINDOWS)

# Effective ILS resolution (χ² minimum). Re-tuned to 0.44 cm⁻¹ once the retrieval
# was moved onto the observation grid (the ILS integrates the model onto the
# measured wavenumbers — no data interpolation); the earlier 0.50 partly
# compensated the smoothing of linear interpolation onto a finer grid.  This is
# the FWHM the *physical* (NB-apodized ME/PE) ILS is built to; a per-band
# ``ils_scale`` nuisance then refines the width in the retrieval.
RES_EFF_CM = 0.44


def build_physical_ils(data_dir, fwhm_cm=RES_EFF_CM, apod="nb_medium",
                       opd_cm=EM27_OPD_CM):
    """The default EM27 ILS: NB-apodized self-apodizing ME/PE kernel from
    ``ils_list.csv`` at ``fwhm_cm``.  Carries an analytic g'(δ) for closed-form
    dispersion / ILS-width Jacobians."""
    row = read_ils_list(Path(data_dir) / "ils_list.csv").iloc[0]
    return ils_physical(opd_cm, float(row["ME1"]), float(row["PE1"]),
                        float(row["ME2"]), float(row["PE2"]),
                        apod=apod, fwhm_cm=fwhm_cm)
GASES = ("co2", "ch4", "co", "n2o", "h2o")
# O2/XAIR is opt-in: pass windows=EM27_WINDOWS and gases=GASES_WITH_O2.  The O2
# 1.27 µm band under-absorbs in the current ABSCO (no CIA), so XAIR is biased
# (~0.54 at res 0.5) and O2 is NOT in the default Xgas retrieval.
GASES_WITH_O2 = ("co2", "ch4", "co", "n2o", "h2o", "o2_1p27")
O2_DRY_MF = 0.2095   # assumed dry-air O2 mole fraction (PROFFAST XAIR convention)


def solar_doppler_ms(date, ut_h, lat, lon, alt_km):
    """Observer heliocentric velocity toward the Sun [m/s] (0.0 if unavailable)."""
    try:
        from astropy.time import Time
        from astropy.coordinates import EarthLocation, get_sun
        from astropy.utils.iers import conf
        import astropy.units as u
        conf.auto_download = False               # avoid per-call network stalls
        loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=alt_km * 1000 * u.m)
        t = Time(f"{date}T00:00:00") + ut_h * u.hour
        return float(get_sun(t).radial_velocity_correction(
            "heliocentric", obstime=t, location=loc).to(u.m / u.s).value)
    except Exception:
        return 0.0


def _build_y_obs(inst, y0, meas_wn, meas_sp, merr, baseline_order=2):
    """Measured spectrum on the GERT radiometric scale (per-window polynomial
    continuum of degree ``baseline_order``), wavelength order, plus diagonal Sy.

    The instrument channels are the observation wavenumbers (``obs_grid``), so
    the measured values are used **directly — no interpolation** — with the ILS
    having integrated the model onto these same points.  ``baseline_order`` sets
    the per-window multiplicative continuum P = (Σ_k c_k xᵏ)·yg fitted to the
    measurement (2 = default; lower orders absorb less broadband structure).
    """
    yobs_parts, sig_parts, off = [], [], 0
    for w in inst.windows:
        n = w.n_channels
        yg = y0[off:off + n][::-1]               # ascending wn (R_band is wavelength order)
        off += n
        wn = w.wn_instrument                      # = measured wavenumbers in this window
        m = (meas_wn >= w.wn_min) & (meas_wn <= w.wn_max)
        ym = meas_sp[m]                           # measured values at the obs grid (no interp)
        if len(ym) != n:                          # obs_grid not used / mismatch -> fall back
            ym = np.interp(wn, meas_wn, meas_sp)
        x = (wn - wn.mean()) / (np.ptp(wn) / 2.0)
        Xp = np.vstack([x ** k for k in range(baseline_order + 1)])   # (order+1, n)
        A = (yg[None, :] * Xp).T
        coef, *_ = np.linalg.lstsq(A, ym, rcond=None)
        P = (coef[:, None] * Xp).sum(0)
        yo = ym / P
        yobs_parts.append(yo[::-1])
        sig_parts.append((np.abs(yo) * merr + 1e-12)[::-1])
    return np.concatenate(yobs_parts), np.concatenate(sig_parts)


def build_eof_basis(resid_npz, n_per_band=1):
    """ACOS-style per-band EOF radiance-correction basis in the retrieval y-layout.

    From an ensemble residual npz (normalized residuals, ascending-wn per window),
    take the top ``n_per_band`` EOFs of each window and place them in the
    concatenated **wavelength-order** y-layout (reversing each window, zero outside
    it).  Returns ``(eof_basis [n_y, n_eof], meta [(label, k), …])`` for
    ``retrieve_spectrum(eof_basis=…)``.  Unit-norm patterns; the retrieval fits a
    coefficient per column.
    """
    z = np.load(resid_npz, allow_pickle=True)
    labels = [str(x) for x in z["win_labels"]]; nch = z["win_nchan"]; R = z["resid"]
    bnd = np.concatenate([[0], np.cumsum(nch)]); n_y = int(bnd[-1])
    npb = [n_per_band] * len(labels) if isinstance(n_per_band, int) else list(n_per_band)
    cols, meta = [], []
    for j, lab in enumerate(labels):
        sl = slice(int(bnd[j]), int(bnd[j + 1]))
        Rb = R[:, sl]
        _U, _S, Vt = np.linalg.svd(Rb - Rb.mean(0), full_matrices=False)
        for k in range(npb[j]):
            pat = np.zeros(n_y)
            pat[sl] = Vt[k][::-1]                 # ascending wn -> wavelength order
            cols.append(pat); meta.append((lab, k))
    return np.array(cols).T, meta                 # (n_y, n_eof)


def build_h2o_shape_profile(atm):
    """Fixed column-neutral H₂O redistribution ramp φ(p) at level resolution.

    ``vmr(p) = vmr_prior(p)·(1 + s + β·φ(p))``: ``s`` (=h2o_scale−1) carries the
    total column, ``β`` (=h2o_shape) the redistribution.  φ is a smooth ramp
    linear in pressure, ``φ_k ∝ (p_k − p̄)``, with the pivot ``p̄`` set to the
    **water-column-weighted mean pressure** so that the first-order column change
    ``Σ_k u_k·vmr_k·φ_k = 0`` (column-neutral), where ``u_k = ½(Δp_{k−1}+Δp_k)``
    is the dry-air partial-column weight around level k.  Normalized to unit
    **water-weighted RMS** (``Σ w_k φ_k² / Σ w_k = 1``) so β is the typical
    fractional water tilt *where the water actually is* (β>0 piles water toward
    the surface); the 0.15 prior then reads as a ~15 % redistribution.
    """
    p = np.asarray(atm.p_levels, dtype=float)          # (n_lev,), surface→top
    vmr = np.asarray(atm.gases["h2o"], dtype=float)
    dp = np.abs(np.diff(p))                             # (n_lay,) layer thickness
    u = np.zeros_like(p)                               # per-level air-mass weight
    u[:-1] += 0.5 * dp
    u[1:]  += 0.5 * dp
    w = u * vmr                                         # water-column weight
    p_bar = float((w * p).sum() / (w.sum() + 1e-30))    # water-weighted mean p
    phi = p - p_bar
    rms = np.sqrt((w * phi**2).sum() / (w.sum() + 1e-30))
    return phi / rms if rms > 0 else phi


def retrieve_spectrum(pick, inv, data_dir, absco, solar, *,
                      res_eff=RES_EFF_CM, merr=0.005, gases=GASES,
                      windows=EM27_WINDOWS[:3], max_iter=12, freeze_gas=None,
                      use_doppler=True, dispersion_order=1, ils=None,
                      retrieve_ils_scale=True, ils_scale_uncert=0.05,
                      eof_basis=None, eof_uncert=1.0,
                      retrieve_h2o_shape=False, h2o_shape_uncert=0.15,
                      baseline_order=2, map_path=None, return_spectra=False):
    """Run the M3 retrieval for one ``.BIN`` spectrum.

    Returns a dict with the retrieved Xgas, PROFFAST Xgas, χ², SZA, time, and
    convergence flag.  ``absco`` / ``solar`` are passed in so they load once.

    ILS: by default a **physical** NB-apodized self-apodizing ME/PE ILS at
    ``res_eff`` cm⁻¹ (from ``ils_list.csv``), with a **per-band ILS width scale**
    (``ils_scale_{b}``) retrieved jointly via analytic Jacobians.  Pass an
    explicit ``ils`` (e.g. ``ils_gaussian(0.44)``) to override the shape, and
    ``retrieve_ils_scale=False`` to freeze the width.

    Ablation knobs (M5): ``use_doppler`` toggles the solar Doppler pre-shift,
    ``dispersion_order=0`` disables the dispersion nuisance.
    """
    data_dir = Path(data_dir)
    row = inv.loc[pick]
    spec = read_spectrum(data_dir / "260406_spectra/cal" / pick)

    # merge the SM (extended) channel for the CO window (~4210 cm-1)
    meas_wn, meas_sp = spec["wn"], spec["spectrum"]
    sm_path = data_dir / "260406_spectra/cal" / pick.replace("SN", "SM")
    if sm_path.exists():
        sm = read_spectrum(sm_path)
        order = np.argsort(np.concatenate([sm["wn"], spec["wn"]]))
        meas_wn = np.concatenate([sm["wn"], spec["wn"]])[order]
        meas_sp = np.concatenate([sm["spectrum"], spec["spectrum"]])[order]

    md = spec["metadata"]
    # By default the nearest 3-hourly .map is used; pass ``map_path`` to hold a
    # single prior atmospheric state fixed across a whole run (surface pressure
    # is still per-scan from gndP, so only the profile *shape* is held fixed).
    if map_path is None:
        map_path = nearest_map(data_dir / "map", md["time_ut_h"], md["date"])
    atm = map_to_atmosphere(map_path, p_surface_pa=float(row["gndP"]) * 100.0)
    prior = {"co2": atm.column_xgas("co2") * 1e6,
             "ch4": atm.column_xgas("ch4") * 1e9,
             "co":  atm.column_xgas("co") * 1e9}

    sd = solar_doppler_ms(md["date"], md["time_ut_h"],
                          float(row["latdeg"]), float(row["londeg"]),
                          md["altitude_km"]) if use_doppler else 0.0

    from gert.geometry import Geometry
    geo = Geometry(sza=spec["sza_deg"], vza=0.0, raa=0.0,
                   observer_altitude=md["altitude_km"])
    ils_use = ils if ils is not None else build_physical_ils(data_dir, fwhm_cm=res_eff)
    inst = build_em27_instrument(windows=windows, obs_wn=meas_wn, ils=ils_use)
    nw = len(inst.windows)
    fm = ForwardModel(atm, absco, inst, geo,
                      solver=TransmissionSolver(jacobians=True), solar_spectrum=solar)
    y0 = fm.run(albedo=np.zeros(nw), solar_doppler=sd).y
    y_obs, sig = _build_y_obs(inst, y0, meas_wn, meas_sp, merr, baseline_order=baseline_order)

    sv = StateVector.transmission_scaling(
        n_bands=nw, gases=list(gases),
        gas_uncerts={"co": 0.30, "n2o": 0.20, "h2o": 0.30},
        co2_uncert=0.10, ch4_uncert=0.20, solar_gain_uncert=0.20,
        include_dispersion=dispersion_order > 0,
        dispersion_order=max(dispersion_order, 1), dispersion_uncert=0.10,
        include_ils_scale=retrieve_ils_scale, ils_scale_uncert=ils_scale_uncert,
        include_eof=eof_basis is not None,
        n_eof=(0 if eof_basis is None else eof_basis.shape[1]), eof_uncert=eof_uncert,
        include_h2o_shape=retrieve_h2o_shape,
        h2o_shape_profile=build_h2o_shape_profile(atm) if retrieve_h2o_shape else None,
        h2o_shape_uncert=h2o_shape_uncert)
    sv.freeze("p_scale")
    if freeze_gas:
        for e in sv.elements:
            mol = e.transform[:-6] if e.transform.endswith("_scale") else None
            if mol in freeze_gas:
                e.prior = freeze_gas[mol]; e.value = e.prior
                sv.freeze(e.name)

    if eof_basis is not None and eof_basis.shape[0] != len(y_obs):
        raise ValueError(f"eof_basis n_y={eof_basis.shape[0]} != retrieval y={len(y_obs)}")
    res = GERTRetrieval(fm, y_obs, np.diag(1.0 / sig**2), sv,
                        prior_albedo=np.zeros(nw), analytical_jacobians=True,
                        solar_doppler=sd, max_iter=max_iter, verbose=False,
                        eof_basis=eof_basis).run()
    eof_ret = {int(e.transform[len("eof_"):]): float(res.x_ret[i])
               for i, e in enumerate(sv.elements) if e.transform.startswith("eof_")}
    scale = {e.transform[:-6]: res.x_ret[i] for i, e in enumerate(sv.elements)
             if e.transform.endswith("_scale") and e.transform != "p_scale"
             and not e.transform.startswith("ils_scale_")}
    ils_scale_ret = {int(e.transform[len("ils_scale_"):]): float(res.x_ret[i])
                     for i, e in enumerate(sv.elements)
                     if e.transform.startswith("ils_scale_")}

    # XAIR (PROFFAST convention): 0.2095·N_dry / N_O2_retrieved.  With the atm
    # built from gndP, N_dry cancels ->  XAIR = 0.2095 / (o2_scale · <o2>_prior).
    xair_gert = np.nan
    if "o2_1p27" in scale:
        o2_prior_mf = atm.column_xgas("o2_1p27")          # ≈ 0.2095 (dry)
        xair_gert = O2_DRY_MF / (scale["o2_1p27"] * o2_prior_mf)

    t_offset = next((float(res.x_ret[i]) for i, e in enumerate(sv.elements)
                     if e.transform == "T_offset"), np.nan)

    out = {
        "spectrum": pick,
        "ut_h": md["time_ut_h"],
        "sza": spec["sza_deg"],
        "airmass": 1.0 / np.cos(np.radians(spec["sza_deg"])),
        "chi2": res.chisq_reduced,
        "t_offset": t_offset,
        "h2o_scale": float(scale.get("h2o", np.nan)),
        "h2o_shape": next((float(res.x_ret[i]) for i, e in enumerate(sv.elements)
                           if e.transform == "h2o_shape"), np.nan),
        "converged": bool(res.converged),
        "xco2_gert": prior["co2"] * scale["co2"] if "co2" in scale else np.nan,
        "xch4_gert": prior["ch4"] * scale["ch4"] if "ch4" in scale else np.nan,
        "xco_gert":  prior["co"]  * scale["co"]  if "co"  in scale else np.nan,
        "xair_gert": xair_gert,
        "xh2o_gert": (atm.column_xgas("h2o") * scale["h2o"] * 1e6
                      if "h2o" in scale else np.nan),          # ppm
        "xco2_proffast": float(row["XCO2"]),
        "xch4_proffast": float(row["XCH4"]) * 1000.0,
        "xco_proffast":  float(row["XCO"]) * 1000.0,
        "xair_proffast": float(row["XAIR"]),
        "xh2o_proffast": float(row["XH2O"]) if "XH2O" in row else np.nan,
    }
    # retrieved per-band ILS width scale (× the physical NB base FWHM ~res_eff)
    for b, s in ils_scale_ret.items():
        out[f"ils_scale_{b}"] = s
        out[f"ils_fwhm_{b}"] = res_eff * s
    for k, cval in eof_ret.items():                       # empirical EOF coefficients
        out[f"eof_{k}"] = cval
    if return_spectra:
        # per-window residual (y_obs − y_ret) on ascending-wavenumber grids
        wn_all, resid_all, yret_all, off = [], [], [], 0
        for w in inst.windows:
            n = w.n_channels
            yo = y_obs[off:off + n][::-1]
            yr = res.y_ret[off:off + n][::-1]
            off += n
            wn_all.append(w.wn_instrument)
            resid_all.append(yo - yr)
            yret_all.append(yr)
        out["wn"] = np.concatenate(wn_all)
        out["resid"] = np.concatenate(resid_all)
        out["yret"] = np.concatenate(yret_all)
        out["win_bounds"] = [(w.label, int(w.n_channels)) for w in inst.windows]
    return out
