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

from .readers import read_spectrum, nearest_map
from .scene import map_to_atmosphere
from .instrument import build_em27_instrument, ils_gaussian, EM27_WINDOWS

# Effective ILS resolution determined in M3 (χ² minimum ≈ EM27 30 mrad FOV + OPD).
RES_EFF_CM = 0.50
GASES = ("co2", "ch4", "co", "n2o", "h2o")


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


def _build_y_obs(inst, y0, meas_wn, meas_sp, merr):
    """Measured spectrum on the GERT radiometric scale (per-window deg-2
    continuum), wavelength order, plus the diagonal Sy."""
    yobs_parts, sig_parts, off = [], [], 0
    for w in inst.windows:
        n = w.n_channels
        yg = y0[off:off + n][::-1]               # ascending wn (R_band is wavelength order)
        off += n
        wn = w.wn_instrument
        ym = np.interp(wn, meas_wn, meas_sp)
        x = (wn - wn.mean()) / (np.ptp(wn) / 2.0)
        A = np.vstack([yg, yg * x, yg * x**2]).T
        coef, *_ = np.linalg.lstsq(A, ym, rcond=None)
        P = coef[0] + coef[1] * x + coef[2] * x**2
        yo = ym / P
        yobs_parts.append(yo[::-1])
        sig_parts.append((np.abs(yo) * merr + 1e-12)[::-1])
    return np.concatenate(yobs_parts), np.concatenate(sig_parts)


def retrieve_spectrum(pick, inv, data_dir, absco, solar, *,
                      res_eff=RES_EFF_CM, merr=0.005, gases=GASES,
                      windows=EM27_WINDOWS[:3], max_iter=12, freeze_gas=None,
                      use_doppler=True, dispersion_order=1, ils=None):
    """Run the M3 retrieval for one ``.BIN`` spectrum.

    Returns a dict with the retrieved Xgas, PROFFAST Xgas, χ², SZA, time, and
    convergence flag.  ``absco`` / ``solar`` are passed in so they load once.

    Ablation knobs (M5): ``use_doppler`` toggles the solar Doppler pre-shift,
    ``dispersion_order=0`` disables the dispersion nuisance, and ``ils`` overrides
    the effective Gaussian ILS (e.g. with ``ils_from_me_pe`` or a bare sinc).
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
    inst = build_em27_instrument(windows=windows,
                                 ils=ils if ils is not None else ils_gaussian(res_eff))
    nw = len(inst.windows)
    fm = ForwardModel(atm, absco, inst, geo,
                      solver=TransmissionSolver(jacobians=True), solar_spectrum=solar)
    y0 = fm.run(albedo=np.zeros(nw), solar_doppler=sd).y
    y_obs, sig = _build_y_obs(inst, y0, meas_wn, meas_sp, merr)

    sv = StateVector.transmission_scaling(
        n_bands=nw, gases=list(gases),
        gas_uncerts={"co": 0.30, "n2o": 0.20, "h2o": 0.30},
        co2_uncert=0.10, ch4_uncert=0.20, solar_gain_uncert=0.20,
        include_dispersion=dispersion_order > 0,
        dispersion_order=max(dispersion_order, 1), dispersion_uncert=0.10)
    sv.freeze("p_scale")
    if freeze_gas:
        for e in sv.elements:
            mol = e.transform[:-6] if e.transform.endswith("_scale") else None
            if mol in freeze_gas:
                e.prior = freeze_gas[mol]; e.value = e.prior
                sv.freeze(e.name)

    res = GERTRetrieval(fm, y_obs, np.diag(1.0 / sig**2), sv,
                        prior_albedo=np.zeros(nw), analytical_jacobians=True,
                        solar_doppler=sd, max_iter=max_iter, verbose=False).run()
    scale = {e.transform[:-6]: res.x_ret[i] for i, e in enumerate(sv.elements)
             if e.transform.endswith("_scale") and e.transform != "p_scale"}

    return {
        "spectrum": pick,
        "ut_h": md["time_ut_h"],
        "sza": spec["sza_deg"],
        "airmass": 1.0 / np.cos(np.radians(spec["sza_deg"])),
        "chi2": res.chisq_reduced,
        "converged": bool(res.converged),
        "xco2_gert": prior["co2"] * scale["co2"],
        "xch4_gert": prior["ch4"] * scale["ch4"],
        "xco_gert":  prior["co"]  * scale["co"],
        "xco2_proffast": float(row["XCO2"]),
        "xch4_proffast": float(row["XCH4"]) * 1000.0,
        "xco_proffast":  float(row["XCO"]) * 1000.0,
    }
