"""EM27/SUN instrument model for GERT.

Defines the COCCON/PROFFAST near-IR microwindows and an ILS derived from the
EM27 maximum optical path difference (OPD).  Only windows whose species are
covered by the available ABSCO table are included — the O2 (1.27 um) window is
intentionally omitted (the table's ``o2`` is the 760 nm A-band), so airmass /
dry-air column comes from the measured ground pressure instead.
"""
from __future__ import annotations

import numpy as np

from gert.instrument import ILS, SpectralWindow
from gert.instrument_config import Instrument


def ils_from_me_pe(
    opd_cm: float,
    me1: float, pe1: float,
    me2: float = None, pe2: float = None,
    n_x: int = 4000, n_off: int = 4000,
) -> ILS:
    """Self-apodizing EM27/PROFFAST ILS from modulation efficiency / phase error.

    PROFFAST parameterizes the ILS by the complex modulation along the optical
    path difference x∈[0, L]: efficiency ME(x) and phase error PE(x), each given
    at x = L/2 (ME1/PE1) and x = L (ME2/PE2), with ME(0)=1, PE(0)=0 and linear
    interpolation between.  The ILS is the cosine transform of that modulation::

        ILS(Δν) = ∫₀ᴸ ME(x)·cos(2π x Δν − PE(x)) dx
    """
    L = float(opd_cm)
    if me2 is None: me2 = me1
    if pe2 is None: pe2 = pe1
    x = np.linspace(0.0, L, n_x)
    ME = np.interp(x, [0.0, L / 2, L], [1.0, me1, me2])
    PE = np.interp(x, [0.0, L / 2, L], [0.0, pe1, pe2])
    dnu0 = 1.0 / (2.0 * L)
    off = np.linspace(-12.0 * dnu0, 12.0 * dnu0, n_off)
    phase = 2.0 * np.pi * np.outer(off, x)                 # (n_off, n_x)
    _trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
    resp = _trapz(ME[None, :] * np.cos(phase - PE[None, :]), x, axis=1)
    return ILS(type="tabulated", wn_offsets=off, response=resp)

# Norton–Beer apodization coefficients A(t)=c0+c1(1−t²)+c2(1−t²)², t=x/L∈[0,1]
# (Norton & Beer 1976; same set gert's ILS.from_mopd uses).
_NB_COEFFS = {
    "nb_weak":   (0.384093, 0.087577, 0.528330),
    "nb_medium": (0.152442, 0.136176, 0.711382),
    "nb_strong": (0.045335, 0.554883, 0.399782),
    "boxcar":    (1.0, 0.0, 0.0),
}


def _fwhm_of(off: np.ndarray, resp: np.ndarray) -> float:
    """FWHM [cm⁻¹] of a peak-normalised, single-peaked kernel by half-max crossings."""
    r = resp / resp.max()
    above = np.where(r >= 0.5)[0]
    if len(above) < 2:
        return float("nan")
    # linear-interpolate the two 0.5 crossings around the central lobe
    iL, iR = above[0], above[-1]
    def cross(i0, i1):
        return off[i0] + (0.5 - r[i0]) * (off[i1] - off[i0]) / (r[i1] - r[i0])
    left  = cross(iL - 1, iL) if iL > 0 else off[iL]
    right = cross(iR, iR + 1) if iR + 1 < len(off) else off[iR]
    return float(right - left)


def ils_physical(
    opd_cm: float,
    me1: float, pe1: float,
    me2: float = None, pe2: float = None,
    apod: str = "nb_medium",
    fwhm_cm: float = 0.44,
    n_x: int = 4000, n_off: int = 4001, half_width_cm: float = 4.0,
) -> ILS:
    """Physical EM27 ILS: OPD-sinc × ME/PE self-apodization × Norton–Beer apodization.

    Combines the PROFFAST self-apodizing modulation ME(x)/phase PE(x) (from
    ``ils_list.csv``) with a Norton–Beer apodization along the optical path,
    giving a *physical* line shape that is also smooth/well-fitting (no bare-sinc
    side-lobes).  The kernel is the cosine transform of the along-OPD modulation::

        g(Δν)  = ∫₀ᴸ M(x)·cos(2π x Δν − PE(x)) dx
        g'(Δν) = −∫₀ᴸ M(x)·2π x·sin(2π x Δν − PE(x)) dx     (analytic derivative)

    with ``M(x) = ME(x)·A_NB(x/L)``.  The offset grid is rescaled so the nominal
    FWHM equals ``fwhm_cm`` (default 0.44 cm⁻¹, the EM27 effective resolution);
    a per-band ``ils_scale`` nuisance refines it in the retrieval.  ``g'`` is
    attached to the returned ILS as ``response_deriv`` for analytic dispersion /
    ILS-width Jacobians (no finite differences).
    """
    L = float(opd_cm)
    if me2 is None: me2 = me1
    if pe2 is None: pe2 = pe1
    c0, c1, c2 = _NB_COEFFS[apod.lower()]

    x = np.linspace(0.0, L, n_x)
    ME = np.interp(x, [0.0, L / 2, L], [1.0, me1, me2])
    PE = np.interp(x, [0.0, L / 2, L], [0.0, pe1, pe2])
    t = x / L
    v = 1.0 - t ** 2
    A_nb = c0 + c1 * v + c2 * v ** 2
    M = ME * A_nb                                          # along-OPD modulation

    off = np.linspace(-half_width_cm, half_width_cm, n_off)
    phase = 2.0 * np.pi * np.outer(off, x)                 # (n_off, n_x)
    _trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
    arg = phase - PE[None, :]
    resp  = _trapz(M[None, :] * np.cos(arg), x, axis=1)
    dresp = _trapz(M[None, :] * (-2.0 * np.pi * x[None, :]) * np.sin(arg), x, axis=1)

    # Rescale the offset axis so the nominal FWHM = fwhm_cm.  Stretching offsets
    # by c leaves g unchanged at the scaled offset but divides g' by c
    # (g_new(cΔ)=g(Δ) ⇒ dg_new/d(cΔ)=g'(Δ)/c).
    fwhm0 = _fwhm_of(off, resp)
    c = float(fwhm_cm) / fwhm0
    off_s = off * c
    dresp_s = dresp / c
    # renormalise to unit peak (keep g and g' consistently scaled)
    peak = resp.max()
    return ILS(type="tabulated", wn_offsets=off_s,
               response=resp / peak, response_deriv=dresp_s / peak)


def ils_gaussian(fwhm_cm: float, half_width_cm: float = 4.0, n: int = 2001) -> ILS:
    """A smooth Gaussian effective ILS of the given FWHM (cm⁻¹).

    For the *calibrated* COCCON L1 spectra the effective line shape is well
    described by a smooth ~0.5 cm⁻¹ kernel rather than a bare OPD sinc: the
    dominant broadening is the EM27/SUN finite field of view (≈30 mrad), which
    the bare ME/PE self-apodization (``ils_from_me_pe``) omits.  The effective
    FWHM is identifiable from the spectrum (a χ² scan over FWHM has a sharp
    minimum near 0.5 cm⁻¹), so it can be treated as a retrievable parameter.
    """
    off = np.linspace(-half_width_cm, half_width_cm, n)
    sigma = float(fwhm_cm) / 2.354820045
    resp = np.exp(-0.5 * (off / sigma) ** 2)
    return ILS(type="tabulated", wn_offsets=off, response=resp,
               response_deriv=-(off / sigma ** 2) * resp)   # g'(δ) = −(δ/σ²)·g


# (label, wn_min, wn_max, molecules) — COCCON/PROFFAST EM27 microwindows
# (Frey et al. 2019).  Molecule names match the ABSCO datasets; the O2 window
# uses the 1.27 um ``o2_1p27`` table (not the 760 nm A-band ``o2``).  All four
# windows have dense 0.01 cm-1 coverage after the EM27 absco_spec extension.
EM27_WINDOWS = [
    ("XCO",  4208.7, 4257.3, ["co", "ch4", "h2o", "n2o"]),
    ("XCH4", 5897.0, 6145.0, ["ch4", "co2", "h2o"]),
    ("XCO2", 6173.0, 6390.0, ["co2", "ch4", "h2o"]),
    ("O2",   7765.0, 8005.0, ["o2_1p27", "h2o"]),
]

EM27_OPD_CM = 1.8  # nominal EM27/SUN max optical path difference


def build_em27_instrument(
    opd_cm: float = EM27_OPD_CM,
    apodization: str = "boxcar",
    snr: float = 300.0,
    channels_per_fwhm: int = 3,
    windows=EM27_WINDOWS,
    ils: ILS = None,
    obs_wn=None,
) -> Instrument:
    """Build a `gert.Instrument` for EM27/SUN.

    The ILS is a sinc from the OPD (M1 first pass).  Pass an explicit ``ils``
    (e.g. the self-apodizing ME/PE ILS from :func:`ils_from_me_pe`, built from
    ``ils_list.csv``) to use the measured instrument line shape instead.

    If ``obs_wn`` (the measured wavenumbers) is given, each window's channels are
    set to the observation grid, so the ILS integrates the high-resolution model
    radiance directly onto the measurement points — the standard observation
    operator, with **no resampling of the data**.
    """
    if ils is None:
        try:
            ils = ILS.from_mopd(opd_cm, apodization=apodization)
        except Exception:
            ils = ILS.from_mopd(opd_cm)  # fall back to default apodization
    obs_wn = None if obs_wn is None else np.asarray(obs_wn, dtype=float)
    win = []
    for (label, wmin, wmax, mols) in windows:
        grid = None if obs_wn is None else obs_wn[(obs_wn >= wmin) & (obs_wn <= wmax)]
        win.append(SpectralWindow(
            wn_min=wmin, wn_max=wmax, ils=ils, molecules=mols, label=label,
            hires_spacing=0.01, channels_per_fwhm=channels_per_fwhm, obs_grid=grid,
        ))
    return Instrument(windows=win, snr=snr)