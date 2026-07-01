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
    return ILS(type="tabulated", wn_offsets=off, response=np.exp(-0.5 * (off / sigma) ** 2))


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
) -> Instrument:
    """Build a `gert.Instrument` for EM27/SUN.

    The ILS is a sinc from the OPD (M1 first pass).  Pass an explicit ``ils``
    (e.g. the self-apodizing ME/PE ILS from :func:`ils_from_me_pe`, built from
    ``ils_list.csv``) to use the measured instrument line shape instead.
    """
    if ils is None:
        try:
            ils = ILS.from_mopd(opd_cm, apodization=apodization)
        except Exception:
            ils = ILS.from_mopd(opd_cm)  # fall back to default apodization
    win = [
        SpectralWindow(
            wn_min=wmin, wn_max=wmax, ils=ils, molecules=mols, label=label,
            hires_spacing=0.01, channels_per_fwhm=channels_per_fwhm,
        )
        for (label, wmin, wmax, mols) in windows
    ]
    return Instrument(windows=win, snr=snr)