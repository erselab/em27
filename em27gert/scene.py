"""Build a `gert.AtmosphericProfile` from a TCCON ``.map`` a-priori file.

Handles the three conversions GERT needs:
  * **ordering** — .map is surface→TOA; GERT wants TOA→surface.
  * **units**    — pressure hPa→Pa; gas columns → mole fraction.
  * **wet→dry**  — .map VMRs are *wet*; GERT gases are dry-air mole fractions,
    ``q_levels`` is specific humidity [kg/kg].
"""
from __future__ import annotations

import numpy as np

from gert.atmosphere import AtmosphericProfile

from .readers import read_map, _MAP_GAS_SCALE

_MW_RATIO = 0.018015 / 0.028964  # M_H2O / M_dry_air ≈ 0.622

# gases we hand to GERT (must exist in the ABSCO table)
_GERT_GASES = ["co2", "ch4", "h2o", "co", "n2o"]


def map_to_atmosphere(map_path, p_surface_pa: float | None = None) -> AtmosphericProfile:
    """Convert a ``.map`` file to a `gert.AtmosphericProfile`.

    Parameters
    ----------
    map_path : path to the ``.map`` file.
    p_surface_pa : optional measured ground pressure [Pa].  If given, the
        bottom level pressure is overridden with this value (the COCCON
        ``gndP``) so the dry-air column matches the real surface.
    """
    df = read_map(map_path)

    # surface→TOA in the file; reverse to TOA→surface for GERT.
    sl = slice(None, None, -1)
    p_hpa = df["pressure"].to_numpy()[sl]
    T = df["temp"].to_numpy()[sl]
    p_pa = p_hpa * 100.0
    if p_surface_pa is not None:
        p_pa[-1] = p_surface_pa  # bottom level = measured ground pressure

    h2o_wet = df["h2o"].to_numpy()[sl] * _MAP_GAS_SCALE["h2o"]
    dry_fac = 1.0 / (1.0 - h2o_wet)  # wet → dry: x_dry = x_wet / (1 - h2o_wet)

    gases: dict[str, np.ndarray] = {}
    for g in _GERT_GASES:
        if g not in df.columns:
            continue
        gases[g] = df[g].to_numpy()[sl] * _MAP_GAS_SCALE[g] * dry_fac

    # specific humidity q [kg/kg] from dry H2O VMR
    w = _MW_RATIO * gases["h2o"]            # mass mixing ratio (kg/kg dry)
    q = w / (1.0 + w)                        # specific humidity

    return AtmosphericProfile(
        p_levels=p_pa, T_levels=T, q_levels=q, gases=gases,
    )