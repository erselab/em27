"""Readers for the EM27/SUN (COCCON / PROFFAST) dataset.

These are *project-specific* loaders — they know the EM27/PROFFAST file formats
and hand clean arrays / DataFrames to the `gert` library.  Nothing here is RT.

Files handled (see ``data/GSFC_SN245_260406/``):
  * ``260406_spectra/cal/*.BIN``  — calibrated solar spectra (via ``proffast_bin``)
  * ``map/*.map``                 — TCCON GINPUT a-priori profiles (wet VMRs)
  * ``comb_invparms_*.csv``       — PROFFAST L2 (XCO2/XCH4/XCO/... per spectrum)
  * ``pressure/*.csv``            — ground-pressure time series
  * ``ils_list.csv``              — per-channel modulation efficiency / phase
"""
from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from .proffast_bin import parse_bin


# ---------------------------------------------------------------------------
# Spectra
# ---------------------------------------------------------------------------
def read_spectrum(path: str | Path) -> dict:
    """Read one PROFFAST ``.BIN`` spectrum.

    Returns the ``parse_bin`` dict (``wn``, ``spectrum``, ``metadata``) plus a
    convenience ``sza_deg`` (= 90 − solar elevation) and the source filename.
    """
    path = Path(path)
    out = parse_bin(path)
    m = out["metadata"]
    out["sza_deg"] = 90.0 - float(m["solar_elevation_deg"])
    out["name"] = path.name
    return out


# ---------------------------------------------------------------------------
# A-priori profiles (.map, TCCON GINPUT)
# ---------------------------------------------------------------------------
# .map columns: Height,Temp,Pressure,Density,h2o,hdo,co2,n2o,co,ch4,hf,o2,gravity
# Units:        km,K,hPa,molec/cm3,parts,parts,ppm,ppb,ppb,ppb,ppt,parts,m/s2
# All gas concentrations are WET mole fractions.
_MAP_GAS_SCALE = {  # → mole fraction (wet)
    "h2o": 1.0, "hdo": 1.0, "o2": 1.0,
    "co2": 1e-6, "n2o": 1e-9, "co": 1e-9, "ch4": 1e-9, "hf": 1e-12,
}
_MAP_TIME_RE = re.compile(r"_(\d{8})(\d{2})Z\.map$")  # ..._YYYYMMDDHHZ.map


def read_map(path: str | Path) -> pd.DataFrame:
    """Read a TCCON ``.map`` file → DataFrame (surface-first, native levels).

    Columns are the raw map columns; gas columns are left in their file units.
    The header row is auto-detected (the line starting with ``Height``).
    """
    path = Path(path)
    lines = path.read_text().splitlines()
    hdr = next(i for i, l in enumerate(lines) if l.strip().startswith("Height"))
    df = pd.read_csv(path, skiprows=hdr, skipinitialspace=True)
    df.columns = [c.strip().lower() for c in df.columns]
    # the row directly under the header holds units (km,K,hPa,...) — drop it
    df = df.iloc[1:].apply(pd.to_numeric, errors="coerce").reset_index(drop=True)
    return df


def map_time(path: str | Path) -> datetime | None:
    """Parse the UT timestamp encoded in a .map filename (``...YYYYMMDDHHZ.map``)."""
    m = _MAP_TIME_RE.search(str(path))
    if not m:
        return None
    return datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H")


def nearest_map(map_dir: str | Path, ut_hour: float, date: str) -> Path:
    """Pick the .map whose UT hour is closest to ``ut_hour`` on ``date`` (ISO)."""
    map_dir = Path(map_dir)
    target = datetime.strptime(date, "%Y-%m-%d").replace(
        hour=int(ut_hour) % 24
    )
    maps = sorted(map_dir.glob("*.map"))
    if not maps:
        raise FileNotFoundError(f"no .map files in {map_dir}")
    return min(maps, key=lambda p: abs(((map_time(p) or target) - target).total_seconds()))


# ---------------------------------------------------------------------------
# PROFFAST L2 (comb_invparms)
# ---------------------------------------------------------------------------
def read_invparms(path: str | Path) -> pd.DataFrame:
    """Read the PROFFAST ``comb_invparms`` L2 file, indexed by spectrum name.

    Columns are whitespace-stripped; ``spectrum`` (the ``.BIN`` filename) is the
    index so a spectrum can be joined directly to its retrieved Xgas.
    """
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    df["spectrum"] = df["spectrum"].astype(str).str.strip()
    return df.set_index("spectrum")


# ---------------------------------------------------------------------------
# Ground pressure + ILS list
# ---------------------------------------------------------------------------
def read_ground_pressure(path: str | Path) -> pd.DataFrame:
    """Read the ground-pressure time series (PROFFAST ``b*_YYYYMMDD.csv``)."""
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    return df


def read_ils_list(path: str | Path) -> pd.DataFrame:
    """Read ``ils_list.csv`` (per-channel modulation efficiency ME / phase PE)."""
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    return df