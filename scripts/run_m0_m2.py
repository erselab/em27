#!/usr/bin/env python3
"""M0–M2 driver: ingest → align → EM27 instrument → forward-vs-real (open loop).

Run from anywhere; pass the gert repo path for the ABSCO/solar tables:
    python scripts/run_m0_m2.py --gert /path/to/gert
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from em27gert.readers import (
    read_spectrum, read_invparms, read_ground_pressure, read_ils_list,
    nearest_map,
)
from em27gert.scene import map_to_atmosphere
from em27gert.instrument import build_em27_instrument

from gert.absco import ABSCOTable
from gert.solar import SolarSpectrum
from gert.geometry import Geometry
from gert.forward_model import ForwardModel
from gert.rt_solver import TransmissionSolver

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data" / "GSFC_SN245_260406"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", required=True, help="path to the gert repo (for absco/solar)")
    args = ap.parse_args()
    gert_root = Path(args.gert)
    absco_path = gert_root / "input/absco/absco.h5"
    solar_path = gert_root / "input/solar/solar.h5"

    # ── M0: ingest ────────────────────────────────────────────────────────
    inv = read_invparms(DATA / "comb_invparms_GSFC_SN245_260406-260406.csv")
    gp = read_ground_pressure(DATA / "pressure/b33_20260406.csv")
    ils = read_ils_list(DATA / "ils_list.csv")
    print(f"[M0] invparms: {len(inv)} spectra | ground-P rows: {len(gp)} | "
          f"ILS ME1={ils['ME1'].iloc[0]}")

    # pick a clean SN (near-IR) sounding: lowest rms among near-noon spectra
    sn = inv[inv.index.str.contains("SN")].copy()
    sn = sn[(sn["appSZA"] < 55)]                       # decent airmass
    pick = sn["job01_rms"].astype(float).idxmin()
    row = inv.loc[pick]
    spec = read_spectrum(DATA / "260406_spectra/cal" / pick)
    print(f"[M0] picked {pick}: SZA={spec['sza_deg']:.2f}  "
          f"rms={row['job01_rms']}  ν=[{spec['wn'].min():.0f},{spec['wn'].max():.0f}]")

    # merge the SM (extended) channel so the XCO window (~4210 cm-1) is covered
    sm_path = DATA / "260406_spectra/cal" / pick.replace("SN", "SM")
    meas_wn, meas_sp = spec["wn"], spec["spectrum"]
    if sm_path.exists():
        sm = read_spectrum(sm_path)
        order = np.argsort(np.concatenate([sm["wn"], spec["wn"]]))
        meas_wn = np.concatenate([sm["wn"], spec["wn"]])[order]
        meas_sp = np.concatenate([sm["spectrum"], spec["spectrum"]])[order]
        print(f"[M0] + SM channel {sm_path.name}: "
              f"ν=[{sm['wn'].min():.0f},{sm['wn'].max():.0f}]")

    # nearest .map prior; surface pressure from PROFFAST gndP (hPa→Pa)
    mp = nearest_map(DATA / "map", spec["metadata"]["time_ut_h"], spec["metadata"]["date"])
    p_sfc = float(row["gndP"]) * 100.0
    atm = map_to_atmosphere(mp, p_surface_pa=p_sfc)
    print(f"[M0] prior map: {mp.name}  ({len(atm.p_levels)} levels)  gndP={row['gndP']} hPa")

    # alignment check: GERT prior column vs PROFFAST retrieved
    xco2_prior = atm.column_xgas("co2") * 1e6
    xch4_prior = atm.column_xgas("ch4") * 1e9
    print(f"[M0] prior  XCO2={xco2_prior:8.2f} ppm   XCH4={xch4_prior:8.1f} ppb")
    print(f"[M0] PROFFAST XCO2={float(row['XCO2']):8.2f} ppm   "
          f"XCH4={float(row['XCH4'])*1000:8.1f} ppb   XCO={float(row['XCO']):.3f}")

    # ── M1: EM27 instrument ───────────────────────────────────────────────
    opd = float(spec["metadata"].get("opd_max_cm", 1.8))
    inst = build_em27_instrument(opd_cm=opd)
    print(f"[M1] instrument: OPD={opd} cm, windows="
          f"{[w.label for w in inst.windows]}")

    # ── M2: forward sim vs real spectrum (open loop) ──────────────────────
    absco = ABSCOTable.load_all(str(absco_path))
    solar = SolarSpectrum.load(str(solar_path))
    geo = Geometry(sza=spec["sza_deg"], vza=0.0, raa=0.0,
                   observer_altitude=spec["metadata"]["altitude_km"])
    fm = ForwardModel(atm, absco, inst, geo,
                      solver=TransmissionSolver(jacobians=False),
                      solar_spectrum=solar)
    res = fm.run(albedo=np.ones(len(inst.windows)))
    y = res.y

    # Per-window residual after removing a low-order continuum/gain (deg-2),
    # matching what PROFFAST does with its baseline (bsl) polynomial and what
    # the GERT retrieval does with solar_gain + slope.  This isolates spectral
    # line-shape fidelity (spectroscopy + ILS + solar) from broadband gain.
    print(f"[M2] forward y: {len(y)} channels; per-window residual "
          f"(continuum/gain removed):")
    off = 0
    for w in inst.windows:
        n = w.n_channels
        wn_c = w.wn_instrument                       # wavenumber order (ascending)
        # ForwardModel.R_band is in *wavelength* order (R_wn[::-1]); reverse the
        # slice back to wavenumber order so it aligns with wn_instrument.
        yg = y[off:off + n][::-1]; off += n
        if wn_c.min() < meas_wn.min() or wn_c.max() > meas_wn.max():
            print(f"   {w.label:5s} [{wn_c.min():.1f}-{wn_c.max():.1f}] "
                  f"not covered by the spectra — skipped")
            continue
        ym = np.interp(wn_c, meas_wn, meas_sp)
        # least-squares fit  ym ≈ (a0 + a1·x + a2·x²)·yg   (multiplicative gain)
        x = (wn_c - wn_c.mean()) / (np.ptp(wn_c) / 2.0)
        A = np.vstack([yg, yg * x, yg * x**2]).T
        coef, *_ = np.linalg.lstsq(A, ym, rcond=None)
        model = A @ coef
        resid = (ym - model) / np.nanmedian(ym)        # fractional
        print(f"   {w.label:5s} [{wn_c.min():.1f}-{wn_c.max():.1f}] "
              f"n={n:4d}  RMS={np.sqrt(np.nanmean(resid**2))*100:5.2f}%  "
              f"max|resid|={np.nanmax(np.abs(resid))*100:5.2f}%")
    print("[M2] done — residual quantifies spectroscopy+ILS+solar fidelity "
          "(target: few %, refined by ILS ME/PE + Doppler in M3/M5).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())