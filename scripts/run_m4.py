"""M4 — loop EM27/SUN SN scans through the GERT retrieval → time series CSV.

Usage:
    PYTHONPATH=. python scripts/run_m4.py --gert ../../gert [--stride 13] [--max N] [--sza-max 80]

Writes ``data/m4_results.csv`` (one row per scan).  Heavy compute lives here so
the notebook only loads/plots the result.
"""
import argparse, sys, time, multiprocessing as mp
from pathlib import Path
import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))

from em27gert.readers import read_invparms
from em27gert.retrieval import retrieve_spectrum

# Populated in the parent before the (fork) pool is created; workers inherit
# these via copy-on-write, so the 2.2 GB ABSCO table is shared, not copied.
_ABSCO = _SOLAR = _INV = _DATA = None


def _worker(pick):
    try:
        return retrieve_spectrum(pick, _INV, _DATA, _ABSCO, _SOLAR)
    except Exception as e:
        return {"spectrum": pick, "error": repr(e)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", default="../../gert", help="path to the gert repo")
    ap.add_argument("--stride", type=int, default=13, help="use every Nth SN scan")
    ap.add_argument("--max", type=int, default=0, help="cap number of scans (0 = no cap)")
    ap.add_argument("--sza-max", type=float, default=80.0)
    ap.add_argument("--workers", type=int, default=8, help="parallel processes (fork)")
    ap.add_argument("--out", default=str(PROJ / "data/m4_results.csv"))
    args = ap.parse_args()

    global _ABSCO, _SOLAR, _INV, _DATA
    gert = Path(args.gert).resolve()
    from gert.absco import ABSCOTable
    from gert.solar import SolarSpectrum
    _ABSCO = ABSCOTable.load_all(str(gert / "input/absco/absco.h5"))
    _SOLAR = SolarSpectrum.load(str(gert / "input/solar/solar.h5"))

    _DATA = PROJ / "data/GSFC_SN245_260406"
    _INV = read_invparms(_DATA / "comb_invparms_GSFC_SN245_260406-260406.csv")
    sn = _INV[_INV.index.str.contains("SN")]
    sn = sn[sn["appSZA"].astype(float) < args.sza_max].sort_values("UTtimeh")
    picks = list(sn.index[::args.stride])
    if args.max:
        picks = picks[:args.max]
    n = len(picks)
    print(f"{n} scans (stride={args.stride}, SZA<{args.sza_max}, workers={args.workers})",
          flush=True)

    rows, t0 = [], time.time()
    ctx = mp.get_context("fork")          # share ABSCO copy-on-write
    with ctx.Pool(args.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(_worker, picks)):
            if "error" in r:
                print(f"[{k+1}/{n}] {r['spectrum']} FAILED: {r['error']}", flush=True)
                continue
            rows.append(r)
            print(f"[{k+1}/{n}] {r['spectrum']} SZA={r['sza']:.1f} chi2={r['chi2']:.2f} "
                  f"dXCO2={r['xco2_gert']-r['xco2_proffast']:+.2f} "
                  f"dXCH4={r['xch4_gert']-r['xch4_proffast']:+.1f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            pd.DataFrame(rows).sort_values("ut_h").to_csv(args.out, index=False)
    print(f"done: {len(rows)} ok -> {args.out} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
