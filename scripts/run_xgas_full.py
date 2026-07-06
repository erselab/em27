"""Full 4-window ensemble (incl. O2) capturing all retrieved Xgases — CO2, CH4,
CO, H2O, XAIR — vs PROFFAST, for the co-variation time series.

Usage: PYTHONPATH=. python scripts/run_xgas_full.py --gert ../../gert --stride 13 --workers 8
"""
import argparse, sys, time, multiprocessing as mp
from pathlib import Path
import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))
from em27gert.readers import read_invparms
from em27gert.retrieval import retrieve_spectrum, GASES_WITH_O2
from em27gert.instrument import EM27_WINDOWS

_ABSCO = _SOLAR = _INV = _DATA = None
COLS = ["spectrum", "ut_h", "airmass", "chi2", "converged",
        "xco2_gert", "xch4_gert", "xco_gert", "xh2o_gert", "xair_gert",
        "xco2_proffast", "xch4_proffast", "xco_proffast", "xh2o_proffast", "xair_proffast"]


def _worker(pick):
    try:
        return retrieve_spectrum(pick, _INV, _DATA, _ABSCO, _SOLAR,
                                 windows=EM27_WINDOWS, gases=GASES_WITH_O2,
                                 res_eff=0.477, retrieve_ils_scale=False)
    except Exception as e:
        return {"spectrum": pick, "error": repr(e)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", default="../../gert")
    ap.add_argument("--stride", type=int, default=13)
    ap.add_argument("--sza-max", type=float, default=80.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=str(PROJ / "data/xgas_full.csv"))
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
    print(f"{len(picks)} scans (4-window incl O2, workers={args.workers})", flush=True)

    rows, t0 = [], time.time()
    with mp.get_context("fork").Pool(args.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(_worker, picks)):
            if "error" in r:
                print(f"[{k+1}] {r['spectrum']} FAIL {r['error']}", flush=True); continue
            rows.append(r)
            if (k + 1) % 10 == 0:
                print(f"[{k+1}/{len(picks)}] {time.time()-t0:.0f}s", flush=True)
    ok = [r for r in rows if r.get("converged")]
    pd.DataFrame(ok)[COLS].sort_values("ut_h").to_csv(args.out, index=False)
    print(f"{len(ok)} ok -> {args.out} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
