"""Recompute the M4 retrievals capturing the per-channel spectral residuals,
and save them for EOF/PCA analysis (-> data/m4_residuals.npz).

Usage: PYTHONPATH=. python scripts/run_m4_residuals.py --gert ../../gert --stride 13 --workers 8
"""
import argparse, sys, time, multiprocessing as mp
from pathlib import Path
import numpy as np

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))
from em27gert.readers import read_invparms
from em27gert.retrieval import retrieve_spectrum

_ABSCO = _SOLAR = _INV = _DATA = None


def _worker(pick):
    try:
        return retrieve_spectrum(pick, _INV, _DATA, _ABSCO, _SOLAR, return_spectra=True)
    except Exception as e:
        return {"spectrum": pick, "error": repr(e)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", default="../../gert")
    ap.add_argument("--stride", type=int, default=13)
    ap.add_argument("--max", type=int, default=0)
    ap.add_argument("--sza-max", type=float, default=80.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=str(PROJ / "data/m4_residuals.npz"))
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
    print(f"{len(picks)} scans, workers={args.workers}", flush=True)

    rows, t0 = [], time.time()
    ctx = mp.get_context("fork")
    with ctx.Pool(args.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(_worker, picks)):
            if "error" in r:
                print(f"[{k+1}] {r['spectrum']} FAIL {r['error']}", flush=True); continue
            rows.append(r)
            print(f"[{k+1}/{len(picks)}] {r['spectrum']} chi2={r['chi2']:.2f} "
                  f"({time.time()-t0:.0f}s)", flush=True)

    rows = [r for r in rows if r.get("converged") and r["chi2"] < 10]
    wn = rows[0]["wn"]
    R = np.vstack([r["resid"] for r in rows])                 # (n_snd, n_chan)
    np.savez(args.out, wn=wn, resid=R,
             win_labels=np.array([b[0] for b in rows[0]["win_bounds"]]),
             win_nchan=np.array([b[1] for b in rows[0]["win_bounds"]]),
             sza=np.array([r["sza"] for r in rows]),
             ut_h=np.array([r["ut_h"] for r in rows]),
             airmass=np.array([r["airmass"] for r in rows]),
             chi2=np.array([r["chi2"] for r in rows]),
             spectrum=np.array([r["spectrum"] for r in rows]))
    print(f"saved {R.shape} residual matrix -> {args.out} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
