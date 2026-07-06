"""Ensemble comparison of two ILS shapes, BOTH with per-band width optimization.

  physical  — NB-apodized self-apodizing ME/PE ILS (new default) + retrieved ils_scale_{b}
  gaussian  — Gaussian ILS + retrieved ils_scale_{b}

For each config, loops the stride-13 SN ensemble, capturing retrieved Xgas, the
retrieved per-band ILS widths, and the per-channel residual (y_obs − y_ret).
Writes data/ils_<config>_results.csv and data/ils_<config>_resid.npz.

Usage: PYTHONPATH=. python scripts/run_ils_experiment.py --gert ../../gert --stride 13 --workers 8
"""
import argparse, sys, time, multiprocessing as mp
from pathlib import Path
import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))
from em27gert.readers import read_invparms
from em27gert.retrieval import retrieve_spectrum
from em27gert.instrument import ils_gaussian

_ABSCO = _SOLAR = _INV = _DATA = _CONFIG = None
_RES_EFF = 0.477            # fixed physical-ILS FWHM for the *_fixed configs


def _worker(pick):
    try:
        if _CONFIG == "physical_fixed":            # physical ILS, width FROZEN at _RES_EFF
            kw = dict(return_spectra=True, retrieve_ils_scale=False, res_eff=_RES_EFF)
        elif _CONFIG == "gaussian":                # Gaussian, width optimized
            kw = dict(return_spectra=True, retrieve_ils_scale=True, ils=ils_gaussian(0.44))
        else:                                      # physical, width optimized (per band)
            kw = dict(return_spectra=True, retrieve_ils_scale=True)
        return retrieve_spectrum(pick, _INV, _DATA, _ABSCO, _SOLAR, **kw)
    except Exception as e:
        return {"spectrum": pick, "error": repr(e)}


def run_config(config, picks, out_csv, out_npz, workers):
    global _CONFIG
    _CONFIG = config
    rows, t0 = [], time.time()
    ctx = mp.get_context("fork")
    with ctx.Pool(workers) as pool:
        for k, r in enumerate(pool.imap_unordered(_worker, picks)):
            if "error" in r:
                print(f"  [{config}] {r['spectrum']} FAIL {r['error']}", flush=True); continue
            rows.append(r)
            w = [r.get(f"ils_fwhm_{b}", np.nan) for b in range(2)]
            print(f"  [{config} {k+1}/{len(picks)}] {r['spectrum']} chi2={r['chi2']:.2f} "
                  f"fwhm=({w[0]:.3f},{w[1]:.3f}) dXCO2={r['xco2_gert']-r['xco2_proffast']:+.1f} "
                  f"dXCH4={r['xch4_gert']-r['xch4_proffast']:+.1f} ({time.time()-t0:.0f}s)", flush=True)

    ok = [r for r in rows if r.get("converged") and r["chi2"] < 10]
    # scalar table
    cols = ["spectrum", "ut_h", "sza", "airmass", "chi2", "t_offset", "h2o_scale",
            "xco2_gert", "xch4_gert", "xco_gert", "xco2_proffast", "xch4_proffast", "xco_proffast"]
    df = pd.DataFrame(ok)[cols + [c for c in ok[0] if c.startswith(("ils_scale_", "ils_fwhm_"))]]
    df.sort_values("ut_h").to_csv(out_csv, index=False)
    # residual matrix — normalized to the per-window continuum so it reads as a
    # fraction of signal (y_obs/y_ret are on the radiance scale ~O(200), not ~1).
    nchan = ok[0]["win_nchan"] if "win_nchan" in ok[0] else \
        np.array([b[1] for b in ok[0]["win_bounds"]])
    bnd = np.concatenate([[0], np.cumsum(nchan)])
    def _frac(r):
        yo = r["resid"] + r["yret"]           # y_obs on radiance scale
        out = np.empty_like(r["resid"])
        for j in range(len(nchan)):
            sl = slice(bnd[j], bnd[j + 1])
            cont = np.median(yo[sl][yo[sl] > np.percentile(yo[sl], 70)])  # continuum level
            out[sl] = r["resid"][sl] / cont
        return out
    R = np.vstack([_frac(r) for r in ok])
    np.savez(out_npz, wn=ok[0]["wn"], resid=R,
             win_labels=np.array([b[0] for b in ok[0]["win_bounds"]]),
             win_nchan=np.array([b[1] for b in ok[0]["win_bounds"]]),
             airmass=np.array([r["airmass"] for r in ok]),
             ut_h=np.array([r["ut_h"] for r in ok]),
             chi2=np.array([r["chi2"] for r in ok]),
             spectrum=np.array([r["spectrum"] for r in ok]))
    print(f"  [{config}] {len(ok)} ok -> {out_csv.name}, {out_npz.name} ({time.time()-t0:.0f}s)\n", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", default="../../gert")
    ap.add_argument("--stride", type=int, default=13)
    ap.add_argument("--max", type=int, default=0)
    ap.add_argument("--sza-max", type=float, default=80.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--configs", nargs="+", default=["physical", "gaussian"],
                    help="physical | physical_fixed | gaussian")
    ap.add_argument("--res-eff", type=float, default=0.477,
                    help="fixed physical-ILS FWHM for the physical_fixed config")
    args = ap.parse_args()

    global _ABSCO, _SOLAR, _INV, _DATA, _RES_EFF
    _RES_EFF = args.res_eff
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
    print(f"{len(picks)} scans, configs={args.configs}, workers={args.workers}", flush=True)

    for cfg in args.configs:
        run_config(cfg, picks, PROJ / f"data/ils_{cfg}_results.csv",
                   PROJ / f"data/ils_{cfg}_resid.npz", args.workers)


if __name__ == "__main__":
    main()
