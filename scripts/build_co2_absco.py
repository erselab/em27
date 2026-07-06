#!/usr/bin/env python3
"""Build a CO2-only ABSCO block with a selectable line shape — Step 0 of the
CO2 line-shape experiment (docs/co2_lineshape_experiment.md).

Speedups vs gert/utils/build_absco.py (all numerically identical to it):
  * builds only the CO2 weak-band sub-range (default 6100-6400), not 4750-6400;
  * parallel over (p,T) slices with fork, so the HITRAN line list is loaded
    ONCE and shared copy-on-write (no per-slice db_begin);
  * reuses the existing absco_co2_p / absco_co2_T grids for a clean drop-in.

Tiers:  --line-shape voigt | sdvoigt | ht   (sdvoigt/ht need SD/HT line params;
        see the spec — CO2 must be re-fetched with those parameter groups.)

Usage:
  PYTHONPATH=. python scripts/build_co2_absco.py --gert ../../gert \
      --line-shape voigt --workers 12 --compare
"""
from __future__ import annotations
import argparse, sys, time, multiprocessing as mp
from pathlib import Path
import numpy as np
import h5py

WN_STEP = 0.01     # cm-1  (must match fetch_hitran / build_absco)
WING    = 25.0     # cm-1  Voigt line-wing cutoff

# globals populated in the parent before forking (children inherit via COW)
_HAPI = _TABLE = _WMIN = _WMAX = _SHAPE = None


def _slice(args):
    """Compute one (p,T) cross-section slice using the pre-loaded HITRAN DB."""
    p_pa, T_k, ip, it = args
    kw = dict(SourceTables=[_TABLE], WavenumberRange=[_WMIN, _WMAX],
              WavenumberStep=WN_STEP, WavenumberWing=WING,
              Environment={"T": T_k, "p": p_pa / 101325.0}, HITRAN_units=True)
    fell_back = False
    try:
        if _SHAPE == "voigt":
            _nu, coef = _HAPI.absorptionCoefficient_Voigt(**kw)
        elif _SHAPE == "sdvoigt":                   # quadratic speed-dependent Voigt
            _nu, coef = _HAPI.absorptionCoefficient_SDVoigt(**kw)
        elif _SHAPE == "sdvoigt_lm":                # qSDV + 1st-order (Rosenkranz) line mixing
            _nu, coef = _HAPI.absorptionCoefficient_SDVoigt(**kw, LineMixingRosen=True)
        else:
            raise ValueError(_SHAPE)
    except IndexError:
        # HAPI SDVoigt fails at near-vacuum (Gamma0→0); there the speed-dependence
        # and line-mixing terms vanish, so plain Voigt is the exact p→0 limit.
        _nu, coef = _HAPI.absorptionCoefficient_Voigt(**kw)
        fell_back = True
    return ip, it, np.asarray(coef, dtype=np.float32), fell_back


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", default="../../gert")
    ap.add_argument("--cache-block", default=None,
                    help="HITRAN cache subdir (table = uppercase). Default: "
                         "co2_5880_6400 for voigt, co2_5880_6400_sdv for the SD tiers "
                         "(which carry the speed-dependence + line-mixing params).")
    ap.add_argument("--line-shape", default="voigt",
                    choices=["voigt", "sdvoigt", "sdvoigt_lm"])
    ap.add_argument("--wn-min", type=float, default=6100.0)
    ap.add_argument("--wn-max", type=float, default=6400.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=None, help="output .npz (default: data/co2_absco_<shape>.npz)")
    ap.add_argument("--compare", action="store_true",
                    help="compare to the existing absco.h5 co2 block on the range")
    args = ap.parse_args()

    gert = Path(args.gert).resolve()
    sys.path.insert(0, str(gert))                       # for `import hapi`
    import certifi                                    # Py3.14 needs an explicit CA bundle
    import os as _os; _os.environ.setdefault("SSL_CERT_FILE", certifi.where())
    global _HAPI, _TABLE, _WMIN, _WMAX, _SHAPE
    import hapi
    _HAPI, _WMIN, _WMAX, _SHAPE = hapi, args.wn_min, args.wn_max, args.line_shape
    cache_block = args.cache_block or ("co2_5880_6400" if args.line_shape == "voigt"
                                       else "co2_5880_6400_sdv")
    _TABLE = cache_block.upper()
    subdir = str(gert / "hitran_cache" / cache_block)

    absco_h5 = gert / "input/absco/absco.h5"
    with h5py.File(absco_h5, "r") as f:                 # reuse existing p/T grid
        p_grid = f["absco_co2_p"][:]
        T_grid = f["absco_co2_T"][:]
    n_p, n_T = T_grid.shape

    hapi.db_begin(subdir)                               # load line list ONCE
    n_wn = round((args.wn_max - args.wn_min) / WN_STEP) + 1
    wn = np.linspace(args.wn_min, args.wn_max, n_wn)
    ds = np.zeros((n_p, n_T, n_wn), dtype=np.float32)
    jobs = [(float(p_grid[ip]), float(T_grid[ip, it]), ip, it)
            for ip in range(n_p) for it in range(n_T)]
    print(f"co2 {args.line_shape}  [{args.wn_min:.0f}-{args.wn_max:.0f}]  "
          f"{n_p}×{n_T}={len(jobs)} slices × {n_wn:,} wn  workers={args.workers}",
          flush=True)

    t0 = time.time()
    n_fb = 0
    ctx = mp.get_context("fork")                        # share DB copy-on-write
    with ctx.Pool(args.workers) as pool:
        for k, (ip, it, sig, fb) in enumerate(pool.imap_unordered(_slice, jobs, chunksize=4)):
            m = min(len(sig), n_wn)
            ds[ip, it, :m] = sig[:m]
            n_fb += fb
            if (k + 1) % 200 == 0 or k + 1 == len(jobs):
                el = time.time() - t0
                print(f"  {k+1}/{len(jobs)}  {el:.0f}s  ETA {el/(k+1)*(len(jobs)-k-1):.0f}s",
                      flush=True)
    print(f"built in {(time.time()-t0)/60:.2f} min  "
          f"({n_fb} near-vacuum slices fell back to Voigt)", flush=True)

    out = Path(args.out) if args.out else Path("data") / f"co2_absco_{args.line_shape}.npz"
    np.savez(out, wn=wn, dataset=ds, p=p_grid, T=T_grid,
             line_shape=args.line_shape, cache_block=args.cache_block)
    print(f"saved -> {out}", flush=True)

    if args.compare:
        with h5py.File(absco_h5, "r") as f:
            wn0 = f["absco_co2_wavenumber"][:]
            m = (wn0 >= args.wn_min - 1e-6) & (wn0 <= args.wn_max + 1e-6)
            ds0 = f["absco_co2_dataset"][:, :, m]
        assert ds0.shape == ds.shape, f"grid mismatch {ds0.shape} vs {ds.shape}"
        denom = np.maximum(np.abs(ds0), 1e-30)
        frac = np.abs(ds - ds0) / denom
        big = ds0 > np.percentile(ds0[ds0 > 0], 50)     # where absorption is significant
        print(f"[compare vs absco.h5]  max|Δ/σ|={frac.max():.2e}  "
              f"mean(where signif.)={frac[big].mean():.2e}  "
              f"RMS(all)={np.sqrt(np.mean((ds-ds0)**2))/ds0.max():.2e}", flush=True)


if __name__ == "__main__":
    main()
