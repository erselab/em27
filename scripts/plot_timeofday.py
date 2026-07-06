"""Diagnose the time-of-day / broadband-tilt residual mode (EOF2).

Correlates the retrieved bulk nuisances (T_offset, H2O scale) against UT hour and
against the per-band PC2 (the tilt amplitude), to see whether the tilt tracks
temperature or water.  Also shows real surface T (PROFFAST gndT) vs UT for context.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent
import sys; sys.path.insert(0, str(PROJ))
from em27gert.readers import read_invparms


def band_pc2(npz):
    """Per-band PC2 (2nd principal-component amplitude), keyed by spectrum name."""
    labels = [str(x) for x in npz["win_labels"]]; nch = npz["win_nchan"]
    bnd = np.concatenate([[0], np.cumsum(nch)])
    specs = [str(s) for s in npz["spectrum"]]
    out = {}
    for j, lab in enumerate(labels):
        Rb = npz["resid"][:, bnd[j]:bnd[j+1]]
        U, S, Vt = np.linalg.svd(Rb - Rb.mean(0), full_matrices=False)
        out[lab] = dict(zip(specs, (U[:, 1] * S[1])))     # PC2 per sounding
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(PROJ / "data/ils_physical_fixed_results.csv"))
    ap.add_argument("--npz", default=str(PROJ / "data/ils_physical_fixed_resid.npz"))
    ap.add_argument("--out", default=str(PROJ / "figures/em27_timeofday.png"))
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    inv = read_invparms(PROJ / "data/GSFC_SN245_260406/comb_invparms_GSFC_SN245_260406-260406.csv")
    df["gndT"] = df["spectrum"].map(inv["gndT"].astype(float))
    pc2 = band_pc2(np.load(args.npz, allow_pickle=True))
    for lab in ("XCH4", "XCO"):
        df[f"pc2_{lab}"] = df["spectrum"].map(pc2[lab])
    df = df.dropna(subset=["t_offset", "h2o_scale", "pc2_XCH4"]).sort_values("ut_h")

    def c(a, b):
        m = np.isfinite(df[a]) & np.isfinite(df[b])
        return np.corrcoef(df[a][m], df[b][m])[0, 1]

    print("correlations:")
    print(f"  T_offset vs UT      : {c('t_offset','ut_h'):+.2f}     H2O_scale vs UT   : {c('h2o_scale','ut_h'):+.2f}")
    print(f"  T_offset vs gndT     : {c('t_offset','gndT'):+.2f}     H2O_scale vs airmass:{c('h2o_scale','airmass'):+.2f}")
    print(f"  PC2(XCH4) vs T_offset: {c('pc2_XCH4','t_offset'):+.2f}   PC2(XCH4) vs H2O  : {c('pc2_XCH4','h2o_scale'):+.2f}   PC2(XCH4) vs UT: {c('pc2_XCH4','ut_h'):+.2f}")
    print(f"  PC2(XCO)  vs T_offset: {c('pc2_XCO','t_offset'):+.2f}   PC2(XCO)  vs H2O  : {c('pc2_XCO','h2o_scale'):+.2f}   PC2(XCO)  vs UT: {c('pc2_XCO','ut_h'):+.2f}")

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    ax[0, 0].scatter(df.ut_h, df.t_offset, c=df.airmass, cmap="plasma", s=16); ax[0, 0].set_title(f"T_offset vs UT (r={c('t_offset','ut_h'):+.2f})"); ax[0, 0].set_ylabel("T_offset (K)")
    ax[0, 1].scatter(df.ut_h, df.h2o_scale, c=df.airmass, cmap="plasma", s=16); ax[0, 1].set_title(f"H2O scale vs UT (r={c('h2o_scale','ut_h'):+.2f})"); ax[0, 1].set_ylabel("H2O scale")
    ax[0, 2].scatter(df.ut_h, df.gndT, c=df.airmass, cmap="plasma", s=16); ax[0, 2].set_title("surface T (gndT) vs UT"); ax[0, 2].set_ylabel("gndT (K)")
    ax[1, 0].scatter(df.t_offset, df.pc2_XCH4, c=df.ut_h, cmap="viridis", s=16); ax[1, 0].set_title(f"PC2(XCH4) vs T_offset (r={c('pc2_XCH4','t_offset'):+.2f})"); ax[1, 0].set_xlabel("T_offset (K)")
    ax[1, 1].scatter(df.h2o_scale, df.pc2_XCH4, c=df.ut_h, cmap="viridis", s=16); ax[1, 1].set_title(f"PC2(XCH4) vs H2O scale (r={c('pc2_XCH4','h2o_scale'):+.2f})"); ax[1, 1].set_xlabel("H2O scale")
    ax[1, 2].scatter(df.ut_h, df.pc2_XCH4, c=df.airmass, cmap="plasma", s=16); ax[1, 2].set_title(f"PC2(XCH4) vs UT (r={c('pc2_XCH4','ut_h'):+.2f})"); ax[1, 2].set_xlabel("UT hour")
    for a in ax.ravel(): a.grid(alpha=0.3)
    fig.suptitle("Time-of-day residual mode: what drives EOF2 (broadband tilt)?")
    fig.tight_layout(); fig.savefig(args.out, dpi=120, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
