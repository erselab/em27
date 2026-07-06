"""Time series + co-variation of the retrieved Xgases (CO2, CH4, CO, H2O, XAIR).

Fig 1: 5 stacked panels, GERT vs PROFFAST vs UT hour (absolute values).
Fig 2: (left) z-scored GERT anomalies overlaid — how they co-vary through the day;
       (right) GERT correlation matrix.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent
GASES = [("xco2", "XCO2", "ppm"), ("xch4", "XCH4", "ppb"), ("xco", "XCO", "ppb"),
         ("xh2o", "XH2O", "ppm"), ("xair", "XAIR", "")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(PROJ / "data/xgas_full.csv"))
    args = ap.parse_args()
    df = pd.read_csv(args.csv).sort_values("ut_h").reset_index(drop=True)
    t = df["ut_h"].to_numpy()

    # ── Fig 1: stacked time series, GERT vs PROFFAST ─────────────────────────
    fig, ax = plt.subplots(len(GASES), 1, figsize=(11, 12), sharex=True)
    for a, (g, lab, u) in zip(ax, GASES):
        pf = f"{g}_proffast"
        if pf in df and df[pf].notna().any():
            a.plot(t, df[pf], "-", color="k", lw=1.1, label="PROFFAST")
        a.plot(t, df[f"{g}_gert"], "o", color="tab:red", ms=3.5, label="GERT")
        a.set_ylabel(f"{lab}" + (f" ({u})" if u else "")); a.set_title(lab, fontsize=9)
        a.legend(fontsize=8, loc="best")
    ax[-1].set_xlabel("UT hour")
    fig.suptitle("Retrieved Xgases vs PROFFAST — GSFC SN245 2026-04-06")
    fig.tight_layout(); fig.savefig(PROJ / "figures/em27_xgas_timeseries_all.png", dpi=120, bbox_inches="tight")
    print("wrote figures/em27_xgas_timeseries_all.png")

    # ── Fig 2: z-scored GERT anomalies + correlation matrix ──────────────────
    Z = {}
    for g, lab, _u in GASES:
        v = df[f"{g}_gert"].to_numpy()
        Z[lab] = (v - np.nanmean(v)) / (np.nanstd(v) + 1e-30)
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.5))
    for lab in Z:
        ax[0].plot(t, Z[lab], "-o", ms=3, lw=0.9, label=lab)
    ax[0].axhline(0, color="k", lw=0.4); ax[0].set_xlabel("UT hour")
    ax[0].set_ylabel("z-scored anomaly"); ax[0].set_title("GERT Xgas anomalies (co-variation)")
    ax[0].legend(fontsize=8, ncol=3)

    labs = [g[1] for g in GASES]
    M = np.vstack([Z[l] for l in labs])
    C = np.corrcoef(M)
    im = ax[1].imshow(C, vmin=-1, vmax=1, cmap="RdBu_r")
    ax[1].set_xticks(range(len(labs))); ax[1].set_xticklabels(labs, rotation=45)
    ax[1].set_yticks(range(len(labs))); ax[1].set_yticklabels(labs)
    for i in range(len(labs)):
        for j in range(len(labs)):
            ax[1].text(j, i, f"{C[i,j]:+.2f}", ha="center", va="center",
                       color="white" if abs(C[i, j]) > 0.5 else "black", fontsize=9)
    ax[1].set_title("GERT Xgas correlation matrix")
    fig.colorbar(im, ax=ax[1], fraction=0.046, label="Pearson r")
    fig.suptitle("Xgas co-variation through the day (GERT)")
    fig.tight_layout(); fig.savefig(PROJ / "figures/em27_xgas_covary.png", dpi=120, bbox_inches="tight")
    print("wrote figures/em27_xgas_covary.png")


if __name__ == "__main__":
    main()
