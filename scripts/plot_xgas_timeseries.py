"""Xgas time series: GERT vs PROFFAST through the day, from an M4-style results CSV."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent
GASES = [("xco2", "XCO2", "ppm"), ("xch4", "XCH4", "ppb"), ("xco", "XCO", "ppb")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(PROJ / "data/ils_physical_fixed_results.csv"))
    ap.add_argument("--out", default=str(PROJ / "figures/em27_alliso_xgas_timeseries.png"))
    ap.add_argument("--title", default="all-iso, physical ILS (fixed 0.477)")
    args = ap.parse_args()

    df = pd.read_csv(args.csv).sort_values("ut_h")
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    for ax, (g, lab, u) in zip(axes, GASES):
        d = df[f"{g}_gert"] - df[f"{g}_proffast"]
        bias, rms = d.mean(), np.sqrt((d**2).mean())
        ax.plot(df["ut_h"], df[f"{g}_proffast"], "-", color="k", lw=1.2, label="PROFFAST")
        ax.plot(df["ut_h"], df[f"{g}_gert"], "o", color="tab:red", ms=4,
                label=f"GERT (bias {bias:+.2f}, RMS {rms:.2f} {u})")
        ax.set_ylabel(f"{lab} ({u})"); ax.set_title(lab, fontsize=9); ax.legend(fontsize=8)
    axes[-1].set_xlabel("UT hour")
    fig.suptitle(f"Xgas vs PROFFAST — {args.title}")
    fig.tight_layout(); fig.savefig(args.out, dpi=120, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
