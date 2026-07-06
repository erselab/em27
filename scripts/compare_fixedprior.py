"""Compare a fixed-single-prior M4 run against the per-scan nearest-map baseline.

For each gas (XCO2/XCH4/XCO): bias, RMS, and airmass slope of (GERT − PROFFAST),
plus how much the single-prior choice moved each scan.

Usage:
    python3 scripts/compare_fixedprior.py \
        --nearest data/m4_results_nearestmap.csv \
        --fixed   data/m4_results_fixedprior.csv
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

GASES = [("xco2", "XCO2", "ppm"), ("xch4", "XCH4", "ppb"), ("xco", "XCO", "ppb")]


def _stats(df, gas):
    d = (df[f"{gas}_gert"] - df[f"{gas}_proffast"]).to_numpy()
    am = df["airmass"].to_numpy()
    ok = np.isfinite(d) & np.isfinite(am)
    d, am = d[ok], am[ok]
    slope, icpt = np.polyfit(am, d, 1)
    r = np.corrcoef(am, d)[0, 1]
    return dict(bias=d.mean(), rms=np.sqrt((d**2).mean()),
                slope=slope, at_am1=slope + icpt, r=r, n=len(d))


def plot_airmass(a, b, out):
    """Δgas vs airmass, nearest-map vs fixed-prior overlaid (M4 Figure 6 style)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, (gas, lab, unit) in zip(axes, GASES):
        for df, name, col in [(a, "nearest-map", "tab:blue"),
                              (b, "fixed-prior", "tab:red")]:
            d = (df[f"{gas}_gert"] - df[f"{gas}_proffast"]).to_numpy()
            am = df["airmass"].to_numpy()
            ok = np.isfinite(d) & np.isfinite(am)
            d, am = d[ok], am[ok]
            slope, icpt = np.polyfit(am, d, 1)
            r = np.corrcoef(am, d)[0, 1]
            xs = np.linspace(am.min(), am.max(), 50)
            ax.scatter(am, d, s=16, color=col, alpha=0.55, edgecolors="none")
            ax.plot(xs, icpt + slope * xs, "-", color=col, lw=1.4,
                    label=f"{name}: {slope:+.2f} {unit}/am, r={r:.2f}")
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_xlabel("airmass (1/cos SZA)")
        ax.set_ylabel(f"Δ{lab} GERT−PROFFAST ({unit})")
        ax.set_title(lab); ax.legend(fontsize=8, loc="best")
    fig.suptitle("Prior-atmosphere sensitivity — bias vs airmass "
                 "(per-scan nearest .map vs one fixed daily prior)")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")


def plot_timeseries(a, b, out):
    """Δgas vs UT hour, nearest-map vs fixed-prior overlaid.  Vertical lines mark
    the 3-hourly .map boundaries where nearest_map switches priors."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a, b = a.sort_values("ut_h"), b.sort_values("ut_h")
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for ax, (gas, lab, unit) in zip(axes, GASES):
        for h in range(0, 24, 3):                     # .map switch boundaries
            ax.axvline(h, color="grey", lw=0.4, ls=":")
        for df, name, col in [(a, "nearest-map", "tab:blue"),
                              (b, "fixed-prior", "tab:red")]:
            d = df[f"{gas}_gert"] - df[f"{gas}_proffast"]
            ax.plot(df["ut_h"], d, "-o", color=col, ms=3, lw=0.8, label=name)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel(f"Δ{lab} ({unit})"); ax.set_title(lab, fontsize=9)
    axes[0].legend(fontsize=8, loc="best")
    axes[-1].set_xlabel("UT hour  (dotted = 3-hourly .map boundaries)")
    fig.suptitle("Prior-atmosphere sensitivity — time series "
                 "(map-switch steps vs one fixed daily prior)")
    fig.tight_layout()
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nearest", default="data/m4_results_nearestmap.csv")
    ap.add_argument("--fixed", default="data/m4_results_fixedprior.csv")
    ap.add_argument("--fig", default="figures/em27_fixedprior_airmass.png")
    ap.add_argument("--fig-ts", default="figures/em27_fixedprior_timeseries.png")
    args = ap.parse_args()

    a = pd.read_csv(args.nearest).query("converged").set_index("spectrum")
    b = pd.read_csv(args.fixed).query("converged").set_index("spectrum")
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]
    print(f"{len(common)} common converged scans "
          f"(airmass {a['airmass'].min():.2f}–{a['airmass'].max():.2f})\n")

    for gas, _lab, unit in GASES:
        sa, sb = _stats(a, gas), _stats(b, gas)
        moved = (b[f"{gas}_gert"] - a[f"{gas}_gert"]).to_numpy()
        print(f"== {gas.upper()} [{unit}] ==")
        print(f"  {'metric':<12}{'nearest-map':>14}{'fixed-prior':>14}{'Δ':>12}")
        for k, lbl in [("bias", "bias"), ("rms", "RMS"),
                       ("slope", "slope/am"), ("at_am1", "@ am=1"), ("r", "r(am)")]:
            print(f"  {lbl:<12}{sa[k]:>14.3f}{sb[k]:>14.3f}{sb[k]-sa[k]:>12.3f}")
        print(f"  per-scan shift from fixing prior: "
              f"mean {np.nanmean(moved):+.3f}, RMS {np.sqrt(np.nanmean(moved**2)):.3f} {unit}\n")

    if args.fig:
        Path(args.fig).parent.mkdir(parents=True, exist_ok=True)
        plot_airmass(a, b, args.fig)
    if args.fig_ts:
        Path(args.fig_ts).parent.mkdir(parents=True, exist_ok=True)
        plot_timeseries(a, b, args.fig_ts)


if __name__ == "__main__":
    main()
