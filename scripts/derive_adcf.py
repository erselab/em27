"""Derive an air-mass-dependent correction (ADCF) for GERT Xgas from one day, the
TCCON way: the airmass artifact is a function of SZA/airmass only (symmetric about
solar noon), while the true Xgas is smooth in time.  Fit
    X_gert(t) = [smooth diurnal poly in time] + beta·(airmass − airmass_mean)
so beta is the airmass artifact — derived from GERT ALONE (no PROFFAST).  Cross-
check beta against the GERT−PROFFAST airmass slope (PROFFAST is already airmass-
corrected, so that slope IS the artifact).  Apply X_corr = X_gert − beta·(am−mean)
and plot before/after.  (ADCF removes the airmass SLOPE; a constant offset remains
— that is the AICF's job, which needs external in-situ data.)
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent
GASES = [("xco2", "XCO2", "ppm"), ("xch4", "XCH4", "ppb")]


def fit_adcf(df, gas, tdeg=2):
    """beta = airmass-artifact slope from GERT alone; smooth time trend of deg tdeg."""
    x = df[f"{gas}_gert"].to_numpy(); am = df["airmass"].to_numpy(); t = df["ut_h"].to_numpy()
    u = (t - t.mean()) / t.std()
    A = np.vstack([u**k for k in range(tdeg + 1)] + [am - am.mean()]).T
    coef, *_ = np.linalg.lstsq(A, x, rcond=None)
    beta = coef[-1]                                        # ppm (or ppb) per airmass
    corr = beta * (am - am.mean())                        # mean-zero artifact
    return beta, corr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(PROJ / "data/ils_physical_fixed_results.csv"))
    args = ap.parse_args()
    df = pd.read_csv(args.csv).sort_values("ut_h").reset_index(drop=True)

    def slope(a, b):
        s, i = np.polyfit(a, b, 1); r = np.corrcoef(a, b)[0, 1]; return s, s * 0 + i, r

    fig1, ax1 = plt.subplots(1, len(GASES), figsize=(6 * len(GASES), 5))
    fig2, ax2 = plt.subplots(len(GASES), 1, figsize=(11, 4 * len(GASES)), sharex=True)
    print(f"{'gas':6}{'β_indep':>10}{'β_vsPROF':>10}{'ADCF %':>9}"
          f"{'slope0':>9}{'slopeC':>9}{'bias0':>8}{'biasC':>8}{'RMS0':>7}{'RMSC':>7}")
    for k, (gas, lab, u) in enumerate(GASES):
        beta, _ = fit_adcf(df, gas)                        # naive single-day fit (degenerate)
        xg = df[f"{gas}_gert"].to_numpy(); xp = df[f"{gas}_proffast"].to_numpy()
        am = df["airmass"].to_numpy(); t = df["ut_h"].to_numpy()
        beta_prof = np.polyfit(am, xg - xp, 1)[0]          # artifact vs PROFFAST (the valid one)
        corr = beta_prof * (am - am.mean())                # airmass correction (mean-zero)
        xc = xg - corr                                     # ADCF-corrected
        d0, dc = xg - xp, xc - xp
        s0 = np.polyfit(am, d0, 1)[0]; sc = np.polyfit(am, dc, 1)[0]
        print(f"{lab:6}{beta:>10.2f}{beta_prof:>10.2f}{100*beta/xg.mean():>9.3f}"
              f"{s0:>9.2f}{sc:>9.2f}{d0.mean():>8.2f}{dc.mean():>8.2f}"
              f"{np.sqrt((d0**2).mean()):>7.2f}{np.sqrt((dc**2).mean()):>7.2f}")

        # Fig 1: ΔXgas vs airmass, before vs after
        xs = np.linspace(am.min(), am.max(), 50)
        ax1[k].scatter(am, d0, s=16, c="tab:red", alpha=.6, label=f"before (slope {s0:+.2f})")
        ax1[k].plot(xs, np.polyval(np.polyfit(am, d0, 1), xs), "r-", lw=1)
        ax1[k].scatter(am, dc, s=16, c="tab:blue", alpha=.6, label=f"after ADCF (slope {sc:+.2f})")
        ax1[k].plot(xs, np.polyval(np.polyfit(am, dc, 1), xs), "b-", lw=1)
        ax1[k].axhline(0, color="k", lw=.5); ax1[k].set_title(f"{lab}: Δ(GERT−PROFFAST) vs airmass")
        ax1[k].set_xlabel("airmass"); ax1[k].set_ylabel(f"Δ{lab} ({u})"); ax1[k].legend(fontsize=8)

        # Fig 2: time series
        ax2[k].plot(t, xp, "k-", lw=1.1, label="PROFFAST")
        ax2[k].plot(t, xg, "o", c="tab:red", ms=3, label="GERT raw")
        ax2[k].plot(t, xc, "o", c="tab:blue", ms=3, label="GERT + ADCF")
        ax2[k].set_ylabel(f"{lab} ({u})"); ax2[k].set_title(lab, fontsize=9); ax2[k].legend(fontsize=8)
    ax2[-1].set_xlabel("UT hour")
    fig1.suptitle("Airmass-dependent correction (ADCF) — Δ vs airmass, before/after")
    fig2.suptitle("Xgas vs PROFFAST — raw vs ADCF-corrected (ADCF derived from GERT alone)")
    fig1.tight_layout(); fig2.tight_layout()
    fig1.savefig(PROJ / "figures/em27_adcf_airmass.png", dpi=120, bbox_inches="tight")
    fig2.savefig(PROJ / "figures/em27_adcf_timeseries.png", dpi=120, bbox_inches="tight")
    print("\nwrote figures/em27_adcf_{airmass,timeseries}.png")


if __name__ == "__main__":
    main()
