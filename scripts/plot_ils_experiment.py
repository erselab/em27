"""Plot the ILS-shape experiment: physical (NB ME/PE) vs Gaussian, both with
per-band width optimization.  Produces:

  figures/em27_ils_timeseries.png  — ΔXgas vs UT hour, both configs
  figures/em27_ils_residuals.png   — mean residual spectrum per window + retrieved
                                      per-band ILS width vs airmass

and prints a bias/RMS/airmass-slope summary per gas and config.

Usage: python3 scripts/plot_ils_experiment.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parent.parent
CFGS = [("physical", "tab:red", "physical NB ME/PE"),
        ("gaussian", "tab:blue", "Gaussian")]
GASES = [("xco2", "XCO2", "ppm"), ("xch4", "XCH4", "ppb"), ("xco", "XCO", "ppb")]


def load(cfg):
    df = pd.read_csv(PROJ / f"data/ils_{cfg}_results.csv")
    npz = np.load(PROJ / f"data/ils_{cfg}_resid.npz", allow_pickle=True)
    return df, npz


def slope(df, gas):
    d = (df[f"{gas}_gert"] - df[f"{gas}_proffast"]).to_numpy()
    am = df["airmass"].to_numpy()
    ok = np.isfinite(d) & np.isfinite(am)
    s, i = np.polyfit(am[ok], d[ok], 1)
    r = np.corrcoef(am[ok], d[ok])[0, 1]
    return d[ok].mean(), np.sqrt((d[ok]**2).mean()), s, s + i, r


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = {cfg: load(cfg) for cfg, _, _ in CFGS}

    # ── summary table ────────────────────────────────────────────────────────
    print(f"{'gas':6}{'config':12}{'bias':>9}{'RMS':>9}{'slope/am':>10}{'@am=1':>9}{'r':>7}"
          f"{'medFWHM0':>10}{'medFWHM1':>10}")
    for gas, lab, unit in GASES:
        for cfg, _, name in CFGS:
            df = data[cfg][0]
            b, rms, sl, a1, r = slope(df, gas)
            w0 = df.get("ils_fwhm_0", pd.Series([np.nan])).median()
            w1 = df.get("ils_fwhm_1", pd.Series([np.nan])).median()
            print(f"{lab:6}{name:12}{b:9.2f}{rms:9.2f}{sl:10.2f}{a1:9.2f}{r:7.2f}{w0:10.3f}{w1:10.3f}")
        print()

    # ── Figure 1: time series ────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for ax, (gas, lab, unit) in zip(axes, GASES):
        for h in range(0, 24, 3):
            ax.axvline(h, color="grey", lw=0.4, ls=":")
        for cfg, col, name in CFGS:
            df = data[cfg][0].sort_values("ut_h")
            d = df[f"{gas}_gert"] - df[f"{gas}_proffast"]
            ax.plot(df["ut_h"], d, "-o", color=col, ms=3, lw=0.8, label=name)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel(f"Δ{lab} ({unit})"); ax.set_title(lab, fontsize=9)
    axes[0].legend(fontsize=8)
    axes[-1].set_xlabel("UT hour")
    fig.suptitle("ILS shape (both width-optimized) — GERT−PROFFAST time series")
    fig.tight_layout(); fig.savefig(PROJ / "figures/em27_ils_timeseries.png", dpi=120, bbox_inches="tight")
    print("wrote figures/em27_ils_timeseries.png")

    # ── Figure 2: residuals + retrieved width ────────────────────────────────
    npz0 = data["physical"][1]
    labels = list(npz0["win_labels"]); nchan = npz0["win_nchan"]
    bnd = np.concatenate([[0], np.cumsum(nchan)])
    nb = len(labels)
    fig, ax = plt.subplots(3, nb, figsize=(5 * nb, 10), squeeze=False)
    for j, lab in enumerate(labels):
        sl = slice(bnd[j], bnd[j + 1])
        ymax = 0.0
        for cfg, col, name in CFGS:
            npz = data[cfg][1]
            wn = npz["wn"][sl]
            mean_res = npz["resid"][:, sl].mean(0) * 100
            # robust RMS: exclude band-edge channels where the continuum poly ≈ 0
            # inflates a handful of points (they carry large Sy, so don't drive χ²).
            absr = np.abs(npz["resid"][:, sl])
            good = absr < np.percentile(absr, 99.0)
            rms_res = np.sqrt((npz["resid"][:, sl][good]**2).mean()) * 100
            ax[0, j].plot(wn, mean_res, lw=0.5, color=col, label=f"{name} (RMS {rms_res:.2f}%)")
            ymax = max(ymax, np.percentile(np.abs(mean_res), 98))
        ax[0, j].set_ylim(-1.5 * ymax, 1.5 * ymax)     # robust limits (ignore edge spikes)
        ax[0, j].axhline(0, color="k", lw=0.4)
        ax[0, j].set_title(f"{lab}: mean residual"); ax[0, j].set_ylabel("resid (%)")
        ax[0, j].set_xlabel("wavenumber cm$^{-1}$"); ax[0, j].legend(fontsize=7)

        # power spectrum of the (edge-trimmed, detrended) mean residual — reveals
        # the oscillation period.  A peak at the line spacing ⇒ per-line
        # spectroscopy; a fixed unrelated period ⇒ instrumental fringe.
        for cfg, col, name in CFGS:
            npz = data[cfg][1]
            wn = npz["wn"][sl]; mr = npz["resid"][:, sl].mean(0)
            e = max(3, int(0.05 * len(wn)))
            wn, mr = wn[e:-e], mr[e:-e]
            wu = np.arange(wn.min(), wn.max(), 0.05); mu = np.interp(wu, wn, mr)
            k = int(15 / 0.05); mu = mu - np.convolve(mu, np.ones(k)/k, mode="same")
            F = np.abs(np.fft.rfft(mu * np.hanning(len(mu)))); fr = np.fft.rfftfreq(len(mu), 0.05)
            per = 1.0 / np.maximum(fr, 1e-9)
            m = (per > 0.8) & (per < 8)
            ax[1, j].plot(per[m], F[m], color=col, lw=0.9, label=name)
        ax[1, j].set_title(f"{lab}: mean-residual power spectrum")
        ax[1, j].set_xlabel("period (cm$^{-1}$)"); ax[1, j].set_ylabel("|FFT|")
        ax[1, j].legend(fontsize=7)

        # retrieved ILS width vs airmass for this band
        for cfg, col, name in CFGS:
            df = data[cfg][0]
            if f"ils_fwhm_{j}" in df:
                ax[2, j].scatter(df["airmass"], df[f"ils_fwhm_{j}"], s=12, color=col, alpha=0.6, label=name)
        ax[2, j].set_title(f"{lab}: retrieved ILS FWHM"); ax[2, j].set_ylabel("FWHM (cm$^{-1}$)")
        ax[2, j].set_xlabel("airmass"); ax[2, j].legend(fontsize=7)
    fig.suptitle("ILS shape experiment — mean residual per window & retrieved width")
    fig.tight_layout(); fig.savefig(PROJ / "figures/em27_ils_residuals.png", dpi=120, bbox_inches="tight")
    print("wrote figures/em27_ils_residuals.png")


if __name__ == "__main__":
    main()
