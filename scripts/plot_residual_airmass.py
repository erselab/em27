"""How the retrieval residual changes with airmass, from a residual npz
(run_ils_experiment.py / run_m4_residuals.py — residuals normalized to continuum).

Top row:    per-sounding residual RMS vs airmass (does the residual grow with path?).
Bottom row: mean residual spectrum, low-airmass vs high-airmass tercile
            (does the spectral pattern intensify with airmass?).

Usage:
    python3 scripts/plot_residual_airmass.py --npz data/ils_physical_fixed_resid.npz \
        --out figures/em27_alliso_residual_airmass.png
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=str(PROJ / "data/ils_physical_fixed_resid.npz"))
    ap.add_argument("--out", default=str(PROJ / "figures/em27_alliso_residual_airmass.png"))
    ap.add_argument("--title", default="all-iso, physical ILS (fixed 0.477)")
    args = ap.parse_args()

    Z = np.load(args.npz, allow_pickle=True)
    wn, R = Z["wn"], Z["resid"]                       # R: (n_snd, n_chan), fractional
    labels = list(Z["win_labels"]); nchan = Z["win_nchan"]
    airmass, ut_h = Z["airmass"], Z["ut_h"]
    bnd = np.concatenate([[0], np.cumsum(nchan)])
    nb = len(labels)
    aml, amh = np.percentile(airmass, 33), np.percentile(airmass, 67)
    lo, hi = airmass <= aml, airmass >= amh
    xs = np.linspace(airmass.min(), airmass.max(), 50)

    fig, ax = plt.subplots(2, nb, figsize=(5 * nb, 8), squeeze=False)
    for j, lab in enumerate(labels):
        sl = slice(bnd[j], bnd[j + 1]); w = wn[sl]; Rb = R[:, sl]
        # robust per-sounding RMS: drop the top-1% |resid| channels (band edges)
        cut = np.percentile(np.abs(Rb), 99, axis=1, keepdims=True)
        rms = np.sqrt(np.nanmean(np.where(np.abs(Rb) < cut, Rb, np.nan) ** 2, axis=1)) * 100

        s, i = np.polyfit(airmass, rms, 1)
        r = np.corrcoef(airmass, rms)[0, 1]
        sc = ax[0, j].scatter(airmass, rms, s=14, c=ut_h, cmap="viridis")
        ax[0, j].plot(xs, i + s * xs, "r-", lw=1.3,
                      label=f"slope {s:+.2f} %/airmass\nr={r:.2f}")
        ax[0, j].set_title(f"{lab}: residual RMS vs airmass")
        ax[0, j].set_xlabel("airmass"); ax[0, j].set_ylabel("residual RMS (%)")
        ax[0, j].legend(fontsize=8)

        mlo, mhi = Rb[lo].mean(0) * 100, Rb[hi].mean(0) * 100
        ymax = 1.5 * max(np.percentile(np.abs(mlo), 98), np.percentile(np.abs(mhi), 98))
        ax[1, j].plot(w, mlo, lw=0.5, color="tab:green", label=f"airmass ≤ {aml:.2f}")
        ax[1, j].plot(w, mhi, lw=0.5, color="tab:purple", label=f"airmass ≥ {amh:.2f}")
        ax[1, j].axhline(0, color="k", lw=0.4); ax[1, j].set_ylim(-ymax, ymax)
        ax[1, j].set_title(f"{lab}: mean residual, low vs high airmass")
        ax[1, j].set_xlabel("wavenumber cm$^{-1}$"); ax[1, j].set_ylabel("resid (%)")
        ax[1, j].legend(fontsize=8)

    fig.colorbar(sc, ax=ax[0, :], label="UT hour", fraction=0.02, pad=0.01)
    fig.suptitle(f"Residual vs airmass — {args.title}")
    fig.savefig(args.out, dpi=120, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
