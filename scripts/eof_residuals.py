"""Per-band EOF/PCA decomposition of retrieval spectral residuals.

Each spectral window is decomposed INDEPENDENTLY — its own SVD on its own
channels — so the modes and the variance fractions are band-specific.  (A single
SVD over the concatenated all-band residual, as before, mixes the windows and
reports one shared variance spectrum, which is not what we want when the bands
carry different physics.)

Reads a residual npz (run_m4_residuals.py / run_ils_experiment.py) and writes
`<prefix>_eofs.png` (mean + first three modes per band) and `<prefix>_pcs.png`
(per-band scree + PC1/PC2 vs airmass).

Usage:
    python3 scripts/eof_residuals.py --npz data/ils_physical_fixed_resid.npz \
        --prefix em27_alliso_physfixed
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent
FIGS = PROJ / "figures"; FIGS.mkdir(exist_ok=True)

ap = argparse.ArgumentParser()
ap.add_argument("--npz", default=str(PROJ / "data/m4_residuals.npz"))
ap.add_argument("--prefix", default="em27_m4_residual")
args = ap.parse_args()

Z = np.load(args.npz, allow_pickle=True)
wn, R = Z["wn"], Z["resid"]
labels = [str(x) for x in Z["win_labels"]]; nchan = Z["win_nchan"]
airmass, ut_h = Z["airmass"], Z["ut_h"]
bnd = np.concatenate([[0], np.cumsum(nchan)])
n_snd = R.shape[0]; nb = len(labels)
print(f"{n_snd} soundings; per-band independent EOFs for {labels}")

# ── independent SVD per band ───────────────────────────────────────────────
band = {}
for j, lab in enumerate(labels):
    sl = slice(bnd[j], bnd[j + 1])
    w = wn[sl]; Rb = R[:, sl]
    mean = Rb.mean(0)
    U, S, Vt = np.linalg.svd(Rb - mean, full_matrices=False)
    var = S**2 / (S**2).sum()
    PC = U * S                                    # (snd, mode)
    band[lab] = dict(w=w, mean=mean, EOF=Vt, PC=PC, var=var)
    r1 = np.corrcoef(airmass, PC[:, 0])[0, 1]
    r2 = np.corrcoef(airmass, PC[:, 1])[0, 1]
    print(f"  {lab:5}: EOF1 {var[0]*100:5.1f}%  EOF2 {var[1]*100:5.1f}%  EOF3 {var[2]*100:5.1f}%"
          f"   corr(PC1,airmass)={r1:+.2f}  corr(PC2,airmass)={r2:+.2f}")

# ── Figure A: mean + first three modes, per band (each from its own SVD) ────
fig, ax = plt.subplots(4, nb, figsize=(4.6 * nb, 9), squeeze=False)
for j, lab in enumerate(labels):
    b = band[lab]; w = b["w"]
    ax[0, j].plot(w, b["mean"] * 100, lw=0.5, color="k")
    ax[0, j].set_title(f"{lab}: mean resid"); ax[0, j].set_ylabel("%")
    for m in range(3):
        sgn = np.sign(np.nanmean(b["EOF"][m]) or 1.0)
        ax[m + 1, j].plot(w, sgn * b["EOF"][m], lw=0.5)
        ax[m + 1, j].set_ylabel(f"EOF{m+1}\n({b['var'][m]*100:.1f}%)")
        if m == 2:
            ax[m + 1, j].set_xlabel("wavenumber cm$^{-1}$")
fig.suptitle(f"{args.prefix}: per-band residual EOFs (independent SVD per window)")
fig.tight_layout(); fig.savefig(FIGS / f"{args.prefix}_eofs.png", dpi=120)

# ── Figure B: per-band scree + PC1/PC2 vs airmass ──────────────────────────
fig, ax = plt.subplots(nb, 3, figsize=(14, 3.3 * nb), squeeze=False)
for j, lab in enumerate(labels):
    b = band[lab]
    ax[j, 0].plot(np.arange(1, 11), b["var"][:10] * 100, "o-")
    ax[j, 0].set_ylabel(f"{lab}\n% variance"); ax[j, 0].set_xlabel("EOF mode")
    ax[j, 0].set_title(f"{lab} scree")
    for col, m in [(1, 0), (2, 1)]:
        r = np.corrcoef(airmass, b["PC"][:, m])[0, 1]
        sc = ax[j, col].scatter(airmass, b["PC"][:, m], c=ut_h, cmap="viridis", s=16)
        ax[j, col].set_xlabel("airmass"); ax[j, col].set_ylabel(f"PC{m+1}")
        ax[j, col].set_title(f"{lab} PC{m+1} vs airmass  (r={r:+.2f})")
fig.colorbar(sc, ax=ax, label="UT hour", fraction=0.015, pad=0.01)
fig.suptitle(f"{args.prefix}: per-band variance & airmass dependence")
fig.savefig(FIGS / f"{args.prefix}_pcs.png", dpi=120)
print(f"saved figures/{args.prefix}_eofs.png and _pcs.png")
