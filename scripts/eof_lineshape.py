"""Diagnose EOF1 of the M4 residuals: is it a wavenumber SHIFT (dispersion),
an ILS width error, or a gas AMOUNT error?

For a modelled transmittance T(ν) the three residual signatures are:
  shift  by δν    :  ΔT ≈ -δν · dT/dν            (antisymmetric at line core)
  ILS broadening  :  ΔT ≈ (σ²/2) · d²T/dν²        (symmetric at core)
  amount scale ε  :  ΔT ≈ ε · T·ln T              (symmetric, ∝ line depth)

We build these three basis vectors per window from a modelled spectrum, then
regress EOF1 onto them and report which explains the mode.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ)); GERT = (PROJ / ".." / ".." / "gert")
FIGS = PROJ / "figures"

# --- per-band EOFs from the saved residual matrix ---
Z = np.load(PROJ / "data/m4_residuals.npz", allow_pickle=True)
wn, R = Z["wn"], Z["resid"]
labels, nchan = Z["win_labels"], Z["win_nchan"]
airmass = Z["airmass"]
bounds = np.concatenate([[0], np.cumsum(nchan)])

def band_eof1(j):
    """Leading EOF (spectral pattern), its %variance, and PC1-airmass corr for band j."""
    sl = slice(bounds[j], bounds[j + 1])
    Ab = R[:, sl] - R[:, sl].mean(0)
    U, Sv, Vt = np.linalg.svd(Ab, full_matrices=False)
    var = Sv[0] ** 2 / (Sv ** 2).sum()
    pc1 = U[:, 0] * Sv[0]
    corr = np.corrcoef(airmass, pc1)[0, 1]
    return Vt[0], var, corr

# --- one modelled spectrum on the same grid (line positions are common) ---
from em27gert.readers import read_invparms
from em27gert.retrieval import retrieve_spectrum
from gert.absco import ABSCOTable; from gert.solar import SolarSpectrum
absco = ABSCOTable.load_all(str((GERT / "input/absco/absco.h5").resolve()))
solar = SolarSpectrum.load(str((GERT / "input/solar/solar.h5").resolve()))
DATA = PROJ / "data/GSFC_SN245_260406"
inv = read_invparms(DATA / "comb_invparms_GSFC_SN245_260406-260406.csv")
sn = inv[inv.index.str.contains("SN")]; sn = sn[sn.appSZA < 55]
PICK = sn.job01_rms.astype(float).idxmin()
r = retrieve_spectrum(PICK, inv, DATA, absco, solar, return_spectra=True)
T = r["yret"]
assert np.allclose(r["wn"], wn), "grid mismatch"

from scipy.ndimage import gaussian_filter1d

def bases(w, t):
    """Physical residual signatures from a transmittance model T, via finite
    perturbations (robust) rather than raw high-order derivatives."""
    t = np.clip(t / np.nanpercentile(t, 99), 1e-4, None)   # continuum ≈ 1
    dwn = np.median(np.diff(w))
    shift = np.gradient(t, w)                                # T(ν+δ)-T ∝ dT/dν
    tb = gaussian_filter1d(t, 0.10 / 2.3548 / dwn)           # broaden ILS by ~0.1 cm-1
    ils = tb - t                                             # symmetric broadening
    amt = t * np.log(t)                                      # amount scale ε
    mask = (1.0 - t) > 0.02                                  # line-core channels only
    return {"shift": shift, "ILS": ils, "amount": amt}, mask

def zscore(v):
    v = v - v.mean(); s = v.std(); return v / s if s > 0 else v

print(f"{'window':6} {'EOF1%':>6} {'r(air)':>7}  {'basis |corr|':40} joint R²")
fig, axes = plt.subplots(1, len(labels), figsize=(5 * len(labels), 4), squeeze=False)
for j, lab in enumerate(labels):
    sl = slice(bounds[j], bounds[j + 1])
    e, var1, corr_air = band_eof1(j)          # this band's OWN leading EOF
    w, t = wn[sl], T[sl]
    B, mask = bases(w, t)
    em = zscore(e[mask])                        # restrict to line-core channels
    corrs = {k: abs(np.corrcoef(zscore(v[mask]), em)[0, 1]) for k, v in B.items()}
    M = np.vstack([zscore(v[mask]) for v in B.values()]).T
    coef, *_ = np.linalg.lstsq(M, em, rcond=None)
    R2 = 1 - np.sum((em - M @ coef) ** 2) / np.sum(em ** 2)
    best = max(corrs, key=corrs.get)
    ranked = ", ".join(f"{k}={corrs[k]:.2f}" for k in
                       sorted(corrs, key=corrs.get, reverse=True))
    # per-line test: does |EOF1| track line depth even if its sign doesn't match a global basis?
    depth = np.clip(1.0 - t / np.nanpercentile(t, 99), 0, None)
    c_mag = np.corrcoef(np.abs(e[mask]), depth[mask])[0, 1]
    print(f"{str(lab):6} {var1*100:5.1f}% {corr_air:+6.2f}  {ranked:34} {R2:.2f}  "
          f"corr(|EOF1|,depth)={c_mag:+.2f}")
    s = np.sign(np.corrcoef(zscore(B[best][mask]), em)[0, 1])
    axes[0, j].plot(w, s * zscore(e), lw=0.5, label="band EOF1")
    axes[0, j].plot(w, zscore(B[best]), lw=0.5, alpha=0.7, label=f"{best} basis")
    axes[0, j].set_title(f"{lab}: EOF1 ({var1*100:.0f}%) vs {best.split()[0]}")
    axes[0, j].legend(fontsize=8); axes[0, j].set_xlabel("wavenumber cm$^{-1}$")
fig.suptitle("Per-band EOF1 line-shape attribution (shift vs ILS-width vs amount)")
fig.tight_layout(); fig.savefig(FIGS / "em27_m4_eof1_lineshape.png", dpi=120)
print("saved figures/em27_m4_eof1_lineshape.png")
