"""Plot representative EM27 ILS shapes at their retrieved (optimized) widths.

  physical NB ME/PE  @ 0.474 cm⁻¹  (ensemble-median retrieved FWHM)
  Gaussian           @ 0.439 cm⁻¹  (ensemble-median retrieved FWHM)
  bare ME/PE sinc                  (unapodized reference — shows the side-lobes)

Left: linear response (main lobe).  Right: |response| on a log axis (wings/side-lobes).
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))
from em27gert.instrument import ils_physical, ils_gaussian, ils_from_me_pe, _fwhm_of

ME1, PE1 = 0.9841, 0.003
curves = [
    ("physical NB ME/PE @ 0.474", ils_physical(1.8, ME1, PE1, apod="nb_medium", fwhm_cm=0.474),
     "tab:red", "-"),
    ("Gaussian @ 0.439", ils_gaussian(0.439), "tab:blue", "-"),
    ("bare ME/PE sinc (unapodized)", ils_from_me_pe(1.8, ME1, PE1), "0.5", "--"),
]

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for lab, ils, col, ls in curves:
    off, r = ils.wn_offsets, ils.response / ils.response.max()
    fw = _fwhm_of(off, r)
    ax[0].plot(off, r, ls, color=col, lw=1.5, label=f"{lab}  (FWHM {fw:.3f})")
    ax[1].plot(off, np.abs(r), ls, color=col, lw=1.3, label=lab)

ax[0].set_xlim(-2, 2); ax[0].axhline(0, color="k", lw=0.4); ax[0].axhline(0.5, color="grey", lw=0.4, ls=":")
ax[0].set_xlabel("wavenumber offset (cm$^{-1}$)"); ax[0].set_ylabel("normalized response")
ax[0].set_title("ILS main lobe (linear)"); ax[0].legend(fontsize=8)

ax[1].set_xlim(-4, 4); ax[1].set_yscale("log"); ax[1].set_ylim(1e-4, 1.3)
ax[1].set_xlabel("wavenumber offset (cm$^{-1}$)"); ax[1].set_ylabel("|normalized response|")
ax[1].set_title("wings / side-lobes (log)"); ax[1].legend(fontsize=8)

fig.suptitle("EM27 ILS shapes at retrieved widths — physical NB ME/PE vs Gaussian (+ bare sinc)")
fig.tight_layout()
out = PROJ / "figures/em27_ils_shapes.png"
fig.savefig(out, dpi=130, bbox_inches="tight")
print("wrote", out)
