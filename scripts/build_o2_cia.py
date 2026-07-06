#!/usr/bin/env python3
"""Build the O2 1.27 µm collision-induced-absorption (CIA) contribution and fold
it into the O2 line cross-section as an effective, pressure-dependent per-O2
cross-section — no gert forward-model change required.

Physics
-------
CIA optical depth is density-squared:

    tau_CIA(nu) = SUM_pairs  B_pair(nu,T) * n(O2) * n(partner) * path

with B the HITRAN binary CIA coefficient [cm^5 molecule^-2].  Factoring out the
O2 column (n_O2 * path) gives an effective per-O2 cross-section [cm^2/molecule]:

    sigma_CIA(nu,p,T) = SUM_pairs  B_pair(nu,T) * VMR_partner * n_air(p,T)

n_air = p / (k_B T)  [molecule/cm^3].  N2/O2 are well mixed so VMR_partner is a
constant to <1%.  sigma_CIA drops straight onto the o2_1p27 (n_p, n_T, n_wn)
grid and adds to the line cross-section: tau = (sigma_line + sigma_CIA) * column.

Inputs
------
HITRAN CIA `.cia` files (Karman et al. 2019) covering the 1.27 µm a¹Δg band
(~7700–8050 cm⁻¹).  Download from https://hitran.org/cia/ .  The collision
partner (and its VMR) is auto-detected from each file's header symbol, e.g.
`O2-O2` (self, VMR 0.2095) and `O2-N2` (VMR 0.7808).

Usage
-----
  # inspect only: build sigma_CIA on the o2_1p27 grid, save npz + diagnostics
  PYTHONPATH=. python scripts/build_o2_cia.py --gert ../../gert \
      --cia O2-O2_2019.cia O2-N2_2019.cia

  # additionally MERGE it into absco.h5 (backs up first)
  PYTHONPATH=. python scripts/build_o2_cia.py --gert ../../gert \
      --cia O2-O2_2019.cia O2-N2_2019.cia --merge
"""
from __future__ import annotations
import argparse, shutil, sys
from pathlib import Path
import numpy as np
import h5py

PROJ = Path(__file__).resolve().parent.parent
K_B = 1.380649e-23            # J/K
# Dry-air volume mixing ratios of the well-mixed collision partners.
VMR = {"N2": 0.780840, "O2": 0.209460, "AR": 0.009340, "CO2": 4.2e-4}
O2_KEY = "o2_1p27"            # absco.h5 dataset for the O2 1.27 µm band


# ---------------------------------------------------------------------------
# HITRAN .cia parsing
# ---------------------------------------------------------------------------

def parse_cia(path: Path):
    """Parse a HITRAN .cia file into a list of (pair, T, wn[], B[]) blocks.

    HITRAN CIA format: each temperature is one block — a header line whose first
    20 chars are the collision-pair symbol followed by
    ``numin numax npoints T maxCIA …``, then ``npoints`` lines of
    ``wavenumber  B``.  A single file can hold many temperatures and bands."""
    lines = Path(path).read_text().splitlines()
    blocks, i = [], 0
    while i < len(lines):
        ln  = lines[i]
        sym = ln[:20].strip()
        rest = ln[20:].split()
        if "-" in sym and len(rest) >= 4:
            try:
                npts = int(rest[2]); T = float(rest[3])
            except ValueError:
                i += 1; continue
            wn, B = [], []
            for d in lines[i + 1: i + 1 + npts]:
                p = d.split()
                if len(p) >= 2:
                    try:
                        wn.append(float(p[0])); B.append(float(p[1]))
                    except ValueError:
                        pass
            if wn:
                blocks.append((sym.upper().replace(" ", ""),
                               T, np.array(wn), np.array(B)))
            i += 1 + npts
        else:
            i += 1
    return blocks


def partner_vmr(pair_sym: str) -> tuple[str, float]:
    """From a pair like 'O2-N2' return (partner, VMR).  O2 must be one member."""
    a, _, b = pair_sym.partition("-")
    partner = b if a == "O2" else a
    if partner not in VMR:
        raise ValueError(f"no VMR for collision partner {partner!r} (pair {pair_sym}); "
                         f"known: {list(VMR)}.  (O2-H2O can't be pre-folded — needs a profile.)")
    return partner, VMR[partner]


class PairCIA:
    """B_pair(nu, T) on a target wavenumber grid, linearly interpolated in T."""
    def __init__(self, blocks, wn_grid):
        # keep only blocks overlapping the target window, resample onto wn_grid
        rows, Ts = [], []
        lo, hi = wn_grid[0], wn_grid[-1]
        for _sym, T, wn, B in blocks:
            if wn.max() < lo or wn.min() > hi:
                continue
            B = np.clip(B, 0.0, None)                      # CIA is non-negative
            rows.append(np.interp(wn_grid, wn, B, left=0.0, right=0.0))
            Ts.append(T)
        if not rows:
            raise ValueError("no CIA temperature blocks overlap the O2 window")
        order = np.argsort(Ts)
        self.T = np.array(Ts)[order]
        self.B = np.array(rows)[order]                     # (n_T_cia, n_wn)

    def at_T(self, Tq: float) -> np.ndarray:
        Tq = float(np.clip(Tq, self.T[0], self.T[-1]))
        j = int(np.clip(np.searchsorted(self.T, Tq) - 1, 0, len(self.T) - 2))
        w = (Tq - self.T[j]) / (self.T[j + 1] - self.T[j])
        return (1.0 - w) * self.B[j] + w * self.B[j + 1]


# ---------------------------------------------------------------------------
# Validation: O2/XAIR retrieval, line-only vs line+CIA (in-memory, no file write)
# ---------------------------------------------------------------------------

def validate(gert: Path, sig_cia: np.ndarray, data_dir: Path, pick: str) -> None:
    """Run the opt-in O2/XAIR retrieval on one clean spectrum before vs after
    adding sigma_CIA to the o2_1p27 table (patched in memory)."""
    sys.path.insert(0, str(PROJ))
    from gert.absco import ABSCOTable
    from gert.solar import SolarSpectrum
    from em27gert.readers import read_invparms
    from em27gert.retrieval import retrieve_spectrum, GASES_WITH_O2
    from em27gert.instrument import EM27_WINDOWS

    absco = ABSCOTable.load_all(str(gert / "input/absco/absco.h5"))
    if absco["o2_1p27"].dataset.shape != sig_cia.shape:
        sys.exit(f"ERROR: o2 grid {absco['o2_1p27'].dataset.shape} != CIA {sig_cia.shape}")
    solar = SolarSpectrum.load(str(gert / "input/solar/solar.h5"))
    inv   = read_invparms(data_dir / "comb_invparms_GSFC_SN245_260406-260406.csv")

    def run():
        return retrieve_spectrum(pick, inv, data_dir, absco, solar,
                                 windows=EM27_WINDOWS, gases=GASES_WITH_O2)

    print(f"\n[validate] O2/XAIR retrieval on {pick} (4 windows incl. O2)…")
    before = run()
    o2 = absco["o2_1p27"]                       # patch line -> line+CIA in memory
    o2.dataset = (o2.dataset + sig_cia).astype(o2.dataset.dtype)
    after = run()

    def o2s(x):                                 # o2_scale ≈ 1/XAIR (o2_prior≈0.2095)
        return float("nan") if not np.isfinite(x) else 1.0 / x
    print(f"\n  {'':14}{'XAIR':>9}{'o2_scale≈':>11}{'χ²':>8}")
    print(f"  {'line only':14}{before['xair_gert']:>9.3f}{o2s(before['xair_gert']):>11.3f}{before['chi2']:>8.2f}")
    print(f"  {'line + CIA':14}{after['xair_gert']:>9.3f}{o2s(after['xair_gert']):>11.3f}{after['chi2']:>8.2f}")
    print(f"  {'PROFFAST':14}{before['xair_proffast']:>9.3f}")
    print(f"\n  XAIR shift: {before['xair_gert']:.3f} -> {after['xair_gert']:.3f} "
          f"(target ~{before['xair_proffast']:.2f})")


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", default="../../gert")
    ap.add_argument("--cia", nargs="+", required=True,
                    help="HITRAN .cia files (e.g. O2-O2 and O2-N2)")
    ap.add_argument("--merge", action="store_true",
                    help="add sigma_CIA into absco.h5 o2_1p27 block (backs up first)")
    ap.add_argument("--validate", action="store_true",
                    help="run the O2/XAIR retrieval line-only vs line+CIA (in memory)")
    ap.add_argument("--pick", default="260406_163455SN.BIN",
                    help="spectrum for --validate (default: the clean low-SZA one)")
    ap.add_argument("--data", default=str(PROJ / "data/GSFC_SN245_260406"),
                    help="EM27 dataset dir for --validate")
    ap.add_argument("--out", default=None, help="output npz (default: data/o2_cia.npz)")
    args = ap.parse_args()

    gert = Path(args.gert).resolve()
    absco_h5 = gert / "input/absco/absco.h5"
    with h5py.File(absco_h5, "r") as f:
        wn      = f[f"absco_{O2_KEY}_wavenumber"][:]
        p_grid  = f[f"absco_{O2_KEY}_p"][:]
        T_grid  = f[f"absco_{O2_KEY}_T"][:]                # (n_p, n_T)
        sig_line = f[f"absco_{O2_KEY}_dataset"][:]         # (n_p, n_T, n_wn)
    n_p, n_T = T_grid.shape
    print(f"o2_1p27 grid: {n_p}×{n_T} (p,T) × {len(wn):,} wn  "
          f"[{wn.min():.1f}-{wn.max():.1f} cm-1]")

    # Parse each .cia file, group by collision pair, build an interpolator each.
    pairs: dict[str, PairCIA] = {}
    vmrs:  dict[str, float]   = {}
    for path in args.cia:
        blocks = parse_cia(Path(path))
        by_pair: dict[str, list] = {}
        for b in blocks:
            by_pair.setdefault(b[0], []).append(b)
        for sym, blks in by_pair.items():
            partner, vmr = partner_vmr(sym)
            try:
                pairs[sym] = PairCIA(blks, wn)
            except ValueError as e:
                print(f"  {sym}: skipped ({e})"); continue
            vmrs[sym] = vmr
            Ts = pairs[sym].T
            print(f"  {sym}: partner {partner} (VMR {vmr:.4f}), "
                  f"{len(Ts)} T blocks {Ts.min():.0f}-{Ts.max():.0f} K in-window")
    if not pairs:
        sys.exit("ERROR: no usable CIA pairs found in the input files.")

    # sigma_CIA(nu,p,T) = n_air(p,T) * SUM_pairs VMR_partner * B_pair(nu,T)
    sig_cia = np.zeros_like(sig_line)
    for ip in range(n_p):
        for it in range(n_T):
            T = float(T_grid[ip, it]); p = float(p_grid[ip])
            n_air = p / (K_B * T) / 1e6                    # molecule/cm^3
            acc = np.zeros(len(wn))
            for sym, pc in pairs.items():
                acc += vmrs[sym] * pc.at_T(T)
            sig_cia[ip, it] = (n_air * acc).astype(np.float32)

    # ---- diagnostics: how much continuum did we add? ----
    ip_sfc = int(np.argmax(p_grid)); it_mid = T_grid.shape[1] // 2
    band = slice(None)
    line_pk = sig_line[ip_sfc, it_mid, band].max()
    cia_pk  = sig_cia[ip_sfc, it_mid, band].max()
    line_int = np.trapezoid(sig_line[ip_sfc, it_mid], wn)
    cia_int  = np.trapezoid(sig_cia[ip_sfc, it_mid],  wn)
    print(f"\nAt surface (p={p_grid[ip_sfc]:.0f} Pa, T={T_grid[ip_sfc,it_mid]:.0f} K):")
    print(f"  peak sigma:  line={line_pk:.3e}  CIA={cia_pk:.3e} cm^2/molec")
    print(f"  band-integrated sigma:  line={line_int:.3e}  CIA={cia_int:.3e}")
    print(f"  CIA / (line+CIA) band-integrated = {cia_int/(line_int+cia_int)*100:.1f}%")

    out = Path(args.out) if args.out else Path("data/o2_cia.npz")
    out.parent.mkdir(exist_ok=True)
    np.savez(out, wn=wn, p=p_grid, T=T_grid, sigma_cia=sig_cia,
             sigma_o2_total=(sig_line + sig_cia).astype(np.float32),
             pairs=np.array(list(pairs)), vmr=np.array([vmrs[s] for s in pairs]))
    print(f"\nsaved sigma_CIA (+ combined) -> {out}")

    if args.validate:                          # before any merge -> clean line-only baseline
        validate(gert, sig_cia, Path(args.data), args.pick)

    if args.merge:
        bkup = absco_h5.with_name("absco.bkup_precia.h5")
        if not bkup.exists():
            shutil.copy2(absco_h5, bkup); print(f"backed up -> {bkup.name}")
        with h5py.File(absco_h5, "a") as f:
            d = f[f"absco_{O2_KEY}_dataset"]
            if d.attrs.get("cia_added"):
                sys.exit("ERROR: o2_1p27 already has CIA merged (cia_added=True). "
                         "Restore from absco.bkup_precia.h5 before re-merging.")
            d[...] = (sig_line + sig_cia).astype(np.float32)
            d.attrs["cia_added"]  = True
            d.attrs["cia_pairs"]  = np.array(list(pairs), dtype=h5py.string_dtype())
        print(f"MERGED CIA into absco.h5 [{O2_KEY}]  (line + CIA)")
    else:
        print("(inspect only — pass --merge to add it into absco.h5)")


if __name__ == "__main__":
    main()
