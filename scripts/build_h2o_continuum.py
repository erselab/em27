#!/usr/bin/env python3
"""Fold the MT_CKD water-vapor continuum into the H2O line cross-section as an
effective, pressure/temperature-dependent per-H2O cross-section — no gert
forward-model change (same trick as build_o2_cia.py).

Physics (MT_CKD; Mlawer et al. 2012)
------------------------------------
Continuum optical depth over an H2O column N_h2o [molec/cm^2]:

  tau(nu) = RADFN(nu,T) · (T0/T) · N_h2o · [
        C_self(nu)·(T0/T)^texp(nu) · (P_h2o/P0)          # self  (∝ H2O²)
      + C_for(nu)                    · (P_dry/P0) ]        # foreign (∝ H2O·dry)

with RADFN(nu,T)=nu·tanh(c2·nu/2T), c2=hc/k=1.4387769 cm·K; C_self=self_absco_ref,
C_for=for_absco_ref, texp=self_texp from the MT_CKD NetCDF; P0=ref_press,T0=ref_temp.
Factoring out N_h2o gives an effective per-H2O cross-section that adds onto the
H2O line block:  sigma_cont(nu,p,T) = tau/N_h2o.  The **foreign** term (∝ P_dry)
is exact and dominant in these dry near-IR windows; the **self** term (∝ P_h2o)
uses a representative H2O profile (a ~10% correction).

Usage
-----
  PYTHONPATH=. python scripts/build_h2o_continuum.py --gert ../../gert \
      --nc ../../gert/hitran_cache/absco-ref_wv-mt-ckd.nc [--validate] [--merge]
"""
from __future__ import annotations
import argparse, shutil, sys
from pathlib import Path
import numpy as np
import h5py

PROJ = Path(__file__).resolve().parent.parent
C2   = 1.4387769           # hc/k  [cm·K]  (radiation-term constant)
H2O_KEY = "h2o"


def read_mtckd(path):
    """Return (wn, self_ref, for_ref, self_texp, T0, P0) from the MT_CKD NetCDF."""
    try:
        import netCDF4 as nc
        d = nc.Dataset(path); g = lambda k: np.asarray(d.variables[k][:], float)
        return (g("wavenumbers"), g("self_absco_ref"), g("for_absco_ref"),
                g("self_texp"), float(g("ref_temp")), float(g("ref_press")))
    except Exception:
        with h5py.File(path) as d:
            g = lambda k: np.asarray(d[k][:], float)
            return (g("wavenumbers"), g("self_absco_ref"), g("for_absco_ref"),
                    g("self_texp"), float(np.ravel(d["ref_temp"][()])[0]),
                    float(np.ravel(d["ref_press"][()])[0]))


def h2o_profile(gert):
    """Representative dry-air H2O VMR vs pressure [Pa] (for the self term)."""
    sys.path.insert(0, str(PROJ))
    from em27gert.readers import read_invparms, read_spectrum, nearest_map
    from em27gert.scene import map_to_atmosphere
    D = PROJ / "data/GSFC_SN245_260406"
    inv = read_invparms(D / "comb_invparms_GSFC_SN245_260406-260406.csv")
    pick = "260406_163455SN.BIN"; row = inv.loc[pick]
    md = read_spectrum(D / "260406_spectra/cal" / pick)["metadata"]
    atm = map_to_atmosphere(nearest_map(D / "map", md["time_ut_h"], md["date"]),
                            p_surface_pa=float(row["gndP"]) * 100.0)
    return np.asarray(atm.p_levels), np.asarray(atm.gases["h2o"])


def sigma_cont(nc_path, gert):
    """Effective H2O continuum cross-section on the H2O ABSCO (p,T,wn) grid."""
    wn0, cs, cf, texp, T0, P0 = read_mtckd(nc_path)
    with h5py.File(gert / "input/absco/absco.h5", "r") as f:
        wn  = f[f"absco_{H2O_KEY}_wavenumber"][:]
        pg  = f[f"absco_{H2O_KEY}_p"][:]                 # (n_p,) Pa
        Tg  = f[f"absco_{H2O_KEY}_T"][:]                 # (n_p, n_T) K
        sig_line = f[f"absco_{H2O_KEY}_dataset"][:]      # (n_p, n_T, n_wn)
    # MT_CKD onto our wn grid (smooth 10 cm-1 -> 0.01; monotonic interp)
    o = np.argsort(wn0)
    Cs = np.interp(wn, wn0[o], cs[o]); Cf = np.interp(wn, wn0[o], cf[o])
    Te = np.interp(wn, wn0[o], texp[o])
    rad_nu = wn                                          # RADFN = nu·tanh(c2·nu/2T)
    p_prof, vmr_prof = h2o_profile(gert)
    n_p, n_T = Tg.shape
    sig = np.zeros_like(sig_line)
    for ip in range(n_p):
        p_mbar = pg[ip] / 100.0
        vmr = float(np.interp(pg[ip], p_prof[::-1], vmr_prof[::-1]))  # p ascending
        for it in range(n_T):
            T = float(Tg[ip, it])
            RADFN = rad_nu * np.tanh(C2 * wn / (2.0 * T))
            self_t = Cs * (T0 / T) ** Te * (vmr * p_mbar / P0)
            for_t  = Cf * ((1.0 - vmr) * p_mbar / P0)
            sig[ip, it] = (RADFN * (T0 / T) * (self_t + for_t)).astype(np.float32)
    return wn, pg, Tg, sig_line, sig


def validate(gert, sig_cont, picks, orders=(0, 1, 2)):
    sys.path.insert(0, str(PROJ))
    from gert.absco import ABSCOTable
    from gert.solar import SolarSpectrum
    from em27gert.readers import read_invparms
    from em27gert.retrieval import retrieve_spectrum
    D = PROJ / "data/GSFC_SN245_260406"
    absco = ABSCOTable.load_all(str(gert / "input/absco/absco.h5"))
    solar = SolarSpectrum.load(str(gert / "input/solar/solar.h5"))
    inv = read_invparms(D / "comb_invparms_GSFC_SN245_260406-260406.csv")
    h2o = absco[H2O_KEY]; base = h2o.dataset.copy()
    labs = ["XCO", "XCH4", "XCO2"]

    def run_set(bo):
        rows = []
        for p in picks:
            r = retrieve_spectrum(p, inv, D, absco, solar, res_eff=0.477,
                                  retrieve_ils_scale=False, baseline_order=bo,
                                  return_spectra=True)
            yo = r["resid"] + r["yret"]; sig = np.abs(yo) * 0.005 + 1e-12
            off = 0; chi = {}
            for lab, n in r["win_bounds"]:
                s = slice(off, off + n); off += n
                chi[lab] = float(np.mean((r["resid"][s] / sig[s]) ** 2))
            rows.append((r["chi2"], r["h2o_scale"], chi))
        return rows

    def med(rs, f): return np.median([f(r) for r in rs])
    print(f"[validate] {len(picks)} scans — median per-window χ², continuum OFF vs ON, by baseline order")
    print(f"  {'order':>6} {'cont':>5}{'χ²tot':>8}{'h2o_sc':>8}" + "".join(f"{'χ²['+l+']':>10}" for l in labs))
    for bo in orders:
        h2o.dataset = base.copy();                    off = run_set(bo)
        h2o.dataset = (base + sig_cont).astype(base.dtype); on = run_set(bo)
        for tag, rs in [("OFF", off), ("ON", on)]:
            row = "".join(f"{med(rs, lambda r: r[2].get(l, np.nan)):>10.2f}" for l in labs)
            print(f"  {bo:>6} {tag:>5}{med(rs,lambda r:r[0]):>8.2f}{med(rs,lambda r:r[1]):>8.3f}{row}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", default="../../gert")
    ap.add_argument("--nc", default="../../gert/hitran_cache/absco-ref_wv-mt-ckd.nc")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--nval", type=int, default=8, help="scans for --validate")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    gert = Path(args.gert).resolve()

    wn, pg, Tg, sig_line, sig = sigma_cont(Path(args.nc), gert)
    ip = int(np.argmax(pg)); it = Tg.shape[1] // 2
    # diagnostics in the CH4/CO2 windows
    for lo, hi in [(4208, 4257), (5897, 6145), (6173, 6390)]:
        m = (wn >= lo) & (wn <= hi)
        if m.sum() == 0: continue
        li = np.trapezoid(sig_line[ip, it, m], wn[m]); co = np.trapezoid(sig[ip, it, m], wn[m])
        print(f"  {lo}-{hi}: continuum/(line+cont) band-int = {100*co/(li+co):.1f}%  "
              f"(surface p={pg[ip]:.0f}Pa)")
    out = Path(args.out) if args.out else PROJ / "data/h2o_continuum.npz"
    np.savez(out, wn=wn, p=pg, T=Tg, sigma_cont=sig)
    print(f"saved sigma_cont -> {out}")

    if args.validate:
        from em27gert.readers import read_invparms
        inv = read_invparms(PROJ / "data/GSFC_SN245_260406/comb_invparms_GSFC_SN245_260406-260406.csv")
        sn = inv[inv.index.str.contains("SN")]; sn = sn[sn["appSZA"].astype(float) < 80].sort_values("UTtimeh")
        picks = list(sn.index[:: max(1, len(sn) // args.nval)])[:args.nval]
        validate(gert, sig, picks)

    if args.merge:
        h5 = gert / "input/absco/absco.h5"; bk = h5.with_name("absco.bkup_precont.h5")
        if not bk.exists(): shutil.copy2(h5, bk); print(f"backed up -> {bk.name}")
        with h5py.File(h5, "a") as f:
            d = f[f"absco_{H2O_KEY}_dataset"]
            if d.attrs.get("h2o_continuum"): sys.exit("ERROR: continuum already merged.")
            d[...] = (sig_line + sig).astype(np.float32)
            d.attrs["h2o_continuum"] = "MT_CKD_4.3"
        print("MERGED MT_CKD H2O continuum into absco.h5 [h2o]")
    elif not args.validate:
        print("(inspect only — pass --validate or --merge)")


if __name__ == "__main__":
    main()
