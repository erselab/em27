"""Step 1 A/C test: does CO2 line-shape (Voigt -> qSDV -> qSDV+LM) shrink the
bandhead residual and let the physical sinc ILS fit?  Splices each CO2 table into
the loaded ABSCO in memory (no 2 GB rewrite) and retrieves the CO2 window.

Usage: PYTHONPATH=. python scripts/run_co2_lineshape_test.py --gert ../../gert
"""
import argparse, sys
from pathlib import Path
import numpy as np

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))
from em27gert.readers import read_invparms
from em27gert.retrieval import retrieve_spectrum
from em27gert.instrument import EM27_WINDOWS
from gert.instrument import ILS

CO2_WIN = EM27_WINDOWS[2:3]                     # ('XCO2', 6173, 6390, [co2,ch4,h2o])
BANDHEADS = lambda wn: ((wn >= 6205) & (wn <= 6250)) | ((wn >= 6325) & (wn <= 6365))


def override_co2(absco, npz_path):
    """Splice a CO2 line-shape table over the matching wn range of absco['co2']."""
    z = np.load(npz_path)
    wn_new, ds_new = z["wn"], z["dataset"]
    t = absco["co2"]
    i0 = int(np.argmin(np.abs(t.wavenumber - wn_new[0])))
    i1 = i0 + len(wn_new)
    assert np.allclose(t.wavenumber[i0:i1], wn_new, atol=1e-6), "wn grid mismatch"
    assert t.dataset[:, :, i0:i1].shape == ds_new.shape, "grid mismatch"
    t.dataset[:, :, i0:i1] = ds_new.astype(t.dataset.dtype)


def bandhead_rms(r):
    return dict(r=r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", default="../../gert")
    args = ap.parse_args()
    gert = Path(args.gert).resolve()
    from gert.absco import ABSCOTable
    from gert.solar import SolarSpectrum
    absco = ABSCOTable.load_all(str(gert / "input/absco/absco.h5"))
    solar = SolarSpectrum.load(str(gert / "input/solar/solar.h5"))
    data = PROJ / "data/GSFC_SN245_260406"
    inv = read_invparms(data / "comb_invparms_GSFC_SN245_260406-260406.csv")
    sn = inv[inv.index.str.contains("SN")]; sn = sn[sn.appSZA < 55]
    PICK = sn.job01_rms.astype(float).idxmin()
    row = inv.loc[PICK]
    sinc = ILS.from_mopd(1.8, apodization="none")     # physical OPD-limited ILS

    def run(tier, ils=None):
        override_co2(absco, PROJ / f"data/co2_absco_{tier}.npz")
        r = retrieve_spectrum(PICK, inv, data, absco, solar,
                              windows=CO2_WIN, gases=("co2", "ch4", "h2o"),
                              ils=ils, return_spectra=True)
        wn, res = r["wn"], r["resid"]
        cont = 1.0  # resid already in model radiance units; normalise by band peak
        bh = BANDHEADS(wn)
        return r, (np.sqrt(np.mean(res[bh] ** 2)), np.sqrt(np.mean(res[~bh] ** 2)))

    def run_res(tier, res):
        override_co2(absco, PROJ / f"data/co2_absco_{tier}.npz")
        r = retrieve_spectrum(PICK, inv, data, absco, solar, res_eff=res,
                              windows=CO2_WIN, gases=("co2", "ch4", "h2o"))
        return r["chi2"], r["xco2_gert"] - float(row.XCO2)

    print(f"PICK={PICK}  PROFFAST XCO2={float(row.XCO2):.2f} ppm\n")
    print("=== A: chi2 vs effective resolution, per tier ===")
    print("(if qSDV+LM is the fix: lower chi2 min, shifted toward narrower/physical res)")
    RES = [0.33, 0.36, 0.40, 0.44, 0.50]
    print(f"{'tier':12}" + "".join(f"{r:>13.2f}" for r in RES) + f"{'  best(res)':>14}")
    for tier in ("voigt", "sdvoigt", "sdvoigt_lm"):
        chis, dcs = zip(*(run_res(tier, r) for r in RES))
        best = int(np.argmin(chis))
        cells = "".join(f"  {c:5.2f}/{d:+5.1f}" for c, d in zip(chis, dcs))
        print(f"{tier:12}{cells}   {chis[best]:.2f}@{RES[best]:.2f}")

    print("\n=== C: physical sinc ILS — does correct spectroscopy let it fit? ===")
    print(f"{'tier':12}{'chi2(sinc)':>12}{'XCO2':>9}{'dPROF':>8}")
    for tier in ("voigt", "sdvoigt", "sdvoigt_lm"):
        r, _ = run(tier, ils=sinc)
        print(f"{tier:12}{r['chi2']:12.2f}{r['xco2_gert']:9.2f}{r['xco2_gert']-float(row.XCO2):+8.2f}")


if __name__ == "__main__":
    main()
