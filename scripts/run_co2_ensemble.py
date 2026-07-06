"""Step 1 B test: does CO2 line-shape (Voigt/qSDV/qSDV+LM) reduce the XCO2
airmass slope?  Runs the ensemble on the CO2 window for each tier (fork-parallel,
ABSCO overridden in the parent and shared copy-on-write).
"""
import argparse, sys, time, multiprocessing as mp
from pathlib import Path
import numpy as np
from scipy.stats import linregress

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ))
from em27gert.readers import read_invparms, read_ils_list
from em27gert.retrieval import retrieve_spectrum
from em27gert.instrument import EM27_WINDOWS, ils_from_me_pe

CO2_WIN = EM27_WINDOWS[2:3]
_ABSCO = _SOLAR = _INV = _DATA = _ILS = None


def build_proffast_ils(data_dir, nu_ref=6281.0, fov_semi_rad=2.36e-3):
    """PROFFAST-style physical ILS: empirical ME/PE self-apodized 1.8 cm OPD-sinc,
    convolved with the finite-FOV (semi-angle) — the same ILS PROFFAST uses."""
    from gert.instrument import ILS
    il = read_ils_list(data_dir / "ils_list.csv").iloc[0]
    base = ils_from_me_pe(1.8, float(il.ME1), float(il.PE1), float(il.ME2), float(il.PE2))
    off = np.asarray(base.wn_offsets); resp = np.asarray(base.response)
    doff = off[1] - off[0]
    width = nu_ref * (1 - np.cos(fov_semi_rad))          # FOV boxcar full width [cm-1]
    k = max(1, int(round(width / doff)))
    if k > 1:
        resp = np.convolve(resp, np.ones(k) / k, mode="same")
    return ILS(type="tabulated", wn_offsets=off, response=resp)


def _worker(pick):
    try:
        r = retrieve_spectrum(pick, _INV, _DATA, _ABSCO, _SOLAR,
                              windows=CO2_WIN, gases=("co2", "ch4", "h2o"), ils=_ILS)
        return (pick, r["airmass"], r["chi2"], r["converged"],
                r["xco2_gert"], r["xco2_proffast"])
    except Exception as e:
        return (pick, None, None, False, None, str(e))


def override_co2(absco, npz_path):
    z = np.load(npz_path); wn_new, ds_new = z["wn"], z["dataset"]
    t = absco["co2"]
    i0 = int(np.argmin(np.abs(t.wavenumber - wn_new[0]))); i1 = i0 + len(wn_new)
    t.dataset[:, :, i0:i1] = ds_new.astype(t.dataset.dtype)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gert", default="../../gert")
    ap.add_argument("--stride", type=int, default=13)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--proffast-ils", action="store_true",
                    help="use PROFFAST's ME/PE + FOV physical ILS instead of the Gaussian")
    ap.add_argument("--apod", default=None,
                    help="use from_mopd(1.8, APOD) apodized ILS (e.g. NB, HN) — "
                         "physical AND smooth/well-fitting")
    args = ap.parse_args()
    global _ABSCO, _SOLAR, _INV, _DATA, _ILS
    gert = Path(args.gert).resolve()
    from gert.absco import ABSCOTable
    from gert.solar import SolarSpectrum
    _ABSCO = ABSCOTable.load_all(str(gert / "input/absco/absco.h5"))
    _SOLAR = SolarSpectrum.load(str(gert / "input/solar/solar.h5"))
    _DATA = PROJ / "data/GSFC_SN245_260406"
    _INV = read_invparms(_DATA / "comb_invparms_GSFC_SN245_260406-260406.csv")
    if args.apod:
        from gert.instrument import ILS
        _ILS = ILS.from_mopd(1.8, apodization=args.apod)
        print(f"ILS = from_mopd(1.8, {args.apod}) apodized")
    elif args.proffast_ils:
        _ILS = build_proffast_ils(_DATA)
        off = np.asarray(_ILS.wn_offsets); r = np.asarray(_ILS.response) / np.max(_ILS.response)
        idx = np.where(r >= 0.5)[0]
        print(f"ILS = PROFFAST ME/PE+FOV  (FWHM~{off[idx[-1]]-off[idx[0]]:.3f} cm-1, "
              f"min={r.min():+.3f})")
    else:
        _ILS = None
        print("ILS = Gaussian effective (res 0.44)")
    sn = _INV[_INV.index.str.contains("SN")]
    sn = sn[sn["appSZA"].astype(float) < 80].sort_values("UTtimeh")
    picks = list(sn.index[::args.stride])
    print(f"{len(picks)} soundings; CO2 window; per-tier airmass slope\n")
    print(f"{'tier':12}{'n':>4}{'bias':>9}{'RMS':>8}{'airmass slope':>16}{'r':>7}{'@am=1':>9}")

    for tier in ("voigt", "sdvoigt", "sdvoigt_lm"):
        override_co2(_ABSCO, PROJ / f"data/co2_absco_{tier}.npz")
        t0 = time.time()
        rows = []
        ctx = mp.get_context("fork")
        with ctx.Pool(args.workers) as pool:
            for res in pool.imap_unordered(_worker, picks):
                if res[4] is not None and np.isfinite(res[4]):   # retrieval produced XCO2
                    rows.append(res)
        am = np.array([r[1] for r in rows]); d = np.array([r[4] - r[5] for r in rows])
        chi = np.median([r[2] for r in rows])
        lr = linregress(am, d)
        print(f"{tier:12}{len(rows):>4}  chi2={chi:6.2f}{d.mean():+9.2f}{np.sqrt((d**2).mean()):8.2f}"
              f"{lr.slope:+12.2f}{'':>4}{lr.rvalue:+7.2f}{lr.intercept:+9.2f}"
              f"   ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
