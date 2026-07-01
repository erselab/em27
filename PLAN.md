# EM27/SUN ↔ GERT: Comparison & Closed-Loop Validation Plan

**Goal.** Use real EM27/SUN (COCCON) measurements as the first downstream
*example project* that consumes `gert` as an installed library, and close the
full loop: real spectra → GERT forward simulation → GERT retrieval → compare
against the PROFFAST L2 product (COCCON `invparms`).

This is GERT's first validation against real data and the template for how an
external project depends on the library (it ships *its own* readers,
instrument definition, and scene; `gert` provides RT + retrieval).

---

## 0. The dataset we have

`data/GSFC_SN245_260406/` — EM27/SUN **SN245** at **NASA GSFC**, **2026-04-06**
(COCCON processing):

| Item | File(s) | Role |
|:-----|:--------|:-----|
| Calibrated solar spectra | `260406_spectra/cal/*.BIN` (2607 files) | **Measurement.** `SM`/`SN` = the two EM27 detector channels per scan (~1303 scans). PROFFAST/OPUS binary: CRLF text header + float32 block. |
| A-priori atmosphere | `map/go_39N_077W_*Z.map` (8 files, 3-hourly) | TCCON **GINPUT** format: `Height, Temp, Pressure, Density, h2o, hdo, co2, n2o, co, ch4, hf, o2, gravity`. **Wet** VMRs. |
| Ground pressure | `pressure/b33_20260406.csv` | Surface-pressure time series → dry-air column / `XAIR`. |
| Instrument line shape | `ils_list.csv` | Per-channel modulation efficiency `ME` and phase error `PE` (PROFFAST ILS parameterization). |
| **PROFFAST L2 (truth)** | `comb_invparms_GSFC_SN245_260406-260406.csv` | Per-spectrum `XH2O, XAIR, XCO2, XCH4, XCO, XCO2_STR, XCH4_S5P`, plus `gndP, gndT, appSZA, azimuth, lat/lon/alt, niter, rms`. |
| COCCON archive | `COCCON_*.nc`, `groundbased_ftir.*.h5` | Standard public L2 (cross-check / metadata). |

**Why this is the right first target:** EM27/SUN is direct-sun FTIR — no
surface, no scattering geometry to model. GERT's `TransmissionSolver`
(direct-sun, line-by-line transmission) is exactly the correct physics, and the
retrieval is a clean column-scaling problem. It exercises the *spectroscopy +
instrument + solar + nuisance* machinery in isolation, before we ever add
scattering.

---

## 1. Project layout

```
em27/                               # separate git repo; depends on gert>=0.1
├── data/GSFC_SN245_260406/         # the dataset above
├── em27gert/                       # project glue (NOT part of core gert)
│   ├── proffast_bin.py             # .BIN reader (seeded from gert/utils)
│   ├── readers.py                  # .map, invparms, ground-P, ils_list
│   ├── instrument.py               # EM27 windows + ILS from ME/PE
│   └── scene.py                    # .map (wet) → gert.AtmosphericProfile (dry)
├── notebooks/
│   ├── em27_realdata.ipynb         # data exploration (from tccon_real_data)
│   └── em27_forward_vs_real.ipynb  # fwd-vs-real + retrieval (from tccon demo)
├── PLAN.md                         # this file
├── README.md
└── pyproject.toml                  # pins gert
```

`em27gert` is deliberately thin: everything physics/RT lives in `gert`; the
project only knows how to *read EM27 data* and *describe the EM27 instrument*.

---

## 2. Milestones

> **STATUS (all done).** M0–M5 complete; `notebooks/em27_realdata.ipynb` runs the
> whole chain. Headline: GERT reproduces PROFFAST **XCO₂ to ≈0.2 ppm (airmass-
> corrected)**, **XCH₄ to ≈2 %**; error budget is **ILS-dominated**. See
> `ONBOARDING.md` §5 for per-milestone results and §7 for open threads.

### M0 — Ingest & alignment *(readers + one smoke test)* — ✅ done
- Finish `em27gert/readers.py`:
  - `read_proffast_bin(path)` → ν grid, measured spectrum, geometry (SZA,
    azimuth, UT time, lat/lon/alt). (Have the converter; wrap it.)
  - `read_map(path)` → DataFrame of z, T, P, density, wet VMRs.
  - `read_invparms(path)` → L2 DataFrame, indexed by `spectrum` filename.
  - `read_ground_pressure`, `read_ils_list`.
- **Smoke test:** load one scan, plot its spectrum, and join it to its
  `invparms` row by `spectrum` filename (the `.BIN` name is the key — confirmed
  present in the CSV). Pick the nearest `.map` by UT hour.
- *Exit:* one spectrum + its prior + its PROFFAST answer, all in memory.

### M1 — EM27 instrument model in GERT — ✅ done
- **Microwindows** (`gert.SpectralWindow`), COCCON/PROFFAST set:
  - XCO₂ ≈ 6173–6390 cm⁻¹ · XCH₄ ≈ 5897–6145 · O₂ ≈ 7765–8005 (airmass)
    · H₂O · XCO ≈ 4208–4257 cm⁻¹ (extended `SM` channel).
- **ILS:** EM27 OPD ≈ 1.8 cm → nominal ~0.5 cm⁻¹ resolution. First pass: a
  sinc ILS at nominal resolution. Refine: fold the `ME`/`PE` from `ils_list.csv`
  into a self-apodized ILS (modulation-efficiency taper + phase). *This is the
  dominant forward-model uncertainty — treat its fidelity as an explicit knob.*
- **ABSCO coverage check (gating risk):** the current table is EMIT/OCO SWIR.
  Verify it spans the CO (≈4210) and O₂-A (≈7765) bands. If O₂-A is missing,
  derive `XAIR`/dry column from `gndP` (PROFFAST `gndP` is in the L2 file) and
  drop the O₂ window for now.
- **Solar:** `gert.SolarSpectrum` (recently extended to full EMIT range) +
  the solar-Doppler pre-shift we added. Confirm it covers 4200–8100 cm⁻¹.

### M2 — Forward vs. real, *open loop* (the physics check) — ✅ done
- One clean spectrum (low SZA, low `invparms` rms): build the scene from its
  `.map` + ground pressure; geometry from the `.BIN`. Run
  `gert.ForwardModel` with `TransmissionSolver` → simulate; apply solar ×
  Doppler × ILS.
- Overplot GERT vs measured per microwindow; inspect residuals.
- *Exit:* residuals spectrally unbiased and < a few % after a fitted continuum
  — validates spectroscopy + ILS + solar end-to-end **before** retrieving.

### M3 — Single-spectrum retrieval vs PROFFAST — ✅ done
- `gert.GERTRetrieval` + `TransmissionSolver`, state vector:
  per-gas column **scaling** (`transmission_scaling`) + continuum/baseline
  + **dispersion** nuisance + **solar Doppler** (all already in `gert`).
- Convert retrieved scaled column → Xgas (dry-air column from O₂ *or* `gndP`,
  same airmass definition PROFFAST uses for `XAIR`).
- Compare to that spectrum's `invparms` row.
- *Exit:* within PROFFAST single-sounding scatter (XCO₂ ≈ 0.5–1 ppm).

### M4 — Full-day time series, *closed loop* (the headline) — ✅ done
- Loop all ~1303 scans → GERT Xgas vs PROFFAST `invparms` across the day.
- Diagnostics: bias / RMS / correlation of (GERT − PROFFAST) for XCO₂, XCH₄,
  XCO; dependence on **SZA/airmass** (the classic FTIR systematic) and on time.
- *Exit:* "GERT reproduces COCCON L2 to within X ppm with airmass slope Y."

### M5 — RT-fidelity & nuisance experiments (GERT's value-add) — ✅ done
- **Solver ladder in direct-sun geometry:** `TransmissionSolver` vs adding
  scattering (`VectorDOSolver`/`XRTMSolver`). Largely a null test confirming
  transmission suffices, but it *quantifies* the Rayleigh/aerosol multiple-
  scatter contribution for any high-AOD scans — the cleanest possible probe of
  the surface/aerosol-MS theme from the GERT roadmap (here with no surface).
- **Nuisance ablations:** ILS (nominal sinc vs ME/PE), solar-Doppler on/off,
  dispersion on/off — show each one's ppm-level impact on retrieved Xgas.
  Directly exercises the machinery we just built.

---

## 3. Risks / dependencies

- **ABSCO band coverage** for CO (4210) and O₂-A (7765) — gating for those
  windows (M1). Fallback: `gndP`-based dry column, CO₂/CH₄ windows only.
- **ILS fidelity** (ME/PE → GERT ILS) is the main forward-model error (M1/M2);
  budget it explicitly via the M5 ablation.
- **Wet→dry VMR**: `.map` profiles are wet mole fractions; apply the conversion
  documented in the `.map` header when building the scene (M0/`scene.py`).
- **Airmass / Xgas definition** must match PROFFAST's (`XAIR`) for a fair M3/M4
  comparison.
- **Solar line list** resolution/coverage across 4200–8100 cm⁻¹ (M1).

---

## 4. Immediate next actions

1. Add `pyproject.toml` pinning `gert` and `em27gert/__init__.py`.
2. Implement `em27gert/readers.py` (M0) + the one-scan smoke test.
3. Adapt `notebooks/em27_realdata.ipynb` (seeded from the old TCCON notebook)
   to drive M0–M2 on this dataset.