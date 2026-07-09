# Changelog

Notable changes to the EM27/SUN ↔ GERT validation project.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). This is
a research project rather than a released library, so entries are **dated
milestones** (M0–M5, see [`PLAN.md`](PLAN.md)) rather than semantic versions.

Detailed reasoning, dead ends and evidence for each finding live in
[`docs/investigation_log.md`](docs/investigation_log.md); this file records
*what changed*, not *how it was worked out*.

This project depends on [`gert`](../../gert) as an installed library and ships
only EM27-specific readers, instrument definition, scene builder and retrieval
glue — no radiative transfer or optimal-estimation code.

---

## [Unreleased]

Nothing pending.

---

## [M5 — nuisance ablations & spectroscopy] — 2026-07-06

Ranked the error budget and tested several candidate spectroscopic fixes for the
CO₂/airmass bias.

### Added
* Nuisance ablations (ILS, dispersion, solar Doppler) and an RT null test
  confirming that scattering into the FOV is ~1e-7 of the direct beam, so
  `TransmissionSolver` is the correct physics for direct-sun FTIR.
* Two-layer H₂O in the state vector.
* Spectroscopy experiments and their build tooling: `scripts/build_co2_absco.py`,
  `scripts/build_h2o_continuum.py`, `scripts/build_o2_cia.py`,
  `scripts/run_co2_lineshape_test.py`, `scripts/run_co2_ensemble.py`,
  `scripts/derive_adcf.py`, and the residual-EOF analyses
  (`scripts/eof_lineshape.py`, `scripts/eof_residuals.py`).
* ILS experiments (`scripts/run_ils_experiment.py`, `plot_ils_*.py`) and a
  physical ILS with an analytically-retrieved per-band width.
* Specs and record: `docs/co2_lineshape_experiment.md`,
  `docs/investigation_log.md`.

### Findings
* **The error budget is ILS-dominated** (ILS ≫ dispersion > Doppler), and the
  residual is largely ILS-*shape*-independent once the width is retrieved.
* **qSDV line shape does not help the EM27**, tested across three ILS choices.
* The XCH₄ offset localizes to the ILS — not the prior and not line strength.
* O₂ / `XAIR` is a **separate spectroscopy gap** from the CO₂ bias.
* Data preparation: interpolating the measured spectrum was contaminating the
  residual; replaced by an observation operator applied inside GERT, so the
  measurement is never resampled.

---

## [M4 — full-day closed loop] — 2026-07-01

First full-day comparison against the PROFFAST L2 product.

### Added
* `scripts/run_m4.py` — batch retrieval over a day of scans (fork-parallel,
  shared ABSCO), writing `data/m4_results.csv`.
* `em27gert/retrieval.py` — reusable retrieval driver.
* Diagnostics: `scripts/run_m4_residuals.py`, `scripts/plot_xgas_timeseries.py`,
  `scripts/plot_residual_airmass.py`, `scripts/plot_timeofday.py`,
  `scripts/plot_xgas_covary.py`, `scripts/run_xgas_full.py`.

### Findings
* **GERT reproduces COCCON/PROFFAST XCO₂ to ≈0.2 ppm once airmass-corrected**,
  and XCH₄ to ≈2 %.
* An **airmass-dependent XCO₂ systematic** appears across the day (the classic
  ground-based FTIR signature), alongside a near-constant XCH₄ offset.
* A forward run held at PROFFAST's retrieved columns fits essentially as well as
  the GERT solution, so the two are consistent to within the spectral
  information content.
* A real ~20 ppm wavenumber-grid / dispersion calibration effect was identified
  and corrected.

---

## [M3 — single-spectrum retrieval] — 2026-06-29

Optimal-estimation retrieval of one spectrum, compared to its PROFFAST row.

### Added
* `gert.GERTRetrieval` + `TransmissionSolver` driver: per-gas column scaling,
  `solar_gain` continuum, first-order dispersion, solar Doppler.
* Retrieval of the **effective ILS resolution** (the O₂ window is dropped;
  airmass is taken from the measured ground pressure `gndP`).

### Findings
* The ILS dominates the retrieved column — the finding that set the direction of
  M5.
* Residual FFT showed power at the line-spacing scale, **not** channel fringes,
  so a higher-order or spline baseline would not help.

---

## [M0–M2 — ingest, instrument, forward model] — 2026-06-29

Data ingest, the EM27 instrument in GERT, and open-loop forward-model validation.

### Added
* `em27gert/readers.py` — PROFFAST `.BIN` spectra, TCCON `.map` priors,
  `comb_invparms` L2, ground pressure, `ils_list.csv`.
* `em27gert/scene.py` — `.map` → `gert.AtmosphericProfile` (wet→dry VMRs,
  hPa→Pa, surface→TOA reversal, surface pressure from `gndP`).
* `em27gert/instrument.py` — the four COCCON/PROFFAST microwindows
  (XCO, XCH4, XCO2, O2) plus a self-apodizing ILS built from the measured
  modulation efficiency / phase error (`ils_from_me_pe`).
* `em27gert/proffast_bin.py` — OPUS/PROFFAST binary reader.
* `scripts/run_m0_m2.py` and `notebooks/em27_realdata.ipynb` (M0–M2 driver plus
  the forward-vs-real and line-position diagnostic figures).
* `PLAN.md`, `README.md`, `ONBOARDING.md`, `pyproject.toml` (pins `gert`).

### Findings
* **M0:** the GERT prior column reproduces PROFFAST L2 to **0.07 ppm XCO₂**
  (431.37 vs 431.44), validating the scene conversions.
* **M1:** the ABSCO table needed extension for the EM27 windows; the O₂ window
  uses the 1.27 µm `o2_1p27` table, *not* the 760 nm A-band `o2`.
* **M2:** open-loop residuals after a continuum fit are **XCH₄ 1.5 %,
  XCO₂ 2.7 %, XCO 4.6 %**, and the high-resolution telluric line positions match
  the measured spectrum exactly.

### Gotchas established (see [`ONBOARDING.md`](ONBOARDING.md) §6)
* `ForwardModel.run().y` is returned in **wavelength order**, while
  `SpectralWindow.wn_instrument` ascends in **wavenumber**. Comparing them
  without reversing looks like a catastrophic forward-model failure and is not.
* `.map` VMRs are wet; GERT gases are dry-air mole fractions.
* Airmass must use PROFFAST's `XAIR` definition for a fair comparison.