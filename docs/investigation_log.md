# EM27/SUN ↔ GERT — investigation log

How we got from "GERT forward model matches real EM27 spectra to a few percent"
to "the residual is missing CO₂/CH₄ line-shape spectroscopy (speed-dependence +
line mixing)." Written as a chain of hypotheses → experiments → findings, so the
reasoning (and the dead ends) are recoverable.

**Dataset:** EM27/SUN SN245, NASA GSFC, 2026-04-06, COCCON L1 (PROFFAST 2.4.1).
**Truth:** PROFFAST L2 `invparms`. **Clean test spectrum:** `260406_163455SN`
(SZA 33°, lowest L2 rms). **Ensemble:** stride-13 subsample, ~100 SN soundings,
airmass 1.2–2.9.

---

## Headline

GERT reproduces PROFFAST **XCO₂ to ~0.2 ppm (airmass-corrected)** and **XCH₄ to
~2 %**. The limiting residual is **not** the ILS, dispersion, continuum, data
prep, or scattering — it is **per-line CO₂/CH₄ line-shape error** (plain-Voigt
ABSCO missing speed-dependence + line mixing), concentrated at the CO₂ R-branch
bandheads and CH₄ Q-branch and scaling with airmass. This matches the spectroscopy
PROFFAST v2 includes (HITRAN2020 + qSDV + LM) and we do not.

---

## The chain

### 1. M3 — single-spectrum retrieval; the ILS dominates the column
- **Approach.** OE retrieval (`GERTRetrieval` + `TransmissionSolver`): per-gas
  column scaling + per-window `solar_gain` + first-order dispersion + solar
  Doppler; O₂ dropped, airmass from `gndP`; `p_scale` frozen.
- **Finding.** Retrieved columns are extremely sensitive to the ILS: apodization
  swung XCO₂ by tens of ppm (BH3 +63 ppm, HN +1.5 ppm). We **retrieved the
  effective ILS resolution** via a χ² scan → minimum ≈ 0.5 cm⁻¹.
- **Correction.** This overturns onboarding gotcha #5 ("ILS is a sub-% effect").
  True for *line positions* (M2), false for the *retrieved column*.

### 2. Residual FFT — channel fringes? No.
- **Hypothesis.** The repeating residual structure could be optical channeling
  (etalon fringes) → a high-order/spline baseline would remove it.
- **Experiment.** FFT of the M2 forward-vs-measured residual.
- **Finding.** Power sits at the **line-spacing scale (1–3 cm⁻¹)**, not at a
  distinct fringe frequency; residual amplitude tracks the band envelope. → It's
  **line-shape/ILS**, not fringes. A spline baseline would not help.

### 3. Forward-at-PROFFAST — are GERT and PROFFAST consistent?
- **Experiment.** Freeze the gas scales at PROFFAST's columns (forward run at the
  "truth"), let only nuisances adjust.
- **Finding.** Δχ² ≈ +0.05 vs the free retrieval. → The spectrum does **not**
  distinguish our values from PROFFAST's; they're consistent within the spectral
  information content. The leftover χ² is structural, not a gas error.

### 4. Dispersion & the wavenumber grid — a real ~20 ppm calibration
- **Observation.** Retrieved dispersion shift is ∝ ν (−2.05×10⁻⁵, constant across
  all three windows) = a single global **wavenumber-scale offset ~20 ppm**.
- **Checks.** Dispersion off → χ²=35, XCO₂ −26 ppm (so it matters); order-1 ≡
  order-2 (purely linear). The `.BIN` grid is internally consistent to 1e-10 —
  **not** a parser bug. It's the FTS HeNe wavenumber calibration, which PROFFAST
  also fits. Dispersion absorbs it cleanly → **not** the residual driver.

### 5. M4 — full-day time series; the airmass systematic appears
- **Method.** Looped the ensemble (fork-parallel, ABSCO shared copy-on-write,
  ~4 min/100 soundings).
- **Finding.** XCO₂ shows a **+3.1 ppm/airmass slope** (r=0.74) but only **+0.2 ppm
  at airmass 1**; XCH₄ a near-constant **+46 ppb (~2.4 %)** offset; XCO
  +16 ppb/airmass. The airmass slope is the classic FTIR systematic and the
  headline symptom to explain.

### 6. M5 — nuisance ablations + RT null test
- **Ablations.** Impact on retrieved XCO₂: ILS/resolution ≫ dispersion > solar
  Doppler (~2.7 ppm). **ILS-dominated error budget.**
- **RT null test.** Clear-sky multiple scattering into the FOV ≈ 1e-7 of the
  direct beam → `TransmissionSolver` is the correct physics (confirmed).

### 7. O₂ / XAIR — a separate spectroscopy gap
- **Experiment.** Re-added the O₂ 1.27 µm window, retrieved XAIR.
- **Finding.** XAIR = 0.54 (o2_scale ≈ 1.85) vs PROFFAST ≈ 1.0 — O₂ band
  **under-absorbs**. XAIR is pathologically sensitive to resolution, and **no
  single resolution fits both O₂ and CO₂/CH₄**. Cause: the `o2_1p27` ABSCO omits
  the **O₂ 1.27 µm collision-induced absorption (CIA)**. (Parked; O₂/XAIR made a
  diagnostic, not part of the Xgas retrieval.)

### 8. Residual EOF — one airmass-scaled mode, per molecule
- **Method.** EOF/PCA of the M4 residual matrix, **per band**.
- **Findings.**
  - One mode = **80–87 %** of each band's residual variance; **PC1 ↔ airmass
    r ≈ −0.98**. A single fixed spectral pattern scaling with path.
  - It is **not** a global wavenumber shift, ILS-width, or amount error (all fit
    R² < 0.07).
  - **Molecule-resolved regression** (residual per unit absorption, β): loads on
    the band's dominant molecule — CH₄ β=0.63 in the CH₄ window, CO₂ β=0.70 in
    the CO₂ window — while the **shared H₂O carries little in either** (β 0.09–
    0.14). An ILS error is molecule-agnostic and cannot do this → **per-molecule
    spectroscopy.**
  - **Per-line ranking:** worst lines are the CO₂ **R-branch bandheads**
    (6210–6245, 6330–6362) and the CH₄ **Q-branch** (~6002) — the textbook
    **line-mixing** regions. XCO residual sits on H₂O/CH₄, not CO (under-determined).

### 9. Jacobians & stationarity — sanity, and why CO is weak
- Gas-scaling Jacobians (∂y/∂scale) localize each species' lines; **CO's Jacobian
  is ~150× weaker** than CO₂'s (σ_post = prior) → CO is under-determined.
- Stationarity check: Newton step ≈ 10⁻⁵ σ_post → genuinely at the minimum despite
  large Jacobians (K is the forward sensitivity, not the cost gradient).

### 10. Data preparation — interpolation was contaminating the residual
- **Concern (raised).** We were resampling the *measured* spectrum onto a finer
  instrument grid (linear interp; measured at ~1.8 points/FWHM — under-sampled).
- **Test.** linear→cubic interpolation moved the line-core residual by ~1.3 % RMS,
  ~¼–⅓ of the ~4 % core residual. **Real contamination.** (Molecule-agnostic, so
  it didn't overturn the spectroscopy conclusion, but it inflated the numbers.)

### 11. Fix: observation-operator in GERT (no data resampling)
- **Change.** Added `obs_grid` to `gert.SpectralWindow`: the ILS integrates the
  high-resolution model **onto the measurement wavenumbers**; the measured values
  are used directly. Dispersion still optimizes on the (non-uniform) obs grid and
  recovers the same ~20 ppm calibration.
- **Consequence.** The χ² minimum re-tuned to **0.44 cm⁻¹, χ²=1.37** (cleaner than
  the old 2.43). Gas biases grew (XCO₂ +8.7, XCH₄ +67) — the honest result once
  interpolation stops masking spectroscopy.

### 12. The physical ILS — dead ends and a decisive file check
- **Hypotheses tried & discarded:**
  - *Bare ME/PE ILS* → too sharp, χ²=102. Rejected.
  - *FOV ~30 mrad* (assumed) → over-broadens (FWHM 0.73), χ²=28. Rejected.
  - *"COCCON L1 is apodized"* → **searched the files: no apodization.**
- **What the `.BIN` header actually says:** OPD 1.8 cm, **semi-FOV = 2.36 mrad**
  (not 30 — FOV broadening is a negligible 0.017 cm⁻¹), ILS = unapodized OPD-sinc
  with ME/PE (0.9841, 0.003). PROFFAST 2.4.1.
- **Test at the correct FOV (2.36 mrad):** physical sinc ILS χ² = **23.5** vs the
  Gaussian's **1.37** — ~17× worse, and worse at *matched* FWHM (shape, not width).
- **Conclusion.** The Gaussian effective ILS is **absorbing a spectroscopy/line-
  shape error** — the ILS↔spectroscopy degeneracy. PROFFAST can use the physical
  sinc because its spectroscopy is internally consistent; our plain-Voigt ABSCO is
  not.

### 13. What PROFFAST uses — the missing physics named
- PROFFAST v2 uses **HITRAN2020** with **line mixing + quadratic speed-dependent
  Voigt (qSDV+LM)** for CO₂ (speed-dependence for O₂, CIA for O₂, LM for CH₄ 2ν₃).
  The literature is explicit that qSDV+LM **reduces airmass dependence**.
- Our ABSCO is **plain Voigt** → missing exactly this. And LM is strongest at
  **bandheads/Q-branches** and scales with **airmass** — i.e. precisely the EOF /
  per-line / airmass-slope signature. The pieces close.

---

### 14. CO₂ line-shape experiment — qSDV does NOT help the EM27 (via three ILS choices)
- **Built** controlled CO₂ ABSCO tiers from HITRAN2020: Voigt / qSDV / qSDV+LM
  (speed-dependence + first-order line-mixing `Y`, both from HITRAN directly —
  `build_co2_absco.py`, gate #0 bit-identical). qSDV concentrates at the R-branch
  bandheads (6241/6361), where the residual peaked — the right region.
- **CO₂ airmass slope (Voigt → qSDV), by ILS:**
  - Gaussian 0.44 (χ²≈1.3, *tuned*): +5.4 → +12.0 — **worse**.
  - bare sinc (χ²≈20, physical but *bad fit*): −2.75 → ~0 — "better" (looked like
    a reversal, but unreliable — the model is far from the data).
  - **Norton–Beer (χ²≈1.6, physical AND well-fitting): +3.0 → +9.1 — worse.**
- **Conclusion:** the two *well-fitting* ILSs agree — **HITRAN2020 qSDV+LM worsens
  the CO₂ airmass slope for the EM27.** The bare-sinc "reversal" was an artifact
  of its poor fit.
- **Reconciliation:** qSDV+LM helps *high-resolution* TCCON (line shape resolved);
  at EM27 low resolution (line ~0.1 ≪ ILS ~0.44) the line shape is unresolved — it
  can't be validated or exploited, and the small airmass-dependent leak-through
  goes the wrong way. Hence COCCON/TCCON correct airmass **empirically**.
- **Bonus finding:** the COCCON L1 spectra are **effectively apodized**
  (Norton–Beer-like) — bare sinc χ²≈11–22 vs smooth χ²≈1.3–1.6 at matched width, a
  *shape* effect. So the Gaussian ILS used in M3–M5 was ~the real effective ILS,
  not merely a spectroscopy fudge.

### 15. Prior-atmosphere sensitivity — one fixed daily prior vs per-scan nearest .map (2026-07-04)
- **Motivation.** The M3–M4 retrieval picks the nearest of 8 3-hourly `.map`
  files per scan (`nearest_map`). Test how much the a-priori *profile* choice
  drives the day's Xgas by holding **one** prior fixed for every scan (surface
  pressure still per-scan from `gndP`, so only the profile *shape* is frozen).
- **Method.** Added a `map_path` override to `retrieve_spectrum` and a `--map`
  flag to `run_m4.py` (`median` picks the map nearest the day's median obs time →
  `go_39N_077W_2026040615Z.map`). Ran the stride-13 ensemble (98 common converged
  soundings, airmass 1.18–2.87) both ways; `scripts/compare_fixedprior.py`.
- **Findings.**
  - **Biases barely move** (XCO₂ Δbias +0.03 ppm, XCH₄ +0.05 ppb). The retrieval
    scales the column and the 8 maps are the same site/day, so the a-priori
    *column* was never the story. **The XCH₄ +67 ppb offset is untouched** →
    confirmed **not** an a-priori artifact.
  - **XCO cleans up materially:** RMS 12.5→8.7 ppb, airmass slope +17.3→+6.4.
    The time series (`em27_fixedprior_timeseries.png`) shows the cause directly:
    the per-scan-map XCO has **step discontinuities at the map-switch boundaries**
    (~13.5 and ~16.5 UT = the midpoints between the 12/15/18Z maps), because the
    XCO a-priori differs step-wise between maps. `nearest_map` injects those steps
    into the series; a single fixed prior removes them. XCO₂/XCH₄ show only tiny
    versions (their a-priori columns barely change between maps).
  - **Airmass slopes all shrink slightly** (XCO₂ 3.96→3.46, XCH₄ 0.12→−2.80, XCO
    17.3→6.4) — an **aliasing** effect, not physics: which map `nearest_map` picks
    is correlated with time-of-day → airmass, so per-scan switching folded a small
    time-varying prior signal onto the airmass axis. Fixing the prior removes the
    confound; the residual slopes are the real instrument/spectroscopy systematic.
- **Conclusion.** For a **single-site, single-day** run, **one fixed daily prior is
  cleaner than nearest-map switching** — it eliminates discontinuity artifacts
  (clearest in XCO) at the cost of a small, correct time-varying prior signal the
  column is insensitive to. The map choice does **not** explain the XCH₄ offset or
  the bulk of the airmass slope. (Open: bound sensitivity to *which* fixed prior
  via the 00Z/21Z extremes.)

### 16. The XCH₄ +67 ppb offset — localized to the ILS, not the prior or line strength (2026-07-04)
Chased the constant XCH₄ offset down the prior → strength → forward-model chain.
- **Prior — faithful (excluded).** GERT's prior XCH₄ (1883.5 ppb, clean spectrum
  `260406_163455SN`) reproduces an independent dry-air pressure-weighted integral
  of the same `.map` to **0.7 ppb** (CO₂ prior to 0.08 ppm). `map_to_atmosphere` /
  wet→dry is not the culprit.
- **Localized to the retrieved scale.** GERT scales CH₄ by **1.060** and CO₂ by
  1.020 to fit the same spectrum → it needs ~6 % more CH₄ / ~2 % more CO₂ than
  PROFFAST. Gas-specific (CH₄ ≫ CO₂) ⇒ not common-mode normalization. Offset lives
  in the **scale**, not the prior.
- **Line strength — identical HITRAN2020 (excluded).** Band-integrated
  ∫σ dν ≡ ΣS (line-shape-independent) over the CH₄ window 5897–6145:
  | source | ΣS (cm/molec) |
  |---|---|
  | GERT ABSCO ∫σ dν | 7.9475e-20 |
  | GERT cache ΣS (iso 1) | 7.9456e-20 (build faithful, 0.02 %) |
  | HITRAN2020 iso 1 | **7.9456e-20 (ratio 1.0000)** |
  | HITRAN2020 all-iso (1+2+3) | 8.0154e-20 (**+0.88 %**) |
  GERT's CH₄ intensities **are** HITRAN2020, bit-identical; the only strength gap
  vs PROFFAST is the missing minor isotopologues, **+0.88 %** — far short of the
  ~3.6 % ensemble offset (+67 ppb) / +6 % single-spectrum scale.
- **Conclusion.** With prior faithful, intensities identical, isotopologues ~0.9 %,
  and line shape already excluded (qSDV worse, §14), the residual ~2.7 % is a
  **forward-model / ILS↔column degeneracy** effect: GERT's empirical Gaussian ILS
  maps the CH₄ 2ν₃ Q-branch absorption to a different column than PROFFAST's
  physical sinc+ME/PE. Band-specificity (CH₄ +6 % vs CO₂ +2 %) is consistent — the
  CH₄ Q-branch is more ILS-sensitive. **The offset is ILS-driven, not spectroscopy.**
  Actionable fix: a physical ME/PE⊛FOV ILS (needs gert FOV support), *not* an ABSCO
  change. (Adding CH₄ minor isos would claw back only ~0.9 %.)

### 17. Physical ILS + analytically-retrieved per-band width — residual is ILS-shape-independent (2026-07-04)
Made the retrieval use a **physical ILS by default** and retrieve a **per-band ILS
width scale** jointly with the gases, then compared ILS *shapes* head-to-head.

- **New default ILS.** `ils_physical` builds an NB-apodized (Norton–Beer medium)
  self-apodizing ME/PE kernel from `ils_list.csv` at 0.44 cm⁻¹ — physical *and*
  smooth/well-fitting (no bare-sinc sidelobes). Gaussian is now opt-in.
- **Per-band `ils_scale_{b}`** added to the state vector (prior 1.0 ± 0.05),
  retrieved jointly. **Analytic Jacobians** (user requirement — no finite
  differences): the ILS carries a closed-form g'(δ) (cosine-transform derivative),
  giving ∂R/∂s = −(1/(s²D))Σ δ g'(δ/s)(I−R). The **dispersion** Jacobian was also
  converted from gert's finite-difference `convolve_center_derivative` to the
  closed form ∂R/∂ν_c = −(1/D)Σ g'(δ)(I−R). Both verified vs FD (dispersion corr
  0.99997; width median fd/an 0.99997 — the FD "misses" were `np.interp` kink
  artifacts that vanish at small step). *Bonus:* the forward's piecewise-linear
  convolve has sub-grid kinks, so the analytic Jacobian is strictly better-behaved
  than FD here. (gert: `ILS.convolve_derivatives`, `SpectralWindow.convolve_width_derivative`,
  `width_scale` through convolve/`ForwardModel.run(ils_scale=)`, `StateVector.
  transmission_scaling(include_ils_scale=)`; em27: `retrieve_spectrum` defaults.)
- **Experiment** (`run_ils_experiment.py`, stride-13, ~96 soundings): physical vs
  Gaussian, **both width-optimized**.
  | gas | physical bias/slope | Gaussian bias/slope | med FWHM |
  |---|---|---|---|
  | XCO₂ | −7.8 / +4.0 | +9.2 / +5.2 | 0.477 vs 0.438 |
  | XCH₄ | +42.8 / −7.3 | +66.4 / −1.0 | 0.477 vs 0.438 |
  | XCO  | −24.0 / +9.5 | −9.6 / +16.1 | 0.477 vs 0.438 |
- **Findings.**
  1. **The residual is ILS-shape-independent.** Correct-scale mean residuals
     (0.26 % XCH₄, 0.40 % XCO₂ RMS) nearly overlap for the two shapes;
     corr(physical, Gaussian mean residual) = 0.68. The **power spectrum peaks at
     the CO₂ line spacing (~1.9 cm⁻¹)** — the oscillatory pattern is **per-line
     spectroscopy**, not a fringe (fixed period) or an ILS artifact (both shapes
     give it).
  2. **EOF is the same mode for both:** EOF1 = 74 % (physical) / 76 % (Gaussian),
     **corr(airmass) = −0.96 / −0.97**, concentrated at the XCO₂ R-branch bandheads
     (6210–6245, 6330–6362) — the §8 airmass-scaled line-mixing signature, now
     shown to **survive an honest physical ILS**.
  3. **The two optimized ILSs collapse to a width difference.** Their kernels are
     nearly identical smooth ~0.44–0.48 cores (physical carries small real
     sidelobes); the retrieved widths differ by ~0.035 cm⁻¹ (0.477 vs 0.438), and
     *that* — via the ILS↔column degeneracy — flips XCO₂ (−7.8 ↔ +9.2) and shifts
     XCH₄ (+42.8 ↔ +66.4).
  4. **The physical ILS cuts the XCH₄ offset +67 → +43 ppb** (width free),
     consistent with §16's "offset is ILS-driven"; it trades ~24 ppb of constant
     offset for a mild XCH₄ airmass slope (−7.3 vs −1.0).
  5. **Shape → column level via line depth (the systematic between the two configs).**
     The Gaussian retrieves *more* Xgas than the physical ILS at **every** time step
     (a clean constant offset in the time series). Mechanism, verified by forwarding
     both ILSs at the **same** prior column in the CO₂ window: the physical NB kernel
     makes **~1.3 % deeper** lines (mean 1−T = 0.0851 vs 0.0840; deeper cores too) →
     needs a **smaller** column to match the measured depths → less Xgas; the Gaussian
     is shallower → more Xgas. **Counterintuitively the physical ILS is *wider* in
     FWHM (0.477 vs 0.438) yet makes deeper lines** — line depth is set by kernel
     *shape near the peak* (the NB sinc-cusp + small sidelobes preserve core depth
     better than a smooth Gaussian bell), **not by FWHM**. The ~1.3 % depth gap maps
     onto the ~1.3 % XCH₄ / larger (saturation-amplified) XCO₂ offset. Convolution
     conserves each line's *area*, so this is a core↔wing redistribution; the fit
     weights the partly-saturated cores and saturation makes depth↔column sub-linear,
     turning a ~1 % depth difference into a several-% column difference.
- **Takeaway (refined).** Two distinct roles for the ILS, and they must not be
  conflated: the ILS *shape* is nearly immaterial to the **residual pattern** (both
  shapes fit equally well, identical EOF1 — the residual is per-line CO₂/CH₄
  bandhead spectroscopy), **but the shape is decisive for the absolute column
  level** — it sets the line depth, the depth sets the retrieved column (ILS↔column
  degeneracy). So getting the ILS *physically right* matters for absolute Xgas /
  the XCH₄ offset even though it does not reduce the residual. Artifacts:
  `figures/em27_ils_{shapes,timeseries,residuals}.png`,
  `em27_ils_{physical,gaussian}_{eofs,pcs}.png`.

### 18. Isotopologues + O₂ CIA — the absolute-level fixes (2026-07-05)
Two spectroscopy *strength* gaps closed, driving the offsets down and separating
them cleanly from the line-shape residual.

- **Generalized the ABSCO tooling.** `gert/utils/fetch_hitran.py` +
  `build_absco.py` now take `iso_ids: [..]|all` and `line_shape: voigt|sdvoigt|
  sdvoigt_lm` per spec block (multi-iso fetch via `fetch_by_ids`; HAPI sums isos
  with correct per-iso mass + Q; SD/LM worker with p→0 Voigt fallback). Variant-
  aware cache paths keep new blocks from colliding with the legacy iso-1 ones.
- **All-isotopologue rebuild** (`utils/absco_spec_iso.yml`, the 6 EM27-window
  blocks). Verified against the backup: grid **bit-identical** (max|Δ|=0, uniform
  0.01 cm⁻¹, no new gaps/resolution change — the "former merge discontinuity" is
  avoided by using identical wn_min/max/step), cross-sections up by the measured
  minor-iso band fractions: **CO₂ +0.96 %, CH₄ +0.88 %, CO +0.46 %, H₂O +10.2 %
  (HDO-dominated)**.
- **O₂ 1.27 µm CIA** (`scripts/build_o2_cia.py`). Key trick: CIA (density-squared)
  folds into an effective p,T-dependent per-O₂ cross-section,
  σ_CIA(ν,p,T)=Σ_pairs B_pair(ν,T)·VMR_partner·n_air, so it drops onto the
  o2_1p27 grid with **no forward-model change**. From HITRAN2024 O₂-O₂ + O₂-N₂
  `.cia`, CIA is **55 % of the band-integrated O₂ absorption** — exactly the
  missing half that made O₂ under-absorb. `--validate` (in-memory, one clean
  spectrum): **XAIR 0.56 → 1.05** (PROFFAST 1.00), o2_scale 1.79→0.97, O₂-window
  χ² 2.85→1.36. Merged into absco.h5 (`cia_added` guard; backup
  `absco.bkup_precia.h5`). Residual +5 % XAIR = O₂ line speed-dependence
  (PROFFAST has it, we don't) — small, orthogonal.
- **Ensemble effect** (stride-13, 97 scans, all-iso + physical ILS *fixed* at
  0.477, width-opt OFF; `run_ils_experiment.py --configs physical_fixed`):
  | gas | bias | slope/am | prior (physical, pre-iso §17) |
  |---|---|---|---|
  | XCH₄ | **+9.0 ppb** | −4.0 | +42.8 |
  | XCO₂ | **−4.1 ppm** | +3.6 | −7.8 |
  χ²med = 0.93. **XCH₄ offset: +67 (Gaussian) → +45 (physical ILS) → +9 (all-iso).**
- **The residual is now cleanly line-shape.** The airmass-scaled per-line mode is
  *unchanged* by the isotopologue rebuild, as expected: isotopologues fix **band
  strength** (offsets), not **line shape** (the airmass residual). Residual RMS is
  **U-shaped in airmass** (min ≈1.4), and PC2 splits off a **late-day (UT 18–19)
  cluster** — a **new time-of-day systematic** at matched airmass (candidates:
  afternoon H₂O, real T evolution vs the 3-hourly `.map`, instrument/solar drift).
  Figures: `em27_alliso_residual_airmass.png`, `em27_alliso_xgas_timeseries.png`.
- **Per-band EOFs (independent SVD per window — `eof_residuals.py` rewritten; a
  joint all-band SVD blended them).** The bands carry *different* variance
  structures and different physics:
  | band | EOF1 | EOF2 | EOF1 character / EOF2 character |
  |---|---|---|---|
  | XCO₂ | **83.7 %** (r=−0.85) | 6.8 % | line-structured at R-branch bandheads 6210–6250/6330–6362 / small |
  | XCH₄ | 68.6 % (r=−0.91) | 11.7 % | line (Q-branch ~6000) / **broadband tilt** |
  | XCO  | 68.2 % (r=+0.96) | **27.4 %** | mixed / **large broadband tilt** |
  Two physically-distinct residuals now separated: a **line-structured airmass
  mode** (per-line spectroscopy / line-mixing — dominant & cleanest in **XCO₂,
  84 %**, at the bandheads) and a **smooth broadband-tilt mode** (continuum /
  temperature / time-of-day — 12 % in XCH₄, 27 % in XCO, ~absent in XCO₂). So the
  spectroscopy problem lives mostly in XCO₂; the continuum/time-of-day problem in
  XCH₄/XCO.
- **Per-window χ² change from the isotopologues** (median, iso-1→all-iso, physical
  fixed 0.477): **XCO 7.18→2.25 (0.31×)**, XCH₄ 0.39→0.35 (0.89×), **XCO₂
  1.31→1.30 (unchanged)**. The huge **XCO** win is the previously-unmodeled **HDO**
  (+8–10 % of H₂O) in the 2.35 µm window; **XCO₂ χ² is unchanged** because its
  residual is the line-shape mode (EOF1), which isotopologues don't touch — the
  strength change just rescales the retrieved column (offset −8→−4 ppm).
- **Xgas time series vs PROFFAST** (`em27_alliso_xgas_timeseries.png`): XCH₄ tracks
  the diurnal rise (1919→1940 ppb) in shape at a near-constant +9 ppb; XCO₂ drifts
  −2→−7 ppm toward solar noon (the airmass slope + late-day effect); XCO ~+28 ppb,
  noisy (under-determined).

### 19. Remaining *shape* residuals are the resolution floor — three hypotheses excluded (2026-07-05)
After the strength fixes (§18), attacked the two residual EOF modes directly.
All three candidate shape-physics were tested and excluded; the residuals are
small and at the EM27 ~0.5 cm⁻¹ resolution floor.

- **CO₂ bandhead (XCO₂ EOF1, 84 %) is NOT line mixing.** Rigorous test: ran the
  retrieval at high airmass with Voigt vs qSDV+LM CO₂ (in-memory tier override,
  `data/co2_absco_*.npz`) and projected the *change in the fitted residual* onto
  the measured EOF1. corr = **−0.10** (qSDV+LM), −0.11 (qSDV), +0.06 (pure LM);
  the physics changes the residual ~20 % RMS but ~orthogonal to EOF1. EOF1 is also
  orthogonal to the amount direction (corr 0.005) and 60 % of its power is in the
  two bandhead clusters. So first-order LM is the wrong shape (as §14), and since
  ECS refines the *same* coupling on the *same* lines, ECS is **not indicated** —
  consistent with the ECS-sizing note (marginal at EM27 res). **ECS engine not
  worth building.**
- **HITRAN line list is already current (per-line-error path exhausted).** A fresh
  hitran.org fetch of the CO₂ bandhead window is **bit-identical** to our cache
  (1026 lines, ΣS ratio 1.000000, max|Δν|=0) — today's iso rebuild already pulled
  the latest HITRAN. So the bandhead residual is not stale line data; any residual
  per-line error is *in* HITRAN2024 and not fixable by a refresh.
- **H₂O tilt (XCH₄ EOF2 12 %, XCO EOF2 27 %) is NOT the water continuum.** First
  identified the tilt as H₂O-driven, not temperature: **PC2↔H₂O_scale = +0.83
  (XCH₄), −0.98 (XCO)**; PC2↔T_offset ≈ 0; `gndT` flat at 285 K all day
  (`plot_timeofday.py`, `em27_timeofday.png`). Built MT_CKD 4.3 continuum as an
  effective σ (`build_h2o_continuum.py`, foreign folds exact + self via a
  reference H₂O profile; continuum = 10.2/6.0/3.9 % of the H₂O band in
  XCO/XCH₄/XCO₂ — ordering matches the tilt). **But it does nothing at any
  baseline order** (Δχ² ≈ 0 at deg-2 *and* deg-0/1; slightly worse at deg-0). The
  big broadband slope is the **solar/instrument continuum** (deg-0→deg-1 χ² 48→0.96,
  a linear term), not water. So the tilt is **H₂O line/profile structure** the bulk
  `h2o_scale` can't fit (small: XCH₄ χ² already 0.35). Continuum **not merged**
  (benign no-op in these dry windows).
- **Bottom line.** The remaining residuals are *shape*, small, and unidentifiable
  by any spectroscopy model these spectra can constrain — the resolution floor
  (§14 theme, now quantitative). The session's real spectroscopy wins were the
  *strength* fixes (§18: isotopologues, O₂ CIA), which are merged. Tooling added:
  `retrieve_spectrum(baseline_order=)`, `_build_y_obs` polynomial-order baseline.

### 20. Empirical corrections (ADCF / AICF / ACOS EOF) — blocked by single-day data, not code (2026-07-06)
With the HITRAN-accessible physics exhausted (§18–19), turned to the empirical
corrections PROFFAST/COCCON apply. All are limited by having **one clear day at
one site**, not by the methods.

- **Current GERT−PROFFAST state** (all-iso+CIA+physical ILS, 97 scans): XCO₂ bias
  −4.1 ppm / slope +3.6 per airmass; XCH₄ +9 ppb; XAIR 1.05; χ²med 0.93. The
  differences decompose into a **constant per-gas offset** (XCO₂ **−0.9 %**, XCH₄
  **+0.5 %** — AICF territory) + a **loose airmass slope** (ADCF) + time-of-day
  (H₂O) scatter. **GERT already agrees to <1 % absolute with zero empirical
  corrections** — for an independent RT/retrieval vs a calibrated product, that
  *is* the validation result.
- **O₂ speed-dependence — dead end.** HITRAN2024 has the SD *columns*
  (`gamma_sdv_*`, `y_sdv_air`) but **all `nan` for O₂ iso-1** (0/497 finite);
  `absorptionCoefficient_SDVoigt`≡Voigt (ratio 1.000). So the XAIR +5 % can't be
  closed via HITRAN SD (needs an external O₂ study, like ECS/CO₂). Not attributable
  to SD anymore; small.
- **ADCF (airmass correction) — not self-derivable single-day.** Fit
  X_gert = smooth-time-poly + β·(airmass): β_indep = **−3.29** vs the true artifact
  (slope of GERT−PROFFAST) **+3.58** — *opposite sign*; "correcting" with it makes
  the slope worse. The airmass artifact and the true diurnal cycle are degenerate
  on one day (TCCON pools many clear days to break this). Even removing the
  PROFFAST-referenced slope only drops RMS 4.56→4.28 (XCO₂) — the difference is
  dominated by the constant AICF offset, not the slope. `derive_adcf.py`,
  `em27_adcf_{airmass,timeseries}.png`.
- **AICF (absolute scale) — needs external in-situ.** COCCON's AICF is
  PROFFAST-specific (non-transferable); tying to PROFFAST is circular. A real GERT
  AICF requires a coincident aircraft/AirCore column over GSFC ≈ 2026-04-06
  (unknown if it exists) or a TCCON tie. Not derivable from this dataset.
- **ACOS EOF radiance correction — capability built, but single-day EOF confounds
  artifact & signal.** Implemented the OCO-2/ACOS empirical orthogonal-function
  correction: `F'(x)=F(x)·(1+E·c)` with retrieved per-band coefficients and
  **analytic** Jacobian `∂F'/∂c_k=F·E_k` (gert: `GERTRetrieval(eof_basis=)`,
  `StateVector` `eof_{k}` + `eof_coeffs()` + `include_eof`; em27:
  `build_eof_basis`, `retrieve_spectrum(eof_basis=)`). Runs, converges, fits
  coefficients. **But a single-day EOF basis worsens XCO₂**: it absorbs ~6 ppm at
  low airmass (XCO₂ slope +2.2 → +5.8), and **no prior tightness fixes it**
  (even eof_uncert 0.01 already degenerate) — the XCO₂ bandhead EOF is degenerate
  with the CO₂ signal on one day. ACOS avoids this by deriving EOFs from a large,
  diverse training set where artifact decorrelates from geophysical signal.
- **Bottom line.** Both self-derived empirical corrections (ADCF, EOF) hit the
  *same* wall — one clear day at one site can't separate the airmass/diurnal
  artifact from real signal. The code (EOF correction) is ready for a multi-day
  archive; on this dataset the honest result is the **<1 % raw physics agreement**,
  with the residual offset being an AICF calibration constant (external in-situ),
  not something derivable here.

### 21. Xgas co-variation + XCO offset — the joint fit and the water driver (2026-07-06)
Traced the residual XCO offset and mapped how the retrieved Xgases move together
through the day.
- **Joint-fit architecture (clarified).** Every window is fit *simultaneously* in
  one OE solve with **global** per-gas scale factors — a gas appearing in >1 window
  (CO₂, CH₄, H₂O) is constrained by all of them at once; CO/N₂O are single-window.
  The gas Jacobians are shared across windows, so an interferer in one window pulls
  a gas that is anchored elsewhere.
- **XCO offset (+25–30 ppb vs PROFFAST) = CH₄-2.3 µm interferer leak.** In the CO
  window the **CH₄ Jacobian is ~44× the CO Jacobian**; CO is a weak absorber sitting
  under strong CH₄ lines, so small CH₄/ILS mismatch in that window re-levels CO. The
  H₂O changes made it *worse* (+28→+30 with the shape test) because they shift the
  CH₄/H₂O balance in the CO window. It's an interferer-coupling floor, not a CO
  spectroscopy error — consistent with the weak-CO Jacobian noted in §9.
- **Xgas co-variation (`run_xgas_full.py` 4-window incl O₂, `plot_xgas_covary.py`).**
  Through the day: **XCO₂↔XH₂O = −0.83**, **XAIR↔XH₂O = +0.86** — water is the
  dominant co-variation driver. XCO₂ drifts −0.9 ppm/hr below PROFFAST as water (and
  airmass) rise. This motivated testing whether a water **profile-shape** DOF (not
  just the bulk `h2o_scale`) absorbs the drift — see §22.
  Artifacts: `data/xgas_full.csv`, `figures/em27_xgas_{timeseries_all,covary}.png`.

### 22. O₂-pinned H₂O, DFS, and a 2-parameter H₂O shape — water profile is NOT the XCO₂ drift (2026-07-06)
Tested the leading XCO₂-drift hypothesis (§21: XCO₂↔XH₂O −0.83) directly, by giving
H₂O its full data-supported profile freedom and asking if the drift moves.
- **O₂ window pinned to surface pressure.** Added the 1.27 µm O₂ window to the fit
  with the O₂ scale **frozen at 1.0** (column fixed by `gndP`), so the band's H₂O
  lines add water information without O₂/H₂O trading off. Clean: O₂ residual
  unchanged (0.52 %); modestly tightens H₂O/CH₄. `freeze_gas={"o2_1p27":1.0}`.
- **DFS check — does the data support a 2nd H₂O DOF?** Built the 4-window analytic
  Jacobian (O₂ pinned), split H₂O into lower/upper layers by finite difference at a
  range of pivots, and computed `DFS = A[lower]+A[upper]` (sanity: the two layer
  Jacobians sum to the full `h2o_scale` Jacobian, corr 1.000). Result: total H₂O
  DFS peaks at **~1.76 with the split at ~850 hPa** (lower 0.96, upper 0.80), falling
  to 1.12 by 500 hPa (all the water is low). **But the lower/upper Jacobians are 99 %
  collinear** — the 2nd DOF is real but spectrally *weak*; in a raw {lower,upper}
  basis the two fight (strong anti-correlation). ⇒ implement as **scale + shape**
  (sum/difference rotation), moderate prior on the shape.
- **2-parameter H₂O implemented (scale + column-neutral smooth ramp).** New opt-in
  `h2o_shape` element: `vmr(p) = vmr_prior(p)·(1 + s + β·φ(p))` with φ a **smooth
  ramp linear in pressure**, pivot = water-column-weighted mean p, normalized to unit
  water-weighted RMS ⇒ **column-neutral** (β=0.05 → ΔXH₂O 6×10⁻⁶) and nearly
  orthogonal to `h2o_scale`. **Analytic Jacobian** (no FD, per standard): folds the
  per-level mid-point chain rule into one layer sum, with the exact transmission
  fallback `K_mol_h·Σ_l τ_lay·w_lay` when the solver gives only column-level K —
  validated **corr(analytic,FD)=0.99998**, ratio 1.0004. gert:
  `StateVector.transmission_scaling(include_h2o_shape=, h2o_shape_profile=)`,
  `apply()` + `_jacobian_mixed` `h2o_shape` block, `_h2o_shape` on the SV; em27:
  `build_h2o_shape_profile(atm)`, `retrieve_spectrum(retrieve_h2o_shape=)`.
- **Result — the shape is real but does NOT fix the drift (hypothesis falsified).**
  10 soundings, O₂ pinned, scale-only vs scale+shape: β retrieved **consistently
  positive** (+0.11…+0.36, mean +0.22 — water piled toward the surface vs the MAP
  prior; some soundings push past the 0.15 prior, so it's data-driven, matching the
  0.8-DOF DFS). **Yet:** XCO₂ drift **−0.92 → −0.96 ppm/hr (unchanged)**, XCO₂↔XH₂O
  −0.83 → −0.78, **χ² 1.195 → 1.199 (no fit improvement)**, XCO +2 ppb worse. The
  χ²-flat, drift-unchanged outcome is exactly what the 99 % collinearity predicted:
  the shape mode is nearly spectrally invisible, so it re-levels a **near-constant**
  prior-shape offset (β is roughly flat through the day, *not* the monotonic
  afternoon boundary-layer growth) rather than the diurnal signal. **⇒ The XCO₂
  diurnal drift is not a water-profile-shape effect** — it tracks **airmass**, back
  to the ADCF (single-day degenerate, §20). Code left in, default **off**.

### 23. `absco.h5` corruption + rebuild — the merged spectroscopy re-validated (2026-07-17)
A later session found `../../gert/input/absco/absco.h5` **reverted** — a Google
Drive restore had rolled it back to a pre-§18 state: `ch4` narrow (5995–6285,
would not serve XCH₄ @ 5897 or XCO₂ @ 6390), iso-1 values, and **no O₂ CIA**. The
§18 all-isotopologue + CIA work was intact only in a backup (`absco.bkup_precia.h5`
held the all-iso base; the CIA lived in `data/o2_cia.npz`).
- **Rebuilt** `absco.h5 = pre-CIA all-iso base + o2_cia.npz`. Integrity verified:
  the npz reproduces the base line to 2.7e-8 (same grid), CIA = 53.4 % of the O₂
  band, floor/peak 1.7e-5 → 6.3e-3, `cia_added=True`. All four EM27 windows load
  on-grid; O₂ A-band untouched.
- **Re-validated end-to-end** (not in-memory) on the clean spectrum
  `260406_163455SN`: **XAIR 1.050** (PROFFAST 1.000), χ² 1.08, XCH₄ +9 ppb —
  reproducing §18/§20. The corrupted file gave ~0.54.
- **Backups consolidated:** kept only `absco.hitran2020_iso1_original.h5` (the
  original HITRAN2020 iso-1 table). **Open provenance question:** the all-iso `ch4`
  band integral is **0.998×** iso-1 (should be ≥1) — likely a HITRAN2020→2024
  line-list revision; diff the two before trusting `ch4` to sub-% absolute. Both
  `.h5` files live only on the Drive mount (no off-Drive copy).

## Synthesis (current understanding)

The airmass-scaled, molecule-specific, bandhead-localized residual — robust across
FFT, EOF, molecule-resolved regression, and per-line ranking, and surviving the
observation-operator fix — is a **per-line CO₂/CH₄ spectroscopy** signature (the
ILS, dispersion, continuum, interpolation, and scattering were each tested and
excluded as its *cause*; the ILS still sets the absolute column level).

**But the specific fix is NOT HITRAN2020 speed-dependence + first-order line
mixing** — the controlled experiment showed qSDV+LM *worsens* the EM27 airmass
slope with any well-fitting ILS. The reason is fundamental: at EM27 resolution the
line shape is **unresolved** (line ~0.1 cm⁻¹ ≪ apodized ILS ~0.44), so the
measurement cannot validate a line-shape model, and imposing one that's imperfect
for these data hurts. Remaining spectroscopic candidates (full ECS line mixing,
per-line intensity/position errors) are **likely not decidable from EM27 spectra
alone** for the same reason. This is precisely why COCCON/TCCON handle the airmass
artifact with an **empirical** AICF + airmass correction rather than spectroscopy
at these resolutions — the correct operational answer for a low-resolution
instrument.

### ECS line mixing — order-of-magnitude estimate (before building it)
Sizing the remaining hypothesis (full ECS relaxation-matrix LM) for CO₂ 1.6 µm at
the pressures the EM27 sees, to set expectations before sourcing/coding it.

- **Nature.** Line mixing is **intensity redistribution**, not broadening. In the
  relaxation-matrix picture the diagonal carries the pressure-broadened
  widths/shifts; LM is the **off-diagonal** part and moves absorption *between*
  lines. It **conserves integrated intensity** (sum rule) → it changes spectral
  *shape*, not total band absorption or a net Xgas bias. (So it is a candidate for
  the airmass-dependent *shape* residual, but structurally **cannot** explain the
  XCH₄ +67 ppb offset.)
- **On line widths (HWHM):** ~**zero to first order** — individual widths stay set
  by γ_air (~0.07 cm⁻¹/atm). Full ECS gives only a slight *apparent* narrowing of
  the band envelope via wing suppression. LM does not show up as a width change.
- **On line depths (~1 atm, ∝ pressure → weighted to the near-surface path):**
  | region | effect | magnitude (surface p) |
  |---|---|---|
  | line cores | slight enhancement | ≲ 0.1–0.5 % |
  | R-branch bandheads (6241, 6361) | intensity piled up → deeper | ~few tenths % to ~1 % |
  | troughs between strong lines / far wings | sub-Lorentzian: suppressed | several %, up to ~10–30 % in deep troughs |
- **Anchored to our data:** first-order Rosenkranz LM here was ~**0.17 % at 0.5 atm**,
  bandhead-localized (§ Step 1, ≈0.35 % at surface); full ECS at bandheads is the
  same order to a few× larger (~0.5–1 %). qSDV alone was 0.2 % RMS in the bandheads.
- **Detectability verdict.** Against PROFFAST rms ~0.35 % and an effective ILS ~0.44
  cm⁻¹ **wider than the line spacing** (so the troughs where LM is largest are
  unresolved/averaged), the exploitable ECS signature (~sub-1 % bandhead, few-%
  unresolved troughs) sits **at or below the EM27 noise/resolution floor**. Worth
  trying as the physically-correct version (first-order `Y` is wrong-signed at
  bandheads), but the magnitudes predict a **marginal** effect — consistent with the
  "unresolved line shape → not decidable from EM27 alone" conclusion above.

## What was ruled out (and how)
| Candidate | Ruled out by |
|---|---|
| Channel fringes | FFT: power at line-spacing, not a fringe frequency |
| Continuum/baseline order | FFT + smooth EOF2 only 8 % |
| Wavenumber-grid / parser bug | grid consistent 1e-10; dispersion absorbs it |
| Global ILS width / shift / amount | EOF attribution R² < 0.07 |
| ILS as *cause of residual* | molecule-resolved test (H₂O loads low); physical-ILS test |
| Data interpolation | obs-grid operator; residual persists, molecule-specific |
| Multiple scattering | RT null test ~1e-7 |
| FOV as broadening source | `.BIN` header: 2.36 mrad → ~0.017 cm⁻¹, negligible |
| qSDV + first-order line mixing | worsens EM27 airmass slope (well-fitting ILS); §14 |
| A-priori profile choice (XCH₄ offset, bulk airmass slope) | one fixed prior vs per-scan map: biases unchanged; §15 |
| A-priori column bookkeeping (XCH₄ offset) | GERT prior = independent .map integral to 0.7 ppb; §16 |
| CH₄ line intensity / HITRAN version (XCH₄ offset) | ∫σ dν = HITRAN2020 ΣS ratio 1.0000; iso gap only +0.88 %; §16 |
| ILS *shape* (physical NB vs Gaussian) as residual cause | width-optimized both: residual corr 0.68, same EOF1 (74/76 %, airmass −0.96); §17 |
| Missing isotopologues as cause of the airmass residual | all-iso rebuild cut offsets (XCH₄ +43→+9) but EOF1 unchanged (73 %, airmass −0.91); §18 |
| O₂ line-only (XAIR under-absorption) | +CIA (55 % of band): XAIR 0.56→1.05; §18 |
| ECS / any line-mixing (XCO₂ bandhead EOF1) | retrieval-projected corr(EOF1, qSDV+LM Δresid)=−0.10; on latest HITRAN; §19 |
| MT_CKD H₂O continuum (H₂O tilt EOF2) | Δχ²≈0 at all baseline orders (deg-0/1/2); §19 |
| O₂ speed-dependence (XAIR +5 %) | HITRAN2024 SD params all-nan for O₂; SDVoigt≡Voigt; §20 |
| Self-derived ADCF / ACOS EOF (airmass slope) | single day confounds artifact & real diurnal signal; both worsen or can't derive; §20 |
| H₂O **profile shape** as the XCO₂ diurnal drift | 2-param scale+shape (DFS 1.76 but 99 % collinear): β data-driven but χ² flat, drift −0.92→−0.96 unchanged; drift is airmass not water; §22 |
| CO spectroscopy (XCO +25–30 ppb) | CH₄-2.3 µm interferer leak — CH₄ Jacobian ~44× CO in the CO window; §21 |

## Corrections to prior assumptions (logged)
- ILS is **not** sub-% for the column (gotcha #5).
- EM27 semi-FOV is **2.36 mrad**, not 30 mrad; FOV broadening is negligible.
- COCCON L1 behaves as **effectively apodized** (Norton–Beer-like, §14); the smooth
  effective ILS is real, not merely a spectroscopy fudge. After per-band width
  fitting a physical NB ME/PE kernel and a Gaussian fit the residual near-equally
  (§17) — the ILS *shape* is immaterial to the **residual**, but shape+width set the
  **line depth**, which sets the absolute **column level** (ILS↔column degeneracy):
  the Gaussian retrieves more Xgas than the physical ILS because it makes ~1.3 %
  shallower lines. FWHM alone does *not* determine line depth — kernel shape does.
- The ~20 ppm wavenumber offset is real FTS calibration, not a bug.

## Open threads / next steps
1. **CO₂ line-shape experiment** (spec: `docs/co2_lineshape_experiment.md`) —
   controlled Voigt → qSDV → qSDV+LM on the CO₂ window. **Step 0 (Voigt) validated
   bit-identically.** Next: qSDV (needs HT/SD re-fetch), then LM (needs the
   Lamouroux/Hartmann package).
2. **XCH₄ offset** — ~RESOLVED. +67 (Gaussian) → +45 (physical ILS §17) → **+9 ppb**
   (all-isotopologue ABSCO §18). Leftover is small.
3. **O₂ CIA / XAIR** — ✅ DONE (§18): HITRAN2024 O₂-O₂+O₂-N₂ CIA folded into an
   effective σ; **XAIR 0.56 → 1.05**. Residual +5 % = O₂ line speed-dependence.
4. **Airmass-scaled per-line residual (XCO₂ EOF1, 84 %, bandheads)** — the dominant
   remaining error, but **shape hypotheses exhausted (§19)**: not line mixing
   (retrieval-projected corr −0.10, ECS not indicated), on latest HITRAN → the
   **resolution floor**. Only untested option is an *empirically-adjusted* line
   list (OCO/GGG), and it's small (XCO₂ ~0.2 ppm airmass-corrected). Likely leave.
5. **H₂O tilt residual (XCH₄ EOF2 12 %, XCO 27 %)** — identified as **H₂O-driven**
   (PC2↔H₂O 0.83/−0.98), **not** temperature, **not** the MT_CKD continuum (§19,
   Δχ²≈0 at all baseline orders). The **2-param H₂O shape** was built and tested
   (§22): the 2nd DOF is real (DFS 1.76, β data-driven) but **spectrally weak
   (99 % collinear, χ² flat)** — it does *not* fix the XCO₂ diurnal drift, which is
   **airmass** (ADCF), not water profile. Available (`retrieve_h2o_shape=`, default
   off) but not a lever on the drift. ~CLOSED.
6. **Spectroscopy engine** — the tiered ABSCO builder (now folded into the general
   `build_absco` with iso + SD/LM support) as a possible standalone repo.
7. **Empirical corrections need multi-day/multi-site data (§20).** The GERT−PROFFAST
   difference is dominated by a **<1 % constant offset (AICF)** + a loose airmass
   slope. Neither is derivable from this one day: ADCF is degenerate with the
   diurnal cycle, the **ACOS EOF correction is built** but a single-day basis eats
   CO₂ signal, and a real **AICF needs a coincident aircraft/AirCore column** (does
   one exist for GSFC ≈ 2026-04-06?). The EOF code is ready for a multi-day archive.
8. Other standing items: full ~1303-scan run for a publication time series;
   O₂ speed-dependence (needs external data — not in HITRAN, §20); the ~20 ppm FTS
   wavenumber calibration applied up front (frees dispersion for physical shifts).

## Key artifacts
- `notebooks/em27_realdata.ipynb` — M0–M5 end-to-end + residual EOF / Jacobians.
- `em27gert/retrieval.py` — `retrieve_spectrum` (obs-grid operator, ablation knobs,
  `map_path` fixed-prior override — §15).
- `scripts/run_m4.py` (`--map` fixed-prior flag), `run_m4_residuals.py`,
  `eof_residuals.py`, `eof_lineshape.py`.
- `scripts/compare_fixedprior.py` — fixed-prior vs nearest-map stats + the two
  `em27_fixedprior_{airmass,timeseries}.png` figures (§15).
- `scripts/build_co2_absco.py` — tiered CO₂ ABSCO builder (Voigt validated).
- **ILS-shape experiment (§17):** `em27gert/instrument.py::ils_physical` (NB ME/PE +
  analytic g'), `scripts/run_ils_experiment.py`, `plot_ils_experiment.py`,
  `plot_ils_shapes.py`; `eof_residuals.py` now takes `--npz/--prefix`.
- **gert additions (§17):** `ILS.response_deriv` + `ILS.convolve_derivatives`
  (closed-form ∂R/∂ν_c and ∂R/∂s), `SpectralWindow.convolve_width_derivative`,
  `width_scale` in convolve + `ForwardModel.run(ils_scale=)`,
  `StateVector.transmission_scaling(include_ils_scale=)` + `ils_scale_params`.
- **Spectroscopy build (§18):** `gert/utils/fetch_hitran.py` + `build_absco.py`
  (iso_ids / line_shape; multi-iso + SD/LM), `gert/utils/absco_spec_iso.yml`,
  `scripts/build_o2_cia.py` (CIA → effective σ, `--validate`/`--merge`).
  The all-iso+CIA result lives in `absco.h5` (rebuilt §23); the only retained
  backup is `absco.hitran2020_iso1_original.h5` (bkup3/bkup_precia deleted §23).
- **§19:** `scripts/build_h2o_continuum.py` (MT_CKD → effective σ, baseline-order
  sweep validate; **not merged** — continuum is a no-op in these dry windows),
  `plot_timeofday.py` (H₂O vs T driver of EOF2), `retrieve_spectrum(baseline_order=)`.
  Data: `data/h2o_continuum.npz`, `data/co2_absco_{voigt,sdvoigt,sdvoigt_lm}.npz`.
- **§20 (empirical corrections):** `scripts/derive_adcf.py` (+`em27_adcf_*.png`).
  ACOS **EOF radiance correction** (opt-in, default off): gert
  `GERTRetrieval(eof_basis=)` + `StateVector` `eof_{k}`/`eof_coeffs`/`include_eof`;
  em27 `build_eof_basis` + `retrieve_spectrum(eof_basis=, eof_uncert=)`. Ready for
  a multi-day archive; single-day basis confounds artifact & signal.
- **§21–22 (Xgas co-variation + 2-param H₂O):** `scripts/run_xgas_full.py`
  (4-window incl O₂ ensemble) + `plot_xgas_covary.py` (`data/xgas_full.csv`,
  `figures/em27_xgas_{timeseries_all,covary}.png`). 2-parameter H₂O (opt-in,
  default off): gert `StateVector.transmission_scaling(include_h2o_shape=,
  h2o_shape_profile=, h2o_shape_uncert=)` + `apply()`/`_jacobian_mixed` `h2o_shape`
  block + `_h2o_shape`; em27 `build_h2o_shape_profile(atm)` +
  `retrieve_spectrum(retrieve_h2o_shape=)`. O₂-pin via `freeze_gas={"o2_1p27":1.0}`.
- `run_ils_experiment.py --configs physical_fixed --res-eff` (now also saves
  `t_offset`/`h2o_scale`), `plot_residual_airmass.py`, `plot_xgas_timeseries.py`,
  `eof_residuals.py` (now **independent per-band SVD**, `--npz/--prefix`).
- `gert/instrument.py` — `SpectralWindow.obs_grid` (observation operator).
- Data: `data/m4_results.csv`, `data/m4_results_{nearestmap,fixedprior}.csv` (§15),
  `data/ils_{physical,gaussian,physical_fixed}_{results.csv,resid.npz}` (§17-18),
  `data/o2_cia.npz` (§18), `data/m4_residuals.npz`, `data/co2_absco_voigt.npz`.
