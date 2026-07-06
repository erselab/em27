# CO₂ line-shape experiment — spec

**Goal.** Decide, with a clean and reusable test, whether the EM27/SUN residual
we localised to the CO₂ 1.6 µm band (airmass-scaled, bandhead-concentrated,
molecule-specific) is caused by **line-shape spectroscopy** — specifically the
absence of **speed dependence (qSDV)** and **line mixing (LM)** in our plain-Voigt
ABSCO — rather than by the ILS or another forward-model term.

This is a controlled, CO₂-only experiment. It is also the **seed / prototype of a
selectable-fidelity ABSCO generator** (a possible standalone repo): the tiers
here (`voigt → sdvoigt → sdvoigt+LM`) are that project in miniature.

---

## 1. Hypothesis & prior evidence

- **H1 (spectroscopy).** GERT's plain-Voigt CO₂ cross-sections omit qSDV+LM.
  These effects are strongest at **bandheads** and **grow with optical path**,
  producing an airmass-scaled residual there.
- Prior evidence assembled this session:
  - Per-band EOF: one mode is 87 % of the CO₂-window residual variance, PC1↔airmass
    r ≈ −0.97.
  - Molecule-resolved regression: residual loads on **CO₂** (β 0.70), not the
    shared H₂O (β 0.09) → molecule-specific, not ILS.
  - Per-line ranking: worst lines are the **R-branch bandheads** (6210–6245,
    6330–6362).
  - Literature: qSDV+LM is documented to reduce airmass dependence for NIR CO₂
    (Mendonça et al.); PROFFAST v2 uses HITRAN2020 **with** CO₂ qSDV+LM, which we
    do not.
  - The physical (sinc) ILS fits our forward model ~17× worse than a Gaussian
    (χ² 23 vs 1.4) — consistent with the Gaussian absorbing a spectroscopy error
    (ILS↔spectroscopy degeneracy).

## 2. Objective & decision criteria

Build CO₂-only ABSCO blocks at three fidelity tiers, drop each into the EM27
retrieval on the **full CO₂ window (6173–6390 cm⁻¹)** across the airmass
ensemble, and decide:

- **GREEN-LIGHT the ABSCO project** if, going Voigt → qSDV+LM:
  1. the bandhead residuals (both clusters) drop materially (target: ≳50 % RMS
     reduction at the ranked worst lines), **and**
  2. the CO₂ **airmass slope** shrinks materially (target: ≳50 %), **and/or**
  3. the **physical sinc ILS** now fits the CO₂ window (χ²_red → ≲ 2) — the
     decisive degeneracy-breaker.
- **STOP / REDIRECT** if qSDV+LM does not move the bandhead residual or the
  airmass slope → the limiter is elsewhere (H₂O continuum, prior profile shape,
  residual ILS, forward-model bug).

Each tier isolates a cause: if **qSDV alone** clears most of it → speed
dependence; if the **bandheads only clear with LM** → line mixing.

## 3. Controlled design (what is held fixed)

Everything except the line shape is identical across tiers:

| Held fixed | Value |
|---|---|
| Line list | same HITRAN CO₂ lines (id 2), same isotopologues |
| Spectral window (hires) | 6173–6390 cm⁻¹ (+ wing pad), spacing 0.01 cm⁻¹ |
| p / T / (H₂O) grid | identical to the current `co2` ABSCO block |
| Broadening | air-broadened; same reference p/T |
| Downstream | same instrument (obs-grid operator), ILS, dispersion, gases, Sy |

Only `line_shape ∈ {voigt, sdvoigt, sdvoigt+LM}` changes.

## 4. Tiers to build

0. **voigt** — reproduce the *current* CO₂ result from a standalone CO₂-only
   table. (Validation gate: must match the existing block to < 0.1 %.)
1. **sdvoigt** — quadratic speed-dependent Voigt; SD width from HITRAN2020
   line-shape parameters (or a nominal SD fraction if absent). HAPI
   `absorptionCoefficient_SDVoigt`.
2. **sdvoigt+LM** — add first-order (Rosenkranz `Y`) line mixing from the
   Lamouroux/Hartmann/Tran CO₂ line-mixing package (the source OCO ABSCO uses).
   Full relaxation-matrix (ECS) is a later option, not in this experiment.

## 5. Format / integration contract

- Output: a CO₂-only `absco.h5`-compatible block — datasets
  `absco_co2_dataset (n_p, n_T, n_wn)`, `absco_co2_p`, `absco_co2_T`,
  `absco_co2_wavenumber` — matching gert's reader exactly.
- Integration: write a copy of `absco.h5` with only the `co2` block swapped
  (back up first), or point gert at a per-molecule override. gert's reader,
  forward model, and retrieval are **unchanged**.
- Units, partition functions, reference p/T, and isotopologue abundances must
  match gert's conventions (validation gate #0 catches mismatches).

## 6. Test protocol

- **A. Single clean spectrum** (`260406_163455SN`): each tier → χ²_red, per-line
  residual at the ranked worst lines, retrieved XCO₂ vs PROFFAST, dispersion
  coeffs (should be unchanged).
- **B. Ensemble** (the M4 stride-13 subsample, ~100 soundings): each tier →
  XCO₂ bias & RMS vs PROFFAST, **airmass slope**, and a per-band residual EOF
  (does the airmass-scaled mode shrink?).
- **C. ILS degeneracy test**: with the best tier, swap the Gaussian ILS for the
  **physical sinc**; does χ² collapse toward the Gaussian's? (Decisive.)

## 7. Metrics (report per tier)

- χ²_red (single + ensemble median).
- RMS residual at the top-N ranked CO₂ lines (both bandhead clusters).
- CO₂ airmass slope [ppm/airmass] and its r.
- XCO₂ − PROFFAST: bias, RMS.
- CO₂-window residual EOF1 variance fraction and PC1↔airmass correlation.
- Gaussian-vs-sinc χ² gap (degeneracy test).

## 8. Confounds & controls

- **Grid/format artifacts** → gate #0 (Voigt tier must reproduce current result).
- **Isolation** → identical line list/grid across tiers; only line shape varies.
- **Dispersion** → verify disp coeffs are unchanged across tiers (a line-shape
  change should not move the wavenumber calibration).
- **CO₂-only** → CH₄/O₂/CO untouched; conclusions limited to CO₂.

## 9. Deliverables

- `scripts/build_co2_absco.py` (tiers: `--line-shape voigt|sdvoigt|sdvoigt_lm`).
- `scripts/run_co2_lineshape_test.py` (A/B/C protocol → CSV/npz + figures).
- A short results note appended here.

## 10. Risks / open questions (resolved by feasibility checks below)

- Does gert's bundled HAPI expose `SDVoigt`/`HT` and the SD parameters?
- Are the CO₂ HITRAN lines cached, and to what wn extent?
- Exact current `co2` p/T/wn grid (for a drop-in block).
- Availability / license of the CO₂ line-mixing package (Lamouroux/Hartmann).

---

## Feasibility check results (checked 2026-07-01)

**Verdict: GO for tiers 0 (Voigt) & 1 (qSDV) now; tier 2 (LM) needs one external package.**

1. **HAPI** — v1.3.0.0 bundled at `gert/hapi.py`; exposes
   `absorptionCoefficient_SDVoigt` **and** `absorptionCoefficient_HT`. Tiers 0/1
   (Voigt, qSDV/HTP) are supported out of the box. ✓
2. **CO₂ line data** — cached (`hitran_cache/CO2_6100_6400`) covering the full
   window, **but only the standard Voigt parameter set** (`gamma_air, n_air,
   delta_air, gamma_self, …`) — **no SD/HT parameters** (`gamma2`/`gamma_SDV`/HT
   groups absent; a `line_mixing_flag` column exists but not the `Y` values).
   → qSDV requires a **re-fetch of CO₂ with the HITRAN2020 HT/SD parameter groups**.
3. **HITRAN key** — available (`HITRAN_KEY` env var set), so the re-fetch is
   possible. (Note: `fetch_hitran.py` reads `HITRAN_API_KEY`/`--key`; pass the key
   explicitly.) ✓
4. **ABSCO grid contract** — `co2` block is `(n_p=64, n_T=17, n_wn)`, **2-D
   p-dependent T grid** `T_grid[ip,it]` (117.5–363.2 K), **air-broadened only
   (no H₂O dimension)**, wn 4750–6400 cm⁻¹ @ 0.01 cm⁻¹, float32 cm²/molecule.
   Reuse the existing `absco_co2_p` / `absco_co2_T` grids and rebuild only the
   weak-band sub-range (~6100–6400). ✓
5. **Line mixing** — **not present** in gert (no relaxation-matrix / `Y` code or
   data). Tier 2 needs the external **Lamouroux/Hartmann/Tran CO₂ line-mixing
   package** (or a self-coded first-order ECS `Y`). This is the only real
   external dependency and the one blocker for tier 2.

### Recommended sequencing
- **Step 0 (now):** `build_co2_absco.py` Voigt tier → reproduce the current CO₂
  result to <0.1 % (validates grid/format/pipeline).
- **Step 1:** re-fetch CO₂ with HT/SD params (key on hand) → qSDV tier → run the
  A/B/C protocol. **qSDV alone may already move the airmass slope** — informative
  before sourcing LM.
- **Step 2:** obtain the CO₂ line-mixing package → LM tier → full comparison.

---

## Results log

### Step 0 — Voigt baseline (validation gate) ✅ PASSED (2026-07-01)
`scripts/build_co2_absco.py --line-shape voigt --workers 12 --compare`
- Reproduces the current `absco.h5` co2 block on 6100–6400 to **max |Δ/σ| = 6.7×10⁻⁸,
  RMS = 0** — bit-identical. Grid/format/line-list/Voigt engine validated.
- Built 1088 (p,T) slices × 30,001 wn in **4.3 min** (12 workers, DB shared COW,
  6100–6400 sub-range) — vs ~1.5–2 h for the naive single-core full-range build.
- 10,066 CO₂ lines (iso 1) from the `CO2_5880_6400` cache block.

### Step 1 — qSDV & qSDV+LM tables built ✅ (2026-07-02)
- **Fetch.** Re-fetched CO₂ 5880–6400 (iso 1) with `ParameterGroups=['par_line',
  'sdvoigt_air','sdvoigt_linemixing_air']` → `hitran_cache/co2_5880_6400_sdv`.
  All 10,066 lines have speed-dependence (`gamma_sdv_2`); **4,861 carry
  first-order line-mixing `y_sdv_air_296`**. So **both qSDV and line mixing come
  straight from HITRAN2020 — the external Lamouroux/Hartmann package is not
  needed** (Step-2 blocker removed). (Py 3.14 SSL fixed via `certifi`.)
- **HAPI wiring.** `absorptionCoefficient_SDVoigt` ignores `Y` by default; line
  mixing is enabled with **`LineMixingRosen=True`** (verified: zeroing `Y` then
  changes the result). Tiers: `sdvoigt` (qSDV), `sdvoigt_lm` (qSDV+LM).
- **Near-vacuum fallback.** HAPI SDVoigt raises at p→0 (Γ₀→0); 17 lowest-p slices
  fall back to Voigt — exact there since SD/LM vanish as p→0.
- **Tables:** `data/co2_absco_{voigt,sdvoigt,sdvoigt_lm}.npz` (64×17×30001, same grid).
- **Physics check (qSDV vs Voigt, surface p, 293 K):** max |Δ|/peak = 1.82 %;
  **RMS in the CO₂ bandheads = 0.199 % vs 0.035 % elsewhere (5.7×)**, largest at
  **6241/6242 and 6361/6362 cm⁻¹** — exactly the R-branch bandheads where the
  retrieval residual peaked. Line mixing adds ~0.17 % at 0.5 atm (∝ p), also
  bandhead-localized.
- **Build time:** ~8 min per SD tier (12 workers, uncontended); SDVoigt's
  `pcqsdhc` is ~20× Voigt per line.

### Step 1 A/B/C — result: **NEGATIVE.** qSDV + first-order LM does *not* fix it (2026-07-03)
In-memory CO₂-table override (`scripts/run_co2_lineshape_test.py`,
`run_co2_ensemble.py`); CO₂ window only.

- **A (χ² vs effective resolution, clean spectrum):** the χ²(res) curve is
  ~identical across tiers — all bottom out at **χ²≈1.25 @ res 0.44**. qSDV+LM does
  **not** lower the χ² minimum nor shift it toward the narrower physical
  resolution. At fixed res it *worsens* XCO₂ (Voigt +7.6 → qSDV +21 ppm).
- **B (ensemble, 99 soundings, XCO₂ airmass slope):**
  | tier | bias | RMS | airmass slope | @am=1 |
  |---|---|---|---|---|
  | voigt | +9.9 | 10.3 | **+5.4** | +1.6 |
  | sdvoigt | +26.3 | 26.9 | **+12.0** | +7.9 |
  | sdvoigt_lm | +25.6 | 26.2 | **+11.8** | +7.5 |
  qSDV+LM **doubles** the airmass slope and worsens the bias — the *opposite* of
  the literature expectation.
- **C (physical sinc ILS):** χ² stays ≈23 for all tiers — qSDV does not rescue
  the sinc fit.
- **Sanity check:** qSDV table is correct — integrated intensity conserved to
  1.0000, core 1.7 % deeper, wings lower (proper qSDV signature). The negative
  result is real, not a build/normalisation bug.

**Conclusion.** The hypothesis "residual = missing qSDV + first-order line
mixing" is **refuted**. Correctly-applied HITRAN2020 qSDV+LM makes CO₂ worse.
The controlled experiment did its job — it ruled out the cheap fix before we
built an engine on it.

### Step 1 — CORRECTION: the negative was an ILS confound; qSDV IS validated (2026-07-03)
The A/B/C tests above used the **Gaussian effective ILS**, which had been *tuned
to make Voigt fit* — and that compensation is airmass-dependent, so it biased the
comparison against qSDV. Re-running the ensemble CO₂ airmass slope with
**PROFFAST's physical ILS** (empirical ME/PE from `ils_list.csv` on the 1.8 cm
OPD-sinc + 2.36 mrad FOV; `run_co2_ensemble.py --proffast-ils`; ILS FWHM 0.333,
first sidelobe −21.4 %):

| tier | airmass slope (ppm/airmass) | r | (Gaussian ILS) |
|---|---|---|---|
| voigt | **−2.75** | −0.47 | (+5.4) |
| sdvoigt | **−0.04** | −0.01 | (+12.0) |
| sdvoigt_lm | **−0.05** | −0.01 | (+11.8) |

**With the physically-correct ILS, qSDV essentially eliminates the CO₂ airmass
slope (−2.75 → ~0)** — the literature result, and the opposite of the Gaussian-ILS
outcome. So **qSDV is the right physics for the airmass dependence**; the earlier
"negative" was entirely the Voigt-tuned Gaussian ILS masking/fighting it. Line
mixing adds little beyond qSDV for the slope.

**Caveat (separate thread):** with the bare sinc the *absolute* bias is large
(−160 ppm) and χ² stays ~20 — the COCCON L1 spectra behave *smoother* than an
unapodized sinc. That is calibration-absorbed (AICF) and orthogonal to the
airmass physics; it is the "effective-ILS shape" puzzle, investigated separately.

**Lesson.** A forward-model component (spectroscopy) can only be fairly judged
with the other components (ILS) physically correct, not tuned.

### Step 1 — FINAL (tie-breaker): the sinc result was unreliable; qSDV does NOT help (2026-07-03)
The bare-sinc reversal above was itself an artifact — the bare sinc fits the
(effectively **apodized**) COCCON L1 spectra badly (χ²≈20), so its airmass slope
is dominated by ILS mismatch, not physics. Disambiguated with a **Norton–Beer
apodized** ILS that is *both physical and well-fitting* (χ²≈1.6):

| ILS (χ²) | Voigt slope | qSDV slope | verdict |
|---|---|---|---|
| Gaussian 0.44 (≈1.3, tuned) | +5.4 | +12.0 | qSDV worse |
| **Norton–Beer (≈1.6, physical)** | **+3.0** | **+9.1** | **qSDV worse** |
| bare sinc (≈20, physical, bad fit) | −2.75 | −0.04 | (unreliable) |

**Both well-fitting ILSs agree: HITRAN2020 qSDV+LM worsens the CO₂ airmass slope
for the EM27.** The "physical-ILS validates qSDV" claim is retracted — it came
only from the badly-fitting bare sinc.

**Reconciliation.** qSDV+LM reduces airmass dependence for **high-resolution
TCCON** (0.02 cm⁻¹), where the line shape is *resolved*. At EM27 low resolution
(line width ~0.1 cm⁻¹ ≪ apodized ILS ~0.44), the line shape is unresolved: the
retrieval can neither validate it nor benefit from it, and the small
airmass-dependent term that leaks through goes the wrong way for these data. This
is why COCCON/TCCON remove the airmass artifact **empirically** (AICF + airmass
correction), not via spectroscopy, at these resolutions.

### Separate finding — the COCCON L1 effective ILS is apodized (not a bare sinc)
At matched width (~0.44 cm⁻¹): bare sinc χ²≈11–22 vs Norton–Beer / Gaussian
χ²≈1.3–1.6. A *shape* effect (apodized fast-decaying wings, no sidelobes), not
width. The `.BIN` "ILS simple" ME/PE flag is the retrieval-side model; the
spectra carry an OPUS/preprocess apodization the header does not advertise. So
the Gaussian effective ILS we used in M3–M5 was, in fact, close to the real
(apodized) effective ILS — not merely a spectroscopy fudge.

### Redirect (remaining candidates for the molecule-specific per-line residual)
1. **Full ECS relaxation-matrix line mixing** — first-order Rosenkranz `Y` is
   invalid/ wrong-signed at bandheads (exactly where our residual lives); the
   proper CO₂ line-mixing (Lamouroux/Hartmann ECS) is the one line-shape path not
   yet tested. This is now the *only* remaining line-shape hypothesis.
2. **Per-line intensity/position errors** in HITRAN2020 — the EOF attribution R²
   was only ~0.35 (rest = line-to-line variation), consistent with specific lines
   being off rather than a uniform shape model. qSDV (a smooth shape change)
   cannot fix per-line intensity errors.
3. Re-examine non-spectroscopic contributors (a-priori CO₂ profile shape, H₂O
   interference) now that the cheap spectroscopy fix is excluded.
