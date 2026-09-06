# Clinical Calculator Bias Audit — Design Spec

**Date:** 2026-09-06
**Owner / director:** Thamer Almutairi
**Implementation:** AI-assisted, under the owner's direction
**Status:** approved design, pre-implementation

---

## 1. Purpose

Quantify how race coefficients in two widely-used clinical calculators change
**patient-level clinical decisions**, using a synthetic cohort. A small,
reproducible methods demonstration that backs the finding of the owner's scoping
review (*AI miscalculation and demographic gaps in healthcare*).

The output is a public GitHub repository: implemented formulas, a seeded
synthetic cohort, an analysis that flips only the race input, committed figures,
and a README that reports the results and compares them to published real-cohort
numbers.

Not a research contribution, not clinical validation. A portfolio artifact that
shows the owner can frame a fairness question, direct its implementation, and
report it honestly.

## 2. Scope

### In

- **eGFR:** CKD-EPI 2009 creatinine equation (with the Black-race multiplier)
  vs CKD-EPI 2021 creatinine equation (race-free).
- **ASCVD:** 2013 ACC/AHA Pooled Cohort Equations (PCE), scoring each synthetic
  patient with the "White" coefficient set and again with the "Black"
  ("African American") coefficient set.
- Seeded synthetic cohort generator.
- Audit that computes, per patient, both versions of each score; the delta; and
  whether the patient crosses a clinical decision boundary.
- Figures (committed PNGs) and a `summary.csv` of reclassification statistics.
- pytest validation of each formula against published reference values.
- GitHub Actions workflow running pytest.
- README with methodology, results, comparison to published cohorts,
  limitations, run instructions, sources.

### Out (non-goals)

- No real patient data.
- No web UI, no interactive dashboard.
- No spirometry (GLI-2012 vs GLI Global 2022) — needs spline lookup tables,
  disproportionate work. README lists it as a planned extension.
- No PyPI packaging.
- **Stretch only, after the two core calculators pass:** add the 2023 AHA
  PREVENT equations as the "race-free replacement" comparison for ASCVD.

## 3. Stack & layout

- Python 3.12. `numpy`, `pandas` for the cohort and audit; `matplotlib` for
  figures; `pytest` for tests. No network, no API keys.
- `requirements.txt` (pinned minor versions). Plain `src/` package, run as
  `python -m src.audit`.
- MIT licence, `LICENSE` attributed to Thamer Almutairi, 2026.
- `.gitignore` for Python.

```
clinical-calculator-bias-audit/
  README.md
  LICENSE
  requirements.txt
  .gitignore
  .github/workflows/ci.yml
  src/
    __init__.py
    calculators/
      __init__.py
      egfr.py        # ckd_epi_2009(scr, age, sex, black), ckd_epi_2021(scr, age, sex)
      ascvd.py       # pce_risk(covariates, race)  -> 10-year risk fraction
      staging.py     # ckd_stage(egfr), statin_recommendation(risk, threshold)
    cohort.py        # generate_cohort(n, seed) -> DataFrame
    audit.py         # run everything, write results/
    figures.py       # build the plots from the audit output
  tests/
    test_egfr.py
    test_ascvd.py
    test_staging.py
  results/           # committed: *.png, summary.csv, headline.json
  docs/superpowers/specs/2026-09-06-clinical-calculator-bias-audit-design.md
```

## 4. Formulas

### 4.1 CKD-EPI 2009 (creatinine) — Levey et al., Ann Intern Med 2009

```
eGFR = 141
     * min(Scr/kappa, 1) ** alpha
     * max(Scr/kappa, 1) ** (-1.209)
     * 0.993 ** age
     * (1.018 if female)
     * (1.159 if Black)
```
`kappa` = 0.7 (female) / 0.9 (male); `alpha` = -0.329 (female) / -0.411 (male).
`Scr` in mg/dL. Units: mL/min/1.73 m^2.

### 4.2 CKD-EPI 2021 (creatinine, race-free) — Inker et al., NEJM 2021

```
eGFR = 142
     * min(Scr/kappa, 1) ** alpha
     * max(Scr/kappa, 1) ** (-1.200)
     * 0.9938 ** age
     * (1.012 if female)
```
`kappa` = 0.7 (female) / 0.9 (male); `alpha` = -0.241 (female) / -0.302 (male).

### 4.3 Pooled Cohort Equations — Goff et al., Circulation 2013 (2013 ACC/AHA guideline, Appendix 7)

Four coefficient sets: {White, African American} x {female, male}. Predictors:
`ln(age)`, `ln(total cholesterol)`, `ln(HDL)`, `ln(treated SBP)` or
`ln(untreated SBP)`, current smoker (0/1), diabetes (0/1), plus the
group-specific interaction terms (`ln(age)*ln(total chol)`,
`ln(age)*ln(HDL)`, `ln(age)*ln(treated SBP)` / `ln(age)*ln(untreated SBP)`,
`ln(age)*smoker`, and for one group `ln(age)^2`).

```
risk = 1 - S0 ** exp( sum(beta*x) - mean_sum )
```
`S0` (baseline 10-year survival) and `mean_sum` are the published per-group
constants. Coefficients transcribed from the 2013 guideline Appendix 7 table and
cross-checked against the ACC ASCVD Risk Estimator Plus.

Race is supplied as an argument (`"white"` / `"black"`); the audit calls the
function twice with the same covariates and both race values.

## 5. Reference values for tests

Each calculator validated against **3-5 published input/output pairs**, tolerance
+/- 1 unit for eGFR (rounding in published tools) and +/- 0.3 percentage points
for PCE risk. Sources for the vectors, pinned during implementation:

- **eGFR:** National Kidney Foundation eGFR calculator worked values; Inker 2021
  supplementary examples.
- **PCE:** worked examples in the 2013 ACC/AHA guideline; ACC ASCVD Risk
  Estimator Plus for a small grid of inputs.

Test vectors are hard-coded in the test files with a comment citing the source
for each.

## 6. Synthetic cohort

`generate_cohort(n=10_000, seed=20260906)` -> `pandas.DataFrame`, one row per
patient. **Race is not generated** — it is the toggle applied later.

| field | distribution (illustrative, not epidemiologically calibrated) |
|---|---|
| `age` | uniform integer 40-79 |
| `sex` | Bernoulli 0.5 -> {"female","male"} |
| `scr_mg_dl` | lognormal, median ~0.9 (female) / ~1.1 (male), clipped to [0.4, 6.0] |
| `total_chol` | normal mean 200 sd 35, clip [120, 320] |
| `hdl` | normal mean 52 sd 15, clip [20, 100] |
| `sbp` | normal mean 130 sd 18, clip [90, 200] |
| `bp_treated` | Bernoulli 0.35 |
| `smoker` | Bernoulli 0.18 |
| `diabetes` | Bernoulli 0.12 |

Distributions are documented in the README as illustrative. The point of the
audit is the *within-patient* change when race flips, which does not depend on
the cohort being a calibrated population sample. A seed makes every number
reproducible.

## 7. Audit & outputs

For every patient:

- **eGFR:** `egfr_2009_black = ckd_epi_2009(..., black=True)`,
  `egfr_2021 = ckd_epi_2021(...)`. Record `d_egfr = egfr_2021 - egfr_2009_black`,
  the CKD stage under each, and whether the stage changed.
- **ASCVD:** `risk_white = pce_risk(..., race="white")`,
  `risk_black = pce_risk(..., race="black")`. Record `d_risk`, the statin
  recommendation under each at thresholds {5%, 7.5%, 20%}, and whether it changed.

### Reported statistics (`results/summary.csv`, one row per metric)

- eGFR: mean and IQR of `d_egfr`; % of cohort crossing each boundary
  (90, 60, 45, 30, 15); count moving to a **more advanced** CKD stage under the
  race-free 2021 equation.
- ASCVD: mean and IQR of `d_risk`; % crossing 7.5% (and 5%, 20%); count
  **losing** and **gaining** a statin recommendation when scored as Black vs White.
- `results/headline.json`: the single sentence figure —
  `"<pct>% of the synthetic cohort receives a different clinical category
  depending on a race coefficient"` (union across both calculators).
- `results/results_table.md`: a Markdown table of the same statistics, written by
  `audit.py` and pasted verbatim into the README so the two cannot drift.

### Figures (`results/*.png`, committed)

1. `egfr_delta_hist.png` — histogram of `d_egfr`.
2. `egfr_stage_heatmap.png` — CKD stage under 2009-Black (rows) vs 2021 (cols),
   cell counts.
3. `egfr_boundary_crossings.png` — bar: % of cohort crossing each eGFR boundary.
4. `ascvd_delta_hist.png` — histogram of `d_risk`.
5. `ascvd_threshold_scatter.png` — `risk_white` vs `risk_black`, 7.5% lines,
   points coloured by "statin decision changed".
6. `ascvd_statin_changes.png` — bar: gained / lost statin recommendation at each
   threshold.

All figures regenerated by `python -m src.audit` (which calls `figures.py`).
No figure or number in the README is hand-typed.

### Plausibility check (README prose, not a test)

Compare the synthetic reclassification rates to published real-cohort numbers
(US Military Health System: +58% Black adults at CKD stage 3-5, 45.8%
reclassified more severe; PCE -> PREVENT: 21.7% of non-Hispanic Black adults
with PCE >= 7.5% fall below the PREVENT 5% threshold). State that the direction
and rough magnitude line up, and that this is a plausibility check, not a
validation of the synthetic cohort.

## 8. README arc

1. What this is / why it matters — framed on Vyas, Eisenstein, Jones,
   *"Hidden in Plain Sight — Reconsidering the Use of Race Correction in
   Clinical Algorithms,"* NEJM 2020.
2. The two calculators and the exact equation versions compared.
3. Method — synthetic cohort (seeded), the race toggle, decision boundaries.
4. Results — embedded figures, filled results table, the headline sentence.
5. Our numbers vs published cohorts (the plausibility check).
6. Limitations —
   - synthetic data, distributions illustrative, not a population sample;
   - not clinical validation, not for clinical use;
   - race is deliberately operationalized as a crude binary toggle (White/Black)
     because that is exactly what the audited equations do; PCE has no other
     categories;
   - the eGFR "% reclassified" describes the affected subgroup (patients a
     2009-era clinic would have recorded as Black), an upper bound on
     population-level impact.
7. How to run.
8. Sources.
9. "How this was built" — designed and directed by Thamer Almutairi;
   implementation with AI assistance.

## 9. Acceptance criteria

- `pip install -r requirements.txt && python -m src.audit` runs end to end,
  deterministic, writes `results/` (figures + `summary.csv` + `headline.json`).
- `pytest` green: every calculator within tolerance of its published reference
  values; staging/threshold helpers unit-tested.
- README renders on GitHub with all six figures embedded and the results table
  populated from `summary.csv`.
- `.github/workflows/ci.yml` passes on push and pull request (Python 3.12,
  install, pytest).
- Re-running the audit reproduces identical `summary.csv` and `headline.json`
  (seed-locked).
- Every quantitative claim in the README traces to a file in `results/`.

## 10. Risks

- **Wrong coefficients** — the PCE table is error-prone to transcribe. Mitigated
  by the reference-value tests against the ACC estimator; if a vector fails, fix
  the coefficient, not the tolerance.
- **Scope creep into PREVENT / spirometry** — held as explicit stretch/out.
- **Over-claiming** — the plausibility check must not be written as validation;
  the limitations section is load-bearing.
