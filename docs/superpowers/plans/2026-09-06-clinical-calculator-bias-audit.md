# Clinical Calculator Bias Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a public Python repo that measures how race coefficients in the CKD-EPI eGFR equations and the Pooled Cohort Equations change patient-level clinical decisions, on a seeded synthetic cohort.

**Architecture:** Pure functions per calculator in `src/calculators/`, a seeded cohort generator, an audit module that scores every synthetic patient both ways and aggregates the reclassification statistics, and a figures module. `python -m src.audit` regenerates every number and figure in `results/`. No network, no API keys.

**Tech Stack:** Python 3.12, numpy, pandas, matplotlib, pytest, GitHub Actions.

## Global Constraints

- Python 3.12. Dependencies limited to: `numpy>=2.0,<3`, `pandas>=2.2,<3`, `matplotlib>=3.9,<4`, `pytest>=8,<9`. No other runtime dependencies. No network calls anywhere in `src/` or `tests/`.
- Package layout `src/` run as a module (`python -m src.audit`). Every module has `__init__.py`.
- Cohort seed is fixed at `20260906`. Re-running the audit must reproduce byte-identical `results/summary.csv`, `results/headline.json`, `results/results_table.md`.
- `results/*.png`, `results/*.csv`, `results/*.json`, `results/*.md` are committed to the repo (they are not build artifacts to ignore).
- Serum creatinine (`scr`) is in mg/dL. eGFR is in mL/min/1.73 m². Risk is a fraction in [0, 1].
- `sex` is the string `"female"` or `"male"`. `race` is the string `"white"` or `"black"`. Race is never a column in the cohort — it is an argument passed twice by the audit.
- LICENSE is MIT, `Copyright (c) 2026 Thamer Almutairi`.
- README must contain a "How this was built" line: designed and directed by Thamer Almutairi; implementation AI-assisted.
- No quantitative claim in the README may be hand-typed — each traces to a file in `results/`.
- Commit after every task. Use `feat:`, `test:`, `docs:`, `chore:` prefixes. End commit messages with `Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>`.

---

## File structure

| File | Responsibility |
|---|---|
| `requirements.txt` | Pinned dependency ranges (Global Constraints). |
| `.gitignore` | Python ignores only. Does **not** ignore `results/`. |
| `LICENSE` | MIT, Thamer Almutairi 2026. |
| `.github/workflows/ci.yml` | Install + `pytest -v` on push and PR, Python 3.12. |
| `src/__init__.py` | Empty package marker. |
| `src/calculators/__init__.py` | Empty package marker. |
| `src/calculators/egfr.py` | `ckd_epi_2009(scr, age, sex, black)`, `ckd_epi_2021(scr, age, sex)`. Pure. |
| `src/calculators/staging.py` | `ckd_stage(egfr)`, `crosses_boundary(a, b, boundaries)`, `statin_recommendation(risk, threshold)`. Pure. |
| `src/calculators/ascvd.py` | `pce_risk(age, total_chol, hdl, sbp, bp_treated, smoker, diabetes, sex, race)`. Pure. |
| `src/cohort.py` | `generate_cohort(n=10_000, seed=20260906) -> pandas.DataFrame`. |
| `src/audit.py` | `run_audit(cohort) -> AuditResult`; `write_results(result, outdir)`; `__main__` glue. |
| `src/figures.py` | `generate_all_figures(result, outdir) -> list[Path]`. Six PNGs. |
| `tests/test_egfr.py` | eGFR formula vs computed and externally-cross-checked values. |
| `tests/test_staging.py` | Staging bins, boundary crossing, statin threshold. |
| `tests/test_ascvd.py` | PCE vs worked examples; race-delta on the guideline example. |
| `tests/test_cohort.py` | Shape, schema, ranges, seed determinism. |
| `tests/test_audit.py` | Aggregation on a hand-built 5-row cohort; output files written. |
| `tests/test_figures.py` | Six non-empty PNGs produced. |
| `README.md` | Full write-up; results table pasted from `results/results_table.md`. |

---

## Task 1: Project skeleton, tooling, CI

**Files:**
- Create: `requirements.txt`, `.gitignore`, `LICENSE`, `.github/workflows/ci.yml`
- Create: `src/__init__.py`, `src/calculators/__init__.py`
- Create: `tests/__init__.py`, `tests/test_smoke.py`

**Interfaces:**
- Consumes: nothing.
- Produces: an importable `src` package and a green `pytest` run.

- [ ] **Step 1: Create `requirements.txt`**

```
numpy>=2.0,<3
pandas>=2.2,<3
matplotlib>=3.9,<4
pytest>=8,<9
```

- [ ] **Step 2: Create `.gitignore`**

```
__pycache__/
*.py[cod]
*.egg-info/
.pytest_cache/
.venv/
venv/
.env
.DS_Store
```

- [ ] **Step 3: Create `LICENSE`** — standard MIT text, `Copyright (c) 2026 Thamer Almutairi`.

- [ ] **Step 4: Create `.github/workflows/ci.yml`**

```yaml
name: CI
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.12"
      - run: pip install -r requirements.txt
      - run: pytest -v
```

- [ ] **Step 5: Create package markers** — `src/__init__.py`, `src/calculators/__init__.py`, `tests/__init__.py`, all empty.

- [ ] **Step 6: Write the smoke test** — `tests/test_smoke.py`

```python
def test_src_package_importable():
    import src  # noqa: F401
```

- [ ] **Step 7: Run it**

Run: `pytest -v`
Expected: PASS (1 test).

- [ ] **Step 8: Create and activate a virtualenv, install deps**

Run: `python -m venv .venv && .venv/Scripts/pip install -r requirements.txt`
(Windows Git Bash path shown; on Linux CI it is `.venv/bin/pip`.)

- [ ] **Step 9: Commit**

```bash
git add requirements.txt .gitignore LICENSE .github src tests
git commit -m "chore: project skeleton, deps, CI"
```

---

## Task 2: eGFR calculators (CKD-EPI 2009 and 2021)

**Files:**
- Create: `src/calculators/egfr.py`
- Test: `tests/test_egfr.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `ckd_epi_2009(scr: float, age: float, sex: str, black: bool) -> float`
  - `ckd_epi_2021(scr: float, age: float, sex: str) -> float`
  - Both return mL/min/1.73 m². `sex` is `"female"`/`"male"`.

- [ ] **Step 1: Write the failing tests** — `tests/test_egfr.py`

```python
import math
import pytest
from src.calculators.egfr import ckd_epi_2009, ckd_epi_2021

# Expected values computed directly from the published equations
# (Levey 2009, Ann Intern Med; Inker 2021, NEJM). Tolerance 1.0 mL/min
# absorbs rounding. Cross-check one row against the NKF calculator at
# kidney.org/professionals/kdoqi/gfr_calculator during review.

def test_2021_male_scr_0_9_age_50():
    # kappa=0.9, ratio=1.0 -> both power terms = 1; 0.9938**50 = 0.732731
    # 142 * 0.732731 = 104.05
    assert ckd_epi_2021(0.9, 50, "male") == pytest.approx(104.0, abs=1.0)

def test_2021_female_scr_0_9_age_50():
    # kappa=0.7, ratio=1.2857; max**-1.2 = 0.739587; 0.9938**50 = 0.732731
    # 142 * 0.739587 * 0.732731 * 1.012 = 77.9
    assert ckd_epi_2021(0.9, 50, "female") == pytest.approx(77.9, abs=1.0)

def test_2009_male_scr_0_9_age_50_not_black():
    # 141 * 1 * 1 * 0.993**50 (=0.703830) = 99.24
    assert ckd_epi_2009(0.9, 50, "male", black=False) == pytest.approx(99.2, abs=1.0)

def test_2009_black_multiplier_is_1_159():
    base = ckd_epi_2009(0.9, 50, "male", black=False)
    black = ckd_epi_2009(0.9, 50, "male", black=True)
    assert black == pytest.approx(base * 1.159, rel=1e-9)

def test_2009_female_scr_1_5_age_60_not_black():
    # 141 * 1 * 2.142857**-1.209 (=0.397945) * 0.993**60 (=0.656082) * 1.018 = 37.5
    assert ckd_epi_2009(1.5, 60, "female", black=False) == pytest.approx(37.5, abs=1.0)

def test_race_free_2021_has_no_black_argument():
    with pytest.raises(TypeError):
        ckd_epi_2021(0.9, 50, "male", black=True)  # type: ignore[call-arg]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_egfr.py -v`
Expected: FAIL — `ModuleNotFoundError: src.calculators.egfr`.

- [ ] **Step 3: Implement** — `src/calculators/egfr.py`

```python
"""CKD-EPI creatinine eGFR equations.

ckd_epi_2009: Levey AS et al. Ann Intern Med. 2009;150(9):604-612.
ckd_epi_2021: Inker LA et al. N Engl J Med. 2021;385:1737-1749 (race-free refit).
Serum creatinine in mg/dL; eGFR in mL/min/1.73 m^2.
"""
from __future__ import annotations


def _kappa_alpha(sex: str, alpha_female: float, alpha_male: float) -> tuple[float, float]:
    if sex == "female":
        return 0.7, alpha_female
    if sex == "male":
        return 0.9, alpha_male
    raise ValueError(f"sex must be 'female' or 'male', got {sex!r}")


def ckd_epi_2009(scr: float, age: float, sex: str, black: bool) -> float:
    kappa, alpha = _kappa_alpha(sex, -0.329, -0.411)
    ratio = scr / kappa
    egfr = (
        141.0
        * min(ratio, 1.0) ** alpha
        * max(ratio, 1.0) ** -1.209
        * 0.993 ** age
    )
    if sex == "female":
        egfr *= 1.018
    if black:
        egfr *= 1.159
    return egfr


def ckd_epi_2021(scr: float, age: float, sex: str) -> float:
    kappa, alpha = _kappa_alpha(sex, -0.241, -0.302)
    ratio = scr / kappa
    egfr = (
        142.0
        * min(ratio, 1.0) ** alpha
        * max(ratio, 1.0) ** -1.200
        * 0.9938 ** age
    )
    if sex == "female":
        egfr *= 1.012
    return egfr
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_egfr.py -v`
Expected: PASS (6 tests). If a numeric test fails outside tolerance, fix the constant in `egfr.py` — do not widen `abs=`.

- [ ] **Step 5: Commit**

```bash
git add src/calculators/egfr.py tests/test_egfr.py
git commit -m "feat: CKD-EPI 2009 and 2021 eGFR equations"
```

---

## Task 3: Decision-boundary helpers

**Files:**
- Create: `src/calculators/staging.py`
- Test: `tests/test_staging.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `ckd_stage(egfr: float) -> str` — one of `"G1","G2","G3a","G3b","G4","G5"`.
  - `EGFR_BOUNDARIES: list[float]` — `[90.0, 60.0, 45.0, 30.0, 15.0]`.
  - `crosses_boundary(a: float, b: float, boundaries: list[float]) -> bool` — True if `a` and `b` fall on opposite sides of at least one boundary.
  - `statin_recommendation(risk: float, threshold: float = 0.075) -> bool`.

- [ ] **Step 1: Write the failing tests** — `tests/test_staging.py`

```python
import pytest
from src.calculators.staging import (
    ckd_stage, crosses_boundary, statin_recommendation, EGFR_BOUNDARIES,
)

@pytest.mark.parametrize("egfr,stage", [
    (120, "G1"), (90, "G1"), (75, "G2"), (60, "G2"),
    (52, "G3a"), (45, "G3a"), (37, "G3b"), (30, "G3b"),
    (20, "G4"), (15, "G4"), (9, "G5"), (0, "G5"),
])
def test_ckd_stage_bins(egfr, stage):
    assert ckd_stage(egfr) == stage

def test_boundaries_constant():
    assert EGFR_BOUNDARIES == [90.0, 60.0, 45.0, 30.0, 15.0]

def test_crosses_boundary_true_when_straddling_60():
    assert crosses_boundary(58.0, 63.0, EGFR_BOUNDARIES) is True

def test_crosses_boundary_false_within_same_band():
    assert crosses_boundary(46.0, 59.0, EGFR_BOUNDARIES) is False

def test_crosses_boundary_order_insensitive():
    assert crosses_boundary(63.0, 58.0, EGFR_BOUNDARIES) is True

def test_statin_default_threshold():
    assert statin_recommendation(0.074) is False
    assert statin_recommendation(0.075) is True

def test_statin_custom_threshold():
    assert statin_recommendation(0.06, threshold=0.05) is True
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_staging.py -v`
Expected: FAIL — module missing.

- [ ] **Step 3: Implement** — `src/calculators/staging.py`

```python
"""Clinical decision boundaries used by the audit."""
from __future__ import annotations

EGFR_BOUNDARIES: list[float] = [90.0, 60.0, 45.0, 30.0, 15.0]


def ckd_stage(egfr: float) -> str:
    if egfr >= 90:
        return "G1"
    if egfr >= 60:
        return "G2"
    if egfr >= 45:
        return "G3a"
    if egfr >= 30:
        return "G3b"
    if egfr >= 15:
        return "G4"
    return "G5"


def crosses_boundary(a: float, b: float, boundaries: list[float]) -> bool:
    lo, hi = (a, b) if a <= b else (b, a)
    return any(lo < boundary <= hi for boundary in boundaries)


def statin_recommendation(risk: float, threshold: float = 0.075) -> bool:
    return risk >= threshold
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_staging.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/calculators/staging.py tests/test_staging.py
git commit -m "feat: CKD staging and statin decision helpers"
```

---

## Task 4: Pooled Cohort Equations (2013 ACC/AHA)

**Files:**
- Create: `src/calculators/ascvd.py`
- Test: `tests/test_ascvd.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `pce_risk(age: float, total_chol: float, hdl: float, sbp: float, bp_treated: bool, smoker: bool, diabetes: bool, sex: str, race: str) -> float` — 10-year hard-ASCVD risk as a fraction. `race` is `"white"` or `"black"`.

- [ ] **Step 1: Write the failing tests** — `tests/test_ascvd.py`

```python
import pytest
from src.calculators.ascvd import pce_risk

# Worked example widely used to validate PCE implementations:
# age 55, total chol 213, HDL 50, SBP 120 untreated, non-smoker, non-diabetic.
# ACC/AHA guideline + ACC ASCVD Risk Estimator Plus give ~5.3% (white male)
# and ~6.1% (black male). Cross-check on tools.acc.org/ascvd-risk-estimator-plus
# during review; if outside tolerance a coefficient is wrong.

BASE = dict(age=55, total_chol=213, hdl=50, sbp=120,
            bp_treated=False, smoker=False, diabetes=False)

def test_white_male_reference():
    r = pce_risk(**BASE, sex="male", race="white")
    assert r == pytest.approx(0.053, abs=0.007)

def test_black_male_reference():
    r = pce_risk(**BASE, sex="male", race="black")
    assert r == pytest.approx(0.061, abs=0.010)

def test_white_female_reference():
    r = pce_risk(**BASE, sex="female", race="white")
    assert r == pytest.approx(0.021, abs=0.007)

def test_risk_is_a_fraction():
    r = pce_risk(**BASE, sex="male", race="white")
    assert 0.0 < r < 1.0

def test_race_changes_the_number_same_covariates():
    w = pce_risk(**BASE, sex="male", race="white")
    b = pce_risk(**BASE, sex="male", race="black")
    assert w != b

def test_treated_bp_raises_risk_vs_untreated():
    lo = pce_risk(**{**BASE, "bp_treated": False}, sex="male", race="white")
    hi = pce_risk(**{**BASE, "sbp": 120, "bp_treated": True}, sex="male", race="white")
    assert hi > lo

def test_invalid_race_raises():
    with pytest.raises(KeyError):
        pce_risk(**BASE, sex="male", race="asian")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_ascvd.py -v`
Expected: FAIL — module missing.

- [ ] **Step 3: Implement** — `src/calculators/ascvd.py`

```python
"""2013 ACC/AHA Pooled Cohort Equations for 10-year hard ASCVD risk.

Coefficients transcribed from Goff DC et al. Circulation. 2014;129(25 Suppl 2)
Appendix 7 (the 2013 ACC/AHA cardiovascular risk guideline). Predictors:
ln(age), ln(total cholesterol), ln(HDL), ln(SBP) split by treated/untreated,
current smoker, diabetes, plus group-specific interaction terms.

REVIEW STEP: verify every constant below against a second published source
(e.g. the R `PooledCohort` package) and the ACC ASCVD Risk Estimator Plus
for the three reference patients in tests/test_ascvd.py before finalizing.
"""
from __future__ import annotations

import math

_MODELS: dict[tuple[str, str], dict] = {
    ("white", "female"): dict(
        s0=0.96652, mean_sum=-29.1817,
        coef=dict(
            ln_age=-29.799, ln_age_sq=4.884,
            ln_tc=13.540, ln_age_ln_tc=-3.114,
            ln_hdl=-13.578, ln_age_ln_hdl=3.149,
            ln_treated_sbp=2.019,
            ln_untreated_sbp=1.957,
            smoker=7.574, ln_age_smoker=-1.665,
            diabetes=0.661,
        ),
    ),
    ("black", "female"): dict(
        s0=0.95334, mean_sum=86.6081,
        coef=dict(
            ln_age=17.114,
            ln_tc=0.940,
            ln_hdl=-18.920, ln_age_ln_hdl=4.475,
            ln_treated_sbp=29.291, ln_age_ln_treated_sbp=-6.432,
            ln_untreated_sbp=27.820, ln_age_ln_untreated_sbp=-6.087,
            smoker=0.691,
            diabetes=0.874,
        ),
    ),
    ("white", "male"): dict(
        s0=0.91436, mean_sum=61.1816,
        coef=dict(
            ln_age=12.344,
            ln_tc=11.853, ln_age_ln_tc=-2.664,
            ln_hdl=-7.990, ln_age_ln_hdl=1.769,
            ln_treated_sbp=1.797,
            ln_untreated_sbp=1.764,
            smoker=7.837, ln_age_smoker=-1.795,
            diabetes=0.658,
        ),
    ),
    ("black", "male"): dict(
        s0=0.89536, mean_sum=19.5425,
        coef=dict(
            ln_age=2.469,
            ln_tc=0.302,
            ln_hdl=-0.307,
            ln_treated_sbp=1.916,
            ln_untreated_sbp=1.809,
            smoker=0.549,
            diabetes=0.645,
        ),
    ),
}


def pce_risk(age: float, total_chol: float, hdl: float, sbp: float,
             bp_treated: bool, smoker: bool, diabetes: bool,
             sex: str, race: str) -> float:
    model = _MODELS[(race, sex)]
    c = model["coef"]

    ln_age = math.log(age)
    ln_tc = math.log(total_chol)
    ln_hdl = math.log(hdl)
    ln_sbp = math.log(sbp)

    s = 0.0
    s += c.get("ln_age", 0.0) * ln_age
    s += c.get("ln_age_sq", 0.0) * ln_age * ln_age
    s += c.get("ln_tc", 0.0) * ln_tc
    s += c.get("ln_age_ln_tc", 0.0) * ln_age * ln_tc
    s += c.get("ln_hdl", 0.0) * ln_hdl
    s += c.get("ln_age_ln_hdl", 0.0) * ln_age * ln_hdl
    if bp_treated:
        s += c.get("ln_treated_sbp", 0.0) * ln_sbp
        s += c.get("ln_age_ln_treated_sbp", 0.0) * ln_age * ln_sbp
    else:
        s += c.get("ln_untreated_sbp", 0.0) * ln_sbp
        s += c.get("ln_age_ln_untreated_sbp", 0.0) * ln_age * ln_sbp
    s += c.get("smoker", 0.0) * (1.0 if smoker else 0.0)
    s += c.get("ln_age_smoker", 0.0) * ln_age * (1.0 if smoker else 0.0)
    s += c.get("diabetes", 0.0) * (1.0 if diabetes else 0.0)

    return 1.0 - model["s0"] ** math.exp(s - model["mean_sum"])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_ascvd.py -v`
Expected: PASS (7 tests). If a reference test fails, the coefficient set for that group is wrong — correct it against the guideline appendix and the ACC estimator; do not widen tolerances.

- [ ] **Step 5: Commit**

```bash
git add src/calculators/ascvd.py tests/test_ascvd.py
git commit -m "feat: 2013 Pooled Cohort Equations (white/black coefficient sets)"
```

---

## Task 5: Synthetic cohort generator

**Files:**
- Create: `src/cohort.py`
- Test: `tests/test_cohort.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `DEFAULT_SEED: int = 20260906`
  - `generate_cohort(n: int = 10_000, seed: int = DEFAULT_SEED) -> pandas.DataFrame` with columns exactly: `age` (int64), `sex` (object, `"female"`/`"male"`), `scr_mg_dl` (float64), `total_chol` (float64), `hdl` (float64), `sbp` (float64), `bp_treated` (bool), `smoker` (bool), `diabetes` (bool).

- [ ] **Step 1: Write the failing tests** — `tests/test_cohort.py`

```python
import pandas as pd
import pytest
from src.cohort import generate_cohort, DEFAULT_SEED

EXPECTED_COLUMNS = [
    "age", "sex", "scr_mg_dl", "total_chol", "hdl", "sbp",
    "bp_treated", "smoker", "diabetes",
]

def test_shape_and_columns():
    df = generate_cohort()
    assert len(df) == 10_000
    assert list(df.columns) == EXPECTED_COLUMNS

def test_no_race_column():
    assert "race" not in generate_cohort().columns

def test_seed_is_deterministic():
    a = generate_cohort(seed=DEFAULT_SEED)
    b = generate_cohort(seed=DEFAULT_SEED)
    pd.testing.assert_frame_equal(a, b)

def test_different_seed_differs():
    a = generate_cohort(seed=1)
    b = generate_cohort(seed=2)
    assert not a.equals(b)

def test_value_ranges():
    df = generate_cohort()
    assert df["age"].between(40, 79).all()
    assert df["scr_mg_dl"].between(0.4, 6.0).all()
    assert df["total_chol"].between(120, 320).all()
    assert df["hdl"].between(20, 100).all()
    assert df["sbp"].between(90, 200).all()
    assert set(df["sex"].unique()) <= {"female", "male"}

def test_dtypes():
    df = generate_cohort()
    assert df["age"].dtype == "int64"
    assert df["bp_treated"].dtype == "bool"
    assert df["smoker"].dtype == "bool"
    assert df["diabetes"].dtype == "bool"

def test_small_n():
    assert len(generate_cohort(n=25)) == 25
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_cohort.py -v`
Expected: FAIL — module missing.

- [ ] **Step 3: Implement** — `src/cohort.py`

```python
"""Seeded synthetic patient cohort.

Distributions are illustrative, not calibrated to any real population. The
audit measures the *within-patient* change when the race input flips, which
does not depend on the marginal distributions being population-accurate.
Race is deliberately not generated.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

DEFAULT_SEED = 20260906


def generate_cohort(n: int = 10_000, seed: int = DEFAULT_SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    sex = rng.choice(["female", "male"], size=n)
    age = rng.integers(40, 80, size=n)  # 40..79 inclusive

    scr = np.where(
        sex == "female",
        rng.lognormal(mean=np.log(0.9), sigma=0.25, size=n),
        rng.lognormal(mean=np.log(1.1), sigma=0.25, size=n),
    )
    scr = np.clip(scr, 0.4, 6.0)

    total_chol = np.clip(rng.normal(200, 35, n), 120, 320)
    hdl = np.clip(rng.normal(52, 15, n), 20, 100)
    sbp = np.clip(rng.normal(130, 18, n), 90, 200)

    bp_treated = rng.random(n) < 0.35
    smoker = rng.random(n) < 0.18
    diabetes = rng.random(n) < 0.12

    return pd.DataFrame(
        {
            "age": age.astype("int64"),
            "sex": sex,
            "scr_mg_dl": scr,
            "total_chol": total_chol,
            "hdl": hdl,
            "sbp": sbp,
            "bp_treated": bp_treated,
            "smoker": smoker,
            "diabetes": diabetes,
        }
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_cohort.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cohort.py tests/test_cohort.py
git commit -m "feat: seeded synthetic cohort generator"
```

---

## Task 6: Audit — score both ways, aggregate, write results

**Files:**
- Create: `src/audit.py`
- Test: `tests/test_audit.py`

**Interfaces:**
- Consumes:
  - `src.calculators.egfr.ckd_epi_2009`, `ckd_epi_2021`
  - `src.calculators.staging.ckd_stage`, `crosses_boundary`, `statin_recommendation`, `EGFR_BOUNDARIES`
  - `src.calculators.ascvd.pce_risk`
  - `src.cohort.generate_cohort`, `DEFAULT_SEED`
- Produces:
  - `STATIN_THRESHOLDS: list[float] = [0.05, 0.075, 0.20]`
  - `run_audit(cohort: pandas.DataFrame) -> AuditResult` where `AuditResult` is a dataclass with:
    - `per_patient: pandas.DataFrame` — original columns plus `egfr_2009_black`, `egfr_2021`, `d_egfr`, `stage_2009_black`, `stage_2021`, `egfr_stage_changed` (bool), `egfr_more_advanced_2021` (bool), `risk_white`, `risk_black`, `d_risk`, and for each threshold `t` two bool columns `statin_white_{t}`, `statin_black_{t}` and `statin_changed_{t}`.
    - `summary: pandas.DataFrame` — columns `metric`, `value`; one row per reported statistic.
    - `headline: dict` — `{"pct_reclassified": float, "sentence": str}`.
  - `write_results(result: AuditResult, outdir: pathlib.Path) -> None` — writes `summary.csv`, `headline.json`, `results_table.md`, and `per_patient.csv`.
  - Threshold column suffix format: `f"{t:g}".replace(".", "_")` → `0_05`, `0_075`, `0_2`.

- [ ] **Step 1: Write the failing tests** — `tests/test_audit.py`

```python
import json
import pandas as pd
import pytest
from src.audit import run_audit, write_results, STATIN_THRESHOLDS

def tiny_cohort():
    # Five hand-built patients spanning the interesting range.
    return pd.DataFrame({
        "age":        [50, 65, 45, 72, 60],
        "sex":        ["male", "female", "male", "female", "male"],
        "scr_mg_dl":  [0.9, 1.3, 1.1, 2.0, 1.0],
        "total_chol": [210.0, 250.0, 180.0, 230.0, 200.0],
        "hdl":        [45.0, 55.0, 40.0, 60.0, 50.0],
        "sbp":        [130.0, 145.0, 120.0, 160.0, 125.0],
        "bp_treated": [False, True, False, True, False],
        "smoker":     [False, True, False, False, True],
        "diabetes":   [False, False, False, True, False],
    })

def test_per_patient_has_expected_columns():
    r = run_audit(tiny_cohort())
    for col in ["egfr_2009_black", "egfr_2021", "d_egfr",
                "stage_2009_black", "stage_2021", "egfr_stage_changed",
                "risk_white", "risk_black", "d_risk"]:
        assert col in r.per_patient.columns
    for t in STATIN_THRESHOLDS:
        suf = f"{t:g}".replace(".", "_")
        assert f"statin_white_{suf}" in r.per_patient.columns
        assert f"statin_black_{suf}" in r.per_patient.columns
        assert f"statin_changed_{suf}" in r.per_patient.columns

def test_d_egfr_is_2021_minus_2009_black():
    r = run_audit(tiny_cohort())
    row = r.per_patient.iloc[0]
    assert row["d_egfr"] == pytest.approx(row["egfr_2021"] - row["egfr_2009_black"])

def test_race_free_2021_lower_than_2009_black_for_all_rows():
    # The 1.159 multiplier makes 2009-as-black strictly higher than race-free 2021
    # whenever the other terms match, so d_egfr should be negative for every row.
    r = run_audit(tiny_cohort())
    assert (r.per_patient["d_egfr"] < 0).all()

def test_summary_is_metric_value_long_form():
    r = run_audit(tiny_cohort())
    assert list(r.summary.columns) == ["metric", "value"]
    assert (r.summary["metric"] == "egfr_pct_stage_changed").any()
    assert (r.summary["metric"] == "ascvd_pct_cross_0_075").any()

def test_headline_sentence_contains_percentage():
    r = run_audit(tiny_cohort())
    assert 0.0 <= r.headline["pct_reclassified"] <= 100.0
    assert "%" in r.headline["sentence"]

def test_write_results_creates_all_files(tmp_path):
    r = run_audit(tiny_cohort())
    write_results(r, tmp_path)
    for name in ["summary.csv", "headline.json", "results_table.md", "per_patient.csv"]:
        assert (tmp_path / name).exists()
    payload = json.loads((tmp_path / "headline.json").read_text())
    assert "sentence" in payload and "pct_reclassified" in payload

def test_run_audit_is_deterministic():
    a = run_audit(tiny_cohort()).summary
    b = run_audit(tiny_cohort()).summary
    pd.testing.assert_frame_equal(a, b)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_audit.py -v`
Expected: FAIL — module missing.

- [ ] **Step 3: Implement** — `src/audit.py`

```python
"""Run both calculators both ways over a cohort and aggregate the effect."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from src.calculators.egfr import ckd_epi_2009, ckd_epi_2021
from src.calculators.ascvd import pce_risk
from src.calculators.staging import (
    ckd_stage, crosses_boundary, statin_recommendation, EGFR_BOUNDARIES,
)
from src.cohort import generate_cohort, DEFAULT_SEED

STATIN_THRESHOLDS = [0.05, 0.075, 0.20]
_STAGE_ORDER = ["G1", "G2", "G3a", "G3b", "G4", "G5"]


def _suffix(t: float) -> str:
    return f"{t:g}".replace(".", "_")


@dataclass
class AuditResult:
    per_patient: pd.DataFrame
    summary: pd.DataFrame
    headline: dict


def run_audit(cohort: pd.DataFrame) -> AuditResult:
    df = cohort.copy().reset_index(drop=True)

    df["egfr_2009_black"] = [
        ckd_epi_2009(s, a, x, black=True)
        for s, a, x in zip(df["scr_mg_dl"], df["age"], df["sex"])
    ]
    df["egfr_2021"] = [
        ckd_epi_2021(s, a, x)
        for s, a, x in zip(df["scr_mg_dl"], df["age"], df["sex"])
    ]
    df["d_egfr"] = df["egfr_2021"] - df["egfr_2009_black"]
    df["stage_2009_black"] = df["egfr_2009_black"].map(ckd_stage)
    df["stage_2021"] = df["egfr_2021"].map(ckd_stage)
    df["egfr_stage_changed"] = df["stage_2009_black"] != df["stage_2021"]
    df["egfr_more_advanced_2021"] = [
        _STAGE_ORDER.index(b) > _STAGE_ORDER.index(a)
        for a, b in zip(df["stage_2009_black"], df["stage_2021"])
    ]
    df["egfr_boundary_crossed"] = [
        crosses_boundary(a, b, EGFR_BOUNDARIES)
        for a, b in zip(df["egfr_2009_black"], df["egfr_2021"])
    ]

    df["risk_white"] = [
        pce_risk(a, tc, h, sbp, bt, sm, di, sex=x, race="white")
        for a, tc, h, sbp, bt, sm, di, x in zip(
            df["age"], df["total_chol"], df["hdl"], df["sbp"],
            df["bp_treated"], df["smoker"], df["diabetes"], df["sex"])
    ]
    df["risk_black"] = [
        pce_risk(a, tc, h, sbp, bt, sm, di, sex=x, race="black")
        for a, tc, h, sbp, bt, sm, di, x in zip(
            df["age"], df["total_chol"], df["hdl"], df["sbp"],
            df["bp_treated"], df["smoker"], df["diabetes"], df["sex"])
    ]
    df["d_risk"] = df["risk_black"] - df["risk_white"]

    for t in STATIN_THRESHOLDS:
        suf = _suffix(t)
        w = df["risk_white"].map(lambda r, t=t: statin_recommendation(r, t))
        b = df["risk_black"].map(lambda r, t=t: statin_recommendation(r, t))
        df[f"statin_white_{suf}"] = w
        df[f"statin_black_{suf}"] = b
        df[f"statin_changed_{suf}"] = w != b

    n = len(df)
    rows: list[tuple[str, float]] = [
        ("n_patients", float(n)),
        ("egfr_mean_d", float(df["d_egfr"].mean())),
        ("egfr_q1_d", float(df["d_egfr"].quantile(0.25))),
        ("egfr_q3_d", float(df["d_egfr"].quantile(0.75))),
        ("egfr_pct_stage_changed", 100.0 * df["egfr_stage_changed"].mean()),
        ("egfr_pct_more_advanced_2021", 100.0 * df["egfr_more_advanced_2021"].mean()),
        ("egfr_pct_boundary_crossed", 100.0 * df["egfr_boundary_crossed"].mean()),
        ("ascvd_mean_d_risk_pp", 100.0 * float(df["d_risk"].mean())),
        ("ascvd_q1_d_risk_pp", 100.0 * float(df["d_risk"].quantile(0.25))),
        ("ascvd_q3_d_risk_pp", 100.0 * float(df["d_risk"].quantile(0.75))),
    ]
    for t in STATIN_THRESHOLDS:
        suf = _suffix(t)
        changed = df[f"statin_changed_{suf}"]
        gained = changed & df[f"statin_black_{suf}"] & ~df[f"statin_white_{suf}"]
        lost = changed & df[f"statin_white_{suf}"] & ~df[f"statin_black_{suf}"]
        rows.append((f"ascvd_pct_cross_{suf}", 100.0 * changed.mean()))
        rows.append((f"ascvd_n_gained_statin_{suf}", float(gained.sum())))
        rows.append((f"ascvd_n_lost_statin_{suf}", float(lost.sum())))

    summary = pd.DataFrame(rows, columns=["metric", "value"])

    any_reclassified = (
        df["egfr_stage_changed"] | df["statin_changed_0_075"]
    )
    pct = round(100.0 * any_reclassified.mean(), 1)
    headline = {
        "pct_reclassified": pct,
        "sentence": (
            f"{pct}% of the synthetic cohort (n={n}) receives a different "
            f"clinical category — a CKD stage or a statin recommendation at the "
            f"7.5% threshold — depending only on a race coefficient."
        ),
    }

    return AuditResult(per_patient=df, summary=summary, headline=headline)


def write_results(result: AuditResult, outdir: Path) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    result.summary.to_csv(outdir / "summary.csv", index=False)
    result.per_patient.to_csv(outdir / "per_patient.csv", index=False)
    (outdir / "headline.json").write_text(
        json.dumps(result.headline, indent=2) + "\n"
    )

    lines = ["| Metric | Value |", "| --- | --- |"]
    for _, r in result.summary.iterrows():
        lines.append(f"| {r['metric']} | {r['value']:.3f} |")
    (outdir / "results_table.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    from src.figures import generate_all_figures

    cohort = generate_cohort(seed=DEFAULT_SEED)
    result = run_audit(cohort)
    results_dir = Path(__file__).resolve().parent.parent / "results"
    write_results(result, results_dir)
    generate_all_figures(result, results_dir)
    print(result.headline["sentence"])


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_audit.py -v`
Expected: PASS (7 tests).

- [ ] **Step 5: Run the full suite**

Run: `pytest -v`
Expected: PASS (all tasks so far).

- [ ] **Step 6: Commit**

```bash
git add src/audit.py tests/test_audit.py
git commit -m "feat: audit runner, aggregation, results writers"
```

---

## Task 7: Figures

**Files:**
- Create: `src/figures.py`
- Test: `tests/test_figures.py`

**Interfaces:**
- Consumes: `src.audit.AuditResult`, `src.audit.run_audit`, `src.calculators.staging.EGFR_BOUNDARIES`.
- Produces:
  - `generate_all_figures(result: AuditResult, outdir: pathlib.Path) -> list[pathlib.Path]` — writes and returns paths for exactly:
    `egfr_delta_hist.png`, `egfr_stage_heatmap.png`, `egfr_boundary_crossings.png`,
    `ascvd_delta_hist.png`, `ascvd_threshold_scatter.png`, `ascvd_statin_changes.png`.
  - Uses the non-interactive Agg backend (`matplotlib.use("Agg")` before pyplot import).

- [ ] **Step 1: Write the failing tests** — `tests/test_figures.py`

```python
from src.audit import run_audit
from src.figures import generate_all_figures
from tests.test_audit import tiny_cohort

EXPECTED = {
    "egfr_delta_hist.png", "egfr_stage_heatmap.png", "egfr_boundary_crossings.png",
    "ascvd_delta_hist.png", "ascvd_threshold_scatter.png", "ascvd_statin_changes.png",
}

def test_generates_six_named_nonempty_pngs(tmp_path):
    result = run_audit(tiny_cohort())
    paths = generate_all_figures(result, tmp_path)
    assert {p.name for p in paths} == EXPECTED
    for p in paths:
        assert p.exists() and p.stat().st_size > 0

def test_is_idempotent(tmp_path):
    result = run_audit(tiny_cohort())
    generate_all_figures(result, tmp_path)
    paths = generate_all_figures(result, tmp_path)
    assert len(paths) == 6
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_figures.py -v`
Expected: FAIL — module missing.

- [ ] **Step 3: Implement** — `src/figures.py`

```python
"""Matplotlib figures for the audit. Headless (Agg)."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from src.audit import AuditResult, STATIN_THRESHOLDS  # noqa: E402
from src.calculators.staging import EGFR_BOUNDARIES  # noqa: E402

_STAGES = ["G1", "G2", "G3a", "G3b", "G4", "G5"]


def _save(fig, outdir: Path, name: str) -> Path:
    path = Path(outdir) / name
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def generate_all_figures(result: AuditResult, outdir: Path) -> list[Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = result.per_patient
    out: list[Path] = []

    # 1. eGFR delta histogram
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["d_egfr"], bins=40, color="#3b6ea5")
    ax.set_xlabel("eGFR(2021 race-free) - eGFR(2009 as Black)  [mL/min/1.73 m^2]")
    ax.set_ylabel("patients")
    ax.set_title("Per-patient eGFR change when the race coefficient is removed")
    out.append(_save(fig, outdir, "egfr_delta_hist.png"))

    # 2. Stage heatmap (2009-as-Black rows x 2021 cols)
    ct = (
        df.groupby(["stage_2009_black", "stage_2021"]).size()
        .reindex(
            [(a, b) for a in _STAGES for b in _STAGES], fill_value=0
        )
        .unstack(fill_value=0)
        .reindex(index=_STAGES, columns=_STAGES, fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(ct.values, cmap="Blues")
    ax.set_xticks(range(len(_STAGES)), _STAGES)
    ax.set_yticks(range(len(_STAGES)), _STAGES)
    ax.set_xlabel("CKD stage under CKD-EPI 2021 (race-free)")
    ax.set_ylabel("CKD stage under CKD-EPI 2009 (as Black)")
    ax.set_title("CKD stage migration")
    for i in range(len(_STAGES)):
        for j in range(len(_STAGES)):
            v = ct.values[i, j]
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.8)
    out.append(_save(fig, outdir, "egfr_stage_heatmap.png"))

    # 3. Boundary crossings bar
    crossings = {
        b: float(
            (
                (df[["egfr_2009_black", "egfr_2021"]].min(axis=1) < b)
                & (df[["egfr_2009_black", "egfr_2021"]].max(axis=1) >= b)
            ).mean()
            * 100.0
        )
        for b in EGFR_BOUNDARIES
    }
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar([str(int(b)) for b in crossings], list(crossings.values()),
           color="#3b6ea5")
    ax.set_xlabel("eGFR boundary [mL/min/1.73 m^2]")
    ax.set_ylabel("% of cohort crossing")
    ax.set_title("Share of the cohort whose eGFR crosses each staging boundary")
    out.append(_save(fig, outdir, "egfr_boundary_crossings.png"))

    # 4. ASCVD delta histogram
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["d_risk"] * 100.0, bins=40, color="#a5643b")
    ax.set_xlabel("10-yr ASCVD risk(as Black) - risk(as White)  [percentage points]")
    ax.set_ylabel("patients")
    ax.set_title("Per-patient ASCVD risk change from the race coefficient")
    out.append(_save(fig, outdir, "ascvd_delta_hist.png"))

    # 5. Threshold scatter at 7.5%
    changed = df["statin_changed_0_075"]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(df.loc[~changed, "risk_white"] * 100,
               df.loc[~changed, "risk_black"] * 100,
               s=6, alpha=0.3, color="#888888", label="decision unchanged")
    ax.scatter(df.loc[changed, "risk_white"] * 100,
               df.loc[changed, "risk_black"] * 100,
               s=10, alpha=0.7, color="#a5643b", label="statin decision flips")
    ax.axvline(7.5, color="k", lw=0.8)
    ax.axhline(7.5, color="k", lw=0.8)
    ax.set_xlabel("10-yr risk scored as White [%]")
    ax.set_ylabel("10-yr risk scored as Black [%]")
    ax.set_title("ASCVD risk by race input, 7.5% statin threshold")
    ax.legend()
    out.append(_save(fig, outdir, "ascvd_threshold_scatter.png"))

    # 6. Statin decision changes bar
    gained = []
    lost = []
    labels = []
    for t in STATIN_THRESHOLDS:
        suf = f"{t:g}".replace(".", "_")
        c = df[f"statin_changed_{suf}"]
        g = c & df[f"statin_black_{suf}"] & ~df[f"statin_white_{suf}"]
        l = c & df[f"statin_white_{suf}"] & ~df[f"statin_black_{suf}"]
        gained.append(int(g.sum()))
        lost.append(int(l.sum()))
        labels.append(f"{t*100:g}%")
    x = range(len(labels))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar([i - 0.18 for i in x], gained, width=0.36, label="gains statin as Black",
           color="#3b6ea5")
    ax.bar([i + 0.18 for i in x], lost, width=0.36, label="loses statin as Black",
           color="#a5643b")
    ax.set_xticks(list(x), labels)
    ax.set_xlabel("statin risk threshold")
    ax.set_ylabel("patients")
    ax.set_title("Statin recommendation changes by race input")
    ax.legend()
    out.append(_save(fig, outdir, "ascvd_statin_changes.png"))

    return out
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_figures.py -v`
Expected: PASS.

- [ ] **Step 5: Run the audit end to end**

Run: `python -m src.audit`
Expected: prints the headline sentence; `results/` now contains 6 PNGs, `summary.csv`, `headline.json`, `results_table.md`, `per_patient.csv`.

- [ ] **Step 6: Sanity-check the numbers**

Open `results/summary.csv`. Confirm: `egfr_mean_d` is negative (race-free 2021 is lower than 2009-as-Black); `egfr_pct_stage_changed` is a plausible double-digit percent; `ascvd_pct_cross_0_075` is > 0. If `d_egfr` is not uniformly negative, the eGFR wiring in Task 6 is reversed — fix before committing.

- [ ] **Step 7: Commit**

```bash
git add src/figures.py tests/test_figures.py results
git commit -m "feat: audit figures; commit generated results"
```

---

## Task 8: README

**Files:**
- Create: `README.md`
- Modify: none.

**Interfaces:**
- Consumes: `results/results_table.md`, `results/headline.json`, `results/*.png`.
- Produces: the repository's front page.

- [ ] **Step 1: Read the generated outputs**

Open `results/headline.json` and `results/results_table.md`. You will paste both verbatim.

- [ ] **Step 2: Write `README.md`** with these sections in order:

1. **Title + one-paragraph summary.** What the repo measures and that it uses a synthetic cohort.
2. **Why this matters.** Two or three sentences framed on Vyas V, Eisenstein LG, Jones DS. *Hidden in Plain Sight — Reconsidering the Use of Race Correction in Clinical Algorithms.* N Engl J Med 2020;383:874-882. Link it.
3. **The two calculators.** eGFR: CKD-EPI 2009 (1.159x multiplier for Black patients) vs CKD-EPI 2021 (race-free). ASCVD: 2013 Pooled Cohort Equations scored with the White vs the Black coefficient set on identical covariates.
4. **Method.** Seeded synthetic cohort (`src/cohort.py`, n=10,000, seed 20260906), distributions stated as illustrative. Every patient scored both ways; a patient "reclassified" if their CKD stage changes or their statin recommendation at 7.5% changes.
5. **Results.** Paste the headline sentence from `headline.json`. Paste the table from `results_table.md`. Embed all six figures with `![...](results/<name>.png)` and a one-line caption each.
6. **How these compare to published cohorts.** State the real-world numbers (US Military Health System: Black adults at CKD stage 3-5 rose 58%, 45.8% reclassified to a more advanced stage; PCE -> PREVENT: 21.7% of non-Hispanic Black adults with PCE >= 7.5% fall below the PREVENT 5% threshold), link the two sources below, and say explicitly: *the synthetic results move in the same direction and rough magnitude; this is a plausibility check, not a validation of the synthetic cohort.*
7. **Limitations.** Bullet list, verbatim intent from spec section 8.6: synthetic data / illustrative distributions; not clinical validation, not for clinical use; race operationalized as a crude binary toggle by design because that is what the audited equations do, and PCE has no other categories; the eGFR "% reclassified" describes the affected subgroup (patients a 2009-era clinic would have recorded as Black), an upper bound on population impact.
8. **Run it.**

   ```bash
   python -m venv .venv && .venv/bin/pip install -r requirements.txt
   .venv/bin/python -m src.audit    # regenerates everything in results/
   .venv/bin/pytest -v
   ```

9. **Sources.** Bulleted links:
   - Vyas et al., NEJM 2020 — https://www.nejm.org/doi/full/10.1056/NEJMms2004740
   - Inker et al., CKD-EPI 2021, NEJM — https://www.nejm.org/doi/full/10.1056/NEJMoa2102953
   - Goff et al., 2013 ACC/AHA risk guideline (PCE), Circulation — https://www.ahajournals.org/doi/10.1161/01.cir.0000437741.48606.98
   - Race-free eGFR impact, US Military Health System — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11295453/
   - Statin eligibility, PCE -> PREVENT — https://www.sciencedirect.com/science/article/pii/S1933287426004319
10. **How this was built.** "Designed and directed by Thamer Almutairi. Implementation was AI-assisted under his direction. Not peer-reviewed; not medical advice."
11. **License.** MIT.

- [ ] **Step 3: Verify no hand-typed numbers**

Re-read the README. Every number in sections 5 is copied from `results/`. Section 6's numbers are the cited external figures only. If you typed a statistic anywhere else, remove it or source it.

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "docs: project README with results, limitations, sources"
```

---

## Task 9: Publish to GitHub

**Files:** none (repo operations only).

**Interfaces:**
- Consumes: a working local repo with a green `pytest`.
- Produces: a public GitHub repo `H4z305/clinical-calculator-bias-audit` with CI passing.

- [ ] **Step 1: Final local verification**

Run: `pytest -v && python -m src.audit`
Expected: all tests pass; audit reprints the headline; `git status` shows a clean tree (results already committed).

- [ ] **Step 2: Confirm the correct gh account**

Run: `GITHUB_TOKEN= gh auth status`
Expected: active account `H4z305`. The ambient `GITHUB_TOKEN` env var is stale and must be blanked for every `gh` call in this task (`GITHUB_TOKEN= gh ...`).

- [ ] **Step 3: Create the remote and push**

```bash
GITHUB_TOKEN= gh repo create clinical-calculator-bias-audit \
  --public --source . --remote origin --push \
  --description "How race coefficients in eGFR and ASCVD calculators change patient-level clinical decisions — a synthetic-cohort audit."
```

- [ ] **Step 4: Verify CI**

Run: `GITHUB_TOKEN= gh run watch --exit-status` (or check the Actions tab).
Expected: the CI workflow passes on the pushed commit.

- [ ] **Step 5: Confirm the README renders**

Run: `GITHUB_TOKEN= gh repo view --web`
Expected: figures display, results table populated, links resolve.

- [ ] **Step 6: Record the URL**

Note the repo URL for the résumé. No commit needed.

---

## Self-Review

**1. Spec coverage**

| Spec section | Task(s) |
|---|---|
| §2 eGFR 2009 vs 2021 | Task 2 |
| §2 ASCVD White vs Black | Task 4 |
| §2 seeded cohort | Task 5 |
| §2 audit flips only race | Task 6 |
| §2 6 figures + summary.csv + headline.json + results_table.md | Tasks 6, 7 |
| §2 pytest vs published references | Tasks 2, 4 (+ 3, 5, 6, 7 unit tests) |
| §2 GitHub Actions | Task 1 |
| §2 README | Task 8 |
| §3 stack / layout / MIT / .gitignore | Task 1 |
| §4 formulas | Tasks 2, 4 |
| §5 reference values + cross-check steps | Tasks 2, 4 (review steps) |
| §6 cohort fields & distributions | Task 5 |
| §7 per-patient columns, summary rows, headline, figures, plausibility check | Tasks 6, 7, 8 |
| §8 README arc incl. "How this was built" | Task 8 |
| §9 acceptance criteria | Tasks 7 (end-to-end), 9 (CI), plus per-task tests |
| §10 risks (wrong coefficients, scope creep, over-claiming) | Tasks 2/4 review steps; PREVENT/spirometry kept out; Task 8 step 3 + limitations |

No spec requirement is left without a task. PREVENT and spirometry are intentionally absent (stretch / out of scope per spec).

**2. Placeholder scan**

- All test bodies contain real assertions with concrete expected values. eGFR expected values are computed from the equations in-line in comments; ASCVD expected values carry an explicit "cross-check on tools.acc.org" review step because transcribed coefficients must be verified against an authority, not trusted from this plan.
- No "TBD", "handle errors appropriately", or "similar to Task N" left in the plan.

**3. Type consistency**

- `sex` is `"female"`/`"male"` everywhere (egfr, ascvd, cohort, audit).
- `race` is `"white"`/`"black"` everywhere; never a cohort column.
- `ckd_epi_2009(scr, age, sex, black)` / `ckd_epi_2021(scr, age, sex)` — same signatures in Task 2 interface, Task 2 tests, and Task 6 consumption.
- `pce_risk(age, total_chol, hdl, sbp, bp_treated, smoker, diabetes, sex, race)` — identical in Task 4 and Task 6.
- Threshold suffix rule `f"{t:g}".replace(".", "_")` → `0_05`, `0_075`, `0_2` — defined in Task 6 interface, used identically in Tasks 6 and 7 and `tests/test_audit.py`.
- `AuditResult` fields (`per_patient`, `summary`, `headline`) — defined in Task 6, consumed by Task 7 and Task 8 exactly as named.
- `generate_all_figures(result, outdir) -> list[Path]` — same in Task 7 interface, tests, and `src/audit.py:main`.

Fixed inline: `tests/test_figures.py` imports `tiny_cohort` from `tests.test_audit`, which requires `tests/__init__.py` — already created in Task 1 Step 5.
