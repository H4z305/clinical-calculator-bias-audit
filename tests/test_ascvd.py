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
