"""2013 ACC/AHA Pooled Cohort Equations for 10-year hard ASCVD risk.

Coefficients transcribed from Goff DC et al. Circulation. 2014;129(25 Suppl 2)
Appendix 7 (the 2013 ACC/AHA cardiovascular risk guideline). Predictors:
ln(age), ln(total cholesterol), ln(HDL), ln(SBP) split by treated/untreated,
current smoker, diabetes, plus group-specific interaction terms.

All coefficients, s0, and mean_sum below were cross-checked against the
`bcjaeger/PooledCohort` R package source (R/predict_risk.R, race_sex_coefs
table) and match Goff et al. 2013 Appendix 7 for all four {race, sex} groups.
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
