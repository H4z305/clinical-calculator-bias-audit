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
