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
