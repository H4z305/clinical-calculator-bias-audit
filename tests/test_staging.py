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
