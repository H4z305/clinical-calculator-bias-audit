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
