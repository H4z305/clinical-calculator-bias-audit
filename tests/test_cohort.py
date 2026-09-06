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
