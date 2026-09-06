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
