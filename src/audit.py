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
