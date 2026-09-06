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
