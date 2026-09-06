# HANDOFF — continue this build in a terminal

**Repo:** `C:\Emily\clinical-calculator-bias-audit`
**Branch:** `impl/calculator-audit` (do NOT work on `main`)
**HEAD:** `e7e19bc` — "feat: audit runner, aggregation, results writers"
**Test state:** `46 passed` (full suite), working tree clean.
**Date handed off:** 2026-09-06

The project was being built with the Superpowers `subagent-driven-development`
skill. Two process restarts interrupted it. Everything through Task 6's
implementation is committed. Task 6's REVIEW never finished. Tasks 7–9 are
not started.

---

## The plan and the ledger

- Plan (all 9 tasks, full code in each): `docs/superpowers/plans/2026-09-06-clinical-calculator-bias-audit.md`
- Design spec: `docs/superpowers/specs/2026-09-06-clinical-calculator-bias-audit-design.md`
- Progress ledger: `.superpowers/sdd/2026-09-06-clinical-calculator-bias-audit/progress.md`
- Per-task briefs + implementer reports: same `.superpowers/sdd/...` folder
  (`task-1-brief.md` … `task-6-brief.md`, `task-1-report.md` … `task-6-report.md`)

## Done (committed + reviewed clean)

| Task | Commit | What |
|---|---|---|
| 1 | `e64ad4b` | skeleton, `requirements.txt`, `.gitignore`, MIT `LICENSE`, `.github/workflows/ci.yml`, smoke test |
| 2 | `5fcc20b` | `src/calculators/egfr.py` — CKD-EPI 2009 (race) + 2021 (race-free) |
| 3 | `8a7b29b` | `src/calculators/staging.py` — `ckd_stage`, `crosses_boundary`, `statin_recommendation` |
| 4 | `7fbfe64` | `src/calculators/ascvd.py` — 2013 Pooled Cohort Equations. All 44 constants independently verified vs the `bcjaeger/PooledCohort` R package, zero transcription errors. Reference patients: white male 5.38%, black male 6.06%, white female 2.05% |
| 5 | `c7cd2e6` | `src/cohort.py` — seeded synthetic cohort (`DEFAULT_SEED = 20260906`) |

## Done (committed, REVIEW NOT FINISHED)

| Task | Commit | What |
|---|---|---|
| 6 | `e7e19bc` | `src/audit.py` — `run_audit`, `write_results`, `main()`. `tests/test_audit.py`. Full suite 46/46 |

**Resume point: review Task 6, then do Tasks 7, 8, 9.**

---

## How to resume

Open a terminal in `C:\Emily\clinical-calculator-bias-audit` and either:

### Option A — hand it back to Claude Code with the skill

Start `claude` in the repo and tell it:

> Continue executing `docs/superpowers/plans/2026-09-06-clinical-calculator-bias-audit.md`
> with the `superpowers:subagent-driven-development` skill. The ledger is at
> `.superpowers/sdd/2026-09-06-clinical-calculator-bias-audit/progress.md`.
> Tasks 1–5 are complete and reviewed. Task 6 is committed at `e7e19bc` but its
> task review never ran — do that review first (BASE `c7cd2e6`, HEAD `e7e19bc`),
> then continue with Tasks 7, 8, 9.

### Option B — finish it by hand

Remaining work, from the plan:

- **Task 6 review** — sanity-check `src/audit.py` against plan §Task 6: `d_egfr = egfr_2021 - egfr_2009_black` (must be negative for every patient), `d_risk = risk_black - risk_white`, `summary` is long-form `["metric","value"]`, threshold suffixes `0_05 / 0_075 / 0_2`, `from src.figures import …` is INSIDE `main()`.
- **Task 7 — `src/figures.py`** (plan has the full code). Six PNGs:
  `egfr_delta_hist`, `egfr_stage_heatmap`, `egfr_boundary_crossings`,
  `ascvd_delta_hist`, `ascvd_threshold_scatter`, `ascvd_statin_changes`.
  `generate_all_figures(result, outdir) -> list[Path]`. Uses Agg backend.
  Then run `python -m src.audit` end-to-end and **commit `results/`**
  (the PNGs + `summary.csv` + `headline.json` + `results_table.md` +
  `per_patient.csv` are tracked, not ignored).
- **Task 8 — `README.md`** (plan §Task 8 has the section-by-section outline).
  Paste the headline from `results/headline.json` and the table from
  `results/results_table.md` verbatim — no hand-typed numbers. Include the
  "How this was built: designed and directed by Thamer Almutairi, implementation
  AI-assisted" line, the limitations section, and the sources list.
- **Task 9 — publish.** `pytest -v && python -m src.audit` clean, then:
  ```bash
  GITHUB_TOKEN= gh repo create clinical-calculator-bias-audit \
    --public --source . --remote origin --push \
    --description "How race coefficients in eGFR and ASCVD calculators change patient-level clinical decisions — a synthetic-cohort audit."
  ```
  Then merge `impl/calculator-audit` into `main`.

---

## Environment gotchas

- **Python:** local venv is `.venv` on **Python 3.14** (3.12 wasn't installed). CI pins 3.12 and is the source of truth. Run tests with `.venv/Scripts/python -m pytest -q` (Git Bash) or `.venv\Scripts\python -m pytest -q` (PowerShell).
- **`gh` / GitHub:** an ambient `GITHUB_TOKEN` env var on this machine is **stale/invalid**. Prefix every `gh` call with `GITHUB_TOKEN=` to force it to use the keyring auth for account **`H4z305`**. Verify with `GITHUB_TOKEN= gh auth status`.
- **Do not run `python -m src.audit` until Task 7 exists** — `main()` imports `src.figures`.
- **`.venv/` and `__pycache__/` are gitignored** — don't commit them.

## Deferred minors (sweep in the final whole-branch review)

- `tests/test_egfr.py` — unused `import math` (came from the plan text).
- `tests/test_cohort.py` — unused `import pytest` (came from the plan text).
- `src/calculators/ascvd.py` — source citation is 3 lines where the plan said "one line"; content is correct.

## After the final review is clean

Delete the scratch workspace: `rm -rf .superpowers/sdd/2026-09-06-clinical-calculator-bias-audit`
(git history is the record). Then use `superpowers:finishing-a-development-branch`.

## Second repo (not started)

The résumé also needs a **Demographic Prompt-Bias Evaluation** repo. It needs an
Anthropic API key (console.anthropic.com) — Thamer said he'd make one. Its
design round hasn't happened yet; brainstorm it separately after this repo ships.
