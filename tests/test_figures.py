from src.audit import run_audit
from src.figures import generate_all_figures
from tests.test_audit import tiny_cohort

EXPECTED = {
    "egfr_delta_hist.png", "egfr_stage_heatmap.png", "egfr_boundary_crossings.png",
    "ascvd_delta_hist.png", "ascvd_threshold_scatter.png", "ascvd_statin_changes.png",
}

def test_generates_six_named_nonempty_pngs(tmp_path):
    result = run_audit(tiny_cohort())
    paths = generate_all_figures(result, tmp_path)
    assert {p.name for p in paths} == EXPECTED
    for p in paths:
        assert p.exists() and p.stat().st_size > 0

def test_is_idempotent(tmp_path):
    result = run_audit(tiny_cohort())
    generate_all_figures(result, tmp_path)
    paths = generate_all_figures(result, tmp_path)
    assert len(paths) == 6
