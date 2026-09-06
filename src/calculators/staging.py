"""Clinical decision boundaries used by the audit."""
from __future__ import annotations

EGFR_BOUNDARIES: list[float] = [90.0, 60.0, 45.0, 30.0, 15.0]


def ckd_stage(egfr: float) -> str:
    if egfr >= 90:
        return "G1"
    if egfr >= 60:
        return "G2"
    if egfr >= 45:
        return "G3a"
    if egfr >= 30:
        return "G3b"
    if egfr >= 15:
        return "G4"
    return "G5"


def crosses_boundary(a: float, b: float, boundaries: list[float]) -> bool:
    lo, hi = (a, b) if a <= b else (b, a)
    return any(lo < boundary <= hi for boundary in boundaries)


def statin_recommendation(risk: float, threshold: float = 0.075) -> bool:
    return risk >= threshold
