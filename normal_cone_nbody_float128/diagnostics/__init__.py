"""
Diagnostics, Symplectic Metrology, Collision Radar, and Publication Visualizer Package.

Exposes:
- CollisionRadar, CollisionRadarReport, PairApproachPrediction, predict_pair_approach
- CelestialKPIs, KPIReport, compute_celestial_kpis, compute_energy_secular_slope
- PublicationVisualizer
"""

from .collision_radar import (
    PairApproachPrediction,
    CollisionRadarReport,
    predict_pair_approach,
    CollisionRadar,
)
from .metrics import (
    CelestialKPIs,
    KPIReport,
    compute_energy_secular_slope,
    compute_celestial_kpis,
    compute_osculating_elements,
    compute_laplace_resonant_angle,
)
from .visualizer import (
    PublicationVisualizer,
)
from .profiler import (
    ComputationalProfiler,
    StageCostMetrics,
    MacroStepProfile,
)

__all__ = [
    "PairApproachPrediction",
    "CollisionRadarReport",
    "predict_pair_approach",
    "CollisionRadar",
    "CelestialKPIs",
    "KPIReport",
    "compute_energy_secular_slope",
    "compute_celestial_kpis",
    "PublicationVisualizer",
    "ComputationalProfiler",
    "StageCostMetrics",
    "MacroStepProfile",
]
