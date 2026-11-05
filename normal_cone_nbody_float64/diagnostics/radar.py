"""
Compatibility Alias: Collision Radar for Singular Boundary Divisors in Q[N].

Exposes CollisionRadar and predictor utilities under the concise module name `radar`.
"""

from .collision_radar import (
    CollisionRadar,
    CollisionRadarReport,
    PairApproachPrediction,
    CollisionWarning,
    compute_pairwise_free_fall_times,
    predict_pair_approach,
    predict_closest_approach,
)

__all__ = [
    "CollisionRadar",
    "CollisionRadarReport",
    "PairApproachPrediction",
    "CollisionWarning",
    "compute_pairwise_free_fall_times",
    "predict_pair_approach",
    "predict_closest_approach",
]

