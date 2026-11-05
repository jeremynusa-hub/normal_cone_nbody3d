"""
Predictive Collision Radar and Exceptional Divisor Proximity Anticipation.

Theoretical Foundations:
- Monograph Reference: Chapter 6 (Exceptional Collision Divisors E_I in partial Q[N]),
  and Section 4.3 (Intrinsic Compactness and Marchal-Sundman Threshold).
- Divisor Proximity Metric:
    d(q, E_I) ~ min_{i < j in I} ||q_i - q_j||.
- Ballistic & Keplerian Time-to-Collision Extrapolation:
    For any pair (i, j), with relative separation r = ||q_i - q_j|| and relative velocity v = v_i - v_j:
    Radial velocity: v_r = (r . v) / r.
    Approach predicate: is_approaching <=> v_r < 0.
    First-order ballistic impact time: t_approach = r / |v_r|.
    Osculating Keplerian periastron: r_peri = h^2 / [mu (1 + e)].
- Normalized Singularity Hazard Score (Convex Multi-Scale Combination):
    sigma = clamp(0.6 * exp(-min_{i<j} r_{ij} / r_ref) + 0.4 * [1 / (1 + t_approach / t_ref)], 0, 1).
- Unsoftened Geometry: Exact collision dynamics (eps = 0, d0 = 0).
- Certified Formal Proofs:
    Lean 4 modules `01_PartitionLattice_Defs.lean` (`collision_diagonal_non_transversality`) and
    `02_LieBargmannExactInvariance_Theorems.lean`.
"""

from __future__ import annotations
import math
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Any, Optional
import numpy as np

from ..core.constants import GRAVITATIONAL_CONSTANT
from ..core.state import ExtendedPhaseState


@dataclass(frozen=True)
class PairApproachPrediction:
    """
    Predictive close encounter indicators for an interacting pair (i, j).
    """
    body_i: int
    body_j: int
    r_current: float
    radial_velocity: float
    is_approaching: bool
    time_to_approach: float
    predicted_r_peri: float
    osculating_eccentricity: float
    osculating_energy: float

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pair": (self.body_i, self.body_j),
            "r_current": float(self.r_current),
            "radial_velocity": float(self.radial_velocity),
            "is_approaching": bool(self.is_approaching),
            "time_to_approach": float(self.time_to_approach),
            "predicted_r_peri": float(self.predicted_r_peri),
            "osculating_eccentricity": float(self.osculating_eccentricity),
            "osculating_energy": float(self.osculating_energy),
        }


@dataclass(frozen=True)
class CollisionRadarReport:
    """
    Global instantaneous collision anticipation report across all N bodies.
    """
    t_current: float
    min_separation: float
    closest_pair: Tuple[int, int]
    imminent_collision: bool
    hazard_score: float
    pair_predictions: List[PairApproachPrediction]
    critical_pair: Optional[PairApproachPrediction] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "t_current": float(self.t_current),
            "min_separation": float(self.min_separation),
            "closest_pair": self.closest_pair,
            "imminent_collision": bool(self.imminent_collision),
            "hazard_score": float(self.hazard_score),
            "critical_pair": self.critical_pair.to_dict() if self.critical_pair else None,
            "n_approaching_pairs": len([p for p in self.pair_predictions if p.is_approaching]),
        }


def predict_pair_approach(
    q_i: np.ndarray,
    q_j: np.ndarray,
    v_i: np.ndarray,
    v_j: np.ndarray,
    m_i: float,
    m_j: float,
    body_i: int = 0,
    body_j: int = 1,
    G: float = GRAVITATIONAL_CONSTANT,
) -> PairApproachPrediction:
    """
    Computes exact osculating two-body encounter metrics between body i and body j.

    Parameters:
        q_i, q_j: Cartesian coordinates (3,).
        v_i, v_j: Cartesian velocities (3,).
        m_i, m_j: Masses of the bodies.
        body_i, body_j: Indices of the bodies.
        G: Gravitational constant.

    Returns:
        PairApproachPrediction instance.
    """
    r_vec = q_i - q_j
    r = float(np.linalg.norm(r_vec))
    v_rel = v_i - v_j
    v_sq = float(np.dot(v_rel, v_rel))

    if r <= 1e-30:
        # Exact coincidence collision
        return PairApproachPrediction(
            body_i=body_i,
            body_j=body_j,
            r_current=0.0,
            radial_velocity=0.0,
            is_approaching=True,
            time_to_approach=0.0,
            predicted_r_peri=0.0,
            osculating_eccentricity=1.0,
            osculating_energy=float('-inf'),
        )

    r_dot_v = float(np.dot(r_vec, v_rel))
    v_r = r_dot_v / r
    is_approaching = (v_r < 0.0)

    # Gravitational parameter mu = G(m_i + m_j)
    mu = G * (m_i + m_j)

    # Specific angular momentum vector h = r x v
    h_vec = np.cross(r_vec, v_rel)
    h_sq = float(np.dot(h_vec, h_vec))

    # Osculating specific energy E = 0.5 * v^2 - mu / r
    energy_spec = 0.5 * v_sq - (mu / r)

    # Osculating eccentricity e = sqrt(max(0, 1 + 2 * E * h^2 / mu^2))
    if mu > 1e-30:
        disc = 1.0 + (2.0 * energy_spec * h_sq) / (mu**2)
        ecc = math.sqrt(max(0.0, disc))
        if h_sq > 1e-25:
            r_peri = h_sq / (mu * (1.0 + ecc))
        else:
            r_peri = 0.0
    else:
        ecc = 1.0
        r_peri = r

    # Time-to-approach calculation
    if is_approaching:
        if abs(v_r) > 1e-12:
            time_to_approach = r / abs(v_r)
        elif mu > 1e-30:
            # Free-fall time from rest under pure gravity
            time_to_approach = (math.pi / (2.0 * math.sqrt(2.0))) * math.sqrt(r**3 / mu)
        else:
            time_to_approach = float('inf')
    else:
        time_to_approach = float('inf')

    return PairApproachPrediction(
        body_i=body_i,
        body_j=body_j,
        r_current=r,
        radial_velocity=v_r,
        is_approaching=is_approaching,
        time_to_approach=time_to_approach,
        predicted_r_peri=r_peri,
        osculating_eccentricity=ecc,
        osculating_energy=energy_spec,
    )


class CollisionRadar:
    """
    Anticipates proximity to Fulton--MacPherson exceptional collision divisors E_I
    and ranks impending singular close encounters across all pairs.
    """

    def __init__(
        self,
        G: float = GRAVITATIONAL_CONSTANT,
        hazard_distance_threshold: float = 0.1,
        hazard_time_threshold: float = 1.0,
    ) -> None:
        """
        Initializes the collision radar.

        Parameters:
            G: Gravitational constant.
            hazard_distance_threshold: Proximity distance scale for high risk score.
            hazard_time_threshold: Time scale for imminent collision alert.
        """
        self.G = float(G)
        self.r_ref = float(hazard_distance_threshold)
        self.t_ref = float(hazard_time_threshold)

    def scan(self, state: ExtendedPhaseState) -> CollisionRadarReport:
        """
        Performs an instantaneous all-to-all pair scan of the system state.

        Parameters:
            state: ExtendedPhaseState instance.

        Returns:
            CollisionRadarReport containing pairwise predictions and risk score.
        """
        N = state.N
        velocities = state.velocities
        predictions: List[PairApproachPrediction] = []

        min_sep = float('inf')
        closest_pair = (0, 1 if N > 1 else 0)
        critical_pred: Optional[PairApproachPrediction] = None
        min_time = float('inf')

        for i in range(N):
            for j in range(i + 1, N):
                pred = predict_pair_approach(
                    q_i=state.q[i],
                    q_j=state.q[j],
                    v_i=velocities[i],
                    v_j=velocities[j],
                    m_i=float(state.masses[i]),
                    m_j=float(state.masses[j]),
                    body_i=i,
                    body_j=j,
                    G=self.G,
                )
                predictions.append(pred)

                if pred.r_current < min_sep:
                    min_sep = pred.r_current
                    closest_pair = (i, j)

                if pred.is_approaching and pred.time_to_approach < min_time:
                    min_time = pred.time_to_approach
                    critical_pred = pred

        if critical_pred is None and predictions:
            # Fallback to closest pair prediction if none are actively approaching
            for p in predictions:
                if (p.body_i, p.body_j) == closest_pair:
                    critical_pred = p
                    break

        # Compute normalized singularity hazard score in [0, 1]
        imminent = (min_sep < self.r_ref) or (min_time < self.t_ref)
        spatial_factor = math.exp(-max(0.0, min_sep) / max(1e-12, self.r_ref))
        temporal_factor = 1.0 / (1.0 + max(0.0, min_time) / max(1e-12, self.t_ref)) if min_time < float('inf') else 0.0
        hazard_score = float(max(0.0, min(1.0, 0.6 * spatial_factor + 0.4 * temporal_factor)))

        return CollisionRadarReport(
            t_current=float(state.t),
            min_separation=min_sep,
            closest_pair=closest_pair,
            imminent_collision=imminent,
            hazard_score=hazard_score,
            pair_predictions=predictions,
            critical_pair=critical_pred,
        )


def compute_pairwise_free_fall_times(
    state: ExtendedPhaseState,
    G: float = GRAVITATIONAL_CONSTANT,
) -> Dict[Tuple[int, int], float]:
    """
    Computes the classical two-body gravitational free-fall time from rest
    t_ff(i, j) = (pi / (2 * sqrt(2))) * sqrt(r_ij^3 / (G * (m_i + m_j)))
    for all pairs (i, j) with i < j.
    """
    N = state.N
    G_val = float(G)
    result: Dict[Tuple[int, int], float] = {}
    prefactor = math.pi / (2.0 * math.sqrt(2.0))
    for i in range(N):
        q_i = np.asarray(state.q[i], dtype=np.float64)
        m_i = float(state.masses[i])
        for j in range(i + 1, N):
            q_j = np.asarray(state.q[j], dtype=np.float64)
            m_j = float(state.masses[j])
            r = float(np.linalg.norm(q_i - q_j))
            mu = G_val * (m_i + m_j)
            if r <= 1e-30:
                result[(i, j)] = 0.0
            elif mu > 1e-30:
                result[(i, j)] = prefactor * math.sqrt((r ** 3) / mu)
            else:
                result[(i, j)] = float("inf")
    return result


# Backward-compatible aliases
CollisionWarning = PairApproachPrediction
predict_closest_approach = predict_pair_approach

