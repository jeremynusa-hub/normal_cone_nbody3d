"""
Unit Tests for Layer 6: Collision Radar and Exceptional Divisor Proximity Anticipation (collision_radar.py).

Theoretical Alignments:
- Monograph: Chapter 6 (Exceptional Collision Divisors E_I in partial Q[N]),
  and Section 4.3 (Intrinsic Compactness and Marchal-Sundman Threshold).
- Invariants: Exact ballistic impact time on head-on collision (< 1e-14),
  receding pair non-alert (is_approaching = False), and hazard score bounding [0, 1].
"""

import math
import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.diagnostics.collision_radar import (
    predict_pair_approach,
    CollisionRadar,
    PairApproachPrediction,
    CollisionRadarReport,
)


def test_head_on_ballistic_collision_prediction():
    """
    Validates exact ballistic time-to-approach calculation on a 1D head-on collision.
    Two unit-mass particles at separation r=2.0 moving toward each other at |v_rel| = 2.0:
    Impact occurs at exactly t_approach = r / |v_r| = 2.0 / 2.0 = 1.0.
    """
    q0 = np.array([-1.0, 0.0, 0.0], dtype=np.float64)
    q1 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    v0 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    v1 = np.array([-1.0, 0.0, 0.0], dtype=np.float64)

    pred = predict_pair_approach(
        q_i=q0,
        q_j=q1,
        v_i=v0,
        v_j=v1,
        m_i=1.0,
        m_j=1.0,
        body_i=0,
        body_j=1,
        G=1.0,
    )

    assert pred.r_current == pytest.approx(2.0, abs=1e-15)
    assert pred.radial_velocity == pytest.approx(-2.0, abs=1e-15)
    assert pred.is_approaching is True
    # Impact time must be machine-exact
    assert pred.time_to_approach == pytest.approx(1.0, abs=1e-14)
    # Head-on collision has zero angular momentum => r_peri == 0
    assert pred.predicted_r_peri == 0.0


def test_receding_pair_no_approach_alert():
    """
    Validates that a pair moving away from each other has is_approaching=False
    and infinite time-to-approach.
    """
    q0 = np.array([-1.0, 0.0, 0.0], dtype=np.float64)
    q1 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    v0 = np.array([-1.0, 0.0, 0.0], dtype=np.float64)
    v1 = np.array([1.0, 0.0, 0.0], dtype=np.float64)

    pred = predict_pair_approach(
        q_i=q0,
        q_j=q1,
        v_i=v0,
        v_j=v1,
        m_i=1.0,
        m_j=1.0,
        body_i=0,
        body_j=1,
        G=1.0,
    )

    assert pred.radial_velocity == pytest.approx(2.0, abs=1e-15)
    assert pred.is_approaching is False
    assert math.isinf(pred.time_to_approach)


def test_circular_orbit_osculating_periastron():
    """
    Validates circular two-body orbit:
    r = 1.0, v_circ = sqrt(G*M/r) = sqrt(2), radial velocity v_r = 0,
    osculating eccentricity e ~ 0, r_peri ~ 1.0.
    """
    q0 = np.array([-0.5, 0.0, 0.0], dtype=np.float64)
    q1 = np.array([0.5, 0.0, 0.0], dtype=np.float64)
    # Total mass M = 2.0, G = 1.0 => v_circ = sqrt(2.0 / 1.0) = sqrt(2)
    # Relative velocity: v_rel = (0, sqrt(2), 0)
    v0 = np.array([0.0, -math.sqrt(2.0) / 2.0, 0.0], dtype=np.float64)
    v1 = np.array([0.0, math.sqrt(2.0) / 2.0, 0.0], dtype=np.float64)

    pred = predict_pair_approach(
        q_i=q0,
        q_j=q1,
        v_i=v0,
        v_j=v1,
        m_i=1.0,
        m_j=1.0,
        body_i=0,
        body_j=1,
        G=1.0,
    )

    assert pred.r_current == pytest.approx(1.0, abs=1e-15)
    assert abs(pred.radial_velocity) < 1e-15
    assert pred.osculating_eccentricity == pytest.approx(0.0, abs=1e-12)
    assert pred.predicted_r_peri == pytest.approx(1.0, abs=1e-12)


def test_collision_radar_scan_and_hazard_score():
    """
    Validates full 3-body system radar scanning and normalized hazard scoring.
    """
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    # Tight approaching pair (0, 1) with separation 0.05 < hazard_ref = 0.1
    q = np.array([
        [-0.025, 0.0, 0.0],
        [0.025, 0.0, 0.0],
        [0.0, 5.0, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.5, 0.0, 0.0],
        [-0.5, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    radar = CollisionRadar(G=1.0, hazard_distance_threshold=0.1, hazard_time_threshold=1.0)
    report = radar.scan(state)

    assert report.closest_pair == (0, 1)
    assert report.min_separation == pytest.approx(0.05, abs=1e-15)
    assert report.imminent_collision is True
    assert 0.0 <= report.hazard_score <= 1.0
    assert report.hazard_score > 0.5  # High hazard due to proximity and approach

    # Test serialization
    d = report.to_dict()
    assert d["closest_pair"] == (0, 1)
    assert d["imminent_collision"] is True


def test_radar_module_compatibility_exports_and_free_fall_times():
    """
    Validates that diagnostics.radar imports cleanly and exposes CollisionRadar,
    CollisionRadarReport, PairApproachPrediction, CollisionWarning,
    predict_closest_approach, and compute_pairwise_free_fall_times (D-146-04).
    """
    from normal_cone_nbody_float64.diagnostics.radar import (
        CollisionRadar as RadarAlias,
        CollisionRadarReport as ReportAlias,
        PairApproachPrediction as PredAlias,
        CollisionWarning,
        predict_closest_approach,
        compute_pairwise_free_fall_times,
    )

    assert RadarAlias is CollisionRadar
    assert ReportAlias is CollisionRadarReport
    assert PredAlias is PairApproachPrediction
    assert CollisionWarning is PairApproachPrediction
    assert predict_closest_approach is predict_pair_approach

    masses = np.array([1.0, 1.0], dtype=np.float64)
    q = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    ff_times = compute_pairwise_free_fall_times(state, G=1.0)
    expected_t_ff = (math.pi / (2.0 * math.sqrt(2.0))) * math.sqrt((2.0 ** 3) / 2.0)
    assert (0, 1) in ff_times
    assert ff_times[(0, 1)] == pytest.approx(expected_t_ff, rel=1e-14)

