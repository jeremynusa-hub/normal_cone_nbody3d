"""
Unit Tests for Layer 5: Trajectory and Invariant History Tracker (history.py).

Theoretical Alignments:
- Monograph: Chapter 2 Section 2.3 & Chapter 5 Section 5.1.
- Invariants: Vectorized property dimensions, drift estimators, and cadence sampling.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.system.history import HistoryTracker


def test_history_recording_and_vectorized_shapes():
    """
    Validates recording snapshots and retrieving vectorized NumPy arrays:
    times (S,), positions (S, N, 3), momenta (S, N, 3), and hamiltonians (S,).
    """
    tracker = HistoryTracker(record_interval=1)
    masses = np.array([1.0, 1.5, 0.5], dtype=np.float64)

    for step in range(5):
        q = np.ones((3, 3), dtype=np.float64) * (step + 1.0)
        p = np.ones((3, 3), dtype=np.float64) * 0.1
        state = ExtendedPhaseState(q=q, p=p, masses=masses, t=step * 0.05)
        tracker.record(state)

    assert len(tracker) == 5
    assert tracker.times.shape == (5,)
    assert np.allclose(tracker.times, [0.0, 0.05, 0.10, 0.15, 0.20])
    assert tracker.positions.shape == (5, 3, 3)
    assert tracker.momenta.shape == (5, 3, 3)
    assert tracker.hamiltonians.shape == (5,)
    assert tracker.linear_momenta.shape == (5, 3)
    assert tracker.angular_momenta.shape == (5, 3)
    assert tracker.boosts.shape == (5, 3)


def test_history_drift_metrics():
    """
    Validates max_energy_drift, max_momentum_drift, max_angular_momentum_drift calculations.
    """
    tracker = HistoryTracker(record_interval=1)
    masses = np.array([1.0, 1.0], dtype=np.float64)

    # Initial state
    q0 = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], dtype=np.float64)
    p0 = np.array([[0.0, 0.5, 0.0], [0.0, -0.5, 0.0]], dtype=np.float64)
    s0 = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=0.0)
    tracker.record(s0)

    # Slightly perturbed state
    q1 = np.array([[1.001, 0.0, 0.0], [-1.0, 0.0, 0.0]], dtype=np.float64)
    p1 = np.array([[0.0, 0.5, 0.0], [0.0, -0.499, 0.0]], dtype=np.float64)
    s1 = ExtendedPhaseState(q=q1, p=p1, masses=masses, t=0.01)
    tracker.record(s1)

    metrics = tracker.summary_metrics()
    assert metrics["total_checkpoints"] == 2.0
    assert metrics["delta_H_max"] > 0.0
    assert metrics["delta_P_max"] > 0.0
    assert tracker.max_energy_drift() == metrics["delta_H_max"]


def test_history_cadence_and_max_records():
    """
    Validates record_interval subsampling and max_records circular capping.
    """
    # Record only every 3 steps, max 4 records
    tracker = HistoryTracker(record_interval=3, max_records=4)
    masses = np.array([1.0, 1.0], dtype=np.float64)
    q = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], dtype=np.float64)
    p = np.array([[0.0, 0.5, 0.0], [0.0, -0.5, 0.0]], dtype=np.float64)

    recorded_count = 0
    for step in range(12):
        s = ExtendedPhaseState(q=q, p=p, masses=masses, t=step * 0.1)
        if tracker.record(s):
            recorded_count += 1

    # Steps 0, 3, 6, 9 -> 4 records
    assert len(tracker) == 4
    # Earliest record should be step 0
    assert tracker.times[0] == 0.0

    # Add 2 more recorded steps (steps 12, 15) -> total capped at 4
    tracker.record(ExtendedPhaseState(q=q, p=p, masses=masses, t=1.2), force=True)
    tracker.record(ExtendedPhaseState(q=q, p=p, masses=masses, t=1.5), force=True)
    assert len(tracker) == 4
    assert tracker.times[-1] == 1.5
