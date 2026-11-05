"""
Unit Tests for Extended Precision (Float128 / Mpmath 35 dps) System History Tracker.

Theoretical Alignments:
- Monograph: Chapter 2 Section 2.3 & Chapter 5 Section 5.1.
- Validates that HistoryTracker preserves 35-decimal mpmath precision on times, positions,
  momenta, and Lie--Bargmann invariants without intermediate float64 truncation.
"""

import numpy as np
import pytest
import mpmath as mp

from normal_cone_nbody_float128.core.state import ExtendedPhaseState
from normal_cone_nbody_float128.core.invariants import compute_invariants
from normal_cone_nbody_float128.system.history import HistoryTracker


def test_history_recording_float128_shapes_and_precision():
    """
    Validates recording snapshots with ExtendedPhaseState(is_f128=True) and retrieving
    exact 35-digit precision without float64 truncation.
    """
    tracker = HistoryTracker(record_interval=1)
    masses = np.array([mp.mpf("1.0"), mp.mpf("1.5"), mp.mpf("0.5")], dtype=object)

    for step in range(5):
        # Precise 35-digit time and coordinates
        t_step = mp.mpf(step) * mp.mpf("0.050000000000000000000000000000000001")
        q = np.array([
            [mp.mpf("1.12345678901234567890123456789012345") + step, mp.mpf("0.0"), mp.mpf("0.0")],
            [mp.mpf("-0.5"), mp.mpf("0.0"), mp.mpf("0.0")],
            [mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")],
        ], dtype=object)
        p = np.array([
            [mp.mpf("0.0"), mp.mpf("0.1"), mp.mpf("0.0")],
            [mp.mpf("0.0"), mp.mpf("-0.1"), mp.mpf("0.0")],
            [mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")],
        ], dtype=object)
        state = ExtendedPhaseState(q=q, p=p, masses=masses, t=t_step)
        tracker.record(state)

    assert len(tracker) == 5
    assert tracker.is_f128 is True
    assert tracker.times.shape == (5,)
    assert tracker.times.dtype == object

    # Verify that the 35th decimal is preserved exactly
    t_first = tracker.times[1]
    assert abs(t_first - mp.mpf("0.050000000000000000000000000000000001")) < mp.mpf("1e-35")

    # Verify positions and momenta retain object dtype
    assert tracker.positions.shape == (5, 3, 3)
    assert tracker.positions.dtype == object
    q00 = tracker.positions[0, 0, 0]
    assert abs(q00 - mp.mpf("1.12345678901234567890123456789012345")) < mp.mpf("1e-35")

    assert tracker.momenta.shape == (5, 3, 3)
    assert tracker.momenta.dtype == object

    assert tracker.hamiltonians.shape == (5,)
    assert tracker.hamiltonians.dtype == object

    # Verify float64 accessor conversion works smoothly
    assert tracker.times_f64.dtype == np.float64
    assert tracker.times_f64.shape == (5,)
    assert tracker.positions_f64.dtype == np.float64
    assert tracker.positions_f64.shape == (5, 3, 3)
    assert tracker.momenta_f64.dtype == np.float64
    assert tracker.hamiltonians_f64.dtype == np.float64


def test_history_drift_metrics_float128():
    """
    Validates max_energy_drift, max_momentum_drift, max_angular_momentum_drift
    calculating infinitesimal variations (~ 1e-30) without float64 underflow.
    """
    tracker = HistoryTracker(record_interval=1)
    masses = np.array([mp.mpf("1.0"), mp.mpf("1.0")], dtype=object)

    # Initial state
    q0 = np.array([[mp.mpf("1.0"), mp.mpf("0.0"), mp.mpf("0.0")],
                   [mp.mpf("-1.0"), mp.mpf("0.0"), mp.mpf("0.0")]], dtype=object)
    p0 = np.array([[mp.mpf("0.0"), mp.mpf("0.5"), mp.mpf("0.0")],
                   [mp.mpf("0.0"), mp.mpf("-0.5"), mp.mpf("0.0")]], dtype=object)
    s0 = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=mp.mpf("0.0"))
    tracker.record(s0)

    # Micro-perturbed state: displacement of 1e-30 (which would be 0 in float64!)
    perturbation = mp.mpf("1e-30")
    q1 = np.array([[mp.mpf("1.0") + perturbation, mp.mpf("0.0"), mp.mpf("0.0")],
                   [mp.mpf("-1.0"), mp.mpf("0.0"), mp.mpf("0.0")]], dtype=object)
    p1 = np.array([[mp.mpf("0.0"), mp.mpf("0.5"), mp.mpf("0.0")],
                   [mp.mpf("0.0"), mp.mpf("-0.5") + perturbation, mp.mpf("0.0")]], dtype=object)
    s1 = ExtendedPhaseState(q=q1, p=p1, masses=masses, t=mp.mpf("0.01"))
    tracker.record(s1)

    delta_H = tracker.max_energy_drift()
    delta_P = tracker.max_momentum_drift()

    # In float64, perturbation 1e-30 would vanish into 0.0
    # In float128 (mpmath 35 dps), it must be strictly positive and around 1e-30!
    assert delta_H > mp.mpf("0.0")
    assert delta_H < mp.mpf("1e-28")
    assert delta_P > mp.mpf("0.0")
    assert delta_P < mp.mpf("1e-28")

    metrics = tracker.summary_metrics()
    assert metrics["total_checkpoints"] == 2
    assert abs(metrics["delta_H_max"] - delta_H) < mp.mpf("1e-35")


def test_history_cadence_and_max_records_float128():
    """
    Validates cadence subsampling and circular capping under float128 mode.
    """
    tracker = HistoryTracker(record_interval=2, max_records=3)
    masses = np.array([mp.mpf("1.0"), mp.mpf("1.0")], dtype=object)
    q = np.array([[mp.mpf("1.0"), mp.mpf("0.0"), mp.mpf("0.0")],
                  [mp.mpf("-1.0"), mp.mpf("0.0"), mp.mpf("0.0")]], dtype=object)
    p = np.array([[mp.mpf("0.0"), mp.mpf("0.5"), mp.mpf("0.0")],
                  [mp.mpf("0.0"), mp.mpf("-0.5"), mp.mpf("0.0")]], dtype=object)

    for step in range(10):
        s = ExtendedPhaseState(q=q, p=p, masses=masses, t=mp.mpf(step) * mp.mpf("0.1"))
        tracker.record(s)

    # 10 steps with interval 2 -> recorded steps 0, 2, 4, 6, 8 (5 steps), capped at max_records=3
    assert len(tracker) == 3
    # The 3 remaining should be steps 4, 6, 8
    assert abs(tracker.times[0] - mp.mpf("0.4")) < mp.mpf("1e-32")
    assert abs(tracker.times[-1] - mp.mpf("0.8")) < mp.mpf("1e-32")
