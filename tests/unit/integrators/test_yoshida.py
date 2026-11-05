"""
Unit Tests for Layer 4: 4th-Order Symmetric Yoshida Integrator S_4.

Theoretical Alignments:
- Monograph: Chapter 5, Section 5.2.
- Formal Lean 4: `yoshida_order4_symmetric_composition`.
- Invariants: Microscopic time reversibility (< 1e-14), machine-precision
  momentum conservation (< 1e-14), and 4th-order convergence slope ~ 16.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.integrators.yoshida import Yoshida4Integrator


def test_yoshida_microscopic_time_reversibility():
    """Validates microscopic time reversibility Phi(-dt) o Phi(dt) = Id (< 1e-14)."""
    masses = np.array([1.0, 1.5, 2.0], dtype=np.float64)
    q0 = np.array([
        [0.8, -0.2, 0.1],
        [-0.4, 0.9, -0.2],
        [0.1, -0.5, 0.5],
    ], dtype=np.float64)
    p0 = np.array([
        [0.15, -0.1, 0.05],
        [-0.2, 0.15, -0.1],
        [0.05, -0.05, 0.05],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    integrator = Yoshida4Integrator(G=1.0)
    dt = 0.04

    s_fwd = integrator.step(state0, tree, dt=dt)
    s_back = integrator.step(s_fwd, tree, dt=-dt)

    assert np.allclose(s_back.q, state0.q, atol=1e-14), f"Reversibility error q: {np.max(np.abs(s_back.q - state0.q))}"
    assert np.allclose(s_back.p, state0.p, atol=1e-14), f"Reversibility error p: {np.max(np.abs(s_back.p - state0.p))}"


def test_yoshida_momentum_conservation():
    """Validates Lie-Bargmann linear and angular momentum conservation."""
    masses = np.array([1.2, 2.2, 1.8], dtype=np.float64)
    q0 = np.array([
        [1.0, 0.0, 0.0],
        [-0.5, 0.8, 0.0],
        [0.2, -0.3, 1.0],
    ], dtype=np.float64)
    p0 = np.array([
        [0.1, -0.2, 0.3],
        [-0.3, 0.1, -0.1],
        [0.2, 0.1, -0.2],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    inv0 = compute_invariants(state0)
    integrator = Yoshida4Integrator(G=1.0)

    curr_state = state0
    dt = 0.02
    for _ in range(5):
        curr_state = integrator.step(curr_state, tree, dt=dt)

    inv_final = compute_invariants(curr_state)
    assert np.allclose(inv_final.P, inv0.P, atol=1e-14)
    assert np.allclose(inv_final.L, inv0.L, atol=1e-13)
    assert np.allclose(inv_final.G, inv0.G, atol=1e-14)


def test_yoshida_fourth_order_convergence():
    """
    Validates 4th-order convergence slope:
        Error(2h) / Error(h) ~= 2^4 = 16 on a 3-body system over fixed interval T.
    """
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q0 = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 8.0, 0.5],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, -0.7, 0.0],
        [0.0, 0.7, 0.0],
        [0.25, 0.0, -0.1],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    integrator = Yoshida4Integrator(G=1.0)
    inv0 = compute_invariants(state0)

    T = 0.2
    # Step with 2h (dt = 0.02)
    s2 = state0
    for _ in range(int(round(T / 0.02))):
        s2 = integrator.step(s2, tree, dt=0.02)
    err2 = abs(compute_invariants(s2).H - inv0.H)

    # Step with h (dt = 0.01)
    s1 = state0
    for _ in range(int(round(T / 0.01))):
        s1 = integrator.step(s1, tree, dt=0.01)
    err1 = abs(compute_invariants(s1).H - inv0.H)

    ratio = err2 / max(err1, 1e-16)
    # Ratio should be close to 16.0 (4th order)
    assert 13.0 <= ratio <= 20.0, f"Expected 4th-order convergence ratio ~ 16, got {ratio}"


def test_yoshida_fsal_evaluations_count(monkeypatch):
    """Validates that Yoshida4Integrator executes exactly 4 force evaluations per step (FSAL fusion)."""
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q0 = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 8.0, 0.5],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, -0.7, 0.0],
        [0.0, 0.7, 0.0],
        [0.25, 0.0, -0.1],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    integrator = Yoshida4Integrator(G=1.0)
    assert integrator.evaluations_per_step == 4

    call_count = 0
    import normal_cone_nbody_float64.integrators.yoshida as yoshida_module
    orig_kick = yoshida_module.step_apex_kick

    def counting_kick(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        return orig_kick(*args, **kwargs)

    monkeypatch.setattr(yoshida_module, "step_apex_kick", counting_kick)

    s_next = integrator.step(state0, tree, dt=0.01)
    assert call_count == 4, f"Expected exactly 4 kick evaluations under FSAL, got {call_count}"


def test_yoshida_registers_current_integrator_stages(monkeypatch):
    """Validates that Yoshida4Integrator.step registers its 3 drift and 4 impulse stages in execution context."""
    from normal_cone_nbody_float64.integrators.common import get_current_integrator_stages
    from normal_cone_nbody_float64.integrators.yoshida import (
        YOSHIDA_DRIFT_WEIGHTS,
        YOSHIDA_IMPULSE_WEIGHTS,
    )
    import normal_cone_nbody_float64.integrators.yoshida as yoshida_module

    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q0 = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 8.0, 0.5]], dtype=np.float64)
    p0 = np.array([[0.0, -0.7, 0.0], [0.0, 0.7, 0.0], [0.25, 0.0, -0.1]], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    observed_stages = []
    orig_drift = yoshida_module.step_drift_and_kepler

    def inspecting_drift(*args, **kwargs):
        observed_stages.append(get_current_integrator_stages())
        return orig_drift(*args, **kwargs)

    monkeypatch.setattr(yoshida_module, "step_drift_and_kepler", inspecting_drift)

    assert get_current_integrator_stages() is None
    integrator = Yoshida4Integrator(G=1.0)
    integrator.step(state0, tree, dt=0.01)
    assert get_current_integrator_stages() is None
    assert len(observed_stages) == 3
    for stages in observed_stages:
        assert stages is not None, "Yoshida4Integrator.step failed to set current_integrator_stages"
        drifts, impulses = stages
        assert len(drifts) == 3
        assert len(impulses) == 4
        assert tuple(drifts) == YOSHIDA_DRIFT_WEIGHTS
        assert tuple(impulses) == YOSHIDA_IMPULSE_WEIGHTS


