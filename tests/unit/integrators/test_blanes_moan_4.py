"""
Unit Tests for Layer 4: Optimal 7-Stage 4th-Order Blanes--Moan Integrator BM_4^{(7)}.

Theoretical Alignments:
- Monograph: Chapter 5, Section 5.3.
- Formal Lean 4: `01_BlanesMoanShadowConfinement_Theorems.lean`.
- Invariants: Microscopic time reversibility (< 1e-14), machine-precision
  momentum conservation (< 1e-14), 4th-order convergence slope, and spectral error superiority over Yoshida.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float64.integrators.yoshida import Yoshida4Integrator


def test_bm4_microscopic_time_reversibility():
    """Validates microscopic time reversibility Phi(-dt) o Phi(dt) = Id (< 1e-14)."""
    masses = np.array([1.0, 1.2, 0.8], dtype=np.float64)
    q0 = np.array([
        [0.5, -0.3, 0.2],
        [-0.4, 0.6, -0.1],
        [0.2, -0.4, 0.3],
    ], dtype=np.float64)
    p0 = np.array([
        [0.1, -0.1, 0.2],
        [-0.2, 0.1, -0.1],
        [0.1, 0.05, -0.05],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    integrator = BlanesMoan4Integrator(G=1.0)
    dt = 0.03

    s_fwd = integrator.step(state0, tree, dt=dt)
    s_back = integrator.step(s_fwd, tree, dt=-dt)

    assert np.allclose(s_back.q, state0.q, atol=1e-15), f"Reversibility error q: {np.max(np.abs(s_back.q - state0.q))}"
    assert np.allclose(s_back.p, state0.p, atol=1e-15), f"Reversibility error p: {np.max(np.abs(s_back.p - state0.p))}"

    inv0 = compute_invariants(state0)
    inv_fwd = compute_invariants(s_fwd)
    assert np.allclose(inv_fwd.P, inv0.P, atol=1e-15)
    assert np.allclose(inv_fwd.L, inv0.L, atol=1e-14)
    assert np.allclose(inv_fwd.G, inv0.G, atol=1e-15)


def test_bm4_fourth_order_convergence():
    """Validates 4th-order convergence slope Error(2h)/Error(h) ~ 16 on a 3-body system over fixed interval T."""
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

    integrator = BlanesMoan4Integrator(G=1.0)
    inv0 = compute_invariants(state0)

    T = 0.2
    # Step with 2h (dt = 0.04)
    s2 = state0
    for _ in range(int(round(T / 0.04))):
        s2 = integrator.step(s2, tree, dt=0.04)
    err2 = abs(compute_invariants(s2).H - inv0.H)

    # Step with h (dt = 0.02)
    s1 = state0
    for _ in range(int(round(T / 0.02))):
        s1 = integrator.step(s1, tree, dt=0.02)
    err1 = abs(compute_invariants(s1).H - inv0.H)

    ratio = err2 / max(err1, 1e-16)
    assert 12.0 <= ratio <= 22.0, f"Expected 4th-order convergence ratio ~ 16, got {ratio}"


def test_bm4_spectral_superiority_over_yoshida():
    """
    Validates Theorem 5.2: The residual error constant of BM4 is significantly smaller
    than Yoshida4 (|C_5^{BM}| << |C_5^{Yoshida}|).
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

    bm4 = BlanesMoan4Integrator(G=1.0)
    y4 = Yoshida4Integrator(G=1.0)
    inv0 = compute_invariants(state0)

    T = 0.2
    dt = 0.02
    s_bm4 = state0
    for _ in range(int(round(T / dt))):
        s_bm4 = bm4.step(s_bm4, tree, dt=dt)

    s_y4 = state0
    for _ in range(int(round(T / dt))):
        s_y4 = y4.step(s_y4, tree, dt=dt)

    err_bm4 = abs(compute_invariants(s_bm4).H - inv0.H)
    err_y4 = abs(compute_invariants(s_y4).H - inv0.H)

    # BM4 error must be significantly smaller than Yoshida4 error (factor > 5x)
    assert err_bm4 < err_y4 * 0.2, f"BM4 error ({err_bm4}) not significantly smaller than Yoshida ({err_y4})"


def test_bm4_fsal_kick_fusion_numerical_equivalence():
    """
    Validates Chapter 5 Section 5.3 & Iteration 146 Phase 3 requirement:
    1. BM4 uses FSAL fused kicks: exactly 8 impulse weights and 7 drift weights.
    2. Sum of drift weights == 1.0 and sum of impulse weights == 1.0 (< 1e-15).
    3. The fused scheme produces numerical results identical to the 7-stage Strang composition
       at machine precision (< 1e-15).
    4. Evaluates exactly 8 apex kicks per isolated step instead of 14.
    """
    from normal_cone_nbody_float64.integrators.blanes_moan_4 import (
        BM4_DRIFT_WEIGHTS,
        BM4_IMPULSE_WEIGHTS,
        BM4_WEIGHTS,
    )
    from normal_cone_nbody_float64.integrators.strang import StrangIntegrator

    bm4 = BlanesMoan4Integrator(G=1.0)
    assert len(bm4.drift_sequence) == 7
    assert len(bm4.impulse_sequence) == 8
    assert np.isclose(sum(bm4.drift_sequence), 1.0, atol=1e-15)
    assert np.isclose(sum(bm4.impulse_sequence), 1.0, atol=1e-15)

    masses = np.array([1.0, 1.2, 0.8], dtype=np.float64)
    q0 = np.array([
        [0.5, -0.3, 0.2],
        [-0.4, 0.6, -0.1],
        [0.2, -0.4, 0.3],
    ], dtype=np.float64)
    p0 = np.array([
        [0.1, -0.1, 0.2],
        [-0.2, 0.1, -0.1],
        [0.1, 0.05, -0.05],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    dt = 0.025

    # 1. Step with fused BM4
    state_fused = bm4.step(state0, tree, dt=dt)

    # 2. Step with explicit Strang composition (14 kicks)
    strang = StrangIntegrator(G=1.0, form="kick_exterior")
    state_strang = state0
    for w in BM4_WEIGHTS:
        state_strang = strang.step(state_strang, tree, dt=w * dt)

    # Exact mathematical equivalence
    assert np.allclose(state_fused.q, state_strang.q, atol=1e-15), (
        f"Position mismatch: {np.max(np.abs(state_fused.q - state_strang.q))}"
    )
    assert np.allclose(state_fused.p, state_strang.p, atol=1e-15), (
        f"Momentum mismatch: {np.max(np.abs(state_fused.p - state_strang.p))}"
    )

    # Count calls to apex forces during step
    call_count = 0
    original_compute = bm4._default_apex.compute_apex_forces

    def counting_compute(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        return original_compute(*args, **kwargs)

    bm4._default_apex.compute_apex_forces = counting_compute
    try:
        bm4.step(state0, tree, dt=dt)
        # Exactly 8 calls for 8 fused impulse stages (down from 14)
        assert call_count == 8, f"Expected 8 apex force calls, got {call_count}"
    finally:
        bm4._default_apex.compute_apex_forces = original_compute

