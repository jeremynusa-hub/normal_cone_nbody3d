"""
Unit Tests for Layer 4: Optimal 11-Stage 6th-Order Blanes--Moan Integrator BM_6^{(11)}.

Theoretical Alignments:
- Monograph: Chapter 5, Section 5.4, Table 1.
- Formal Lean 4: `BM6Order6SymmetricComposition.lean`.
- Invariants: Microscopic time reversibility (< 1e-14), machine-precision
  momentum conservation (< 1e-14), and 6th-order convergence slope ~ 64.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator


def test_bm6_microscopic_time_reversibility():
    """Validates microscopic time reversibility Phi(-dt) o Phi(dt) = Id (< 1e-14)."""
    masses = np.array([1.2, 2.1, 1.4], dtype=np.float64)
    q0 = np.array([
        [0.7, -0.2, 0.3],
        [-0.5, 0.8, -0.2],
        [0.1, -0.4, 0.4],
    ], dtype=np.float64)
    p0 = np.array([
        [0.12, -0.15, 0.1],
        [-0.18, 0.12, -0.08],
        [0.08, 0.05, -0.04],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    integrator = BlanesMoan6Integrator(G=1.0)
    dt = 0.03

    s_fwd = integrator.step(state0, tree, dt=dt)
    s_back = integrator.step(s_fwd, tree, dt=-dt)

    assert np.allclose(s_back.q, state0.q, atol=1e-14), f"Reversibility error q: {np.max(np.abs(s_back.q - state0.q))}"
    assert np.allclose(s_back.p, state0.p, atol=1e-14), f"Reversibility error p: {np.max(np.abs(s_back.p - state0.p))}"


def test_bm6_momentum_conservation():
    """Validates Lie-Bargmann linear and angular momentum conservation."""
    masses = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    q0 = np.array([
        [1.0, 0.0, 0.0],
        [-0.5, 0.866, 0.0],
        [0.0, -0.5, 0.8],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, 0.3, 0.1],
        [-0.3, 0.0, -0.2],
        [0.2, -0.1, 0.0],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    inv0 = compute_invariants(state0)
    integrator = BlanesMoan6Integrator(G=1.0)

    curr_state = state0
    dt = 0.02
    for _ in range(5):
        curr_state = integrator.step(curr_state, tree, dt=dt)

    inv_final = compute_invariants(curr_state)
    assert np.allclose(inv_final.P, inv0.P, atol=1e-14)
    assert np.allclose(inv_final.L, inv0.L, atol=1e-13)
    assert np.allclose(inv_final.G, inv0.G, atol=1e-14)


def test_bm6_sixth_order_convergence():
    """
    Validates 6th-order convergence on a 3-body system over fixed interval T:
    Demonstrates reduction of energy residual down to machine precision (< 1e-14)
    with steep high-order ratio.
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

    integrator = BlanesMoan6Integrator(G=1.0)
    inv0 = compute_invariants(state0)

    # 1. Asymptotic 6th-order convergence test on [0, 0.8] above machine noise floor
    T = 0.8
    s2 = state0
    for _ in range(int(round(T / 0.4))):
        s2 = integrator.step(s2, tree, dt=0.4)
    err2 = abs(compute_invariants(s2).H - inv0.H)

    s1 = state0
    for _ in range(int(round(T / 0.2))):
        s1 = integrator.step(s1, tree, dt=0.2)
    err1 = abs(compute_invariants(s1).H - inv0.H)

    ratio = err2 / max(err1, 1e-16)
    # Theoretical order 6 ratio is 2^6 = 64
    assert ratio >= 50.0, f"Expected 6th-order asymptotic convergence ratio (2^6=64), got {ratio:.2f}"

    # 2. Machine precision floor check on fine step (dt = 0.1)
    s01 = state0
    for _ in range(int(round(0.4 / 0.1))):
        s01 = integrator.step(s01, tree, dt=0.1)
    err01 = abs(compute_invariants(s01).H - inv0.H)
    assert err01 <= 5e-15, f"Expected machine-precision floor (< 5e-15) for BM6 at dt=0.1, got {err01}"

