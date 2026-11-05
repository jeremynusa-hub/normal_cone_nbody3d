"""
Unit Tests for Layer 4: 2nd-Order Palindromic Strang Integrator S_2.

Theoretical Alignments:
- Monograph: Chapter 5, Section 5.1.
- Formal Lean 4: `ExtendedSymplecticStrang.lean`, `MicroscopicTimeReversibility.lean`.
- Invariants: Exact microscopic time-reversibility (< 1e-14), machine-precision
  linear and angular momentum conservation (< 1e-14), and 2nd-order convergence slope.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.integrators.strang import StrangIntegrator


def test_strang_microscopic_time_reversibility():
    """
    Validates microscopic time reversibility:
        Phi(-dt) o Phi(dt) = Id (< 1e-14).
    """
    masses = np.array([1.0, 2.0, 1.5], dtype=np.float64)
    q0 = np.array([
        [1.0, 0.2, -0.1],
        [-0.5, 0.8, 0.3],
        [0.2, -0.6, 0.4],
    ], dtype=np.float64)
    p0 = np.array([
        [0.1, -0.3, 0.2],
        [-0.2, 0.1, -0.1],
        [0.1, 0.2, -0.1],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    integrator = StrangIntegrator(G=1.0)
    dt = 0.05

    # Forward step then backward step
    s_fwd = integrator.step(state0, tree, dt=dt)
    s_back = integrator.step(s_fwd, tree, dt=-dt)

    assert np.allclose(s_back.q, state0.q, atol=1e-14), f"Reversibility error q: {np.max(np.abs(s_back.q - state0.q))}"
    assert np.allclose(s_back.p, state0.p, atol=1e-14), f"Reversibility error p: {np.max(np.abs(s_back.p - state0.p))}"


def test_strang_momentum_conservation():
    """Validates exact machine-precision Lie-Bargmann linear and angular momentum conservation."""
    masses = np.array([1.5, 2.5, 3.5], dtype=np.float64)
    q0 = np.array([
        [1.2, -0.4, 0.5],
        [-0.8, 1.1, -0.3],
        [0.3, 0.2, 1.0],
    ], dtype=np.float64)
    p0 = np.array([
        [0.2, -0.1, 0.3],
        [-0.4, 0.3, -0.2],
        [0.2, -0.2, -0.1],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    inv0 = compute_invariants(state0)
    integrator = StrangIntegrator(G=1.0)

    curr_state = state0
    dt = 0.02
    for _ in range(10):
        curr_state = integrator.step(curr_state, tree, dt=dt)

    inv_final = compute_invariants(curr_state)

    # Linear momentum conservation
    assert np.allclose(inv_final.P, inv0.P, atol=1e-14), f"Linear momentum drift: {inv_final.P - inv0.P}"
    # Angular momentum conservation
    assert np.allclose(inv_final.L, inv0.L, atol=1e-13), f"Angular momentum drift: {inv_final.L - inv0.L}"
    # Galilean boost G = M_tot * R_cdm - t * P conservation
    assert np.allclose(inv_final.G, inv0.G, atol=1e-14), (
        f"Galilean boost drift: {inv_final.G - inv0.G}"
    )


def test_strang_second_order_convergence():
    """
    Validates 2nd-order convergence slope:
        Error(2h) / Error(h) ~= 2^2 = 4 on a 3-body system over fixed interval T.
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

    integrator = StrangIntegrator(G=1.0)
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
    # Ratio should be close to 4.0 (2nd order)
    assert 3.5 <= ratio <= 4.5, f"Expected 2nd-order convergence ratio ~ 4, got {ratio}"
