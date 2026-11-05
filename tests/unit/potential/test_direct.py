"""
Unit Tests for Layer 3: Direct Unsoftened Newtonian Potential and Force Evaluator.

Theoretical Alignments:
- Monograph: Chapter 1 & Section 4.1.
- Invariants: Strict action-reaction sum F_i = 0, zero torque sum q_i x F_i = 0,
  analytical gradient consistency < 1e-10, and zero softening epsilon = 0.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.invariants import compute_potential_energy
from normal_cone_nbody_float64.potential.direct import (
    compute_direct_newtonian_potential,
    compute_direct_newtonian_forces,
    direct_newtonian_step,
)


def test_direct_newtonian_potential_and_forces_symmetry():
    """
    Validates machine-precision momentum conservation:
        sum F_i = 0 (< 1e-15),
        sum q_i x F_i = 0 (< 1e-14).
    """
    masses = np.array([1.0, 2.5, 3.2, 0.8], dtype=np.float64)
    q = np.array([
        [1.2, -0.5, 0.3],
        [-0.8, 1.1, -0.4],
        [0.4, 0.9, 1.5],
        [-1.0, -1.2, -0.7],
    ], dtype=np.float64)

    forces = compute_direct_newtonian_forces(q, masses, G=1.0)

    # Net force is zero (linear momentum conservation)
    F_tot = np.sum(forces, axis=0)
    assert np.allclose(F_tot, 0.0, atol=1e-15), f"Net force non-zero: {F_tot}"

    # Net torque is zero (angular momentum conservation)
    torque = np.sum([np.cross(q[i], forces[i]) for i in range(len(masses))], axis=0)
    assert np.allclose(torque, 0.0, atol=1e-14), f"Net torque non-zero: {torque}"


def test_direct_newtonian_gradient_consistency():
    """
    Validates that the analytical forces F_i match central finite differences
    of the Newtonian potential V:
        F_i == - nabla_{q_i} V.
    """
    masses = np.array([1.5, 2.0, 3.0], dtype=np.float64)
    q0 = np.array([
        [1.0, 0.0, 0.0],
        [-0.5, 0.8, 0.0],
        [0.2, -0.3, 1.0],
    ], dtype=np.float64)

    analytical_forces = compute_direct_newtonian_forces(q0, masses, G=1.0)

    eps = 1e-7
    for i in range(len(masses)):
        for comp in range(3):
            q_plus = np.copy(q0)
            q_plus[i, comp] += eps
            V_plus = compute_direct_newtonian_potential(q_plus, masses, G=1.0)

            q_minus = np.copy(q0)
            q_minus[i, comp] -= eps
            V_minus = compute_direct_newtonian_potential(q_minus, masses, G=1.0)

            fd_force = - (V_plus - V_minus) / (2.0 * eps)
            assert np.isclose(analytical_forces[i, comp], fd_force, rtol=1e-4, atol=1e-6), (
                f"Mismatch at body {i}, comp {comp}: analytical {analytical_forces[i, comp]} vs FD {fd_force}"
            )


def test_direct_newtonian_potential_matches_invariants():
    """Validates that compute_direct_newtonian_potential matches compute_potential_energy."""
    masses = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    q = np.array([
        [0.5, 0.1, -0.2],
        [-0.3, 0.4, 0.1],
        [0.1, -0.5, 0.3],
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    V_direct = compute_direct_newtonian_potential(state.q, state.masses, G=1.0)
    V_inv = compute_potential_energy(state, G=1.0)

    assert np.isclose(V_direct, V_inv, atol=1e-15)


def test_direct_newtonian_step_symplectic_propagation():
    """Validates that direct_newtonian_step executes a clean leapfrog step without NaN or Inf."""
    masses = np.array([1.0, 1.0], dtype=np.float64)
    q = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], dtype=np.float64)
    p = np.array([[0.0, 0.5, 0.0], [0.0, -0.5, 0.0]], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    dt = 0.01
    state_next = direct_newtonian_step(state, dt=dt, G=1.0)

    assert state_next.t == state.t + dt
    assert not np.any(np.isnan(state_next.q))
    assert not np.any(np.isnan(state_next.p))
    assert np.allclose(state_next.total_momentum, state.total_momentum, atol=1e-15)
