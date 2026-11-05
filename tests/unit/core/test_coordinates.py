"""
Unit Tests for Layer 1 Coordinates (normal_cone_nbody_float64.core.coordinates).

Mathematical Verification Objectives:
- Symplectic isomorphism: (q, p) <-> (r_v, w_v) bijective round-trip error < 10^-15.
- Darboux 2-form invariance: sum <q_i, p_i> = <R_cm, P_tot> + sum <r_v, w_v> (omega0_skew, omega0At_nondeg).
- Kinetic energy diagonal decomposition: T_cart == T_cm + T_jacobi.
- Rees deformation limit: lim_{t -> 0} ev_t(q_parallel, n) = q_parallel (deformation_limit_zero).
"""

import pytest
import numpy as np

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.core.coordinates import (
    cartesian_to_jacobi,
    jacobi_to_cartesian,
    compute_jacobi_kinetic_energy,
    darboux_pairing_error,
    rees_deformation_retraction,
)
from normal_cone_nbody_float64.core.invariants import compute_kinetic_energy


@pytest.fixture
def random_system():
    rng = np.random.default_rng(999)
    N = 5
    q = rng.standard_normal((N, 3))
    p = rng.standard_normal((N, 3))
    masses = rng.uniform(0.5, 5.0, size=N)
    state = ExtendedPhaseState(q=q, p=p, masses=masses, t=1.0, E0=-3.0)
    tree = build_canonical_tree(masses, topology="comb")
    return state, tree


def test_jacobi_roundtrip_bijective_inversion(random_system):
    state, tree = random_system
    coords = cartesian_to_jacobi(state, tree)

    # Reconstruct Cartesian state
    reconstructed = jacobi_to_cartesian(
        coords=coords,
        tree=tree,
        masses=state.masses,
        t=state.t,
        E=state.E,
        E0=state.E0,
        center_of_mass=state.center_of_mass,
        total_momentum=state.total_momentum,
    )

    # Bijective precision to machine roundoff
    np.testing.assert_allclose(reconstructed.q, state.q, atol=1e-15)
    np.testing.assert_allclose(reconstructed.p, state.p, atol=1e-15)
    assert reconstructed.t == state.t
    assert reconstructed.E == state.E


def test_darboux_pairing_conservation(random_system):
    state, tree = random_system
    coords = cartesian_to_jacobi(state, tree)

    error = darboux_pairing_error(state, coords, tree)
    assert error < 1e-14


def test_kinetic_energy_diagonal_decoupling(random_system):
    state, tree = random_system
    coords = cartesian_to_jacobi(state, tree)

    t_cart = compute_kinetic_energy(state)
    t_jacobi = compute_jacobi_kinetic_energy(coords, tree)

    # Center-of-mass kinetic energy: (1/2) * ||P_tot||^2 / M_tot
    p_tot = state.total_momentum
    t_cm = 0.5 * float(np.dot(p_tot, p_tot)) / state.total_mass

    assert t_cart == pytest.approx(t_cm + t_jacobi, rel=1e-14)


def test_rees_deformation_retraction(random_system):
    state, tree = random_system
    coords = cartesian_to_jacobi(state, tree)

    # Apply scaling t = 0.1
    scaled_coords = rees_deformation_retraction(coords, tree, t_param=0.1)
    for node_id in coords.r:
        np.testing.assert_allclose(scaled_coords.r[node_id], 0.1 * coords.r[node_id], atol=1e-15)
        np.testing.assert_allclose(scaled_coords.w[node_id], coords.w[node_id] / 0.1, atol=1e-15)

    # At t -> 0, transverse distances retract onto the collision locus (deformation_limit_zero)
    zero_coords = rees_deformation_retraction(coords, tree, t_param=0.0)
    for node_id in zero_coords.r:
        np.testing.assert_allclose(zero_coords.r[node_id], np.zeros(3), atol=1e-15)
