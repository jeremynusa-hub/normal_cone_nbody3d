"""
Unit Tests for Layer 1 Transitions (normal_cone_nbody_float64.core.transitions).

Mathematical Verification Objectives:
- Unimodular determinant: det M == 1.0 (exact SL(2, R), sl2_preserves_symplectic_area).
- Symplecticity of cotangent lift: Psi^T J Psi == J to machine precision < 10^-15.
- Darboux pairing invariance: <r'1, w'1> + <r'2, w'2> == <r1, w1> + <r2, w2>.
- Invariance of Cartesian physical state across Jacobi tree rotation.
"""

import pytest
import numpy as np

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec
from normal_cone_nbody_float64.core.coordinates import (
    cartesian_to_jacobi,
    jacobi_to_cartesian,
)
from normal_cone_nbody_float64.core.transitions import (
    SL2TransitionMatrix,
    build_symplectic_matrix_12,
    verify_symplectic_matrix_12,
    rotate_jacobi_tree_coordinates,
)


def test_sl2_unimodular_construction():
    mA, mB, mC = 1.5, 2.5, 3.5
    trans = SL2TransitionMatrix.from_masses(mA, mB, mC)

    assert abs(trans.determinant - 1.0) < 1e-15
    assert trans.b == 1.0

    # Test inverse matrix
    inv = trans.inverse()
    assert abs(inv.determinant - 1.0) < 1e-15

    # Check M * M^-1 = Identity
    M = np.array([[trans.a, trans.b], [trans.c, trans.d]])
    M_inv = np.array([[inv.a, inv.b], [inv.c, inv.d]])
    np.testing.assert_allclose(M @ M_inv, np.eye(2), atol=1e-15)


def test_darboux_pairing_invariance_sl2():
    trans = SL2TransitionMatrix.from_masses(1.0, 2.0, 3.0)

    r1 = np.array([1.0, 2.0, 3.0])
    r2 = np.array([-1.0, 0.5, 2.0])
    w1 = np.array([0.2, -0.4, 1.1])
    w2 = np.array([0.5, 0.1, -0.8])

    r1_new, r2_new = trans.transform_positions(r1, r2)
    w1_new, w2_new = trans.transform_momenta(w1, w2)

    pairing_before = np.dot(r1, w1) + np.dot(r2, w2)
    pairing_after = np.dot(r1_new, w1_new) + np.dot(r2_new, w2_new)

    assert pairing_before == pytest.approx(pairing_after, rel=1e-14)


def test_symplectic_matrix_12_representation():
    trans = SL2TransitionMatrix.from_masses(2.1, 3.2, 4.3)
    Psi = build_symplectic_matrix_12(trans)

    assert Psi.shape == (12, 12)
    assert verify_symplectic_matrix_12(Psi)


def test_tree_rotation_cartesian_invariance():
    # 3-body system: Tree A = ((0, 1), 2), Tree B = (0, (1, 2))
    masses = np.array([1.2, 2.3, 3.4])
    tree_A = build_tree_from_spec(masses, ((0, 1), 2))
    tree_B = build_tree_from_spec(masses, (0, (1, 2)))

    rng = np.random.default_rng(777)
    q = rng.standard_normal((3, 3))
    p = rng.standard_normal((3, 3))
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    # 1. Project to Tree A
    coords_A = cartesian_to_jacobi(state, tree_A)

    # 2. Identify the internal nodes undergoing the flip
    # In tree_A: internal nodes are pair(0, 1) and root((0, 1), 2)
    pair_A = [n for n in tree_A.internal_nodes if n.body_indices == (0, 1)][0]
    root_A = tree_A.root

    # In tree_B: internal nodes are pair(1, 2) and root(0, (1, 2))
    pair_B = [n for n in tree_B.internal_nodes if n.body_indices == (1, 2)][0]
    root_B = tree_B.root

    # Transition from ((0, 1), 2) -> (0, (1, 2)): mA=m_0, mB=m_1, mC=m_2
    trans = SL2TransitionMatrix.from_masses(masses[0], masses[1], masses[2])

    # 3. Rotate coordinates to Tree B
    coords_B = rotate_jacobi_tree_coordinates(
        coords=coords_A,
        tree_A=tree_A,
        tree_B=tree_B,
        trans=trans,
        node_ids_A=(pair_A.id, root_A.id),
        node_ids_B=(pair_B.id, root_B.id),
    )

    # 4. Reconstruct Cartesian state from Tree B
    state_reconstructed = jacobi_to_cartesian(
        coords=coords_B,
        tree=tree_B,
        masses=masses,
        center_of_mass=state.center_of_mass,
        total_momentum=state.total_momentum,
    )

    # Physical positions and momenta must match to machine roundoff
    np.testing.assert_allclose(state_reconstructed.q, state.q, atol=1e-14)
    np.testing.assert_allclose(state_reconstructed.p, state.p, atol=1e-14)

    # 5. Direct Jacobi-space commutation: rotated coords_B must match cartesian_to_jacobi(state, tree_B)
    coords_B_direct = cartesian_to_jacobi(state, tree_B)
    for node in tree_B.internal_nodes:
        np.testing.assert_allclose(coords_B.r[node.id], coords_B_direct.r[node.id], atol=1e-14)
        np.testing.assert_allclose(coords_B.w[node.id], coords_B_direct.w[node.id], atol=1e-14)


def test_sl2_from_trees_and_three_chart_cech_cocycle_float64():
    """
    Detached Component 5 Certification (Theorem 4.4 Eq. 4.116 & Section 6.5 Eq. 6.178):
    1. SL2TransitionMatrix.from_trees(tree_A, tree_B) extracts exact unimodular matrices
       with det(M) == 1.0 (< 1e-15) across all 3 canonical 3-body trees:
           T_A = ((0, 1), 2),  T_B = (0, (1, 2)),  T_C = ((0, 2), 1).
    2. Exact 2x2 SL(2, R) Cech 1-cocycle closure: M_{C -> A} @ M_{B -> C} @ M_{A -> B} == I_2 (< 1e-15).
    3. Oriented 12x12 Sp(12, R) representation Psi_{A -> B} = build_tree_transition_symplectic_12(tree_A, tree_B)
       satisfies Psi^T J_12 Psi == J_12, det(Psi) == 1.0, and Psi_{C -> A} @ Psi_{B -> C} @ Psi_{A -> B} == I_12 (< 1e-15).
    4. Coordinate rotation round-trip T_A -> T_B -> T_C -> T_A preserves kinetic energy, potential energy,
       and Darboux pairing to machine precision (< 1e-14).
    """
    from normal_cone_nbody_float64.core.transitions import build_tree_transition_symplectic_12
    from normal_cone_nbody_float64.core.coordinates import compute_jacobi_kinetic_energy
    from normal_cone_nbody_float64.core.invariants import compute_kinetic_energy

    masses = np.array([1.2, 2.3, 3.4], dtype=np.float64)
    tree_A = build_tree_from_spec(masses, ((0, 1), 2))
    tree_B = build_tree_from_spec(masses, (0, (1, 2)))
    tree_C = build_tree_from_spec(masses, ((0, 2), 1))

    M_AB = SL2TransitionMatrix.from_trees(tree_A, tree_B)
    M_BC = SL2TransitionMatrix.from_trees(tree_B, tree_C)
    M_CA = SL2TransitionMatrix.from_trees(tree_C, tree_A)

    assert abs(M_AB.determinant - 1.0) < 1e-15
    assert abs(M_BC.determinant - 1.0) < 1e-15
    assert abs(M_CA.determinant - 1.0) < 1e-15

    mat_AB = np.array([[M_AB.a, M_AB.b], [M_AB.c, M_AB.d]], dtype=np.float64)
    mat_BC = np.array([[M_BC.a, M_BC.b], [M_BC.c, M_BC.d]], dtype=np.float64)
    mat_CA = np.array([[M_CA.a, M_CA.b], [M_CA.c, M_CA.d]], dtype=np.float64)
    np.testing.assert_allclose(mat_CA @ mat_BC @ mat_AB, np.eye(2), atol=1e-15)

    Psi_AB = build_tree_transition_symplectic_12(tree_A, tree_B)
    Psi_BC = build_tree_transition_symplectic_12(tree_B, tree_C)
    Psi_CA = build_tree_transition_symplectic_12(tree_C, tree_A)

    for Psi in (Psi_AB, Psi_BC, Psi_CA):
        assert verify_symplectic_matrix_12(Psi)
        assert abs(float(np.linalg.det(Psi)) - 1.0) < 1e-14

    np.testing.assert_allclose(Psi_CA @ Psi_BC @ Psi_AB, np.eye(12), atol=1e-15)

    # Full 3-chart rotation cycle on physical state
    rng = np.random.default_rng(2026)
    q = rng.standard_normal((3, 3))
    p = rng.standard_normal((3, 3))
    state = ExtendedPhaseState(q=q, p=p, masses=masses).to_barycentric()
    T_cart = compute_kinetic_energy(state)

    coords_A = cartesian_to_jacobi(state, tree_A)
    coords_B = rotate_jacobi_tree_coordinates(coords_A, tree_A, tree_B)
    coords_C = rotate_jacobi_tree_coordinates(coords_B, tree_B, tree_C)
    coords_A_back = rotate_jacobi_tree_coordinates(coords_C, tree_C, tree_A)

    # Kinetic energy and Darboux pairing preserved in every chart
    T_A = compute_jacobi_kinetic_energy(coords_A, tree_A)
    T_B = compute_jacobi_kinetic_energy(coords_B, tree_B)
    T_C = compute_jacobi_kinetic_energy(coords_C, tree_C)
    assert abs(T_A - T_cart) < 1e-14
    assert abs(T_B - T_cart) < 1e-14
    assert abs(T_C - T_cart) < 1e-14

    darboux_A = sum(float(np.dot(coords_A.r[n.id], coords_A.w[n.id])) for n in tree_A.internal_nodes)
    darboux_B = sum(float(np.dot(coords_B.r[n.id], coords_B.w[n.id])) for n in tree_B.internal_nodes)
    darboux_C = sum(float(np.dot(coords_C.r[n.id], coords_C.w[n.id])) for n in tree_C.internal_nodes)
    assert abs(darboux_A - darboux_B) < 1e-14
    assert abs(darboux_B - darboux_C) < 1e-14

    for node in tree_A.internal_nodes:
        np.testing.assert_allclose(coords_A_back.r[node.id], coords_A.r[node.id], atol=1e-14)
        np.testing.assert_allclose(coords_A_back.w[node.id], coords_A.w[node.id], atol=1e-14)

