"""
Unit Tests for PBM4 Symplectic Gauge Transition Operator T_{A -> B}^{(PBM4)} in Quadruple Precision (Float128 / mpmath 35 dps).

Mathematical Monograph & Certified Lean 4 Alignment:
- Chapter 5: Section 5.5 (Conormal Symplectic Processing and Commutators).
- Proof Document: docs/pbm4_symplectic_gauge_connection_proof.md
- Certified Lean 4 Module: 06_ConormalGaugeTransitionPBM4_Order6.lean
- Axioms Audited:
    1. Exact Symplecticity Sp(6N-6, R)
    2. Strict Position Stationarity Delta q == 0
    3. Galilean Lie--Bargmann 10 Invariants (< 1e-34 in f128)
    4. Exact Baker--Campbell--Hausdorff Collapse (Zero Higher Remainder)
    5. Cech 1-Cocycle Transitivity T_{B -> C} o T_{A -> B} = T_{A -> C}
    6. Time-Reversibility / Involution T_{B -> A} o T_{A -> B} = Id
"""

from __future__ import annotations
import pytest
import numpy as np
import mpmath as mp

from normal_cone_nbody_float128.core.constants import to_f128
from normal_cone_nbody_float128.core.state import ExtendedPhaseState
from normal_cone_nbody_float128.core.jacobi import JacobiTree, JacobiTopologyNode
from normal_cone_nbody_float128.core.coordinates import cartesian_to_jacobi
from normal_cone_nbody_float128.core.transitions import apply_pbm4_gauge_transition_jacobi
from normal_cone_nbody_float128.integrators.correctors import CorrectorC4
from normal_cone_nbody_float128.integrators.processed import apply_gauge_connection


def _create_three_body_trees(state: ExtendedPhaseState) -> tuple[JacobiTree, JacobiTree, JacobiTree]:
    """Builds the three canonical non-isomorphic Jacobi trees on 3 bodies:
    Tree A: ((0, 1), 2)
    Tree B: (0, (1, 2))
    Tree C: ((0, 2), 1)
    """
    m0 = to_f128(state.masses[0])
    m1 = to_f128(state.masses[1])
    m2 = to_f128(state.masses[2])
    zero = mp.mpf('0.0')

    # Leaf nodes
    l0 = JacobiTopologyNode(id=0, is_leaf=True, body_index=0, M=m0, body_indices=(0,))
    l1 = JacobiTopologyNode(id=1, is_leaf=True, body_index=1, M=m1, body_indices=(1,))
    l2 = JacobiTopologyNode(id=2, is_leaf=True, body_index=2, M=m2, body_indices=(2,))

    # Tree A: ((0, 1), 2)
    m01 = m0 + m1
    mu01 = (m0 * m1) / m01
    p01_A = JacobiTopologyNode(id=10, is_leaf=False, left=l0, right=l1, body_indices=(0, 1), M=m01, mu=mu01, c_left=m0/m01, c_right=m1/m01)
    m_tot = m01 + m2
    mu_tot_A = (m01 * m2) / m_tot
    root_A = JacobiTopologyNode(id=11, is_leaf=False, left=p01_A, right=l2, body_indices=(0, 1, 2), M=m_tot, mu=mu_tot_A, c_left=m01/m_tot, c_right=m2/m_tot)
    tree_A = JacobiTree(
        root=root_A,
        masses=state.masses,
        node_map={0: l0, 1: l1, 2: l2, 10: p01_A, 11: root_A},
        internal_nodes=(p01_A, root_A),
        leaf_pairs=(p01_A,),
        composite_nodes=(root_A,),
    )

    # Tree B: (0, (1, 2))
    m12 = m1 + m2
    mu12 = (m1 * m2) / m12
    p12_B = JacobiTopologyNode(id=20, is_leaf=False, left=l1, right=l2, body_indices=(1, 2), M=m12, mu=mu12, c_left=m1/m12, c_right=m2/m12)
    mu_tot_B = (m0 * m12) / m_tot
    root_B = JacobiTopologyNode(id=21, is_leaf=False, left=l0, right=p12_B, body_indices=(0, 1, 2), M=m_tot, mu=mu_tot_B, c_left=m0/m_tot, c_right=m12/m_tot)
    tree_B = JacobiTree(
        root=root_B,
        masses=state.masses,
        node_map={0: l0, 1: l1, 2: l2, 20: p12_B, 21: root_B},
        internal_nodes=(p12_B, root_B),
        leaf_pairs=(p12_B,),
        composite_nodes=(root_B,),
    )

    # Tree C: ((0, 2), 1)
    m02 = m0 + m2
    mu02 = (m0 * m2) / m02
    p02_C = JacobiTopologyNode(id=30, is_leaf=False, left=l0, right=l2, body_indices=(0, 2), M=m02, mu=mu02, c_left=m0/m02, c_right=m2/m02)
    mu_tot_C = (m02 * m1) / m_tot
    root_C = JacobiTopologyNode(id=31, is_leaf=False, left=p02_C, right=l1, body_indices=(0, 1, 2), M=m_tot, mu=mu_tot_C, c_left=m02/m_tot, c_right=m1/m_tot)
    tree_C = JacobiTree(
        root=root_C,
        masses=state.masses,
        node_map={0: l0, 1: l1, 2: l2, 30: p02_C, 31: root_C},
        internal_nodes=(p02_C, root_C),
        leaf_pairs=(p02_C,),
        composite_nodes=(root_C,),
    )

    return tree_A, tree_B, tree_C


@pytest.fixture
def f128_three_body_state() -> ExtendedPhaseState:
    """Fixture providing a hierarchical 3-body system in quadruple precision."""
    with mp.workdps(35):
        masses = np.array([mp.mpf('1.0'), mp.mpf('1.0'), mp.mpf('0.5')], dtype=object)
        q = np.array([
            [mp.mpf('0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('-0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.0'), mp.mpf('5.0'), mp.mpf('0.2')],
        ], dtype=object)
        p = np.array([
            [mp.mpf('0.0'), mp.mpf('0.7'), mp.mpf('0.1')],
            [mp.mpf('0.0'), mp.mpf('-0.7'), mp.mpf('-0.1')],
            [mp.mpf('0.1'), mp.mpf('0.0'), mp.mpf('0.0')],
        ], dtype=object)
        return ExtendedPhaseState(q=q, p=p, masses=masses, t=mp.mpf('0.0'))


def test_pbm4_gauge_position_stationarity(f128_three_body_state):
    """Axiom 2: Strict spatial stationarity Delta q == 0 for unipotent P_0 mode (while default is homological)."""
    tree_A, tree_B, _ = _create_three_body_trees(f128_three_body_state)
    assert CorrectorC4(G=mp.mpf('1.0')).mode == "homological"
    corrector = CorrectorC4(G=mp.mpf('1.0'), mode="global")
    dt = mp.mpf('0.05')

    st_trans = apply_gauge_connection(f128_three_body_state, tree_A, tree_B, dt=dt, corrector=corrector)

    # Positions must be bit-exact identical
    for i in range(3):
        for d in range(3):
            assert st_trans.q[i, d] == f128_three_body_state.q[i, d]


def test_pbm4_gauge_involutivity_and_reversibility(f128_three_body_state):
    """Axiom 6: Time-reversibility and involution T_{B -> A} o T_{A -> B} = Id (< 1e-34)."""
    tree_A, tree_B, _ = _create_three_body_trees(f128_three_body_state)
    corrector = CorrectorC4(G=mp.mpf('1.0'), mode="global")
    dt = mp.mpf('0.05')

    st_B = apply_gauge_connection(f128_three_body_state, tree_A, tree_B, dt=dt, corrector=corrector)
    st_roundtrip = apply_gauge_connection(st_B, tree_B, tree_A, dt=dt, corrector=corrector)

    for i in range(3):
        for d in range(3):
            diff = abs(st_roundtrip.p[i, d] - f128_three_body_state.p[i, d])
            assert diff < mp.mpf('1e-34'), f"Roundtrip momentum error {diff} exceeds 1e-34 at ({i}, {d})"


def test_pbm4_gauge_lie_bargmann_invariants(f128_three_body_state):
    """Axiom 3: Exact conservation of linear and angular momentum (< 1e-34)."""
    tree_A, tree_B, _ = _create_three_body_trees(f128_three_body_state)
    corrector = CorrectorC4(G=mp.mpf('1.0'), mode="global")
    dt = mp.mpf('0.05')

    st_B = apply_gauge_connection(f128_three_body_state, tree_A, tree_B, dt=dt, corrector=corrector)

    # 1. Total linear momentum conservation: sum_i Delta p_i == 0
    dp_tot = [mp.mpf('0.0') for _ in range(3)]
    for i in range(3):
        for d in range(3):
            dp_tot[d] += (st_B.p[i, d] - f128_three_body_state.p[i, d])
    for d in range(3):
        assert abs(dp_tot[d]) < mp.mpf('1e-34'), f"Total linear momentum jump {dp_tot[d]} exceeds 1e-34"

    # 2. Total angular momentum conservation: sum_i q_i x Delta p_i == 0
    dL_tot = [mp.mpf('0.0') for _ in range(3)]
    for i in range(3):
        q_i = f128_three_body_state.q[i]
        dp_i = [st_B.p[i, d] - f128_three_body_state.p[i, d] for d in range(3)]
        dL_tot[0] += (q_i[1] * dp_i[2] - q_i[2] * dp_i[1])
        dL_tot[1] += (q_i[2] * dp_i[0] - q_i[0] * dp_i[2])
        dL_tot[2] += (q_i[0] * dp_i[1] - q_i[1] * dp_i[0])
    for d in range(3):
        assert abs(dL_tot[d]) < mp.mpf('1e-34'), f"Total angular momentum jump {dL_tot[d]} exceeds 1e-34"


def test_pbm4_gauge_cech_cocycle_transitivity(f128_three_body_state):
    """Axiom 5: Cech 1-cocycle transitivity T_{B -> C} o T_{A -> B} = T_{A -> C} (< 1e-34)."""
    tree_A, tree_B, tree_C = _create_three_body_trees(f128_three_body_state)
    corrector = CorrectorC4(G=mp.mpf('1.0'), mode="global")
    dt = mp.mpf('0.05')

    # Two-step transition: A -> B -> C
    st_B = apply_gauge_connection(f128_three_body_state, tree_A, tree_B, dt=dt, corrector=corrector)
    st_C_twostep = apply_gauge_connection(st_B, tree_B, tree_C, dt=dt, corrector=corrector)

    # Direct transition: A -> C
    st_C_direct = apply_gauge_connection(f128_three_body_state, tree_A, tree_C, dt=dt, corrector=corrector)

    for i in range(3):
        for d in range(3):
            diff = abs(st_C_twostep.p[i, d] - st_C_direct.p[i, d])
            assert diff < mp.mpf('1e-34'), f"Cech cocycle error {diff} exceeds 1e-34 at ({i}, {d})"


def test_pbm4_gauge_momentum_independence(f128_three_body_state):
    """Axiom 4: Momentum independence grad_w K_4 == 0 => jump is independent of initial momentum."""
    tree_A, tree_B, _ = _create_three_body_trees(f128_three_body_state)
    corrector = CorrectorC4(G=mp.mpf('1.0'), mode="global")
    dt = mp.mpf('0.05')

    # State 1 with initial momentum p
    st_B_1 = apply_gauge_connection(f128_three_body_state, tree_A, tree_B, dt=dt, corrector=corrector)
    dp_1 = st_B_1.p - f128_three_body_state.p

    # State 2 with scaled initial momentum 5.0 * p
    p_scaled = np.array([
        [mp.mpf('5.0') * f128_three_body_state.p[i, d] for d in range(3)]
        for i in range(3)
    ], dtype=object)
    st_scaled = f128_three_body_state.with_updates(p=p_scaled)
    st_B_2 = apply_gauge_connection(st_scaled, tree_A, tree_B, dt=dt, corrector=corrector)
    dp_2 = st_B_2.p - st_scaled.p

    for i in range(3):
        for d in range(3):
            diff = abs(dp_1[i, d] - dp_2[i, d])
            assert diff < mp.mpf('1e-34'), f"Momentum jump dependent on initial momentum: {diff}"


def test_pbm4_gauge_jacobi_coordinates_isomorphism(f128_three_body_state):
    """Tests the Jacobi-coordinate representation apply_pbm4_gauge_transition_jacobi."""
    tree_A, tree_B, _ = _create_three_body_trees(f128_three_body_state)
    corrector = CorrectorC4(G=mp.mpf('1.0'), mode="global")
    dt = mp.mpf('0.05')

    # Cartesian reference
    st_B_cart = apply_gauge_connection(f128_three_body_state, tree_A, tree_B, dt=dt, corrector=corrector)
    coords_B_cart = cartesian_to_jacobi(st_B_cart, tree_B)

    # Jacobi coordinate direct transition
    coords_A = cartesian_to_jacobi(f128_three_body_state, tree_A)
    coords_B_jac = apply_pbm4_gauge_transition_jacobi(
        coords_A, tree_A, tree_B, f128_three_body_state, dt=dt, G=mp.mpf('1.0')
    )

    # Verify positions and momenta match across all nodes in Tree B
    for node_id in tree_B.internal_nodes:
        nid = node_id.id
        r_cart = coords_B_cart.r[nid]
        r_jac = coords_B_jac.r[nid]
        w_cart = coords_B_cart.w[nid]
        w_jac = coords_B_jac.w[nid]
        for d in range(3):
            assert abs(r_cart[d] - r_jac[d]) < mp.mpf('1e-34')
            assert abs(w_cart[d] - w_jac[d]) < mp.mpf('1e-34')


def test_pbm4_gauge_five_strata_supernode_masking_and_screening(f128_three_body_state):
    """Validates Theorem 5.16 & Theorem 5.17: apply_gauge_connection and apply_pbm4_gauge_transition_jacobi
    respect active supernode leaf-pair masking and screening attenuation in 35-digit arithmetic."""
    from normal_cone_nbody_float128.potential.screened_apex import ScreenedApexPotential

    tree_A, tree_B, _ = _create_three_body_trees(f128_three_body_state)
    corrector = CorrectorC4(G=mp.mpf('1.0'), mode="global")
    dt = mp.mpf('0.05')

    # When both leaf pairs are inside an active supernode {0, 1, 2}, the gauge jump vanishes identically
    st_masked = apply_gauge_connection(
        f128_three_body_state, tree_A, tree_B, dt=dt, corrector=corrector, supernode_bodies={0, 1, 2}
    )
    for i in range(3):
        for d in range(3):
            assert abs(st_masked.p[i, d] - f128_three_body_state.p[i, d]) == mp.mpf('0.0')

    # Equivalence with explicit sandwich C_B^{-1} o C_A when apex screening is active
    apex = ScreenedApexPotential(G=mp.mpf('1.0'), screening_mode="c6")
    st_unipotent = apply_gauge_connection(
        f128_three_body_state, tree_A, tree_B, dt=dt, corrector=corrector, apex=apex, supernode_bodies={0, 1}
    )
    st_sandwich = corrector.apply_inverse(
        corrector.apply(f128_three_body_state, tree_A, dt, inverse=False, apex=apex, supernode_bodies={0, 1}),
        tree_B, dt, apex=apex, supernode_bodies={0, 1}
    )
    for i in range(3):
        for d in range(3):
            assert abs(st_unipotent.p[i, d] - st_sandwich.p[i, d]) < mp.mpf('1e-34')


def test_sl2_and_homological_gauge_cech_cocycle_and_dynamic_jump_float128():
    """
    Detached Component 5 Certification in 35-digit Float128 (Theorems 4.4, 5.17, 5.18 & Section 6.5):
    1. Exact SL(2, R) unimodularity |det(M) - 1| < 1e-34 and 3-chart Sp(12, R) Cech 1-cocycle closure
       ||Psi_{C -> A} Psi_{B -> C} Psi_{A -> B} - I_12||_inf < 1e-34.
    2. Exact non-commutative gauge 1-cocycle closure with default mode='homological' (Z_4^{BM4} in P_1)
       across (T_A, h_A) -> (T_B, h_B) -> (T_C, h_C) -> (T_A, h_A) to < 1e-33, and conservation of
       all 9 Lie--Bargmann generators (P, L, G) to < 1e-33.
    3. Dynamic multi-step interfacial jump suppression ([CRIT-7.1.2-01]):
       On the 3D eccentric hierarchical triple in Float128, integrating M_1 = 8 steps at h_A = 0.05,
       switching to h_B = 0.025 via T_{(A, h_A) -> (A, h_B)}^{(PBM4)}, and integrating M_2 = 16 steps at h_B
       reduces the post-transition energy error by > 500x compared to a naive transition or the obstructed P_0 kick.
    """
    from normal_cone_nbody_float128.core.transitions import (
        SL2TransitionMatrix,
        build_tree_transition_symplectic_12,
        verify_symplectic_matrix_12,
        rotate_jacobi_tree_coordinates,
    )
    from normal_cone_nbody_float128.core.invariants import compute_invariants
    from normal_cone_nbody_float128.integrators.processed import (
        PBM4GaugeTransitionOperator,
        create_pbm4_integrator,
    )

    with mp.workdps(35):
        masses = np.array([mp.mpf('0.5'), mp.mpf('0.5'), mp.mpf('1.0')], dtype=object)
        a_in = mp.mpf('1.0')
        e_in = mp.mpf('0.30')
        inc = mp.mpf('28.0') * mp.pi / mp.mpf('180.0')
        r_p = a_in * (mp.mpf('1.0') - e_in)
        v_p = mp.sqrt((mp.mpf('1.0') + e_in) / (a_in * (mp.mpf('1.0') - e_in)))
        R_out = mp.mpf('30.0')
        v_out = mp.sqrt(mp.mpf('2.0') / R_out)

        q_rel_in = np.array([r_p, mp.mpf('0.0'), mp.mpf('0.0')], dtype=object)
        v_rel_in = np.array([mp.mpf('0.0'), v_p * mp.cos(inc), v_p * mp.sin(inc)], dtype=object)
        q0 = np.array([
            [-mp.mpf('0.5') * q_rel_in[0], -mp.mpf('0.5') * q_rel_in[1], -mp.mpf('0.5') * q_rel_in[2]],
            [+mp.mpf('0.5') * q_rel_in[0], +mp.mpf('0.5') * q_rel_in[1], +mp.mpf('0.5') * q_rel_in[2]],
            [R_out * mp.mpf('0.8'), R_out * mp.mpf('0.6'), mp.mpf('0.0')],
        ], dtype=object)
        v0 = np.array([
            [-mp.mpf('0.5') * v_rel_in[0], -mp.mpf('0.5') * v_rel_in[1], -mp.mpf('0.5') * v_rel_in[2]],
            [+mp.mpf('0.5') * v_rel_in[0], +mp.mpf('0.5') * v_rel_in[1], +mp.mpf('0.5') * v_rel_in[2]],
            [-mp.mpf('0.6') * v_out, mp.mpf('0.8') * v_out, mp.mpf('0.0')],
        ], dtype=object)
        p0 = np.array([[masses[i] * v0[i, d] for d in range(3)] for i in range(3)], dtype=object)
        state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=mp.mpf('0.0')).to_barycentric()

        tree_A, tree_B, tree_C = _create_three_body_trees(state0)

        # 1. Exact SL(2, R) unimodularity and 12x12 Sp(12, R) Cech 1-cocycle closure in Float128
        M_AB = SL2TransitionMatrix.from_trees(tree_A, tree_B)
        M_BC = SL2TransitionMatrix.from_trees(tree_B, tree_C)
        M_CA = SL2TransitionMatrix.from_trees(tree_C, tree_A)
        assert M_AB.is_f128 and abs(M_AB.determinant - mp.mpf('1.0')) < mp.mpf('1e-34')
        assert M_BC.is_f128 and abs(M_BC.determinant - mp.mpf('1.0')) < mp.mpf('1e-34')
        assert M_CA.is_f128 and abs(M_CA.determinant - mp.mpf('1.0')) < mp.mpf('1e-34')

        Psi_AB = build_tree_transition_symplectic_12(tree_A, tree_B)
        Psi_BC = build_tree_transition_symplectic_12(tree_B, tree_C)
        Psi_CA = build_tree_transition_symplectic_12(tree_C, tree_A)
        for Psi in (Psi_AB, Psi_BC, Psi_CA):
            assert verify_symplectic_matrix_12(Psi, tol=mp.mpf('1e-33'))

        Psi_cycle = Psi_CA @ Psi_BC @ Psi_AB
        for i in range(12):
            for j in range(12):
                expected = mp.mpf('1.0') if i == j else mp.mpf('0.0')
                assert abs(Psi_cycle[i, j] - expected) < mp.mpf('1e-34')

        coords_A = cartesian_to_jacobi(state0, tree_A)
        coords_B = rotate_jacobi_tree_coordinates(coords_A, tree_A, tree_B)
        coords_C = rotate_jacobi_tree_coordinates(coords_B, tree_B, tree_C)
        coords_A_back = rotate_jacobi_tree_coordinates(coords_C, tree_C, tree_A)
        for node in tree_A.internal_nodes:
            for d in range(3):
                assert abs(coords_A_back.r[node.id][d] - coords_A.r[node.id][d]) < mp.mpf('1e-33')
                assert abs(coords_A_back.w[node.id][d] - coords_A.w[node.id][d]) < mp.mpf('1e-33')

        # 2. Non-commutative 3-chart + 3-scale Cech 1-cocycle closure in Float128 (mode="homological")
        gauge_op = PBM4GaugeTransitionOperator(G=mp.mpf('1.0'))
        assert getattr(gauge_op.corrector, "mode", None) == "homological"

        cocycle_err = gauge_op.verify_cech_cocycle(
            state0,
            tree_A,
            tree_B,
            tree_C,
            dt_A=mp.mpf('0.04'),
            dt_B=mp.mpf('0.02'),
            dt_C=mp.mpf('0.01'),
        )
        assert cocycle_err < mp.mpf('1e-33'), f"Float128 non-commutative Cech cocycle error: {cocycle_err}"

        inv0 = compute_invariants(state0)
        s_AB = apply_gauge_connection(state0, tree_A, tree_B, dt=mp.mpf('0.04'), dt_new=mp.mpf('0.02'))
        inv_AB = compute_invariants(s_AB)
        for d in range(3):
            assert abs(inv_AB.P[d] - inv0.P[d]) < mp.mpf('1e-33')
            assert abs(inv_AB.L[d] - inv0.L[d]) < mp.mpf('1e-33')
            assert abs(inv_AB.G[d] - inv0.G[d]) < mp.mpf('1e-33')

        # 3. Dynamic multi-step interfacial jump suppression in Float128
        pbm4_hom = create_pbm4_integrator(G=mp.mpf('1.0'), mode="homological")
        pbm4_p0 = create_pbm4_integrator(G=mp.mpf('1.0'), mode="global")

        h_A = mp.mpf('0.05')
        h_B = mp.mpf('0.025')
        M_1 = 8
        M_2 = 16
        E0 = inv0.H

        sh_hom_1 = pbm4_hom.enter_shadow(state0, tree_A, h_A)
        sh_p0_1 = pbm4_p0.enter_shadow(state0, tree_A, h_A)
        for _ in range(M_1):
            sh_hom_1 = pbm4_hom.step_shadow(sh_hom_1, tree_A, h_A)
            sh_p0_1 = pbm4_p0.step_shadow(sh_p0_1, tree_A, h_A)

        sh_naive = sh_hom_1
        sh_p0_trans = apply_gauge_connection(
            sh_p0_1, tree_A, tree_A, dt=h_A, dt_new=h_B, corrector=pbm4_p0.corrector
        )
        sh_hom_trans = apply_gauge_connection(sh_hom_1, tree_A, tree_A, dt=h_A, dt_new=h_B)

        max_err_naive = mp.mpf('0.0')
        max_err_p0 = mp.mpf('0.0')
        max_err_hom = mp.mpf('0.0')

        for _ in range(M_2):
            sh_naive = pbm4_hom.step_shadow(sh_naive, tree_A, h_B)
            sh_p0_trans = pbm4_p0.step_shadow(sh_p0_trans, tree_A, h_B)
            sh_hom_trans = pbm4_hom.step_shadow(sh_hom_trans, tree_A, h_B)

            phys_naive = pbm4_hom.observe_physical(sh_naive, tree_A, h_B)
            phys_p0 = pbm4_p0.observe_physical(sh_p0_trans, tree_A, h_B)
            phys_hom = pbm4_hom.observe_physical(sh_hom_trans, tree_A, h_B)

            max_err_naive = max(max_err_naive, abs(compute_invariants(phys_naive).H - E0) / abs(E0))
            max_err_p0 = max(max_err_p0, abs(compute_invariants(phys_p0).H - E0) / abs(E0))
            max_err_hom = max(max_err_hom, abs(compute_invariants(phys_hom).H - E0) / abs(E0))

        suppression_vs_naive = float(max_err_naive / max_err_hom)
        suppression_vs_p0 = float(max_err_p0 / max_err_hom)
        assert suppression_vs_naive > 500.0, f"Expected > 500x jump suppression vs naive in Float128, got {suppression_vs_naive:.1f}x"
        assert suppression_vs_p0 > 500.0, f"Expected > 500x jump suppression vs P_0 in Float128, got {suppression_vs_p0:.1f}x"

