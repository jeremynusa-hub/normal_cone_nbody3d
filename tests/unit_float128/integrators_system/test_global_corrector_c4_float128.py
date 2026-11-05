"""
Unit Tests for Global Conormal Symplectic Corrector GlobalCorrectorC4 in Quadruple Precision (Float128 / mpmath 35 dps).

Mathematical Monograph & Certified Lean 4 Alignment:
- Monograph: Chapter 5, Section 5.5, Theorem 5.8 & Chapter 6, Section 6.5.
- Formal Proof: docs/global_conormal_corrector_pbm4_proof.md
- Specification: docs/well_posed_global_conormal_corrector_pbm4.md
- Certified Lean 4 Module: 07_GlobalConormalCorrectorPBM4_Order6.lean (10 theorems certified without sorry).
- 7 Admissibility Axioms:
    1. Strict Position Stationarity: nabla_p K_4 == 0 <=> Delta q == 0.
    2. Pericenter Conormal Regularity: Closed form regular at r -> 0.
    3. Keplerian Tidal Commutator Cancellation: homological coefficient 1/480.
    4. Composite Apex Curvature Decoupling: zero Cartesian drift.
    5. Exact Lie--Bargmann Galilean Invariance: 10/10 invariants conserved (< 1e-34).
    6. Time-Reversal Parity: C_4(-h) = C_4(h)^{-1}.
    7. Gauge Functor Compatibility: Cech 1-cocycle without residual drift.
"""

from __future__ import annotations
import pytest
import numpy as np
import mpmath as mp

from normal_cone_nbody_float128.core.constants import to_f128, G_DEFAULT
from normal_cone_nbody_float128.core.state import ExtendedPhaseState
from normal_cone_nbody_float128.core.invariants import compute_invariants
from normal_cone_nbody_float128.core.jacobi import JacobiTree, JacobiTopologyNode
from normal_cone_nbody_float128.integrators.correctors import (
    GlobalCorrectorC4,
    CorrectorC4,
    compute_global_conormal_c4_kick,
)


@pytest.fixture(autouse=True)
def setup_mpmath():
    mp.mp.dps = 35
    yield
    mp.mp.dps = 35


def _create_hierarchical_3body_state() -> tuple[ExtendedPhaseState, JacobiTree]:
    """
    Creates a hierarchical 3-body system in Float128:
    Inner binary: bodies 0 and 1 (r = 1.0)
    Outer body: body 2 (R = 8.0)
    """
    m0 = mp.mpf('1.0')
    m1 = mp.mpf('1.0')
    m2 = mp.mpf('1.0')
    masses = np.array([m0, m1, m2], dtype=object)

    # Positions
    q = np.array([
        [mp.mpf('-0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
        [mp.mpf('0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
        [mp.mpf('0.0'), mp.mpf('8.0'), mp.mpf('0.0')],
    ], dtype=object)

    # Momenta (circular-like velocities)
    # v_bin = sqrt(G * M_bin / r) / 2 = sqrt(2.0) / 2 ~= 0.7071
    v_bin = mp.sqrt(mp.mpf('2.0')) / mp.mpf('2.0')
    # v_out = sqrt(G * M_tot / R) ~= sqrt(3.0 / 8.0)
    v_out = mp.sqrt(mp.mpf('3.0') / mp.mpf('8.0'))

    p = np.array([
        [mp.mpf('0.0'), -v_bin, mp.mpf('0.0')],
        [mp.mpf('0.0'), v_bin, mp.mpf('0.0')],
        [v_out, mp.mpf('0.0'), mp.mpf('0.0')],
    ], dtype=object)

    state = ExtendedPhaseState(q=q, p=p, masses=masses, t=mp.mpf('0.0'))

    # Jacobi tree: ((0, 1), 2)
    l0 = JacobiTopologyNode(id=0, is_leaf=True, body_index=0, M=m0, body_indices=(0,))
    l1 = JacobiTopologyNode(id=1, is_leaf=True, body_index=1, M=m1, body_indices=(1,))
    l2 = JacobiTopologyNode(id=2, is_leaf=True, body_index=2, M=m2, body_indices=(2,))

    m01 = m0 + m1
    mu01 = (m0 * m1) / m01
    p01 = JacobiTopologyNode(id=10, is_leaf=False, left=l0, right=l1, body_indices=(0, 1), M=m01, mu=mu01, c_left=m0/m01, c_right=m1/m01)

    m_tot = m01 + m2
    mu_tot = (m01 * m2) / m_tot
    root = JacobiTopologyNode(id=11, is_leaf=False, left=p01, right=l2, body_indices=(0, 1, 2), M=m_tot, mu=mu_tot, c_left=m01/m_tot, c_right=m2/m_tot)

    tree = JacobiTree(
        root=root,
        masses=masses,
        node_map={0: l0, 1: l1, 2: l2, 10: p01, 11: root},
        internal_nodes=(p01, root),
        leaf_pairs=(p01,),
        composite_nodes=(root,),
    )
    return state, tree


def test_algebraic_inversion_exactness():
    """
    Test 1: Exact algebraic inversion C_4^{-1} o C_4 = Id at < 1e-32 in Float128.
    """
    state, tree = _create_hierarchical_3body_state()
    corrector = GlobalCorrectorC4()
    h = mp.mpf('0.05')

    # Forward then inverse
    s_fwd = corrector.apply(state, tree, h, inverse=False)
    s_inv = corrector.apply(s_fwd, tree, h, inverse=True)

    # Position roundtrip error
    max_dq = max(abs(s_inv.q[i, d] - state.q[i, d]) for i in range(3) for d in range(3))
    # Momentum roundtrip error
    max_dp = max(abs(s_inv.p[i, d] - state.p[i, d]) for i in range(3) for d in range(3))

    assert max_dq == mp.mpf('0.0'), f"Positions must not change: max_dq = {max_dq}"
    assert max_dp < mp.mpf('1e-32'), f"Momentum inversion error must be < 1e-32: max_dp = {max_dp}"


def test_spatial_stationarity_absolute():
    """
    Test 2: Strict spatial stationarity nabla_p K_4 == 0 <=> Delta q == 0.0.
    """
    state, tree = _create_hierarchical_3body_state()
    corrector = GlobalCorrectorC4()
    h = mp.mpf('0.05')

    s_applied = corrector.apply(state, tree, h, inverse=False)

    for i in range(3):
        for d in range(3):
            diff = s_applied.q[i, d] - state.q[i, d]
            assert diff == mp.mpf('0.0'), f"Body {i} coordinate {d} moved by {diff}"


def test_total_linear_momentum_conservation():
    """
    Test 3: Strict conservation of total linear momentum: sum Delta p_i == 0 (< 1e-34).
    """
    state, tree = _create_hierarchical_3body_state()
    corrector = GlobalCorrectorC4()
    h = mp.mpf('0.05')

    dp = compute_global_conormal_c4_kick(state, tree, h)

    sum_dp = [sum(dp[i, d] for i in range(3)) for d in range(3)]
    norm_sum_dp = mp.sqrt(sum(c**2 for c in sum_dp))

    assert norm_sum_dp < mp.mpf('1e-34'), f"Total linear momentum kick sum must vanish: {norm_sum_dp}"


def test_center_of_mass_position_stationarity():
    """
    Test 4: Strict stationarity of spatial center of mass: Delta R_cm == 0.0.
    """
    state, tree = _create_hierarchical_3body_state()
    corrector = GlobalCorrectorC4()
    h = mp.mpf('0.05')

    s_applied = corrector.apply(state, tree, h)
    inv0 = compute_invariants(state)
    inv1 = compute_invariants(s_applied)

    dR_cm = mp.sqrt(sum((inv1.R_cm[d] - inv0.R_cm[d])**2 for d in range(3)))
    assert dR_cm == mp.mpf('0.0'), f"Center of mass shifted by {dR_cm}"


def test_total_angular_momentum_conservation():
    """
    Test 5: Conservation of total angular momentum: sum q_i x Delta p_i == 0 (< 1e-34).
    """
    state, tree = _create_hierarchical_3body_state()
    corrector = GlobalCorrectorC4()
    h = mp.mpf('0.05')

    dp = compute_global_conormal_c4_kick(state, tree, h)

    # Compute sum_i (q_i x Delta p_i)
    torque_tot = [mp.mpf('0.0') for _ in range(3)]
    for i in range(3):
        qx, qy, qz = state.q[i, 0], state.q[i, 1], state.q[i, 2]
        dpx, dpy, dpz = dp[i, 0], dp[i, 1], dp[i, 2]
        torque_tot[0] += qy * dpz - qz * dpy
        torque_tot[1] += qz * dpx - qx * dpz
        torque_tot[2] += qx * dpy - qy * dpx

    norm_torque = mp.sqrt(sum(t**2 for t in torque_tot))
    assert norm_torque < mp.mpf('1e-34'), f"Total torque from kick must vanish: {norm_torque}"


def test_time_reversal_parity_and_c4_delegation():
    """
    Test 6: Inversion operator exactness and CorrectorC4 wrapper delegation parity.
    """
    state, tree = _create_hierarchical_3body_state()
    global_c4 = GlobalCorrectorC4()
    wrapper_c4 = CorrectorC4(mode="global")
    h = mp.mpf('0.05')

    # 1. Direct apply_inverse cancelation: C_4(h)^{-1} (C_4(h)(x)) == x
    s_fwd = global_c4.apply(state, tree, h)
    s_rec = global_c4.apply_inverse(s_fwd, tree, h)
    max_rec_err = max(abs(s_rec.p[i, d] - state.p[i, d]) for i in range(3) for d in range(3))
    assert max_rec_err < mp.mpf('1e-32'), f"Recovery error too large: {max_rec_err}"

    # 2. Wrapper forward delegation: CorrectorC4(mode='global') vs GlobalCorrectorC4
    s_wrapper_fwd = wrapper_c4.apply(state, tree, h)
    s_direct_fwd = global_c4.apply(state, tree, h)
    max_diff_fwd = max(abs(s_wrapper_fwd.p[i, d] - s_direct_fwd.p[i, d]) for i in range(3) for d in range(3))
    assert max_diff_fwd == mp.mpf('0.0'), f"Forward wrapper delegation mismatch: {max_diff_fwd}"

    # 3. Wrapper inverse delegation: CorrectorC4(mode='global') vs GlobalCorrectorC4
    s_wrapper_inv = wrapper_c4.apply_inverse(s_fwd, tree, h)
    s_direct_inv = global_c4.apply_inverse(s_fwd, tree, h)
    max_diff_inv = max(abs(s_wrapper_inv.p[i, d] - s_direct_inv.p[i, d]) for i in range(3) for d in range(3))
    assert max_diff_inv == mp.mpf('0.0'), f"Inverse wrapper delegation mismatch: {max_diff_inv}"


def test_five_strata_supernode_masking_and_screening_float128():
    """
    Test 7: Validates Theorem 5.16 (5-Strata Partition) in Float128:
    - Democratic Supernode Stratum D_super: active supernode leaf pairs produce zero kick.
    - Screening Transition Stratum D_trans: spectator body in u in (1, 2) is smoothly attenuated by Theta(u)
      while conserving total linear and angular momentum to < 1e-34.
    """
    from normal_cone_nbody_float128.potential.screened_apex import ScreenedApexPotential

    state, tree = _create_hierarchical_3body_state()
    h = mp.mpf('0.05')

    # 1. Supernode masking
    dp_masked = compute_global_conormal_c4_kick(state, tree, h, supernode_bodies={0, 1, 2})
    for i in range(3):
        for d in range(3):
            assert dp_masked[i, d] == mp.mpf('0.0'), f"Masked supernode kick must be 0.0, got {dp_masked[i, d]}"

    # 2. Screening transition stratum: move body 2 to R = 0.75 (u = 0.75 / 0.5 = 1.5 in (1, 2))
    q_trans = np.array([
        [mp.mpf('-0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
        [mp.mpf('0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
        [mp.mpf('0.1'), mp.mpf('0.74'), mp.mpf('0.05')],
    ], dtype=object)
    state_trans = state.with_updates(q=q_trans)
    apex_screened = ScreenedApexPotential(G=mp.mpf('1.0'), screening_mode="c6")

    dp_unscreened = compute_global_conormal_c4_kick(state_trans, tree, h, apex=None)
    dp_screened = compute_global_conormal_c4_kick(state_trans, tree, h, apex=apex_screened)

    norm_unscreened = mp.sqrt(sum(dp_unscreened[i, d]**2 for i in range(3) for d in range(3)))
    norm_screened = mp.sqrt(sum(dp_screened[i, d]**2 for i in range(3) for d in range(3)))
    assert mp.mpf('0.0') < norm_screened < norm_unscreened

    # Conservation of P and L to < 1e-34
    sum_dp = [sum(dp_screened[i, d] for i in range(3)) for d in range(3)]
    assert mp.sqrt(sum(c**2 for c in sum_dp)) < mp.mpf('1e-34')

    torque_tot = [mp.mpf('0.0') for _ in range(3)]
    for i in range(3):
        qx, qy, qz = state_trans.q[i, 0], state_trans.q[i, 1], state_trans.q[i, 2]
        dpx, dpy, dpz = dp_screened[i, 0], dp_screened[i, 1], dp_screened[i, 2]
        torque_tot[0] += qy * dpz - qz * dpy
        torque_tot[1] += qz * dpx - qx * dpz
        torque_tot[2] += qx * dpy - qy * dpx
    assert mp.sqrt(sum(t**2 for t in torque_tot)) < mp.mpf('1e-34')


def test_normal_cone_canonical_exports_float128():
    """
    Validates D-143-08 in Float128:
    `potential.screened_apex` and `integrators.correctors` expose the canonical
    `compute_normal_cone_*` functions, and `NBodySystem` uses the `normal-cone geometry` warning phrasing.
    """
    import warnings
    import normal_cone_nbody_float128.potential.screened_apex as apex_mod
    import normal_cone_nbody_float128.integrators.correctors as corr_mod
    from normal_cone_nbody_float128.system.nbody_system import NBodySystem

    assert hasattr(apex_mod, "compute_normal_cone_tidal_tensor")
    assert hasattr(corr_mod, "compute_normal_cone_tidal_tensor")
    assert hasattr(corr_mod, "compute_normal_cone_corrector_kick")
    assert hasattr(corr_mod, "compute_global_normal_cone_c4_kick")

    state, _ = _create_hierarchical_3body_state()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _ = NBodySystem(state=state, auto_supernodes=False)
    assert len(caught) == 1
    assert "normal-cone geometry" in str(caught[0].message)


def test_homological_and_cotangent_lift_float128():
    """
    Validates CotangentLiftCorrectorC4 and HomologicalLieCorrectorC4 in 35-digit Float128:
    1. Exact group inversion C_4^{-1}(h) o C_4(h) = Id (< 1e-32).
    2. Exact conservation of Lie-Bargmann invariants P, L, G (< 1e-32).
    """
    from normal_cone_nbody_float128.integrators.correctors import (
        CotangentLiftCorrectorC4,
        HomologicalLieCorrectorC4,
    )
    from normal_cone_nbody_float128.core.invariants import compute_invariants

    state0, tree = _create_hierarchical_3body_state()
    h = mp.mpf('0.10')
    inv0 = compute_invariants(state0)

    for corr in (CotangentLiftCorrectorC4(G=mp.mpf('1.0')), HomologicalLieCorrectorC4(G=mp.mpf('1.0'))):
        st_fwd = corr.apply(state0, tree, h, inverse=False)
        st_back = corr.apply(st_fwd, tree, h, inverse=True)
        for i in range(3):
            for d in range(3):
                assert abs(st_back.q[i, d] - state0.q[i, d]) < mp.mpf('1e-32')
                assert abs(st_back.p[i, d] - state0.p[i, d]) < mp.mpf('1e-32')
        inv_fwd = compute_invariants(st_fwd)
        for d in range(3):
            assert abs(inv_fwd.P[d] - inv0.P[d]) < mp.mpf('1e-32')
            assert abs(inv_fwd.L[d] - inv0.L[d]) < mp.mpf('1e-32')
            assert abs(inv_fwd.G[d] - inv0.G[d]) < mp.mpf('1e-32')



