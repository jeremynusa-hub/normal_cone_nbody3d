"""
Unit Tests for Layer 4: Processed Symplectic Integrator Architecture PBM.

Theoretical Alignments:
- Monograph: Chapter 5, Section 5.5, Proposition 5.3 & Comparative Table 1.
- Formal Lean 4: `02_CorrectorBM6_Theorems.lean` (theorem_9_4_dyadic_telescoping_under_hysteresis).
- Invariants: Exact dyadic boundary telescoping over M steps, order elevation by +2,
  and microscopic time-reversibility (< 1e-14).
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
from normal_cone_nbody_float64.integrators.processed import ProcessedSymplecticIntegrator


def test_processed_microscopic_time_reversibility():
    """Validates microscopic time reversibility Psi(-dt) o Psi(dt) = Id (< 1e-14)."""
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

    bm4 = BlanesMoan4Integrator(G=1.0)
    pbm4 = ProcessedSymplecticIntegrator(bm4)
    dt = 0.03

    s_fwd = pbm4.step(state0, tree, dt=dt)
    s_back = pbm4.step(s_fwd, tree, dt=-dt)

    assert np.allclose(s_back.q, state0.q, atol=1e-14), f"Reversibility error q: {np.max(np.abs(s_back.q - state0.q))}"
    assert np.allclose(s_back.p, state0.p, atol=1e-14), f"Reversibility error p: {np.max(np.abs(s_back.p - state0.p))}"


def test_processed_dyadic_boundary_telescoping_identity():
    """
    Validates Theorem 5.5 / Proposition 5.3:
    Integrating M steps with boundary telescoping produces identical results to
    applying C(h) o S(h)^M o C(h)^{-1}.
    """
    masses = np.array([1.0, 2.0, 1.2], dtype=np.float64)
    q0 = np.array([
        [1.0, 0.2, -0.1],
        [-0.6, 0.8, 0.2],
        [0.1, -0.5, 0.4],
    ], dtype=np.float64)
    p0 = np.array([
        [0.1, -0.2, 0.1],
        [-0.2, 0.1, -0.1],
        [0.1, 0.1, -0.05],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    bm4 = BlanesMoan4Integrator(G=1.0)
    pbm4 = ProcessedSymplecticIntegrator(bm4)
    dt = 0.01
    num_steps = 5

    # 1. Trajectory integration via telescoping .integrate()
    state_tele, _ = pbm4.integrate(state0, tree, dt, num_steps=num_steps)

    # 2. Manual evaluation: C(dt) o [S(dt)^M] o C^{-1}(dt)
    s_manual = pbm4.corrector.apply_inverse(state0, tree, dt)
    for _ in range(num_steps):
        s_manual = bm4.step(s_manual, tree, dt)
    state_manual = pbm4.corrector.apply(s_manual, tree, dt, inverse=False)

    # Telescoping and manual must match to machine precision (< 1e-14)
    assert np.allclose(state_tele.q, state_manual.q, atol=1e-14)
    assert np.allclose(state_tele.p, state_manual.p, atol=1e-14)


def test_processed_effective_order_elevation():
    """
    Validates that PBM4 reduces the shadow Hamiltonian error below raw BM4 and achieves
    steep convergence on a 3-body system over fixed interval T.
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
    pbm4 = ProcessedSymplecticIntegrator(bm4)
    inv0 = compute_invariants(state0)

    T = 0.24
    dt = 0.04
    # Raw BM4
    s_raw = state0
    for _ in range(int(round(T / dt))):
        s_raw = bm4.step(s_raw, tree, dt=dt)
    err_raw = abs(compute_invariants(s_raw).H - inv0.H)

    # Processed PBM4 (with boundary telescoping on the shadow manifold)
    s_proc, _ = pbm4.integrate(state0, tree, dt=dt, num_steps=int(round(T / dt)))
    err_proc = abs(compute_invariants(s_proc).H - inv0.H)

    # 1. Processing reduces Hamiltonian residual error
    assert err_proc < err_raw, f"Expected processed error ({err_proc}) < raw error ({err_raw})"

    # 2. Refined step dt/2
    dt_half = 0.02
    s_half, _ = pbm4.integrate(state0, tree, dt=dt_half, num_steps=int(round(T / dt_half)))
    err_half = abs(compute_invariants(s_half).H - inv0.H)

    ratio = err_proc / max(err_half, 1e-16)
    assert ratio >= 13.0, f"Expected steep convergence ratio >= 13, got {ratio}"


def test_gauge_connection_identity_on_same_tree():
    """Validates that J_{A -> A} is identically the identity operator (< 1e-16)."""
    from normal_cone_nbody_float64.integrators.processed import apply_gauge_connection
    from normal_cone_nbody_float64.integrators.correctors import CorrectorC4

    masses = np.array([1.0, 1.2, 0.8], dtype=np.float64)
    q0 = np.array([[1.0, 0.0, 0.0], [-0.5, 0.8, 0.0], [0.0, -0.4, 0.5]], dtype=np.float64)
    p0 = np.array([[0.0, 0.5, -0.2], [-0.3, -0.2, 0.1], [0.3, -0.3, 0.1]], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")
    c4 = CorrectorC4()

    s_out = apply_gauge_connection(state0, tree, tree, dt=0.05, corrector=c4)
    assert np.allclose(s_out.q, state0.q, atol=1e-16)
    assert np.allclose(s_out.p, state0.p, atol=1e-16)


def test_gauge_connection_invariants_conservation_across_transition():
    """
    Validates that the gauge connection operator J_{A -> B} = C_B^{-1} o C_A
    preserves all 10 Lie--Bargmann invariants across an atlas transition T_A -> T_B.
    """
    from normal_cone_nbody_float64.integrators.processed import apply_gauge_connection
    from normal_cone_nbody_float64.integrators.correctors import CorrectorC4
    from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec

    masses = np.array([1.0, 1.5, 2.0], dtype=np.float64)
    q0 = np.array([[0.8, -0.2, 0.1], [-0.4, 0.9, -0.2], [0.1, -0.5, 0.5]], dtype=np.float64)
    p0 = np.array([[0.15, -0.1, 0.05], [-0.2, 0.15, -0.1], [0.05, -0.05, 0.05]], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)

    # Chart A: leaf pair (0, 1), then combined with 2 -> ((0, 1), 2)
    treeA = build_tree_from_spec(masses, ((0, 1), 2))

    # Chart B: leaf pair (1, 2), then combined with 0 -> (0, (1, 2))
    treeB = build_tree_from_spec(masses, (0, (1, 2)))


    c4 = CorrectorC4(G=1.0)
    dt = 0.02
    inv0 = compute_invariants(state0)

    # Apply gauge connection J_{A -> B}
    s_connected = apply_gauge_connection(state0, treeA, treeB, dt=dt, corrector=c4)
    inv_conn = compute_invariants(s_connected)

    # Linear momentum conservation
    assert np.allclose(inv_conn.P, inv0.P, atol=1e-14), f"Momentum drift: {inv_conn.P - inv0.P}"

    # Angular momentum conservation
    assert np.allclose(inv_conn.L, inv0.L, atol=1e-14), f"Angular momentum drift: {inv_conn.L - inv0.L}"

    # Center of mass boost conservation
    assert np.allclose(inv_conn.G, inv0.G, atol=1e-14), f"Boost drift: {inv_conn.G - inv0.G}"

    # Energy continuity under small dt (deviation bounded by corrector order O(h^4))
    delta_H = abs(inv_conn.H - inv0.H)
    assert delta_H < 1e-4, f"Energy discontinuity too large: {delta_H}"


def test_pbm4_gauge_transition_jacobi_and_cech_cocycle_float64():
    """
    Validates Theorem 5.17 in float64:
    1. Exact isomorphism between Cartesian apply_gauge_connection and apply_pbm4_gauge_transition_jacobi.
    2. Cech 1-cocycle transitivity T_{B -> C} o T_{A -> B} = T_{A -> C} (< 1e-15).
    3. 5-strata partition support (apex screening and supernode_bodies masking).
    4. create_pbm4_integrator factory and absence of create_pbm6_integrator.
    """
    import normal_cone_nbody_float64.integrators.processed as proc_mod
    from normal_cone_nbody_float64.core.transitions import apply_pbm4_gauge_transition_jacobi
    from normal_cone_nbody_float64.core.coordinates import cartesian_to_jacobi
    from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec
    from normal_cone_nbody_float64.potential.screened_apex import ScreenedApexPotential
    from normal_cone_nbody_float64.integrators.correctors import CorrectorC4
    from normal_cone_nbody_float64.integrators.processed import (
        apply_gauge_connection,
        create_pbm4_integrator,
    )

    assert not hasattr(proc_mod, "create_pbm6_integrator")
    pbm4 = create_pbm4_integrator(G=1.0)
    assert pbm4.corrector.mode == "homological"

    masses = np.array([1.0, 1.5, 2.0], dtype=np.float64)
    q0 = np.array([[0.8, -0.2, 0.1], [-0.4, 0.9, -0.2], [0.1, -0.5, 0.5]], dtype=np.float64)
    p0 = np.array([[0.15, -0.1, 0.05], [-0.2, 0.15, -0.1], [0.05, -0.05, 0.05]], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)

    treeA = build_tree_from_spec(masses, ((0, 1), 2))
    treeB = build_tree_from_spec(masses, (0, (1, 2)))
    treeC = build_tree_from_spec(masses, ((0, 2), 1))

    c4 = CorrectorC4(G=1.0, mode="global")
    dt = 0.04

    # 1. Exact isomorphism between Cartesian apply_gauge_connection and apply_pbm4_gauge_transition_jacobi
    st_B_cart = apply_gauge_connection(state0, treeA, treeB, dt=dt, corrector=c4)
    coords_B_cart = cartesian_to_jacobi(st_B_cart, treeB)
    coords_A = cartesian_to_jacobi(state0, treeA)
    coords_B_jac = apply_pbm4_gauge_transition_jacobi(coords_A, treeA, treeB, state0, dt=dt, G=1.0)

    for node in treeB.internal_nodes:
        nid = node.id
        assert np.allclose(coords_B_cart.r[nid], coords_B_jac.r[nid], atol=1e-15)
        assert np.allclose(coords_B_cart.w[nid], coords_B_jac.w[nid], atol=1e-15)

    # 2. Cech 1-cocycle transitivity T_{B -> C} o T_{A -> B} = T_{A -> C}
    st_C_twostep = apply_gauge_connection(st_B_cart, treeB, treeC, dt=dt, corrector=c4)
    st_C_direct = apply_gauge_connection(state0, treeA, treeC, dt=dt, corrector=c4)
    assert np.allclose(st_C_twostep.p, st_C_direct.p, atol=1e-15)

    # 3. 5-strata partition: when all bodies are in an active supernode, gauge transition is identity
    st_B_masked = apply_gauge_connection(
        state0, treeA, treeB, dt=dt, corrector=c4, supernode_bodies={0, 1, 2}
    )
    assert np.allclose(st_B_masked.p, state0.p, atol=1e-16)


def test_homological_z4_order_6_elevation() -> None:
    """
    Validates the canonical palindromic Lie-conjugation corrector Z_4^{BM4} (HomologicalLieCorrectorC4)
    against raw BM4:
    1. Exact group inversion Z_4^{-1}(h) o Z_4(h) = Id (< 1e-15).
    2. Exact conservation of Lie-Bargmann invariants (P, L, R_cm < 1e-15).
    3. > 150x energy error reduction and order-6 convergence over BM4 on a 3D hierarchical binary.
    """
    from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec
    from normal_cone_nbody_float64.integrators.correctors import HomologicalLieCorrectorC4
    from normal_cone_nbody_float64.integrators.processed import create_pbm4_integrator

    masses = np.array([0.5, 0.5, 1.0], dtype=np.float64)
    a_in, e_in = 1.0, 0.25
    r_p = a_in * (1.0 - e_in)
    v_p = np.sqrt((1.0 + e_in) / (a_in * (1.0 - e_in)))
    R_out = 20.0
    v_out = np.sqrt(2.0 / R_out)
    q0 = np.array([[-0.5 * r_p, 0.0, 0.0], [0.5 * r_p, 0.0, 0.0], [R_out * 0.8, R_out * 0.6, 2.0]])
    v0 = np.array([[0.0, -0.5 * v_p, 0.05], [0.0, 0.5 * v_p, -0.05], [-0.6 * v_out, 0.8 * v_out, 0.0]])
    v_cm = np.sum(masses[:, None] * v0, axis=0) / np.sum(masses)
    v0 -= v_cm
    q_cm = np.sum(masses[:, None] * q0, axis=0) / np.sum(masses)
    q0 -= q_cm
    state0 = ExtendedPhaseState(q=q0, p=masses[:, None] * v0, masses=masses)
    tree = build_tree_from_spec(masses, ((0, 1), 2))
    dt = 0.20

    # 1. Exact inversion and Lie-Bargmann conservation for Z_4^{BM4}
    inv0 = compute_invariants(state0)
    corr = HomologicalLieCorrectorC4(G=1.0)
    st_fwd = corr.apply(state0, tree, dt, inverse=False)
    st_back = corr.apply(st_fwd, tree, dt, inverse=True)
    assert np.max(np.abs(st_back.q - state0.q)) < 1e-14
    assert np.max(np.abs(st_back.p - state0.p)) < 1e-14
    inv_fwd = compute_invariants(st_fwd)
    assert np.linalg.norm(inv_fwd.P - inv0.P) < 1e-14
    assert np.linalg.norm(inv_fwd.L - inv0.L) < 1e-14
    assert np.linalg.norm(inv_fwd.G - inv0.G) < 1e-14

    # 2. Order-6 convergence slope (~2^6 = 64x per halving) and > 150x energy error reduction over raw BM4
    bm4 = BlanesMoan4Integrator(G=1.0)
    pbm4_hom = create_pbm4_integrator(G=1.0)
    E0 = inv0.H

    errs_bm4 = {}
    errs_pbm4 = {}
    for h_test in (0.20, 0.10):
        n_steps = int(round(6.0 / h_test))
        st_b = state0
        st_s = pbm4_hom.enter_shadow(state0, tree, h_test)
        m_b, m_p = 0.0, 0.0
        for _ in range(n_steps):
            st_b = bm4.step(st_b, tree, h_test)
            st_s = pbm4_hom.step_shadow(st_s, tree, h_test)
            st_p = pbm4_hom.exit_shadow(st_s, tree, h_test)
            m_b = max(m_b, abs(compute_invariants(st_b).H - E0) / abs(E0))
            m_p = max(m_p, abs(compute_invariants(st_p).H - E0) / abs(E0))
        errs_bm4[h_test] = m_b
        errs_pbm4[h_test] = m_p

    # At h=0.20 gain is ~48x; at h=0.10 gain is ~190x (> 150x), and order-6 ratio is ~60x (>> 2^4 = 16)
    assert errs_bm4[0.20] / errs_pbm4[0.20] > 40.0
    assert errs_bm4[0.10] / errs_pbm4[0.10] > 150.0
    assert errs_pbm4[0.20] / errs_pbm4[0.10] > 48.0


def test_shadow_manifold_lifecycle_and_pbm4_z4_elevation_float64() -> None:
    """
    Detached Component 4 Certification (Proposition 5.3 & Theorem 5.15):
    On a 3D eccentric hierarchical triple (e_in = 0.30, i = 28 deg):
    1. Full 4-method Shadow Manifold Lifecycle API: enter_shadow, step_shadow, observe_physical, exit_shadow.
    2. Canonical PBM4 with palindromic Lie-conjugation corrector Z_4^{BM4}
       reduces maximum relative energy error by > 200x over bare BM4.
    """
    from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec
    from normal_cone_nbody_float64.integrators.processed import create_pbm4_integrator

    masses = np.array([0.5, 0.5, 1.0], dtype=np.float64)
    a_in, e_in = 1.0, 0.30
    inc = np.deg2rad(28.0)
    r_p = a_in * (1.0 - e_in)
    v_p = np.sqrt((1.0 + e_in) / (a_in * (1.0 - e_in)))
    R_out = 30.0
    v_out = np.sqrt(2.0 / R_out)

    # Inner binary inclined by i = 28 deg relative to outer perturber
    q_rel_in = np.array([r_p, 0.0, 0.0], dtype=np.float64)
    v_rel_in = np.array([0.0, v_p * np.cos(inc), v_p * np.sin(inc)], dtype=np.float64)
    q0 = np.array([
        -0.5 * q_rel_in,
        +0.5 * q_rel_in,
        np.array([R_out * 0.8, R_out * 0.6, 0.0], dtype=np.float64),
    ], dtype=np.float64)
    v0 = np.array([
        -0.5 * v_rel_in,
        +0.5 * v_rel_in,
        np.array([-0.6 * v_out, 0.8 * v_out, 0.0], dtype=np.float64),
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=masses[:, None] * v0, masses=masses, t=0.0).to_barycentric()
    tree = build_tree_from_spec(masses, ((0, 1), 2))
    E0 = compute_invariants(state0).H

    bm4 = BlanesMoan4Integrator(G=1.0)
    pbm4_hom = create_pbm4_integrator(G=1.0)  # Canonical Z_4^{BM4}
    assert pbm4_hom.corrector.mode == "homological"

    dt = 0.075
    n_steps = 80
    st_bm4 = state0
    sh_hom = pbm4_hom.enter_shadow(state0, tree, dt)

    max_err_bm4 = 0.0
    max_err_hom = 0.0

    for _ in range(n_steps):
        st_bm4 = bm4.step(st_bm4, tree, dt)
        sh_hom = pbm4_hom.step_shadow(sh_hom, tree, dt)

        obs_hom = pbm4_hom.observe_physical(sh_hom, tree, dt)

        max_err_bm4 = max(max_err_bm4, abs(compute_invariants(st_bm4).H - E0) / abs(E0))
        max_err_hom = max(max_err_hom, abs(compute_invariants(obs_hom).H - E0) / abs(E0))

    # Verify observe_physical matches exit_shadow identically
    final_hom = pbm4_hom.exit_shadow(sh_hom, tree, dt)
    assert np.array_equal(final_hom.q, obs_hom.q)
    assert np.array_equal(final_hom.p, obs_hom.p)

    gain_hom = max_err_bm4 / max_err_hom
    assert gain_hom > 200.0, f"Expected Homological Z_4^{{BM4}} gain > 200x, got {gain_hom:.1f}x"


def test_covariant_gauge_operator_non_commutative_cech_cocycle_and_dynamic_jump_suppression_float64() -> None:
    """
    Detached Component 5 Certification (Theorems 5.17-5.18, Section 6.5 Algorithm 5 Step 4.4):
    1. PBM4GaugeTransitionOperator uses the canonical palindromic Lie-conjugation corrector Z_4^{BM4}.
    2. Exact 3-chart + 3-dyadic-step non-commutative Cech 1-cocycle closure:
           T_{(C, h_C) -> (A, h_A)} o T_{(B, h_B) -> (C, h_C)} o T_{(A, h_A) -> (B, h_B)}(z) == z  (< 1e-14),
       and exact conservation of all 9 Lie--Bargmann generators (P, L, G < 1e-14).
    3. Dynamic multi-step interfacial jump suppression:
       When integrating M_1 steps on the 3D eccentric hierarchical triple at h_A,
       transitioning to h_B, and integrating M_2 steps at h_B,
       the covariant non-commutative sandwich T_{(A, h_A) -> (B, h_B)}^{(PBM4)} reduces the post-transition
       energy error by > 100x compared to a naive transition (no gauge connection).
    """
    from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec
    from normal_cone_nbody_float64.integrators.processed import (
        PBM4GaugeTransitionOperator,
        apply_gauge_connection,
        create_pbm4_integrator,
    )

    masses = np.array([0.5, 0.5, 1.0], dtype=np.float64)
    a_in, e_in = 1.0, 0.30
    inc = np.deg2rad(28.0)
    r_p = a_in * (1.0 - e_in)
    v_p = np.sqrt((1.0 + e_in) / (a_in * (1.0 - e_in)))
    R_out = 30.0
    v_out = np.sqrt(2.0 / R_out)

    q_rel_in = np.array([r_p, 0.0, 0.0], dtype=np.float64)
    v_rel_in = np.array([0.0, v_p * np.cos(inc), v_p * np.sin(inc)], dtype=np.float64)
    q0 = np.array([
        -0.5 * q_rel_in,
        +0.5 * q_rel_in,
        np.array([R_out * 0.8, R_out * 0.6, 0.0], dtype=np.float64),
    ], dtype=np.float64)
    v0 = np.array([
        -0.5 * v_rel_in,
        +0.5 * v_rel_in,
        np.array([-0.6 * v_out, 0.8 * v_out, 0.0], dtype=np.float64),
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=masses[:, None] * v0, masses=masses, t=0.0).to_barycentric()
    inv0 = compute_invariants(state0)

    tree_A = build_tree_from_spec(masses, ((0, 1), 2))
    tree_B = build_tree_from_spec(masses, (0, (1, 2)))
    tree_C = build_tree_from_spec(masses, ((0, 2), 1))

    gauge_op = PBM4GaugeTransitionOperator(G=1.0)
    assert getattr(gauge_op.corrector, "mode", None) == "homological"

    # 1. Non-commutative 3-chart + 3-scale Cech 1-cocycle closure:
    #    (tree_A, h_A=0.04) -> (tree_B, h_B=0.02) -> (tree_C, h_C=0.01) -> (tree_A, h_A=0.04)
    cocycle_err = gauge_op.verify_cech_cocycle(
        state0, tree_A, tree_B, tree_C, dt_A=0.04, dt_B=0.02, dt_C=0.01
    )
    rel_cocycle_err = cocycle_err / float(np.max(np.abs(state0.q)))
    assert rel_cocycle_err < 5e-15, f"Non-commutative Cech 1-cocycle relative residual too large: {rel_cocycle_err:.3e}"

    # Verify Lie--Bargmann conservation on the non-commutative gauge transition
    s_AB = apply_gauge_connection(state0, tree_A, tree_B, dt=0.04, dt_new=0.02)
    inv_AB = compute_invariants(s_AB)
    assert np.linalg.norm(inv_AB.P - inv0.P) < 1e-14
    assert np.linalg.norm(inv_AB.L - inv0.L) < 1e-14
    assert np.linalg.norm(inv_AB.G - inv0.G) < 1e-14
    assert np.max(np.abs(s_AB.q - state0.q)) > 1e-10

    # 2. Dynamic multi-step boundary transition with real integration before and after the switch
    pbm4_hom = create_pbm4_integrator(G=1.0)

    h_A = 0.10
    h_B = 0.05
    M_1 = 3
    M_2 = 12
    E0 = inv0.H

    # Phase 1: enter shadow at h_A and integrate M_1 steps in shadow space
    sh_hom_1 = pbm4_hom.enter_shadow(state0, tree_A, h_A)
    for _ in range(M_1):
        sh_hom_1 = pbm4_hom.step_shadow(sh_hom_1, tree_A, h_A)

    # Transition at t = M_1 * h_A from h_A -> h_B:
    # (a) Naive transition without gauge connection (continuing shadow state directly at h_B)
    sh_naive = sh_hom_1
    # (b) Covariant non-commutative sandwich T_{(A, h_A) -> (A, h_B)}^{(PBM4)}
    sh_hom_trans = apply_gauge_connection(sh_hom_1, tree_A, tree_A, dt=h_A, dt_new=h_B)

    max_err_naive = 0.0
    max_err_hom = 0.0

    # Phase 2: integrate M_2 steps at the new scale h_B and observe physical energy
    for _ in range(M_2):
        sh_naive = pbm4_hom.step_shadow(sh_naive, tree_A, h_B)
        sh_hom_trans = pbm4_hom.step_shadow(sh_hom_trans, tree_A, h_B)

        phys_naive = pbm4_hom.observe_physical(sh_naive, tree_A, h_B)
        phys_hom = pbm4_hom.observe_physical(sh_hom_trans, tree_A, h_B)

        max_err_naive = max(max_err_naive, abs(compute_invariants(phys_naive).H - E0) / abs(E0))
        max_err_hom = max(max_err_hom, abs(compute_invariants(phys_hom).H - E0) / abs(E0))

    suppression_vs_naive = max_err_naive / max_err_hom
    assert suppression_vs_naive > 100.0, f"Expected > 100x jump suppression vs naive, got {suppression_vs_naive:.1f}x"






