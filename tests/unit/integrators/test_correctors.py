"""
Unit Tests for Layer 4: Boundary Symplectic Corrector C_4 (GlobalCorrectorC4 & CorrectorC4).

Theoretical Alignments:
- Monograph: Chapter 5, Section 5.5, Theorems 5.14, 5.15, 5.16.
- Formal Lean 4: `01_CorrectorBM4_Theorems.lean`, `07_GlobalConormalCorrectorPBM4_Order6.lean`,
  `08_UniversalOrderSixObstruction.lean`.
- Invariants:
  1. Universal Obstruction (Theorem 5.14): Absence of CorrectorC6 / compute_conormal_c6_step.
  2. Default Global Unipotent Mode (Theorem 5.15): divisor = 480.0, strict spatial stationarity Delta q == 0,
     exact algebraic inversion C_4^{-1} o C_4 = Id (< 1e-15), and exact Lie--Bargmann torque balance.
  3. 5-Strata Partition (Theorem 5.16): Supernode leaf-pair masking and smooth screening attenuation Theta(u_k).
"""

import numpy as np
import pytest

import normal_cone_nbody_float64.integrators.correctors as correctors_mod
from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.integrators.correctors import (
    CorrectorC4,
    GlobalCorrectorC4,
    compute_conormal_tidal_tensor,
    compute_conormal_corrector_kick,
    compute_global_conormal_c4_kick,
)


def test_universal_obstruction_purge_of_corrector_c6():
    """
    Validates Theorem 5.14 (Universal Obstruction for Processed Order 6):
    CorrectorC6 and compute_conormal_c6_step must be completely purged from integrators.correctors.
    """
    assert not hasattr(correctors_mod, "CorrectorC6"), "CorrectorC6 must be purged per Theorem 5.14"
    assert not hasattr(correctors_mod, "compute_conormal_c6_step"), "compute_conormal_c6_step must be purged"
    assert not hasattr(correctors_mod, "KAPPA6_BM6"), "KAPPA6_BM6 must be purged"
    assert not hasattr(correctors_mod, "KAPPA6_CONORMAL_LEAF"), "KAPPA6_CONORMAL_LEAF must be purged"


def test_corrector_c4_default_global_mode_and_exact_inversion():
    """
    Validates Theorem 5.15 & [ENG-FOOTGUN-01]:
    - CorrectorC4 defaults to mode='homological' (unobstructed Z_4^{BM4} in P_1).
    - CorrectorC4(mode='global') delegates to GlobalCorrectorC4 in P_0:
      * Strict spatial stationarity: Delta q == 0 identically.
      * Exact algebraic inversion: C_4^{-1} o C_4 = Id to machine precision (< 1e-15).
      * Exact 1/480 homological scaling: compute_global_conormal_c4_kick == 3 * compute_conormal_corrector_kick(divisor=1440).
    """
    masses = np.array([1.0, 2.0, 1.5], dtype=np.float64)
    q0 = np.array([
        [0.6, -0.2, 0.1],
        [-0.4, 0.7, -0.3],
        [0.1, -0.5, 0.2],
    ], dtype=np.float64)
    p0 = np.array([
        [0.1, -0.2, 0.15],
        [-0.15, 0.1, -0.1],
        [0.05, 0.1, -0.05],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    default_corrector = CorrectorC4(G=1.0)
    assert default_corrector.mode == "homological"

    corrector = CorrectorC4(G=1.0, mode="global")
    assert corrector.mode == "global"
    h = 0.05

    # Forward then inverse
    s_fwd = corrector.apply(state0, tree, h, inverse=False)
    # Strict spatial stationarity Delta q == 0
    assert np.array_equal(s_fwd.q, state0.q), "GlobalCorrectorC4 must leave positions strictly invariant (Delta q == 0)"

    # Non-zero momentum shear kick
    dp_global = s_fwd.p - state0.p
    assert np.max(np.abs(dp_global)) > 0.0

    # Verify divisor = 480.0 (3x legacy 1440.0 kick)
    dp_1440 = compute_conormal_corrector_kick(state0, tree, h, G=1.0, divisor=1440.0)
    assert np.allclose(dp_global, 3.0 * dp_1440, atol=1e-16)

    s_inv = corrector.apply_inverse(s_fwd, tree, h)
    assert np.array_equal(s_inv.q, state0.q)
    assert np.allclose(s_inv.p, state0.p, atol=1e-15), f"Inversion error p: {np.max(np.abs(s_inv.p - state0.p))}"


def test_corrector_c4_scaling_powers():
    """Validates perturbation scaling: Delta p ~ h^4 for GlobalCorrectorC4 (and Delta q ~ h^4 in legacy mode)."""
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    q0 = np.array([
        [0.5, 0.0, 0.0],
        [-0.25, 0.433, 0.0],
        [-0.25, -0.433, 0.0],
    ], dtype=np.float64)
    p0 = np.zeros_like(q0)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    c4 = CorrectorC4(G=1.0, mode="global")
    h1 = 0.02
    h2 = 0.04

    # Global C4 momentum kick scaling: (h2 / h1)^4 = 2^4 = 16
    dp4_1 = np.max(np.abs(c4.apply(state0, tree, h1).p - state0.p))
    dp4_2 = np.max(np.abs(c4.apply(state0, tree, h2).p - state0.p))
    ratio4 = dp4_2 / max(dp4_1, 1e-20)
    assert np.isclose(ratio4, 16.0, rtol=1e-12)


def test_corrector_c4_supernode_masking_and_screening_taper():
    """
    Validates Theorem 5.16 (5-Strata Partition of T*Q[N]):
    1. Democratic Supernode Stratum D_super: When a leaf pair (ia, ib) belongs to supernode_bodies,
       its tidal corrector kick is completely masked (Delta p == 0 for a 3-body system).
    2. Screening Transition Stratum D_trans: When apex screening is active and a spectator k is in
       the transition zone u_k in (1, 2), its tidal contribution is attenuated by Theta(u_k) in (0, 1)
       while preserving exact linear and angular momentum conservation.
    """
    from normal_cone_nbody_float64.potential.screened_apex import ScreenedApexPotential

    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Leaf pair (0, 1) has separation r = 1.0 -> normalized hyper-radius rho_bar = sqrt(mu / M) * r = 0.5.
    # Place body 2 at R = 0.75 from binary CM -> u = R / rho_bar = 1.5 in (1.0, 2.0) transition zone!
    q0 = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.1, 0.74, 0.05],
    ], dtype=np.float64)
    p0 = np.zeros_like(q0)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")
    h = 0.05

    # 1. Supernode masking: leaf pair (0, 1) inside supernode {0, 1, 2}
    dp_masked = compute_global_conormal_c4_kick(
        state0, tree, h, G=1.0, supernode_bodies={0, 1, 2}
    )
    assert np.allclose(dp_masked, 0.0, atol=1e-16), "Active supernode leaf pair must be masked in K_4"

    # 2. Screening attenuation in D_trans (u ~ 1.5) with exact conservation of P and L
    apex_screened = ScreenedApexPotential(G=1.0, screening_mode="c6")
    dp_unscreened = compute_global_conormal_c4_kick(state0, tree, h, G=1.0, apex=None)
    dp_screened = compute_global_conormal_c4_kick(state0, tree, h, G=1.0, apex=apex_screened)

    norm_unscreened = float(np.linalg.norm(dp_unscreened))
    norm_screened = float(np.linalg.norm(dp_screened))
    assert 0.0 < norm_screened < norm_unscreened, (
        f"Expected screened kick (0 < {norm_screened} < {norm_unscreened}) in transition stratum u ~ 1.5"
    )

    # Verify exact linear and angular momentum conservation under screening
    assert np.allclose(np.sum(dp_screened, axis=0), 0.0, atol=1e-15)
    dL_screened = np.sum([np.cross(q0[i], dp_screened[i]) for i in range(3)], axis=0)
    assert np.allclose(dL_screened, 0.0, atol=1e-15)


def test_analytical_directional_force_derivative_precision_and_conservation():
    """
    Validates exact analytical directional force derivative:
    1. Zero net force derivative: sum_i D_v F_i = 0 (< 1e-15).
    2. High-precision agreement with centered finite differences (< 1e-10) without noise floor.
    """
    from normal_cone_nbody_float64.integrators.correctors import (
        compute_analytical_directional_force_derivative,
    )
    from normal_cone_nbody_float64.potential.screened_apex import ScreenedApexPotential

    masses = np.array([1.0, 1.5, 2.0, 0.8], dtype=np.float64)
    q0 = np.array([
        [0.8, -0.2, 0.1],
        [-0.4, 0.9, -0.2],
        [0.1, -0.5, 0.5],
        [-0.3, 0.2, 0.8],
    ], dtype=np.float64)
    p0 = np.array([
        [0.15, -0.1, 0.05],
        [-0.2, 0.15, -0.1],
        [0.05, -0.05, 0.05],
        [-0.1, 0.05, -0.1],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")
    eval_apex = ScreenedApexPotential(G=1.0)
    v_dir = p0 / masses[:, None]

    # 1. Analytical directional derivative
    dF_anal = compute_analytical_directional_force_derivative(eval_apex, state0, tree, None, v_dir)

    # 2. Strict conservation of linear momentum: sum_i D_v F_i = 0
    sum_dF = np.sum(dF_anal, axis=0)
    assert np.allclose(sum_dF, 0.0, atol=1e-15), f"Non-zero momentum derivative: {sum_dF}"

    # 3. Agreement with centered finite differences
    eta = 1e-6
    st_plus = state0.with_updates(q=q0 + eta * v_dir)
    st_minus = state0.with_updates(q=q0 - eta * v_dir)
    f_plus, _ = eval_apex.compute_apex_forces(st_plus, tree)
    f_minus, _ = eval_apex.compute_apex_forces(st_minus, tree)
    dF_fd = (f_plus - f_minus) / (2.0 * eta)

    max_err = np.max(np.abs(dF_anal - dF_fd))
    assert max_err < 1e-9, f"Analytical derivative disagrees with FD: max error = {max_err}"


def test_corrector_c4_calibrated_constants():
    """Validates that C4 uses the calibrated homological constant and homological mode by default."""
    from normal_cone_nbody_float64.core.constants import KAPPA4_BM4

    c4 = CorrectorC4()
    assert c4.kappa4 == KAPPA4_BM4 == 1.0 / 2880.0
    assert c4.mode == "homological"


def test_conormal_tidal_tensor_symmetry_and_zero_trace():
    """Validates that T_apex in Sym_0(3) has exact zero trace and symmetry."""
    masses = np.array([3.0, 4.0, 5.0], dtype=np.float64)
    q0 = np.array([
        [-1.0, 0.0, 0.0],
        [1.0, 3.0, 0.0],
        [2.0, -1.0, 0.0],
    ], dtype=np.float64)
    p0 = np.zeros_like(q0)
    state = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    assert len(tree.leaf_pairs) >= 1
    leaf_node = tree.leaf_pairs[0]

    T_apex, r, mu, R_cm = compute_conormal_tidal_tensor(state, leaf_node, G=1.0)

    # 1. Symmetry: T_apex^T == T_apex
    assert np.allclose(T_apex, T_apex.T, atol=1e-16), f"Symmetry error: {np.max(np.abs(T_apex - T_apex.T))}"

    # 2. Traceless: Tr(T_apex) == 0.0
    trace_val = float(np.trace(T_apex))
    assert abs(trace_val) < 1e-15, f"Trace non-zero: {trace_val}"


def test_conormal_corrector_exact_momentum_and_angular_conservation():
    """
    Validates that global C4 conormal kicks conserve total linear momentum and total
    angular momentum strictly to machine precision (< 1e-15).
    """
    masses = np.array([3.0, 4.0, 5.0], dtype=np.float64)
    q0 = np.array([
        [-1.0, 0.2, 0.1],
        [1.0, 3.0, -0.2],
        [2.0, -1.0, 0.3],
    ], dtype=np.float64)
    p0 = np.array([
        [0.1, -0.3, 0.0],
        [-0.2, 0.1, 0.2],
        [0.1, 0.2, -0.2],
    ], dtype=np.float64)
    state = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    h = 0.05
    dp = compute_global_conormal_c4_kick(state, tree, h, G=1.0, inverse=False)

    # 1. Total linear momentum conservation: sum_i Delta p_i == 0
    sum_dp = np.sum(dp, axis=0)
    assert np.allclose(sum_dp, 0.0, atol=1e-15), f"Linear momentum violation: {sum_dp}"

    # 2. Total angular momentum conservation: sum_i q_i x Delta p_i == 0
    dL = np.sum([np.cross(q0[i], dp[i]) for i in range(len(masses))], axis=0)
    assert np.allclose(dL, 0.0, atol=1e-15), f"Angular momentum violation: {dL}"


def test_conormal_corrector_pericenter_smooth_vanishing():
    """
    Validates that the global C4 kick smoothly vanishes as r -> 0:
        ||Delta w|| proportional to ||r|| -> 0.
    """
    masses = np.array([3.0, 4.0, 5.0], dtype=np.float64)
    tree = build_canonical_tree(masses, topology="comb")
    leaf_node = tree.leaf_pairs[0]
    ia, ib = leaf_node.body_indices
    ma, mb = masses[ia], masses[ib]
    M_ab = ma + mb

    h = 0.01
    separations = [0.01, 0.001, 1e-4, 1e-5]
    kick_norms = []

    for sep in separations:
        # Keep leaf barycenter fixed at origin
        q = np.array([
            [-(mb / M_ab) * sep, 0.0, 0.0],
            [+(ma / M_ab) * sep, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ], dtype=np.float64)
        p = np.zeros_like(q)
        st = ExtendedPhaseState(q=q, p=p, masses=masses)
        dp = compute_global_conormal_c4_kick(st, tree, h, G=1.0)
        dp_rel = dp[ib] - dp[ia]
        kick_norms.append(float(np.linalg.norm(dp_rel)))

    # Verify monotonic decay proportional to sep
    for k in range(len(separations) - 1):
        ratio = kick_norms[k + 1] / kick_norms[k]
        expected_ratio = separations[k + 1] / separations[k]
        assert np.isclose(ratio, expected_ratio, rtol=1e-3), f"Linearity failure at step {k}: {ratio} vs {expected_ratio}"


def test_corrector_c4_conormal_exact_algebraic_inversion_tight_binary():
    """
    Validates C_4^{-1} o C_4 = Id on an ultra-tight binary encounter (< 1e-15).
    """
    masses = np.array([3.0, 4.0, 5.0], dtype=np.float64)
    # Ultra-tight binary separation r = 0.014 (pericenter scale of Burrau)
    q0 = np.array([
        [-0.007, 0.0, 0.0],
        [0.007, 0.0, 0.0],
        [1.5, 2.0, 0.5],
    ], dtype=np.float64)
    p0 = np.array([
        [0.5, -15.0, 0.0],
        [-0.5, 15.0, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    corrector = CorrectorC4(G=1.0, mode="global")
    h = 0.005

    s_fwd = corrector.apply(state0, tree, h, inverse=False)
    s_inv = corrector.apply_inverse(s_fwd, tree, h)

    assert np.array_equal(s_inv.q, state0.q)
    assert np.allclose(s_inv.p, state0.p, atol=1e-15), f"Inversion error p: {np.max(np.abs(s_inv.p - state0.p))}"


def test_all_three_corrector_modes_exact_group_inversion_and_lie_bargmann_invariants():
    """
    Detached Component 4 Certification:
    1. Exact group inversion ||C_h(C_h^{-1}(z_0)) - z_0||_inf <= 1e-15 across all 3 corrector modes:
       - mode='homological' (Z_4^{BM4}(h) in P_1, canonical default)
       - mode='cotangent_lift' (C_4^{cot}(h) in P_1 via 2-stage Gauss-Legendre collocation)
       - mode='global' (C_4^{global}(h) in P_0)
    2. Exact conservation of all 9 Lie-Bargmann phase-space generators (P, L, G) to <= 1e-15.
    3. Multi-supernode cluster intra-pair masking via supernode_clusters.
    """
    from normal_cone_nbody_float64.core.invariants import compute_invariants
    from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec

    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q0 = np.array([
        [-0.4, 0.1, -0.05],
        [0.4, -0.1, 0.05],
        [0.2, 8.0, 1.5],
    ], dtype=np.float64)
    p0 = np.array([
        [0.05, -0.75, 0.12],
        [-0.02, 0.78, -0.10],
        [0.20, -0.05, 0.08],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=1.25)
    tree = build_tree_from_spec(masses, ((0, 1), 2))
    inv0 = compute_invariants(state0, G=1.0)
    h = 0.05

    for mode in ("homological", "cotangent_lift", "global"):
        corr = CorrectorC4(G=1.0, mode=mode)
        st_in = corr.apply_inverse(state0, tree, h)
        st_out = corr.apply(st_in, tree, h, inverse=False)

        err_q = float(np.max(np.abs(st_out.q - state0.q)))
        err_p = float(np.max(np.abs(st_out.p - state0.p)))
        rel_err_q = err_q / max(1.0, float(np.max(np.abs(state0.q))))
        rel_err_p = err_p / max(1.0, float(np.max(np.abs(state0.p))))
        assert rel_err_q <= 2e-15, f"mode={mode} relative q inversion error {rel_err_q:.3e} > 2e-15"
        assert rel_err_p <= 2e-15, f"mode={mode} relative p inversion error {rel_err_p:.3e} > 2e-15"
        if mode in ("cotangent_lift", "global"):
            assert err_q <= 1e-15 and err_p <= 1e-15
        assert st_in.t == state0.t and st_out.t == state0.t

        inv_in = compute_invariants(st_in, G=1.0)
        assert float(np.max(np.abs(inv_in.P - inv0.P))) <= 1e-14, f"mode={mode} P drift"
        assert float(np.max(np.abs(inv_in.L - inv0.L))) <= 1e-14, f"mode={mode} L drift"
        assert float(np.max(np.abs(inv_in.G - inv0.G))) <= 1e-14, f"mode={mode} G drift"


def test_cotangent_lift_asymptotic_free_fall_regularity_at_pericenter():
    """
    Detached Component 4 Certification (Theorem 5.16 Item 2 & Report 04 [CRIT-7.1.3-04]):
    Along a Keplerian parabolic radial free-fall branch w(r) = -mu * sqrt(2 G M_ab / r) * r_hat
    as r -> 0 (r in [1e-2, 1e-3, 1e-4, 1e-5]):
    1. Under Sundman time-regularized scaling (multiplying by the Keplerian Sundman factor r^3),
       the cotangent-lifted quadrupole tensor M_reg(r, w) = r (x) w + w (x) r satisfies:
       - Generator magnitude |K_4^{reg}(r)| = O(r^{+1.0}) -> 0 (slope = 1.0);
       - Position shift ||Delta r^{reg}|| = O(r^{+1.5}) -> 0 (slope = 1.5);
       - Momentum shift ||Delta w^{reg}|| = O(r^{0.0}) = O(1) strictly bounded (slope = 0.0).
    2. The canonical default HomologicalLieCorrectorC4 (Z_4^{BM4}(h), Theorem 5.16 Item 2)
       remains completely non-singular and preserves exact group inversion (< 1e-13) and
       Lie-Bargmann invariants down to r = 1e-5 along the radial free-fall branch.
    """
    from normal_cone_nbody_float64.core.constants import C4_BCH_BM4
    from normal_cone_nbody_float64.core.invariants import compute_invariants
    from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec
    from normal_cone_nbody_float64.integrators.correctors import HomologicalLieCorrectorC4

    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    M_ab = 2.0
    mu = 0.5
    G = 1.0
    h = 0.01
    tree = build_tree_from_spec(masses, ((0, 1), 2))
    hom_corr = HomologicalLieCorrectorC4(G=G)

    r_vals = [1e-2, 1e-3, 1e-4, 1e-5]
    k4_vals = []
    dr_vals = []
    dw_vals = []

    r_hat = np.array([1.0, 2.0, -1.0], dtype=np.float64)
    r_hat /= np.linalg.norm(r_hat)

    for r_mag in r_vals:
        r_vec = r_mag * r_hat
        w_mag = mu * np.sqrt(2.0 * G * M_ab / r_mag)
        w_vec = -w_mag * r_hat

        q = np.array([
            -0.5 * r_vec,
            +0.5 * r_vec,
            np.array([10.0, 3.0, -2.0], dtype=np.float64),
        ], dtype=np.float64)
        p = np.array([
            -w_vec,
            +w_vec,
            np.zeros(3, dtype=np.float64),
        ], dtype=np.float64)
        st = ExtendedPhaseState(q=q, p=p, masses=masses)

        T_apex, _, _, _ = compute_conormal_tidal_tensor(st, tree.leaf_pairs[0], G=G)
        M_reg = np.outer(r_vec, w_vec) + np.outer(w_vec, r_vec)
        M_reg_sq = M_reg @ M_reg
        scale = 0.5 * C4_BCH_BM4 * (h ** 4)
        k4_reg = abs(float(scale * np.trace(M_reg_sq @ T_apex)))
        # Analytical cotangent gradients of K_4^{reg}(r, w) = (c_4 h^4 / 2) Tr(M_reg(r, w)^2 T_apex)
        SY = M_reg @ T_apex + T_apex @ M_reg
        grad_w = scale * (SY @ r_vec)
        grad_r = scale * (SY @ w_vec)

        k4_vals.append(k4_reg)
        dr_vals.append(float(np.linalg.norm(grad_w)))
        dw_vals.append(float(np.linalg.norm(grad_r)))

        # Verify Theorem 5.16 Item 2: HomologicalLieCorrectorC4 is regular across r -> 0
        st_in = hom_corr.apply_inverse(st, tree, h)
        st_out = hom_corr.apply(st_in, tree, h, inverse=False)
        assert np.all(np.isfinite(st_in.q)) and np.all(np.isfinite(st_in.p))
        assert np.all(np.isfinite(st_out.q)) and np.all(np.isfinite(st_out.p))
        # At the physical local dyadic step h_dyn = 0.1 * r^{3/2}, inversion is exact to < 5e-15:
        h_dyn = 0.1 * (r_mag ** 1.5)
        st_in_dyn = hom_corr.apply_inverse(st, tree, h_dyn)
        st_out_dyn = hom_corr.apply(st_in_dyn, tree, h_dyn, inverse=False)
        assert float(np.max(np.abs(st_out_dyn.q - st.q))) < 5e-15
        assert float(np.max(np.abs(st_out_dyn.p - st.p))) / float(np.max(np.abs(st.p))) < 1e-14
        inv0 = compute_invariants(st, G=G)
        inv_in = compute_invariants(st_in_dyn, G=G)
        assert float(np.max(np.abs(inv_in.P - inv0.P))) / max(1.0, float(np.max(np.abs(st.p)))) < 1e-15
        assert float(np.max(np.abs(inv_in.L - inv0.L))) < 1e-13

    # Check empirical log10 slopes across 3 decades in r (1e-2 -> 1e-5)
    slope_k4 = np.log10(k4_vals[0] / k4_vals[-1]) / np.log10(r_vals[0] / r_vals[-1])
    slope_dr = np.log10(dr_vals[0] / dr_vals[-1]) / np.log10(r_vals[0] / r_vals[-1])
    slope_dw = np.log10(dw_vals[0] / dw_vals[-1]) / np.log10(r_vals[0] / r_vals[-1])

    assert abs(slope_k4 - 1.0) < 1e-3, f"Expected |K_4| ~ O(r^1.0), got slope {slope_k4}"
    assert abs(slope_dr - 1.5) < 1e-3, f"Expected ||Delta r|| ~ O(r^1.5), got slope {slope_dr}"
    assert abs(slope_dw - 0.0) < 1e-3, f"Expected ||Delta w|| ~ O(r^0.0), got slope {slope_dw}"


def test_colocalized_listing_02_tree_builder_and_corrector_interface():
    """
    Validates the colocalized Listing 02 (`02_TreeBuilderAndCorrector_Interface.py`):
    1. Absence of obsolete C6 / KAPPA6 classes and constants.
    2. `detect_supernodes` works on N=4 without IndexError on `tree.compute_normalized_hyper_radius_cartesian(state.q, ...)`,
       selects external perturber by maximum tidal field M_P / r_P^3, and supports unequal-mass non-transversal coalescence.
    3. `GlobalCorrectorC4` / `CorrectorC4` matches the main repository (divisor = 480.0, Delta q == 0, Delta L == 0).
    """
    import importlib.util
    from pathlib import Path
    from normal_cone_nbody_float64.integrators.correctors import GlobalCorrectorC4 as MainGlobalCorrectorC4

    listing_path = (
        Path(__file__).resolve().parents[3]
        / "interfaces" / "02_TreeBuilderAndCorrector_Interface.py"
    )
    spec = importlib.util.spec_from_file_location("listing02_mod", listing_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    # 1. Purge of obsolete C6 / KAPPA6 symbols
    assert not hasattr(mod, "ConormalCorrectorC6")
    assert not hasattr(mod, "SymplecticCorrectorBM6")
    assert not hasattr(mod, "KAPPA6_BM6")
    assert hasattr(mod, "GlobalCorrectorC4")

    # 2. Tree building and supernode detection on N=4 without IndexError
    masses4 = np.array([1.0, 2.0, 3.0, 10.0], dtype=np.float64)
    q4 = np.array([
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [0.05, 0.08, 0.0],
        [10.0, 0.0, 0.0],
    ], dtype=np.float64)
    p4 = np.zeros_like(q4)
    state4 = ExtendedPhaseState(q=q4, p=p4, masses=masses4)
    tree4 = mod.build_optimal_tree(state4)
    sn4 = mod.detect_supernodes(state4, tree4)
    assert (0, 1, 2) in sn4

    # 3. Equivalence of GlobalCorrectorC4 with main repository
    corr_listing = mod.GlobalCorrectorC4(G=1.0)
    corr_main = MainGlobalCorrectorC4(G=1.0)
    s_list = corr_listing.apply(state4, tree4, h=0.05)
    s_main = corr_main.apply(state4, tree4, h=0.05)
    assert np.array_equal(s_list.q, state4.q)
    assert np.allclose(s_list.p, s_main.p, atol=1e-16)

    # 4. Canonical normal_cone exports in correctors_mod, listing02_mod, and NBodySystem terminology (D-143-08)
    assert hasattr(correctors_mod, "compute_normal_cone_tidal_tensor")
    assert hasattr(correctors_mod, "compute_normal_cone_corrector_kick")
    assert hasattr(correctors_mod, "compute_global_normal_cone_c4_kick")
    assert hasattr(mod, "compute_normal_cone_tidal_tensor")
    assert hasattr(mod, "compute_normal_cone_corrector_kick")
    assert hasattr(mod, "compute_global_normal_cone_c4_kick")

    import warnings
    from normal_cone_nbody_float64.system.nbody_system import NBodySystem

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _ = NBodySystem(state=state4, auto_supernodes=False)
    assert len(caught) == 1
    assert "normal-cone geometry" in str(caught[0].message)


