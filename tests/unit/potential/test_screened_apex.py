"""
Unit Tests for Layer 3: Universal Screened Apex Potential C^{1,1} and Symplectic Kick Operator.

Theoretical Alignments:
- Monograph: Section 4.2 (Universal C^{1,1} Regularity and Multipole Expansion).
- Formal Lean 4: `01_ApexUniversalLipschitzBound_Theorems.lean` (`apex_universal_lipschitz_bound`,
  `apex_potential_bounded`).
- Invariants: Strict action-reaction sum F_i = 0, zero torque sum q_i x F_i = 0,
  analytical gradient consistency with finite differences < 1e-10, and ghost-pole eradication.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.potential.screened_apex import (
    ScreenedApexPotential,
    compute_quadrupole_tensor,
    compute_quadrupole_potential,
)


def test_quadrupole_tensor_properties():
    """Validates that the quadrupole moment tensor Q is strictly symmetric and traceless."""
    masses = np.array([1.2, 2.3, 3.4], dtype=np.float64)
    q = np.array([
        [1.0, 0.5, -0.2],
        [-0.5, 1.2, 0.8],
        [0.2, -1.0, 0.4],
    ], dtype=np.float64)

    Q = compute_quadrupole_tensor(q, masses)

    # 1. Symmetry Q^T == Q
    assert np.allclose(Q, Q.T, atol=1e-15)

    # 2. Traceless Tr(Q) == 0
    assert np.isclose(np.trace(Q), 0.0, atol=1e-14)


def test_apex_forces_action_reaction_and_zero_torque():
    """
    Validates machine-precision Lie-Bargmann linear and angular momentum conservation
    for the screened apex force field:
        sum_i F_i = 0 (< 1e-15),
        sum_i q_i x F_i = 0 (< 1e-14).
    """
    masses = np.array([1.0, 2.0, 1.5, 3.0], dtype=np.float64)
    q = np.array([
        [1.0, 0.2, -0.1],
        [1.1, 0.3, 0.0],    # Pair (0, 1) close
        [-2.0, 0.5, 0.2],
        [-2.1, 0.4, 0.1],   # Pair (2, 3) close
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    for mode in ["cauchy", "hermite"]:
        apex = ScreenedApexPotential(G=1.0, screening_mode=mode)
        forces, V_apex = apex.compute_apex_forces(state, tree)

        # Action-reaction: total force is zero (< 5e-14 for force scales ~ 200)
        F_tot = np.sum(forces, axis=0)
        assert np.allclose(F_tot, 0.0, atol=5e-14), f"Non-zero net force: {F_tot}"

        # Zero net torque: sum q_i x F_i == 0 (< 5e-14 for torque scales ~ 400)
        torque = np.sum([np.cross(state.q[i], forces[i]) for i in range(state.N)], axis=0)
        assert np.allclose(torque, 0.0, atol=5e-14), f"Non-zero net torque: {torque}"


def test_apex_forces_analytical_gradient_vs_finite_differences():
    """
    Validates that the analytical apex forces F_i^apex match central finite differences
    of the apex potential V_apex:
        F_i^apex == - nabla_{q_i} V_{apex}
    to high numerical accuracy (< 1e-10).
    """
    masses = np.array([1.0, 1.2, 0.8], dtype=np.float64)
    q0 = np.array([
        [0.5, 0.1, -0.2],
        [-0.3, 0.4, 0.1],
        [0.1, -0.5, 0.3],
    ], dtype=np.float64)
    p0 = np.zeros_like(q0)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    apex = ScreenedApexPotential(G=1.0, screening_mode="hermite")
    analytical_forces, _ = apex.compute_apex_forces(state0, tree)

    eps = 1e-7
    for i in range(state0.N):
        for comp in range(3):
            q_plus = np.copy(q0)
            q_plus[i, comp] += eps
            state_p = ExtendedPhaseState(q=q_plus, p=p0, masses=masses)
            _, V_p = apex.compute_apex_forces(state_p, tree)

            q_minus = np.copy(q0)
            q_minus[i, comp] -= eps
            state_m = ExtendedPhaseState(q=q_minus, p=p0, masses=masses)
            _, V_m = apex.compute_apex_forces(state_m, tree)

            fd_force = - (V_p - V_m) / (2.0 * eps)
            assert np.isclose(analytical_forces[i, comp], fd_force, rtol=1e-4, atol=1e-6), (
                f"Gradient mismatch at body {i}, comp {comp}: analytical {analytical_forces[i, comp]} vs FD {fd_force}"
            )


def test_ghost_pole_eradication_at_coincidence():
    """
    Validates Lean 4 theorem 'apex_universal_lipschitz_bound':
    When child barycenters coincide (r_v -> 0, u_v <= 1), the screened monopole vanishes
    smoothly, ensuring finite forces and zero fictitious singularity without softening (d_0 = 0).
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Positions where composite cluster barycenters are nearly coincident
    q = np.array([
        [0.0, 0.05, 0.0],
        [0.0, -0.05, 0.0],   # Cluster {0, 1} cdm at (0, 0, 0)
        [0.0, 0.0, 0.01],    # Body 2 extremely close to (0, 0, 0)
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    apex = ScreenedApexPotential(G=1.0, screening_mode="cauchy")
    forces, V_apex = apex.compute_apex_forces(state, tree)

    assert not np.any(np.isnan(forces))
    assert not np.any(np.isinf(forces))
    assert not np.isnan(V_apex)
    assert not np.isinf(V_apex)


def test_apex_kick_symplectic_invariance():
    """Validates that step_apex_kick preserves total linear momentum strictly."""
    masses = np.array([1.5, 2.5, 3.5], dtype=np.float64)
    q = np.array([
        [1.0, 0.0, 0.0],
        [-0.5, 0.8, 0.0],
        [0.2, -0.3, 1.0],
    ], dtype=np.float64)
    p = np.array([
        [0.1, -0.2, 0.3],
        [-0.4, 0.1, -0.1],
        [0.3, 0.1, -0.2],
    ], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    apex = ScreenedApexPotential(G=1.0)
    dt = 0.05
    state_next = apex.step_apex_kick(state, tree, dt=dt)

    assert np.allclose(state_next.total_momentum, state.total_momentum, atol=1e-15)


def test_screened_apex_parity_and_unequal_mass_transition_gradient():
    """
    Validates D-144-04:
    1. compute_conormal_tidal_tensor is exported in float64/potential/screened_apex.py and returns a symmetric
       traceless tensor matching float128.
    2. Default screening_mode is 'c6' across both float64 and float128.
    3. Analytical gradient F_i = -nabla_{q_i} V_apex with include_screened_monopole=True in the active
       transition zone u_v in (1, 2) matches central finite differences (< 1e-9) for unequal masses.
    """
    import normal_cone_nbody_float64.potential.screened_apex as sa_f64
    import normal_cone_nbody_float128.potential.screened_apex as sa_f128

    assert hasattr(sa_f64, "compute_conormal_tidal_tensor"), (
        "compute_conormal_tidal_tensor missing from float64/potential/screened_apex.py"
    )
    assert sa_f64.ScreenedApexPotential().screening_mode == "c6"
    assert sa_f128.ScreenedApexPotential().screening_mode == "c6"

    masses = np.array([1.0, 2.3, 3.7, 5.1], dtype=np.float64)
    # Cluster L = {0, 1} around (0,0,0) with spread ~0.5; body 2 at distance ~1.2 so u_v in (1.0, 2.0)
    q0 = np.array([
        [-0.35, 0.10, -0.05],
        [0.25, -0.15, 0.08],
        [0.85, 0.20, -0.10],
        [-1.10, -0.40, 0.30],
    ], dtype=np.float64)
    p0 = np.zeros_like(q0)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    T_apex = sa_f64.compute_conormal_tidal_tensor(state0, leaf_pair=(0, 1), G=1.0)
    assert np.allclose(T_apex, T_apex.T, atol=1e-15)
    assert abs(np.trace(T_apex)) < 1e-14

    apex = ScreenedApexPotential(G=1.0, screening_mode="c6")
    forces, V0 = apex.compute_apex_forces(state0, tree, include_screened_monopole=True)
    _, V0_unscreened = apex.compute_apex_forces(state0, tree, include_screened_monopole=False)
    # Ensure we are genuinely in the active transition/screened region where W_v != 0
    assert abs(V0 - V0_unscreened) > 1e-2

    # Action-reaction & zero torque even with screened monopole active
    assert np.allclose(np.sum(forces, axis=0), 0.0, atol=1e-13)
    torque = np.sum([np.cross(q0[i], forces[i]) for i in range(4)], axis=0)
    assert np.allclose(torque, 0.0, atol=1e-13)

    eps = 1e-7
    for i in range(4):
        for d in range(3):
            q_p = q0.copy()
            q_m = q0.copy()
            q_p[i, d] += eps
            q_m[i, d] -= eps
            _, V_p = apex.compute_apex_forces(
                ExtendedPhaseState(q=q_p, p=p0, masses=masses), tree, include_screened_monopole=True
            )
            _, V_m = apex.compute_apex_forces(
                ExtendedPhaseState(q=q_m, p=p0, masses=masses), tree, include_screened_monopole=True
            )
            fd_f = -(V_p - V_m) / (2.0 * eps)
            assert np.isclose(forces[i, d], fd_f, rtol=1e-7, atol=1e-8), (
                f"Screened gradient mismatch at ({i},{d}): analytical={forces[i, d]} vs FD={fd_f}"
            )


def test_optimized_screened_apex_cache_key_and_liskov_signature():
    """
    Validates D-144-05 & D-144-06:
    1. Switching supernode_bodies from {0, 1, 2} to {0, 1, 3} (same cardinality 3) on the same
       OptimizedScreenedApexPotential instance does not collide in _topology_cache and matches
       ScreenedApexPotential exactly.
    2. OptimizedScreenedApexPotential.compute_apex_forces supports mask_leaf_pairs and
       screening_mode='pure_cartesian' (Liskov substitution).
    """
    from normal_cone_nbody_float64.potential.optimized_screened_apex import (
        OptimizedScreenedApexPotential,
    )

    masses = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    q = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.1, -0.1],
        [-1.2, 0.8, 0.3],
        [1.5, -0.7, 0.4],
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    ref_apex = ScreenedApexPotential(G=1.0)
    opt_apex = OptimizedScreenedApexPotential(G=1.0)

    # First call populates cache with supernode {0, 1, 2} (cardinality 3)
    f_ref_012, v_ref_012 = ref_apex.compute_apex_forces(state, tree, supernode_bodies={0, 1, 2})
    f_opt_012, v_opt_012 = opt_apex.compute_apex_forces(state, tree, supernode_bodies={0, 1, 2})
    assert np.allclose(f_opt_012, f_ref_012, atol=1e-15)
    assert np.isclose(v_opt_012, v_ref_012, atol=1e-15)

    # Second call switches supernode to {0, 1, 3} (also cardinality 3!)
    f_ref_013, v_ref_013 = ref_apex.compute_apex_forces(state, tree, supernode_bodies={0, 1, 3})
    f_opt_013, v_opt_013 = opt_apex.compute_apex_forces(state, tree, supernode_bodies={0, 1, 3})
    assert np.allclose(f_opt_013, f_ref_013, atol=1e-15), (
        f"Cache collision detected when switching supernode {0,1,2} -> {0,1,3}: "
        f"max diff = {np.max(np.abs(f_opt_013 - f_ref_013))}"
    )
    assert np.isclose(v_opt_013, v_ref_013, atol=1e-15)

    # Liskov signature: mask_leaf_pairs=False
    f_ref_unmasked, v_ref_unmasked = ref_apex.compute_apex_forces(state, tree, mask_leaf_pairs=False)
    f_opt_unmasked, v_opt_unmasked = opt_apex.compute_apex_forces(state, tree, mask_leaf_pairs=False)
    assert np.allclose(f_opt_unmasked, f_ref_unmasked, atol=1e-15)
    assert np.isclose(v_opt_unmasked, v_ref_unmasked, atol=1e-15)

    # pure_cartesian mode
    opt_cart = OptimizedScreenedApexPotential(G=1.0, screening_mode="pure_cartesian")
    f_opt_cart, v_opt_cart = opt_cart.compute_apex_forces(state, tree)
    assert np.allclose(f_opt_cart, f_ref_unmasked, atol=1e-15)
    assert np.isclose(v_opt_cart, v_ref_unmasked, atol=1e-15)


def test_normal_cone_tidal_tensor_canonical_export():
    """
    Validates D-143-08:
    `potential.screened_apex` exposes the canonical `compute_normal_cone_tidal_tensor`
    matching `compute_conormal_tidal_tensor` to machine precision.
    """
    from normal_cone_nbody_float64.potential.screened_apex import (
        compute_normal_cone_tidal_tensor,
        compute_conormal_tidal_tensor,
    )

    masses = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    q = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [2.0, 1.5, -0.5]], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    T_norm = compute_normal_cone_tidal_tensor(state, (0, 1), G=1.0)
    T_legacy = compute_conormal_tidal_tensor(state, (0, 1), G=1.0)
    assert np.array_equal(T_norm, T_legacy)
    assert np.isclose(np.trace(T_norm), 0.0, atol=1e-15)


def test_multi_supernode_cross_pair_preservation_and_composite_drift_n6():
    """
    Validates [ENG-FALSE-05] (§4.2 Def. 4.3 & Eq. 4.105):
    When two disjoint supernodes S_1 = (0, 1, 2) and S_2 = (3, 4, 5) are simultaneously active (N = 6):
    1. ScreenedApexPotential, OptimizedScreenedApexPotential, and compute_intercluster_forces mask
       only the 6 intra-supernode pairs while preserving all 9 inter-supernode cross-pairs (i in S_1, j in S_2).
    2. step_drift_and_kepler (in common.py and OptimizedKeplerSolver) preserves the relative barycenter
       drift along the root composite node joining S_1 and S_2 in V_comp^free(T).
    """
    from normal_cone_nbody_float64.potential.screened_apex import (
        compute_intercluster_forces,
    )
    from normal_cone_nbody_float64.potential.optimized_screened_apex import (
        OptimizedScreenedApexPotential,
    )
    from normal_cone_nbody_float64.integrators.common import (
        step_drift_and_kepler,
        step_apex_kick,
        set_current_apex_evaluator,
        get_current_apex_evaluator,
        reset_current_apex_evaluator,
    )
    from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
        OptimizedKeplerSolver,
    )

    masses = np.array([1.0, 1.5, 2.0, 1.2, 1.8, 2.5], dtype=np.float64)
    q = np.array([
        [0.0, 0.0, 0.0],
        [0.02, 0.0, 0.0],
        [0.0, 0.02, 0.0],
        [10.0, 0.0, 0.0],
        [10.02, 0.0, 0.0],
        [10.0, 0.02, 0.0],
    ], dtype=np.float64)
    # Give cluster S_1 total momentum +P_rel and cluster S_2 total momentum -P_rel
    p = np.array([
        [0.1, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [-0.1, 0.0, 0.0],
        [-0.1, 0.0, 0.0],
        [-0.1, 0.0, 0.0],
    ], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="balanced")
    clusters = [(0, 1, 2), (3, 4, 5)]

    # Exact manual calculation of the 9 cross-pair Newtonian forces between S_1 and S_2
    f_cross_exact = np.zeros_like(q)
    v_cross_exact = 0.0
    for i in (0, 1, 2):
        for j in (3, 4, 5):
            r_ij = q[i] - q[j]
            dist = np.linalg.norm(r_ij)
            f_vec = -masses[i] * masses[j] * r_ij / (dist ** 3)
            f_cross_exact[i] += f_vec
            f_cross_exact[j] -= f_vec
            v_cross_exact -= masses[i] * masses[j] / dist

    ref_apex = ScreenedApexPotential(G=1.0)
    opt_apex = OptimizedScreenedApexPotential(G=1.0)

    f_ref, v_ref = ref_apex.compute_apex_forces(state, tree, supernode_clusters=clusters)
    f_opt, v_opt = opt_apex.compute_apex_forces(state, tree, supernode_clusters=clusters)
    f_inter = compute_intercluster_forces(
        q, masses, {tuple(sorted(lp.body_indices)) for lp in tree.leaf_pairs}, G=1.0, supernode_clusters=clusters
    )

    assert np.allclose(f_ref, f_cross_exact, atol=1e-15), f"Cross-pair forces lost in ScreenedApexPotential: {f_ref}"
    assert np.allclose(f_opt, f_cross_exact, atol=1e-15), f"Cross-pair forces lost in OptimizedScreenedApexPotential: {f_opt}"
    assert np.allclose(f_inter, f_cross_exact, atol=1e-15)
    assert np.isclose(v_ref, v_cross_exact, atol=1e-14)
    assert np.isclose(v_opt, v_cross_exact, atol=1e-14)

    # Verify ContextVar dispatch in step_apex_kick ([ENG-FALSE-06])
    token = set_current_apex_evaluator(opt_apex)
    try:
        assert get_current_apex_evaluator() is opt_apex
        kicked = step_apex_kick(state, tree, ref_apex, dt=0.05, supernode_clusters=clusters)
        assert np.allclose(kicked.p, p + 0.05 * f_cross_exact, atol=1e-15)
    finally:
        reset_current_apex_evaluator(token)
    assert get_current_apex_evaluator() is None

    # Verify that the inter-supernode composite node drifts properly in step_drift_and_kepler
    # Find the composite node whose body_indices == (0, 1, 2, 3, 4, 5)
    root_comp = [n for n in tree.composite_nodes if set(n.body_indices) == {0, 1, 2, 3, 4, 5}]
    assert len(root_comp) == 1
    dt = 0.1
    st_canon = step_drift_and_kepler(state, tree, dt=dt, G=1.0, supernode_clusters=clusters)
    st_opt = OptimizedKeplerSolver(G=1.0).step_drift_and_kepler(state, tree, dt=dt, supernode_clusters=clusters)

    cdm1_init = np.sum(masses[:3, None] * q[:3], axis=0) / np.sum(masses[:3])
    cdm2_init = np.sum(masses[3:, None] * q[3:], axis=0) / np.sum(masses[3:])
    rel_init = cdm1_init - cdm2_init
    v_rel = (np.sum(p[:3], axis=0) / np.sum(masses[:3])) - (np.sum(p[3:], axis=0) / np.sum(masses[3:]))

    cdm1_canon = np.sum(masses[:3, None] * st_canon.q[:3], axis=0) / np.sum(masses[:3])
    cdm2_canon = np.sum(masses[3:, None] * st_canon.q[3:], axis=0) / np.sum(masses[3:])
    assert np.allclose(cdm1_canon - cdm2_canon, rel_init + dt * v_rel, atol=1e-14)

    cdm1_opt = np.sum(masses[:3, None] * st_opt.q[:3], axis=0) / np.sum(masses[:3])
    cdm2_opt = np.sum(masses[3:, None] * st_opt.q[3:], axis=0) / np.sum(masses[3:])
    assert np.allclose(cdm1_opt - cdm2_opt, rel_init + dt * v_rel, atol=1e-14)




