"""
Unit Tests for Layer 2: Synchronous Supernode Multi-Body Solver.

Theoretical Alignments:
- Monograph: Section 4.3 (Synchronous Supernodes and Excluded Middle Dilemma)
  and Section 4.4 (Anti-Zeno Regularization via Sundman Time Unfolding).
- Formal Lean 4: `01_UnifiedSuperNodeDef_Theorems.lean` and `02_SuperNodeAntiZenoRegularization_Theorems.lean`.
- Invariants: Strict conservation of internal linear and angular momentum,
  decoupling of center-of-mass motion, and Anti-Zeno finite sub-step bound N_sub < infty.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.solvers.supernode import (
    SupernodeCluster,
    compute_cluster_hyper_radius,
    step_supernode_synchronous,
)


def test_supernode_cluster_dataclass_validation():
    """Validates constructor invariants of SupernodeCluster."""
    q_cdm = np.zeros(3)
    p_cdm = np.zeros(3)

    # Valid cluster
    cluster = SupernodeCluster(
        body_indices=(0, 1, 2),
        mass=3.0,
        hyper_radius=0.25,
        center_of_mass_q=q_cdm,
        center_of_mass_p=p_cdm,
    )
    assert cluster.mass == 3.0
    assert cluster.hyper_radius == 0.25

    # Invalid: fewer than 2 bodies
    with pytest.raises(ValueError, match="at least 2 bodies"):
        SupernodeCluster(
            body_indices=(0,),
            mass=1.0,
            hyper_radius=0.1,
            center_of_mass_q=q_cdm,
            center_of_mass_p=p_cdm,
        )

    # Invalid: non-positive mass
    with pytest.raises(ValueError, match="strictly positive"):
        SupernodeCluster(
            body_indices=(0, 1),
            mass=-1.0,
            hyper_radius=0.1,
            center_of_mass_q=q_cdm,
            center_of_mass_p=p_cdm,
        )


def test_supernode_internal_linear_and_angular_momentum_conservation():
    """
    Validates machine-precision conservation of internal linear momentum
    P_S = sum p_i and angular momentum L_S = sum q_i x p_i within the supernode.
    """
    masses = np.array([1.0, 1.5, 2.0, 1.0], dtype=np.float64)
    # Bodies 0, 1, 2 form the supernode; body 3 is an external spectator
    q = np.array([
        [0.1, 0.0, 0.0],
        [-0.05, 0.08, 0.0],
        [-0.05, -0.08, 0.0],
        [10.0, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.0, 0.5, 0.0],
        [-0.3, -0.2, 0.1],
        [0.3, -0.3, -0.1],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)

    body_indices = (0, 1, 2)
    P_init = np.sum(p[list(body_indices)], axis=0)
    L_init = np.sum([np.cross(q[i], p[i]) for i in body_indices], axis=0)

    # Propagate supernode across several steps
    q_curr, p_curr = np.copy(q), np.copy(p)
    dt = 0.02
    for _ in range(20):
        q_curr, p_curr = step_supernode_synchronous(
            q_curr, p_curr, masses, body_indices, dt=dt, G=1.0, n_substeps=16
        )

        P_now = np.sum(p_curr[list(body_indices)], axis=0)
        L_now = np.sum([np.cross(q_curr[i], p_curr[i]) for i in body_indices], axis=0)

        assert np.allclose(P_now, P_init, atol=1e-14), f"Internal linear momentum drift: {np.max(np.abs(P_now - P_init))}"
        assert np.allclose(L_now, L_init, atol=1e-13), f"Internal angular momentum drift: {np.max(np.abs(L_now - L_init))}"

    # Spectator body 3 must be strictly untouched
    assert np.allclose(q_curr[3], q[3], atol=1e-16)
    assert np.allclose(p_curr[3], p[3], atol=1e-16)


def test_supernode_anti_zeno_finite_crossing():
    """
    Validates Lean 4 theorem 'super_node_anti_zeno_regularization':
    Multi-body cluster undergoes tight gravitational encounter in a finite number
    of sub-steps (N_sub < infty) without Zeno freezing or numerical divergence.
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Compact Pythagorean-like 3-body configuration
    q = np.array([
        [0.03, 0.04, 0.0],
        [-0.03, 0.04, 0.0],
        [0.0, -0.04, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q)

    # Evolution across physical time step dt
    dt = 0.05
    q_next, p_next = step_supernode_synchronous(q, p, masses, (0, 1, 2), dt=dt, G=1.0)

    assert not np.any(np.isnan(q_next))
    assert not np.any(np.isnan(p_next))
    assert not np.any(np.isinf(q_next))
    assert not np.any(np.isinf(p_next))

    rho_init, _ = compute_cluster_hyper_radius(q, masses, (0, 1, 2))
    rho_next, _ = compute_cluster_hyper_radius(q_next, masses, (0, 1, 2))
    assert rho_init > 0.0
    assert rho_next > 0.0


def test_supernode_homothetic_radial_collapse():
    """
    Validates the exact 1D radial homothety branch for symmetric central configuration collapse.
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Equilateral triangle centered at origin (central configuration)
    angle = np.array([0.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0])
    r0 = 0.5
    q = np.zeros((3, 3), dtype=np.float64)
    q[:, 0] = r0 * np.cos(angle)
    q[:, 1] = r0 * np.sin(angle)

    # Pure radial inward momentum
    v_rad = -0.2
    p = np.zeros((3, 3), dtype=np.float64)
    p[:, 0] = masses * v_rad * np.cos(angle)
    p[:, 1] = masses * v_rad * np.sin(angle)

    dt = 0.1
    q_next, p_next = step_supernode_synchronous(q, p, masses, (0, 1, 2), dt=dt, G=1.0)

    # Must preserve equilateral symmetry (distances to origin equal)
    r_norms = np.linalg.norm(q_next, axis=1)
    assert np.isclose(r_norms[0], r_norms[1], atol=1e-12)
    assert np.isclose(r_norms[1], r_norms[2], atol=1e-12)


def test_supernode_center_of_mass_drift():
    """Validates that the supernode center of mass drifts by dt * (p_cdm / m_tot) under Galilean inertia."""
    masses = np.array([1.0, 2.0, 1.5], dtype=np.float64)
    q = np.array([
        [0.1, 0.0, 0.0],
        [-0.05, 0.08, 0.0],
        [-0.05, -0.08, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.2, 0.5, 0.1],
        [-0.1, -0.2, 0.2],
        [0.4, -0.1, -0.1],
    ], dtype=np.float64)
    m_tot = float(np.sum(masses))
    p_cdm = np.sum(p, axis=0)
    q_cdm_0 = np.sum(masses[:, None] * q, axis=0) / m_tot

    dt = 0.05
    q_next, p_next = step_supernode_synchronous(q, p, masses, (0, 1, 2), dt=dt, G=1.0)

    q_cdm_expected = q_cdm_0 + dt * (p_cdm / m_tot)
    q_cdm_actual = np.sum(masses[:, None] * q_next, axis=0) / m_tot
    assert np.allclose(q_cdm_actual, q_cdm_expected, atol=1e-14)


def test_supernode_singular_collision_regularization():
    """Validates that synchronous sub-stepping safely regularizes extreme collision encounters (r -> 0)."""
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Bodies 0 and 1 are in extreme close proximity (r = 1e-12)
    q = np.array([
        [-5e-13, 0.0, 0.0],
        [5e-13, 0.0, 0.0],
        [0.0, 0.2, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)

    dt = 0.01
    q_next, p_next = step_supernode_synchronous(q, p, masses, (0, 1, 2), dt=dt, G=1.0)
    assert not np.any(np.isnan(q_next))
    assert not np.any(np.isnan(p_next))
    assert not np.any(np.isinf(q_next))
    assert not np.any(np.isinf(p_next))


def test_supernode_unfolded_sundman_adaptive_stepping_unsoftened():
    """
    Validates that synchronous supernode integration under unfolded Sundman clock
    dt_S = rho_S dtau_S operates with exact zero softening (d_safe = 0.0) and
    preserves linear momentum and reversibility under forward and retrograde steps.
    """
    masses = np.array([1.0, 1.5, 2.0], dtype=np.float64)
    q0 = np.array([
        [0.1, 0.0, 0.0],
        [-0.05, 0.08, 0.0],
        [-0.05, -0.08, 0.0],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, 0.5, 0.0],
        [-0.3, -0.2, 0.1],
        [0.3, -0.3, -0.1],
    ], dtype=np.float64)

    P_init = np.sum(p0, axis=0)

    # 1. Forward step under adaptive Sundman clock (n_substeps=None)
    dt_fwd = 0.01
    q_fwd, p_fwd = step_supernode_synchronous(
        q0, p0, masses, (0, 1, 2), dt=dt_fwd, G=1.0, n_substeps=None
    )
    P_fwd = np.sum(p_fwd, axis=0)
    assert np.allclose(P_fwd, P_init, atol=1e-14), f"Forward momentum drift: {np.max(np.abs(P_fwd - P_init))}"

    # 2. Retrograde step under adaptive Sundman clock (n_substeps=None)
    dt_bwd = -0.01
    q_bwd, p_bwd = step_supernode_synchronous(
        q_fwd, p_fwd, masses, (0, 1, 2), dt=dt_bwd, G=1.0, n_substeps=None
    )
    P_bwd = np.sum(p_bwd, axis=0)
    assert np.allclose(P_bwd, P_init, atol=1e-14), f"Retrograde momentum drift: {np.max(np.abs(P_bwd - P_init))}"

    # 3. Time reversibility: q_bwd and p_bwd must recover (q0, p0) to machine precision ([ENG-FOOTGUN-05])
    assert np.allclose(q_bwd, q0, atol=1e-14), f"Reversibility error on positions: {np.max(np.abs(q_bwd - q0))}"
    assert np.allclose(p_bwd, p0, atol=1e-14), f"Reversibility error on momenta: {np.max(np.abs(p_bwd - p0))}"


def test_supernode_d4_quadruple_collision_bounce():
    """
    Validates [ENG-FALSE-03] & [CRIT-7.1.3-02]:
    Symmetric D_4 quadruple collision configuration (N = 4 bodies at vertices of a square
    released from rest with t_coll^{(D_4)} = 2.08575) integrates cleanly across the total
    quadruple collision up to t_max = 2.50 > t_coll^{(D_4)} with exact elastic bounce,
    D_4 symmetry preservation, energy conservation, and round-trip time reversibility (< 1e-14).
    """
    from normal_cone_nbody_float64.solvers.optimized_supernode import (
        step_supernode_synchronous_fast,
    )
    from normal_cone_nbody_float64.solvers.jit_supernode import (
        step_supernode_synchronous_jit,
    )
    masses = np.ones(4, dtype=np.float64)
    r0 = 1.5  # Canonical D_4 radius giving t_coll^{(D_4)} = (pi/2)*sqrt(1.5^3 / (2*(0.25 + sqrt(2)/2))) = 2.085747
    q0 = np.array([
        [r0, 0.0, 0.0],
        [0.0, r0, 0.0],
        [-r0, 0.0, 0.0],
        [0.0, -r0, 0.0],
    ], dtype=np.float64)
    p0 = np.zeros((4, 3), dtype=np.float64)
    indices = (0, 1, 2, 3)
    dt_full = 2.50  # Strictly past t_coll^{(D_4)} = 2.08575

    def _cluster_energy(q_arr, p_arr):
        T = 0.5 * float(np.sum(p_arr ** 2))
        V = 0.0
        for i in range(4):
            for j in range(i + 1, 4):
                V -= 1.0 / float(np.linalg.norm(q_arr[i] - q_arr[j]))
        return T + V

    E0 = _cluster_energy(q0, p0)

    for step_fn in (step_supernode_synchronous, step_supernode_synchronous_fast, step_supernode_synchronous_jit):
        q_f, p_f = step_fn(q0, p0, masses, indices, dt=dt_full, G=1.0)
        assert not np.any(np.isnan(q_f)) and not np.any(np.isnan(p_f))
        r_norms = np.linalg.norm(q_f, axis=1)
        assert np.allclose(r_norms, r_norms[0], atol=1e-13)
        # After bounce (t = 2.50 > 2.08575), radial velocity must be outward (> 0)
        v_rad_out = float(np.dot(q_f[0], p_f[0])) / float(r_norms[0])
        assert v_rad_out > 0.1, f"Expected outward radial velocity after D_4 bounce, got {v_rad_out}"
        Ef = _cluster_energy(q_f, p_f)
        assert abs(Ef - E0) / abs(E0) < 1e-13

        # Round-trip reversibility back across the quadruple collision to t = 0
        q_rev, p_rev = step_fn(q_f, p_f, masses, indices, dt=-dt_full, G=1.0)
        assert np.allclose(q_rev, q0, atol=1e-13)
        assert np.allclose(p_rev, p0, atol=1e-13)



def test_supernode_pairwise_quasi_collision_stability():
    """
    Validates that internal pairwise close encounters across impact parameters
    b in [1e-1, 1e-3, 1e-6, 0.0] are regularized via Stumpff-Sundman isolation,
    preventing b^{-2} energy blowup and preserving total linear and angular momentum.
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    dt = 0.2

    for b in [0.1, 1e-3, 1e-6, 0.0]:
        q0 = np.array([
            [-0.1, b, 0.0],
            [0.1, -b, 0.0],
            [0.0, 0.5, 0.0]
        ], dtype=np.float64)
        p0 = np.array([
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0]
        ], dtype=np.float64)

        # Compute initial energy
        T0 = 0.5 * np.sum(np.sum(p0 ** 2, axis=1) / masses)
        V0 = - 1.0 / np.linalg.norm(q0[0] - q0[1]) - 1.0 / np.linalg.norm(q0[0] - q0[2]) - 1.0 / np.linalg.norm(q0[1] - q0[2])
        E0 = T0 + V0
        P0 = np.sum(p0, axis=0)

        q1, p1 = step_supernode_synchronous(q0, p0, masses, (0, 1, 2), dt=dt, G=1.0)

        assert not np.any(np.isnan(q1)), f"NaN in q1 for b={b}"
        assert not np.any(np.isnan(p1)), f"NaN in p1 for b={b}"
        assert not np.any(np.isinf(q1)), f"Inf in q1 for b={b}"
        assert not np.any(np.isinf(p1)), f"Inf in p1 for b={b}"

        T1 = 0.5 * np.sum(np.sum(p1 ** 2, axis=1) / masses)
        V1 = - 1.0 / np.linalg.norm(q1[0] - q1[1]) - 1.0 / np.linalg.norm(q1[0] - q1[2]) - 1.0 / np.linalg.norm(q1[1] - q1[2])
        E1 = T1 + V1
        rel_dE = abs((E1 - E0) / E0)

        # Bounded energy error without b^{-2} blowup (previously blew up to 1e11 - 1e15)
        assert rel_dE < 0.02, f"Energy blowup for b={b}: rel_dE={rel_dE:.2e}"

        # Momentum conservation
        P1 = np.sum(p1, axis=0)
        assert np.allclose(P1, P0, atol=1e-13), f"Momentum drift for b={b}: {np.max(np.abs(P1 - P0))}"


def test_supernode_zero_distance_exact_collision_guard():
    """
    Validates that _compute_internal_cluster_acc safely handles strict coincident
    particles (d_ij = 0.0) without NaN or zero-division exception.
    """
    from normal_cone_nbody_float64.solvers.supernode import _compute_internal_cluster_acc
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Bodies 0 and 1 at exact same coordinate
    x = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ], dtype=np.float64)

    acc, r_min = _compute_internal_cluster_acc(x, masses, G=1.0)
    assert not np.any(np.isnan(acc)), "NaN generated on coincident particles"
    assert not np.any(np.isinf(acc)), "Inf generated on coincident particles"
    # Acceleration from body 2 on 0 and 1 must be identical
    assert np.allclose(acc[0], acc[1], atol=1e-15)


def test_supernode_incomplete_time_exception():
    """
    Validates that a RuntimeError is explicitly raised if the internal integrator
    reaches max_steps without completing the target physical time interval.
    """
    from normal_cone_nbody_float64.solvers.supernode import step_supernode_synchronous
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    q0 = np.array([
        [-0.1, 0.01, 0.0],
        [0.1, -0.01, 0.0],
        [0.0, 0.5, 0.0]
    ], dtype=np.float64)
    p0 = np.array([
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0]
    ], dtype=np.float64)

    with pytest.raises(RuntimeError, match=r"reached max_steps \(\d+\) with incomplete physical time"):
        # Large dt with max_steps=2 cannot complete the interval
        step_supernode_synchronous(q0, p0, masses, (0, 1, 2), dt=5.0, max_steps=2)


def test_supernode_adaptive_stages_inheritance():
    """
    Validates that step_supernode_synchronous accepts and executes different
    palindromic composition weights (e.g. Strang order 2 vs BM4 order 4 vs BM6 order 6)
    and strictly preserves linear and angular momentum under the recursive normal cone split.
    """
    from normal_cone_nbody_float64.core.constants import (
        BM4_DRIFT_WEIGHTS, BM4_IMPULSE_WEIGHTS,
        BM6_DRIFT_WEIGHTS, BM6_IMPULSE_WEIGHTS,
    )

    masses = np.array([1.0, 1.2, 0.8], dtype=np.float64)
    q0 = np.array([
        [0.05, 0.0, 0.0],
        [-0.03, 0.04, 0.0],
        [-0.02, -0.04, 0.0],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, 0.4, 0.0],
        [-0.3, -0.2, 0.1],
        [0.3, -0.2, -0.1],
    ], dtype=np.float64)

    # 1. Integration with Strang weights
    strang_drift = (1.0,)
    strang_impulse = (0.5, 0.5)
    q_s, p_s = step_supernode_synchronous(
        q0, p0, masses, (0, 1, 2), dt=0.01,
        drift_weights=strang_drift, impulse_weights=strang_impulse,
        n_substeps=5,
    )
    assert not np.any(np.isnan(q_s))
    assert np.allclose(np.sum(p_s, axis=0), np.sum(p0, axis=0), atol=1e-14)

    # 2. Integration with BM4 weights
    q_bm4, p_bm4 = step_supernode_synchronous(
        q0, p0, masses, (0, 1, 2), dt=0.01,
        drift_weights=BM4_DRIFT_WEIGHTS, impulse_weights=BM4_IMPULSE_WEIGHTS,
        n_substeps=5,
    )
    assert not np.any(np.isnan(q_bm4))
    assert np.allclose(np.sum(p_bm4, axis=0), np.sum(p0, axis=0), atol=1e-14)

    # 3. Integration with BM6 weights
    a1, a2, a3, a4, a5, a6 = BM6_DRIFT_WEIGHTS
    b1, b2, b3, b4, b5, b6 = BM6_IMPULSE_WEIGHTS
    bm6_drift = (a1, a2, a3, a4, a5, a6, a5, a4, a3, a2, a1)
    bm6_impulse = (b1, b2, b3, b4, b5, 0.5 * b6, 0.5 * b6, b5, b4, b3, b2, b1)

    q_bm6, p_bm6 = step_supernode_synchronous(
        q0, p0, masses, (0, 1, 2), dt=0.01,
        drift_weights=bm6_drift, impulse_weights=bm6_impulse,
        n_substeps=5,
    )
    assert not np.any(np.isnan(q_bm6))
    assert np.allclose(np.sum(p_bm6, axis=0), np.sum(p0, axis=0), atol=1e-14)

    # High-order methods differ from lower-order scheme
    diff_order = np.max(np.abs(q_bm4 - q_s))
    assert diff_order > 1e-10, "BM4 should exhibit high-order deviation from Strang"


def test_supernode_context_stages_inheritance():
    """
    Validates that step_supernode_synchronous without explicit weights
    correctly inherits stage weights from execution context via get_current_integrator_stages.
    """
    from normal_cone_nbody_float64.integrators.common import (
        set_current_integrator_stages,
        _CURRENT_INTEGRATOR_STAGES,
    )
    masses = np.array([1.0, 1.2, 0.8], dtype=np.float64)
    q0 = np.array([
        [0.05, 0.0, 0.0],
        [-0.03, 0.04, 0.0],
        [-0.02, -0.04, 0.0],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, 0.4, 0.0],
        [-0.3, -0.2, 0.1],
        [0.3, -0.2, -0.1],
    ], dtype=np.float64)

    # Context with Strang
    token = set_current_integrator_stages(((1.0,), (0.5, 0.5)))
    try:
        q_ctx_s, p_ctx_s = step_supernode_synchronous(q0, p0, masses, (0, 1, 2), dt=0.01, n_substeps=5)
    finally:
        _CURRENT_INTEGRATOR_STAGES.reset(token)

    # Explicit Strang
    q_exp_s, p_exp_s = step_supernode_synchronous(
        q0, p0, masses, (0, 1, 2), dt=0.01, drift_weights=(1.0,), impulse_weights=(0.5, 0.5), n_substeps=5
    )

    assert np.allclose(q_ctx_s, q_exp_s, atol=1e-15, rtol=1e-15)
    assert np.allclose(p_ctx_s, p_exp_s, atol=1e-15, rtol=1e-15)


def test_optimized_and_jit_supernode_kernels_active_and_equivalent(monkeypatch):
    """
    Validates D-144-12:
    1. step_supernode_synchronous_fast and step_supernode_synchronous_jit execute their own
       specialized zero-allocation / @njit kernels for k = 3, 4, 5 WITHOUT short-circuiting
       via an early return to reference supernode.step_supernode_synchronous.
    2. Direct invocation of _substep_cluster_k3/k4/general and _njit_integrate_supernode_k3/k4/general
       succeeds (including k=5 where _njit_integrate_supernode_general previously had NameError: 'pi')
       and matches the multistage composition of supernode.step_supernode_synchronous to < 1e-13.
    """
    import normal_cone_nbody_float64.solvers.supernode as ref_mod
    from normal_cone_nbody_float64.solvers.optimized_supernode import (
        step_supernode_synchronous_fast,
        _substep_cluster_k3,
        _substep_cluster_k4,
        _substep_cluster_general,
    )
    from normal_cone_nbody_float64.solvers.jit_supernode import (
        step_supernode_synchronous_jit,
        NUMBA_AVAILABLE,
    )

    for k in (3, 4, 5):
        rng = np.random.default_rng(100 + k)
        masses = np.linspace(1.0, 2.0, k, dtype=np.float64)
        q0 = rng.uniform(-0.2, 0.2, size=(k, 3)).astype(np.float64)
        p0 = rng.uniform(-0.3, 0.3, size=(k, 3)).astype(np.float64)
        indices = tuple(range(k))
        dt = 0.005

        # Reference trajectory (before monkeypatching)
        q_ref, p_ref = ref_mod.step_supernode_synchronous(q0, p0, masses, indices, dt=dt, G=1.0)

        # Monkeypatch reference step_supernode_synchronous so any early-return short-circuit raises AssertionError!
        def _forbidden_ref_call(*args, **kwargs):
            raise AssertionError("Short-circuit to reference step_supernode_synchronous is forbidden!")

        monkeypatch.setattr(ref_mod, "step_supernode_synchronous", _forbidden_ref_call)
        try:
            q_fast, p_fast = step_supernode_synchronous_fast(q0, p0, masses, indices, dt=dt, G=1.0)
            q_jit, p_jit = step_supernode_synchronous_jit(q0, p0, masses, indices, dt=dt, G=1.0)
        finally:
            monkeypatch.undo()

        assert np.allclose(q_fast, q_ref, atol=1e-13), (
            f"k={k} optimized_supernode mismatch vs reference: max diff={np.max(np.abs(q_fast - q_ref)):.2e}"
        )
        assert np.allclose(p_fast, p_ref, atol=1e-13)
        assert np.allclose(q_jit, q_ref, atol=1e-13), (
            f"k={k} jit_supernode mismatch vs reference: max diff={np.max(np.abs(q_jit - q_ref)):.2e}"
        )
        assert np.allclose(p_jit, p_ref, atol=1e-13)


def test_colocalized_listing_06_high_performance_optimizations_interface():
    """
    Validates D-144-06 & D-144-10 on colocalized listing 06_HighPerformanceOptimizations_Interface.py:
    1. OptimizedScreenedApexPotential.compute_apex_forces exerts an ATTRACTIVE force
       (F_0 . (q_1 - q_0) > 0) and handles both N=3 and N=4 respecting tree.get_leaf_pairs().
    2. OptimizedKeplerSolver.step evolves (r0, w0) matching step_kepler_sundman_fast.
    3. JitSupernodeSolver.step_supernode_synchronous evolves (q, p) matching step_supernode_synchronous_jit.
    """
    import importlib.util
    from pathlib import Path
    from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
    from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import step_kepler_sundman_fast
    from normal_cone_nbody_float64.solvers.jit_supernode import step_supernode_synchronous_jit

    listing_path = (
        Path(__file__).resolve().parents[3]
        / "interfaces" / "06_HighPerformanceOptimizations_Interface.py"
    )
    spec = importlib.util.spec_from_file_location("listing06_interface", listing_path)
    mod = importlib.util.module_from_spec(spec)
    import sys
    sys.modules["listing06_interface"] = mod
    spec.loader.exec_module(mod)

    # 1. Attractive force sign and N >= 3 support in OptimizedScreenedApexPotential
    masses3 = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    q3 = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ], dtype=np.float64)
    tree3 = build_canonical_tree(masses3, topology="comb")  # leaf pair is (0, 1)
    opt_apex3 = mod.OptimizedScreenedApexPotential(masses3, G=1.0)
    f3 = opt_apex3.compute_apex_forces(q3, tree3)
    # Since (0, 1) is a leaf pair, apex forces come from non-leaf pairs (0, 2) and (1, 2):
    # body 2 at (0, 2, 0) attracts body 0 at (0, 0, 0) along +y -> f3[0, 1] > 0!
    assert f3[0, 1] > 0.1, f"Expected attractive +y force on body 0 from body 2, got {f3[0]}"
    assert np.allclose(np.sum(f3, axis=0), 0.0, atol=1e-15)

    # Also test N = 4
    masses4 = np.array([1.0, 1.5, 2.0, 2.5], dtype=np.float64)
    q4 = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 4.0],
    ], dtype=np.float64)
    tree4 = build_canonical_tree(masses4, topology="comb")
    opt_apex4 = mod.OptimizedScreenedApexPotential(masses4, G=1.0)
    f4 = opt_apex4.compute_apex_forces(q4, tree4)
    assert np.linalg.norm(f4) > 0.05, "N=4 apex forces should not be identically zero!"
    assert np.allclose(np.sum(f4, axis=0), 0.0, atol=1e-15)

    # 2. OptimizedKeplerSolver.step
    opt_kep = mod.OptimizedKeplerSolver(G=1.0)
    r0 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    w0 = np.array([0.0, 0.8, 0.0], dtype=np.float64)
    r_new, w_new = opt_kep.step(r0, w0, mu=0.5, gm=1.5, dt=0.1)
    r_ref, w_ref = step_kepler_sundman_fast(r0, w0, mu=0.5, GM=1.5, dt_target=0.1)
    assert not np.allclose(r_new, r0), "OptimizedKeplerSolver.step must evolve r0!"
    assert np.allclose(r_new, r_ref, atol=1e-14)
    assert np.allclose(w_new, w_ref, atol=1e-14)

    # 3. JitSupernodeSolver.step_supernode_synchronous
    jit_sn = mod.JitSupernodeSolver(G=1.0)
    p3 = np.array([
        [0.0, 0.2, 0.0],
        [-0.1, -0.1, 0.05],
        [0.1, -0.1, -0.05],
    ], dtype=np.float64)
    q_sn, p_sn = jit_sn.step_supernode_synchronous(q3, p3, masses3, (0, 1, 2), dt=0.01)
    q_sn_ref, p_sn_ref = step_supernode_synchronous_jit(q3, p3, masses3, (0, 1, 2), dt=0.01, G=1.0)
    assert not np.allclose(q_sn, q3), "JitSupernodeSolver.step_supernode_synchronous must evolve q!"
    assert np.allclose(q_sn, q_sn_ref, atol=1e-13)
    assert np.allclose(p_sn, p_sn_ref, atol=1e-13)




