"""
Unit Tests for Optimized Kepler-Sundman Diffeomorphic Leaf Solver & Fast Drift.

Validates:
1. Bit-exact equivalence between stumpff_c0123_scalar and stumpff_c[:4].
2. Bit-exact concordance between step_kepler_sundman_fast and step_kepler_sundman (< 1e-14).
3. Exact time-reversibility S(-dt) o S(dt) = Id (< 1e-14).
4. Direct N=3 fast drift concordance against canonical step_drift_and_kepler (< 1e-14).
5. Seamless plug-and-play compatibility in NBodySystem.
"""

import math
import numpy as np
import pytest

from normal_cone_nbody_float64.solvers.stumpff import stumpff_c
from normal_cone_nbody_float64.solvers.kepler_sundman import step_kepler_sundman
from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
    stumpff_c0123_scalar,
    step_kepler_sundman_fast,
    OptimizedKeplerSolver,
)
from normal_cone_nbody_float64.solvers.jit_supernode import _njit_step_kepler_sundman
from normal_cone_nbody_float64.integrators.common import step_drift_and_kepler
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
from normal_cone_nbody_float64.integrators.processed import create_pbm4_integrator
from normal_cone_nbody_float64.benchmarks.registry import (
    create_eccentric_binary_benchmark,
    create_radial_collision_benchmark,
)
from normal_cone_nbody_float64.diagnostics.metrics import compute_forward_backward_reversibility


def test_stumpff_c0123_scalar_exact_match():
    """Validates that stumpff_c0123_scalar matches stumpff_c[:4] bit-exactly (0.0 ULP) across all regimes."""
    test_z_values = [
        -100.0, -10.0, -2.5, -0.5, -0.15, -0.09, -0.05, -0.01, -1e-8,
        0.0,
        1e-8, 0.01, 0.05, 0.09, 0.15, 0.5, 2.5, 10.0, 100.0,
    ]
    for z in test_z_values:
        c0_f, c1_f, c2_f, c3_f, _, _ = stumpff_c(z)
        c0_o, c1_o, c2_o, c3_o = stumpff_c0123_scalar(z)
        assert c0_f == c0_o, f"Mismatch in c0 at z={z}: {c0_f} vs {c0_o}"
        assert c1_f == c1_o, f"Mismatch in c1 at z={z}: {c1_f} vs {c1_o}"
        assert c2_f == c2_o, f"Mismatch in c2 at z={z}: {c2_f} vs {c2_o}"
        assert c3_f == c3_o, f"Mismatch in c3 at z={z}: {c3_f} vs {c3_o}"


def test_step_kepler_sundman_fast_agreement():
    """Validates 0.0 ULP bit-exact agreement between step_kepler_sundman, step_kepler_sundman_fast, and _njit_step_kepler_sundman."""
    mu = 12.0 / 7.0
    GM = 1.0

    test_states = [
        # Circular orbit
        (np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.8, 0.0]), 0.05),
        # Highly eccentric elliptic orbit (e ~ 0.95)
        (np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.15, 0.0]), 0.05),
        # Hyperbolic orbit
        (np.array([0.5, 1.0, 0.0]), np.array([1.5, 1.0, 0.5]), 0.05),
        # Retrograde step (negative dt)
        (np.array([1.0, 0.5, -0.2]), np.array([-0.2, 0.6, 0.4]), -0.04),
        # Small step
        (np.array([2.0, -1.0, 0.5]), np.array([0.1, -0.3, 0.2]), 1e-4),
    ]

    for r0, w0, dt in test_states:
        r_base, w_base = step_kepler_sundman(r0, w0, mu, GM, dt)
        r_fast, w_fast = step_kepler_sundman_fast(r0, w0, mu, GM, dt)
        rx_j, ry_j, rz_j, wx_j, wy_j, wz_j = _njit_step_kepler_sundman(
            float(r0[0]), float(r0[1]), float(r0[2]),
            float(w0[0]), float(w0[1]), float(w0[2]),
            mu, GM, dt
        )
        r_jit = np.array([rx_j, ry_j, rz_j], dtype=np.float64)
        w_jit = np.array([wx_j, wy_j, wz_j], dtype=np.float64)

        assert np.array_equal(r_base, r_fast), f"r mismatch (fast) for dt={dt}: {r_base - r_fast}"
        assert np.array_equal(w_base, w_fast), f"w mismatch (fast) for dt={dt}: {w_base - w_fast}"
        assert np.array_equal(r_base, r_jit), f"r mismatch (jit) for dt={dt}: {r_base - r_jit}"
        assert np.array_equal(w_base, w_jit), f"w mismatch (jit) for dt={dt}: {w_base - w_jit}"


def test_extreme_pericenter_traversal_1e_minus_12_zero_ulp():
    """
    Validates 0.0 ULP bit-exact isomorphism across the extreme pericenter branch
    (1 - e = 1e-12, rho_raw < 0.25 * rho0) and across a full 129-step L1 orbit.
    """
    e_ext = 1.0 - 1e-12
    L_val = math.sqrt(1.0 - e_ext * e_ext)
    vy_in = L_val
    vx_in = -math.sqrt(max(1.0 - vy_in * vy_in, 0.0))
    r0_in_3d = np.array([0.8, 0.6, 0.0], dtype=np.float64)
    w0_in_3d = np.array([0.8 * vx_in - 0.6 * vy_in, 0.6 * vx_in + 0.8 * vy_in, 1e-7], dtype=np.float64)
    dt_peri = math.pi / (2.0 * math.sqrt(2.0))

    r_base, w_base = step_kepler_sundman(r0_in_3d, w0_in_3d, 1.0, 1.0, dt_peri)
    r_fast, w_fast = step_kepler_sundman_fast(r0_in_3d, w0_in_3d, 1.0, 1.0, dt_peri)
    rx_j, ry_j, rz_j, wx_j, wy_j, wz_j = _njit_step_kepler_sundman(
        float(r0_in_3d[0]), float(r0_in_3d[1]), float(r0_in_3d[2]),
        float(w0_in_3d[0]), float(w0_in_3d[1]), float(w0_in_3d[2]),
        1.0, 1.0, dt_peri
    )
    r_jit = np.array([rx_j, ry_j, rz_j], dtype=np.float64)
    w_jit = np.array([wx_j, wy_j, wz_j], dtype=np.float64)

    assert np.array_equal(r_base, r_fast)
    assert np.array_equal(w_base, w_fast)
    assert np.array_equal(r_base, r_jit)
    assert np.array_equal(w_base, w_jit)

    # Full 129-step orbit across pericenter at 1 - e = 1e-12
    st_l1 = create_eccentric_binary_benchmark(e=e_ext, a=1.0, m0=1.0, m1=1.0, G=1.0)
    dt_l1 = (2.0 * math.pi / math.sqrt(2.0)) / 129.0
    sys_base = NBodySystem(state=st_l1.with_updates(), integrator=create_pbm4_integrator(G=1.0), G=1.0, use_optimized_kernels=False)
    sys_opt = NBodySystem(state=st_l1.with_updates(), integrator=create_pbm4_integrator(G=1.0), G=1.0, use_optimized_kernels=True)
    for _ in range(129):
        sys_base.step(dt_l1)
        sys_opt.step(dt_l1)
    assert np.array_equal(sys_base.state.q, sys_opt.state.q)
    assert np.array_equal(sys_base.state.p, sys_opt.state.p)


def test_optimized_kepler_time_reversibility():
    """Validates that S(-dt) o S(dt) = Id to float64 machine precision (< 1e-15)."""
    mu = 1.2
    GM = 2.5
    r0 = np.array([1.5, -0.8, 0.4])
    w0 = np.array([0.3, 1.1, -0.6])
    dt = 0.08

    r_fwd, w_fwd = step_kepler_sundman_fast(r0, w0, mu, GM, dt)
    r_rev, w_rev = step_kepler_sundman_fast(r_fwd, w_fwd, mu, GM, -dt)

    assert np.max(np.abs(r_rev - r0)) < 1e-15
    assert np.max(np.abs(w_rev - w0)) < 1e-15


def test_radial_collision_concordance():
    """Validates 0.0 ULP concordance on head-on radial encounters (e = 1, L2.1a)."""
    mu = 0.5
    GM = 2.0
    r0 = np.array([1.0, 0.0, 0.0])
    w0 = np.array([0.0, 0.0, 0.0])
    dt = 1.5

    r_base, w_base = step_kepler_sundman(r0, w0, mu, GM, dt)
    r_fast, w_fast = step_kepler_sundman_fast(r0, w0, mu, GM, dt)
    rx_j, ry_j, rz_j, wx_j, wy_j, wz_j = _njit_step_kepler_sundman(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, mu, GM, dt)

    assert np.array_equal(r_base, r_fast)
    assert np.array_equal(w_base, w_fast)
    assert np.array_equal(r_base, np.array([rx_j, ry_j, rz_j]))
    assert np.array_equal(w_base, np.array([wx_j, wy_j, wz_j]))

    # L2.1a 80-step forward/backward round-trip through d_0 = 0
    rad_st = create_radial_collision_benchmark(r0=2.0, m0=1.0, m1=1.0, G=1.0)
    fb = compute_forward_backward_reversibility(rad_st, create_pbm4_integrator(G=1.0), dt=0.025, n_steps=80, G=1.0)
    assert fb["delta_E_rel"] < 2e-13


def test_angular_momentum_and_energy_conservation():
    """Validates that orbital energy and angular momentum are conserved to < 1e-15."""
    mu = 1.0
    GM = 1.5
    r0 = np.array([2.0, 0.5, -0.3])
    w0 = np.array([0.2, 0.9, 0.4])

    v0 = w0 / mu
    E0 = 0.5 * np.sum(v0**2) - GM / np.linalg.norm(r0)
    L0 = np.cross(r0, w0)

    dt = 0.25
    r_next, w_next = step_kepler_sundman_fast(r0, w0, mu, GM, dt)

    v_next = w_next / mu
    E_next = 0.5 * np.sum(v_next**2) - GM / np.linalg.norm(r_next)
    L_next = np.cross(r_next, w_next)

    assert abs(E_next - E0) / abs(E0) < 1e-15
    assert np.max(np.abs(L_next - L0)) < 1e-15


def test_optimized_kepler_solver_n3_drift_agreement():
    """Validates 0.0 ULP bit-exact agreement of N=3 drift including non-zero total momentum and large R_out."""
    masses = np.array([1.0, 1e-3, 1e-3], dtype=np.float64)
    tree = build_canonical_tree(masses, topology="comb")

    q0 = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [100.0, 50.0, -25.0],
    ], dtype=np.float64)
    v0 = np.array([
        [0.01, -0.02, 0.005],
        [0.01, 0.98, 0.005],
        [0.01 + 1e-4, -0.02 + 2e-4, 0.005 - 1e-4],
    ], dtype=np.float64)
    p0 = v0 * masses[:, None]
    state = ExtendedPhaseState(q=q0, p=p0, masses=masses)

    solver = OptimizedKeplerSolver(G=1.0)
    dt = 1e-3

    state_base = step_drift_and_kepler(state, tree, dt, G=1.0, kepler_solver=step_kepler_sundman)
    state_opt = solver.step_drift_and_kepler(state, tree, dt, G=1.0)

    assert np.array_equal(state_base.q, state_opt.q), f"Drift q mismatch: {state_base.q - state_opt.q}"
    assert np.array_equal(state_base.p, state_opt.p), f"Drift p mismatch: {state_base.p - state_opt.p}"


def test_nbody_system_with_kepler_evaluator():
    """Validates 0.0 ULP bit-exact NBodySystem simulation with OptimizedKeplerSolver."""
    masses = np.array([3.0, 4.0, 5.0])
    q0 = np.array([
        [1.0, 3.0, 0.0],
        [-2.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
    ])
    p0 = np.zeros((3, 3), dtype=np.float64)
    state = ExtendedPhaseState(q=q0, p=p0, masses=masses)

    integrator_base = BlanesMoan6Integrator(G=1.0)
    integrator_opt = BlanesMoan6Integrator(G=1.0)

    system_base = NBodySystem(
        state=state,
        integrator=integrator_base,
        G=1.0,
    )
    system_opt = NBodySystem(
        state=state,
        integrator=integrator_opt,
        G=1.0,
        kepler_evaluator=OptimizedKeplerSolver(G=1.0),
    )

    for _ in range(10):
        system_base.step(dt=0.01)
        system_opt.step(dt=0.01)

    assert np.array_equal(system_base.state.q, system_opt.state.q)
    assert np.array_equal(system_base.state.p, system_opt.state.p)
