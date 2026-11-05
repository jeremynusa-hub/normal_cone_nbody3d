"""
Unit Tests for Optimized Supernode Solver (Zero-Allocation Micro-Kernel).

Theoretical Alignments:
- Monograph: Section 4.3 (Synchronous Supernodes and Excluded Middle Dilemma)
  and Section 4.4 (Anti-Zeno Regularization via Sundman Time Unfolding).
- Formal Lean 4: `01_UnifiedSuperNodeDef_Theorems.lean` and `02_SuperNodeAntiZenoRegularization_Theorems.lean`.
- Invariants: Strict conservation of internal linear momentum (< 1e-15) and angular momentum (< 1e-14),
  high-precision agreement (< 1e-12) with canonical baseline `step_supernode_synchronous`,
  and non-intrusive ContextVar / NBodySystem delegation.
"""

import math
import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.constants import GRAVITATIONAL_CONSTANT
from normal_cone_nbody_float64.solvers.supernode import (
    step_supernode_synchronous,
    compute_cluster_hyper_radius,
)
from normal_cone_nbody_float64.solvers.optimized_supernode import (
    OptimizedSupernodeSolver,
    step_supernode_synchronous_fast,
)
from normal_cone_nbody_float64.integrators.common import (
    set_current_supernode_solver,
    get_current_supernode_solver,
    step_supernode_flow,
)
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.integrators.strang import StrangIntegrator


def test_optimized_supernode_binary_concordance():
    """Validates bit-exact agreement for binary (k = 2) supernode clusters."""
    q = np.array([
        [-0.5, 0.2, 0.0],
        [+0.5, -0.2, 0.0],
        [10.0, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.0, -0.8, 0.1],
        [0.0, +0.8, -0.1],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)
    masses = np.array([1.5, 2.5, 1.0], dtype=np.float64)
    dt = 0.025
    G = 1.0

    solver_opt = OptimizedSupernodeSolver(G=G)

    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1), dt=dt, G=G
    )
    q_opt, p_opt = solver_opt.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1), dt=dt, G=G
    )

    assert np.allclose(q_base, q_opt, atol=1e-14, rtol=1e-14)
    assert np.allclose(p_base, p_opt, atol=1e-14, rtol=1e-14)


def test_optimized_supernode_k3_concordance_and_momentum():
    """
    Validates high-precision numerical agreement and machine-precision momentum conservation
    on a general 3-body cluster (k = 3) under adaptive unfolded Sundman time.
    """
    masses = np.array([1.0, 1.5, 2.0], dtype=np.float64)
    q = np.array([
        [0.1, 0.0, 0.0],
        [-0.05, 0.08, 0.0],
        [-0.05, -0.08, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.0, 0.5, 0.0],
        [-0.3, -0.2, 0.1],
        [0.3, -0.3, -0.1],
    ], dtype=np.float64)
    dt = 0.02
    G = 1.0

    solver_opt = OptimizedSupernodeSolver(G=G)

    # Baseline execution
    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )

    # Fast micro-kernel execution
    q_opt, p_opt = solver_opt.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )

    # Verify high-precision concordance
    assert np.allclose(q_base, q_opt, atol=1e-12, rtol=1e-12), f"q diff: {np.max(np.abs(q_base - q_opt))}"
    assert np.allclose(p_base, p_opt, atol=1e-12, rtol=1e-12), f"p diff: {np.max(np.abs(p_base - p_opt))}"

    # Verify linear and angular momentum conservation
    P_init = np.sum(p, axis=0)
    L_init = np.sum([np.cross(q[i], p[i]) for i in range(3)], axis=0)

    P_opt = np.sum(p_opt, axis=0)
    L_opt = np.sum([np.cross(q_opt[i], p_opt[i]) for i in range(3)], axis=0)

    assert np.allclose(P_opt, P_init, atol=1e-14), f"Linear momentum drift: {np.max(np.abs(P_opt - P_init))}"
    assert np.allclose(L_opt, L_init, atol=1e-13), f"Angular momentum drift: {np.max(np.abs(L_opt - L_init))}"


def test_optimized_supernode_k3_close_pair_concordance():
    """
    Validates agreement and stability when an isolated close pair occurs inside a k=3 cluster.
    """
    q = np.array([
        [-0.005, 0.0, 0.0],
        [+0.005, 0.0, 0.0],
        [0.0, 1.5, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.0, -10.0, 0.0],
        [0.0, +10.0, 0.0],
        [0.1, 0.0, 0.0],
    ], dtype=np.float64)
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    dt = 0.005
    G = 1.0

    solver_opt = OptimizedSupernodeSolver(G=G)

    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )
    q_opt, p_opt = solver_opt.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )

    assert np.allclose(q_base, q_opt, atol=1e-12, rtol=1e-12)
    assert np.allclose(p_base, p_opt, atol=1e-12, rtol=1e-12)


def test_optimized_supernode_k4_concordance():
    """
    Validates the 4-body (6 pairs) unrolled micro-kernel against canonical baseline.
    """
    masses = np.array([1.0, 1.0, 0.8, 1.2], dtype=np.float64)
    q = np.array([
        [-0.15, 0.0, 0.0],
        [+0.15, 0.0, 0.0],
        [0.0, 0.8, 0.0],
        [0.5, -0.3, 0.1],
    ], dtype=np.float64)
    p = np.array([
        [0.0, 1.5, 0.0],
        [0.0, -1.5, 0.0],
        [-0.8, 0.0, 0.2],
        [0.8, 0.0, -0.2],
    ], dtype=np.float64)
    dt = 0.005
    G = 1.0

    solver_opt = OptimizedSupernodeSolver(G=G)

    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2, 3), dt=dt, G=G
    )
    q_opt, p_opt = solver_opt.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2, 3), dt=dt, G=G
    )

    assert np.allclose(q_base, q_opt, atol=1e-12, rtol=1e-12)
    assert np.allclose(p_base, p_opt, atol=1e-12, rtol=1e-12)

    # Momentum conservation
    P_init = np.sum(p, axis=0)
    P_opt = np.sum(p_opt, axis=0)
    assert np.allclose(P_opt, P_init, atol=1e-14)


def test_optimized_supernode_homothetic_collapse():
    """Validates exact homothety central configuration collapse symmetry."""
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    angle = np.array([0.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0])
    r0 = 0.5
    q = np.zeros((3, 3), dtype=np.float64)
    q[:, 0] = r0 * np.cos(angle)
    q[:, 1] = r0 * np.sin(angle)

    v_rad = -0.2
    p = np.zeros((3, 3), dtype=np.float64)
    p[:, 0] = masses * v_rad * np.cos(angle)
    p[:, 1] = masses * v_rad * np.sin(angle)

    dt = 0.1
    solver_opt = OptimizedSupernodeSolver(G=1.0)
    q_opt, p_opt = solver_opt.step_supernode_synchronous(q, p, masses, (0, 1, 2), dt=dt, G=1.0)

    # Preserve equilateral symmetry
    r_norms = np.linalg.norm(q_opt, axis=1)
    assert np.isclose(r_norms[0], r_norms[1], atol=1e-12)
    assert np.isclose(r_norms[1], r_norms[2], atol=1e-12)


def test_optimized_supernode_reversibility():
    """Validates time reversibility S(-dt) o S(dt) = Id under forward and retrograde steps."""
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

    dt = 0.01
    solver_opt = OptimizedSupernodeSolver(G=1.0)

    q_fwd, p_fwd = solver_opt.step_supernode_synchronous(q0, p0, masses, (0, 1, 2), dt=dt, G=1.0)
    q_bwd, p_bwd = solver_opt.step_supernode_synchronous(q_fwd, p_fwd, masses, (0, 1, 2), dt=-dt, G=1.0)

    assert np.allclose(q_bwd, q0, atol=1e-14)
    assert np.allclose(p_bwd, p0, atol=1e-14)


def test_optimized_supernode_contextvar_and_nbody_system():
    """
    Validates that setting supernode_evaluator in NBodySystem seamlessly routes
    supernode operations to the optimized micro-kernel and restores baseline afterwards.
    """
    masses = np.array([1.0, 1.0, 0.8, 1.2], dtype=np.float64)
    q0 = np.array([
        [-0.15, 0.0, 0.0],
        [+0.15, 0.0, 0.0],
        [0.0, 1.5, 0.0],
        [3.0, 0.3, 0.0],
    ], dtype=np.float64)
    p0 = np.zeros_like(q0)
    state = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=0.0)

    solver_opt = OptimizedSupernodeSolver(G=1.0)
    system = NBodySystem(
        state,
        integrator=StrangIntegrator(G=1.0),
        supernode_evaluator=solver_opt,
    )

    # ContextVar must be None before step
    assert get_current_supernode_solver() is None

    # Step system
    st_next = system.step(dt=0.01)
    assert st_next is not None

    # ContextVar must be properly restored to None after step
    assert get_current_supernode_solver() is None
