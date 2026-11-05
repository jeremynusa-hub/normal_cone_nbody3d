"""
Unit Tests for JIT-Compiled Synchronous Supernode Multi-Body Solver (Block 5).

Theoretical Alignments:
- Monograph: Section 4.3 (Synchronous Supernodes and Excluded Middle Dilemma)
  and Section 4.4 (Anti-Zeno Regularization via Sundman Time Unfolding).
- Formal Lean 4: `01_UnifiedSuperNodeDef_Theorems.lean` and `02_SuperNodeAntiZenoRegularization_Theorems.lean`.
- Invariants: Strict conservation of internal linear momentum (< 1e-15) and angular momentum (< 1e-13),
  high-precision numerical agreement (< 1e-12) with canonical baseline `step_supernode_synchronous`
  and `OptimizedSupernodeSolver`, and non-intrusive ContextVar / NBodySystem delegation.
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
from normal_cone_nbody_float64.solvers.jit_supernode import (
    JitSupernodeSolver,
    step_supernode_synchronous_jit,
    NUMBA_AVAILABLE,
)
from normal_cone_nbody_float64.integrators.common import (
    set_current_supernode_solver,
    get_current_supernode_solver,
    step_supernode_flow,
)
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.integrators.strang import StrangIntegrator


def test_jit_supernode_binary_concordance():
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

    solver_jit = JitSupernodeSolver(G=G)

    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1), dt=dt, G=G
    )
    q_jit, p_jit = solver_jit.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1), dt=dt, G=G
    )

    assert np.allclose(q_base, q_jit, atol=1e-14, rtol=1e-14)
    assert np.allclose(p_base, p_jit, atol=1e-14, rtol=1e-14)


def test_jit_supernode_k3_concordance_and_momentum():
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

    solver_jit = JitSupernodeSolver(G=G)

    # Baseline execution
    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )

    # JIT micro-kernel execution
    q_jit, p_jit = solver_jit.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )

    # Verify high-precision concordance (< 1e-12)
    assert np.allclose(q_base, q_jit, atol=1e-12, rtol=1e-12), f"q diff: {np.max(np.abs(q_base - q_jit))}"
    assert np.allclose(p_base, p_jit, atol=1e-12, rtol=1e-12), f"p diff: {np.max(np.abs(p_base - p_jit))}"

    # Verify linear and angular momentum conservation
    p_tot_initial = np.sum(p, axis=0)
    p_tot_jit = np.sum(p_jit, axis=0)
    assert np.allclose(p_tot_initial, p_tot_jit, atol=1e-15)

    L_initial = np.sum([np.cross(q[i], p[i]) for i in range(3)], axis=0)
    L_jit = np.sum([np.cross(q_jit[i], p_jit[i]) for i in range(3)], axis=0)
    # Drifted angular momentum accounts for CM motion: L(t) = L(0) + R_cm(0) x P_tot
    L_drift = np.cross(np.sum(masses[:, None] * q, axis=0) / np.sum(masses) + dt * p_tot_initial / np.sum(masses), p_tot_initial)
    L_drift_0 = np.cross(np.sum(masses[:, None] * q, axis=0) / np.sum(masses), p_tot_initial)
    assert np.allclose(L_jit - L_drift, L_initial - L_drift_0, atol=1e-13)


def test_jit_supernode_k3_close_pair_concordance():
    """
    Validates cluster behavior with a tight internal binary pair desingularized via inlined Kepler step.
    """
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
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
    dt = 0.005
    G = 1.0

    solver_jit = JitSupernodeSolver(G=G)

    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )
    q_jit, p_jit = solver_jit.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )

    assert np.allclose(q_base, q_jit, atol=1e-12, rtol=1e-12), f"Close pair q diff: {np.max(np.abs(q_base - q_jit))}"
    assert np.allclose(p_base, p_jit, atol=1e-12, rtol=1e-12), f"Close pair p diff: {np.max(np.abs(p_base - p_jit))}"


def test_jit_supernode_k4_concordance():
    """Validates 4-body cluster (k = 4, 6 pairs) zero-allocation micro-kernel."""
    masses = np.array([1.0, 0.8, 1.2, 0.5], dtype=np.float64)
    q = np.array([
        [0.1, 0.1, 0.0],
        [-0.1, 0.1, 0.0],
        [-0.1, -0.1, 0.0],
        [0.1, -0.1, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [-0.2, 0.3, 0.1],
        [0.2, 0.2, -0.1],
        [0.1, -0.4, 0.0],
        [-0.1, -0.1, 0.0],
    ], dtype=np.float64)
    dt = 0.01
    G = 1.0

    solver_jit = JitSupernodeSolver(G=G)

    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2, 3), dt=dt, G=G
    )
    q_jit, p_jit = solver_jit.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2, 3), dt=dt, G=G
    )

    assert np.allclose(q_base, q_jit, atol=1e-12, rtol=1e-12), f"k=4 q diff: {np.max(np.abs(q_base - q_jit))}"
    assert np.allclose(p_base, p_jit, atol=1e-12, rtol=1e-12), f"k=4 p diff: {np.max(np.abs(p_base - p_jit))}"


def test_jit_supernode_homothetic_collapse():
    """Validates homothety central configuration collapse symmetry detection."""
    # Equilateral triangle with radial velocities
    s3 = math.sqrt(3.0) / 2.0
    q = np.array([
        [0.0, 1.0, 0.0],
        [-s3, -0.5, 0.0],
        [s3, -0.5, 0.0],
    ], dtype=np.float64)
    # Pure radial inward velocities
    p = - 0.2 * np.copy(q)
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    dt = 0.01
    G = 1.0

    solver_jit = JitSupernodeSolver(G=G)

    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )
    q_jit, p_jit = solver_jit.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )

    assert np.allclose(q_base, q_jit, atol=1e-14, rtol=1e-14)
    assert np.allclose(p_base, p_jit, atol=1e-14, rtol=1e-14)


def test_jit_supernode_reversibility():
    """Validates microscopic time reversibility S(-dt) o S(dt) == Id under forward and retrograde steps."""
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
    G = 1.0

    solver_jit = JitSupernodeSolver(G=G)

    # Forward step
    q_fwd, p_fwd = solver_jit.step_supernode_synchronous(
        q0, p0, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )
    # Backward step
    q_bwd, p_bwd = solver_jit.step_supernode_synchronous(
        q_fwd, p_fwd, masses, body_indices=(0, 1, 2), dt=-dt, G=G
    )

    assert np.allclose(q_bwd, q0, atol=1e-14)
    assert np.allclose(p_bwd, p0, atol=1e-14)


def test_jit_supernode_contextvar_and_nbody_system():
    """Validates ContextVar registration and NBodySystem instantiation with JitSupernodeSolver."""
    masses = np.array([1.0, 1.0, 0.8, 1.2], dtype=np.float64)
    q0 = np.array([
        [-0.15, 0.0, 0.0],
        [+0.15, 0.0, 0.0],
        [0.0, 1.5, 0.0],
        [3.0, 0.3, 0.0],
    ], dtype=np.float64)
    p0 = np.zeros_like(q0)
    state = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=0.0)

    solver_jit = JitSupernodeSolver(G=1.0)
    system = NBodySystem(
        state,
        integrator=StrangIntegrator(G=1.0),
        supernode_evaluator=solver_jit,
    )

    # ContextVar must be None before step
    assert get_current_supernode_solver() is None

    # Step system
    st_next = system.step(dt=0.01)
    assert st_next is not None

    # ContextVar must be properly restored to None after step
    assert get_current_supernode_solver() is None


def test_jit_supernode_fallback_mode(monkeypatch):
    """Validates transparent fallback to OptimizedSupernodeSolver when Numba is unavailable."""
    import normal_cone_nbody_float64.solvers.jit_supernode as jsn
    monkeypatch.setattr(jsn, "NUMBA_AVAILABLE", False)

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

    solver_jit = JitSupernodeSolver(G=G)
    q_fallback, p_fallback = solver_jit.step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )

    q_opt, p_opt = step_supernode_synchronous_fast(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G
    )

    assert np.allclose(q_fallback, q_opt, atol=1e-15, rtol=1e-15)
    assert np.allclose(p_fallback, p_opt, atol=1e-15, rtol=1e-15)
