"""
Unit Tests for Optimized Kepler Solver Propagation to Synchronous Supernodes.

Validates that when the Kepler solver optimizer flag or ContextVar is active,
the synchronous supernode cluster solver seamlessly routes binary pairs,
close-pair regularizations, and homothetic collapses through `step_kepler_sundman_fast`,
achieving bit-exact physical agreement (< 1e-14) with canonical baseline.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.solvers.kepler_sundman import step_kepler_sundman
from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
    OptimizedKeplerSolver,
    step_kepler_sundman_fast,
)
from normal_cone_nbody_float64.solvers.supernode import (
    step_supernode_synchronous,
    _resolve_kepler_step,
)
from normal_cone_nbody_float64.solvers.optimized_supernode import (
    step_supernode_synchronous_fast,
    OptimizedSupernodeSolver,
)
from normal_cone_nbody_float64.solvers.jit_supernode import (
    step_supernode_synchronous_jit,
    JitSupernodeSolver,
)
from normal_cone_nbody_float64.integrators.common import (
    set_current_kepler_solver,
    get_current_kepler_solver,
    set_current_supernode_solver,
    step_supernode_flow,
)
from normal_cone_nbody_float64.benchmarks.registry import (
    create_burrau_pythagorean_benchmark,
    create_figure8_choreography_benchmark,
)


def test_resolve_kepler_step_routing():
    """Validates dynamic routing between baseline and optimized Kepler solvers."""
    fn_default = _resolve_kepler_step(None)
    assert fn_default is step_kepler_sundman

    solver_opt = OptimizedKeplerSolver()
    fn_explicit = _resolve_kepler_step(solver_opt)
    assert fn_explicit is step_kepler_sundman_fast

    token = set_current_kepler_solver(solver_opt)
    try:
        fn_ctx = _resolve_kepler_step(None)
        assert fn_ctx is step_kepler_sundman_fast
    finally:
        set_current_kepler_solver(None)

    assert _resolve_kepler_step(None) is step_kepler_sundman


def test_supernode_binary_concordance():
    """Validates 0.0 ULP bit-exact agreement for binary (k = 2) supernode clusters."""
    q = np.array([
        [-0.5, 0.2, 0.0],
        [+0.5, -0.2, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.0, -0.8, 0.1],
        [0.0, +0.8, -0.1],
    ], dtype=np.float64)
    masses = np.array([1.5, 2.5], dtype=np.float64)
    dt = 0.025
    G = 1.0

    solver_opt = OptimizedKeplerSolver()

    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1), dt=dt, G=G, kepler_solver=None
    )
    q_opt, p_opt = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1), dt=dt, G=G, kepler_solver=solver_opt
    )
    q_fast, p_fast = step_supernode_synchronous_fast(
        q, p, masses, body_indices=(0, 1), dt=dt, G=G, kepler_solver=solver_opt
    )
    q_jit, p_jit = step_supernode_synchronous_jit(
        q, p, masses, body_indices=(0, 1), dt=dt, G=G, kepler_solver=solver_opt
    )

    assert np.array_equal(q_base, q_opt)
    assert np.array_equal(p_base, p_opt)
    assert np.array_equal(q_base, q_fast)
    assert np.array_equal(p_base, p_fast)
    assert np.array_equal(q_base, q_jit)
    assert np.array_equal(p_base, p_jit)


def test_supernode_close_pair_regularization_concordance():
    """
    Validates 0.0 ULP bit-exact agreement on multi-body clusters (k = 3) across
    baseline, fast, and compiled @njit supernode solvers (including Burrau and Figure-8).
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
    dt = 0.002
    G = 1.0

    solver_opt = OptimizedKeplerSolver()

    q_base, p_base = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G, kepler_solver=None
    )
    q_opt, p_opt = step_supernode_synchronous(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G, kepler_solver=solver_opt
    )
    q_fast, p_fast = step_supernode_synchronous_fast(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G, kepler_solver=solver_opt
    )
    q_jit, p_jit = step_supernode_synchronous_jit(
        q, p, masses, body_indices=(0, 1, 2), dt=dt, G=G, kepler_solver=solver_opt
    )

    assert np.array_equal(q_base, q_opt)
    assert np.array_equal(p_base, p_opt)
    assert np.array_equal(q_base, q_fast)
    assert np.array_equal(p_base, p_fast)
    assert np.array_equal(q_base, q_jit)
    assert np.array_equal(p_base, p_jit)

    # Also test Figure-8 and Burrau 3-4-5 at 0.0 ULP
    for st in (create_figure8_choreography_benchmark(), create_burrau_pythagorean_benchmark(G=1.0)):
        qb, pb = step_supernode_synchronous(st.q, st.p, st.masses, (0, 1, 2), dt=0.05, G=1.0)
        qf, pf = step_supernode_synchronous_fast(st.q, st.p, st.masses, (0, 1, 2), dt=0.05, G=1.0)
        qj, pj = step_supernode_synchronous_jit(st.q, st.p, st.masses, (0, 1, 2), dt=0.05, G=1.0, kepler_solver=solver_opt)
        assert np.array_equal(qb, qf)
        assert np.array_equal(pb, pf)
        assert np.array_equal(qb, qj)
        assert np.array_equal(pb, pj)


def test_step_supernode_flow_propagation():
    """
    Validates 0.0 ULP bit-exact propagation in step_supernode_flow across baseline,
    OptimizedSupernodeSolver, and JitSupernodeSolver without silent @njit bailout.
    """
    q = np.array([
        [-0.005, 0.0, 0.0],
        [+0.005, 0.0, 0.0],
        [0.0, 1.5, 0.0],
        [2.0, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.0, -10.0, 0.0],
        [0.0, +10.0, 0.0],
        [0.1, 0.0, 0.0],
        [-0.1, 0.0, 0.0],
    ], dtype=np.float64)
    masses = np.array([1.0, 1.0, 0.5, 1.0], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0)

    solver_opt = OptimizedKeplerSolver()

    st_base = step_supernode_flow(
        state,
        dt=0.001,
        G=1.0,
        supernode_clusters=[(0, 1, 2)],
        kepler_solver=None,
    )
    st_opt = step_supernode_flow(
        state,
        dt=0.001,
        G=1.0,
        supernode_clusters=[(0, 1, 2)],
        kepler_solver=solver_opt,
    )
    st_fast_cls = OptimizedSupernodeSolver(G=1.0, kepler_solver=solver_opt).step_supernode_flow(
        state,
        dt=0.001,
        G=1.0,
        supernode_clusters=[(0, 1, 2)],
    )
    st_jit_cls = JitSupernodeSolver(G=1.0, kepler_solver=solver_opt).step_supernode_flow(
        state,
        dt=0.001,
        G=1.0,
        supernode_clusters=[(0, 1, 2)],
    )

    assert np.array_equal(st_base.q, st_opt.q)
    assert np.array_equal(st_base.p, st_opt.p)
    assert np.array_equal(st_base.q, st_fast_cls.q)
    assert np.array_equal(st_base.p, st_fast_cls.p)
    assert np.array_equal(st_base.q, st_jit_cls.q)
    assert np.array_equal(st_base.p, st_jit_cls.p)
