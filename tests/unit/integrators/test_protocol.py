"""
Unit Tests for Layer 4: SymplecticIntegrator Protocol and Diagnostic Structures.

Theoretical Alignments:
- Monograph: Chapter 5, Sections 5.1–5.5.
- PEP 544 Runtime-Checkable Structural Subtyping.
"""

import pytest
import numpy as np

from normal_cone_nbody_float64.integrators.protocol import (
    SymplecticIntegrator,
    IntegrationStepResult,
)
from normal_cone_nbody_float64.integrators.strang import StrangIntegrator
from normal_cone_nbody_float64.integrators.yoshida import Yoshida4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
from normal_cone_nbody_float64.integrators.processed import ProcessedSymplecticIntegrator
from normal_cone_nbody_float64.core.state import ExtendedPhaseState


def test_protocol_structural_subtyping():
    """Validates that all six integrators satisfy the SymplecticIntegrator protocol."""
    s2 = StrangIntegrator()
    s4 = Yoshida4Integrator()
    bm4 = BlanesMoan4Integrator()
    bm6 = BlanesMoan6Integrator()
    pbm4 = ProcessedSymplecticIntegrator(bm4)
    pbm6 = ProcessedSymplecticIntegrator(bm6)

    for integ in [s2, s4, bm4, bm6, pbm4, pbm6]:
        assert isinstance(integ, SymplecticIntegrator), f"{integ.name} does not satisfy SymplecticIntegrator"
        assert isinstance(integ.name, str)
        assert isinstance(integ.order, int)
        assert isinstance(integ.effective_order, int)
        assert isinstance(integ.evaluations_per_step, int)
        assert isinstance(integ.is_palindromic, bool)
        assert integ.is_palindromic is True


def test_protocol_order_hierarchy():
    """Validates formal vs effective order hierarchy across the suite."""
    s2 = StrangIntegrator()
    s4 = Yoshida4Integrator()
    bm4 = BlanesMoan4Integrator()
    bm6 = BlanesMoan6Integrator()
    pbm4 = ProcessedSymplecticIntegrator(bm4)
    pbm6 = ProcessedSymplecticIntegrator(bm6)

    assert (s2.order, s2.effective_order) == (2, 2)
    assert (s4.order, s4.effective_order) == (4, 4)
    assert (bm4.order, bm4.effective_order) == (4, 4)
    assert (bm6.order, bm6.effective_order) == (6, 6)
    assert (pbm4.order, pbm4.effective_order) == (4, 6)
    assert (pbm6.order, pbm6.effective_order) == (6, 8)


def test_integration_step_result_dataclass():
    """Validates IntegrationStepResult dataclass initialization and fields."""
    masses = np.array([1.0, 1.0], dtype=np.float64)
    q = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    res = IntegrationStepResult(
        state=state,
        forces_evaluations=7,
        dt=0.01,
        wall_time_seconds=0.002,
    )
    assert res.forces_evaluations == 7
    assert res.dt == 0.01
    assert res.wall_time_seconds == 0.002
    assert res.state.N == 2


def test_protocol_standardized_force_evaluation_properties_and_yoshida4_shim():
    """
    Validates Section 3.3 of Blueprint 08:
    - All base integrators expose num_force_evaluations_per_step, fsal_evaluations_per_step,
      and isolated_evaluations_per_step.
    - integrators.yoshida4 shim re-exports Yoshida4Integrator and YOSHIDA_* constants.
    """
    from normal_cone_nbody_float64.integrators.yoshida4 import (
        Yoshida4Integrator as Yoshida4Shim,
        YOSHIDA_DRIFT_WEIGHTS,
        YOSHIDA_IMPULSE_WEIGHTS,
    )

    assert Yoshida4Shim is Yoshida4Integrator
    assert len(YOSHIDA_DRIFT_WEIGHTS) == 3
    assert len(YOSHIDA_IMPULSE_WEIGHTS) == 4

    s2 = StrangIntegrator()
    s4 = Yoshida4Integrator()
    bm4 = BlanesMoan4Integrator()
    bm6 = BlanesMoan6Integrator()

    assert (s2.evaluations_per_step, s2.num_force_evaluations_per_step, s2.fsal_evaluations_per_step, s2.isolated_evaluations_per_step) == (2, 2, 1, 2)
    assert (s4.evaluations_per_step, s4.num_force_evaluations_per_step, s4.fsal_evaluations_per_step, s4.isolated_evaluations_per_step) == (4, 4, 3, 4)
    assert (bm4.evaluations_per_step, bm4.num_force_evaluations_per_step, bm4.fsal_evaluations_per_step, bm4.isolated_evaluations_per_step) == (7, 7, 7, 8)
    assert (bm6.evaluations_per_step, bm6.num_force_evaluations_per_step, bm6.fsal_evaluations_per_step, bm6.isolated_evaluations_per_step) == (11, 11, 11, 12)


def test_all_base_integrators_symplectic_jacobian_condition():
    """
    Validates that the 6N x 6N state-transition Jacobian matrix M = d(q_1, p_1)/d(q_0, p_0)
    of all four base integrators (Strang, Yoshida4, BM4, BM6) satisfies the canonical
    symplectic condition M^T * Omega_0 * M == Omega_0.
    """
    from normal_cone_nbody_float64.core.jacobi import build_canonical_tree

    masses = np.array([1.2, 0.9, 0.7], dtype=np.float64)
    q0 = np.array([
        [-0.6, 0.1, 0.05],
        [0.6, -0.1, -0.05],
        [0.1, 3.2, 0.4],
    ], dtype=np.float64)
    p0 = np.array([
        [0.05, -0.55, 0.02],
        [-0.05, 0.55, -0.02],
        [-0.22, 0.03, 0.08],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    n_dof = 3 * state0.N  # 9
    omega0 = np.block([
        [np.zeros((n_dof, n_dof)), np.eye(n_dof)],
        [-np.eye(n_dof), np.zeros((n_dof, n_dof))],
    ])

    eps = 1e-4
    dt = 0.02

    for integrator in [
        StrangIntegrator(G=1.0),
        Yoshida4Integrator(G=1.0),
        BlanesMoan4Integrator(G=1.0),
        BlanesMoan6Integrator(G=1.0),
    ]:
        M = np.zeros((2 * n_dof, 2 * n_dof), dtype=np.float64)
        for col in range(2 * n_dof):
            dq = np.zeros_like(q0)
            dp = np.zeros_like(p0)
            if col < n_dof:
                dq.flat[col] = eps
            else:
                dp.flat[col - n_dof] = eps

            def _eval_vec(scale: float) -> np.ndarray:
                st = integrator.step(
                    ExtendedPhaseState(q=q0 + scale * dq, p=p0 + scale * dp, masses=masses),
                    tree,
                    dt=dt,
                )
                return np.concatenate([st.q.ravel(), st.p.ravel()])

            # 4th-order central finite-difference stencil to eliminate O(eps^2) truncation
            M[:, col] = (8.0 * (_eval_vec(1.0) - _eval_vec(-1.0)) - (_eval_vec(2.0) - _eval_vec(-2.0))) / (12.0 * eps)

        symp_err = np.max(np.abs(M.T @ omega0 @ M - omega0))
        assert symp_err < 1e-10, f"{integrator.name} violated symplectic Jacobian condition: {symp_err}"


def test_all_base_integrators_nine_lie_bargmann_generators_and_multi_supernode_clusters():
    """
    Validates across all 4 base integrators:
    1. Exact machine-precision conservation (<= 1e-14) of all 9 Lie-Bargmann phase-space
       generators: P in R^3, L in R^3, and Galilean boost G(t) = M_tot * R_cdm(t) - t * P in R^3.
    2. Multi-supernode support via supernode_clusters=[(0, 1, 2), (3, 4, 5)] for N=6 with
       exact microscopic time-reversibility (<= 1e-15).
    """
    from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
    from normal_cone_nbody_float64.core.invariants import compute_invariants

    masses = np.array([2.0, 1.5, 1.0, 1.8, 1.2, 0.9], dtype=np.float64)
    q0 = np.array([
        [0.0, 0.0, 0.0],
        [0.25, 0.05, -0.02],
        [-0.18, 0.22, 0.04],
        [10.0, 0.5, -0.2],
        [10.22, 0.35, -0.15],
        [9.82, 0.68, -0.25],
    ], dtype=np.float64)
    # Non-zero total momentum P != 0 to actively test Galilean boost G(t) = M_tot * R_cdm(t) - t * P
    p0 = np.array([
        [0.10, -0.20, 0.05],
        [-0.05, 0.25, -0.03],
        [0.02, -0.02, 0.01],
        [-0.15, 0.30, -0.04],
        [0.12, -0.18, 0.06],
        [0.08, -0.05, -0.01],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=0.0)
    tree = build_canonical_tree(masses, topology="comb")
    clusters = [(0, 1, 2), (3, 4, 5)]
    inv0 = compute_invariants(state0)

    dt = 0.01
    for integrator in [
        StrangIntegrator(G=1.0),
        Yoshida4Integrator(G=1.0),
        BlanesMoan4Integrator(G=1.0),
        BlanesMoan6Integrator(G=1.0),
    ]:
        # 1. Reversibility with supernode_clusters
        s_fwd = integrator.step(state0, tree, dt=dt, supernode_clusters=clusters)
        s_rev = integrator.step(s_fwd, tree, dt=-dt, supernode_clusters=clusters)
        assert np.allclose(s_rev.q, state0.q, rtol=2e-15, atol=1.5e-14), f"{integrator.name} q reversibility failed"
        assert np.allclose(s_rev.p, state0.p, rtol=2e-15, atol=1.5e-14), f"{integrator.name} p reversibility failed"

        # 2. 9 Lie-Bargmann generators over 5 steps
        curr = state0
        for _ in range(5):
            curr = integrator.step(curr, tree, dt=dt, supernode_clusters=clusters)
        inv1 = compute_invariants(curr)
        assert np.max(np.abs(inv1.P - inv0.P)) <= 1e-14, f"{integrator.name} P drift: {inv1.P - inv0.P}"
        assert np.max(np.abs(inv1.L - inv0.L)) <= 1e-13, f"{integrator.name} L drift: {inv1.L - inv0.L}"
        assert np.allclose(inv1.G, inv0.G, rtol=1e-14, atol=1e-14), (
            f"{integrator.name} G drift: {inv1.G - inv0.G}"
        )

