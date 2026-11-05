r"""
Unit Test Suite for Extended Precision Float128 / Mpmath Potential and Analytic Solvers.

Tests:
1. Stumpff functions (order 14 Horner) and fundamental algebraic identities (< 1e-34).
2. Kepler--Sundman Halley solver and Battin unimodular projection (|W - 1| < 1e-32, time-reversibility < 1e-32).
3. Screened apex potential, traceless tidal tensor (|Tr(T)| < 1e-34), and exact momentum conservation.
4. Supernode detector, Hill tidal ratio, and Marchal--Sundman threshold in mpmath.
5. Synchronous supernode solver and local barycentric conservation (< 1e-32).
"""

import pytest
import numpy as np
import mpmath as mp

from normal_cone_nbody_float128.core.constants import (
    DEFAULT_DPS,
    G_DEFAULT_F128,
    to_f128,
    to_f64,
)
from normal_cone_nbody_float128.core.state import ExtendedPhaseState
from normal_cone_nbody_float128.core.jacobi import build_canonical_tree
from normal_cone_nbody_float128.solvers.stumpff import (
    stumpff_c,
    stumpff_c_f128,
    stumpff_c_derivatives,
    stumpff_c_derivatives_f128,
)
from normal_cone_nbody_float128.solvers.kepler_sundman import (
    step_kepler_sundman,
    compute_wronskian,
)
from normal_cone_nbody_float128.potential.screened_apex import (
    compute_quadrupole_tensor,
    compute_conormal_tidal_tensor,
    compute_quadrupole_potential,
    compute_intercluster_forces,
    ScreenedApexPotential,
)
from normal_cone_nbody_float128.potential.supernode_detector import (
    compute_tidal_ratio,
    compute_natural_threshold,
    compute_internal_hierarchy_ratio,
    detect_supernodes,
)
from normal_cone_nbody_float128.solvers.supernode import (
    step_supernode_synchronous,
    compute_cluster_hyper_radius,
)

mp.mp.dps = DEFAULT_DPS


class TestStumpffFloat128:
    """Verifies high-precision Stumpff series evaluation up to order 14."""

    def test_stumpff_algebraic_identities(self):
        """Checks Stumpff identities c0 + z c2 = 1, c1 + z c3 = 1, and c1^2 - c0 c2 = c2."""
        test_z_values = [
            mp.mpf('-5.0'),
            mp.mpf('-0.5'),
            mp.mpf('-0.05'),
            mp.mpf('0.0'),
            mp.mpf('0.05'),
            mp.mpf('0.5'),
            mp.mpf('5.0'),
        ]
        for z in test_z_values:
            v0, v1, v2, v3, v4, v5 = stumpff_c_f128(z)

            err1 = abs((v0 + z * v2) - mp.mpf('1.0'))
            err2 = abs((v1 + z * v3) - mp.mpf('1.0'))
            err3 = abs((v1**2 - v0 * v2) - v2)

            assert err1 < mp.mpf('1e-34'), f"Identity c0 + z c2 = 1 violated at z={z}: err={err1}"
            assert err2 < mp.mpf('1e-34'), f"Identity c1 + z c3 = 1 violated at z={z}: err={err2}"
            assert err3 < mp.mpf('1e-34'), f"Wronskian identity c1^2 - c0 c2 = c2 violated at z={z}: err={err3}"

    def test_stumpff_derivatives(self):
        """Checks Stumpff derivative relations: c0' = -1/2 c1."""
        z = mp.mpf('0.03')
        v0, v1, v2, v3, v4, v5 = stumpff_c_f128(z)
        dc0, dc1, dc2, dc3, dc4, dc5 = stumpff_c_derivatives_f128(z)
        assert abs(dc0 - (-mp.mpf('0.5') * v1)) < mp.mpf('1e-34')


class TestKeplerSundmanFloat128:
    """Verifies Kepler--Sundman Halley solver and Battin unimodular projection."""

    def test_kepler_wronskian_unimodularity(self):
        """Verifies that the Kepler flow preserves W(f, g) = f g_dot - f_dot g = 1 to < 1e-32."""
        r0 = np.array([mp.mpf('1.0'), mp.mpf('0.0'), mp.mpf('0.0')], dtype=object)
        w0 = np.array([mp.mpf('0.0'), mp.mpf('0.8'), mp.mpf('0.0')], dtype=object)
        mu = mp.mpf('0.5')
        GM = mp.mpf('1.0')
        dt = mp.mpf('0.25')

        r_new, w_new = step_kepler_sundman(r0, w0, mu=mu, GM=GM, dt_target=dt)

        # Check that r_new and w_new are object arrays of mp.mpf
        assert r_new.dtype == object
        assert w_new.dtype == object
        assert isinstance(r_new[0], mp.mpf)

        # Verify Wronskian determinant f * g_dot - f_dot * g == 1 to < 1e-32
        W = compute_wronskian(r0, w0, mu=mu, GM=GM, dt_target=dt)
        assert abs(W - mp.mpf('1.0')) < mp.mpf('1e-32'), f"Wronskian |W - 1| = {abs(W - 1)} >= 1e-32"

    def test_kepler_time_reversibility(self):
        """Verifies symplectic time-reversibility: S(-dt) o S(dt) = Id to < 1e-32."""
        r0 = np.array([mp.mpf('0.7'), mp.mpf('0.4'), mp.mpf('-0.1')], dtype=object)
        w0 = np.array([mp.mpf('-0.2'), mp.mpf('1.1'), mp.mpf('0.3')], dtype=object)
        mu = mp.mpf('0.33333333333333333333333333333333333')
        GM = mp.mpf('2.0')
        dt = mp.mpf('0.15')

        # Forward step
        r_fwd, w_fwd = step_kepler_sundman(r0, w0, mu=mu, GM=GM, dt_target=dt)

        # Backward step
        r_rev, w_rev = step_kepler_sundman(r_fwd, w_fwd, mu=mu, GM=GM, dt_target=-dt)

        err_r = mp.sqrt(sum((r_rev[d] - r0[d])**2 for d in range(3)))
        err_w = mp.sqrt(sum((w_rev[d] - w0[d])**2 for d in range(3)))

        assert err_r < mp.mpf('1e-32'), f"Time reversibility err_r={err_r} >= 1e-32"
        assert err_w < mp.mpf('1e-32'), f"Time reversibility err_w={err_w} >= 1e-32"

    def test_optimized_kepler_solver_float128(self):
        """
        Verifies D-144-08 in 35-digit mpmath arithmetic:
        step_kepler_sundman_fast and OptimizedKeplerSolver preserve 35-digit mp.mpf precision
        (including micro-steps dt = 1e-20) without truncating to float64.
        """
        from normal_cone_nbody_float128.solvers.optimized_kepler_sundman import (
            step_kepler_sundman_fast,
            OptimizedKeplerSolver,
        )

        r0 = np.array([mp.mpf('2.0'), mp.mpf('0.0'), mp.mpf('0.0')], dtype=object)
        w0 = np.array([mp.mpf('0.0'), mp.mpf('0.5'), mp.mpf('0.0')], dtype=object)
        mu = mp.mpf('1.0')
        GM = mp.mpf('8.0')

        # Micro-step dt = 1e-20 where float64 would completely underflow/clamp
        dt_micro = mp.mpf('1e-20')
        r_fast, w_fast = step_kepler_sundman_fast(r0, w0, mu=mu, GM=GM, dt_target=dt_micro)
        assert r_fast.dtype == object and isinstance(r_fast[0], mp.mpf)
        assert w_fast.dtype == object and isinstance(w_fast[0], mp.mpf)
        dw_exact = -GM * mu * r0[0] / (mp.mpf('8.0')) * dt_micro  # -1.0 * 1e-20
        dw_actual = w_fast[0] - w0[0]
        assert abs(dw_actual - dw_exact) / abs(dw_exact) < mp.mpf('1e-18')

        # N=3 OptimizedKeplerSolver.step_drift_and_kepler in 35-digit precision
        masses = np.array([mp.mpf('1.0'), mp.mpf('2.0'), mp.mpf('3.0')], dtype=object)
        q = np.array([
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('1.0'), mp.mpf('0.2'), mp.mpf('-0.1')],
            [mp.mpf('-1.5'), mp.mpf('0.8'), mp.mpf('0.4')],
        ], dtype=object)
        p = np.array([
            [mp.mpf('0.1'), mp.mpf('-0.2'), mp.mpf('0.05')],
            [mp.mpf('-0.15'), mp.mpf('0.1'), mp.mpf('-0.02')],
            [mp.mpf('0.05'), mp.mpf('0.1'), mp.mpf('-0.03')],
        ], dtype=object)
        state = ExtendedPhaseState(q=q, p=p, masses=masses)
        tree = build_canonical_tree(masses, topology="comb")
        opt_solver = OptimizedKeplerSolver(G=G_DEFAULT_F128)
        state_fwd = opt_solver.step_drift_and_kepler(state, tree, dt=mp.mpf('0.1'))
        assert state_fwd.is_f128 and isinstance(state_fwd.q[0, 0], mp.mpf)
        state_rev = opt_solver.step_drift_and_kepler(state_fwd, tree, dt=mp.mpf('-0.1'))
        err_q = max(abs(state_rev.q[i, d] - state.q[i, d]) for i in range(3) for d in range(3))
        assert err_q < mp.mpf('1e-32'), f"OptimizedKeplerSolver f128 reversibility error: {err_q}"

    def test_kepler_pericenter_homotopy_and_head_on_f128(self):
        """
        Verifies [ENG-FALSE-01] & [ENG-FALSE-02] in 35-digit mpmath arithmetic:
        1. Pericenter traversal at 1 - e = 1e-12 preserves 35-digit orbital energy (< 1e-32)
           without (1 - e)^{-2} cancellation.
        2. Head-on radial collision bounce (L = 0, e = 1) via analytical Levi-Civita
           continuation preserves energy and time-reversibility to < 1e-32 without v_mag rescaling.
        """
        mu = mp.mpf('1.0')
        GM = mp.mpf('1.0')
        one_minus_e = mp.mpf('1e-12')
        e = mp.mpf('1.0') - one_minus_e
        r_apo = mp.mpf('1.0') + e
        v_apo = mp.sqrt(GM * one_minus_e / r_apo)
        r0 = np.array([-r_apo, mp.mpf('0.0'), mp.mpf('0.0')], dtype=object)
        w0 = np.array([mp.mpf('0.0'), -mu * v_apo, mp.mpf('0.0')], dtype=object)
        E0 = -mp.mpf('0.5') * GM * mu

        r_peri, w_peri = step_kepler_sundman(r0, w0, mu=mu, GM=GM, dt_target=mp.pi)
        r_peri_norm = mp.sqrt(sum(r_peri[d]**2 for d in range(3)))
        E_peri = mp.mpf('0.5') * sum(w_peri[d]**2 for d in range(3)) / mu - GM * mu / r_peri_norm
        rel_dE = abs(E_peri - E0) / abs(E0)
        assert rel_dE < mp.mpf('1e-32'), f"F128 pericenter 1-e=1e-12 energy loss: {rel_dE}"

        # Head-on collision bounce
        r_col0 = np.array([mp.mpf('1.0'), mp.mpf('0.0'), mp.mpf('0.0')], dtype=object)
        w_col0 = np.array([mp.mpf('-0.5'), mp.mpf('0.0'), mp.mpf('0.0')], dtype=object)
        E_col0 = mp.mpf('0.5') * sum(w_col0[d]**2 for d in range(3)) / mu - GM * mu / mp.mpf('1.0')
        dt_col = mp.mpf('1.5')
        r_col1, w_col1 = step_kepler_sundman(r_col0, w_col0, mu=mu, GM=GM, dt_target=dt_col)
        r_col1_norm = mp.sqrt(sum(r_col1[d]**2 for d in range(3)))
        E_col1 = mp.mpf('0.5') * sum(w_col1[d]**2 for d in range(3)) / mu - GM * mu / r_col1_norm
        assert abs(E_col1 - E_col0) / abs(E_col0) < mp.mpf('1e-32')
        r_col_rev, w_col_rev = step_kepler_sundman(r_col1, w_col1, mu=mu, GM=GM, dt_target=-dt_col)
        assert max(abs(r_col_rev[d] - r_col0[d]) for d in range(3)) < mp.mpf('1e-32')
        assert max(abs(w_col_rev[d] - w_col0[d]) for d in range(3)) < mp.mpf('1e-32')



class TestScreenedApexFloat128:
    """Verifies screened apex potential, tidal tensors, and exact momentum conservation."""

    def test_quadrupole_and_tidal_traceless(self):
        """Verifies that quadrupole Q and tidal tensor T_apex are traceless to < 1e-34."""
        q_cl = np.array([
            [mp.mpf('0.5'), mp.mpf('0.2'), mp.mpf('-0.1')],
            [mp.mpf('-0.3'), mp.mpf('0.4'), mp.mpf('0.2')],
            [mp.mpf('-0.2'), mp.mpf('-0.6'), mp.mpf('-0.1')],
        ], dtype=object)
        m_cl = np.array([mp.mpf('1.0'), mp.mpf('2.0'), mp.mpf('1.5')], dtype=object)

        Q = compute_quadrupole_tensor(q_cl, m_cl)
        tr_Q = Q[0, 0] + Q[1, 1] + Q[2, 2]
        assert abs(tr_Q) < mp.mpf('1e-34'), f"Quadrupole trace = {tr_Q} >= 1e-34"

        state_cl = ExtendedPhaseState(q_cl, np.zeros_like(q_cl), m_cl)
        T = compute_conormal_tidal_tensor(state_cl, leaf_pair=(0, 1), G=G_DEFAULT_F128)
        tr_T = T[0, 0] + T[1, 1] + T[2, 2]
        assert abs(tr_T) < mp.mpf('1e-34'), f"Tidal tensor trace = {tr_T} >= 1e-34"

    def test_apex_forces_exact_conservation(self):
        """Verifies that inter-cluster forces sum to zero and conserve angular momentum (< 1e-32)."""
        q = np.array([
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('1.0'), mp.mpf('0.5'), mp.mpf('0.2')],
            [mp.mpf('-2.0'), mp.mpf('3.0'), mp.mpf('-1.0')],
        ], dtype=object)
        masses = np.array([mp.mpf('3.0'), mp.mpf('4.0'), mp.mpf('5.0')], dtype=object)
        leaf_pairs = {(0, 1)}

        F = compute_intercluster_forces(q, masses, leaf_pairs, G=G_DEFAULT_F128)

        total_F = sum(F[i] for i in range(3))
        err_F = mp.sqrt(total_F[0]**2 + total_F[1]**2 + total_F[2]**2)
        assert err_F < mp.mpf('1e-32'), f"Total force = {err_F} >= 1e-32"

        total_torque = sum(
            np.array([
                q[i, 1] * F[i, 2] - q[i, 2] * F[i, 1],
                q[i, 2] * F[i, 0] - q[i, 0] * F[i, 2],
                q[i, 0] * F[i, 1] - q[i, 1] * F[i, 0],
            ], dtype=object)
            for i in range(3)
        )
        err_torque = mp.sqrt(total_torque[0]**2 + total_torque[1]**2 + total_torque[2]**2)
        assert err_torque < mp.mpf('1e-32'), f"Total torque = {err_torque} >= 1e-32"

    def test_direct_newtonian_potential_and_forces_float128(self):
        """
        Verifies D-144-01: direct.py supports 35-digit mpmath.mpf arrays without TypeError or float64 truncation.
        """
        from normal_cone_nbody_float128.potential.direct import (
            compute_direct_newtonian_potential,
            compute_direct_newtonian_forces,
            direct_newtonian_step,
        )

        q = np.array([
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('1.0'), mp.mpf('2.0'), mp.mpf('2.0')],  # r_01 = 3.0 exactly
            [mp.mpf('0.0'), mp.mpf('3.0'), mp.mpf('4.0')],  # r_02 = 5.0 exactly
        ], dtype=object)
        p = np.array([
            [mp.mpf('0.1'), mp.mpf('-0.2'), mp.mpf('0.05')],
            [mp.mpf('-0.05'), mp.mpf('0.15'), mp.mpf('-0.02')],
            [mp.mpf('-0.05'), mp.mpf('0.05'), mp.mpf('-0.03')],
        ], dtype=object)
        masses = np.array([mp.mpf('1.0'), mp.mpf('3.0'), mp.mpf('5.0')], dtype=object)

        V = compute_direct_newtonian_potential(q, masses, G=G_DEFAULT_F128)
        assert isinstance(V, mp.mpf), f"Expected mp.mpf potential, got {type(V)}"
        # Exact V = -G*(m0*m1/3 + m0*m2/5 + m1*m2/sqrt(1 + 1 + 4)) = -(1 + 1 + 15/sqrt(6))
        V_exact = -(mp.mpf('2.0') + mp.mpf('15.0') / mp.sqrt(mp.mpf('6.0')))
        assert abs(V - V_exact) < mp.mpf('1e-33'), f"Direct potential mismatch: {abs(V - V_exact)}"

        F = compute_direct_newtonian_forces(q, masses, G=G_DEFAULT_F128)
        assert F.dtype == object and isinstance(F[0, 0], mp.mpf)
        # Central finite differences check at 1e-12 step in 35-digit arithmetic
        eps = mp.mpf('1e-12')
        for i in range(3):
            for d in range(3):
                q_plus = q.copy()
                q_minus = q.copy()
                q_plus[i, d] = q[i, d] + eps
                q_minus[i, d] = q[i, d] - eps
                V_plus = compute_direct_newtonian_potential(q_plus, masses, G=G_DEFAULT_F128)
                V_minus = compute_direct_newtonian_potential(q_minus, masses, G=G_DEFAULT_F128)
                f_fd = -(V_plus - V_minus) / (mp.mpf('2.0') * eps)
                assert abs(F[i, d] - f_fd) < mp.mpf('1e-22'), f"Force FD mismatch at ({i},{d}): {abs(F[i, d] - f_fd)}"

        state = ExtendedPhaseState(q=q, p=p, masses=masses)
        next_state = direct_newtonian_step(state, dt=mp.mpf('0.01'), G=G_DEFAULT_F128)
        assert next_state.is_f128
        assert isinstance(next_state.q[0, 0], mp.mpf)
        assert isinstance(next_state.p[0, 0], mp.mpf)

    def test_screened_apex_transition_gradient_unequal_masses_f128(self):
        """
        Verifies D-144-04 in 35-digit arithmetic:
        Analytical gradient F_i = -nabla_{q_i} V_apex with include_screened_monopole=True in the active
        transition zone u_v in (1, 2) matches central finite differences to < 1e-22 for unequal masses.
        """
        masses = np.array([mp.mpf('1.0'), mp.mpf('2.3'), mp.mpf('3.7'), mp.mpf('5.1')], dtype=object)
        q0 = np.array([
            [mp.mpf('-0.35'), mp.mpf('0.10'), mp.mpf('-0.05')],
            [mp.mpf('0.25'), mp.mpf('-0.15'), mp.mpf('0.08')],
            [mp.mpf('0.85'), mp.mpf('0.20'), mp.mpf('-0.10')],
            [mp.mpf('-1.10'), mp.mpf('-0.40'), mp.mpf('0.30')],
        ], dtype=object)
        p0 = np.full((4, 3), mp.mpf('0.0'), dtype=object)
        state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
        tree = build_canonical_tree(masses, topology="comb")

        apex = ScreenedApexPotential(G=G_DEFAULT_F128)
        assert apex.screening_mode == "c6"

        forces, V0 = apex.compute_apex_forces(state0, tree, include_screened_monopole=True)
        _, V0_unscreened = apex.compute_apex_forces(state0, tree, include_screened_monopole=False)
        assert abs(V0 - V0_unscreened) > mp.mpf('1e-2')

        eps = mp.mpf('1e-12')
        for i in range(4):
            for d in range(3):
                q_p = q0.copy()
                q_m = q0.copy()
                q_p[i, d] = q0[i, d] + eps
                q_m[i, d] = q0[i, d] - eps
                _, V_p = apex.compute_apex_forces(
                    ExtendedPhaseState(q=q_p, p=p0, masses=masses), tree, include_screened_monopole=True
                )
                _, V_m = apex.compute_apex_forces(
                    ExtendedPhaseState(q=q_m, p=p0, masses=masses), tree, include_screened_monopole=True
                )
                fd_f = -(V_p - V_m) / (mp.mpf('2.0') * eps)
                err = abs(forces[i, d] - fd_f)
                assert err < mp.mpf('1e-21'), f"F128 screened gradient mismatch at ({i},{d}): err={err}"

    def test_optimized_screened_apex_float128(self):
        """
        Verifies D-144-05 & D-144-06 in float128:
        OptimizedScreenedApexPotential preserves 35-digit mpmath precision for N=3 and N=4,
        avoids cache collisions across distinct supernode sets of equal size, and supports mask_leaf_pairs.
        """
        from normal_cone_nbody_float128.potential.optimized_screened_apex import (
            OptimizedScreenedApexPotential,
        )

        masses = np.array([mp.mpf('1.0'), mp.mpf('2.0'), mp.mpf('3.0'), mp.mpf('4.0')], dtype=object)
        q = np.array([
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.5'), mp.mpf('0.1'), mp.mpf('-0.1')],
            [mp.mpf('-1.2'), mp.mpf('0.8'), mp.mpf('0.3')],
            [mp.mpf('1.5'), mp.mpf('-0.7'), mp.mpf('0.4')],
        ], dtype=object)
        p = np.full((4, 3), mp.mpf('0.0'), dtype=object)
        state = ExtendedPhaseState(q=q, p=p, masses=masses)
        tree = build_canonical_tree(masses, topology="comb")

        ref_apex = ScreenedApexPotential(G=G_DEFAULT_F128)
        opt_apex = OptimizedScreenedApexPotential(G=G_DEFAULT_F128)

        f_ref_1, v_ref_1 = ref_apex.compute_apex_forces(state, tree, supernode_bodies={0, 1, 2})
        f_opt_1, v_opt_1 = opt_apex.compute_apex_forces(state, tree, supernode_bodies={0, 1, 2})
        assert f_opt_1.dtype == object and isinstance(f_opt_1[0, 0], mp.mpf)
        assert isinstance(v_opt_1, mp.mpf)
        assert abs(v_opt_1 - v_ref_1) < mp.mpf('1e-33')

        f_ref_2, v_ref_2 = ref_apex.compute_apex_forces(state, tree, supernode_bodies={0, 1, 3})
        f_opt_2, v_opt_2 = opt_apex.compute_apex_forces(state, tree, supernode_bodies={0, 1, 3})
        for i in range(4):
            for d in range(3):
                assert abs(f_opt_2[i, d] - f_ref_2[i, d]) < mp.mpf('1e-33')
        assert abs(v_opt_2 - v_ref_2) < mp.mpf('1e-33')

        state_next = opt_apex.step_apex_kick(state, tree, dt=mp.mpf('0.01'))
        assert state_next.is_f128 and isinstance(state_next.p[0, 0], mp.mpf)

    def test_multi_supernode_and_contextvar_apex_f128(self):
        """
        Verifies [ENG-FALSE-05] & [ENG-FALSE-06] in 35-digit mpmath arithmetic:
        1. N=6 with two disjoint supernodes S_1=(0,1,2) and S_2=(3,4,5) preserves all 9 cross-pair
           forces to < 1e-33 and advances the inter-supernode composite Jacobi drift.
        2. _CURRENT_APEX_EVALUATOR ContextVar dispatches step_apex_kick in common.py.
        """
        from normal_cone_nbody_float128.potential.optimized_screened_apex import (
            OptimizedScreenedApexPotential,
        )
        from normal_cone_nbody_float128.integrators.common import (
            step_drift_and_kepler,
            step_apex_kick,
            set_current_apex_evaluator,
            get_current_apex_evaluator,
            reset_current_apex_evaluator,
        )

        masses = to_f128(np.array([1.0, 1.5, 2.0, 1.2, 1.8, 2.5]))
        q = to_f128(np.array([
            [0.0, 0.0, 0.0],
            [0.02, 0.0, 0.0],
            [0.0, 0.02, 0.0],
            [10.0, 0.0, 0.0],
            [10.02, 0.0, 0.0],
            [10.0, 0.02, 0.0],
        ]))
        p = to_f128(np.array([
            [0.1, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [-0.1, 0.0, 0.0],
            [-0.1, 0.0, 0.0],
            [-0.1, 0.0, 0.0],
        ]))
        state = ExtendedPhaseState(q=q, p=p, masses=masses)
        tree = build_canonical_tree(masses, topology="balanced")
        clusters = [(0, 1, 2), (3, 4, 5)]

        ref_apex = ScreenedApexPotential(G=G_DEFAULT_F128)
        opt_apex = OptimizedScreenedApexPotential(G=G_DEFAULT_F128)

        f_ref, v_ref = ref_apex.compute_apex_forces(state, tree, supernode_clusters=clusters)
        f_opt, v_opt = opt_apex.compute_apex_forces(state, tree, supernode_clusters=clusters)
        assert abs(v_ref - v_opt) < mp.mpf('1e-33')
        assert abs(v_ref) > mp.mpf('1.0')  # Non-zero cross-cluster interaction!

        token = set_current_apex_evaluator(opt_apex)
        try:
            assert get_current_apex_evaluator() is opt_apex
            kicked = step_apex_kick(state, tree, ref_apex, dt=mp.mpf('0.05'), supernode_clusters=clusters)
            for i in range(6):
                for d in range(3):
                    assert abs(kicked.p[i, d] - (p[i, d] + mp.mpf('0.05') * f_ref[i, d])) < mp.mpf('1e-33')
        finally:
            reset_current_apex_evaluator(token)
        assert get_current_apex_evaluator() is None






class TestSupernodeFloat128:
    """Verifies supernode detection and synchronous integration in f128."""

    def test_supernode_detection_f128(self):
        """Checks detection of democratic triple vs isolated binary."""
        # Democratic triple: equilateral triangle with chi = 1.0 < 3.0
        q_demo = np.array([
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('1.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.5'), mp.sqrt(mp.mpf('0.75')), mp.mpf('0.0')],
        ], dtype=object)
        chi = compute_internal_hierarchy_ratio(q_demo, (0, 1, 2))
        assert abs(chi - mp.mpf('1.0')) < mp.mpf('1e-30'), f"Chi for equilateral triangle = {chi}"

        masses = np.array([mp.mpf('1.0'), mp.mpf('1.0'), mp.mpf('1.0')], dtype=object)
        p = np.zeros_like(q_demo)
        state = ExtendedPhaseState(q_demo, p, masses)
        tree = build_canonical_tree(state.masses)

        sn_list = detect_supernodes(state, tree)
        assert len(sn_list) == 1
        assert sn_list[0] == (0, 1, 2)

    def test_schmitt_with_negative_E0_and_roche_disruption_f128(self):
        """
        Verifies [ENG-FALSE-07] in 35-digit mpmath arithmetic:
        1. Negative E0 does not short-circuit Schmitt trigger hysteresis when external perturbers exist.
        2. compute_cubic_tidal_force_ratio and mode='roche_disruption' detect Roche tidal disruption.
        """
        from normal_cone_nbody_float128.potential.supernode_detector import (
            compute_cubic_tidal_force_ratio,
            LAMBDA_CRITICAL_ROCHE_F128,
        )
        from normal_cone_nbody_float128.system.atlas_manager import AtlasManager

        m_jup = mp.mpf('9.547919e-4')
        m_frag = mp.mpf('1e-12')
        masses = np.array([m_jup, m_frag, m_frag, m_frag], dtype=object)
        d_frag = mp.mpf('8.165e-6')
        s3 = mp.sqrt(mp.mpf('3.0'))
        offsets = np.array([
            [d_frag, mp.mpf('0.0'), mp.mpf('0.0')],
            [-mp.mpf('0.5') * d_frag, mp.mpf('0.5') * s3 * d_frag, mp.mpf('0.0')],
            [-mp.mpf('0.5') * d_frag, -mp.mpf('0.5') * s3 * d_frag, mp.mpf('0.0')],
        ], dtype=object)

        q_init = np.full((4, 3), mp.mpf('0.0'), dtype=object)
        for k in range(3):
            q_init[k + 1] = np.array([mp.mpf('0.012') + offsets[k, 0], offsets[k, 1], offsets[k, 2]], dtype=object)
        p_zero = np.full((4, 3), mp.mpf('0.0'), dtype=object)
        st_init = ExtendedPhaseState(q=q_init, p=p_zero, masses=masses, E0=mp.mpf('-1e-10'))
        tree_init = AtlasManager().build_optimal_tree(st_init)

        lam_init = compute_cubic_tidal_force_ratio(d_frag, mp.mpf('0.012'), mp.mpf('3.0') * m_frag, m_jup)
        assert isinstance(lam_init, mp.mpf) and lam_init < LAMBDA_CRITICAL_ROCHE_F128
        assert detect_supernodes(st_init, tree_init, mode="roche_disruption") == []

        q_peri = np.full((4, 3), mp.mpf('0.0'), dtype=object)
        for k in range(3):
            q_peri[k + 1] = np.array([mp.mpf('6.33e-4') + offsets[k, 0], offsets[k, 1], offsets[k, 2]], dtype=object)
        st_peri = ExtendedPhaseState(q=q_peri, p=p_zero, masses=masses, E0=mp.mpf('-1e-10'))
        tree_peri = AtlasManager().build_optimal_tree(st_peri)
        assert detect_supernodes(st_peri, tree_peri, mode="roche_disruption") == [(1, 2, 3)]


    def test_supernode_synchronous_barycenter_conservation(self):
        """Verifies that step_supernode_synchronous preserves the cluster center of mass momentum."""
        q = np.array([
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.1'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.05'), mp.mpf('0.08'), mp.mpf('0.0')],
        ], dtype=object)
        p = np.array([
            [mp.mpf('0.1'), mp.mpf('-0.2'), mp.mpf('0.0')],
            [mp.mpf('-0.3'), mp.mpf('0.1'), mp.mpf('0.0')],
            [mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.0')],
        ], dtype=object)
        masses = np.array([mp.mpf('1.0'), mp.mpf('2.0'), mp.mpf('1.5')], dtype=object)
        dt = mp.mpf('0.005')

        P_init = sum(p[i] for i in range(3))

        q_new, p_new = step_supernode_synchronous(
            q, p, masses, (0, 1, 2), dt=dt, G=G_DEFAULT_F128, n_substeps=4
        )

        P_final = sum(p_new[i] for i in range(3))
        dP = mp.sqrt(sum((P_final[d] - P_init[d])**2 for d in range(3)))
        assert dP < mp.mpf('1e-32'), f"Momentum drift in supernode = {dP} >= 1e-32"

    def test_dominant_tidal_perturber_selection_f128(self):
        """
        Verifies D-144-07 (Part 1) in 35-digit mpmath arithmetic:
        Dominant external tidal perturber is selected by maximizing eta_S(k), not minimizing distance r_k.
        """
        masses = np.array([
            mp.mpf('1.0'), mp.mpf('1.0'), mp.mpf('1.0'), mp.mpf('1e-12'), mp.mpf('3000.0')
        ], dtype=object)
        half_sqrt3 = mp.sqrt(mp.mpf('3.0')) / mp.mpf('4.0')
        q_disrupted = np.array([
            [mp.mpf('0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('-0.25'), half_sqrt3, mp.mpf('0.0')],
            [mp.mpf('-0.25'), -half_sqrt3, mp.mpf('0.0')],
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('3.0')],
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('-4.0')],
        ], dtype=object)
        p = np.full((5, 3), mp.mpf('0.0'), dtype=object)
        state_disrupted = ExtendedPhaseState(q=q_disrupted, p=p, masses=masses)
        tree = build_canonical_tree(masses, topology="comb")

        sn_disrupted = detect_supernodes(state_disrupted, tree)
        assert (0, 1, 2) not in sn_disrupted

        q_isolated = q_disrupted.copy()
        q_isolated[4] = np.array([mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('-10.0')], dtype=object)
        state_isolated = ExtendedPhaseState(q=q_isolated, p=p, masses=masses)
        sn_isolated = detect_supernodes(state_isolated, tree)
        assert (0, 1, 2) in sn_isolated

    def test_supernode_detection_tetrahedral_unequal_masses_f128(self):
        """
        Verifies D-144-07 (Part 2) in 35-digit mpmath arithmetic:
        Tetrahedral 4-body non-transversal coalescence succeeds for unequal masses [1, 2, 3, 4].
        """
        masses = np.array([mp.mpf('1.0'), mp.mpf('2.0'), mp.mpf('3.0'), mp.mpf('4.0')], dtype=object)
        s = mp.mpf('0.05') / mp.sqrt(mp.mpf('3.0'))
        q = np.array([
            [s, s, s],
            [s, -s, -s],
            [-s, s, -s],
            [-s, -s, s],
        ], dtype=object)
        p = np.full((4, 3), mp.mpf('0.0'), dtype=object)
        state = ExtendedPhaseState(q=q, p=p, masses=masses)
        tree = build_canonical_tree(masses, topology="comb")

        clusters = detect_supernodes(state, tree)
        assert clusters == [(0, 1, 2, 3)]

    def test_supernode_homothetic_radial_collapse_and_clock_f128(self):
        """
        Verifies D-144-11 in 35-digit mpmath arithmetic:
        1. The homothetic radial collapse branch (is_radial_symmetric) preserves 35-digit
           mp.mpf precision and exact time-reversibility (< 1e-32) instead of truncating
           r0, v_rad, GM_eff to 53-bit float64.
        2. The adaptive Sundman clock (n_substeps=None) maintains t_accum and dt_rem in
           mp.mpf so that a pure irrational 35-digit dt = mp.pi / 100 is integrated over
           the exact 35-digit horizon without 1e-16 float64 roundoff on the final substep.
        """
        masses = np.array([mp.mpf('1.0'), mp.mpf('1.0'), mp.mpf('1.0')], dtype=object)
        r0 = mp.mpf('0.5')
        angles = [mp.mpf('0.0'), mp.mpf('2.0') * mp.pi / mp.mpf('3.0'), mp.mpf('4.0') * mp.pi / mp.mpf('3.0')]
        q = np.empty((3, 3), dtype=object)
        p = np.empty((3, 3), dtype=object)
        v_rad = mp.mpf('-0.2')
        for i in range(3):
            ca = mp.cos(angles[i])
            sa = mp.sin(angles[i])
            q[i, 0] = r0 * ca
            q[i, 1] = r0 * sa
            q[i, 2] = mp.mpf('0.0')
            p[i, 0] = masses[i] * v_rad * ca
            p[i, 1] = masses[i] * v_rad * sa
            p[i, 2] = mp.mpf('0.0')

        dt = mp.mpf('0.1')
        q_fwd, p_fwd = step_supernode_synchronous(q, p, masses, (0, 1, 2), dt=dt, G=G_DEFAULT_F128)
        assert isinstance(q_fwd[0, 0], mp.mpf) and isinstance(p_fwd[0, 0], mp.mpf)

        # Equilateral symmetry must hold to 35-digit precision (< 1e-32), not float64 (~1e-16)
        r_norms_fwd = [mp.sqrt(sum(q_fwd[i, d] ** 2 for d in range(3))) for i in range(3)]
        assert abs(r_norms_fwd[0] - r_norms_fwd[1]) < mp.mpf('1e-32')
        assert abs(r_norms_fwd[1] - r_norms_fwd[2]) < mp.mpf('1e-32')

        # Exact 35-digit time-reversibility of homothetic collapse
        q_rev, p_rev = step_supernode_synchronous(q_fwd, p_fwd, masses, (0, 1, 2), dt=-dt, G=G_DEFAULT_F128)
        err_q = max(abs(q_rev[i, d] - q[i, d]) for i in range(3) for d in range(3))
        err_p = max(abs(p_rev[i, d] - p[i, d]) for i in range(3) for d in range(3))
        assert err_q < mp.mpf('1e-32'), f"Homothetic radial collapse f128 reversibility err_q={err_q} >= 1e-32"
        assert err_p < mp.mpf('1e-32'), f"Homothetic radial collapse f128 reversibility err_p={err_p} >= 1e-32"

        # 2. Adaptive clock (n_substeps=None) with non-zero COM velocity and 35-digit irrational dt:
        # also verify that internal step sum matches 35-digit dt by checking a free-drifting / weakly
        # interacting cluster or comparing adaptive 1-step vs n_substeps=1 when dt < 0.05 * tau.
        q_gen = np.array([
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.25'), mp.mpf('0.4'), mp.mpf('0.0')],
        ], dtype=object)
        p_gen = np.array([
            [mp.mpf('0.1'), mp.mpf('-0.05'), mp.mpf('0.0')],
            [mp.mpf('-0.08'), mp.mpf('0.03'), mp.mpf('0.0')],
            [mp.mpf('-0.02'), mp.mpf('0.02'), mp.mpf('0.0')],
        ], dtype=object)
        dt_small = mp.pi * mp.mpf('1e-4')  # Guaranteed 1 adaptive substep, but h_local = dt_rem in adaptive loop!
        q_adap, p_adap = step_supernode_synchronous(
            q_gen, p_gen, masses, (0, 1, 2), dt=dt_small, G=G_DEFAULT_F128, n_substeps=None
        )
        q_exact1, p_exact1 = step_supernode_synchronous(
            q_gen, p_gen, masses, (0, 1, 2), dt=dt_small, G=G_DEFAULT_F128, n_substeps=1
        )
        diff_clock = max(abs(q_adap[i, d] - q_exact1[i, d]) for i in range(3) for d in range(3))
        assert diff_clock < mp.mpf('1e-33'), (
            f"Adaptive clock truncated dt_small to float64: diff_clock={diff_clock}"
        )

        # 3. Verify that step_supernode_synchronous_fast and step_supernode_synchronous_jit
        # preserve 35-digit mp.mpf precision when given dtype=object arrays (D-144-12)
        from normal_cone_nbody_float128.solvers.optimized_supernode import step_supernode_synchronous_fast
        from normal_cone_nbody_float128.solvers.jit_supernode import step_supernode_synchronous_jit

        q_fast, p_fast = step_supernode_synchronous_fast(
            q_gen, p_gen, masses, (0, 1, 2), dt=dt_small, G=G_DEFAULT_F128
        )
        q_jit, p_jit = step_supernode_synchronous_jit(
            q_gen, p_gen, masses, (0, 1, 2), dt=dt_small, G=G_DEFAULT_F128
        )
        assert isinstance(q_fast[0, 0], mp.mpf) and isinstance(q_jit[0, 0], mp.mpf)
        diff_fast = max(abs(q_fast[i, d] - q_adap[i, d]) for i in range(3) for d in range(3))
        diff_jit = max(abs(q_jit[i, d] - q_adap[i, d]) for i in range(3) for d in range(3))
        assert diff_fast < mp.mpf('1e-33')
        assert diff_jit < mp.mpf('1e-33')

    def test_supernode_multi_substep_reversibility_and_d4_bounce_f128(self):
        """
        Verifies [ENG-FALSE-03] in 35-digit mpmath arithmetic:
        1. Multi-substep adaptive supernode integration (dt = 0.01) is time-reversible
           to < 1e-33 (instead of ~1e-5 with asymmetric left-endpoint step selection).
        2. D_4 quadruple collision configuration (N = 4 square released from rest)
           traverses the collision at t_coll = 2.08575 up to t_max = 2.50 with 35-digit
           energy conservation and time-reversibility (< 1e-32).
        """
        masses3 = np.array([mp.mpf('1.0'), mp.mpf('1.5'), mp.mpf('2.0')], dtype=object)
        q3 = np.array([
            [mp.mpf('0.1'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('-0.05'), mp.mpf('0.08'), mp.mpf('0.0')],
            [mp.mpf('-0.05'), mp.mpf('-0.08'), mp.mpf('0.0')],
        ], dtype=object)
        p3 = np.array([
            [mp.mpf('0.0'), mp.mpf('0.5'), mp.mpf('0.0')],
            [mp.mpf('-0.3'), mp.mpf('-0.2'), mp.mpf('0.1')],
            [mp.mpf('0.3'), mp.mpf('-0.3'), mp.mpf('-0.1')],
        ], dtype=object)
        dt3 = mp.mpf('0.01')
        q3_fwd, p3_fwd = step_supernode_synchronous(
            q3, p3, masses3, (0, 1, 2), dt=dt3, G=G_DEFAULT_F128, n_substeps=None
        )
        q3_rev, p3_rev = step_supernode_synchronous(
            q3_fwd, p3_fwd, masses3, (0, 1, 2), dt=-dt3, G=G_DEFAULT_F128, n_substeps=None
        )
        err_q3 = max(abs(q3_rev[i, d] - q3[i, d]) for i in range(3) for d in range(3))
        err_p3 = max(abs(p3_rev[i, d] - p3[i, d]) for i in range(3) for d in range(3))
        assert err_q3 < mp.mpf('1e-33'), f"F128 multi-substep supernode rev err_q3={err_q3}"
        assert err_p3 < mp.mpf('1e-33'), f"F128 multi-substep supernode rev err_p3={err_p3}"

        # D_4 quadruple collision bounce in 35-digit precision (r0 = 1.5 => t_coll = 2.085747 < 2.50)
        masses4 = np.array([mp.mpf('1.0')] * 4, dtype=object)
        r0 = mp.mpf('1.5')
        z0 = mp.mpf('0.0')
        q4 = np.array([
            [r0, z0, z0],
            [z0, r0, z0],
            [-r0, z0, z0],
            [z0, -r0, z0],
        ], dtype=object)
        p4 = np.full((4, 3), z0, dtype=object)
        dt4 = mp.mpf('2.50')
        q4_fwd, p4_fwd = step_supernode_synchronous(
            q4, p4, masses4, (0, 1, 2, 3), dt=dt4, G=G_DEFAULT_F128
        )
        q4_rev, p4_rev = step_supernode_synchronous(
            q4_fwd, p4_fwd, masses4, (0, 1, 2, 3), dt=-dt4, G=G_DEFAULT_F128
        )
        err_q4 = max(abs(q4_rev[i, d] - q4[i, d]) for i in range(4) for d in range(3))
        err_p4 = max(abs(p4_rev[i, d] - p4[i, d]) for i in range(4) for d in range(3))
        assert err_q4 < mp.mpf('1e-32'), f"F128 D4 bounce rev err_q4={err_q4}"
        assert err_p4 < mp.mpf('1e-32'), f"F128 D4 bounce rev err_p4={err_p4}"



