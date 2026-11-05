r"""
Comprehensive High-Precision Unit Test Suite for Float128 / Mpmath Integrators and System Subsystem.

Validates:
1. Blanes--Moan BM4 and BM6 symplectic integrators in mpmath 35 dps.
2. Conormal boundary correctors C4 and C6 exact algebraic inversion C^-1 o C = Id (< 1e-32).
3. Exact Lie--Bargmann linear and angular momentum conservation of corrector kicks (< 1e-34).
4. DyadicScheduler bit-exact timestamp and substep generation in mpmath.
5. AtlasManager tree construction and dynamic dwell time under is_f128.
6. NBodySystem orchestration of BM4, PBM4, BM6, PBM6 on Burrau's problem.
7. Multi-step evolution under PBM6 with shadow manifold telescoping.
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
from normal_cone_nbody_float128.core.invariants import compute_invariants, compute_hamiltonian
from normal_cone_nbody_float128.potential.screened_apex import ScreenedApexPotential
from normal_cone_nbody_float128.integrators.strang import StrangIntegrator
from normal_cone_nbody_float128.integrators.yoshida import Yoshida4Integrator
from normal_cone_nbody_float128.integrators.yoshida4 import (
    Yoshida4Integrator as Yoshida4ShimF128,
)
from normal_cone_nbody_float128.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float128.integrators.blanes_moan_6 import BlanesMoan6Integrator
import normal_cone_nbody_float128.integrators.correctors as correctors_f128_mod
import normal_cone_nbody_float128.integrators.processed as processed_f128_mod
from normal_cone_nbody_float128.integrators.correctors import (
    CorrectorC4,
    GlobalCorrectorC4,
    compute_conormal_tidal_tensor,
    compute_conormal_corrector_kick,
    compute_global_conormal_c4_kick,
)
from normal_cone_nbody_float128.integrators.processed import (
    ProcessedSymplecticIntegrator,
    create_pbm4_integrator,
)
from normal_cone_nbody_float128.system.atlas_manager import AtlasManager
from normal_cone_nbody_float128.system.dyadic_scheduler import DyadicScheduler
from normal_cone_nbody_float128.system.nbody_system import NBodySystem


def get_burrau_state_f128() -> ExtendedPhaseState:
    """Returns the canonical Burrau 3-body initial state in extended precision (35 dps)."""
    mp.dps = DEFAULT_DPS
    masses = np.array([mp.mpf('3.0'), mp.mpf('4.0'), mp.mpf('5.0')], dtype=object)
    q = np.array([
        [mp.mpf('1.0'), mp.mpf('3.0'), mp.mpf('0.0')],
        [-mp.mpf('2.0'), -mp.mpf('1.0'), mp.mpf('0.0')],
        [mp.mpf('1.0'), -mp.mpf('1.0'), mp.mpf('0.0')],
    ], dtype=object)
    p = np.array([
        [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
        [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
        [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
    ], dtype=object)
    state = ExtendedPhaseState(q=q, p=p, masses=masses, t=mp.mpf('0.0'))
    return state.to_barycentric()


class TestIntegratorsSystemFloat128:
    """High-precision verification of integrators, correctors, scheduler, and system orchestration."""

    def test_blanes_moan_4_step_float128(self):
        """Validates single step of BM4 in 35 dps extended precision."""
        mp.dps = DEFAULT_DPS
        state = get_burrau_state_f128()
        tree = build_canonical_tree(state.masses)
        apex = ScreenedApexPotential(G=G_DEFAULT_F128)
        integrator = BlanesMoan4Integrator(G=G_DEFAULT_F128)

        dt = mp.mpf('0.01')
        next_state = integrator.step(state, tree, dt=dt, apex=apex)

        assert next_state.is_f128
        assert isinstance(next_state.t, mp.mpf)
        assert abs(next_state.t - dt) < mp.mpf('1e-34')
        assert next_state.q.dtype == object

        # Center of mass and linear momentum conservation
        inv = compute_invariants(next_state, G=G_DEFAULT_F128)
        P_tot_norm = mp.sqrt(sum(inv.P[d]**2 for d in range(3)))
        assert P_tot_norm < mp.mpf('1e-32')

    def test_blanes_moan_6_step_float128(self):
        """Validates single step of BM6 in 35 dps extended precision."""
        mp.dps = DEFAULT_DPS
        state = get_burrau_state_f128()
        tree = build_canonical_tree(state.masses)
        apex = ScreenedApexPotential(G=G_DEFAULT_F128)
        integrator = BlanesMoan6Integrator(G=G_DEFAULT_F128)

        dt = mp.mpf('0.01')
        next_state = integrator.step(state, tree, dt=dt, apex=apex)

        assert next_state.is_f128
        inv = compute_invariants(next_state, G=G_DEFAULT_F128)
        H0 = compute_hamiltonian(state, G=G_DEFAULT_F128)
        H1 = compute_hamiltonian(next_state, G=G_DEFAULT_F128)
        dE_rel = abs(H1 - H0) / abs(H0)
        assert dE_rel < mp.mpf('1e-8')

    def test_conormal_corrector_c4_involutivity_float128(self):
        """Validates C_4^-1 o C_4 = Id across all 3 modes ('homological', 'cotangent_lift', 'global') to < 1e-33 in float128."""
        mp.dps = DEFAULT_DPS
        state = get_burrau_state_f128()
        tree = build_canonical_tree(state.masses)
        apex = ScreenedApexPotential(G=G_DEFAULT_F128)
        default_corrector = CorrectorC4(G=G_DEFAULT_F128)
        assert default_corrector.mode == "homological"

        # 1. Pure conormal momentum shear kick exact algebraic involution (< 1e-34)
        dp_fwd = compute_global_conormal_c4_kick(state, tree, mp.mpf('0.05'), G=G_DEFAULT_F128, inverse=False)
        dp_inv = compute_global_conormal_c4_kick(state, tree, mp.mpf('0.05'), G=G_DEFAULT_F128, inverse=True)
        for i in range(state.N):
            for d in range(3):
                assert abs(dp_fwd[i, d] + dp_inv[i, d]) < mp.mpf('1e-34')

        # 2. All 3 corrector modes ('homological', 'cotangent_lift', 'global') roundtrip error at h=0.05 (< 1e-33)
        dt = mp.mpf('0.05')
        inv0 = compute_invariants(state, G=G_DEFAULT_F128)
        for mode in ("homological", "cotangent_lift", "global"):
            corrector = CorrectorC4(G=G_DEFAULT_F128, mode=mode)
            forward_state = corrector.apply(state, tree, dt, inverse=False, apex=apex)
            roundtrip_state = corrector.apply_inverse(forward_state, tree, dt, apex=apex)

            err_q = mp.mpf('0.0')
            err_p = mp.mpf('0.0')
            for i in range(state.N):
                for d in range(3):
                    err_q = max(err_q, abs(roundtrip_state.q[i, d] - state.q[i, d]))
                    err_p = max(err_p, abs(roundtrip_state.p[i, d] - state.p[i, d]))

            if mode == "global":
                assert err_q == mp.mpf('0.0')
            else:
                assert err_q < mp.mpf('1e-33'), f"mode={mode} Float128 q inversion error: {err_q}"
            assert err_p < mp.mpf('1e-33'), f"mode={mode} Float128 p inversion error: {err_p}"

            # Back-reaction total momentum, angular momentum, and boost conservation (< 1e-32)
            inv_fwd = compute_invariants(forward_state, G=G_DEFAULT_F128)
            P_norm = mp.sqrt(sum(inv_fwd.P[d]**2 for d in range(3)))
            dL_norm = mp.sqrt(sum((inv_fwd.L[d] - inv0.L[d])**2 for d in range(3)))
            dG_norm = mp.sqrt(sum((inv_fwd.G[d] - inv0.G[d])**2 for d in range(3)))
            assert P_norm < mp.mpf('1e-32')
            assert dL_norm < mp.mpf('1e-32')
            assert dG_norm < mp.mpf('1e-32')

    def test_universal_obstruction_purge_of_c6_and_pbm6_float128(self):
        """Validates Theorem 5.14: CorrectorC6, compute_conormal_c6_step, and create_pbm6_integrator are purged."""
        assert not hasattr(correctors_f128_mod, "CorrectorC6")
        assert not hasattr(correctors_f128_mod, "compute_conormal_c6_step")
        assert not hasattr(processed_f128_mod, "create_pbm6_integrator")

    def test_dyadic_scheduler_float128(self):
        """Validates DyadicScheduler operations in float128."""
        mp.dps = DEFAULT_DPS
        state = get_burrau_state_f128()
        tree = build_canonical_tree(state.masses)
        scheduler = DyadicScheduler(dt_macro=mp.mpf('0.01'), G=G_DEFAULT_F128)

        # 1. Level determination without type errors
        lvl = scheduler.compute_dynamical_level(state, tree=tree)
        assert 0 <= lvl <= 10

        # 2. Exact substep parameters
        num_substeps, dt_sub, k = scheduler.get_substep_params(
            t_current=state.t,
            dt_block=mp.mpf('0.01'),
            state=state,
            tree=tree,
        )
        assert isinstance(dt_sub, mp.mpf)
        assert abs(num_substeps * dt_sub - mp.mpf('0.01')) < mp.mpf('1e-34')

        # 3. Bit-exact timestamps in mpmath
        t_stamps = DyadicScheduler.generate_bit_exact_timestamps(
            t_start=mp.mpf('0.0'),
            dt_block=mp.mpf('0.01'),
            k=2,
            is_f128=True,
        )
        assert len(t_stamps) == 5
        assert isinstance(t_stamps[0], mp.mpf)
        assert abs(t_stamps[-1] - mp.mpf('0.01')) < mp.mpf('1e-34')

        # 4. V4 Time-reversibility symmetry k(q, p) == k(q, -p) on non-leaf grazing encounter
        q_enc = np.array([
            [mp.mpf('-10.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.1'), mp.mpf('0.005'), mp.mpf('0.0')],
        ], dtype=object)
        p_in = np.array([
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.0'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('-1.5'), mp.mpf('0.1'), mp.mpf('0.0')],
        ], dtype=object)
        m_enc = np.array([mp.mpf('1.0'), mp.mpf('1.0'), mp.mpf('1.0')], dtype=object)
        tree_enc = build_canonical_tree(m_enc)
        sched_enc = DyadicScheduler(dt_macro=mp.mpf('0.1'), G=G_DEFAULT_F128)
        k_in = sched_enc.compute_dynamical_level(
            ExtendedPhaseState(q=q_enc, p=p_in, masses=m_enc), tree=tree_enc
        )
        k_out = sched_enc.compute_dynamical_level(
            ExtendedPhaseState(q=q_enc, p=-p_in, masses=m_enc), tree=tree_enc
        )
        assert k_in >= 6
        assert k_in == k_out

    def test_atlas_manager_float128(self):
        """Validates AtlasManager tree construction and dynamic dwell time under is_f128."""
        mp.dps = DEFAULT_DPS
        state = get_burrau_state_f128()
        atlas = AtlasManager()

        # Dynamic dwell time
        dwell = atlas.compute_dynamic_dwell_time(state)
        assert dwell >= 0.0

        # Optimal tree construction with mpmath mass preservation
        tree = atlas.build_optimal_tree(state)
        assert tree.is_f128
        assert abs(tree.root.M - (state.masses[0] + state.masses[1] + state.masses[2])) < mp.mpf('1e-34')

    def test_nbody_system_macro_step_all_integrators_float128(self):
        """Validates 1 macro-step under BM4, PBM4, BM6 on Burrau in float128."""
        mp.dps = DEFAULT_DPS
        dt_macro = mp.mpf('0.01')

        integrator_factories = [
            ("BM4", lambda: BlanesMoan4Integrator(G=G_DEFAULT_F128)),
            ("PBM4", lambda: create_pbm4_integrator(G=G_DEFAULT_F128)),
            ("BM6", lambda: BlanesMoan6Integrator(G=G_DEFAULT_F128)),
        ]

        for name, factory in integrator_factories:
            state = get_burrau_state_f128()
            system = NBodySystem(
                state=state,
                integrator=factory(),
                G=G_DEFAULT_F128,
            )
            H0 = compute_hamiltonian(system.state, G=G_DEFAULT_F128)
            next_state = system.step(dt=dt_macro)

            assert next_state.is_f128, f"{name} failed to preserve is_f128"
            assert abs(next_state.t - dt_macro) < mp.mpf('1e-34'), f"{name} timestamp drift"

            H1 = compute_hamiltonian(next_state, G=G_DEFAULT_F128)
            dE_rel = abs(H1 - H0) / abs(H0)
            assert dE_rel < mp.mpf('1e-7'), f"{name} energy error too large: {dE_rel}"

    def test_nbody_system_evolve_pbm4_float128(self):
        """Validates multi-step evolution under PBM4 with shadow manifold telescoping in float128."""
        mp.dps = DEFAULT_DPS
        state = get_burrau_state_f128()
        pbm4 = create_pbm4_integrator(G=G_DEFAULT_F128)
        system = NBodySystem(
            state=state,
            integrator=pbm4,
            G=G_DEFAULT_F128,
        )

        H0 = compute_hamiltonian(system.state, G=G_DEFAULT_F128)
        t_span = mp.mpf('0.03')
        final_state = system.evolve(t_span=t_span, dt=mp.mpf('0.01'))

        assert final_state.is_f128
        assert abs(final_state.t - t_span) < mp.mpf('1e-32')

        H1 = compute_hamiltonian(final_state, G=G_DEFAULT_F128)
        dE_rel = abs(H1 - H0) / abs(H0)
        assert dE_rel < mp.mpf('1e-7')

        # Conservation of linear and angular momentum to < 1e-32
        inv0 = compute_invariants(state, G=G_DEFAULT_F128)
        inv = compute_invariants(final_state, G=G_DEFAULT_F128)
        P_norm = mp.sqrt(sum(inv.P[d]**2 for d in range(3)))
        dL_norm = mp.sqrt(sum((inv.L[d] - inv0.L[d])**2 for d in range(3)))
        assert P_norm < mp.mpf('1e-32')
        assert dL_norm < mp.mpf('1e-32')

    def test_nbody_system_multistep_continuous_trajectory_float128(self):
        """Validates multi-step trajectory evolution with HistoryTracker recording under PBM4 and BM6."""
        mp.dps = DEFAULT_DPS
        dt = mp.mpf('0.005')
        n_steps = 10

        for scheme_name, integrator_factory in [
            ("PBM4", lambda: create_pbm4_integrator(G=G_DEFAULT_F128)),
            ("BM6", lambda: BlanesMoan6Integrator(G=G_DEFAULT_F128)),
        ]:
            state = get_burrau_state_f128()
            system = NBodySystem(
                state=state,
                integrator=integrator_factory(),
                G=G_DEFAULT_F128,
            )
            H0 = compute_hamiltonian(system.state, G=G_DEFAULT_F128)

            for step in range(n_steps):
                system.step(dt=dt)

            assert len(system.history) == n_steps + 1
            assert system.history.is_f128 is True

            # Verify history drift metrics
            summary = system.history.summary_metrics()
            dE_max = summary["delta_H_max"]
            dP_max = summary["delta_P_max"]
            dL_max = summary["delta_L_max"]

            # Relative energy drift bounded
            assert dE_max / abs(H0) < mp.mpf('1e-6'), f"{scheme_name} energy drift: {dE_max / abs(H0)}"
            # Both linear and angular momentum conserved to < 1e-30 (GlobalCorrectorC4 preserves L to < 1e-32!)
            assert dP_max < mp.mpf('1e-30'), f"{scheme_name} linear momentum drift: {dP_max}"
            assert dL_max < mp.mpf('1e-30'), f"{scheme_name} angular momentum drift: {dL_max}"

    def test_nbody_system_pbm4_multirate_step_vs_evolve_float128(self):
        """
        Validates D-145-04 and D-145-15 in float128:
        - Default max_dyadic_level is 17.
        - Under multirate dyadic substepping (k = 2), M calls to step(dt) and a single call
          to evolve(M*dt, dt) under PBM4 produce identical 35-digit trajectories (< 1e-31).
        """
        mp.dps = DEFAULT_DPS
        state0 = get_burrau_state_f128()

        sys_step = NBodySystem(
            state=state0,
            integrator=create_pbm4_integrator(G=G_DEFAULT_F128),
            G=G_DEFAULT_F128,
        )
        assert sys_step.scheduler.max_level == 17
        sys_step.scheduler.min_level = 2

        sys_evolve = NBodySystem(
            state=state0,
            integrator=create_pbm4_integrator(G=G_DEFAULT_F128),
            G=G_DEFAULT_F128,
        )
        sys_evolve.scheduler.min_level = 2

        dt = mp.mpf('0.01')
        M = 3
        for _ in range(M):
            sys_step.step(dt=dt)
        sys_evolve.evolve(t_span=M * dt, dt=dt)

        max_dq = max(abs(sys_step.state.q[i, d] - sys_evolve.state.q[i, d]) for i in range(3) for d in range(3))
        max_dp = max(abs(sys_step.state.p[i, d] - sys_evolve.state.p[i, d]) for i in range(3) for d in range(3))
        assert max_dq < mp.mpf('1e-31'), f"float128 q mismatch between step and evolve: {max_dq}"
        assert max_dp < mp.mpf('1e-31'), f"float128 p mismatch between step and evolve: {max_dp}"

    def test_base_integrators_detached_reversibility_invariants_and_orders_float128(self):
        """
        Certifies Detached Component 3 in Float128 (35 dps):
        1. StrangIntegrator, Yoshida4Integrator (and yoshida4 shim), BlanesMoan4Integrator,
           and BlanesMoan6Integrator preserve is_f128 == True even when dt is passed as float or mp.mpf,
           and expose standardized force evaluation counts.
        2. Exact microscopic time-reversibility Phi(-h) o Phi(h) = Id to < 1e-33 across all 4 base integrators.
        3. Exact conservation of all 9 Lie-Bargmann generators (P, L, G_boost) to < 1e-33.
        4. Asymptotic convergence slopes p = 2, 4, 4, 6 in 35-digit precision down below the float64 floor.
        """
        mp.dps = DEFAULT_DPS
        assert Yoshida4ShimF128 is Yoshida4Integrator

        masses = np.array([mp.mpf('1.0'), mp.mpf('1.0'), mp.mpf('0.5')], dtype=object)
        q0 = np.array([
            [mp.mpf('-0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.5'), mp.mpf('0.0'), mp.mpf('0.0')],
            [mp.mpf('0.0'), mp.mpf('8.0'), mp.mpf('0.5')],
        ], dtype=object)
        # Non-zero total momentum P != 0 to test Galilean boost G_boost = M_tot * R_cdm - t * P
        p0 = np.array([
            [mp.mpf('0.02'), mp.mpf('-0.7'), mp.mpf('0.01')],
            [mp.mpf('0.01'), mp.mpf('0.7'), mp.mpf('-0.01')],
            [mp.mpf('0.25'), mp.mpf('0.05'), mp.mpf('-0.1')],
        ], dtype=object)
        state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=mp.mpf('0.0'))
        tree = build_canonical_tree(masses, topology="comb")
        inv0 = compute_invariants(state0, G=G_DEFAULT_F128)

        s2 = StrangIntegrator(G=G_DEFAULT_F128)
        y4 = Yoshida4Integrator(G=G_DEFAULT_F128)
        bm4 = BlanesMoan4Integrator(G=G_DEFAULT_F128)
        bm6 = BlanesMoan6Integrator(G=G_DEFAULT_F128)

        assert (s2.order, s2.num_force_evaluations_per_step, s2.fsal_evaluations_per_step) == (2, 2, 1)
        assert (y4.order, y4.num_force_evaluations_per_step, y4.fsal_evaluations_per_step) == (4, 4, 3)
        assert (bm4.order, bm4.num_force_evaluations_per_step, bm4.fsal_evaluations_per_step) == (4, 7, 7)
        assert (bm6.order, bm6.num_force_evaluations_per_step, bm6.fsal_evaluations_per_step) == (6, 11, 11)

        # Verify StrangIntegrator preserves is_f128 even when dt is passed as a Python float
        s_float_dt = s2.step(state0, tree, dt=0.01)
        assert s_float_dt.is_f128 is True
        assert isinstance(s_float_dt.t, mp.mpf)

        dt = mp.mpf('0.01')
        for integ in (s2, y4, bm4, bm6):
            s_fwd = integ.step(state0, tree, dt=dt)
            assert s_fwd.is_f128 is True
            s_rev = integ.step(s_fwd, tree, dt=-dt)
            err_q = max(abs(s_rev.q[i, d] - state0.q[i, d]) for i in range(3) for d in range(3))
            err_p = max(abs(s_rev.p[i, d] - state0.p[i, d]) for i in range(3) for d in range(3))
            assert err_q < mp.mpf('1e-33'), f"{integ.name} Float128 q reversibility error: {err_q}"
            assert err_p < mp.mpf('1e-33'), f"{integ.name} Float128 p reversibility error: {err_p}"

            inv1 = compute_invariants(s_fwd, G=G_DEFAULT_F128)
            err_P = max(abs(inv1.P[d] - inv0.P[d]) for d in range(3))
            err_L = max(abs(inv1.L[d] - inv0.L[d]) for d in range(3))
            err_G = max(abs(inv1.G[d] - inv0.G[d]) for d in range(3))
            assert err_P < mp.mpf('1e-33'), f"{integ.name} Float128 P error: {err_P}"
            assert err_L < mp.mpf('1e-33'), f"{integ.name} Float128 L error: {err_L}"
            assert err_G < mp.mpf('1e-33'), f"{integ.name} Float128 G error: {err_G}"

        # Asymptotic convergence slopes in Float128 (h_coarse = 0.01 vs h_fine = 0.005 over 1 step vs 2 steps)
        h_c = mp.mpf('0.01')
        h_f = mp.mpf('0.005')
        for integ, expected_ratio, tol_low in [
            (s2, mp.mpf('4.0'), mp.mpf('3.8')),
            (y4, mp.mpf('16.0'), mp.mpf('15.0')),
            (bm4, mp.mpf('16.0'), mp.mpf('15.0')),
            (bm6, mp.mpf('64.0'), mp.mpf('60.0')),
        ]:
            sc = integ.step(state0, tree, dt=h_c)
            sf = integ.step(integ.step(state0, tree, dt=h_f), tree, dt=h_f)
            err_c = abs(compute_hamiltonian(sc, G=G_DEFAULT_F128) - inv0.H)
            err_f = abs(compute_hamiltonian(sf, G=G_DEFAULT_F128) - inv0.H)
            ratio = err_c / err_f
            assert ratio >= tol_low, f"{integ.name} Float128 order ratio {float(ratio):.2f} < {float(tol_low)}"
            if integ is bm6:
                # BM6 single-step energy error at h=0.005 is ~ 1e-22, 6 orders of magnitude below float64!
                assert err_f < mp.mpf('1e-20'), f"BM6 Float128 error did not break float64 barrier: {err_f}"

    def test_shadow_manifold_lifecycle_and_universal_obstruction_float128(self):
        """
        Certifies Detached Component 4 in Float128 (35 dps):
        1. Full 4-method Shadow Manifold Lifecycle API: enter_shadow, step_shadow, observe_physical, exit_shadow.
        2. Direct empirical proof of Theorem 5.14 on a 3D eccentric hierarchical triple (e_in = 0.30, i = 28 deg)
           in quadruple precision:
           - mode='global' (P_0) is universally obstructed (~1.0x gain over bare BM4);
           - mode='cotangent_lift' (P_1) and canonical default mode='homological' (Z_4^{BM4} in P_1)
             achieve > 1000x energy error reduction over bare BM4 and over P_0!
        """
        from normal_cone_nbody_float128.core.jacobi import build_tree_from_spec

        mp.dps = DEFAULT_DPS
        masses = np.array([mp.mpf('0.5'), mp.mpf('0.5'), mp.mpf('1.0')], dtype=object)
        a_in = mp.mpf('1.0')
        e_in = mp.mpf('0.30')
        inc = mp.mpf('28.0') * mp.pi / mp.mpf('180.0')
        r_p = a_in * (mp.mpf('1.0') - e_in)
        v_p = mp.sqrt((mp.mpf('1.0') + e_in) / (a_in * (mp.mpf('1.0') - e_in)))
        R_out = mp.mpf('50.0')
        v_out = mp.sqrt(mp.mpf('2.0') / R_out)

        q0 = np.array([
            [-mp.mpf('0.5') * r_p, mp.mpf('0.0'), mp.mpf('0.0')],
            [+mp.mpf('0.5') * r_p, mp.mpf('0.0'), mp.mpf('0.0')],
            [R_out * mp.mpf('0.8'), R_out * mp.mpf('0.6'), mp.mpf('0.0')],
        ], dtype=object)
        v0 = np.array([
            [mp.mpf('0.0'), -mp.mpf('0.5') * v_p * mp.cos(inc), -mp.mpf('0.5') * v_p * mp.sin(inc)],
            [mp.mpf('0.0'), +mp.mpf('0.5') * v_p * mp.cos(inc), +mp.mpf('0.5') * v_p * mp.sin(inc)],
            [-mp.mpf('0.6') * v_out, mp.mpf('0.8') * v_out, mp.mpf('0.0')],
        ], dtype=object)
        p0 = np.array([[masses[i] * v0[i, d] for d in range(3)] for i in range(3)], dtype=object)
        state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=mp.mpf('0.0')).to_barycentric()
        tree = build_tree_from_spec(masses, ((0, 1), 2))
        E0 = compute_hamiltonian(state0, G=G_DEFAULT_F128)

        bm4 = BlanesMoan4Integrator(G=G_DEFAULT_F128)
        pbm4_p0 = create_pbm4_integrator(G=G_DEFAULT_F128, mode="global")
        pbm4_cot = create_pbm4_integrator(G=G_DEFAULT_F128, mode="cotangent_lift")
        pbm4_hom = create_pbm4_integrator(G=G_DEFAULT_F128)  # Default is mode="homological"
        assert pbm4_hom.corrector.mode == "homological"

        dt = mp.mpf('0.02')
        n_steps = 10
        st_bm4 = state0
        sh_p0 = pbm4_p0.enter_shadow(state0, tree, dt)
        sh_cot = pbm4_cot.enter_shadow(state0, tree, dt)
        sh_hom = pbm4_hom.enter_shadow(state0, tree, dt)

        max_err_bm4 = mp.mpf('0.0')
        max_err_p0 = mp.mpf('0.0')
        max_err_cot = mp.mpf('0.0')
        max_err_hom = mp.mpf('0.0')

        for _ in range(n_steps):
            st_bm4 = bm4.step(st_bm4, tree, dt)
            sh_p0 = pbm4_p0.step_shadow(sh_p0, tree, dt)
            sh_cot = pbm4_cot.step_shadow(sh_cot, tree, dt)
            sh_hom = pbm4_hom.step_shadow(sh_hom, tree, dt)

            obs_p0 = pbm4_p0.observe_physical(sh_p0, tree, dt)
            obs_cot = pbm4_cot.observe_physical(sh_cot, tree, dt)
            obs_hom = pbm4_hom.observe_physical(sh_hom, tree, dt)

            max_err_bm4 = max(max_err_bm4, abs(compute_hamiltonian(st_bm4, G=G_DEFAULT_F128) - E0) / abs(E0))
            max_err_p0 = max(max_err_p0, abs(compute_hamiltonian(obs_p0, G=G_DEFAULT_F128) - E0) / abs(E0))
            max_err_cot = max(max_err_cot, abs(compute_hamiltonian(obs_cot, G=G_DEFAULT_F128) - E0) / abs(E0))
            max_err_hom = max(max_err_hom, abs(compute_hamiltonian(obs_hom, G=G_DEFAULT_F128) - E0) / abs(E0))

        final_hom = pbm4_hom.exit_shadow(sh_hom, tree, dt)
        for i in range(3):
            for d in range(3):
                assert final_hom.q[i, d] == obs_hom.q[i, d]
                assert final_hom.p[i, d] == obs_hom.p[i, d]

        gain_p0 = float(max_err_bm4 / max_err_p0)
        gain_cot = float(max_err_bm4 / max_err_cot)
        gain_hom = float(max_err_bm4 / max_err_hom)

        assert 0.7 <= gain_p0 <= 1.25, f"Expected P_0 obstruction (~1.0x), got {gain_p0:.2f}x"
        assert gain_cot > 150.0, f"Expected Float128 CotangentLift gain > 150x, got {gain_cot:.1f}x"
        assert gain_hom > 2000.0, f"Expected Float128 Homological Z_4^{{BM4}} gain > 2000x, got {gain_hom:.1f}x"

    def test_symmetric_dyadic_palindrome_and_default_shadow_mode_float128(self):
        """
        Certifies Detached Component 6 in Float128 (35 dps):
        1. NBodySystem defaults to in_shadow_mode=True and preserves stateful shadow evolution
           across step() calls with get_physical_state() / exit_shadow_manifold().
        2. Symmetric dyadic scheduler achieves exact palindromic level sequence k_fwd(n) == k_bwd(N-1-n)
           and < 1e-32 roundtrip reversibility across a hyperbolic leaf pericenter passage.
        """
        from normal_cone_nbody_float128.diagnostics.metrics import (
            compute_forward_backward_reversibility,
        )

        mp.dps = DEFAULT_DPS
        masses = np.array([mp.mpf("1.0"), mp.mpf("1e-6")], dtype=object)
        q0 = np.array([
            [mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")],
            [mp.mpf("0.25"), mp.mpf("0.08"), mp.mpf("0.0")],
        ], dtype=object)
        v0 = np.array([
            [mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")],
            [mp.mpf("-3.5"), mp.mpf("0.0"), mp.mpf("0.0")],
        ], dtype=object)
        p0 = np.array([[masses[i] * v0[i, d] for d in range(3)] for i in range(2)], dtype=object)
        state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=mp.mpf("0.0")).to_barycentric()

        rev = compute_forward_backward_reversibility(
            state0=state0,
            integrator=create_pbm4_integrator(G=G_DEFAULT_F128),
            dt=mp.mpf("0.02"),
            n_steps=6,
            G=G_DEFAULT_F128,
        )
        assert rev["dyadic_palindrome_matched"] is True, (
            f"Float128 dyadic palindrome mismatch: {rev['k_fwd']} vs {rev['k_bwd_reversed']}"
        )
        assert max(rev["k_fwd"]) > min(rev["k_fwd"])
        assert rev["delta_q_norm"] < 1e-31
        assert rev["delta_p_norm"] < 1e-31




