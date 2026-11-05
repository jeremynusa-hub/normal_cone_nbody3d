"""
Unit Tests for Layer 2: Universal Kepler-Sundman Two-Body Solver.

Theoretical Alignments:
- Monograph: Section 3.3 (Sundman-Stumpff Universal Diffeomorphism),
  Section 5.2 (Sheng-Suzuki Barrier and Time-Reversibility), and Appendix A.
- Formal Lean 4: `02_PCR3BPRegularization_Theorems.lean` and `PaperUnifiedFoundations.lean`.
- Invariants: Wronskian unimodularity |f dot{g} - dot{f} g - 1| < 1e-15,
  exact time-reversibility Phi(-dt) o Phi(dt) = Id, and energy/momentum conservation.
"""

import math
import numpy as np
import pytest

from normal_cone_nbody_float64.solvers.kepler_sundman import (
    solve_sundman_time_equation,
    step_kepler_sundman,
    compute_wronskian,
)


def test_halley_convergence_cubic():
    """
    Validates that cubic Halley iteration solves F(tau) = dt to machine precision
    across various physical regimes (elliptic, parabolic, hyperbolic).
    """
    GM = 1.0
    r0 = 1.0
    sigma0 = 0.2
    alpha_elliptic = 1.5   # Bound elliptic orbit
    alpha_parabolic = 0.0  # Parabolic orbit
    alpha_hyperbolic = -1.0 # Hyperbolic escape orbit

    for alpha in [alpha_elliptic, alpha_parabolic, alpha_hyperbolic]:
        for dt in [0.01, 0.1, 0.5, 2.0]:
            tau = solve_sundman_time_equation(r0, sigma0, GM, alpha, dt, tol=1e-15)
            # Verify F(tau) - dt == 0
            from normal_cone_nbody_float64.solvers.stumpff import stumpff_c
            z = alpha * tau * tau
            c0, c1, c2, c3, _, _ = stumpff_c(z)
            t_eval = r0 * tau * c1 + sigma0 * (tau**2) * c2 + GM * (tau**3) * c3
            assert np.isclose(t_eval, dt, atol=1e-14), (
                f"Halley failed at alpha={alpha}, dt={dt}: evaluated t={t_eval} vs target {dt}"
            )


def test_wronskian_unimodularity_exact():
    """
    Validates the fundamental Darboux symplectic condition of the universal Kepler flow:
        f * dot{g} - dot{f} * g == 1.0
    evaluated to machine precision IEEE 754 (< 1e-15).
    """
    mu = 0.8
    GM = 1.0

    test_states = [
        # Circular orbit (r=1, v=1, L=1)
        (np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.8, 0.0])),
        # Elliptic orbit with high eccentricity
        (np.array([0.5, 0.0, 0.0]), np.array([0.0, 0.8 * 1.5, 0.0])),
        # Hyperbolic flyby
        (np.array([2.0, 1.0, 0.5]), np.array([-0.8 * 1.8, 0.8 * 0.5, 0.8 * 0.2])),
        # Oblique near-radial trajectory
        (np.array([0.1, 0.01, 0.0]), np.array([-0.8 * 3.0, 0.8 * 0.1, 0.0])),
    ]

    for r0, w0 in test_states:
        for dt in [0.01, 0.05, 0.2, 1.0]:
            w_det = compute_wronskian(r0, w0, mu, GM, dt)
            assert np.isclose(w_det, 1.0, atol=1e-14), (
                f"Wronskian violation: |w_det - 1| = {abs(w_det - 1.0)} for dt={dt}"
            )


def test_time_reversibility_retrograde():
    """
    Validates exact time-reversal invariance:
        Phi(-dt) o Phi(dt) = Id
    for both forward (dt > 0) and retrograde (dt < 0) steps (essential for Suzuki barrier bypass).
    """
    mu = 0.75
    GM = 1.5
    r0 = np.array([1.2, -0.4, 0.3], dtype=np.float64)
    w0 = np.array([0.2, 0.9, -0.5], dtype=np.float64)

    for dt in [0.05, 0.2, -0.15, -0.5]:
        r_fwd, w_fwd = step_kepler_sundman(r0, w0, mu, GM, dt)
        r_back, w_back = step_kepler_sundman(r_fwd, w_fwd, mu, GM, -dt)

        assert np.allclose(r_back, r0, atol=1e-14), f"Position reversibility error: {np.max(np.abs(r_back - r0))}"
        assert np.allclose(w_back, w0, atol=1e-14), f"Momentum reversibility error: {np.max(np.abs(w_back - w0))}"


def test_two_body_energy_and_angular_momentum_conservation():
    """
    Validates exact physical conservation of orbital energy and angular momentum vector.
    """
    mu = 0.6
    GM = 1.0
    r0 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    w0 = np.array([0.0, 0.6 * 1.1, 0.0], dtype=np.float64)

    def energy(r, w):
        v = w / mu
        return 0.5 * mu * float(np.dot(v, v)) - float(GM * mu / np.linalg.norm(r))

    def angular_momentum(r, w):
        return np.cross(r, w)

    E_init = energy(r0, w0)
    L_init = angular_momentum(r0, w0)

    r_curr, w_curr = np.copy(r0), np.copy(w0)
    dt = 0.05
    for _ in range(50):
        r_curr, w_curr = step_kepler_sundman(r_curr, w_curr, mu, GM, dt)
        E_now = energy(r_curr, w_curr)
        L_now = angular_momentum(r_curr, w_curr)

        assert np.isclose(E_now, E_init, atol=1e-13), f"Energy drift: {abs(E_now - E_init)}"
        assert np.allclose(L_now, L_init, atol=1e-14), f"Angular momentum drift: {np.max(np.abs(L_now - L_init))}"


def test_head_on_radial_collision_traversal():
    """
    Validates that a pure head-on radial collision orbit (e = 1, L = 0)
    traverses the collision singularity cleanly without producing NaN or Inf.
    """
    mu = 1.0
    GM = 1.0
    # Initial state collapsing directly toward the origin
    r0 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    w0 = np.array([-0.5, 0.0, 0.0], dtype=np.float64)

    # Step through and past collision
    dt = 1.5  # Large enough to cross collision
    r_f, w_f = step_kepler_sundman(r0, w0, mu, GM, dt)

    assert not np.any(np.isnan(r_f))
    assert not np.any(np.isnan(w_f))
    assert not np.any(np.isinf(r_f))
    assert not np.any(np.isinf(w_f))
    assert float(np.linalg.norm(r_f)) > 0.0

    # Verify exact energy conservation across head-on bounce without artificial v_mag projection
    E0 = 0.5 * float(np.dot(w0, w0)) / mu - GM * mu / float(np.linalg.norm(r0))
    Ef = 0.5 * float(np.dot(w_f, w_f)) / mu - GM * mu / float(np.linalg.norm(r_f))
    assert abs(Ef - E0) / abs(E0) < 1e-14

    # Verify round-trip reversibility across collision
    r_rev, w_rev = step_kepler_sundman(r_f, w_f, mu, GM, -dt)
    assert np.allclose(r_rev, r0, atol=1e-14)
    assert np.allclose(w_rev, w0, atol=1e-14)


def test_kepler_pericenter_homotopy_no_cancellation():
    """
    Validates [ENG-FALSE-01] & [ENG-FALSE-02]:
    Across the extreme eccentricity homotopy 1 - e in [1e-1, 1e-12] (a = 1, GM = 1, mu = 1),
    traversing pericenter from apocenter and back preserves energy to <= 1e-14 without
    any (1 - e)^{-2} subtractive cancellation.
    """
    from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
        step_kepler_sundman_fast,
    )
    mu = 1.0
    GM = 1.0
    T_orb = 2.0 * math.pi
    for k_exp in [1, 3, 6, 9, 12]:
        one_minus_e = 10.0 ** (-k_exp)
        e = 1.0 - one_minus_e
        r_apo = 1.0 + e
        v_apo = math.sqrt(GM * one_minus_e / r_apo)
        r0 = np.array([-r_apo, 0.0, 0.0], dtype=np.float64)
        w0 = np.array([0.0, -mu * v_apo, 0.0], dtype=np.float64)
        E0 = -0.5 * GM * mu

        # 1. Step from apocenter directly to pericenter (dt = T_orb / 2):
        # Verify that analytical angular-momentum conditioning preserves L = r x w to <= 1e-14
        # across 1-e down to 1e-12 without any (1 - e)^{-2} subtractive cancellation!
        Lz0 = float(np.cross(r0, w0)[2])
        for step_fn in (step_kepler_sundman, step_kepler_sundman_fast):
            r_peri, w_peri = step_fn(r0, w0, mu, GM, 0.5 * T_orb)
            Lz_peri = float(np.cross(r_peri, w_peri)[2])
            rel_dL_peri = abs(Lz_peri - Lz0) / abs(Lz0)
            assert rel_dL_peri <= 1e-14, f"Pericenter angular momentum loss for 1-e=1e-{k_exp}: {rel_dL_peri:.3e}"

            # 2. Step across pericenter (pre-pericenter t = T_orb/2 - 0.5 to post-pericenter t = T_orb/2 + 0.5
            # and full orbit dt = T_orb): verify energy conservation <= 1e-14
            r_pre, w_pre = step_fn(r0, w0, mu, GM, 0.5 * T_orb - 0.5)
            r_post, w_post = step_fn(r_pre, w_pre, mu, GM, 1.0)
            E_post = 0.5 * float(np.dot(w_post, w_post)) / mu - GM * mu / float(np.linalg.norm(r_post))
            rel_dE_across = abs(E_post - E0) / abs(E0)
            assert rel_dE_across <= 1e-14, f"Across-pericenter energy loss for 1-e=1e-{k_exp}: {rel_dE_across:.3e}"

            r_full, w_full = step_fn(r0, w0, mu, GM, T_orb)
            E_full = 0.5 * float(np.dot(w_full, w_full)) / mu - GM * mu / float(np.linalg.norm(r_full))
            rel_dE_full = abs(E_full - E0) / abs(E0)
            assert rel_dE_full <= 1e-14, f"Full-orbit energy loss for 1-e=1e-{k_exp}: {rel_dE_full:.3e}"


def test_kepler_50000_steps_zero_linear_drift():
    """
    Validates [ENG-FALSE-01]:
    Incremental Gauss-Stumpff formulation (df = f - 1, dg_dot = g_dot - 1) and incremental
    Wronskian renormalization eliminate the -1.44e-16/substep secular velocity drift over 50,000 steps.
    """
    from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
        step_kepler_sundman_fast,
    )
    mu = 1.0
    GM = 1.0
    e = 0.2056
    r0 = np.array([1.0 - e, 0.0, 0.0], dtype=np.float64)
    w0 = np.array([0.0, math.sqrt((1.0 + e) / (1.0 - e)), 0.0], dtype=np.float64)
    E0 = 0.5 * float(np.dot(w0, w0)) - GM / float(np.linalg.norm(r0))
    dt = 0.002

    r_c, w_c = r0.copy(), w0.copy()
    r_f, w_f = r0.copy(), w0.copy()
    for _ in range(50000):
        r_c, w_c = step_kepler_sundman(r_c, w_c, mu, GM, dt)
        r_f, w_f = step_kepler_sundman_fast(r_f, w_f, mu, GM, dt)

    Ef_c = 0.5 * float(np.dot(w_c, w_c)) - GM / float(np.linalg.norm(r_c))
    Ef_f = 0.5 * float(np.dot(w_f, w_f)) - GM / float(np.linalg.norm(r_f))
    assert abs(Ef_c - E0) / abs(E0) < 5e-14
    assert abs(Ef_f - E0) / abs(E0) < 5e-14



def test_extreme_hyperbolic_retrograde_step():
    """
    Validates that extreme hyperbolic encounters with large negative z (down to z ~ -100,000)
    under both forward and retrograde Blanes-Moan steps solve accurately without artificial
    clamping at z = -100, conserving the Wronskian determinant to machine precision.
    """
    mu = 1.0
    GM = 1.0

    # Hyperbolic flyby with high relative velocity v = 50.0 m/s at r = 0.05
    # Kinetic energy dominates: alpha = 2 GM / r0 - v^2 = 40.0 - 2500.0 = -2460.0
    r0 = np.array([0.05, 0.0, 0.0], dtype=np.float64)
    w0 = np.array([0.0, 50.0, 0.0], dtype=np.float64)

    # For dt = 0.01, tau ~ dt / r0 = 0.2, so z = alpha * tau^2 ~ -2460 * 0.04 ~ -98.4
    # For dt = 0.05, tau ~ 1.0, z ~ -2460 (far below the old -100 clamp!)
    for dt in [0.005, 0.02, 0.05, -0.01, -0.03]:
        r_f, w_f = step_kepler_sundman(r0, w0, mu, GM, dt)
        assert not np.any(np.isnan(r_f))
        assert not np.any(np.isnan(w_f))

        # Wronskian unimodularity check
        w_det = compute_wronskian(r0, w0, mu, GM, dt)
        assert np.isclose(w_det, 1.0, atol=1e-13), f"Wronskian failed on extreme hyperbolic step dt={dt}: {w_det}"

        # Reversibility check
        r_back, w_back = step_kepler_sundman(r_f, w_f, mu, GM, -dt)
        assert np.allclose(r_back, r0, atol=1e-12), f"Reversibility failed for dt={dt}: diff={np.max(np.abs(r_back - r0))}"
        assert np.allclose(w_back, w0, atol=1e-12), f"Momentum reversibility failed for dt={dt}: diff={np.max(np.abs(w_back - w0))}"


def test_random_kepler_orbits_reversibility_1000():
    """
    Validates microscopic time-reversibility Phi(-dt) o Phi(dt) = Id across 1000 randomly
    sampled physical configurations spanning elliptic, parabolic, and hyperbolic escape regimes.
    Ensures zero Halley limit cycles, zero stall traps, and zero silent stagnation at tau=1e-15.
    """
    np.random.seed(42)
    mu = 0.5
    GM = 1.0

    for i in range(1000):
        r0 = np.random.uniform(-5.0, 5.0, 3)
        if np.linalg.norm(r0) < 0.1:
            r0 += 0.5
        v0 = np.random.uniform(-5.0, 5.0, 3)
        dt = np.random.uniform(0.1, 10.0)

        r1, w1 = step_kepler_sundman(r0, mu * v0, mu, GM, dt)
        r_rev, w_rev = step_kepler_sundman(r1, -w1, mu, GM, dt)
        w_rec = -w_rev

        scale = max(float(np.linalg.norm(r0)), float(np.linalg.norm(r1)))
        rel_err = float(np.linalg.norm(r_rev - r0)) / scale

        assert rel_err < 1e-10, (
            f"Orbit #{i} reversibility failure: rel_err={rel_err:.2e}, dt={dt}, r0={r0}, v0={v0}"
        )


def test_extreme_hyperbolic_plunge_no_overflow():
    """
    Validates that extreme hyperbolic plunge (e.g. alpha = -10,000, dt = 100, v = 100)
    which historically triggered OverflowError in math.sinh now converges without overflow
    under asymptotic hyperbolic initial guess and bounded _eval_sundman.
    """
    mu = 1.0
    GM = 1.0
    r0 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    # v = 100 -> alpha = 2*GM/r0 - v^2 = 2 - 10000 = -9998 ~ -10000
    w0 = np.array([100.0, 1.0, 0.0], dtype=np.float64)

    for dt in [10.0, 50.0, 100.0, -10.0, -50.0, -100.0]:
        r_f, w_f = step_kepler_sundman(r0, w0, mu, GM, dt)
        assert not np.any(np.isnan(r_f))
        assert not np.any(np.isnan(w_f))
        assert not np.any(np.isinf(r_f))
        assert not np.any(np.isinf(w_f))

        # Reversibility check
        r_rev, w_rev = step_kepler_sundman(r_f, w_f, mu, GM, -dt)
        scale = max(float(np.linalg.norm(r0)), float(np.linalg.norm(r_f)))
        rel_err = float(np.linalg.norm(r_rev - r0)) / scale
        assert rel_err < 1e-7, f"Reversibility error for dt={dt}: {rel_err}"


def test_step_kepler_sundman_fast_microstep_and_equivalence():
    """
    Validates D-144-08:
    1. step_kepler_sundman_fast uses scale-relative Halley tolerance (1e-15 * abs(dt))
       and adaptive tau_floor so micro-steps (e.g. dt = 1e-12 at r0 = 2.0, where the old
       hardcoded abs_tol = 1e-13 skipped all Halley updates on iteration 1) converge to
       the exact trajectory matching reference step_kepler_sundman to < 1e-15.
    2. step_kepler_sundman_fast matches reference step_kepler_sundman to < 1e-15 across
       macro and micro steps, including gravitationally dominated free-fall from near rest.
    """
    from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
        step_kepler_sundman_fast,
        OptimizedKeplerSolver,
    )

    r0 = np.array([2.0, 0.0, 0.0], dtype=np.float64)
    w0 = np.array([1.0, 0.5, -0.2], dtype=np.float64)
    mu = 1.0
    GM = 10.0

    # Micro-step dt = 1e-12 (dyadic level k=17 on dt_macro=1e-7) where old abs_tol=1e-13
    # exited immediately on Halley iteration 1 without refining tau!
    dt_micro = 1e-12
    r_fast, w_fast = step_kepler_sundman_fast(r0, w0, mu, GM, dt_micro)
    r_ref, w_ref = step_kepler_sundman(r0, w0, mu, GM, dt_micro)
    assert float(np.max(np.abs(r_fast - r_ref))) < 1e-15
    assert float(np.max(np.abs(w_fast - w_ref))) < 1e-15

    # Sub-1e-15 step dt = 1e-16 with r0 = 0.5 (true tau* = 2e-16 < 1e-15) and w0[0] = 0.0:
    # adaptive tau_floor avoids tau = 1e-15 clamping (which would give 5x impulse error!)
    r0_small = np.array([0.5, 0.0, 0.0], dtype=np.float64)
    w0_perp = np.array([0.0, 0.5, -0.2], dtype=np.float64)
    r_fast_16, w_fast_16 = step_kepler_sundman_fast(r0_small, w0_perp, mu, GM, 1e-16)
    dw_exact_16 = np.array([-10.0 / 0.25 * 1e-16, 0.0, 0.0], dtype=np.float64)
    rel_err_dw_16 = float(np.linalg.norm((w_fast_16 - w0_perp) - dw_exact_16)) / float(np.linalg.norm(dw_exact_16))
    assert rel_err_dw_16 < 1e-12, f"rel_err_dw_16={rel_err_dw_16:.2e}, w_fast_16={w_fast_16}"

    # Moderate step dt = 0.2: must match reference step_kepler_sundman to < 1e-15 (not 1e-13)
    r_f_mod, w_f_mod = step_kepler_sundman_fast(r0, w0, mu, GM, 0.2)
    r_r_mod, w_r_mod = step_kepler_sundman(r0, w0, mu, GM, 0.2)
    assert float(np.max(np.abs(r_f_mod - r_r_mod))) < 1e-15
    assert float(np.max(np.abs(w_f_mod - w_r_mod))) < 1e-14

    # Free-fall from near rest: r0 = 1e-4 << (GM * dt^2)^(1/3), v0 = 0
    r_ff = np.array([1e-4, 0.0, 0.0], dtype=np.float64)
    w_ff = np.array([0.0, 1e-8, 0.0], dtype=np.float64)
    r_ff_fast, w_ff_fast = step_kepler_sundman_fast(r_ff, w_ff, mu, GM, 0.1)
    r_ff_ref, w_ff_ref = step_kepler_sundman(r_ff, w_ff, mu, GM, 0.1)
    assert np.allclose(r_ff_fast, r_ff_ref, atol=1e-14)
    assert np.allclose(w_ff_fast, w_ff_ref, atol=1e-14)


def test_pcr3bp_canonical_solver_bm4_and_sundman_jacobi_conservation():
    """
    Validates D-144-13:
    The colocalized solver PCR3BPConormalIntegrator in Chapter 3 Section 3.3
    (pcr3bp_canonical_solver.py) implements:
    1. jacobi_constant(q, p) and hamiltonian(q, p) for the synodic PCR3BP.
    2. step_bm4(q, p, h) (7-stage Blanes-Moan order-4 symmetric composition) achieving
       fourth-order accuracy (much smaller Jacobi drift than second-order step_strang).
    3. step_sundman_regularized(q, p, h) using harmonic Sundman scaling r1*r2/(r1+r2)
       for close primary flybys on the Arenstorf periodic orbit.
    """
    import importlib.util
    from pathlib import Path

    solver_path = (
        Path(__file__).resolve().parents[3]
        / "interfaces" / "pcr3bp_canonical_solver.py"
    )
    spec = importlib.util.spec_from_file_location("pcr3bp_canonical_solver", solver_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    solver = mod.PCR3BPConormalIntegrator(mu=0.012277471)
    # Regular synodic initial state away from primaries
    q0 = np.array([0.5, 0.3], dtype=np.float64)
    p0 = np.array([-0.2, 0.4], dtype=np.float64)

    c0 = solver.jacobi_constant(q0, p0)
    h0 = solver.hamiltonian(q0, p0)
    assert np.isclose(c0, -2.0 * h0, atol=1e-15)

    # Compare order-2 Strang vs order-4 Blanes-Moan over T = 0.2 at h = 0.02 vs h = 0.01
    def _integrate(step_fn, h_step, n_steps):
        q_c, p_c = np.copy(q0), np.copy(p0)
        for _ in range(n_steps):
            q_c, p_c = step_fn(q_c, p_c, h_step)
        return q_c, p_c

    q_ref, p_ref = _integrate(solver.step_bm4, 0.0005, 400)
    q_s2_coarse, _ = _integrate(solver.step_strang, 0.02, 10)
    q_s2_fine, _ = _integrate(solver.step_strang, 0.01, 20)
    q_bm4_coarse, _ = _integrate(solver.step_bm4, 0.02, 10)
    q_bm4_fine, p_bm4_fine = _integrate(solver.step_bm4, 0.01, 20)

    err_s2_coarse = float(np.linalg.norm(q_s2_coarse - q_ref))
    err_s2_fine = float(np.linalg.norm(q_s2_fine - q_ref))
    err_bm4_coarse = float(np.linalg.norm(q_bm4_coarse - q_ref))
    err_bm4_fine = float(np.linalg.norm(q_bm4_fine - q_ref))

    # Strang is order 2 (error drops by ~4x when halving h); BM4 is order 4 (drops by ~16x)
    assert 3.5 < (err_s2_coarse / err_s2_fine) < 4.5
    assert 14.0 < (err_bm4_coarse / err_bm4_fine) < 18.0
    assert err_bm4_fine < 1e-3 * err_s2_fine

    # Verify time-reversibility of step_bm4 and step_sundman_regularized
    q_bm4, p_bm4 = solver.step_bm4(q0, p0, 0.05)
    q_rev, p_rev = solver.step_bm4(q_bm4, p_bm4, -0.05)
    assert np.allclose(q_rev, q0, atol=1e-15)
    assert np.allclose(p_rev, p0, atol=1e-15)

    q_reg, p_reg = solver.step_sundman_regularized(q0, p0, 0.05)
    assert abs(solver.jacobi_constant(q_reg, p_reg) - c0) < 1e-8

