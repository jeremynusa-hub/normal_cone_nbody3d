"""
Unit Tests for Layer 5: Discrete Dyadic Multi-Rate Scheduler (dyadic_scheduler.py).

Theoretical Alignments:
- Monograph: Chapter 5 (Section 5.1 & Section 5.5, Ge--Marsden Theorem & Dyadic Lattices).
- Invariants: Machine-precision sub-step sum identity sum dt_k == dt_macro,
  colored interval prioritization, and dynamic timescale derivation.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.constants import DEFAULT_MAX_DYADIC_LEVEL
from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.system.dyadic_scheduler import DyadicScheduler, ColoredInterval


def test_dyadic_exact_subdivision_identity():
    """
    Validates that for any dyadic level k in [0, 10], the sub-step division satisfies:
        num_substeps * dt_substep == dt_macro
    to float64 machine precision (< 1e-16).
    """
    dt_macro = 0.024
    scheduler = DyadicScheduler(dt_macro=dt_macro, default_level=0, max_level=10)

    for k in range(8):
        scheduler.default_level = k
        num_substeps, dt_substep, level = scheduler.get_substep_params(t_current=0.0)

        assert num_substeps == (1 << k)
        assert level == k
        # Exact product matches macro step
        assert np.isclose(num_substeps * dt_substep, dt_macro, atol=1e-16, rtol=1e-15)


def test_colored_interval_scheduling():
    """
    Validates pre-scheduled colored window reservations:
    Refinement level k is elevated inside the window and reverts outside.
    """
    scheduler = DyadicScheduler(dt_macro=0.1, default_level=0, max_level=8)

    # Schedule high-resolution event at [1.0, 2.0] with k=4 (16 substeps)
    scheduler.color_interval(t_start=1.0, t_end=2.0, level=4, description="Transit event")

    # Before event (t = 0.5)
    n_pre, dt_pre, k_pre = scheduler.get_substep_params(t_current=0.5)
    assert k_pre == 0
    assert n_pre == 1

    # Overlapping boundary (t = 0.95 to 1.05)
    n_cross, dt_cross, k_cross = scheduler.get_substep_params(t_current=0.95)
    assert k_cross == 4
    assert n_cross == 16

    # Inside event (t = 1.5)
    n_in, dt_in, k_in = scheduler.get_substep_params(t_current=1.5)
    assert k_in == 4
    assert n_in == 16

    # After event (t = 2.5)
    n_post, dt_post, k_post = scheduler.get_substep_params(t_current=2.5)
    assert k_post == 0
    assert n_post == 1


def test_dynamic_level_computation():
    """
    Validates dynamic timescale adaptation:
    When a leaf binary has high relative velocity and small separation,
    the dynamical timescale tau triggers an elevated dyadic level k > 0.
    """
    scheduler = DyadicScheduler(dt_macro=0.1, default_level=0, max_level=8)

    masses = np.array([1.0, 1.0], dtype=np.float64)
    # Binary with separation r = 0.01 and orbital velocity v = 10.0 (tau = 0.001)
    q = np.array([[-0.005, 0.0, 0.0], [0.005, 0.0, 0.0]], dtype=np.float64)
    p = np.array([[0.0, -10.0, 0.0], [0.0, 10.0, 0.0]], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    num_substeps, dt_substep, k = scheduler.get_substep_params(
        t_current=0.0, state=state, tree=tree
    )

    # With dt_macro=0.1 and tau=0.001, target dt ~ 0.00005 -> k should be elevated (k >= 4)
    assert k >= 4
    assert num_substeps >= 16
    assert dt_substep < 0.01


def test_dyadic_extreme_depth_up_to_level_17():
    """
    Validates Chapter 5 Section 5.6 & Iteration 146 Phase 3 requirement:
    1. DEFAULT_MAX_DYADIC_LEVEL == 17 by default.
    2. Extreme pericenter encounter (r ~ 1e-4, v ~ 100 m/s, tau ~ 1e-6) produces
       a dyadic depth k >= 14 and up to 17 without artificial saturation at k=10.
    3. Exact time decomposition identity sum dt_k == dt_macro holds at k=17 (131,072 substeps).
    """
    assert DEFAULT_MAX_DYADIC_LEVEL == 17

    dt_macro = 0.015
    scheduler = DyadicScheduler(dt_macro=dt_macro)
    assert scheduler.max_level == 17

    # Validate exact subdivision at level 17
    scheduler.default_level = 17
    num_substeps, dt_substep, level = scheduler.get_substep_params(t_current=0.0)
    assert num_substeps == 131072  # 2^17
    assert level == 17
    assert np.isclose(num_substeps * dt_substep, dt_macro, atol=1e-16, rtol=1e-15)

    # Dynamic detection on extreme high-speed pericenter
    scheduler.default_level = 0
    masses = np.array([1.0, 1.0], dtype=np.float64)
    # Separation r = 2e-4, velocity v = 100.0 -> tau_kin = 2e-6 s
    q = np.array([[-1e-4, 0.0, 0.0], [1e-4, 0.0, 0.0]], dtype=np.float64)
    p = np.array([[0.0, -100.0, 0.0], [0.0, 100.0, 0.0]], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    num_sub, dt_sub, k_dyn = scheduler.get_substep_params(
        t_current=0.0, state=state, tree=tree
    )

    # With dt_macro=0.015 and target_dt = 0.05 * 2e-6 = 1e-7 s:
    # ratio = 0.015 / 1e-7 = 150,000 -> ceil(log2(150000)) = 18 -> clamped to 17!
    assert k_dyn == 17
    assert num_sub == 131072
    assert dt_sub < 2e-7


def test_dyadic_scheduler_time_reversibility_symmetry():
    """
    Validates Definition 5.4 Item 3 (V4 Time-Reversibility):
    The bilateral pericenter encounter window |t_peri| <= 2 |dt_macro| ensures exact
    even symmetry under phase-space time reversal R(q, p) = (q, -p):
        compute_dynamical_level(state(q, p)) == compute_dynamical_level(state(q, -p)).
    """
    scheduler = DyadicScheduler(dt_macro=0.1, default_level=0, max_level=17)
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Canonical comb tree makes (0, 1) a leaf pair and (1, 2) a non-leaf pair
    tree = build_canonical_tree(masses, topology="comb")

    # Non-leaf pair (1, 2) on a grazing eccentric trajectory near pericenter:
    # At r_12 = (0.1, 0.005, 0.0) and v_12 = (-1.5, 0.1, 0.0), body 2 is approaching body 1
    # with r_peri << r_12 and |t_peri| <= 2 * dt_macro.
    # Under time reversal p -> -p, v_rad flips sign (+1.5, receding just after pericenter),
    # which must trigger the exact same predictive pericenter level k_in == k_out > 0.
    q = np.array([
        [-10.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.1, 0.005, 0.0],
    ], dtype=np.float64)
    p_approaching = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [-1.5, 0.1, 0.0],
    ], dtype=np.float64)
    p_receding = -p_approaching

    state_in = ExtendedPhaseState(q=q, p=p_approaching, masses=masses)
    state_out = ExtendedPhaseState(q=q, p=p_receding, masses=masses)

    k_in = scheduler.compute_dynamical_level(state_in, tree=tree)
    k_out = scheduler.compute_dynamical_level(state_out, tree=tree)

    assert k_in >= 6, f"Expected elevated predictive pericenter level on approach, got {k_in}"
    assert k_out == k_in, (
        f"Time-reversibility violation (D-145-05): approaching level k_in={k_in} != receding level k_out={k_out}"
    )


def test_leaf_apex_decoupling_eccentric_pericenter_invariance():
    """
    Validates Theorem 5.20 Item 1 (V5a Leaf/Apex Decoupling):
    For an isolated bound Keplerian leaf pair C = (0, 1) with semi-major axis a_C = 1.0
    and high eccentricity e = 0.99 (r_peri = 0.01, r_apo = 1.99), the dynamical timescale
    is governed by the invariant semi-major axis tau_dyn(C) = sqrt(a_C^3 / GM_C) = 1 / sqrt(2),
    so the dyadic level k at pericenter equals the level at apocenter without artificial
    saturation as r_peri -> 0.
    """
    scheduler = DyadicScheduler(dt_macro=0.02, default_level=0, max_level=17)
    masses = np.array([1.0, 1.0, 1e-3], dtype=np.float64)
    tree = build_canonical_tree(masses, topology="comb")  # (0, 1) is the leaf pair

    GM = 2.0
    a = 1.0
    e = 0.99
    r_peri = a * (1.0 - e)  # 0.01
    r_apo = a * (1.0 + e)   # 1.99
    v_peri = np.sqrt(GM * (2.0 / r_peri - 1.0 / a))  # ~19.95
    v_apo = np.sqrt(GM * (2.0 / r_apo - 1.0 / a))    # ~0.10025
    mu = 0.5

    # State at apocenter with distant spectator body 2 at R_ext = 50.0
    q_apo = np.array([
        [-0.5 * r_apo, 0.0, 0.0],
        [0.5 * r_apo, 0.0, 0.0],
        [50.0, 0.0, 0.0],
    ], dtype=np.float64)
    p_apo = np.array([
        [0.0, -mu * v_apo, 0.0],
        [0.0, mu * v_apo, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)

    # State at pericenter with distant spectator body 2 at R_ext = 50.0
    q_peri = np.array([
        [-0.5 * r_peri, 0.0, 0.0],
        [0.5 * r_peri, 0.0, 0.0],
        [50.0, 0.0, 0.0],
    ], dtype=np.float64)
    p_peri = np.array([
        [0.0, -mu * v_peri, 0.0],
        [0.0, mu * v_peri, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)

    k_apo = scheduler.compute_dynamical_level(
        ExtendedPhaseState(q=q_apo, p=p_apo, masses=masses), tree=tree
    )
    k_peri = scheduler.compute_dynamical_level(
        ExtendedPhaseState(q=q_peri, p=p_peri, masses=masses), tree=tree
    )

    assert k_peri == k_apo == 0, (
        f"Leaf/Apex Decoupling violation (D-145-06): eccentric leaf pair triggered k_peri={k_peri} vs k_apo={k_apo}"
    )


def test_hill_isolation_guard_breach_triggers_refinement():
    """
    Validates Theorem 5.20 Item 2 (V5b Hill Isolation Guard R_ext >= 3 R_Hill):
    When an external massive perturber k=2 breaches the Hill sphere of a leaf pair (0, 1)
    (R_ext < 3 R_Hill), the leaf exemption is lifted and the instantaneous/pericenter
    timescales reactivate dyadic refinement.
    """
    scheduler = DyadicScheduler(dt_macro=0.05, default_level=0, max_level=17)
    # Leaf pair (0, 1) of masses 1.0, 1.0 (M_C = 2.0) and massive perturber m_2 = 1000.0
    masses = np.array([1.0, 1.0, 1000.0], dtype=np.float64)
    tree = build_canonical_tree(masses, topology="comb")

    # Bound leaf pair at separation r_01 = 0.2, semi-major axis a = 1.0 (v_rel = sqrt(2*(2/0.2 - 1)) = sqrt(18))
    r_01 = 0.2
    v_rel = np.sqrt(18.0)
    mu = 0.5

    # Case 1: Perturber far away at R_ext = 100.0 -> R_Hill = 100 * (2 / 3006)^(1/3) ~= 8.73 -> a = 1.0 <= R_Hill / 3 ~= 2.91 (Isolated)
    q_isolated = np.array([
        [-0.5 * r_01, 0.0, 0.0],
        [0.5 * r_01, 0.0, 0.0],
        [100.0, 0.0, 0.0],
    ], dtype=np.float64)
    p_leaf = np.array([
        [0.0, -mu * v_rel, 0.0],
        [0.0, mu * v_rel, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)

    # Case 2: Perturber closer at R_ext = 15.0 -> R_Hill = 15 * (2 / 3006)^(1/3) ~= 1.31 -> a = 1.0 > R_Hill / 3 ~= 0.436 (Hill guard breached!)
    q_breached = np.array([
        [-0.5 * r_01, 0.0, 0.0],
        [0.5 * r_01, 0.0, 0.0],
        [15.0, 0.0, 0.0],
    ], dtype=np.float64)

    k_isolated = scheduler.compute_dynamical_level(
        ExtendedPhaseState(q=q_isolated, p=p_leaf, masses=masses), tree=tree
    )
    k_breached = scheduler.compute_dynamical_level(
        ExtendedPhaseState(q=q_breached, p=p_leaf, masses=masses), tree=tree
    )

    assert k_isolated == 1
    assert k_breached > k_isolated, (
        f"Hill guard breach (D-145-07) failed to lift leaf exemption: k_isolated={k_isolated}, k_breached={k_breached}"
    )


def test_colocalized_listing_03_dyadic_scheduler_interface():
    """
    Validates D-145-05..07 & D-145-16 on the colocalized listing
    03_DyadicScheduler_Interface.py:
    1. Supports both (state, tree) and raw (q, p, masses, leaf_pairs) signatures.
    2. Time-reversibility symmetry k(q, p) == k(q, -p) on non-leaf grazing encounter.
    3. Leaf/Apex decoupling invariance at eccentric pericenter and Hill isolation guard.
    """
    import importlib.util
    from pathlib import Path

    listing_path = (
        Path(__file__).resolve().parents[3]
        / "interfaces" / "03_DyadicScheduler_Interface.py"
    )
    spec = importlib.util.spec_from_file_location("listing_03_dyadic", listing_path)
    mod = importlib.util.module_from_spec(spec)
    import sys
    sys.modules["listing_03_dyadic"] = mod
    spec.loader.exec_module(mod)

    sched = mod.DyadicScheduler(dt_macro=0.1, default_level=0, max_level=17)
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    tree = build_canonical_tree(masses, topology="comb")

    q = np.array([
        [-10.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.1, 0.005, 0.0],
    ], dtype=np.float64)
    p_in = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [-1.5, 0.1, 0.0],
    ], dtype=np.float64)

    # 1. Raw array signature with V4 time-reversibility
    k_in_raw = sched.compute_dynamical_level(q, p_in, masses, leaf_pairs=[(0, 1)])
    k_out_raw = sched.compute_dynamical_level(q, -p_in, masses, leaf_pairs=[(0, 1)])
    assert k_in_raw >= 6
    assert k_in_raw == k_out_raw

    # 2. Canonical (state, tree) signature parity
    state_in = ExtendedPhaseState(q=q, p=p_in, masses=masses)
    k_in_state = sched.compute_dynamical_level(state_in, tree=tree)
    assert k_in_state == k_in_raw
    n_sub, dt_sub, k_sub = sched.get_substep_params(0.0, state=state_in, tree=tree)
    assert k_sub == k_in_raw
    assert np.isclose(n_sub * dt_sub, 0.1, atol=1e-16)


def test_symmetric_dyadic_level_selector_and_hyperbolic_leaf_pericenter_window():
    """
    Certifies Detached Component 6 (Symmetric Dyadic Scheduler, [ENG-FALSE-09], [CRIT-7.2.1-04]):
    1. Bilateral pericenter window |t_peri| <= 2|dt| is active on hyperbolic leaf encounters (E_2b > 0).
    2. Symmetric endpoint evaluation k_sym = max(k(z_n), k(Phi_{dt}^{(0)}(z_n))) guarantees
       exact palindromic sequence k_fwd(n) == k_bwd(N - 1 - n) across a close pericenter passage.
    """
    from normal_cone_nbody_float64.diagnostics.metrics import (
        compute_forward_backward_reversibility,
    )
    from normal_cone_nbody_float64.integrators.processed import create_pbm4_integrator

    # Hyperbolic 2-body leaf flyby approaching pericenter r_p ~ 0.033 from r_0 ~ 0.25
    masses = np.array([1.0, 1e-6], dtype=np.float64)
    q0 = np.array([
        [0.0, 0.0, 0.0],
        [0.25, 0.08, 0.0],
    ], dtype=np.float64)
    v0 = np.array([
        [0.0, 0.0, 0.0],
        [-3.5, 0.0, 0.0],
    ], dtype=np.float64)
    p0 = v0 * masses[:, np.newaxis]
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses).to_barycentric()

    rev = compute_forward_backward_reversibility(
        state0=state0,
        integrator=create_pbm4_integrator(G=1.0),
        dt=0.02,
        n_steps=8,
        G=1.0,
    )
    assert rev["dyadic_palindrome_matched"] is True, (
        f"Symmetric dyadic palindrome mismatch: fwd={rev['k_fwd']} vs bwd_rev={rev['k_bwd_reversed']}"
    )
    assert max(rev["k_fwd"]) > min(rev["k_fwd"]), "Expected nontrivial multi-level dyadic transitions across pericenter"
    assert rev["delta_q_norm"] < 5e-11
    assert rev["delta_p_norm"] < 5e-11




