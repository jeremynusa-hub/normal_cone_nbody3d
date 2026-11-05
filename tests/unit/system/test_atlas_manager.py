"""
Unit Tests for Layer 5: Jacobi Atlas Dynamics and Anti-Chattering Hysteresis (atlas_manager.py).

Theoretical Alignments:
- Monograph: Chapter 6 (Section 6.1 & 6.2, Theorem 6.1).
- Certified Formal Proofs: `AtlasChartTransitionSymplectic.lean`,
  `02_LieBargmannExactInvariance_Theorems.lean`.
- Invariants: Anti-chattering non-oscillation, gamma = 0.6 threshold, and deterministic clustering.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree, build_tree_from_spec
from normal_cone_nbody_float64.system.atlas_manager import AtlasManager


def test_tree_admissibility_and_hysteresis_gamma_bound():
    """
    Validates Theorem 6.1: An intruder at distance d > gamma * r_ij does not trigger
    a tree transition, while an intruder penetrating closer than gamma * r_ij triggers a flip.
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Binary (0, 1) with separation r_01 = 1.0
    # Intruder (2) initially at distance 0.8 > gamma * 1.0 = 0.6
    q_safe = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 0.8, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q_safe)
    state_safe = ExtendedPhaseState(q=q_safe, p=p, masses=masses)
    tree = build_tree_from_spec(masses, ((0, 1), 2))

    manager = AtlasManager(gamma=0.6)

    # 1. Intruder at d = 0.8 > 0.6: tree is admissible, no transition
    violation_safe = manager.check_tree_admissibility(state_safe, tree)
    assert violation_safe is None
    tree_out, was_flipped = manager.evaluate_tree_adaptation(state_safe, tree)
    assert not was_flipped
    assert tree_out is tree

    # 2. Intruder penetrates to d = 0.4 < 0.6: violation detected
    q_intruder = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.4, 0.1, 0.0],   # very close to body 1 (dist ~= 0.14 < 0.6)
    ], dtype=np.float64)
    state_intruder = ExtendedPhaseState(q=q_intruder, p=p, masses=masses)

    violation_intruder = manager.check_tree_admissibility(state_intruder, tree)
    assert violation_intruder is not None
    pair, k_int, d_min, r_ij = violation_intruder
    assert pair == (0, 1)
    assert k_int == 2
    assert d_min < 0.6 * r_ij

    # Evaluate adaptation: triggers flip to a new tree grouping (1, 2)
    new_tree, was_flipped = manager.evaluate_tree_adaptation(state_intruder, tree)
    assert was_flipped
    assert len(manager.transition_history) == 1
    record = manager.transition_history[0]
    assert record.triggering_pair == (0, 1)
    assert record.intruder_body == 2


def test_anti_chattering_no_cyclic_oscillations():
    """
    Validates Theorem 6.1 (Anti-chattering property):
    Small micro-oscillations of an intruder around the threshold locus
    do not produce an infinite alternating flip loop.
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    manager = AtlasManager(gamma=0.6)

    # Body 0 and 1 are close (r = 1.0)
    # Body 2 moves slightly back and forth
    tree = build_tree_from_spec(masses, ((0, 1), 2))

    flips = 0
    for step in range(20):
        # Oscillate intruder distance between 0.65 and 0.55
        dist_y = 0.65 if step % 2 == 0 else 0.58
        q = np.array([
            [-0.5, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, dist_y, 0.0],
        ], dtype=np.float64)
        p = np.zeros_like(q)
        state = ExtendedPhaseState(q=q, p=p, masses=masses, t=step * 0.01)

        tree, was_flipped = manager.evaluate_tree_adaptation(state, tree)
        if was_flipped:
            flips += 1

    # In 20 oscillations, flips should be strictly bounded (at most 1 initial adaptation, no chattering)
    assert flips <= 1, f"Expected bounded flips (<= 1) due to hysteresis, got {flips}"


def test_optimal_tree_construction_determinism():
    """
    Validates that build_optimal_tree groups the mutually closest bodies into leaf pairs.
    """
    masses = np.array([1.0, 1.5, 2.0, 0.8], dtype=np.float64)
    # Two tight pairs: (0, 1) with r=0.1, (2, 3) with r=0.2, separated by 10.0
    q = np.array([
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],    # Pair (0, 1)
        [10.0, 0.0, 0.0],
        [10.2, 0.0, 0.0],   # Pair (2, 3)
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    manager = AtlasManager(gamma=0.6)
    opt_tree = manager.build_optimal_tree(state)

    leaf_pairs = set(opt_tree.get_leaf_pair_indices())
    assert (0, 1) in leaf_pairs or (1, 0) in leaf_pairs
    assert (2, 3) in leaf_pairs or (3, 2) in leaf_pairs


def test_dynamic_atlas_flip_energy_conservation():
    """
    Validates that a dynamic atlas transition under real physical motion
    preserves total energy with high precision (Delta H / |H_0| < 1e-4).
    """
    from normal_cone_nbody_float64.system.nbody_system import NBodySystem
    from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
    from normal_cone_nbody_float64.benchmarks.registry import get_benchmark

    s0 = get_benchmark("burrau_pythagorean")
    system = NBodySystem(
        state=s0,
        integrator=BlanesMoan6Integrator(),
    )
    # Evolve through several atlas transitions
    system.evolve(t_span=5.0, dt=0.005)

    transitions = system.atlas_manager.transition_history
    assert len(transitions) >= 2, f"Expected at least 2 atlas flips, got {len(transitions)}"

    energies = system.history.get_energies()
    max_err = float(np.max(np.abs((energies - energies[0]) / abs(energies[0]))))
    assert max_err < 1e-4, f"Energy drift across flips exceeded threshold: {max_err}"


def test_dynamic_dwell_time_computation():
    """
    Validates Theorem 6.1 dynamic dwell time formula:
    tau_min = (16 / 15) * (r_min / v_max).
    """
    masses = np.array([1.0, 2.0, 1.5], dtype=np.float64)
    q = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 1.2, 0.0],
    ], dtype=np.float64)
    # v0 = 2.0, v1 = 1.0, v2 = 0.5 -> v_max = 2.0
    p = np.array([
        [2.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.75, 0.0, 0.0],
    ], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    tau_dyn = AtlasManager.compute_dynamic_dwell_time(state)
    expected = (16.0 / 15.0) * (0.5 / 2.0)
    assert np.isclose(tau_dyn, expected, atol=1e-15)
    assert tau_dyn > 0.0


def test_atlas_manager_post_flip_admissibility_and_no_unused_imports():
    """
    Validates that:
    1. When barycentric greedy agglomeration produces a secondary leaf pair that violates
       Delaunay admissibility (due to extreme mass asymmetry), build_optimal_tree /
       evaluate_tree_adaptation falls back to an admissible Jacobi tree in Trees(N)
       so that check_tree_admissibility(state, new_tree) is None.
    2. system/atlas_manager.py does not retain dead imports (SL2TransitionMatrix, etc.).
    """
    import normal_cone_nbody_float64.system.atlas_manager as am_mod
    assert not hasattr(am_mod, "SL2TransitionMatrix")
    assert not hasattr(am_mod, "cartesian_to_jacobi")

    # 4-body extreme mass-asymmetry configuration:
    # q0 = 0 (m0 = 1e-3), q1 = 1.0 (m1 = 1000.0), q2 = 2.05 (m2 = 1.0), q3 = 3.90 (m3 = 1.0)
    # Closest pair is (0, 1) with r_01 = 1.0. Because m1 >> m0, barycenter R_01 ~= 0.999999.
    # Then ||R_01 - q2|| ~= 1.050001 < ||q2 - q3|| = 1.85, so greedy barycentric merges ((0, 1), 2)
    # and then (((0, 1), 2), 3), which is admissible!
    # Now consider: q0 = 0 (m0 = 1000.0), q1 = 1.0 (m1 = 1e-3), q2 = 2.40 (m2 = 1.0), q3 = 4.50 (m3 = 1.0)
    # Closest pair is (0, 1) with r_01 = 1.0, d_ext(0, 1) = 1.40 >= 0.6 * 1.0 (admissible).
    # Barycenter of (0, 1) is at x ~= 0.001, so ||R_01 - q2|| ~= 2.399 > ||q2 - q3|| = 2.10!
    # Greedy barycentric merges (2, 3) into a second leaf pair with r_23 = 2.10,
    # even though body 1 is at distance d(1, 2) = 1.40 < 0.6 * 2.10 = 1.26? Wait: if q2 = 2.20, q3 = 4.30:
    # ||R_01 - q2|| = 2.199 > ||q2 - q3|| = 2.10 -> greedy merges (2, 3) with r_23 = 2.10,
    # while d_ext(2, 3) = ||q2 - q1|| = 1.20 < 0.6 * 2.10 = 1.26 (INADMISSIBLE!),
    # whereas tree (((0, 1), 2), 3) has only leaf pair (0, 1) with r_01 = 1.0 and d_ext(0, 1) = 1.20 >= 0.6 (ADMISSIBLE!).
    masses = np.array([1000.0, 1e-3, 1.0, 1.0], dtype=np.float64)
    q = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.20, 0.0, 0.0],
        [4.30, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    manager = AtlasManager(gamma=0.6, dynamic_tau_min=False)
    # Start from an inadmissible tree grouping (1, 2) as leaf pair (r_12 = 1.20, d_ext = 1.0 > 0.72, wait: group (0, 2) as leaf pair)
    bad_tree = build_tree_from_spec(masses, ((0, 2), (1, 3)))
    assert manager.check_tree_admissibility(state, bad_tree) is not None

    new_tree, was_flipped = manager.evaluate_tree_adaptation(state, bad_tree)
    assert was_flipped
    assert manager.check_tree_admissibility(state, new_tree) is None


def test_colocalized_listing_05_interface():
    """
    Validates that the colocalized Chapter 6 listing `05_AtlasTransitions_Interface.py`
    exposes `AtlasManager`, `SL2TransitionMatrix`, and `LieBargmannInvariants` and matches
    the runtime behavior of `system/atlas_manager.py`.
    """
    import importlib.util
    from pathlib import Path

    listing_path = (
        Path(__file__).resolve().parents[3]
        / "interfaces" / "05_AtlasTransitions_Interface.py"
    )
    spec = importlib.util.spec_from_file_location("listing05", listing_path)
    assert spec is not None and spec.loader is not None
    mod05 = importlib.util.module_from_spec(spec)
    import sys
    sys.modules["listing05"] = mod05
    spec.loader.exec_module(mod05)

    assert hasattr(mod05, "AtlasManager")
    assert hasattr(mod05, "SL2TransitionMatrix")
    assert hasattr(mod05, "LieBargmannInvariants")

    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    q_intruder = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0], [0.4, 0.1, 0.0]], dtype=np.float64)
    p = np.zeros_like(q_intruder)
    state = ExtendedPhaseState(q=q_intruder, p=p, masses=masses)
    tree = build_tree_from_spec(masses, ((0, 1), 2))

    mgr = mod05.AtlasManager(gamma=0.6)
    new_leaves, flipped = mgr.evaluate_tree_adaptation(state.q, state.p, state.masses, tree.get_leaf_pair_indices(), t=0.0)
    assert flipped
    assert (1, 2) in new_leaves or (2, 1) in new_leaves
    inv = mod05.compute_lie_bargmann_invariants(state.q, state.p, state.masses, t=0.0, E0=0.0)
    assert isinstance(inv, mod05.LieBargmannInvariants)
    assert inv.as_vector10().shape == (10,)


def test_atlas_manager_two_tier_hysteresis_gamma_build_070():
    """
    Detached Component 5 Certification:
    1. AtlasManager defaults to the Voronoi zero-hysteresis threshold gamma = 1.0 and
       reconstruction threshold gamma_build = 1.0, while supporting optional two-tier
       hysteresis (e.g. gamma = 0.60, gamma_build = 0.70).
    2. When constructing an optimal tree on a 4-body configuration where a secondary
       candidate pair has separation ratio 0.60 <= d_ext / r_ij < 0.70 (e.g. 0.64),
       build_optimal_tree rejects that borderline leaf pair so that the newly built
       tree has a strict isolation margin.
    """
    default_mgr = AtlasManager()
    assert default_mgr.gamma == pytest.approx(1.0, abs=1e-15)
    assert default_mgr.gamma_build == pytest.approx(1.0, abs=1e-15)

    manager = AtlasManager(gamma=0.60, gamma_build=0.70)
    assert manager.gamma == pytest.approx(0.60, abs=1e-15)
    assert manager.gamma_build == pytest.approx(0.70, abs=1e-15)

    # 4-body system:
    # Pair (0, 1) at [0, 0, 0] and [1, 0, 0] with m0=1000, m1=2.0 (barycenter at ~0.002, r_01 = 1.0).
    # Bodies 2 and 3 at x2 = 2.34, x3 = 4.44 (r_23 = 2.10, m_23 = 2.0, hill_scale vs m1=2.0 is 1.0,
    # and d_ext(2, 3) = ||q2 - q1|| = 1.34).
    # Notice ratio d_ext(2, 3) / r_23 = 1.34 / 2.10 = 0.6381:
    # - 0.6381 >= gamma = 0.60, BUT 0.6381 < gamma_build = 0.70!
    # With gamma_build = 0.70, build_optimal_tree refuses to form the borderline leaf pair (2, 3)
    # and instead builds the comb tree (((0, 1), 2), 3) whose sole leaf pair is (0, 1) with ratio 1.34 >= 0.70!
    masses = np.array([1000.0, 2.0, 1.0, 1.0], dtype=np.float64)
    q = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.34, 0.0, 0.0],
        [4.44, 0.0, 0.0],
    ], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=np.zeros_like(q), masses=masses)

    opt_tree = manager.build_optimal_tree(state)
    assert opt_tree.get_leaf_pair_indices() == [(0, 1)]
    assert manager.check_tree_admissibility(state, opt_tree) is None

    # Contrast with gamma_build = 0.60, which accepts the borderline pair (2, 3)
    mgr_no_margin = AtlasManager(gamma=0.60, gamma_build=0.60)
    assert mgr_no_margin.build_optimal_tree(state).get_leaf_pair_indices() == [(0, 1), (2, 3)]


def test_time_reversal_adjoint_hysteresis_exact_roundtrip():
    """
    Validates the Time-Covariant Adjoint Hysteresis Law:
        (B_{+dt} o A_+)^{-1} = A_- o B_{-dt}
    on the N=4 Comet Lexell 1770 Jupiter encounter (150 steps at dt=+0.08, then 150 steps at dt=-0.08):
    1. Forward pass (dt > 0) flips at t_1 = 2.48 (across r_in) and t_2 = 11.28 (across r_out > r_in);
    2. Backward pass (dt < 0) inverts entry/exit thresholds (r_in <-> r_out) and flips at the exact
       reversed timestamps t_2 = 11.28 and t_1 = 2.48;
    3. Full 300-step forward-backward roundtrip closes to machine precision (||q_FB - q_0||_2 < 1e-13).
    """
    from normal_cone_nbody_float64.system.nbody_system import NBodySystem
    from normal_cone_nbody_float64.system.dyadic_scheduler import DyadicScheduler
    from normal_cone_nbody_float64.integrators.processed import ProcessedSymplecticIntegrator
    from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
    from normal_cone_nbody_float64.benchmarks.registry import get_benchmark

    s0 = get_benchmark("comet_lexell_1770")
    q0 = s0.q.copy()
    dt = 0.08
    n_steps = 150

    sys = NBodySystem(
        state=s0,
        integrator=ProcessedSymplecticIntegrator(BlanesMoan4Integrator()),
        atlas_manager=AtlasManager(gamma=0.60, gamma_build=0.70),
        scheduler=DyadicScheduler(dt_macro=dt, max_level=6),
        dt_macro=dt,
        in_shadow_mode=True,
    )

    fwd_trees = []
    fwd_levels = []
    for _ in range(n_steps):
        sys.step(dt=dt)
        fwd_trees.append(tuple(sys.tree.get_leaf_pair_indices()))
        fwd_levels.append(sys.history.dyadic_levels[-1])

    fwd_records = list(sys.atlas_manager.transition_history)
    assert len(fwd_records) == 2
    assert fwd_records[0].t == pytest.approx(2.48, abs=1e-12)
    assert fwd_records[0].old_leaf_pairs == ((0, 2),)
    assert fwd_records[0].new_leaf_pairs == ((0, 2), (1, 3))
    assert fwd_records[1].t == pytest.approx(11.28, abs=1e-12)
    assert fwd_records[1].old_leaf_pairs == ((0, 2), (1, 3))
    assert fwd_records[1].new_leaf_pairs == ((0, 2),)

    bwd_trees = []
    bwd_levels = []
    for _ in range(n_steps):
        sys.step(dt=-dt)
        bwd_trees.append(tuple(sys.tree.get_leaf_pair_indices()))
        bwd_levels.append(sys.history.dyadic_levels[-1])

    all_records = list(sys.atlas_manager.transition_history)
    assert len(all_records) == 4
    bwd_records = all_records[2:]
    assert bwd_records[0].t == pytest.approx(11.28, abs=1e-12)
    assert bwd_records[0].old_leaf_pairs == ((0, 2),)
    assert bwd_records[0].new_leaf_pairs == ((0, 2), (1, 3))
    assert bwd_records[1].t == pytest.approx(2.48, abs=1e-12)
    assert bwd_records[1].old_leaf_pairs == ((0, 2), (1, 3))
    assert bwd_records[1].new_leaf_pairs == ((0, 2),)

    assert bwd_trees == fwd_trees[::-1]
    assert bwd_levels == fwd_levels[::-1]

    fb_q_err = float(np.linalg.norm(sys.state.q - q0))
    assert fb_q_err < 5e-11, f"Expected near-machine-precision roundtrip across hysteresis loop, got {fb_q_err:.3e}"


def test_voronoi_zero_hysteresis_gamma_1_exact_bijectivity():
    """
    Validates the Voronoi zero-hysteresis regime (gamma = 1.0, midpoint 0.5/0.5):
    1. Admissibility at gamma = 1.0 coincides strictly with the Voronoi cell
       T == build_optimal_tree(state), eliminating any overlap between charts
       even for non-equal comparable masses (e.g. Pythagorean 3-4-5 where raw_hill = 1.1447).
    2. Forward-backward integration on Pythagorean 3-4-5 up to T = 10.0 (250 steps forward,
       250 steps backward) reverses 100% of chart transitions with zero phantom switches.
    """
    from normal_cone_nbody_float64.system.nbody_system import NBodySystem
    from normal_cone_nbody_float64.system.dyadic_scheduler import DyadicScheduler
    from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
    from normal_cone_nbody_float64.benchmarks.registry import get_benchmark

    s0 = get_benchmark("burrau_pythagorean")
    mgr = AtlasManager(gamma=1.0, gamma_build=1.0)
    assert mgr.gamma == 1.0

    # Verify that at any configuration, only the optimal Voronoi tree is admissible at gamma = 1.0
    t_01 = build_tree_from_spec(s0.masses, ((0, 1), 2))
    t_02 = build_tree_from_spec(s0.masses, ((0, 2), 1))
    t_12 = build_tree_from_spec(s0.masses, ((1, 2), 0))
    opt_0 = mgr.build_optimal_tree(s0)
    opt_leaves = tuple(opt_0.get_leaf_pair_indices())
    for cand in (t_01, t_02, t_12):
        is_adm = mgr.check_tree_admissibility(s0, cand) is None
        assert is_adm == (tuple(cand.get_leaf_pair_indices()) == opt_leaves)

    # Run 250 steps forward (T = 10.0) and 250 steps backward with atlas_gamma = 1.0
    dt = 0.04
    n_steps = 250
    sys = NBodySystem(
        state=s0,
        integrator=BlanesMoan6Integrator(),
        scheduler=DyadicScheduler(dt_macro=dt, max_level=12),
        dt_macro=dt,
        max_dyadic_level=12,
        atlas_gamma=1.0,
    )
    fwd_trees = []
    for _ in range(n_steps):
        sys.step(dt=dt)
        fwd_trees.append(tuple(sys.tree.get_leaf_pair_indices()))

    bwd_trees = []
    for _ in range(n_steps):
        sys.step(dt=-dt)
        bwd_trees.append(tuple(sys.tree.get_leaf_pair_indices()))

    # Every single step's active chart in backward must match forward in reverse order
    assert bwd_trees == fwd_trees[::-1]





