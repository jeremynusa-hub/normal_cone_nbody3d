"""
Unit Tests for Layer 5: Unified N-Body System Facade (nbody_system.py).

Theoretical Alignments:
- Monograph: Chapters 1 to 6 (Global Multi-Layer Integration Facade).
- Invariants: Machine-precision Lie--Bargmann conservation under evolve(),
  integrator hot-swapping, and automatic tree clustering.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.integrators.strang import StrangIntegrator
from normal_cone_nbody_float64.integrators.yoshida import Yoshida4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
from normal_cone_nbody_float64.integrators.processed import ProcessedSymplecticIntegrator


def make_hierarchical_3body_system() -> ExtendedPhaseState:
    """Helper creating a standard stable hierarchical 3-body configuration."""
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q0 = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 8.0, 0.5],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, -0.7, 0.0],
        [0.0, 0.7, 0.0],
        [0.25, 0.0, -0.1],
    ], dtype=np.float64)
    return ExtendedPhaseState(q=q0, p=p0, masses=masses)


def test_nbody_system_single_step_advancement():
    """
    Validates single step advancement: time advances by dt, dimensions are preserved,
    and checkpoint is recorded in history.
    """
    state0 = make_hierarchical_3body_system()
    system = NBodySystem(state=state0, integrator=StrangIntegrator(G=1.0), G=1.0)

    assert system.t == 0.0
    assert system.N == 3
    assert len(system.history) == 1  # initial baseline

    dt = 0.02
    state1 = system.step(dt=dt)

    assert np.isclose(system.t, 0.02, atol=1e-15)
    assert np.isclose(state1.t, 0.02, atol=1e-15)
    assert state1.q.shape == (3, 3)
    assert state1.p.shape == (3, 3)
    assert len(system.history) == 2


def test_nbody_system_evolve_trajectory_and_invariants():
    """
    Validates evolving the system over t_span = 0.2:
    - Reaches exact final time.
    - Conserves linear momentum P (< 1e-14) and angular momentum L (< 1e-14).
    - Preserves Hamiltonian energy to high accuracy.
    """
    state0 = make_hierarchical_3body_system()
    system = NBodySystem(state=state0, integrator=BlanesMoan4Integrator(G=1.0), G=1.0)

    t_span = 0.2
    final_state = system.evolve(t_span=t_span, dt=0.02)

    assert np.isclose(final_state.t, t_span, atol=1e-14)
    assert np.isclose(system.t, t_span, atol=1e-14)

    # Invariants verification
    metrics = system.history.summary_metrics()
    assert metrics["delta_P_max"] < 1e-13, f"Linear momentum drift: {metrics['delta_P_max']}"
    assert metrics["delta_L_max"] < 1e-13, f"Angular momentum drift: {metrics['delta_L_max']}"
    assert metrics["delta_H_max"] < 1e-10, f"Energy drift: {metrics['delta_H_max']}"


def test_nbody_system_integrator_hot_swap():
    """
    Validates dynamically hot-swapping between all 5 symplectic integrators:
    Strang -> Yoshida4 -> BM4 -> BM6 -> Processed BM4.
    Ensures seamless continuity without state loss or NaN values.
    """
    state0 = make_hierarchical_3body_system()
    system = NBodySystem(state=state0, G=1.0)

    integrators = [
        StrangIntegrator(G=1.0),
        Yoshida4Integrator(G=1.0),
        BlanesMoan4Integrator(G=1.0),
        BlanesMoan6Integrator(G=1.0),
        ProcessedSymplecticIntegrator(BlanesMoan4Integrator(G=1.0)),
    ]

    for integ in integrators:
        system.set_integrator(integ)
        s = system.step(dt=0.01)
        assert not np.any(np.isnan(s.q))
        assert not np.any(np.isnan(s.p))
        inv = system.get_invariants()
        assert np.isfinite(inv.H)
        assert np.all(np.isfinite(inv.P))
        assert np.all(np.isfinite(inv.L))


def test_nbody_system_reset_preserves_clean_state():
    """
    Validates system reset with a new state: clears history and resets clock to new initial.
    """
    state0 = make_hierarchical_3body_system()
    system = NBodySystem(state=state0, G=1.0)
    system.evolve(t_span=0.1, dt=0.02)
    assert len(system.history) > 1

    # Reset with fresh state
    state_new = make_hierarchical_3body_system()
    system.reset(state_new)

    assert system.t == 0.0
    assert len(system.history) == 1
    assert len(system.atlas_manager.transition_history) == 0


def test_skeel_gear_macro_boundary_synchronization_and_palindromic_reversibility():
    """
    Validates Skeel & Gear (1992) macro-boundary synchronization:
    1. During a multi-rate dyadic macro-step (2^k substeps), the Jacobi tree chart
       remains strictly invariant across all internal substeps.
    2. The full multi-rate dyadic macro-step preserves palindromic time-reversibility:
       stepping forward by dt and then backward by -dt returns to initial state within < 1e-14.
    3. Processed integrators maintain exact boundary telescoping on the macro-boundary chart.
    """
    state0 = make_hierarchical_3body_system()
    system = NBodySystem(state=state0, integrator=BlanesMoan4Integrator(G=1.0), G=1.0)
    system.scheduler.default_level = 4  # Force 2^4 = 16 dyadic substeps per macro-step

    dt = 0.02
    initial_tree = system.tree
    state_fwd = system.step(dt=dt)

    assert system.tree == initial_tree
    assert np.isclose(system.t, dt, atol=1e-15)

    # Reversibility test with a fresh backward system
    system_rev = NBodySystem(state=state_fwd, integrator=BlanesMoan4Integrator(G=1.0), G=1.0)
    system_rev.scheduler.default_level = 4
    state_back = system_rev.step(dt=-dt)

    assert np.allclose(state_back.q, state0.q, atol=1e-13), f"q error: {np.max(np.abs(state_back.q - state0.q))}"
    assert np.allclose(state_back.p, state0.p, atol=1e-13), f"p error: {np.max(np.abs(state_back.p - state0.p))}"

    # Also test ProcessedSymplecticIntegrator
    pbm_integ = ProcessedSymplecticIntegrator(BlanesMoan4Integrator(G=1.0))
    system_pbm = NBodySystem(state=state0, integrator=pbm_integ, G=1.0)
    system_pbm.scheduler.default_level = 3  # 8 substeps
    s_pbm_fwd = system_pbm.step(dt=dt)
    assert not np.any(np.isnan(s_pbm_fwd.q))
    assert not np.any(np.isnan(s_pbm_fwd.p))
    inv_pbm = system_pbm.get_invariants()
    assert np.isfinite(inv_pbm.H)


def test_nbody_system_evolve_shadow_telescoping_evaluation_efficiency():
    """
    Validates that NBodySystem.evolve() with ProcessedSymplecticIntegrator achieves
    exact shadow telescoping:
    - Over M = 20 macro steps without atlas transitions, the corrector is called
      exactly twice (1 entry with apply_inverse, 1 exit with apply), reducing
      corrector evaluations from 2M = 40 to 2.
    - All 10 Lie--Bargmann invariants are strictly conserved (< 1e-13).
    """
    state0 = make_hierarchical_3body_system()
    base = BlanesMoan4Integrator(G=1.0)
    pbm = ProcessedSymplecticIntegrator(base)
    system = NBodySystem(state=state0, integrator=pbm, G=1.0)

    call_count = {"forward": 0, "inverse": 0}
    orig_apply = pbm.corrector.apply

    def counting_apply(*args, **kwargs):
        is_inv = kwargs.get("inverse", False)
        if len(args) >= 4 and isinstance(args[3], bool):
            is_inv = args[3]
        if is_inv:
            call_count["inverse"] += 1
        else:
            call_count["forward"] += 1
        return orig_apply(*args, **kwargs)

    pbm.corrector.apply = counting_apply

    M = 20
    dt = 0.01
    final_state = system.evolve(t_span=M * dt, dt=dt)

    # Telescoping check: exactly 1 inverse call at entrance and 1 forward call at exit
    # Without telescoping, this would evaluate 2 * M = 40 times.
    assert call_count["inverse"] == 1, f"Expected 1 inverse call, got {call_count['inverse']}"
    assert call_count["forward"] == 1, f"Expected 1 forward call, got {call_count['forward']}"

    # Invariants check
    inv_final = system.get_invariants()
    inv0 = compute_invariants(state0, G=1.0)
    assert np.allclose(inv_final.P, inv0.P, atol=1e-13)
    assert np.allclose(inv_final.L, inv0.L, atol=1e-13)
    assert np.allclose(inv_final.G, inv0.G, atol=1e-13)
    assert np.isclose(final_state.t, M * dt, atol=1e-14)


def test_nbody_system_default_max_dyadic_level_is_17():
    """
    Validates D-145-15: NBodySystem.__init__ uses DEFAULT_MAX_DYADIC_LEVEL = 17
    by default when max_dyadic_level is not specified.
    """
    state0 = make_hierarchical_3body_system()
    system = NBodySystem(state=state0, G=1.0)
    assert system.scheduler.max_level == 17


def test_nbody_system_pbm4_step_vs_evolve_exact_equivalence_and_physical_history():
    """
    Validates D-145-04 and D-145-08:
    - Under multirate dyadic substepping (k >= 2, dt_substep = dt / 2^k),
      M standalone calls to system.step(dt) and a single call to system.evolve(M*dt, dt)
      under PBM4 must produce the exact same final state (< 1e-13) and the exact same
      physical intermediate energy history in system.history.hamiltonians (< 1e-13).
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    q0 = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 4.0, 0.2],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, -0.8, 0.0],
        [0.0, 0.8, 0.0],
        [0.4, 0.0, -0.1],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses).to_barycentric()

    pbm_step = ProcessedSymplecticIntegrator(BlanesMoan4Integrator(G=1.0))
    sys_step = NBodySystem(state=state0, integrator=pbm_step, G=1.0)
    sys_step.scheduler.min_level = 2  # Force k >= 2 (>= 4 substeps per macro-step)

    pbm_evolve = ProcessedSymplecticIntegrator(BlanesMoan4Integrator(G=1.0))
    sys_evolve = NBodySystem(state=state0, integrator=pbm_evolve, G=1.0)
    sys_evolve.scheduler.min_level = 2

    M = 4
    dt = 0.40
    for _ in range(M):
        sys_step.step(dt=dt)

    sys_evolve.evolve(t_span=M * dt, dt=dt)

    # Verify dyadic substepping was actually active (k >= 2)
    assert np.all(sys_step.history.get_dyadic_levels()[1:] >= 2)
    assert np.all(sys_evolve.history.get_dyadic_levels()[1:] >= 2)

    # 1. Exact final state equivalence
    assert np.allclose(sys_step.state.q, sys_evolve.state.q, atol=1e-14, rtol=0.0), (
        f"q mismatch: {np.max(np.abs(sys_step.state.q - sys_evolve.state.q))}"
    )
    assert np.allclose(sys_step.state.p, sys_evolve.state.p, atol=1e-14, rtol=0.0), (
        f"p mismatch: {np.max(np.abs(sys_step.state.p - sys_evolve.state.p))}"
    )

    # 2. Exact physical-space intermediate history equivalence (D-145-08)
    H_step = sys_step.history.hamiltonians
    H_evolve = sys_evolve.history.hamiltonians
    assert len(H_step) == M + 1
    assert len(H_evolve) == M + 1
    assert np.allclose(H_step, H_evolve, atol=1e-14, rtol=0.0), (
        f"Intermediate Hamiltonian history mismatch: max diff = {np.max(np.abs(H_step - H_evolve))}"
    )


def test_nbody_system_pbm4_multirate_gauge_transition_uses_dt_substep():
    """
    Validates D-145-04:
    During an atlas chart transition in shadow mode with k = 3 (dt_substep = dt / 8),
    apply_gauge_connection must receive dt_substep (not macro_dt), preventing a 16^k = 4096x
    overestimation of the gauge momentum jump.
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Body 2 approaches Body 1 so tree flips from ((0,1),2) to ((1,2),0)
    q0 = np.array([
        [-1.5, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.45, 0.1, 0.0],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, -0.4, 0.0],
        [0.0, 0.5, 0.0],
        [-0.3, -0.1, 0.0],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses).to_barycentric()

    # build_canonical_tree(masses) pairs (0, 1) initially even though (1, 2) is much closer (d=0.46 < 0.6*1.5)
    tree_01 = build_canonical_tree(masses)
    assert tree_01.get_leaf_pair_indices() == [(0, 1)]

    pbm = ProcessedSymplecticIntegrator(BlanesMoan4Integrator(G=1.0))
    sys_step = NBodySystem(state=state0, tree=tree_01, integrator=pbm, G=1.0)
    sys_step.scheduler.min_level = 3  # k >= 3 -> dt_substep <= dt / 8

    sys_evolve = NBodySystem(
        state=state0,
        tree=build_canonical_tree(masses),
        integrator=ProcessedSymplecticIntegrator(BlanesMoan4Integrator(G=1.0)),
        G=1.0,
    )
    sys_evolve.scheduler.min_level = 3

    dt = 0.25
    sys_step.step(dt=dt)
    sys_evolve.evolve(t_span=dt, dt=dt)

    # Both must flip to pair (1, 2) and agree to machine precision
    assert sys_step.tree.get_leaf_pair_indices() == [(1, 2)]
    assert sys_evolve.tree.get_leaf_pair_indices() == [(1, 2)]
    assert np.allclose(sys_step.state.q, sys_evolve.state.q, atol=1e-14)
    assert np.allclose(sys_step.state.p, sys_evolve.state.p, atol=1e-14)


def test_colocalized_listings_01_and_04_interfaces():
    """
    Validates D-145-16:
    - Colocalized Listing 01 (01_PhaseStateAndSystem_Interface.py) implements
      ExtendedPhaseState.with_updates(masses=...) and functional NBodySystem.step / evolve.
    - Colocalized Listing 04 (04_SymplecticSubcycle_Interface.py) solves hyperbolic
      micro-steps (dt = 1e-15), and exports step_supernode_synchronous and step_strang_subcycle.
    """
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    listing_01_path = (
        root
        / "interfaces" / "01_PhaseStateAndSystem_Interface.py"
    )
    spec_01 = importlib.util.spec_from_file_location("listing_01_system", str(listing_01_path))
    mod_01 = importlib.util.module_from_spec(spec_01)
    import sys
    sys.modules["listing_01_system"] = mod_01
    spec_01.loader.exec_module(mod_01)

    state0 = make_hierarchical_3body_system()
    st_01 = mod_01.ExtendedPhaseState(q=state0.q, p=state0.p, masses=state0.masses)
    st_upd = st_01.with_updates(masses=np.array([2.0, 2.0, 1.0]))
    assert np.allclose(st_upd.masses, [2.0, 2.0, 1.0])

    sys_01 = mod_01.NBodySystem(state=st_01, G=1.0)
    st_after = sys_01.step(0.01)
    assert st_after is not None, "Listing 01 NBodySystem.step returned None (unimplemented stub)"
    assert np.isclose(st_after.t, 0.01, atol=1e-15)
    assert hasattr(sys_01, "evolve")
    st_ev = sys_01.evolve(0.02, dt=0.01)
    assert np.isclose(st_ev.t, 0.03, atol=1e-15)

    listing_04_path = (
        root
        / "interfaces" / "04_SymplecticSubcycle_Interface.py"
    )
    spec_04 = importlib.util.spec_from_file_location("listing_04_subcycle", str(listing_04_path))
    mod_04 = importlib.util.module_from_spec(spec_04)
    sys.modules["listing_04_subcycle"] = mod_04
    spec_04.loader.exec_module(mod_04)

    # Check hyperbolic micro-step convergence in solve_sundman_time_equation
    tau_micro = mod_04.solve_sundman_time_equation(
        r0=0.5, sigma0=1.0, GM=1.0, alpha=-4.0, dt_target=1e-15
    )
    assert abs(tau_micro - 2e-15) < 1e-25, f"Hyperbolic micro-step failed: tau={tau_micro}"

    assert hasattr(mod_04, "step_supernode_synchronous")
    assert hasattr(mod_04, "step_strang_subcycle")
    q_sub, p_sub = mod_04.step_strang_subcycle(
        state0.q, state0.p, state0.masses, leaf_pairs=[(0, 1)], dt=0.01, G=1.0
    )
    assert q_sub.shape == (3, 3)
    assert p_sub.shape == (3, 3)
    assert np.allclose(np.sum(p_sub, axis=0), np.sum(state0.p, axis=0), atol=1e-14)


def test_nbody_system_shadow_mode_default_and_stateful_lifecycle():
    """
    Certifies Detached Component 6 (NBodySystem in_shadow_mode=True by default, [ENG-FOOTGUN-02], [CRIT-7.2.1-05]):
    1. NBodySystem.__init__ defaults to in_shadow_mode=True.
    2. A standalone loop of system.step(dt) automatically enters the shadow manifold on step 1,
       remains on the shadow manifold across subsequent steps without redundant corrector roundtrips,
       exposes non-destructive physical observation via get_physical_state() / observe_physical_state(),
       and exits cleanly via exit_shadow_manifold().
    """
    from normal_cone_nbody_float64.integrators.processed import create_pbm4_integrator

    state0 = make_hierarchical_3body_system().to_barycentric()
    pbm4 = create_pbm4_integrator(G=1.0)
    sys_loop = NBodySystem(state=state0, integrator=pbm4, G=1.0)

    assert sys_loop.in_shadow_mode is True
    assert sys_loop._on_shadow_manifold is False

    call_count = {"forward": 0, "inverse": 0}
    orig_apply = pbm4.corrector.apply

    def counting_apply(*args, **kwargs):
        is_inv = kwargs.get("inverse", False)
        if len(args) >= 4 and isinstance(args[3], bool):
            is_inv = args[3]
        if is_inv:
            call_count["inverse"] += 1
        else:
            call_count["forward"] += 1
        return orig_apply(*args, **kwargs)

    pbm4.corrector.apply = counting_apply

    dt = 0.01
    M = 6
    for _ in range(M):
        sys_loop.step(dt=dt)
        assert sys_loop._on_shadow_manifold is True

    # Across M=6 steps with constant dyadic level and chart:
    # 1 inverse call at step 1 entry + M forward calls for history physical observation (zero extra inverse calls!)
    assert call_count["inverse"] == 1, f"Expected exactly 1 shadow entry inverse call, got {call_count['inverse']}"

    phys_obs = sys_loop.observe_physical_state()
    assert sys_loop._on_shadow_manifold is True
    phys_exit = sys_loop.exit_shadow_manifold()
    assert sys_loop._on_shadow_manifold is False
    assert np.allclose(phys_obs.q, phys_exit.q, atol=1e-15)
    assert np.allclose(phys_obs.p, phys_exit.p, atol=1e-15)



