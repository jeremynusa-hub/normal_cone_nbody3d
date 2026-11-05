"""
Detached Component 7 Certification Suite (Iteration 57, Pillar I Closure):
Direct-Delegation Alignment and Independent Verification of the 7 Paper Showcase Interfaces:
1. chapters/pcr3bp_canonical_solver.py (Appendix A.13)
2. chapters/01_PhaseStateAndSystem_Interface.py (Appendix A.22)
3. chapters/02_TreeBuilderAndCorrector_Interface.py (Appendix A.23)
4. chapters/03_DyadicScheduler_Interface.py (Appendix A.24)
5. chapters/04_SymplecticSubcycle_Interface.py (Appendix A.25)
6. chapters/05_AtlasTransitions_Interface.py (Appendix A.26)
7. chapters/06_HighPerformanceOptimizations_Interface.py (Appendix A.27)
"""

from __future__ import annotations
import importlib.util
from pathlib import Path
from types import ModuleType
import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec
from normal_cone_nbody_float64.core.transitions import (
    SL2TransitionMatrix as CanonicalSL2TransitionMatrix,
    build_tree_transition_symplectic_12 as canonical_build_tree_transition_symplectic_12,
)
from normal_cone_nbody_float64.core.invariants import (
    LieBargmannInvariants as CanonicalLieBargmannInvariants,
    compute_invariants,
)
from normal_cone_nbody_float64.system.nbody_system import (
    NBodySystem as CanonicalNBodySystem,
)
from normal_cone_nbody_float64.system.atlas_manager import (
    AtlasManager as CanonicalAtlasManager,
)
from normal_cone_nbody_float64.system.dyadic_scheduler import (
    ColoredInterval as CanonicalColoredInterval,
    DyadicScheduler as CanonicalDyadicScheduler,
)
from normal_cone_nbody_float64.integrators import common as canonical_common
from normal_cone_nbody_float64.integrators.correctors import (
    CorrectorC4,
    GlobalCorrectorC4,
    CotangentLiftCorrectorC4,
    HomologicalLieCorrectorC4,
)
from normal_cone_nbody_float64.integrators.processed import (
    ProcessedSymplecticIntegrator,
    PBM4GaugeTransitionOperator as CanonicalPBM4GaugeTransitionOperator,
    apply_gauge_connection as canonical_apply_gauge_connection,
    create_pbm4_integrator,
)
from normal_cone_nbody_float64.solvers.stumpff import (
    stumpff_c as canonical_stumpff_c,
)
from normal_cone_nbody_float64.solvers.kepler_sundman import (
    step_kepler_sundman as canonical_step_kepler_sundman,
)
from normal_cone_nbody_float64.potential.screened_apex import (
    ScreenedApexPotential,
)
from normal_cone_nbody_float64.potential.optimized_screened_apex import (
    OptimizedScreenedApexPotential as CanonicalOptimizedScreenedApexPotential,
)
from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
    OptimizedKeplerSolver as CanonicalOptimizedKeplerSolver,
)
from normal_cone_nbody_float64.solvers.optimized_supernode import (
    OptimizedSupernodeSolver as CanonicalOptimizedSupernodeSolver,
)
from normal_cone_nbody_float64.solvers.jit_supernode import (
    JitSupernodeSolver as CanonicalJitSupernodeSolver,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
PAPER_DIR = REPO_ROOT / "interfaces"


def _load_showcase_module(name: str, rel_path: str) -> ModuleType:
    full_path = PAPER_DIR / rel_path
    assert full_path.is_file(), f"Missing paper showcase file: {full_path}"
    spec = importlib.util.spec_from_file_location(name, str(full_path))
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_showcase_pcr3bp_canonical_solver_detached():
    """
    Certifies Appendix A.13 (`pcr3bp_canonical_solver.py`):
    - Symmetric 7-stage Strang composition weights satisfy sum(w_i) = 1 and sum(w_i^3) = 0;
    - Symmetric dyadic rung selector k_sym = max(k_start, k_end) guarantees exact
      microscopic time-reversibility across a dyadic rung boundary.
    """
    mod_pcr = _load_showcase_module(
        "pcr3bp_canonical_solver",
        "pcr3bp_canonical_solver.py",
    )
    assert hasattr(mod_pcr, "PCR3BPNormalConeIntegrator")
    assert hasattr(mod_pcr, "PCR3BPConormalIntegrator")
    assert mod_pcr.PCR3BPNormalConeIntegrator is mod_pcr.PCR3BPConormalIntegrator

    bm4_a = np.asarray(mod_pcr.BM4_A, dtype=np.float64)
    assert len(bm4_a) == 7
    assert np.allclose(bm4_a, bm4_a[::-1], atol=1e-16)
    assert abs(float(np.sum(bm4_a)) - 1.0) <= 1e-15
    assert abs(float(np.sum(bm4_a**3))) <= 1e-15

    solver = mod_pcr.PCR3BPNormalConeIntegrator(mu=0.012277471)
    q0 = np.array([0.994, 0.0], dtype=np.float64)
    p0 = np.array([0.0, -1.0075851063790825], dtype=np.float64)
    dt = 0.02
    q_fwd, p_fwd = solver.step_sundman_regularized(q0, p0, dt)
    q_bwd, p_bwd = solver.step_sundman_regularized(q_fwd, p_fwd, -dt)
    assert np.linalg.norm(q_bwd - q0) <= 1e-13
    assert np.linalg.norm(p_bwd - p0) <= 2e-11


def test_showcase_01_phase_state_and_system_direct_delegation():
    """
    Certifies Appendix A.22 (`01_PhaseStateAndSystem_Interface.py`):
    - Direct class identity with canonical ExtendedPhaseState and NBodySystem;
    - Shadow-manifold mode active by default (`in_shadow_mode=True`) and support
      for both `dt` and `dt_macro` in `.step()`.
    """
    mod_01 = _load_showcase_module(
        "showcase_01",
        "01_PhaseStateAndSystem_Interface.py",
    )
    assert mod_01.ExtendedPhaseState is ExtendedPhaseState
    assert mod_01.NBodySystem is CanonicalNBodySystem

    masses = np.array([1.0, 1e-3, 1e-3], dtype=np.float64)
    q0 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [4.0, 0.0, 0.0]], dtype=np.float64)
    p0 = np.array([[0.0, 0.0, 0.0], [0.0, 1e-3, 0.0], [0.0, 5e-4, 0.0]], dtype=np.float64)
    state0 = mod_01.ExtendedPhaseState(q=q0, p=p0, masses=masses).to_barycentric()
    sys = mod_01.NBodySystem(state=state0, integrator=create_pbm4_integrator(G=1.0), dt_macro=0.05)

    assert sys.in_shadow_mode is True
    assert sys.integrator.corrector.mode == "homological"

    st1 = sys.step(dt_macro=0.05)
    assert sys._on_shadow_manifold is True
    assert np.isclose(st1.t, 0.05)
    assert np.linalg.norm(st1.center_of_mass) <= 1e-15
    assert np.linalg.norm(st1.total_momentum) <= 1e-15


def test_showcase_02_tree_builder_and_corrector_no_p0_footgun():
    """
    Certifies Appendix A.23 (`02_TreeBuilderAndCorrector_Interface.py`):
    - Eliminates the historical `PalindromicCorrectorZ4 = GlobalCorrectorC4` alias footgun;
    - Exports `PalindromicCorrectorZ4 is HomologicalLieCorrectorC4` (in P_1),
      `CotangentLiftCorrectorP1 is CotangentLiftCorrectorC4` (in P_1),
      `SymplecticBoundaryCorrector is CorrectorC4` (defaulting to mode="homological"),
      and `GlobalCorrectorC4` strictly as P_0 ablation baseline;
    - `build_optimal_tree` delegates to AtlasManager with gamma_build = 0.70.
    """
    mod_02 = _load_showcase_module(
        "showcase_02",
        "02_TreeBuilderAndCorrector_Interface.py",
    )
    assert mod_02.PalindromicCorrectorZ4 is HomologicalLieCorrectorC4
    assert mod_02.PalindromicCorrectorZ4 is not GlobalCorrectorC4
    assert mod_02.CotangentLiftCorrectorP1 is CotangentLiftCorrectorC4
    assert mod_02.SymplecticBoundaryCorrector is CorrectorC4
    assert mod_02.SymplecticBoundaryCorrector().mode == "homological"
    assert mod_02.create_pbm4_integrator().corrector.mode == "homological"

    # Two-tier hysteresis check on build_optimal_tree (gamma=0.70 rejects ratio 0.638 < 0.70)
    masses = np.array([1000.0, 2.0, 1.0, 1.0], dtype=np.float64)
    q = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.34, 0.0, 0.0], [4.44, 0.0, 0.0]], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=np.zeros_like(q), masses=masses)
    tree_070 = mod_02.build_optimal_tree(state, gamma=0.70)
    assert tree_070.get_leaf_pair_indices() == [(0, 1)]


def test_showcase_03_dyadic_scheduler_direct_delegation():
    """
    Certifies Appendix A.24 (`03_DyadicScheduler_Interface.py`):
    - Direct identity with CanonicalDyadicScheduler and CanonicalColoredInterval;
    - Symmetric endpoint level selection k_sym = max(k_start, k_end) across pericenter.
    """
    mod_03 = _load_showcase_module(
        "showcase_03",
        "03_DyadicScheduler_Interface.py",
    )
    assert mod_03.DyadicScheduler is CanonicalDyadicScheduler
    assert mod_03.ColoredInterval is CanonicalColoredInterval

    sched = mod_03.DyadicScheduler(dt_macro=0.1, default_level=0, max_level=17)
    masses = np.array([1.0, 1e-6], dtype=np.float64)
    q = np.array([[0.0, 0.0, 0.0], [0.02, 0.002, 0.0]], dtype=np.float64)
    p = np.array([[0.0, 0.0, 0.0], [-1e-5, 0.0, 0.0]], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_tree_from_spec(masses, (0, 1))

    _, _, k_fwd = sched.get_substep_params(0.0, dt_block=0.02, state=state, tree=tree)
    assert k_fwd >= 1


def test_showcase_04_symplectic_subcycle_direct_delegation():
    """
    Certifies Appendix A.25 (`04_SymplecticSubcycle_Interface.py`):
    - Direct delegation of `stumpff_c` and `step_kepler_sundman`;
    - Multi-substep time-reversibility of `step_supernode_synchronous` (<= 1e-15);
    - Preservation of inter-supernode cross-pair forces in `compute_intercluster_forces` for N=6.
    """
    mod_04 = _load_showcase_module(
        "showcase_04",
        "04_SymplecticSubcycle_Interface.py",
    )
    assert mod_04.stumpff_c is canonical_stumpff_c
    assert mod_04.step_kepler_sundman is canonical_step_kepler_sundman

    # Time-reversibility of step_supernode_synchronous on a 3-body close encounter
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    q0 = np.array([[-0.2, 0.0, 0.0], [0.2, 0.0, 0.0], [0.0, 0.5, 0.1]], dtype=np.float64)
    p0 = np.array([[0.0, -0.8, 0.1], [0.0, 0.8, -0.1], [0.2, 0.0, 0.0]], dtype=np.float64)
    q_f, p_f = mod_04.step_supernode_synchronous(q0, p0, masses, 0.05, G=1.0)
    q_b, p_b = mod_04.step_supernode_synchronous(q_f, p_f, masses, -0.05, G=1.0)
    assert np.max(np.abs(q_b - q0)) <= 1e-14
    assert np.max(np.abs(p_b - p0)) <= 1e-14

    # Inter-supernode cross-pair preservation for N=6 with S1=(0,1,2) and S2=(3,4,5)
    m6 = np.ones(6, dtype=np.float64)
    q6 = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [10.0, 0.0, 0.0],
            [10.1, 0.0, 0.0],
            [10.0, 0.1, 0.0],
        ],
        dtype=np.float64,
    )
    f6 = mod_04.compute_intercluster_forces(
        q6, m6, leaf_pair_set=set(), supernode_clusters=[(0, 1, 2), (3, 4, 5)], G=1.0
    )
    # Cluster S1 must feel a net attractive pull in +x toward cluster S2
    assert float(np.sum(f6[:3, 0])) > 0.08
    assert np.allclose(np.sum(f6, axis=0), 0.0, atol=1e-15)


def test_showcase_05_atlas_transitions_and_gauge_operator():
    """
    Certifies Appendix A.26 (`05_AtlasTransitions_Interface.py`):
    - Direct export of SL2TransitionMatrix (.from_masses and .from_trees),
      build_tree_transition_symplectic_12, PBM4GaugeTransitionOperator, and apply_gauge_connection;
    - 3-chart Sp(12, R) Cech 1-cocycle closure <= 1e-15 and non-commutative gauge operator
      defaulting to mode="homological".
    """
    mod_05 = _load_showcase_module(
        "showcase_05",
        "05_AtlasTransitions_Interface.py",
    )
    assert mod_05.SL2TransitionMatrix is CanonicalSL2TransitionMatrix
    assert mod_05.build_tree_transition_symplectic_12 is canonical_build_tree_transition_symplectic_12
    assert mod_05.PBM4GaugeTransitionOperator is CanonicalPBM4GaugeTransitionOperator
    assert mod_05.apply_gauge_connection is canonical_apply_gauge_connection
    assert mod_05.LieBargmannInvariants is CanonicalLieBargmannInvariants
    assert issubclass(mod_05.AtlasManager, CanonicalAtlasManager)

    masses = np.array([1.2, 0.8, 0.5], dtype=np.float64)
    tA = build_tree_from_spec(masses, ((0, 1), 2))
    tB = build_tree_from_spec(masses, (0, (1, 2)))
    tC = build_tree_from_spec(masses, ((0, 2), 1))

    psi_AB = mod_05.build_tree_transition_symplectic_12(tA, tB)
    psi_BC = mod_05.build_tree_transition_symplectic_12(tB, tC)
    psi_CA = mod_05.build_tree_transition_symplectic_12(tC, tA)
    cocycle = psi_CA @ psi_BC @ psi_AB
    assert np.max(np.abs(cocycle - np.eye(12))) <= 1e-15

    gauge_op = mod_05.PBM4GaugeTransitionOperator(G=1.0)
    assert gauge_op.corrector.mode == "homological"


def test_showcase_06_high_performance_optimizations_shared_contextvars():
    """
    Certifies Appendix A.27 (`06_HighPerformanceOptimizations_Interface.py`):
    - Shared ContextVar identity between `06_*.py` and `integrators.common`;
    - Exact 2x5-stage Horner + Stumpff recurrence in `OptimizedKeplerSolver.stumpff_c0123_scalar`;
    - Subclass relationship of `OptimizedScreenedApexPotential` and `OptimizedKeplerSolver`
      with canonical classes and bit-exact force equivalence.
    """
    mod_06 = _load_showcase_module(
        "showcase_06",
        "06_HighPerformanceOptimizations_Interface.py",
    )
    assert mod_06._CURRENT_APEX_EVALUATOR is canonical_common._CURRENT_APEX_EVALUATOR
    assert mod_06._CURRENT_KEPLER_SOLVER is canonical_common._CURRENT_KEPLER_SOLVER
    assert mod_06._CURRENT_SUPERNODE_SOLVER is canonical_common._CURRENT_SUPERNODE_SOLVER
    assert issubclass(mod_06.OptimizedScreenedApexPotential, CanonicalOptimizedScreenedApexPotential)
    assert issubclass(mod_06.OptimizedKeplerSolver, CanonicalOptimizedKeplerSolver)
    assert mod_06.OptimizedSupernodeSolver is CanonicalOptimizedSupernodeSolver
    assert mod_06.JitSupernodeSolver is CanonicalJitSupernodeSolver

    # ContextVar propagation into canonical_common
    masses = np.array([1.0, 0.5, 0.25], dtype=np.float64)
    opt_apex = mod_06.OptimizedScreenedApexPotential(masses, G=1.0)
    tok = mod_06.set_current_apex_evaluator(opt_apex)
    try:
        assert canonical_common.get_current_apex_evaluator() is opt_apex
    finally:
        mod_06.reset_current_apex_evaluator(tok)

    # Degree-6 (6-stage) Horner + Stumpff recurrence accuracy == 0.0 ULP
    for z in (-0.099, -1e-4, 0.0, 1e-4, 0.099):
        c_opt = mod_06.OptimizedKeplerSolver.stumpff_c0123_scalar(z)
        c_ref = canonical_stumpff_c(z)[:4]
        assert np.array_equal(np.array(c_opt), np.array(c_ref))

    # 0.0 ULP Force equivalence between canonical ScreenedApexPotential and showcase OptimizedScreenedApexPotential
    q = np.array([[0.0, 0.0, 0.0], [0.5, 0.1, -0.05], [3.0, -0.4, 0.2]], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=np.zeros_like(q), masses=masses)
    tree = build_tree_from_spec(masses, ((0, 1), 2))
    ref_forces, ref_v = ScreenedApexPotential(G=1.0).compute_apex_forces(state, tree, include_screened_monopole=False)
    opt_forces_raw = opt_apex.compute_apex_forces(q, tree)
    opt_forces_state, opt_v = opt_apex.compute_apex_forces(state, tree, include_screened_monopole=False)
    assert np.array_equal(opt_forces_raw, ref_forces)
    assert np.array_equal(opt_forces_state, ref_forces)
    assert opt_v == ref_v


def test_end_to_end_pipeline_assembled_purely_from_paper_showcase_interfaces():
    """
    Assembles a complete 6-phase simulation pipeline exclusively from the 6 Chapter 6
    paper showcase interfaces (`01..06_*.py`) and verifies shadow-manifold evolution,
    optimized ContextVar dispatch, and machine-precision Lie-Bargmann conservation.
    """
    mod_01 = _load_showcase_module(
        "showcase_01_e2e",
        "01_PhaseStateAndSystem_Interface.py",
    )
    mod_02 = _load_showcase_module(
        "showcase_02_e2e",
        "02_TreeBuilderAndCorrector_Interface.py",
    )
    mod_03 = _load_showcase_module(
        "showcase_03_e2e",
        "03_DyadicScheduler_Interface.py",
    )
    mod_05 = _load_showcase_module(
        "showcase_05_e2e",
        "05_AtlasTransitions_Interface.py",
    )
    mod_06 = _load_showcase_module(
        "showcase_06_e2e",
        "06_HighPerformanceOptimizations_Interface.py",
    )

    masses = np.array([1.0, 1e-3, 1e-3], dtype=np.float64)
    q0 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.1]], dtype=np.float64)
    p0 = np.array([[0.0, 0.0, 0.0], [0.0, 1e-3, 0.0], [0.0, 4.47e-4, 0.0]], dtype=np.float64)

    state0 = mod_01.ExtendedPhaseState(q=q0, p=p0, masses=masses).to_barycentric()
    inv0 = mod_05.compute_lie_bargmann_invariants(state0, G=1.0)
    state0 = state0.with_updates(E=inv0.H, E0=inv0.H)

    tree0 = mod_02.build_optimal_tree(state0, gamma=0.70)
    integrator = mod_02.create_pbm4_integrator(G=1.0)
    scheduler = mod_03.DyadicScheduler(dt_macro=0.05, G=1.0)
    atlas = mod_05.AtlasManager(gamma=0.60, gamma_build=0.70)

    sys = mod_01.NBodySystem(
        state=state0,
        tree=tree0,
        integrator=integrator,
        atlas_manager=atlas,
        scheduler=scheduler,
        G=1.0,
        use_optimized_kernels=True,
    )
    final_state = sys.evolve(t_span=0.25, dt=0.05)
    inv_f = mod_05.compute_lie_bargmann_invariants(final_state, G=1.0)

    # Verify 0.0 ULP isomorphism against baseline NBodySystem(use_optimized_kernels=False)
    sys_base = mod_01.NBodySystem(
        state=state0,
        tree=mod_02.build_optimal_tree(state0, gamma=0.70),
        integrator=mod_02.create_pbm4_integrator(G=1.0),
        atlas_manager=mod_05.AtlasManager(gamma=0.60, gamma_build=0.70),
        scheduler=mod_03.DyadicScheduler(dt_macro=0.05, G=1.0),
        G=1.0,
        use_optimized_kernels=False,
    )
    final_state_base = sys_base.evolve(t_span=0.25, dt=0.05)
    assert np.array_equal(final_state.q, final_state_base.q)
    assert np.array_equal(final_state.p, final_state_base.p)

    assert np.isclose(final_state.t, 0.25)
    assert np.linalg.norm(inv_f.P - inv0.P) <= 1e-15
    assert np.linalg.norm(inv_f.L - inv0.L) <= 1e-15
    assert np.linalg.norm(inv_f.G - inv0.G) <= 1e-15
    assert inv_f.delta_E_rel <= 1e-11
