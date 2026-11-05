"""
Unified Normal-Cone N-Body Integration: Initialization and Shadow Pre-Processing Interface.
Mathematical Monograph Alignment: Chapter 6, Section 6.2 (Phases 0 and 1, Appendix A.23).
Reference Codebase: normal_cone_nbody_float64/system/atlas_manager.py
                   normal_cone_nbody_float64/potential/supernode_detector.py
                   normal_cone_nbody_float64/integrators/correctors.py
                   normal_cone_nbody_float64/integrators/processed.py

Reference Implementation Interface:
1. build_optimal_tree: Delegates directly to AtlasManager(gamma_build=gamma).build_optimal_tree(state).
2. detect_supernodes, compute_natural_threshold, compute_internal_hierarchy_ratio,
   compute_cubic_tidal_force_ratio: Delegate directly to potential.supernode_detector.
3. Symplectic Corrector Hierarchy (Theorems 5.4, 5.5 & 5.6):
   - PalindromicCorrectorZ4 = HomologicalLieCorrectorC4 (unobstructed 4-shell palindromic
     Lie conjugator Z_4^{BM4}(h) in P_1, canonical default of CorrectorC4);
   - CotangentLiftCorrectorP1 = CotangentLiftCorrectorC4 (unobstructed cotangent-lifted
     flow exp(-X_{K_4}) in P_1 with 2-stage Gauss-Legendre symplectic collocation);
   - SymplecticBoundaryCorrector = CorrectorC4 (unified corrector dispatcher defaulting
     to mode="homological");
   - GlobalCorrectorC4 (strictly labeled configuration-only P_0 ablation baseline,
     obstructed on eccentric/3D orbits per Theorem 5.4).
"""

from __future__ import annotations
import sys
from pathlib import Path
from typing import Optional, List, Tuple, Set, Union

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import (
    JacobiTree,
    JacobiTopologyNode,
)
from normal_cone_nbody_float64.core.constants import (
    ETA_CRITICAL_HILL,
    ETA_ENTER_HILL,
    ETA_EXIT_HILL,
    LAMBDA_CRITICAL_ROCHE,
    KAPPA4_BM4,
)
from normal_cone_nbody_float64.system.atlas_manager import AtlasManager
from normal_cone_nbody_float64.potential.supernode_detector import (
    compute_natural_threshold,
    compute_internal_hierarchy_ratio,
    compute_cubic_tidal_force_ratio,
    detect_supernodes as _canonical_detect_supernodes,
)
from normal_cone_nbody_float64.integrators.correctors import (
    CorrectorC4,
    GlobalCorrectorC4,
    CotangentLiftCorrectorC4,
    HomologicalLieCorrectorC4,
    compute_conormal_tidal_tensor as compute_normal_cone_tidal_tensor,
    compute_conormal_corrector_kick as compute_normal_cone_corrector_kick,
    compute_global_conormal_c4_kick as compute_global_normal_cone_c4_kick,
    compute_cotangent_lift_c4_flow,
)
compute_cotangent_generator_K4 = compute_cotangent_lift_c4_flow
from normal_cone_nbody_float64.integrators.processed import (
    ProcessedSymplecticIntegrator,
    create_pbm4_integrator,
)


def build_optimal_tree(state: ExtendedPhaseState, gamma: float = 1.0) -> JacobiTree:
    """
    Constructs an optimal binary Jacobi tree T_0 in Trees(N) via agglomerative
    hierarchical clustering on mutual Cartesian barycenters (Phase 0.2), equipped
    with the unilateral Delaunay-Hill isolation guard (gamma_build = gamma).
    """
    eff_gamma_monitor = 1.0 if float(gamma) >= 1.0 - 1e-15 else min(0.60, float(gamma))
    manager = AtlasManager(gamma=eff_gamma_monitor, gamma_build=float(gamma))
    return manager.build_optimal_tree(state)


def detect_supernodes(
    state: ExtendedPhaseState,
    tree: JacobiTree,
    eta_enter: float = ETA_ENTER_HILL,
    eta_exit: float = ETA_EXIT_HILL,
    active_supernodes: Optional[Union[Set[Tuple[int, ...]], List[Tuple[int, ...]]]] = None,
    G: float = 1.0,
    mode: str = "canonical",
    lambda_crit: float = LAMBDA_CRITICAL_ROCHE,
    enforce_non_hierarchical: bool = True,
) -> List[Tuple[int, ...]]:
    """
    Proactively detects compact multi-body clusters (k >= 3) via the canonical
    Schmitt trigger band [eta_enter, eta_exit] = [2^(-1/3), 2^(-1/6)] in Hill
    bound mode or cubic tidal ratio Lambda_S >= Lambda_crit in Roche disruption mode.
    """
    eff_mode = "hill_bound" if mode == "hill_isolation" else mode
    return _canonical_detect_supernodes(
        state=state,
        tree=tree,
        eta_enter=eta_enter,
        eta_exit=eta_exit,
        active_supernodes=active_supernodes,
        G=G,
        mode=eff_mode,
        lambda_crit=lambda_crit,
        enforce_non_hierarchical=enforce_non_hierarchical,
    )


# Canonical Cotangent-Lifted Symplectic Corrector Exports (P_1):
PalindromicCorrectorZ4 = HomologicalLieCorrectorC4
CotangentLiftCorrectorP1 = CotangentLiftCorrectorC4
SymplecticBoundaryCorrector = CorrectorC4

__all__ = [
    "ExtendedPhaseState",
    "JacobiTree",
    "JacobiTopologyNode",
    "ETA_CRITICAL_HILL",
    "ETA_ENTER_HILL",
    "ETA_EXIT_HILL",
    "LAMBDA_CRITICAL_ROCHE",
    "KAPPA4_BM4",
    "compute_natural_threshold",
    "compute_internal_hierarchy_ratio",
    "compute_cubic_tidal_force_ratio",
    "build_optimal_tree",
    "detect_supernodes",
    "PalindromicCorrectorZ4",
    "CotangentLiftCorrectorP1",
    "SymplecticBoundaryCorrector",
    "CorrectorC4",
    "HomologicalLieCorrectorC4",
    "CotangentLiftCorrectorC4",
    "GlobalCorrectorC4",
    "ProcessedSymplecticIntegrator",
    "create_pbm4_integrator",
    "compute_normal_cone_tidal_tensor",
    "compute_normal_cone_corrector_kick",
    "compute_global_normal_cone_c4_kick",
    "compute_cotangent_lift_c4_flow",
    "compute_cotangent_generator_K4",
]


