"""
Section 6.6 Standalone Interface Specification (Appendix A.27):
High-Performance Computational Optimizations and Modular Micro-Kernels

Reference Implementation Interface:
1. Unified Thread-Local ContextVar Dispatch:
   Re-exports _CURRENT_APEX_EVALUATOR, _CURRENT_KEPLER_SOLVER, and
   _CURRENT_SUPERNODE_SOLVER (and their set/get/reset helpers) directly from
   normal_cone_nbody_float64.integrators.common, guaranteeing
   that setting an optimized kernel via this interface is immediately active
   inside step_apex_kick, step_drift_and_kepler, and step_supernode_flow.
2. OptimizedScreenedApexPotential:
   Subclasses canonical potential.optimized_screened_apex.OptimizedScreenedApexPotential
   (supporting supernode_clusters and r^2 > 0 guards) while accepting both
   OptimizedScreenedApexPotential(G=1.0) and OptimizedScreenedApexPotential(masses, G=1.0),
   and both compute_apex_forces(state, tree) and compute_apex_forces(q, tree).
3. OptimizedKeplerSolver / OptimizedKeplerSundmanSolver:
   Delegates directly to solvers.optimized_kepler_sundman (degree-6 (6-stage) Horner on c_2, c_3
   with exact Stumpff recurrence c_0 = 1 - z*c_2, c_1 = 1 - z*c_3 and incremental
   Gauss-Stumpff conditioning).
4. OptimizedSupernodeSolver / JitSupernodeSolver:
   Re-exports canonical solvers.optimized_supernode and solvers.jit_supernode.
"""

from __future__ import annotations
import sys
from pathlib import Path
from typing import Any, FrozenSet, List, Optional, Sequence, Set, Tuple, Union
import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import JacobiTree
from normal_cone_nbody_float64.integrators.common import (
    _CURRENT_APEX_EVALUATOR,
    _CURRENT_KEPLER_SOLVER,
    _CURRENT_SUPERNODE_SOLVER,
    _CURRENT_SUPERNODE_N_SUBSTEPS,
    set_current_apex_evaluator,
    get_current_apex_evaluator,
    reset_current_apex_evaluator,
    set_current_kepler_solver,
    get_current_kepler_solver,
    reset_current_kepler_solver,
    set_current_supernode_solver,
    get_current_supernode_solver,
    reset_current_supernode_solver,
    set_current_supernode_n_substeps,
    get_current_supernode_n_substeps,
    reset_current_supernode_n_substeps,
)
from normal_cone_nbody_float64.potential.screened_apex import (
    _normalize_supernode_clusters,
)
from normal_cone_nbody_float64.potential.optimized_screened_apex import (
    OptimizedScreenedApexPotential as _CanonicalOptimizedScreenedApexPotential,
    _HAS_NUMBA,
    _njit_accumulate_apex_pairs,
)
from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
    OptimizedKeplerSolver as _CanonicalOptimizedKeplerSolver,
    stumpff_c0123_scalar,
    step_kepler_sundman_fast,
)
from normal_cone_nbody_float64.solvers.optimized_supernode import (
    OptimizedSupernodeSolver,
    step_supernode_synchronous_fast,
)
from normal_cone_nbody_float64.solvers.jit_supernode import (
    JitSupernodeSolver,
    step_supernode_synchronous_jit,
)


class OptimizedScreenedApexPotential(_CanonicalOptimizedScreenedApexPotential):
    """
    Tensor-free screened apex potential subclassing the canonical
    OptimizedScreenedApexPotential. Supports both:
    - Canonical initialization: OptimizedScreenedApexPotential(G=1.0, screening_mode=None, order=None)
    - Pre-bound mass initialization: OptimizedScreenedApexPotential(masses, G=1.0)
    and both:
    - Canonical call: compute_apex_forces(state, tree, ...) -> (forces, V_apex)
    - Raw array call: compute_apex_forces(q, tree, ...) -> forces
    """

    def __init__(
        self,
        masses_or_G: Union[np.ndarray, List[float], float] = 1.0,
        G: float = 1.0,
        screening_mode: Optional[str] = None,
        order: Optional[int] = None,
    ) -> None:
        if isinstance(masses_or_G, (np.ndarray, list, tuple)):
            eff_G = float(G)
            self.masses: Optional[np.ndarray] = np.asarray(masses_or_G, dtype=np.float64)
            self.N: int = len(self.masses)
            self.G_mass_matrix: Optional[np.ndarray] = eff_G * np.outer(self.masses, self.masses)
            np.fill_diagonal(self.G_mass_matrix, 0.0)
        else:
            eff_G = float(masses_or_G)
            self.masses = None
            self.N = 0
            self.G_mass_matrix = None
        super().__init__(G=eff_G, screening_mode=screening_mode, order=order)
        self._cached_tree_ref: Any = None
        self._cached_tree_id: Optional[int] = None
        self._cached_sn_key: Optional[Tuple[Tuple[int, ...], ...]] = None
        self._cached_non_leaf_idx_a: Optional[np.ndarray] = None
        self._cached_non_leaf_idx_b: Optional[np.ndarray] = None

    def invalidate_cache(self) -> None:
        """Invalidate topological pair caches upon atlas chart flip."""
        self._topology_cache.clear()
        self._cached_tree_ref = None
        self._cached_tree_id = None
        self._cached_sn_key = None
        self._cached_non_leaf_idx_a = None
        self._cached_non_leaf_idx_b = None

    def _get_non_leaf_pairs(
        self,
        tree: Any,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
        n_bodies: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Cache and return 1D index arrays (idx_a, idx_b) of non-leaf, non-intra-supernode pairs."""
        N = self.N if n_bodies is None else int(n_bodies)
        tree_id = id(tree)
        clusters = _normalize_supernode_clusters(supernode_bodies, supernode_clusters)
        sn_key = tuple(sorted(clusters))
        if (
            self._cached_tree_ref is tree
            and self._cached_tree_id == tree_id
            and self._cached_sn_key == sn_key
            and self._cached_non_leaf_idx_a is not None
        ):
            return self._cached_non_leaf_idx_a, self._cached_non_leaf_idx_b  # type: ignore

        leaf_pairs = set()
        if tree is not None:
            if hasattr(tree, "get_leaf_pair_indices"):
                for i, j in tree.get_leaf_pair_indices():
                    leaf_pairs.add((min(i, j), max(i, j)))
            elif hasattr(tree, "get_leaf_pairs"):
                for i, j in tree.get_leaf_pairs():
                    leaf_pairs.add((min(i, j), max(i, j)))
            elif isinstance(tree, (list, set, tuple)):
                for i, j in tree:
                    leaf_pairs.add((min(i, j), max(i, j)))

        cluster_sets = [set(c) for c in clusters]
        idx_a_list: List[int] = []
        idx_b_list: List[int] = []
        for i in range(N):
            for j in range(i + 1, N):
                if (i, j) in leaf_pairs:
                    continue
                if cluster_sets and any(i in c_set and j in c_set for c_set in cluster_sets):
                    continue
                idx_a_list.append(i)
                idx_b_list.append(j)

        self._cached_tree_ref = tree
        self._cached_tree_id = tree_id
        self._cached_sn_key = sn_key
        self._cached_non_leaf_idx_a = np.array(idx_a_list, dtype=np.intp)
        self._cached_non_leaf_idx_b = np.array(idx_b_list, dtype=np.intp)
        return self._cached_non_leaf_idx_a, self._cached_non_leaf_idx_b

    def compute_apex_forces(  # type: ignore[override]
        self,
        state_or_q: Union[ExtendedPhaseState, np.ndarray],
        tree: Any = None,
        supernode_bodies: Optional[Set[int]] = None,
        include_screened_monopole: bool = False,
        mask_leaf_pairs: Optional[bool] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> Union[Tuple[np.ndarray, float], np.ndarray]:
        """
        Dispatches to canonical super().compute_apex_forces when given an ExtendedPhaseState,
        or evaluates raw array forces when given q: np.ndarray.
        """
        if isinstance(state_or_q, ExtendedPhaseState):
            return super().compute_apex_forces(
                state_or_q,
                tree,
                supernode_bodies=supernode_bodies,
                include_screened_monopole=include_screened_monopole,
                mask_leaf_pairs=mask_leaf_pairs,
                supernode_clusters=supernode_clusters,
            )

        q = np.asarray(state_or_q, dtype=np.float64)
        N = q.shape[0]
        forces = np.zeros_like(q, dtype=np.float64)
        if self.G_mass_matrix is None:
            raise ValueError("Raw array call compute_apex_forces(q, ...) requires masses at initialization.")
        idx_a, idx_b = self._get_non_leaf_pairs(
            tree,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
            n_bodies=N,
        )
        n_pairs = len(idx_a)
        if n_pairs == 0:
            return forces

        g_mm = self.G_mass_matrix[idx_a, idx_b]
        if _HAS_NUMBA:
            forces, _ = _njit_accumulate_apex_pairs(
                q,
                np.asarray(idx_a, dtype=np.int64),
                np.asarray(idx_b, dtype=np.int64),
                np.asarray(g_mm, dtype=np.float64),
                N,
            )
            return forces

        dr = q[idx_b] - q[idx_a]
        r2 = np.einsum("ij,ij->i", dr, dr)
        valid = r2 > 0.0
        inv_r3 = np.zeros_like(r2)
        inv_r3[valid] = 1.0 / (r2[valid] * np.sqrt(r2[valid]))
        f_pair = (g_mm * inv_r3)[:, None] * dr
        np.add.at(forces, idx_a, f_pair)
        np.add.at(forces, idx_b, -f_pair)
        return forces


class OptimizedKeplerSolver(_CanonicalOptimizedKeplerSolver):
    """
    Inlined Stumpff-Sundman universal leaf propagator subclassing the canonical
    OptimizedKeplerSolver. Exposes both:
    - OptimizedKeplerSolver.stumpff_c0123_scalar(z) (degree-6 (6-stage) Horner + Stumpff recurrence)
    - solver.step(r0, w0, mu, gm, dt) / solver.step_kepler(r0, w0, mu, GM, dt_target)
    - solver.step_drift_and_kepler(state, tree, dt, ...) for ContextVar injection.
    """

    def __init__(self, G: float = 1.0, max_iter: int = 50, tol: Optional[float] = None) -> None:
        super().__init__(G=float(G), max_iter=max_iter, tol=tol)

    @staticmethod
    def stumpff_c0123_scalar(z: float) -> Tuple[float, float, float, float]:
        """Delegates directly to canonical solvers.optimized_kepler_sundman.stumpff_c0123_scalar."""
        return stumpff_c0123_scalar(float(z))

    def step(
        self,
        r0: np.ndarray,
        w0: np.ndarray,
        mu: float,
        gm: float,
        dt: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Delegates directly to canonical step_kepler_sundman_fast."""
        return step_kepler_sundman_fast(
            r0=r0,
            w0=w0,
            mu=float(mu),
            GM=float(gm),
            dt_target=float(dt),
            abs_tol=self.tol,
            max_iter=self.max_iter,
        )


OptimizedKeplerSundmanSolver = OptimizedKeplerSolver


__all__ = [
    "_CURRENT_APEX_EVALUATOR",
    "_CURRENT_KEPLER_SOLVER",
    "_CURRENT_SUPERNODE_SOLVER",
    "_CURRENT_SUPERNODE_N_SUBSTEPS",
    "set_current_apex_evaluator",
    "get_current_apex_evaluator",
    "reset_current_apex_evaluator",
    "set_current_kepler_solver",
    "get_current_kepler_solver",
    "reset_current_kepler_solver",
    "set_current_supernode_solver",
    "get_current_supernode_solver",
    "reset_current_supernode_solver",
    "set_current_supernode_n_substeps",
    "get_current_supernode_n_substeps",
    "reset_current_supernode_n_substeps",
    "OptimizedScreenedApexPotential",
    "OptimizedKeplerSolver",
    "OptimizedKeplerSundmanSolver",
    "stumpff_c0123_scalar",
    "step_kepler_sundman_fast",
    "OptimizedSupernodeSolver",
    "step_supernode_synchronous_fast",
    "JitSupernodeSolver",
    "step_supernode_synchronous_jit",
]


