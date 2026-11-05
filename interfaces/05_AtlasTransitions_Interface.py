"""
Standalone interface for Section 6.5 (Appendix A.26): Atlas Transitions,
Gauge Connection, and Symplectic Post-Processing (Phases 4 & 5).

Reference Implementation Interface:
Re-exports and delegates directly to the canonical reference classes:
1. SL2TransitionMatrix (including .from_masses and .from_trees with intrinsic
   node orientation signs closing the 3-chart Cech 1-cocycle to < 1e-15),
   build_tree_transition_symplectic_12, verify_wronskian_unimodularity
   from core.transitions.
2. AtlasManager (gamma_monitor=0.60, gamma_build=0.70 two-tier hysteresis and
   dwell-time guard tau_min) from system.atlas_manager.
3. PBM4GaugeTransitionOperator and apply_gauge_connection: true non-commutative
   symplectic sandwich T_{(A,h_A)->(B,h_B)}^{(PBM4)} = (Z_4^{BM4,(B)}(h_B))^{-1} o M_{A->B} o Z_4^{BM4,(A)}(h_A)
   with default CorrectorC4(mode="homological") from integrators.processed.
4. LieBargmannInvariants, compute_invariants, compute_lie_bargmann_invariants
   from core.invariants.
"""

from __future__ import annotations
import sys
from pathlib import Path
from typing import Tuple, List, Optional, Union
import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import JacobiTree
from normal_cone_nbody_float64.core.transitions import (
    SL2TransitionMatrix,
    build_tree_transition_symplectic_12,
)


def verify_wronskian_unimodularity(f: float, g: float, f_dot: float, g_dot: float) -> float:
    r"""
    Computes absolute deviation from unimodularity:
        |f * g_dot - f_dot * g - 1.0|.
    """
    return float(abs(f * g_dot - f_dot * g - 1.0))
from normal_cone_nbody_float64.system.atlas_manager import (
    AtlasManager as _CanonicalAtlasManager,
)
from normal_cone_nbody_float64.integrators.processed import (
    PBM4GaugeTransitionOperator,
    apply_gauge_connection,
)
from normal_cone_nbody_float64.core.invariants import (
    LieBargmannInvariants,
    compute_invariants,
    compute_truncation_residual,
)


def compute_dynamic_dwell_time(q: np.ndarray, p: np.ndarray, masses: np.ndarray) -> float:
    r"""
    Computes the theoretical minimal dwell time tau_min:
        tau_min = (16 / 15) * (r_min / v_max).
    Delegates directly to _CanonicalAtlasManager.compute_dynamic_dwell_time.
    """
    state = ExtendedPhaseState(
        q=np.asarray(q, dtype=np.float64),
        p=np.asarray(p, dtype=np.float64),
        masses=np.asarray(masses, dtype=np.float64),
    )
    return _CanonicalAtlasManager.compute_dynamic_dwell_time(state)


def check_tree_admissibility(
    q: np.ndarray,
    leaf_pairs: List[Tuple[int, int]],
    gamma: float = 1.0,
    masses: Optional[np.ndarray] = None,
) -> Optional[Tuple[Tuple[int, int], int, float, float]]:
    r"""
    Evaluates unilateral Hill-weighted Delaunay--Jacobi admissibility (Phase 4.1).
    Returns None if all pairs satisfy d_ext^{Hill} >= gamma * r_ij,
    or the worst violation ((i, j), k_intruder, d_ext, r_ij).
    """
    N = len(q)
    if N <= 2:
        return None

    worst_violation = None
    min_ratio = float("inf")

    for i, j in leaf_pairs:
        r_ij = float(np.linalg.norm(q[i] - q[j]))
        if r_ij <= 0.0:
            continue

        d_min = float("inf")
        k_intruder = -1
        m_ij = float(masses[i] + masses[j]) if masses is not None else 0.0
        is_hier = (
            (masses is not None)
            and (min(float(masses[i]), float(masses[j])) < 0.01 * max(float(masses[i]), float(masses[j])))
        )

        for k in range(N):
            if k == i or k == j:
                continue
            d_ik = float(np.linalg.norm(q[i] - q[k]))
            d_jk = float(np.linalg.norm(q[j] - q[k]))
            d_raw = min(d_ik, d_jk)
            if masses is not None and m_ij > 0.0:
                m_k = max(float(masses[k]), 1e-30)
                raw_hill = (m_ij / (2.0 * m_k)) ** (1.0 / 3.0)
                if is_hier:
                    hill_scale = raw_hill
                elif gamma >= 1.0 - 1e-15:
                    hill_scale = 1.0
                elif gamma > 0.60:
                    w_hyst = (1.0 - gamma) / 0.40
                    hill_scale = 1.0 + w_hyst * (max(1.0, raw_hill) - 1.0)
                else:
                    hill_scale = max(1.0, raw_hill)
            else:
                hill_scale = 1.0
            d_cand = d_raw * hill_scale
            if d_cand < d_min:
                d_min = d_cand
                k_intruder = k

        ratio = d_min / r_ij
        if ratio < gamma and ratio < min_ratio:
            min_ratio = ratio
            worst_violation = ((i, j), k_intruder, d_min, r_ij)

    return worst_violation


class AtlasManager(_CanonicalAtlasManager):
    r"""
    Canonical dynamical Jacobi atlas manager subclassing _CanonicalAtlasManager
    (gamma = 0.60, gamma_build = 0.70, dwell-time guard tau_min), supporting both:
    - Canonical signature: evaluate_tree_adaptation(state, current_tree) -> (JacobiTree, bool)
    - Raw array signature: evaluate_tree_adaptation(q, p, masses, current_leaf_pairs, t=0.0)
      -> (List[Tuple[int, int]], bool)
    """

    @property
    def last_transition_time(self) -> Optional[float]:
        if "_last_transition_time" in self.__dict__:
            return self.__dict__["_last_transition_time"]
        return self.transition_history[-1].t if self.transition_history else None

    @last_transition_time.setter
    def last_transition_time(self, value: Optional[float]) -> None:
        self.__dict__["_last_transition_time"] = value

    def build_optimal_leaf_pairs(
        self,
        q: np.ndarray,
        masses: np.ndarray,
    ) -> List[Tuple[int, int]]:
        """Builds the optimal Jacobi tree via super().build_optimal_tree and returns leaf pair indices."""
        state = ExtendedPhaseState(
            q=np.asarray(q, dtype=np.float64),
            p=np.zeros_like(q, dtype=np.float64),
            masses=np.asarray(masses, dtype=np.float64),
        )
        tree = self.build_optimal_tree(state)
        return tree.get_leaf_pair_indices()

    def evaluate_tree_adaptation(  # type: ignore[override]
        self,
        state_or_q: Union[ExtendedPhaseState, np.ndarray],
        tree_or_p: Union[JacobiTree, np.ndarray],
        masses: Optional[np.ndarray] = None,
        current_leaf_pairs: Optional[List[Tuple[int, int]]] = None,
        t: float = 0.0,
    ) -> Union[Tuple[JacobiTree, bool], Tuple[List[Tuple[int, int]], bool]]:
        """Dispatches between canonical (state, current_tree) and raw (q, p, masses, leaf_pairs, t)."""
        if isinstance(state_or_q, ExtendedPhaseState) and isinstance(tree_or_p, JacobiTree):
            return super().evaluate_tree_adaptation(state_or_q, tree_or_p)

        q_arr = np.asarray(state_or_q, dtype=np.float64)
        p_arr = np.asarray(tree_or_p, dtype=np.float64)
        m_arr = np.asarray(masses, dtype=np.float64)
        leaves = list(current_leaf_pairs or [])

        eff_tau = self.tau_min
        if eff_tau == 0.0 and self.dynamic_tau_min and self.last_transition_time is not None:
            eff_tau = compute_dynamic_dwell_time(q_arr, p_arr, m_arr)
        if eff_tau > 0.0 and self.last_transition_time is not None:
            if (t - self.last_transition_time) < eff_tau - 1e-15:
                return leaves, False

        violation = check_tree_admissibility(q_arr, leaves, gamma=self.gamma, masses=m_arr)
        if violation is None:
            return leaves, False

        new_leaf_pairs = self.build_optimal_leaf_pairs(q_arr, m_arr)
        if set(new_leaf_pairs) == set(leaves):
            return leaves, False
        self.last_transition_time = float(t)
        return new_leaf_pairs, True


def compute_lie_bargmann_invariants(
    q_or_state: Union[np.ndarray, ExtendedPhaseState],
    p: Optional[np.ndarray] = None,
    masses: Optional[np.ndarray] = None,
    t: float = 0.0,
    G: float = 1.0,
    E0: float = 0.0,
) -> LieBargmannInvariants:
    r"""
    Computes the 10 canonical generators of the Galilean Lie--Bargmann group
    by direct delegation to core.invariants.compute_invariants.
    Supports both (state, G=1.0) and (q, p, masses, t=0.0, G=1.0, E0=0.0).
    """
    if isinstance(q_or_state, ExtendedPhaseState):
        return compute_invariants(q_or_state, G=float(G))
    state = ExtendedPhaseState(
        q=np.asarray(q_or_state, dtype=np.float64),
        p=np.asarray(p, dtype=np.float64),
        masses=np.asarray(masses, dtype=np.float64),
        t=float(t),
        E=float(E0),
        E0=float(E0),
    )
    return compute_invariants(state, G=float(G))


__all__ = [
    "SL2TransitionMatrix",
    "build_tree_transition_symplectic_12",
    "verify_wronskian_unimodularity",
    "compute_dynamic_dwell_time",
    "check_tree_admissibility",
    "AtlasManager",
    "PBM4GaugeTransitionOperator",
    "apply_gauge_connection",
    "LieBargmannInvariants",
    "compute_invariants",
    "compute_lie_bargmann_invariants",
    "compute_truncation_residual",
]

