"""
Discrete Dyadic Multi-Rate Scheduler and Palindromic Time Lattice.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.1 & Section 5.5 (Multi-Rate Symplectic Scheduling and Dyadic Lattices),
  Equations (5.106)-(5.132).
- Theoretical Foundations:
  1. Ge--Marsden Theorem (1988/1993): Arbitrary continuous step-size adaptation dt(t) in
     geometric numerical integration destroys the existence of an asymptotic shadow Hamiltonian H_tilde,
     introducing secular energy drift.
  2. Dyadic Lattice Invariance: Subdividing the global macro-step dt_macro into discrete dyadic powers
     dt_k = dt_macro / 2^k (with k in N_0) preserves exact symplecticity, with block transitions
     occurring exclusively at macro-block boundaries t_n = n * dt_macro.
  3. Palindromic Symmetry Preservation: Inside each dyadic macro-block, sub-steps are executed
     symmetrically to preserve microscopic time reversibility Phi(-h) o Phi(h) = Id.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Callable, Set, Union
import math
import numpy as np
import mpmath as mp

from ..core.constants import DEFAULT_MAX_DYADIC_LEVEL, to_f128
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree


@dataclass(frozen=True)
class ColoredInterval:
    """
    Time interval colored with a dyadic refinement level k in N_0.

    Attributes:
        t_start: Initial time bound.
        t_end: Final time bound.
        level: Dyadic refinement exponent k >= 0 (substep count = 2^k).
        description: Informative label (e.g. 'Binary periastron passage').
    """
    t_start: float
    t_end: float
    level: int
    description: str = ""
    t_min: float = field(init=False)
    t_max: float = field(init=False)

    def __post_init__(self) -> None:
        if self.level < 0:
            raise ValueError(f"Dyadic level must be non-negative, got {self.level}")
        t_min = min(float(self.t_start), float(self.t_end))
        t_max = max(float(self.t_start), float(self.t_end))
        object.__setattr__(self, "t_min", t_min)
        object.__setattr__(self, "t_max", t_max)

    def contains(self, t: float) -> bool:
        """Returns True if t is contained in the closed interval [t_min, t_max]."""
        return self.t_min <= t <= self.t_max

    def intersects(self, t1: float, t2: float) -> bool:
        """Returns True if [min(t1, t2), max(t1, t2)] overlaps this colored window."""
        b_min = min(t1, t2)
        b_max = max(t1, t2)
        return not (self.t_max < b_min or self.t_min > b_max)


class DyadicScheduler:
    """
    Symplectic multi-rate time-step scheduler operating on a discrete dyadic lattice:
        dt_k = dt_macro / 2^k.

    Features:
    1. Static / Pre-scheduled interval coloring (e.g. close encounters, transit events).
    2. Dynamic timescale adaptation based on local cluster velocities and curvatures.
    3. Exact dyadic sub-step subdivision satisfying sum dt_substep == dt_macro to machine precision.
    4. Palindromic sub-step ordering preserving microscopic time reversibility.
    """

    def __init__(
        self,
        dt_macro: float,
        default_level: int = 0,
        max_level: int = DEFAULT_MAX_DYADIC_LEVEL,
        G: float = 1.0,
        warn_on_saturation: bool = False,
    ) -> None:
        if dt_macro == 0.0:
            raise ValueError("dt_macro cannot be zero.")
        if default_level < 0:
            raise ValueError(f"default_level must be >= 0, got {default_level}")
        if max_level < default_level:
            raise ValueError(f"max_level ({max_level}) must be >= default_level ({default_level})")

        self.dt_macro = float(dt_macro)
        self.default_level = int(default_level)
        self.min_level = int(default_level)
        self.max_level = int(max_level)
        self.G = float(G)
        self.warn_on_saturation = bool(warn_on_saturation)
        self.last_ideal_level: int = int(default_level)
        self.saturation_events: List[dict] = []
        self.intervals: List[ColoredInterval] = []

    def color_interval(
        self,
        t_start: float,
        t_end: float,
        level: int,
        description: str = "",
    ) -> DyadicScheduler:
        """Registers a colored temporal window with dyadic refinement level k."""
        if level < 0:
            raise ValueError(f"Dyadic level must be >= 0, got {level}")
        clamped_level = min(level, self.max_level)
        self.intervals.append(
            ColoredInterval(
                t_start=t_start,
                t_end=t_end,
                level=clamped_level,
                description=description,
            )
        )
        return self

    def compute_dynamical_level(
        self,
        state: ExtendedPhaseState | np.ndarray,
        p: Optional[np.ndarray] = None,
        masses: Optional[np.ndarray] = None,
        leaf_pairs: Optional[List[Tuple[int, int]]] = None,
        tree: Optional[JacobiTree] = None,
        dt_macro: Optional[float] = None,
        safety_factor: float = 0.05,
        G: Optional[float] = None,
    ) -> int:
        r"""
        Dynamically derives the optimal dyadic refinement exponent k from local
        orbital velocities, mutual separations across all pairs, and proximity thresholds:
            \tau_{\min} = \min_{i < j} (r_{ij} / ||v_i - v_j||),
            k_{\tau} = \lceil\log_2(\Delta t_{\mathrm{macro}} / (f_{\mathrm{safety}} \cdot \tau_{\min}))\rceil,
            k_{\mathrm{prox}} = \text{conormal proximity lattice level from } r_{\min}.

        Supports both the canonical (state, tree=...) signature and raw arrays
        (q, p, masses, leaf_pairs=...).
        """
        if isinstance(p, JacobiTree) or (p is not None and hasattr(p, "leaf_pairs")):
            tree = p  # type: ignore[assignment]
            p = None

        dt = float(self.dt_macro if dt_macro is None else dt_macro)
        eff_G = self.G if G is None else float(G)
        min_lvl = max(self.default_level, int(getattr(self, "min_level", 0)))

        is_f128 = False
        if hasattr(state, "q") and hasattr(state, "p") and hasattr(state, "masses"):
            is_f128 = bool(getattr(state, "is_f128", False))
            N = int(state.N)
            if N < 2:
                return min_lvl
            if is_f128:
                q_flt = np.array([[float(state.q[i, d]) for d in range(3)] for i in range(N)], dtype=np.float64)
                p_flt = np.array([[float(state.p[i, d]) for d in range(3)] for i in range(N)], dtype=np.float64)
                m_flt = np.array([float(state.masses[i]) for i in range(N)], dtype=np.float64)
            else:
                q_flt = np.asarray(state.q, dtype=np.float64)
                p_flt = np.asarray(state.p, dtype=np.float64)
                m_flt = np.asarray(state.masses, dtype=np.float64)
        else:
            q_flt = np.asarray(state, dtype=np.float64)
            p_flt = np.asarray(p, dtype=np.float64)
            m_flt = np.asarray(masses, dtype=np.float64)
            N = int(len(m_flt))
            if N < 2:
                return min_lvl

        leaf_pair_set: Set[Tuple[int, int]] = set()
        if tree is not None and hasattr(tree, "leaf_pairs"):
            for node in tree.leaf_pairs:
                if node.body_indices and len(node.body_indices) == 2:
                    i_a, i_b = sorted(node.body_indices)
                    leaf_pair_set.add((i_a, i_b))
        if leaf_pairs is not None:
            for pair in leaf_pairs:
                if len(pair) == 2:
                    leaf_pair_set.add((min(pair[0], pair[1]), max(pair[0], pair[1])))

        min_tau = float("inf")
        limiting_pair = (0, 1)
        v_rel_floor = 1e-30 if is_f128 else 1e-14
        gm_floor = 1e-100 if is_f128 else 1e-50

        for i in range(N):
            v_i = p_flt[i] / m_flt[i]
            for j in range(i + 1, N):
                r_vec = q_flt[i] - q_flt[j]
                r_ij = float(np.linalg.norm(r_vec))
                if r_ij <= 0.0:
                    continue

                v_vec = v_i - (p_flt[j] / m_flt[j])
                v_rel = float(np.linalg.norm(v_vec))
                GM = eff_G * (m_flt[i] + m_flt[j])

                # 1. Instantaneous kinematic and dynamical free-fall timescales
                tau_kin = r_ij / max(v_rel, v_rel_floor)
                tau_dyn = math.sqrt((r_ij ** 3) / max(GM, gm_floor))
                tau_loc = min(tau_kin, tau_dyn)
                r_loc = r_ij

                # Two-body Keplerian energy and semi-major axis
                E_2b = 0.5 * (v_rel ** 2) - GM / r_ij
                is_leaf_pair = (i, j) in leaf_pair_set

                # 2. Definition 5.4 & Definition 6.2 Item 3 (V4): Universal Keplerian Pericenter Timescale
                if v_rel > v_rel_floor and not (is_leaf_pair and N == 2 and E_2b < -1e-12):
                    L_vec = np.cross(r_vec, v_vec)
                    L = float(np.linalg.norm(L_vec))
                    e_sq = 1.0 + 2.0 * E_2b * (L ** 2) / max(GM ** 2, 1e-100)
                    if e_sq >= 0.0:
                        e = math.sqrt(e_sq)
                        r_peri = (L ** 2) / max(GM * (1.0 + e), gm_floor)
                        if 0.0 < r_peri < r_ij:
                            rv = abs(float(np.dot(r_vec, v_vec)))
                            energy_dimless = 2.0 * abs(E_2b) * r_ij / max(GM, gm_floor)
                            if energy_dimless < 1e-4 or e < 1e-12:
                                t_kepler = (
                                    math.sqrt(2.0 * max(0.0, r_ij - r_peri))
                                    / math.sqrt(max(GM, gm_floor))
                                ) * ((r_ij + 2.0 * r_peri) / 3.0)
                            elif E_2b < 0.0:
                                a_orb = -GM / (2.0 * E_2b)
                                n_mean = math.sqrt(max(GM, gm_floor) / max(a_orb ** 3, 1e-100))
                                e_cos_E = 1.0 - r_ij / max(a_orb, gm_floor)
                                e_sin_E = rv / math.sqrt(max(GM * a_orb, 1e-100))
                                E_anom = math.atan2(e_sin_E, e_cos_E)
                                M_anom = max(0.0, E_anom - e * math.sin(E_anom))
                                t_kepler = M_anom / max(n_mean, gm_floor)
                            else:
                                a_hyp = GM / (2.0 * E_2b)
                                n_hyp = math.sqrt(max(GM, gm_floor) / max(a_hyp ** 3, 1e-100))
                                e_sinh_H = rv / math.sqrt(max(GM * a_hyp, 1e-100))
                                H_anom = math.asinh(e_sinh_H / max(e, 1.0))
                                M_hyp = max(0.0, e * math.sinh(H_anom) - H_anom)
                                t_kepler = M_hyp / max(n_hyp, gm_floor)

                            if t_kepler <= 2.0 * abs(dt):
                                dt_rem = max(0.0, t_kepler - abs(dt))
                                r_eff = min(
                                    r_ij,
                                    (r_peri ** 1.5 + 1.5 * math.sqrt(2.0 * max(GM, gm_floor)) * dt_rem) ** (2.0 / 3.0),
                                )
                                v_eff = math.sqrt(max(1e-14, v_rel ** 2 + 2.0 * GM * (1.0 / max(r_eff, gm_floor) - 1.0 / r_ij)))
                                tau_peri_kin = r_eff / v_eff
                                tau_peri_dyn = math.sqrt((r_eff ** 3) / max(GM, gm_floor))
                                tau_loc = min(tau_loc, tau_peri_kin, tau_peri_dyn)
                                r_loc = min(r_ij, r_eff)

                # 3. Theorem 5.20 & Theorem 6.6 (V5a & V5b): Two-Regime Leaf/Apex Decoupling and Hill Isolation Guard
                tau_pair = tau_loc
                if is_leaf_pair and E_2b < -1e-12:
                    M_C = float(m_flt[i] + m_flt[j])
                    R_cm = (m_flt[i] * q_flt[i] + m_flt[j] * q_flt[j]) / M_C
                    a_C = -GM / (2.0 * E_2b)
                    tau_dyn_inv = math.sqrt((a_C ** 3) / max(GM, gm_floor))
                    if N == 2:
                        tau_pair = tau_dyn_inv
                    else:
                        is_isolated_leaf = True
                        min_r_hill_third = float("inf")
                        for k_ext in range(N):
                            if k_ext == i or k_ext == j:
                                continue
                            R_ext = float(np.linalg.norm(q_flt[k_ext] - R_cm))
                            m_ext = float(m_flt[k_ext])
                            R_hill = R_ext * ((M_C / max(3.0 * (M_C + m_ext), gm_floor)) ** (1.0 / 3.0))
                            r_hill_third = R_hill / 3.0
                            if r_hill_third < min_r_hill_third:
                                min_r_hill_third = r_hill_third
                            if R_ext < 3.0 * r_ij or r_ij > r_hill_third:
                                is_isolated_leaf = False
                                break

                        if is_isolated_leaf:
                            a_eff = min(a_C, min_r_hill_third)
                            char_size = max(r_ij, a_eff)
                            tau_dyn_eff = math.sqrt((a_eff ** 3) / max(GM, gm_floor))
                            eps_tidal = 0.0
                            for k_ext in range(N):
                                if k_ext == i or k_ext == j:
                                    continue
                                R_ext = float(np.linalg.norm(q_flt[k_ext] - R_cm))
                                m_ext = float(m_flt[k_ext])
                                eps_k = (m_ext / max(M_C, gm_floor)) * ((char_size / max(R_ext, 1e-30)) ** 3)
                                if eps_k > eps_tidal:
                                    eps_tidal = eps_k

                            if eps_tidal <= 1e-30:
                                tau_pair = tau_dyn_inv
                            else:
                                quad_factor = (min(r_loc, char_size) / max(char_size, gm_floor)) ** 2
                                denom_reg1 = max((eps_tidal * quad_factor) ** 0.25, 1e-30)
                                tau_reg1 = tau_loc / denom_reg1
                                tau_reg2 = tau_dyn_eff * (0.1 * (safety_factor ** 3)) / max(eps_tidal ** 0.75, 1e-30)
                                tau_pair = min(tau_dyn_inv, max(tau_reg1, tau_reg2))

                if tau_pair < min_tau:
                    min_tau = tau_pair
                    limiting_pair = (i, j)

        # Continuous scale-invariant dyadic refinement level:
        # k = ceil(log2(Delta t_macro / (eta * tau_enc)))
        dynamic_k = 0
        if not math.isinf(min_tau) and min_tau > 0.0:
            target_dt = max(min_tau * safety_factor, 1e-15)
            ratio = abs(dt) / target_dt
            if ratio > 1.0:
                dynamic_k = math.ceil(math.log2(ratio))

        ideal_k = max(dynamic_k, min_lvl)
        self.last_ideal_level = ideal_k
        if ideal_k > self.max_level:
            t_val = float(state.t) if isinstance(state, ExtendedPhaseState) else float("nan")
            event = {
                "t": t_val,
                "ideal_k": ideal_k,
                "clamped_k": self.max_level,
                "min_tau": min_tau,
                "limiting_pair": limiting_pair,
            }
            self.saturation_events.append(event)
            if self.warn_on_saturation:
                print(
                    f"[DyadicScheduler SATURATION] t={t_val:.4f}: ideal k={ideal_k} clamped to "
                    f"max_level={self.max_level} (min_tau={min_tau:.4e}, pair={limiting_pair})"
                )

        return min(ideal_k, self.max_level)

    def get_level_for_block(
        self,
        t_current: float,
        dt_block: Optional[float] = None,
        state: Optional[ExtendedPhaseState] = None,
        tree: Optional[JacobiTree] = None,
        state_end: Optional[ExtendedPhaseState] = None,
        symmetric: bool = False,
        q: Optional[np.ndarray] = None,
        p: Optional[np.ndarray] = None,
        masses: Optional[np.ndarray] = None,
        leaf_pairs: Optional[List[Tuple[int, int]]] = None,
    ) -> int:
        """
        Determines the dyadic refinement level k for macro-block [t_current, t_current + dt_block],
        taking the maximum between pre-scheduled colored windows and dynamic timescale estimates.
        When state_end is provided, evaluates the symmetric block level:
            k_sym = max(k(z_start), k(z_end)).
        """
        dt = self.dt_macro if dt_block is None else float(dt_block)
        min_lvl = max(self.default_level, int(getattr(self, "min_level", 0)))

        # 1. Check pre-scheduled intervals
        matching_levels = [
            inter.level for inter in self.intervals
            if inter.intersects(t_current, t_current + dt)
        ]
        scheduled_k = max(matching_levels) if matching_levels else min_lvl

        # 2. Check dynamic state if provided
        dynamic_k = 0
        if state is not None:
            dynamic_k = self.compute_dynamical_level(
                state, tree=tree, leaf_pairs=leaf_pairs, dt_macro=dt
            )
            if state_end is not None:
                k_end = self.compute_dynamical_level(
                    state_end, tree=tree, leaf_pairs=leaf_pairs, dt_macro=dt
                )
                dynamic_k = max(dynamic_k, k_end)
        elif q is not None and p is not None and masses is not None:
            dynamic_k = self.compute_dynamical_level(
                q, p=p, masses=masses, leaf_pairs=leaf_pairs, tree=tree, dt_macro=dt
            )

        return min(max(scheduled_k, dynamic_k, min_lvl), self.max_level)

    def get_substep_params(
        self,
        t_current: Union[float, mp.mpf],
        q: Optional[np.ndarray] = None,
        p: Optional[np.ndarray] = None,
        masses: Optional[np.ndarray] = None,
        leaf_pairs: Optional[List[Tuple[int, int]]] = None,
        dt_block: Optional[Union[float, mp.mpf]] = None,
        state: Optional[ExtendedPhaseState] = None,
        tree: Optional[JacobiTree] = None,
        state_end: Optional[ExtendedPhaseState] = None,
        symmetric: bool = False,
    ) -> Tuple[int, Union[float, mp.mpf], int]:
        """
        Computes exact subdivision parameters for the macro-block:
            num_substeps = 2^k,
            dt_substep = dt_block / 2^k.

        Exact identity:
            num_substeps * dt_substep == dt_block (to float machine precision).
        """
        if isinstance(q, (int, float, mp.mpf)) and dt_block is None:
            dt_block = q
            q = None

        dt = self.dt_macro if dt_block is None else dt_block
        k = self.get_level_for_block(
            float(t_current),
            dt_block=float(dt),
            state=state,
            tree=tree,
            state_end=state_end,
            symmetric=symmetric,
            q=q,
            p=p,
            masses=masses,
            leaf_pairs=leaf_pairs,
        )
        num_substeps = 1 << k
        if (state is not None and state.is_f128) or isinstance(dt, mp.mpf):
            dt_substep = to_f128(dt) / num_substeps
        else:
            dt_substep = float(dt) / float(num_substeps)
        return num_substeps, dt_substep, k

    @staticmethod
    def generate_bit_exact_timestamps(
        t_start: Union[float, mp.mpf],
        dt_block: Union[float, mp.mpf],
        k: int,
        is_f128: bool = False,
    ) -> np.ndarray:
        """
        Generates bit-exact intermediate timestamps:
            t^{(m)} = t_start + m * (dt_block / 2^k),  m in {0, ..., 2^k}.

        Eliminates recursive accumulation noise sum(dt_substep) != dt_block.
        """
        num_substeps = 1 << k
        if is_f128 or isinstance(t_start, mp.mpf) or isinstance(dt_block, mp.mpf):
            t_s = to_f128(t_start)
            dt_b = to_f128(dt_block)
            step_frac = dt_b / num_substeps
            return np.array([t_s + m * step_frac for m in range(num_substeps + 1)], dtype=object)
        step_fraction = float(dt_block) / float(num_substeps)
        indices = np.arange(num_substeps + 1, dtype=np.float64)
        return float(t_start) + indices * step_fraction

    def clear_intervals(self) -> None:
        """Clears all colored interval reservations."""
        self.intervals.clear()
