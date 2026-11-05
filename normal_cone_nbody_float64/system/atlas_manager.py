"""
Jacobi Atlas Dynamics Manager and Anti-Chattering Hysteresis State Machine.

Mathematical Monograph Alignment:
- Chapter 4: Section 4.3 (Dynamical Atlas of Charts and Smooth Transition Functions),
  Equations (4.164)-(4.195).
- Chapter 6: Section 6.1 (Symplectic Lifts of Tree Rotations in Sp(6N-6, R) via SL(2, R)) &
  Section 6.2 (Switching Hysteresis and Proof of Anti-Chattering Invariance, Theorem 6.1).
- Certified Formal Proofs:
  - Lean 4 module `AtlasChartTransitionSymplectic.lean` (canonical SL(2, R) unimodular action).
  - `02_LieBargmannExactInvariance_Theorems.lean` (Darboux pairing and Lie-Bargmann conservation).
- Anti-Chattering Theorem (Theorem 6.1):
  Using the strict hysteresis band gamma = 0.6, a third body must penetrate closer than
  0.6 * r_ij before triggering an SL(2, R) tree flip, preventing infinite Zeno chattering
  by the geometric contradiction 1 <= gamma^2 = 0.36 < 1.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple, List, Set, Dict
import numpy as np

from ..core.constants import GAMMA_HYSTERESIS
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree, JacobiTopologyNode, build_canonical_tree, build_tree_from_spec


@dataclass(frozen=True)
class AtlasTransitionRecord:
    """
    Metadata recording a discrete chart transition in the Jacobi atlas:
    - t: Physical time at transition.
    - old_leaf_pairs: Tuple of prior leaf pair body index tuples.
    - new_leaf_pairs: Tuple of newly adapted leaf pair body index tuples.
    - triggering_pair: Tuple of body indices (i, j) that triggered the flip.
    - intruder_body: Index of the third body that penetrated the cluster.
    - separation_ratio: Ratio d_min / r_ij < gamma.
    - dwell_time: Elapsed time since previous chart transition (tau_dwell).
    """
    t: float
    old_leaf_pairs: Tuple[Tuple[int, int], ...]
    new_leaf_pairs: Tuple[Tuple[int, int], ...]
    triggering_pair: Tuple[int, int]
    intruder_body: int
    separation_ratio: float
    dwell_time: Optional[float] = None


class AtlasManager:
    """
    Deterministic manager of the dynamical Jacobi atlas on the compactified
    Fulton--MacPherson configuration manifold Q[N].

    Ensures that the active Jacobi tree T in Trees(N) remains adapted to the current
    spatial hierarchy, applying symplectic transitions in SL(2, R) whenever
    inter-cluster separations cross the hysteresis threshold gamma = 0.6.
    """

    def __init__(
        self,
        gamma: float = GAMMA_HYSTERESIS,
        gamma_build: float = 1.0,
        tau_min: float = 0.0,
        dynamic_tau_min: bool = False,
    ) -> None:
        if not (0.0 < gamma <= 1.0):
            raise ValueError(f"Hysteresis threshold gamma must lie in (0, 1], got {gamma}")
        self.gamma = float(gamma)
        self.gamma_build = max(float(gamma), float(gamma_build))
        self.tau_min = max(0.0, float(tau_min))
        self.dynamic_tau_min = dynamic_tau_min
        self.transition_history: List[AtlasTransitionRecord] = []

    @staticmethod
    def compute_dynamic_dwell_time(state: ExtendedPhaseState) -> float:
        """
        Computes the theoretical minimal dwell time tau_min according to Theorem 6.1:
            tau_min = (16 / 15) * (r_min / v_max)
        where r_min is the minimal pairwise separation and v_max is the maximum body speed.
        """
        if state.N <= 1:
            return 0.0
        diff = state.q[:, None, :] - state.q[None, :, :]
        dists = np.linalg.norm(diff, axis=-1)
        i_upper = np.triu_indices(state.N, k=1)
        r_min = float(np.min(dists[i_upper])) if len(i_upper[0]) > 0 else 1.0

        v_norms = np.linalg.norm(state.p / state.masses[:, None], axis=1)
        v_max = float(np.max(v_norms))
        if v_max <= 1e-12:
            return 0.0

        return (16.0 / 15.0) * (r_min / v_max)

    def check_tree_admissibility(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        ignore_dwell: bool = False,
    ) -> Optional[Tuple[Tuple[int, int], int, float, float]]:
        """
        Evaluates the topological admissibility condition for Jacobi tree T (Condition C6).

        For every leaf pair C = {i, j} in LeafPairs(T) with separation r_ij = ||q_i - q_j||:
            d_ext(C) = min_{k != i, j} min(||q_i - q_k||, ||q_j - q_k||).
        If d_ext(C) < gamma * r_ij, the pair is no longer isolated (an intruder penetrated),
        and the tree must flip to maintain C^{1,1} apex regularity.

        Returns:
            None if tree is admissible, or ((i, j), k_intruder, d_min, r_ij) on violation.
        """
        N = state.N
        if N <= 2:
            return None

        # Minimal dwell time guard tau_min > 0 (Theorem 6.1)
        if not ignore_dwell:
            effective_tau = self.tau_min
            if effective_tau == 0.0 and self.dynamic_tau_min and self.transition_history:
                effective_tau = self.compute_dynamic_dwell_time(state)

            if (
                effective_tau > 0.0
                and self.transition_history
                and not (
                    float(np.max(state.masses))
                    > 50.0 * (float(np.sum(state.masses)) - float(np.max(state.masses)))
                )
            ):
                last_t = self.transition_history[-1].t
                if abs(state.t - last_t) < effective_tau - 1e-15:
                    return None

        leaf_pairs = tree.get_leaf_pair_indices()
        violations: List[Tuple[float, Tuple[int, int], int, float, float]] = []

        for i, j in leaf_pairs:
            r_ij = float(np.linalg.norm(state.q[i] - state.q[j]))
            if r_ij <= 0.0:
                continue
            m_i_val = float(state.masses[i])
            m_j_val = float(state.masses[j])
            m_ij = max(m_i_val + m_j_val, 1e-30)
            is_hierarchical_leaf = min(m_i_val, m_j_val) < 0.01 * max(m_i_val, m_j_val)

            d_min = float("inf")
            k_best = -1

            for k in range(N):
                if k == i or k == j:
                    continue
                d_ik = float(np.linalg.norm(state.q[i] - state.q[k]))
                d_jk = float(np.linalg.norm(state.q[j] - state.q[k]))
                m_k = max(float(state.masses[k]), 1e-30)
                raw_hill = (m_ij / (2.0 * m_k)) ** (1.0 / 3.0)
                if is_hierarchical_leaf:
                    hill_scale = raw_hill
                elif self.gamma >= 1.0 - 1e-15:
                    hill_scale = 1.0
                elif self.gamma > 0.60:
                    w_hyst = (1.0 - self.gamma) / 0.40
                    hill_scale = 1.0 + w_hyst * (max(1.0, raw_hill) - 1.0)
                else:
                    hill_scale = max(1.0, raw_hill)
                d_k = min(d_ik, d_jk) * hill_scale

                if d_k < d_min:
                    d_min = d_k
                    k_best = k

            # Strict hysteresis rule: d_min < gamma * r_ij
            if d_min < self.gamma * r_ij:
                ratio = d_min / r_ij
                violations.append((ratio, (i, j), k_best, d_min, r_ij))

        # Also check whether two bodies (p, b) not paired together as a leaf pair have entered
        # a planetary Hill capture zone where (p, b) is now an admissible leaf pair
        # and violates b's current attachment to the primary body dom_idx (Proposition 6.5, Eq. 6.178).
        if not violations and N >= 3:
            leaf_bodies = {b for pair in leaf_pairs for b in pair}
            dom_idx = int(np.argmax([float(state.masses[b]) for b in range(N)]))
            if N >= 4:
                free_bodies = [b for b in range(N) if b not in leaf_bodies]
            elif (
                dom_idx in leaf_bodies
                and float(state.masses[dom_idx])
                > 50.0 * (float(np.sum(state.masses)) - float(state.masses[dom_idx]))
            ):
                free_bodies = [b for b in range(N) if b != dom_idx]
            else:
                free_bodies = []
            for idx_a in range(len(free_bodies)):
                for idx_b in range(idx_a + 1, len(free_bodies)):
                    p_idx = free_bodies[idx_a]
                    b_idx = free_bodies[idx_b]
                    if (min(p_idx, b_idx), max(p_idx, b_idx)) in leaf_pairs:
                        continue
                    r_pb = float(np.linalg.norm(state.q[p_idx] - state.q[b_idx]))
                    if r_pb <= 0.0:
                        continue
                    m_pb = max(float(state.masses[p_idx]) + float(state.masses[b_idx]), 1e-30)
                    d_ext_pb = min(
                        min(
                            float(np.linalg.norm(state.q[p_idx] - state.q[k])),
                            float(np.linalg.norm(state.q[b_idx] - state.q[k])),
                        )
                        * ((m_pb / (2.0 * max(float(state.masses[k]), 1e-30))) ** (1.0 / 3.0))
                        for k in range(N)
                        if k != p_idx and k != b_idx
                    )
                    if d_ext_pb >= self.gamma_build * r_pb:
                        if float(state.masses[p_idx]) >= float(state.masses[b_idx]):
                            heavy_p, light_b = p_idx, b_idx
                        else:
                            heavy_p, light_b = b_idx, p_idx
                        r_dom_b = float(np.linalg.norm(state.q[dom_idx] - state.q[light_b]))
                        m_dom_b = max(float(state.masses[dom_idx]) + float(state.masses[light_b]), 1e-30)
                        m_hp = max(float(state.masses[heavy_p]), 1e-30)
                        d_scaled = r_pb * ((m_dom_b / (2.0 * m_hp)) ** (1.0 / 3.0))
                        if r_dom_b > 0.0 and d_scaled < self.gamma * r_dom_b:
                            pair_dom = (min(dom_idx, light_b), max(dom_idx, light_b))
                            violations.append((d_scaled / r_dom_b, pair_dom, heavy_p, d_scaled, r_dom_b))

        if not violations:
            return None

        # Deterministic tie-breaking: sort by most acute violation (smallest ratio)
        violations.sort(key=lambda item: (item[0], item[1], item[2]))
        _, best_pair, best_k, best_d, best_r = violations[0]
        return best_pair, best_k, best_d, best_r

    def _build_optimal_tree_impl(
        self,
        state: ExtendedPhaseState,
        gamma_build_override: Optional[float] = None,
        forbidden_leaf_pairs: Optional[Set[Tuple[int, int]]] = None,
    ) -> JacobiTree:
        """
        Internal agglomerative builder supporting threshold override (gamma_build_override)
        and forbidden leaf-pair constraints for time-reversed adjoint hysteresis.
        """
        N = state.N
        eff_gamma_build = self.gamma_build if gamma_build_override is None else float(gamma_build_override)
        node_id_counter = 0

        def _next_id() -> int:
            nonlocal node_id_counter
            cur = node_id_counter
            node_id_counter += 1
            return cur

        # 1. Base initialization: each body forms an elementary leaf node
        active_nodes: List[JacobiTopologyNode] = []
        active_barycenters: List[np.ndarray] = []

        for i in range(N):
            m = float(state.masses[i])
            leaf = JacobiTopologyNode(
                id=_next_id(),
                is_leaf=True,
                body_index=i,
                M=m,
                mu=0.0,
                c_left=0.0,
                c_right=0.0,
                body_indices=(i,),
            )
            active_nodes.append(leaf)
            active_barycenters.append(np.copy(state.q[i]))

        # 2. Agglomerative pairwise merging with Delaunay isolation guard on leaf pairs
        while len(active_nodes) > 1:
            min_dist = float("inf")
            min_rdot = float("inf")
            best_pair = (0, 1)
            found_admissible = False

            for i in range(len(active_nodes)):
                for j in range(i + 1, len(active_nodes)):
                    dist = float(np.linalg.norm(active_barycenters[i] - active_barycenters[j]))
                    is_cand_leaf_pair = active_nodes[i].is_leaf and active_nodes[j].is_leaf
                    cand_admissible = True
                    rdot = 0.0
                    if is_cand_leaf_pair and dist > 0.0:
                        b_i_vel = active_nodes[i].body_index
                        b_j_vel = active_nodes[j].body_index
                        if b_i_vel is not None and b_j_vel is not None:
                            v_bi = state.p[b_i_vel] / max(float(state.masses[b_i_vel]), 1e-30)
                            v_bj = state.p[b_j_vel] / max(float(state.masses[b_j_vel]), 1e-30)
                            rdot = float(np.dot(state.q[b_i_vel] - state.q[b_j_vel], v_bi - v_bj)) / dist
                    if is_cand_leaf_pair and N > 2 and dist > 0.0:
                        b_i = active_nodes[i].body_index
                        b_j = active_nodes[j].body_index
                        assert b_i is not None and b_j is not None
                        pair_key = (min(b_i, b_j), max(b_i, b_j))
                        if forbidden_leaf_pairs is not None and pair_key in forbidden_leaf_pairs:
                            cand_admissible = False
                        else:
                            m_bi = float(state.masses[b_i])
                            m_bj = float(state.masses[b_j])
                            m_ij = max(m_bi + m_bj, 1e-30)
                            is_hier_cand = min(m_bi, m_bj) < 0.01 * max(m_bi, m_bj)
                            d_ext = float("inf")
                            for k in range(N):
                                if k == b_i or k == b_j:
                                    continue
                                d_ik = float(np.linalg.norm(state.q[b_i] - state.q[k]))
                                d_jk = float(np.linalg.norm(state.q[b_j] - state.q[k]))
                                m_k = max(float(state.masses[k]), 1e-30)
                                raw_hill = (m_ij / (2.0 * m_k)) ** (1.0 / 3.0)
                                if is_hier_cand:
                                    hill_scale = raw_hill
                                elif eff_gamma_build >= 1.0 - 1e-15:
                                    hill_scale = 1.0
                                elif eff_gamma_build > 0.70:
                                    w_hyst = (1.0 - eff_gamma_build) / 0.30
                                    hill_scale = 1.0 + w_hyst * (max(1.0, raw_hill) - 1.0)
                                else:
                                    hill_scale = max(1.0, raw_hill)
                                d_ext = min(d_ext, min(d_ik, d_jk) * hill_scale)
                            if d_ext < eff_gamma_build * dist:
                                cand_admissible = False

                    if cand_admissible and not found_admissible:
                        found_admissible = True
                        min_dist = dist
                        min_rdot = rdot
                        best_pair = (i, j)
                    elif cand_admissible == found_admissible:
                        if dist < min_dist * (1.0 - 1e-14):
                            min_dist = dist
                            min_rdot = rdot
                            best_pair = (i, j)
                        elif abs(dist - min_dist) <= 1e-14 * max(dist, min_dist, 1e-30) and rdot < min_rdot - 1e-14:
                            min_dist = dist
                            min_rdot = rdot
                            best_pair = (i, j)

            i, j = best_pair
            node_i = active_nodes[i]
            node_j = active_nodes[j]

            # Deterministic canonical ordering by minimal body index
            if node_i.body_indices[0] <= node_j.body_indices[0]:
                left_child, right_child = node_i, node_j
                R_l, R_r = active_barycenters[i], active_barycenters[j]
            else:
                left_child, right_child = node_j, node_i
                R_l, R_r = active_barycenters[j], active_barycenters[i]

            M_new = left_child.M + right_child.M
            mu_new = (left_child.M * right_child.M) / M_new
            c_l = left_child.M / M_new
            c_r = right_child.M / M_new
            R_bary = c_l * R_l + c_r * R_r

            parent_node = JacobiTopologyNode(
                id=_next_id(),
                is_leaf=False,
                left=left_child,
                right=right_child,
                M=M_new,
                mu=mu_new,
                c_left=c_l,
                c_right=c_r,
                body_indices=tuple(sorted(left_child.body_indices + right_child.body_indices)),
            )

            # Remove merged children (higher index first to preserve indexing)
            idx_max = max(i, j)
            idx_min = min(i, j)
            del active_nodes[idx_max]
            del active_barycenters[idx_max]
            del active_nodes[idx_min]
            del active_barycenters[idx_min]

            active_nodes.append(parent_node)
            active_barycenters.append(R_bary)

        root = active_nodes[0]
        all_nodes = root.collect_internal_nodes() + root.collect_leaves()
        node_map = {n.id: n for n in all_nodes}
        internal_nodes = tuple(root.collect_internal_nodes())
        leaf_pairs = tuple(n for n in internal_nodes if n.is_leaf_pair())
        composite_nodes = tuple(n for n in internal_nodes if n.is_composite_node())

        return JacobiTree(
            root=root,
            masses=state.masses,
            node_map=node_map,
            internal_nodes=internal_nodes,
            leaf_pairs=leaf_pairs,
            composite_nodes=composite_nodes,
        )

    def build_optimal_tree(self, state: ExtendedPhaseState) -> JacobiTree:
        """
        Constructs an optimal binary Jacobi tree for the current Cartesian configuration
        using agglomerative hierarchical clustering on cluster barycenters.

        Guarantees that the closest pairs are clustered at the bottom as leaf pairs,
        maximizing the apex isolation ratio.
        """
        return self._build_optimal_tree_impl(state)

    def evaluate_reverse_tree_adaptation(
        self,
        state_start: ExtendedPhaseState,
        state_end: ExtendedPhaseState,
        current_tree: JacobiTree,
    ) -> Tuple[JacobiTree, bool]:
        """
        Adjoint time-reversed atlas transition operator A_- under dt -> -dt,
        satisfying the group inversion identity (B_{+dt} o A_+)^{-1} = A_- o B_{-dt}.

        Detects whether the backward step [state_start -> state_end] crosses a forward
        hysteresis boundary in reverse (from state_start outside Adm(T_prev) to
        state_end inside Adm(T_prev), with inverted entry/exit thresholds r_in <-> r_out).
        """
        N = state_start.N
        if N <= 2:
            return current_tree, False

        cur_leaves = tuple(current_tree.get_leaf_pair_indices())

        candidates: List[JacobiTree] = []

        # 1. Candidate from outer threshold r_out (re-entry into leaf pair across gamma)
        candidates.append(self._build_optimal_tree_impl(state_end, gamma_build_override=self.gamma))
        # 2. Standard optimal tree at state_end
        candidates.append(self._build_optimal_tree_impl(state_end))

        # 3. Inverted inner threshold r_in exit: if current_tree has any secondary leaf pair
        # (p, b) that is outside r_in at state_end, build agglomerative candidate forbidding (p, b)
        dom_idx = int(np.argmax([float(state_end.masses[b]) for b in range(N)]))
        forbidden_exited: Set[Tuple[int, int]] = set()
        for p_idx, b_idx in cur_leaves:
            if p_idx == dom_idx or b_idx == dom_idx:
                continue
            m_p_val = float(state_end.masses[p_idx])
            m_b_val = float(state_end.masses[b_idx])
            if m_p_val >= m_b_val:
                heavy_p, light_b = p_idx, b_idx
            else:
                heavy_p, light_b = b_idx, p_idx
            r_pb = float(np.linalg.norm(state_end.q[heavy_p] - state_end.q[light_b]))
            r_dom_b = float(np.linalg.norm(state_end.q[dom_idx] - state_end.q[light_b]))
            m_dom_b = max(float(state_end.masses[dom_idx]) + float(state_end.masses[light_b]), 1e-30)
            m_hp = max(float(state_end.masses[heavy_p]), 1e-30)
            d_scaled = r_pb * ((m_dom_b / (2.0 * m_hp)) ** (1.0 / 3.0))
            if r_dom_b > 0.0 and d_scaled >= self.gamma * r_dom_b:
                pair_key = (min(p_idx, b_idx), max(p_idx, b_idx))
                forbidden_exited.add(pair_key)
                candidates.append(
                    self._build_optimal_tree_impl(state_end, forbidden_leaf_pairs={pair_key})
                )
        if len(forbidden_exited) > 1:
            candidates.append(
                self._build_optimal_tree_impl(state_end, forbidden_leaf_pairs=forbidden_exited)
            )

        # 4. General comparable-mass prior-chart candidates: forbid active leaf pair(s)
        for p_idx, b_idx in cur_leaves:
            pair_key = (min(p_idx, b_idx), max(p_idx, b_idx))
            cand_alt = self._build_optimal_tree_impl(
                state_end, gamma_build_override=self.gamma, forbidden_leaf_pairs={pair_key}
            )
            candidates.append(cand_alt)
            for alt_p, alt_b in cand_alt.get_leaf_pair_indices():
                alt_key = (min(alt_p, alt_b), max(alt_p, alt_b))
                if alt_key != pair_key:
                    candidates.append(
                        self._build_optimal_tree_impl(
                            state_end,
                            gamma_build_override=self.gamma,
                            forbidden_leaf_pairs={pair_key, alt_key},
                        )
                    )

        seen_leaf_sets: Set[Tuple[Tuple[int, int], ...]] = {cur_leaves}
        for cand_tree in candidates:
            cand_leaves = tuple(cand_tree.get_leaf_pair_indices())
            if cand_leaves in seen_leaf_sets:
                continue
            seen_leaf_sets.add(cand_leaves)

            # Condition A: T_cand must be admissible at state_end (t_{n-1})
            if self.check_tree_admissibility(state_end, cand_tree, ignore_dwell=True) is not None:
                continue

            # Condition B: T_cand must be inadmissible at state_start (t_n)
            viol_start = self.check_tree_admissibility(state_start, cand_tree, ignore_dwell=True)
            if viol_start is None:
                continue

            # Condition C: Forward adaptation at state_start must map to current_tree
            fwd_opt = self.build_optimal_tree(state_start)
            if tuple(fwd_opt.get_leaf_pair_indices()) != cur_leaves:
                continue

            best_pair, best_k, best_d, best_r = viol_start
            dwell = abs(state_start.t - self.transition_history[-1].t) if self.transition_history else None
            record = AtlasTransitionRecord(
                t=state_start.t,
                old_leaf_pairs=cur_leaves,
                new_leaf_pairs=cand_leaves,
                triggering_pair=best_pair,
                intruder_body=best_k,
                separation_ratio=best_d / max(best_r, 1e-15),
                dwell_time=dwell,
            )
            self.transition_history.append(record)
            return cand_tree, True

        # Mandatory admissibility safeguard: if current_tree itself has become invalid at state_end,
        # we must switch to a valid chart (preferring a candidate that just became admissible at state_end,
        # or falling back to build_optimal_tree(state_end)).
        viol_cur = self.check_tree_admissibility(state_end, current_tree)
        if viol_cur is not None:
            fallback_tree: Optional[JacobiTree] = None
            for cand_tree in candidates:
                cand_leaves = tuple(cand_tree.get_leaf_pair_indices())
                if cand_leaves == cur_leaves:
                    continue
                if self.check_tree_admissibility(state_end, cand_tree, ignore_dwell=True) is None:
                    if self.check_tree_admissibility(state_start, cand_tree, ignore_dwell=True) is not None:
                        fallback_tree = cand_tree
                        break
            if fallback_tree is None:
                fallback_tree = self.build_optimal_tree(state_end)
            new_leaves = tuple(fallback_tree.get_leaf_pair_indices())
            if new_leaves != cur_leaves:
                best_pair, best_k, best_d, best_r = viol_cur
                dwell = abs(state_end.t - self.transition_history[-1].t) if self.transition_history else None
                record = AtlasTransitionRecord(
                    t=state_end.t,
                    old_leaf_pairs=cur_leaves,
                    new_leaf_pairs=new_leaves,
                    triggering_pair=best_pair,
                    intruder_body=best_k,
                    separation_ratio=best_d / max(best_r, 1e-15),
                    dwell_time=dwell,
                )
                self.transition_history.append(record)
                return fallback_tree, True

        return current_tree, False

    def evaluate_tree_adaptation(
        self,
        state: ExtendedPhaseState,
        current_tree: JacobiTree,
        reverse_time: bool = False,
        dt_step: float = -1e-3,
    ) -> Tuple[JacobiTree, bool]:
        """
        Evaluates topological admissibility and executes a symplectic chart transition if needed.

        Returns:
            (adapted_tree, was_flipped)
        """
        if reverse_time:
            dt_probe = -abs(float(dt_step)) if dt_step != 0.0 else -1e-3
            v = state.p / state.masses[:, None]
            state_end_probe = state.with_updates(q=state.q + dt_probe * v, t=state.t + dt_probe)
            return self.evaluate_reverse_tree_adaptation(state, state_end_probe, current_tree)

        violation = self.check_tree_admissibility(state, current_tree)
        if violation is None:
            return current_tree, False

        best_pair, best_k, best_d, best_r = violation
        new_tree = self.build_optimal_tree(state)
        if tuple(new_tree.get_leaf_pair_indices()) == tuple(current_tree.get_leaf_pair_indices()):
            return current_tree, False

        # Log transition record with dwell time tracking
        dwell = abs(state.t - self.transition_history[-1].t) if self.transition_history else None
        record = AtlasTransitionRecord(
            t=state.t,
            old_leaf_pairs=tuple(current_tree.get_leaf_pair_indices()),
            new_leaf_pairs=tuple(new_tree.get_leaf_pair_indices()),
            triggering_pair=best_pair,
            intruder_body=best_k,
            separation_ratio=best_d / max(best_r, 1e-15),
            dwell_time=dwell,
        )
        self.transition_history.append(record)

        return new_tree, True

    def reset_history(self) -> None:
        """Clears the transition record history."""
        self.transition_history.clear()
