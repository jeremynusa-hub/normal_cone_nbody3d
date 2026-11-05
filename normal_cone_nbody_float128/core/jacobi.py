r"""
Immutable Combinatorial Jacobi Binary Trees, Mass Metrics & Hyper-Radii.

Mathematical Monograph Alignment:
- Chapter 1: Section 1.4 (\label{sec:jacobi_trees}) — Canonical Jacobi coordinates and kinetic diagonal cometric.
- Chapter 2: Section 2.1 & Section 2.2 — Jacobi cometric and coordinate trees.
- Chapter 4: Section 4.1 & Section 4.3, Equations (4.13)-(4.22) — Maximal flag Jacobi trees T in Trees(N),
  internal nodes, leaf pairs, and dimensionally normalized hyper-radii bar_rho_C = rho_C / sqrt(M_C) in [L].
- Proposition 4.1: Disjoint cotangent commutation {H_A, H_B} = 0 for disjoint subtrees.
- Certified Lean 4 Alignment: PaperUnifiedFoundations.lean
  (jacobi_node_defines_hierarchy, hyper_radius_sq_pos, hyper_radius_sq_zero_iff,
   disjoint_cotangent_commutation_local).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Dict, Tuple, List, Union
import math
import numpy as np
import mpmath as mp

from .constants import to_f128, to_f64


@dataclass(frozen=True)
class JacobiTopologyNode:
    """
    Immutable topological node of a complete binary Jacobi tree T in Trees(N).

    Strictly decoupled from dynamic phase coordinates (positions r and momenta w are
    stored externally in JacobiCoordinates).

    Attributes:
        id: Unique positive integer identifier for the node within tree T.
        is_leaf: True if the node represents a single physical point mass.
        body_index: Physical particle index in {0, ..., N-1} (None for internal nodes).
        left: Left child sub-tree (None for leaves).
        right: Right child sub-tree (None for leaves).
        body_indices: Ordered tuple of particle indices belonging to this cluster.
        M: Total cluster mass M_v = sum_{i in body_indices} m_i.
        mu: Reduced mass mu_v = (M_left * M_right) / M_v (0.0 for leaves).
        c_left: Mass fraction c_left = M_left / M_v (0.0 for leaves).
        c_right: Mass fraction c_right = M_right / M_v (0.0 for leaves).
    """
    id: int
    is_leaf: bool
    body_index: Optional[int] = None
    left: Optional[JacobiTopologyNode] = None
    right: Optional[JacobiTopologyNode] = None
    body_indices: Tuple[int, ...] = field(default_factory=tuple)
    M: Union[float, mp.mpf] = 0.0
    mu: Union[float, mp.mpf] = 0.0
    c_left: Union[float, mp.mpf] = 0.0
    c_right: Union[float, mp.mpf] = 0.0

    def is_leaf_pair(self) -> bool:
        """Returns True if both children are leaves (represents a physical binary pair)."""
        return (not self.is_leaf) and (self.left is not None and self.left.is_leaf) and (self.right is not None and self.right.is_leaf)

    def is_composite_node(self) -> bool:
        """Returns True if the internal node is not a leaf pair (merges clusters)."""
        return (not self.is_leaf) and (not self.is_leaf_pair())

    def collect_internal_nodes(self) -> List[JacobiTopologyNode]:
        """Collects all internal nodes in post-order traversal (bottom-up)."""
        if self.is_leaf:
            return []
        nodes: List[JacobiTopologyNode] = []
        if self.left is not None:
            nodes.extend(self.left.collect_internal_nodes())
        if self.right is not None:
            nodes.extend(self.right.collect_internal_nodes())
        nodes.append(self)
        return nodes

    def collect_leaves(self) -> List[JacobiTopologyNode]:
        """Collects all leaf nodes."""
        if self.is_leaf:
            return [self]
        leaves: List[JacobiTopologyNode] = []
        if self.left is not None:
            leaves.extend(self.left.collect_leaves())
        if self.right is not None:
            leaves.extend(self.right.collect_leaves())
        return leaves


@dataclass(frozen=True)
class JacobiTree:
    """
    Complete immutable binary Jacobi tree T in Trees(N) equipped with the mass metric.

    Contains N leaves and N - 1 internal nodes, defining an adapted cotangent chart
    on T*Q[N].
    """
    root: JacobiTopologyNode
    masses: np.ndarray
    node_map: Dict[int, JacobiTopologyNode]
    internal_nodes: Tuple[JacobiTopologyNode, ...]
    leaf_pairs: Tuple[JacobiTopologyNode, ...]
    composite_nodes: Tuple[JacobiTopologyNode, ...]

    @property
    def is_f128(self) -> bool:
        """Returns True if the tree masses and nodes are evaluated in extended precision."""
        return self.masses.dtype == object or isinstance(self.root.M, mp.mpf)

    @property
    def N(self) -> int:
        """Number of point masses N."""
        return len(self.masses)

    @property
    def num_internal(self) -> int:
        """Number of internal Jacobi nodes: N - 1."""
        return len(self.internal_nodes)

    def get_node(self, node_id: int) -> JacobiTopologyNode:
        """Retrieves node by ID."""
        return self.node_map[node_id]

    def get_leaf_pair_indices(self) -> List[Tuple[int, int]]:
        """Returns the list of body index pairs (i, j) for all leaf pairs."""
        return [(node.body_indices[0], node.body_indices[1]) for node in self.leaf_pairs]

    def compute_hyper_radius(self, r_dict: Dict[int, np.ndarray]) -> Union[float, mp.mpf]:
        """
        Computes the global cluster hyper-radius rho_C (dimension [M]^(1/2) [L]):
            rho_C = sqrt(sum_{v in V_int(T)} mu_v ||r_v||^2).
        """
        if self.is_f128:
            rho_sq = mp.mpf('0.0')
            for node in self.internal_nodes:
                r_v = r_dict[node.id]
                dot_v = r_v[0]**2 + r_v[1]**2 + r_v[2]**2
                rho_sq += node.mu * dot_v
            return mp.sqrt(max(mp.mpf('0.0'), rho_sq))
        else:
            rho_sq = 0.0
            for node in self.internal_nodes:
                r_v = r_dict[node.id]
                rho_sq += node.mu * float(np.dot(r_v, r_v))
            return float(math.sqrt(max(0.0, rho_sq)))

    def compute_normalized_hyper_radius(self, r_dict: Dict[int, np.ndarray]) -> Union[float, mp.mpf]:
        """
        Computes the normalized hyper-radius bar_rho_C (dimension [L]):
            bar_rho_C = rho_C / sqrt(M_tot) in [L].
        Guarantees that the scale ratio u_v = ||r_v|| / (bar_rho_L + bar_rho_R) is strictly dimensionless.
        """
        rho_c = self.compute_hyper_radius(r_dict)
        if self.is_f128:
            return rho_c / mp.sqrt(self.root.M)
        return float(rho_c / math.sqrt(self.root.M))

    def compute_hyper_radius_cartesian(self, q: np.ndarray) -> Union[float, mp.mpf]:
        """
        Computes global hyper-radius directly from Cartesian positions:
            rho_C = sqrt((1 / M_tot) * sum_{i < j} m_i m_j ||q_i - q_j||^2).
        """
        if self.is_f128:
            N = self.N
            total_sq = mp.mpf('0.0')
            for i in range(N):
                for j in range(i + 1, N):
                    diff = q[i] - q[j]
                    d2 = diff[0]**2 + diff[1]**2 + diff[2]**2
                    total_sq += self.masses[i] * self.masses[j] * d2
            return mp.sqrt(max(mp.mpf('0.0'), total_sq / self.root.M))
        else:
            N = self.N
            diff = q[:, np.newaxis, :] - q[np.newaxis, :, :]
            dist_sq = np.sum(diff ** 2, axis=-1)
            mass_matrix = np.outer(self.masses, self.masses)
            i_idx, j_idx = np.triu_indices(N, k=1)
            total_sq = np.sum(mass_matrix[i_idx, j_idx] * dist_sq[i_idx, j_idx])
            return float(math.sqrt(max(0.0, total_sq / self.root.M)))

    def compute_normalized_hyper_radius_cartesian(self, q: np.ndarray) -> Union[float, mp.mpf]:
        """Computes normalized hyper-radius bar_rho_C in [L] from Cartesian positions."""
        rho_c = self.compute_hyper_radius_cartesian(q)
        if self.is_f128:
            return rho_c / mp.sqrt(self.root.M)
        return float(rho_c / math.sqrt(self.root.M))


def _build_node_recursive(
    spec: Union[int, tuple],
    masses: np.ndarray,
    id_counter: List[int],
    node_map: Dict[int, JacobiTopologyNode],
) -> JacobiTopologyNode:
    """Internal recursive constructor of JacobiTopologyNode from index tuple."""
    node_id = id_counter[0]
    id_counter[0] += 1

    is_f128 = masses.dtype == object or (len(masses) > 0 and isinstance(masses[0], mp.mpf))

    if isinstance(spec, (int, np.integer)):
        body_idx = int(spec)
        m = masses[body_idx] if is_f128 else float(masses[body_idx])
        zero_val = mp.mpf('0.0') if is_f128 else 0.0
        node = JacobiTopologyNode(
            id=node_id,
            is_leaf=True,
            body_index=body_idx,
            body_indices=(body_idx,),
            M=m,
            mu=zero_val,
            c_left=zero_val,
            c_right=zero_val,
        )
        node_map[node_id] = node
        return node

    if not isinstance(spec, (tuple, list)) or len(spec) != 2:
        raise ValueError(f"Tree specification node must be an int or a pair (left, right), got {spec}.")

    left_child = _build_node_recursive(spec[0], masses, id_counter, node_map)
    right_child = _build_node_recursive(spec[1], masses, id_counter, node_map)

    m_tot = left_child.M + right_child.M
    mu_val = (left_child.M * right_child.M) / m_tot
    c_left = left_child.M / m_tot
    c_right = right_child.M / m_tot
    body_indices = tuple(sorted(left_child.body_indices + right_child.body_indices))

    node = JacobiTopologyNode(
        id=node_id,
        is_leaf=False,
        body_index=None,
        left=left_child,
        right=right_child,
        body_indices=body_indices,
        M=m_tot,
        mu=mu_val,
        c_left=c_left,
        c_right=c_right,
    )
    node_map[node_id] = node
    return node


def build_tree_from_spec(masses: np.ndarray, spec: Union[int, tuple]) -> JacobiTree:
    """
    Constructs a complete JacobiTree from a recursive nested tuple specification.
    Example for N=3: ((0, 1), 2)
    Example for N=4: (((0, 1), 2), 3) or ((0, 1), (2, 3))
    """
    is_f128 = (isinstance(masses, np.ndarray) and masses.dtype == object) or (len(masses) > 0 and isinstance(masses[0], mp.mpf))
    if not is_f128:
        masses = np.asarray(masses, dtype=np.float64)
    else:
        masses = np.asarray(masses, dtype=object)

    if masses.ndim != 1 or len(masses) < 2:
        raise ValueError("Masses must be a 1D array of at least 2 particles.")

    id_counter = [0]
    node_map: Dict[int, JacobiTopologyNode] = {}
    root = _build_node_recursive(spec, masses, id_counter, node_map)

    internal_nodes = tuple(root.collect_internal_nodes())
    leaf_pairs = tuple(n for n in internal_nodes if n.is_leaf_pair())
    composite_nodes = tuple(n for n in internal_nodes if n.is_composite_node())

    return JacobiTree(
        root=root,
        masses=masses,
        node_map=node_map,
        internal_nodes=internal_nodes,
        leaf_pairs=leaf_pairs,
        composite_nodes=composite_nodes,
    )


def build_canonical_tree(
    masses: np.ndarray,
    topology: str = "comb",
) -> JacobiTree:
    """
    Builds standard canonical binary tree topologies:
    - 'comb' / 'linear': (((0, 1), 2), 3)...
    - 'balanced': Recursive median-bipartition tree.
    """
    N = len(masses)
    if N < 2:
        raise ValueError(f"Jacobi tree requires N >= 2, got N={N}.")

    if topology in ("comb", "linear"):
        spec: Union[int, tuple] = 0
        for i in range(1, N):
            spec = (spec, i)
        return build_tree_from_spec(masses, spec)

    if topology == "balanced":
        def _build_balanced_spec(indices: List[int]) -> Union[int, tuple]:
            if len(indices) == 1:
                return indices[0]
            mid = len(indices) // 2
            return (_build_balanced_spec(indices[:mid]), _build_balanced_spec(indices[mid:]))

        spec = _build_balanced_spec(list(range(N)))
        return build_tree_from_spec(masses, spec)

    raise ValueError(f"Unknown topology '{topology}'. Choose 'comb' or 'balanced'.")


def enumerate_all_jacobi_trees(masses: np.ndarray) -> List[JacobiTree]:
    """
    Combinatorially enumerates all (2N - 3)!! distinct binary Jacobi trees with N labeled leaves.
    Counts:
    - N = 2: 1 tree ((0, 1))
    - N = 3: 3 trees (((0, 1), 2), ((0, 2), 1), ((1, 2), 0))
    - N = 4: 15 trees
    - N = 5: 105 trees
    """
    N = len(masses)
    if N < 2:
        return []

    def _all_trees_for_set(leaves: Tuple[int, ...]) -> List[tuple]:
        if len(leaves) == 1:
            return [leaves[0]]
        if len(leaves) == 2:
            return [(leaves[0], leaves[1])]

        results: List[tuple] = []
        # Partition leaves into two non-empty subsets (order irrelevant within pair)
        # Fix the first element to left subset to prevent symmetry duplication A <-> B
        first = leaves[0]
        rest = leaves[1:]
        n_rest = len(rest)

        # Iterate over all non-empty proper subsets for the left side
        for mask in range(0, (1 << n_rest) - 1):
            left_subset = [first]
            right_subset = []
            for idx in range(n_rest):
                if (mask >> idx) & 1:
                    left_subset.append(rest[idx])
                else:
                    right_subset.append(rest[idx])

            left_tuple = tuple(sorted(left_subset))
            right_tuple = tuple(sorted(right_subset))

            for l_tree in _all_trees_for_set(left_tuple):
                for r_tree in _all_trees_for_set(right_tuple):
                    results.append((l_tree, r_tree))

        return results

    all_specs = _all_trees_for_set(tuple(range(N)))
    return [build_tree_from_spec(masses, spec) for spec in all_specs]
