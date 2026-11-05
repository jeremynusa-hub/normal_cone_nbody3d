r"""
Symplectic Lifts of Jacobi Tree Rotations in Sp(6N-6, R) via SL(2, R).

Mathematical Monograph Alignment:
- Chapter 3: Section 3.4, Equations (3.22)-(3.28) — Symplectic chart transitions between adjacent Jacobi trees.
- Chapter 4: Section 4.5 — Symplectic atlas structure and transition mappings.
- Chapter 6: Section 6.1, Proposition 6.1 — Exact unimodular SL(2, R) representation and Darboux invariance.
- Certified Lean 4 Alignment: PaperUnifiedFoundations.lean
  (sl2_preserves_symplectic_area, AtlasChartTransitionSymplectic.lean).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, Dict, Optional, Set, Sequence, TYPE_CHECKING
import numpy as np

from .coordinates import JacobiCoordinates, cartesian_to_jacobi, jacobi_to_cartesian
from .jacobi import JacobiTree, JacobiTopologyNode
from .state import ExtendedPhaseState

if TYPE_CHECKING:
    from ..potential.screened_apex import ScreenedApexPotential
    from ..integrators.correctors import SymplecticCorrector


def _identify_tree_rotation_clusters(
    tree_A: JacobiTree,
    tree_B: JacobiTree,
) -> Tuple[
    JacobiTopologyNode,
    JacobiTopologyNode,
    JacobiTopologyNode,
    JacobiTopologyNode,
    Set[int],
    Set[int],
    Set[int],
]:
    """
    Identifies the internal nodes (node_A1, node_A2) in tree_A, (node_B1, node_B2) in tree_B,
    and the three disjoint body-index clusters (cluster_A, cluster_B, cluster_C) undergoing
    an elementary Jacobi tree rotation:
        ((A, B), C) <---> (A, (B, C)),
    where:
    - cluster_B is the pivot cluster shared by both inner nodes node_A1 and node_B1,
    - cluster_A is the sibling of cluster_B in node_A1 (which becomes outer in node_B2),
    - cluster_C is the sibling of cluster_B in node_B1 (which was outer in node_A2).
    """
    map_A = {n.body_indices: n for n in tree_A.internal_nodes}
    map_B = {n.body_indices: n for n in tree_B.internal_nodes}

    diff_A = [n for idx, n in map_A.items() if idx not in map_B]
    diff_B = [n for idx, n in map_B.items() if idx not in map_A]

    if not diff_A or not diff_B:
        # Identical tree topology: fall back to the first two internal nodes if available
        if len(tree_A.internal_nodes) >= 2 and len(tree_B.internal_nodes) >= 2:
            nA1, nA2 = tree_A.internal_nodes[0], tree_A.internal_nodes[-1]
            nB1, nB2 = tree_B.internal_nodes[0], tree_B.internal_nodes[-1]
            assert nA1.left is not None and nA1.right is not None
            cA = set(nA1.left.body_indices)
            cB = set(nA1.right.body_indices)
            cC = set(nA2.body_indices) - (cA | cB)
            return nA1, nA2, nB1, nB2, cA, cB, cC
        raise ValueError("Trees must have at least 2 internal nodes (N >= 3) for SL(2, R) rotation.")

    node_A1 = min(diff_A, key=lambda n: len(n.body_indices))
    node_B1 = min(diff_B, key=lambda n: len(n.body_indices))

    # Find parent node_A2 of node_A1 in tree_A and parent node_B2 of node_B1 in tree_B
    node_A2 = tree_A.root
    for n in tree_A.internal_nodes:
        if (n.left is not None and n.left.id == node_A1.id) or (n.right is not None and n.right.id == node_A1.id):
            node_A2 = n
            break

    node_B2 = tree_B.root
    for n in tree_B.internal_nodes:
        if (n.left is not None and n.left.id == node_B1.id) or (n.right is not None and n.right.id == node_B1.id):
            node_B2 = n
            break

    assert node_A1.left is not None and node_A1.right is not None
    assert node_B1.left is not None and node_B1.right is not None

    children_A1 = [set(node_A1.left.body_indices), set(node_A1.right.body_indices)]
    children_B1 = [set(node_B1.left.body_indices), set(node_B1.right.body_indices)]

    cluster_B: Optional[Set[int]] = None
    cluster_A: Optional[Set[int]] = None
    for c_cand in children_A1:
        if c_cand in children_B1:
            cluster_B = c_cand
        else:
            cluster_A = c_cand

    if cluster_B is None or cluster_A is None:
        cluster_A = children_A1[0]
        cluster_B = children_A1[1]

    cluster_C = (set(node_B1.body_indices) - cluster_B) or (set(node_A2.body_indices) - (cluster_A | cluster_B))
    return node_A1, node_A2, node_B1, node_B2, cluster_A, cluster_B, cluster_C


@dataclass(frozen=True)
class SL2TransitionMatrix:
    """
    Elementary canonical Jacobi chart transition matrix in SL(2, R) governing
    the local tree rotation uniting three clusters (A, B, C):
        ((A, B), C) <---> (A, (B, C)).

    Unimodular condition:
        det(M) = a * d - b * c == 1.0.

    Symplectic action on positions:
        r'_1 = a * r_1 + b * r_2
        r'_2 = c * r_1 + d * r_2

    Contragredient dual action on canonical momenta:
        w'_1 =  d * w_1 - c * w_2
        w'_2 = -b * w_1 + a * w_2

    Preserves the canonical Darboux 2-form:
        <r'_1, w'_1> + <r'_2, w'_2> == <r_1, w_1> + <r_2, w_2>.
    """
    a: float
    b: float
    c: float
    d: float

    def __post_init__(self) -> None:
        det = self.a * self.d - self.b * self.c
        if abs(det - 1.0) > 1e-14:
            raise ValueError(f"Transition matrix must strictly belong to SL(2, R), got det={det:.16e}.")

    @property
    def determinant(self) -> float:
        """Unimodular determinant: ad - bc = 1.0."""
        return float(self.a * self.d - self.b * self.c)

    def transform_positions(self, r1: np.ndarray, r2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Transforms spatial relative Jacobi vectors."""
        r1_new = self.a * r1 + self.b * r2
        r2_new = self.c * r1 + self.d * r2
        return r1_new, r2_new

    def transform_momenta(self, w1: np.ndarray, w2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Transforms conjugate relative momentum vectors."""
        w1_new = self.d * w1 - self.c * w2
        w2_new = -self.b * w1 + self.a * w2
        return w1_new, w2_new

    def inverse(self) -> SL2TransitionMatrix:
        """Returns the exact inverse transition matrix in SL(2, R): [[d, -b], [-c, a]]."""
        return SL2TransitionMatrix(a=self.d, b=-self.b, c=-self.c, d=self.a)

    @classmethod
    def from_masses(cls, mA: float, mB: float, mC: float) -> SL2TransitionMatrix:
        """
        Computes exact mass-dependent coefficients for the canonical rotation:
            ((A, B), C)  |-->  (A, (B, C))

        Aligned with Theorem 4.4 Eq. (4.116) and Proposition 6.5 Eq. (6.178) of the monograph:
            M_coord = [[ -mA / (mA + mB),               1.0             ],
                       [ -(mB * mABC) / (mAB * mBC),   -mC / (mB + mC)  ]]

        Note on orientation:
            In `coordinates.py`, each internal node v defines r_v = X_L(v) - X_R(v).
            For T_A = ((A, B), C) with r_1 = X_A - X_B and r_2 = X_AB - X_C, applying M_coord yields:
                r'_1 = a * r_1 + b * r_2 = X_B - X_C
                r'_2 = c * r_1 + d * r_2 = -(X_A - X_BC) = X_BC - X_A
            In `rotate_jacobi_tree_coordinates`, `sign1 = +1` and `sign2 = -1` restore
            r_{B,2} = X_A - X_BC for the root node of T_B = (A, (B, C)).

        Formulas:
            a = -mA / (mA + mB)
            b = 1.0
            c = -(mB * mABC) / ((mA + mB) * (mB + mC))
            d = -mC / (mB + mC)

        Proof:
            det(M) = a * d - b * c
                   = (mA * mC + mB * mABC) / ((mA + mB)(mB + mC))
                   = (mA * mC + mB * (mA + mB + mC)) / (mA * mB + mA * mC + mB^2 + mB * mC)
                   = 1.0  (Exact algebraic identity in SL(2, R)).
        """
        if mA <= 0.0 or mB <= 0.0 or mC <= 0.0:
            raise ValueError("Sub-cluster masses mA, mB, mC must be strictly positive.")

        mAB = mA + mB
        mBC = mB + mC
        mABC = mA + mB + mC

        a = -mA / mAB
        b = 1.0
        c = -(mB * mABC) / (mAB * mBC)
        d = -mC / mBC

        return cls(a=a, b=b, c=c, d=d)

    @classmethod
    def from_trees(
        cls,
        tree_A: JacobiTree,
        tree_B: JacobiTree,
        masses: Optional[np.ndarray] = None,
    ) -> SL2TransitionMatrix:
        """
        Extracts the exact unimodular SL(2, R) transition matrix M_{A -> B} directly
        from two adjacent Jacobi trees tree_A and tree_B uniting three clusters (A, B, C).

        If tree_A and tree_B share identical internal topologies, returns the identity matrix I_2.
        """
        if tree_A is tree_B or {n.body_indices for n in tree_A.internal_nodes} == {
            n.body_indices for n in tree_B.internal_nodes
        }:
            return cls(a=1.0, b=0.0, c=0.0, d=1.0)

        _, _, _, _, cluster_A, cluster_B, cluster_C = _identify_tree_rotation_clusters(tree_A, tree_B)
        mass_arr = masses if masses is not None else tree_A.masses

        mA = float(sum(float(mass_arr[i]) for i in cluster_A))
        mB = float(sum(float(mass_arr[i]) for i in cluster_B))
        mC = float(sum(float(mass_arr[i]) for i in cluster_C))

        return cls.from_masses(mA=mA, mB=mB, mC=mC)


def build_symplectic_matrix_12(trans: SL2TransitionMatrix) -> np.ndarray:
    """
    Constructs the full 12x12 symplectic representation Psi in Sp(12, R) acting on
    the 2 active 3D pairs: (r_1, r_2, w_1, w_2) in R^12.
    """
    M_pos = np.array([[trans.a, trans.b], [trans.c, trans.d]], dtype=np.float64)
    # Contragredient momentum block: M_mom = (M_pos^-1)^T = [[d, -c], [-b, a]]
    M_mom = np.array([[trans.d, -trans.c], [-trans.b, trans.a]], dtype=np.float64)

    # Kronecker product with 3x3 identity for 3D components
    I3 = np.eye(3, dtype=np.float64)
    block_pos = np.kron(M_pos, I3)   # 6x6
    block_mom = np.kron(M_mom, I3)   # 6x6

    Psi = np.zeros((12, 12), dtype=np.float64)
    Psi[0:6, 0:6] = block_pos
    Psi[6:12, 6:12] = block_mom
    return Psi


def build_tree_transition_symplectic_12(
    tree_A: JacobiTree,
    tree_B: JacobiTree,
    masses: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Constructs the oriented 12x12 symplectic transition matrix Psi_{A -> B} in Sp(12, R)
    mapping the stored internal node coordinates (r_{1,A}, r_{2,A}, w_{1,A}, w_{2,A}) of tree_A
    directly to the stored internal node coordinates (r_{1,B}, r_{2,B}, w_{1,B}, w_{2,B}) of tree_B.

    Incorporates the intrinsic node orientation signs S_A, S_B in {+1, -1}^2 so that:
    - Psi_{A -> B} in Sp(12, R) with det(Psi_{A -> B}) == +1.0,
    - Psi_{C -> A} @ Psi_{B -> C} @ Psi_{A -> B} == I_12 (exact Cech 1-cocycle closure).
    """
    if tree_A is tree_B or {n.body_indices for n in tree_A.internal_nodes} == {
        n.body_indices for n in tree_B.internal_nodes
    }:
        return np.eye(12, dtype=np.float64)

    node_A1, node_A2, node_B1, node_B2, cluster_A, cluster_B, cluster_C = _identify_tree_rotation_clusters(
        tree_A, tree_B
    )
    trans = SL2TransitionMatrix.from_trees(tree_A, tree_B, masses=masses)

    assert node_A1.left is not None and node_A2.left is not None
    assert node_B1.left is not None and node_B2.left is not None

    s1_A = +1.0 if not cluster_A.isdisjoint(node_A1.left.body_indices) else -1.0
    s2_A = +1.0 if not (cluster_A | cluster_B).isdisjoint(node_A2.left.body_indices) else -1.0
    s1_B = +1.0 if not cluster_B.isdisjoint(node_B1.left.body_indices) else -1.0
    s2_B = +1.0 if not (cluster_B | cluster_C).isdisjoint(node_B2.left.body_indices) else -1.0

    S_A = np.diag([s1_A, s2_A])
    S_B = np.diag([s1_B, s2_B])
    M_base = np.array([[trans.a, trans.b], [trans.c, trans.d]], dtype=np.float64)
    M_pos = S_B @ M_base @ S_A
    # Since S_A^-1 = S_A and S_B^-1 = S_B, contragredient block is S_B @ M_mom_base @ S_A
    M_mom_base = np.array([[trans.d, -trans.c], [-trans.b, trans.a]], dtype=np.float64)
    M_mom = S_B @ M_mom_base @ S_A

    I3 = np.eye(3, dtype=np.float64)
    Psi = np.zeros((12, 12), dtype=np.float64)
    Psi[0:6, 0:6] = np.kron(M_pos, I3)
    Psi[6:12, 6:12] = np.kron(M_mom, I3)
    return Psi


def verify_symplectic_matrix_12(Psi: np.ndarray) -> bool:
    """
    Certifies that a 12x12 matrix Psi satisfies the symplectic condition:
        Psi^T * J_12 * Psi == J_12,
    where J_12 = [[0, I_6], [-I_6, 0]].
    """
    if Psi.shape != (12, 12):
        return False

    I6 = np.eye(6, dtype=np.float64)
    J12 = np.zeros((12, 12), dtype=np.float64)
    J12[0:6, 6:12] = I6
    J12[6:12, 0:6] = -I6

    diff = Psi.T @ J12 @ Psi - J12
    return bool(np.max(np.abs(diff)) < 1e-14)


def rotate_jacobi_tree_coordinates(
    coords: JacobiCoordinates,
    tree_A: JacobiTree,
    tree_B: JacobiTree,
    trans: Optional[SL2TransitionMatrix] = None,
    node_ids_A: Optional[Tuple[int, int]] = None,
    node_ids_B: Optional[Tuple[int, int]] = None,
) -> JacobiCoordinates:
    """
    Applies the SL(2, R) tree transition to rotate coordinates between two adjacent trees:
    - node_ids_A = (v1_A, v2_A): the two internal nodes in tree_A undergoing the flip.
    - node_ids_B = (v1_B, v2_B): the corresponding transformed nodes in tree_B.
    If `trans`, `node_ids_A`, or `node_ids_B` are omitted, they are automatically inferred
    from `(tree_A, tree_B)`. All other disjoint Jacobi coordinates are carried over directly.
    """
    nA1, nA2, nB1, nB2, cluster_A, cluster_B, cluster_C = _identify_tree_rotation_clusters(tree_A, tree_B)

    if node_ids_A is None:
        v1_A, v2_A = nA1.id, nA2.id
    else:
        v1_A, v2_A = node_ids_A

    if node_ids_B is None:
        v1_B, v2_B = nB1.id, nB2.id
    else:
        v1_B, v2_B = node_ids_B

    if trans is None:
        trans = SL2TransitionMatrix.from_trees(tree_A, tree_B)

    node_A1 = tree_A.get_node(v1_A)
    node_A2 = tree_A.get_node(v2_A)
    node_B1 = tree_B.get_node(v1_B)
    node_B2 = tree_B.get_node(v2_B)

    assert node_A1.left is not None and node_A2.left is not None
    assert node_B1.left is not None and node_B2.left is not None

    if trans.b >= 0.0:
        # Canonical forward SL(2, R) form from_masses(mA, mB, mC):
        # maps (X_A - X_B, X_AB - X_C) |--> (X_B - X_C, X_BC - X_A)
        s1_in = +1.0 if not cluster_A.isdisjoint(node_A1.left.body_indices) else -1.0
        s2_in = +1.0 if not (cluster_A | cluster_B).isdisjoint(node_A2.left.body_indices) else -1.0
        s1_out = +1.0 if not cluster_B.isdisjoint(node_B1.left.body_indices) else -1.0
        s2_out = +1.0 if not (cluster_B | cluster_C).isdisjoint(node_B2.left.body_indices) else -1.0
    else:
        # Explicit inverse SL(2, R) form from_masses(mC, mB, mA).inverse():
        # maps (X_B - X_C, X_BC - X_A) |--> (X_A - X_B, X_AB - X_C)
        s1_in = +1.0 if not cluster_B.isdisjoint(node_A1.left.body_indices) else -1.0
        s2_in = +1.0 if not (cluster_A | cluster_B).isdisjoint(node_A2.left.body_indices) else -1.0
        s1_out = +1.0 if not cluster_A.isdisjoint(node_B1.left.body_indices) else -1.0
        s2_out = +1.0 if not (cluster_B | cluster_C).isdisjoint(node_B2.left.body_indices) else -1.0

    r1_A = s1_in * coords.r[v1_A]
    r2_A = s2_in * coords.r[v2_A]
    w1_A = s1_in * coords.w[v1_A]
    w2_A = s2_in * coords.w[v2_A]

    r1_B, r2_B = trans.transform_positions(r1_A, r2_A)
    w1_B, w2_B = trans.transform_momenta(w1_A, w2_A)

    r_out: Dict[int, np.ndarray] = {}
    w_out: Dict[int, np.ndarray] = {}

    # Assign transformed nodes with correct spatial orientation
    r_out[v1_B] = s1_out * r1_B
    r_out[v2_B] = s2_out * r2_B
    w_out[v1_B] = s1_out * w1_B
    w_out[v2_B] = s2_out * w2_B

    # Carry over remaining nodes that have identical cluster definitions
    for node_A in tree_A.internal_nodes:
        if node_A.id not in (v1_A, v2_A):
            for node_B in tree_B.internal_nodes:
                if node_B.id not in (v1_B, v2_B) and node_B.body_indices == node_A.body_indices:
                    r_out[node_B.id] = np.copy(coords.r[node_A.id])
                    w_out[node_B.id] = np.copy(coords.w[node_A.id])
                    break

    return JacobiCoordinates(r=r_out, w=w_out)


def apply_pbm4_gauge_transition_jacobi(
    coords: JacobiCoordinates,
    tree_A: JacobiTree,
    tree_B: JacobiTree,
    state: ExtendedPhaseState,
    dt: float,
    G: float = 1.0,
    divisor: float = 480.0,
    apex: Optional["ScreenedApexPotential"] = None,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    mode: str = "global",
    corrector: Optional["SymplecticCorrector"] = None,
    dt_new: Optional[float] = None,
) -> JacobiCoordinates:
    """
    Applies the PBM4 symplectic gauge transition operator in Jacobi coordinates:
        T_{(A, h_A) -> (B, h_B)}^{(PBM4)}: (r^{(A)}, w^{(A)}) |--> (r^{(B)}, w^{(B)})

    Supports:
    - Non-commutative cotangent/homological sandwich (`mode="homological"` or `mode="cotangent_lift"`,
      or when an explicit `corrector` is provided):
          T_{(A, h_A) -> (B, h_B)}^{(PBM4)} = (C_{h_B}^{(B)})^{-1} o M_{A -> B} o C_{h_A}^{(A)}
    - Unipotent configuration-only ablation kick (`mode="global"`):
          Delta w = (M^{-1})^T Delta w_A - Delta w_B.
    """
    from ..integrators.correctors import CorrectorC4, compute_conormal_corrector_kick
    from ..integrators.processed import apply_gauge_connection

    h_old = dt
    h_new = dt if dt_new is None else dt_new

    if (
        (tree_A is tree_B or tree_A.get_leaf_pair_indices() == tree_B.get_leaf_pair_indices())
        and h_old == h_new
    ) or (h_old == 0.0 and h_new == 0.0):
        return coords

    if corrector is not None or mode not in ("global", "legacy"):
        corr = corrector if corrector is not None else CorrectorC4(G=G, mode=mode)
        st_A = jacobi_to_cartesian(
            coords=coords,
            tree=tree_A,
            masses=state.masses,
            center_of_mass=state.center_of_mass,
            total_momentum=state.total_momentum,
            t=state.t,
        )
        st_B = apply_gauge_connection(
            state=st_A,
            old_tree=tree_A,
            new_tree=tree_B,
            dt=h_old,
            corrector=corr,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
            dt_new=h_new,
        )
        return cartesian_to_jacobi(st_B, tree_B)

    dp_A = compute_conormal_corrector_kick(
        state, tree_A, h_old, G=G, inverse=False, divisor=divisor,
        apex=apex, supernode_bodies=supernode_bodies, supernode_clusters=supernode_clusters,
    )
    dp_B = compute_conormal_corrector_kick(
        state, tree_B, h_new, G=G, inverse=False, divisor=divisor,
        apex=apex, supernode_bodies=supernode_bodies, supernode_clusters=supernode_clusters,
    )

    delta_p = dp_A - dp_B

    # Project the momentum jump to Jacobi coordinates in tree_B
    st_kick = state.with_updates(p=delta_p)
    coords_kick = cartesian_to_jacobi(st_kick, tree_B)

    # Base Jacobi transformation from tree_A to tree_B
    coords_B = cartesian_to_jacobi(state, tree_B)

    w_new: Dict[int, np.ndarray] = {}
    for node_id, w_b in coords_B.w.items():
        dw = coords_kick.w.get(node_id)
        if dw is not None:
            w_new[node_id] = w_b + dw
        else:
            w_new[node_id] = w_b

    return JacobiCoordinates(r=coords_B.r, w=w_new)
