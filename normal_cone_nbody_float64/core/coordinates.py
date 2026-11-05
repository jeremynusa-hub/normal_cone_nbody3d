r"""
Symplectic Coordinate Transformations and Kinetic Cometric Decoupling.

Mathematical Monograph Alignment:
- Chapter 1: Section 1.4 (\label{sec:jacobi_trees}) — Canonical cotangent lift and cometric tensor.
- Chapter 2: Section 2.1, Equations (2.5)-(2.9) — Normal cone symplectic splitting and Rees limit.
- Chapter 3: Section 3.1 & 3.2, Equations (3.6)-(3.12) — Canonical symplectic diffeomorphism
  (q, p) <-> (r_v, w_v) and diagonal kinetic decomposition.
- Chapter 4: Section 4.1, Equations (4.13)-(4.15) — Barycentric cluster recursions and cometric tensor.
- Certified Lean 4 Alignment: PaperUnifiedFoundations.lean
  (omega0_skew, omega0At_nondeg, deformation_limit_zero).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple, Optional
import numpy as np

from .state import ExtendedPhaseState
from .jacobi import JacobiTree, JacobiTopologyNode


@dataclass(frozen=True)
class JacobiCoordinates:
    """
    Decoupled dynamic phase coordinates in the adapted Jacobi cotangent chart:
    - r: Dictionary mapping internal node ID -> relative position vector r_v in R^3.
    - w: Dictionary mapping internal node ID -> canonical momentum vector w_v in R^3.
    """
    r: Dict[int, np.ndarray]
    w: Dict[int, np.ndarray]

    def get_r(self, node_id: int) -> np.ndarray:
        """Retrieves relative position vector r_v."""
        return self.r[node_id]

    def get_w(self, node_id: int) -> np.ndarray:
        """Retrieves conjugate relative momentum w_v."""
        return self.w[node_id]


def cartesian_to_jacobi(
    state: ExtendedPhaseState,
    tree: JacobiTree,
) -> JacobiCoordinates:
    """
    Performs the canonical symplectic transformation from Cartesian phase space T*M
    to the adapted Jacobi cotangent chart T*Q[N]:
        (q_i, p_i)_{i=1}^N  |-->  (r_v, w_v)_{v in V_int(T)}.

    Formulas:
        r_v = X_L(v) - X_R(v) = (1/M_L sum_{i in L} m_i q_i) - (1/M_R sum_{j in R} m_j q_j)
        w_v = mu_v (P_L(v) / M_L - P_R(v) / M_R)
    """
    r_dict: Dict[int, np.ndarray] = {}
    w_dict: Dict[int, np.ndarray] = {}

    # Recursive helper computing (X_cluster, P_cluster) bottom-up
    def _compute_cluster_kinematics(node: JacobiTopologyNode) -> Tuple[np.ndarray, np.ndarray]:
        if node.is_leaf:
            assert node.body_index is not None
            idx = node.body_index
            x_leaf = np.copy(state.q[idx])
            p_leaf = np.copy(state.p[idx])
            return x_leaf, p_leaf

        assert node.left is not None and node.right is not None
        x_left, p_left = _compute_cluster_kinematics(node.left)
        x_right, p_right = _compute_cluster_kinematics(node.right)

        # Jacobi relative coordinates
        r_v = x_left - x_right
        # w_v = mu_v * (v_left - v_right) = mu_v * (P_left / M_left - P_right / M_right)
        w_v = node.mu * ((p_left / node.left.M) - (p_right / node.right.M))

        r_dict[node.id] = r_v
        w_dict[node.id] = w_v

        # Aggregate cluster kinematics
        x_cluster = node.c_left * x_left + node.c_right * x_right
        p_cluster = p_left + p_right
        return x_cluster, p_cluster

    _compute_cluster_kinematics(tree.root)
    return JacobiCoordinates(r=r_dict, w=w_dict)


def jacobi_to_cartesian(
    coords: JacobiCoordinates,
    tree: JacobiTree,
    masses: np.ndarray,
    t: float = 0.0,
    E: float = 0.0,
    E0: float = 0.0,
    center_of_mass: Optional[np.ndarray] = None,
    total_momentum: Optional[np.ndarray] = None,
) -> ExtendedPhaseState:
    """
    Performs the exact inverse symplectic transformation from the Jacobi chart back
    to physical Cartesian coordinates:
        (r_v, w_v)_{v in V_int(T)}  |-->  (q_i, p_i)_{i=1}^N.

    Formulas:
        X_left  = X_v + c_right * r_v
        X_right = X_v - c_left  * r_v
        P_left  = c_left  * P_v + w_v
        P_right = c_right * P_v - w_v
    """
    N = len(masses)
    q_out = np.zeros((N, 3), dtype=np.float64)
    p_out = np.zeros((N, 3), dtype=np.float64)

    r_cm = np.zeros(3, dtype=np.float64) if center_of_mass is None else np.asarray(center_of_mass, dtype=np.float64)
    p_tot = np.zeros(3, dtype=np.float64) if total_momentum is None else np.asarray(total_momentum, dtype=np.float64)

    # Top-down recursive reconstruction
    def _distribute_kinematics(
        node: JacobiTopologyNode,
        x_parent: np.ndarray,
        p_parent: np.ndarray,
    ) -> None:
        if node.is_leaf:
            assert node.body_index is not None
            idx = node.body_index
            q_out[idx] = x_parent
            p_out[idx] = p_parent
            return

        assert node.left is not None and node.right is not None
        r_v = coords.r[node.id]
        w_v = coords.w[node.id]

        x_left = x_parent + node.c_right * r_v
        x_right = x_parent - node.c_left * r_v

        p_left = node.c_left * p_parent + w_v
        p_right = node.c_right * p_parent - w_v

        _distribute_kinematics(node.left, x_left, p_left)
        _distribute_kinematics(node.right, x_right, p_right)

    _distribute_kinematics(tree.root, r_cm, p_tot)

    return ExtendedPhaseState(
        q=q_out,
        p=p_out,
        masses=masses,
        t=t,
        E=E,
        E0=E0,
    )


def compute_jacobi_kinetic_energy(
    coords: JacobiCoordinates,
    tree: JacobiTree,
) -> float:
    """
    Computes relative internal kinetic energy diagonalized over Jacobi nodes:
        T_rel = (1/2) * sum_{v in V_int(T)} (||w_v||^2 / mu_v).
    """
    t_kin = 0.0
    for node in tree.internal_nodes:
        w_v = coords.w[node.id]
        t_kin += 0.5 * float(np.dot(w_v, w_v)) / node.mu
    return float(t_kin)


def darboux_pairing_error(
    state: ExtendedPhaseState,
    coords: JacobiCoordinates,
    tree: JacobiTree,
) -> float:
    """
    Verifies exact symplectic pairing conservation:
        sum_{i=1}^N <q_i, p_i> == <R_cm, P_tot> + sum_{v in V_int(T)} <r_v, w_v>.
    Returns the absolute difference |Cartesian_pairing - Jacobi_pairing|.
    """
    cartesian_sum = float(np.sum(state.q * state.p))
    r_cm = state.center_of_mass
    p_tot = state.total_momentum
    cm_sum = float(np.dot(r_cm, p_tot))

    jacobi_sum = 0.0
    for node in tree.internal_nodes:
        jacobi_sum += float(np.dot(coords.r[node.id], coords.w[node.id]))

    return abs(cartesian_sum - (cm_sum + jacobi_sum))


def rees_deformation_retraction(
    coords: JacobiCoordinates,
    tree: JacobiTree,
    t_param: float,
) -> JacobiCoordinates:
    """
    Instantiates the graded Rees algebra deformation morphism:
        ev_t(r_v, w_v) = (t * r_v, w_v / t) for t > 0.
    Certified Lean 4 Theorem: deformation_limit_zero.
    As t -> 0, transverse distances scale linearly to zero, retracting onto the singular locus Z.
    """
    r_scaled: Dict[int, np.ndarray] = {}
    w_scaled: Dict[int, np.ndarray] = {}
    t_val = float(t_param)

    for node_id, r_vec in coords.r.items():
        r_scaled[node_id] = t_val * r_vec
        w_scaled[node_id] = np.copy(coords.w[node_id]) if t_val == 0.0 else coords.w[node_id] / t_val

    return JacobiCoordinates(r=r_scaled, w=w_scaled)
