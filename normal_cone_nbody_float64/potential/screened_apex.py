"""
Universal Screened Apex Potential C^{1,1} and Symplectic Kick Operator.

Theoretical Foundations:
- Monograph Reference: Section 4.2 (Universal C^{1,1} Regularity and Multipole Expansion),
  Equations (4.73)-(4.163).
- Universal C^{1,1} Lipschitz Regularity:
    ||nabla V_{apex}^{hier}|| <= M_{apex} < infty,
    ||Hess V_{apex}^{hier}|| <= L_{apex} < infty (strict zero softening d_0 = 0).
- Pure Quadrupole Asymptotics (u_v >= 2):
    V_{apex, v}(q) = - G / ||r_v||^3 * Q_v(s) + O(r_v^{-4}),
    where Q_v(s) = (1/2) [M_R (s^T Q_L s) + M_L (s^T Q_R s)].
- Ghost-Pole Eradication:
    When ||r_v|| -> 0, u_v <= 1.0 => Theta(u_v) = 0.0 => screened monopole vanishes,
    leaving only smooth inter-cluster Newtonian interactions without fictitious center-of-mass poles.
- Certified Formal Proofs:
    Lean 4 module `01_ApexUniversalLipschitzBound_Theorems.lean` (theorems `apex_universal_lipschitz_bound`,
    `apex_potential_bounded`).
"""

import math
from typing import Tuple, List, Optional, Set, Sequence, Union
import numpy as np

from ..core.constants import GRAVITATIONAL_CONSTANT
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree, JacobiTopologyNode
from .screening import eval_screening_theta, get_screening_mode_for_order


def _normalize_supernode_clusters(
    supernode_bodies: Optional[Union[Set[int], Sequence[Tuple[int, ...]], Set[Tuple[int, ...]]]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
) -> List[Tuple[int, ...]]:
    """
    Normalizes supernode specifications into a canonical list of disjoint sorted body-index tuples.
    Supports both explicit `supernode_clusters=[(0, 1, 2), (3, 4, 5)]` and legacy
    `supernode_bodies={0, 1, 2}` (or `supernode_bodies=[(0, 1, 2), (3, 4, 5)]`).
    """
    clusters: List[Tuple[int, ...]] = []
    if supernode_clusters is not None:
        for c in supernode_clusters:
            if c and len(c) >= 2:
                clusters.append(tuple(sorted(int(x) for x in c)))
        return clusters
    if supernode_bodies:
        first_elem = next(iter(supernode_bodies))
        if isinstance(first_elem, (int, np.integer)):
            t = tuple(sorted(int(x) for x in supernode_bodies))
            if len(t) >= 2:
                clusters.append(t)
        else:
            for c in supernode_bodies:
                if c and len(c) >= 2:
                    clusters.append(tuple(sorted(int(x) for x in c)))
    return clusters


def compute_quadrupole_tensor(
    q_cluster: np.ndarray,
    masses_cluster: np.ndarray,
) -> np.ndarray:
    """
    Computes the traceless quadrupole moment tensor Q of a cluster:
        Q = sum_{i in C} m_i (3 x_i x_i^T - ||x_i||^2 I_3),
    where x_i = q_i - X_cdm.
    """
    m_tot = float(np.sum(masses_cluster))
    if m_tot <= 0.0:
        return np.zeros((3, 3), dtype=np.float64)

    cdm = np.sum(masses_cluster[:, None] * q_cluster, axis=0) / m_tot
    x = q_cluster - cdm  # shape (k, 3)

    Q = np.zeros((3, 3), dtype=np.float64)
    I3 = np.eye(3, dtype=np.float64)

    for i in range(len(masses_cluster)):
        xi = x[i]
        xi_sq = float(np.dot(xi, xi))
        Q += masses_cluster[i] * (3.0 * np.outer(xi, xi) - xi_sq * I3)

    return Q


def compute_conormal_tidal_tensor(
    state: ExtendedPhaseState,
    leaf_pair: Tuple[int, int],
    G: float = GRAVITATIONAL_CONSTANT,
) -> np.ndarray:
    """
    Computes the conormal tidal apex tensor T_apex in Sym_0(3) of a Keplerian leaf:
        T_apex = sum_{k not in {a, b}} G * m_k * (3 * R_k (x) R_k / R_k^5 - I_3 / R_k^3),
    where R_k = q_k - R_{ab}^{cdm}.
    Guarantees tracelessness Tr(T_apex) == 0.0 to machine precision (< 1e-15).
    """
    ia, ib = leaf_pair
    ma = float(state.masses[ia])
    mb = float(state.masses[ib])
    qa = state.q[ia]
    qb = state.q[ib]
    m_ab = ma + mb
    R_ab = (ma * qa + mb * qb) / m_ab
    T = np.zeros((3, 3), dtype=np.float64)
    I3 = np.eye(3, dtype=np.float64)
    G_val = float(G)
    for k in range(state.N):
        if k == ia or k == ib:
            continue
        mk = float(state.masses[k])
        R_k = state.q[k] - R_ab
        R_k_norm = float(np.linalg.norm(R_k))
        if R_k_norm <= 0.0:
            continue
        R_k_inv3 = 1.0 / (R_k_norm ** 3)
        R_k_inv5 = R_k_inv3 / (R_k_norm ** 2)
        Tk = (G_val * mk) * (3.0 * R_k_inv5 * np.outer(R_k, R_k) - R_k_inv3 * I3)
        T += Tk
    return T


compute_normal_cone_tidal_tensor = compute_conormal_tidal_tensor



def compute_quadrupole_potential(
    Q_L: np.ndarray,
    Q_R: np.ndarray,
    M_L: float,
    M_R: float,
    r_vec: np.ndarray,
    G: float = GRAVITATIONAL_CONSTANT,
) -> float:
    """
    Evaluates the asymptotic quadrupole interaction potential:
        V_quad = - G / ||r_v||^3 * Q_v(s),
    where Q_v(s) = (1/2) [M_R (s^T Q_L s) + M_L (s^T Q_R s)].
    """
    r = float(np.linalg.norm(r_vec))
    if r < 1e-15:
        return 0.0

    s = r_vec / r
    s_Q_L_s = float(s @ Q_L @ s)
    s_Q_R_s = float(s @ Q_R @ s)
    Q_v = 0.5 * (M_R * s_Q_L_s + M_L * s_Q_R_s)

    return float(- G * Q_v / (r**3))


class ScreenedApexPotential:
    """
    Evaluator of the regularized C^{1,1} hierarchical apex potential and kick propagator.
    """
    def __init__(
        self,
        G: float = GRAVITATIONAL_CONSTANT,
        screening_mode: Optional[str] = None,
        order: Optional[int] = None,
    ):
        self.G = float(G)
        if screening_mode is not None:
            self.screening_mode = screening_mode
        elif order is not None:
            self.screening_mode = get_screening_mode_for_order(order)
        else:
            self.screening_mode = "c6"

    def compute_screened_composite_potential(
        self,
        node: JacobiTopologyNode,
        state: ExtendedPhaseState,
    ) -> Tuple[float, np.ndarray]:
        """
        Computes the screened monopole potential and its relative gradient for a composite node v:
            V_v^screened(r_v) = - (G M_L M_R / ||r_v||) * Theta(||r_v|| / (bar_rho_L + bar_rho_R)).

        Mathematical Scope Note:
            This helper routine evaluates the partial derivative with respect to the inter-cluster
            relative vector r_v (partial V_v / partial r_v), keeping hyper-radii fixed.
            In contrast, compute_apex_forces evaluates the complete physical gradient on Cartesian
            coordinates, including internal hyper-radius deformation derivatives (nabla_{q_i} bar_rho_L)
            as established by Fulton-MacPherson conormal blow-up (Chapter 4, Theorem 4.1).

        Returns:
            (V_screened, grad_r): Potential scalar and relative gradient vector in R^3.
        """
        if node.is_leaf or node.left is None or node.right is None:
            return 0.0, np.zeros(3, dtype=np.float64)

        # Barycenters and mass-normalized hyper-radii
        idx_L = list(node.left.body_indices)
        idx_R = list(node.right.body_indices)

        m_L = state.masses[idx_L]
        m_R = state.masses[idx_R]
        q_L = state.q[idx_L]
        q_R = state.q[idx_R]

        M_L = float(np.sum(m_L))
        M_R = float(np.sum(m_R))

        cdm_L = np.sum(m_L[:, None] * q_L, axis=0) / M_L
        cdm_R = np.sum(m_R[:, None] * q_R, axis=0) / M_R

        r_vec = cdm_L - cdm_R
        r_norm = float(np.linalg.norm(r_vec))

        if r_norm < 1e-15:
            return 0.0, np.zeros(3, dtype=np.float64)

        # Mass-normalized hyper-radii bar_rho in [L]
        rho_L = float(np.sqrt(max(0.0, np.sum(m_L * np.sum((q_L - cdm_L)**2, axis=1)))))
        rho_R = float(np.sqrt(max(0.0, np.sum(m_R * np.sum((q_R - cdm_R)**2, axis=1)))))

        bar_rho_L = rho_L / math.sqrt(M_L) if M_L > 0.0 else 0.0
        bar_rho_R = rho_R / math.sqrt(M_R) if M_R > 0.0 else 0.0

        R_scale = max(bar_rho_L + bar_rho_R, 1e-15)
        u_v = r_norm / R_scale

        theta, d_theta = eval_screening_theta(u_v, mode=self.screening_mode)
        if theta == 0.0 and d_theta == 0.0:
            return 0.0, np.zeros(3, dtype=np.float64)

        GM = self.G * M_L * M_R
        V_screened = - (GM / r_norm) * theta

        # Analytical gradient with respect to relative vector r_vec:
        # d/dr [ - GM/r * Theta(r/R) ] = GM/r^2 * Theta - GM/(r*R) * Theta'
        d_val = (GM / (r_norm * r_norm)) * theta - (GM / (r_norm * R_scale)) * d_theta
        grad_r = (d_val / r_norm) * r_vec

        return V_screened, grad_r

    def compute_apex_forces(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        supernode_bodies: Optional[Set[int]] = None,
        include_screened_monopole: bool = False,
        mask_leaf_pairs: Optional[bool] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> Tuple[np.ndarray, float]:
        """
        Computes exact Cartesian apex forces F_i^apex = - nabla_{q_i} V_{apex}^{hier} acting on each body:
            V_{apex}^{hier} = V_{Newton} - sum_{LeafPairs} V_{Kepler} - sum_{Comp} V_{screened} - sum_{Super} V_{int}.

        When self.screening_mode == 'pure_cartesian' or mask_leaf_pairs is False, computes
        the exact unmasked classical Newtonian gravitational forces without conormal splitting.
        """
        N = state.N
        if N <= 1:
            return np.zeros((N, 3), dtype=np.float64), 0.0

        is_pure_cartesian = (self.screening_mode == "pure_cartesian")
        eff_mask_leaves = (not is_pure_cartesian) if mask_leaf_pairs is None else bool(mask_leaf_pairs)
        eff_screened_monopole = False if is_pure_cartesian else include_screened_monopole

        clusters = _normalize_supernode_clusters(supernode_bodies, supernode_clusters)
        cluster_sets = [set(c) for c in clusters]

        # Vectorized pairwise relative displacement: diff[i, j] = q[j] - q[i]
        diff = state.q[None, :, :] - state.q[:, None, :]  # shape (N, N, 3)
        dist_sq = np.sum(diff**2, axis=-1)                # shape (N, N)
        np.fill_diagonal(dist_sq, 1.0)
        dist = np.sqrt(dist_sq)
        with np.errstate(divide='ignore', invalid='ignore'):
            inv_dist3 = np.where(dist_sq > 0, 1.0 / (dist_sq * dist), 0.0)
            inv_dist = np.where(dist > 0, 1.0 / dist, 0.0)
        np.fill_diagonal(inv_dist3, 0.0)
        np.fill_diagonal(inv_dist, 0.0)

        # Mask leaf pairs (handled analytically by Stumpff-Sundman Kepler solver in conormal mode)
        if eff_mask_leaves:
            for i, j in tree.get_leaf_pair_indices():
                inv_dist3[i, j] = 0.0
                inv_dist3[j, i] = 0.0
                inv_dist[i, j] = 0.0
                inv_dist[j, i] = 0.0

        # Mask intra-supernode pairs if specified (handled synchronously by supernode solver)
        for cluster in clusters:
            for idx_a in range(len(cluster)):
                ia = cluster[idx_a]
                for idx_b in range(idx_a + 1, len(cluster)):
                    ib = cluster[idx_b]
                    inv_dist3[ia, ib] = 0.0
                    inv_dist3[ib, ia] = 0.0
                    inv_dist[ia, ib] = 0.0
                    inv_dist[ib, ia] = 0.0

        # Pairwise mass coupling matrix G m_i m_j
        m_matrix = self.G * (state.masses[:, None] * state.masses[None, :])

        # Base inter-cluster Newtonian forces (sequential column reduction matching pair order)
        f_matrix = m_matrix * inv_dist3
        term_matrix = f_matrix[:, :, None] * diff
        forces = np.zeros((N, 3), dtype=np.float64)
        for j in range(N):
            forces += term_matrix[:, j, :]

        # Subtract screened monopole contributions from free composite internal nodes
        # V_{apex} = V_{inter} + sum_{v in V_{comp}^{free}} Theta(u_v) (G M_L M_R / r_v)
        pot_matrix = m_matrix * inv_dist
        V_apex = 0.0
        for i in range(N):
            for j in range(i + 1, N):
                if pot_matrix[i, j] != 0.0:
                    V_apex -= float(pot_matrix[i, j])

        if eff_screened_monopole:
            for node in tree.internal_nodes:
                if node.is_composite_node():
                    if node.is_leaf or node.left is None or node.right is None:
                        continue
                    if cluster_sets and any(set(node.body_indices).issubset(c_set) for c_set in cluster_sets):
                        continue

                    idx_L = list(node.left.body_indices)
                    idx_R = list(node.right.body_indices)

                    m_L = state.masses[idx_L]
                    m_R = state.masses[idx_R]
                    q_L = state.q[idx_L]
                    q_R = state.q[idx_R]

                    M_L = float(np.sum(m_L))
                    M_R = float(np.sum(m_R))

                    cdm_L = np.sum(m_L[:, None] * q_L, axis=0) / M_L
                    cdm_R = np.sum(m_R[:, None] * q_R, axis=0) / M_R

                    r_vec = cdm_L - cdm_R
                    r_norm = float(np.linalg.norm(r_vec))

                    if r_norm < 1e-15:
                        continue

                    # Mass-normalized hyper-radii bar_rho in [L]
                    rho_L = float(np.sqrt(max(0.0, np.sum(m_L * np.sum((q_L - cdm_L)**2, axis=1)))))
                    rho_R = float(np.sqrt(max(0.0, np.sum(m_R * np.sum((q_R - cdm_R)**2, axis=1)))))

                    bar_rho_L = rho_L / math.sqrt(M_L) if M_L > 0.0 else 0.0
                    bar_rho_R = rho_R / math.sqrt(M_R) if M_R > 0.0 else 0.0

                    R_scale = max(bar_rho_L + bar_rho_R, 1e-15)
                    u_v = r_norm / R_scale

                    theta, d_theta = eval_screening_theta(u_v, mode=self.screening_mode)
                    if theta == 0.0 and d_theta == 0.0:
                        continue

                    GM = self.G * M_L * M_R
                    W_v = (GM / r_norm) * theta
                    V_apex += W_v

                    # Exact analytical gradient of W_v:
                    # grad_W_r = d W_v / d r_vec
                    # g_W_rho = d W_v / d R_scale
                    grad_W_r = ((-GM / (r_norm * r_norm)) * theta + (GM / (r_norm * R_scale)) * d_theta) * (r_vec / r_norm)
                    g_W_rho = -(GM / (R_scale * R_scale)) * d_theta

                    # Force contributions: F_k = - grad_{q_k} W_v
                    # For i in L: d r_vec / d q_i = m_i / M_L, d bar_rho_L / d q_i = m_i (q_i - cdm_L) / (M_L * bar_rho_L)
                    for i in idx_L:
                        F_i_W = - (state.masses[i] / M_L) * grad_W_r
                        if bar_rho_L > 1e-14:
                            F_i_W -= g_W_rho * (state.masses[i] * (state.q[i] - cdm_L)) / (M_L * bar_rho_L)
                        forces[i] += F_i_W

                    # For j in R: d r_vec / d q_j = - m_j / M_R, d bar_rho_R / d q_j = m_j (q_j - cdm_R) / (M_R * bar_rho_R)
                    for j in idx_R:
                        F_j_W = + (state.masses[j] / M_R) * grad_W_r
                        if bar_rho_R > 1e-14:
                            F_j_W -= g_W_rho * (state.masses[j] * (state.q[j] - cdm_R)) / (M_R * bar_rho_R)
                        forces[j] += F_j_W

        return forces, V_apex

    def step_apex_kick(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        dt: float,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """
        Applies a symplectic kick step p_i <- p_i + F_i^apex * dt.

        Analytical Monopole Cancellation Identity:
            In the conormal Hamiltonian splitting (Chapter 4, Eq. 4.41-4.48):
                H = H_base + H_leaf + V_apex^hier + H_super^sync.
            The physical inter-cluster potential satisfies:
                V_inter = sum_{v in V_int \\ V_leaf} W_v(q).
            When decomposing W_v into screened monopole and apex remainder:
                W_v = (1 - Theta) V_v^raw + Theta (V_v^raw - G M_L M_R / r_v) + Theta (G M_L M_R / r_v).
            The regularized apex potential V_apex^hier subtracts + Theta G M_L M_R / r_v.
            In the operator-splitting flow:
                V_kick = V_apex^hier + V_kepler^screened
                       = (V_inter - Theta G M_L M_R / r_v) + (Theta G M_L M_R / r_v)
                       = V_inter.
            The + and - screened monopoles cancel analytically to machine zero (+- Theta G M_L M_R / r_v = 0).
            Therefore, calling compute_apex_forces with include_screened_monopole=False directly
            evaluates the exact physical kick force of V_inter, eliminating spurious virtual kicks
            without uncompensated barycentric reactions.
        """
        forces, _ = self.compute_apex_forces(
            state,
            tree,
            supernode_bodies=supernode_bodies,
            include_screened_monopole=False,
            supernode_clusters=supernode_clusters,
        )
        new_p = state.p + dt * forces
        return state.with_updates(p=new_p)


def compute_intercluster_forces(
    q: np.ndarray,
    masses: np.ndarray,
    leaf_pair_set: Set[Tuple[int, int]],
    supernode_bodies: Optional[Set[int]] = None,
    G: float = 1.0,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
) -> np.ndarray:
    """
    Computes exact inter-cluster Newtonian forces acting on Cartesian positions:
        F_i = - nabla_{q_i} V_inter,
    where V_inter excludes pairs belonging to terminal leaves and synchronous supernodes.
    Guarantees sum F_i == 0 and sum q_i x F_i == 0 to IEEE 754 precision (< 1e-15).
    """
    N = len(masses)
    clusters = _normalize_supernode_clusters(supernode_bodies, supernode_clusters)
    diff = q[None, :, :] - q[:, None, :]
    dist_sq = np.sum(diff**2, axis=-1)
    np.fill_diagonal(dist_sq, 1.0)
    dist = np.sqrt(dist_sq)

    with np.errstate(divide='ignore', invalid='ignore'):
        inv_dist3 = np.where(dist_sq > 0, 1.0 / (dist_sq * dist), 0.0)
    np.fill_diagonal(inv_dist3, 0.0)

    for i, j in leaf_pair_set:
        inv_dist3[i, j] = 0.0
        inv_dist3[j, i] = 0.0

    for cluster in clusters:
        for idx_a in range(len(cluster)):
            ia = cluster[idx_a]
            for idx_b in range(idx_a + 1, len(cluster)):
                ib = cluster[idx_b]
                inv_dist3[ia, ib] = 0.0
                inv_dist3[ib, ia] = 0.0

    m_matrix = G * (masses[:, None] * masses[None, :])
    f_matrix = m_matrix * inv_dist3
    return np.sum(f_matrix[:, :, None] * diff, axis=1)
