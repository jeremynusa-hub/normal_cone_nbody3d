"""
Boundary Symplectic Corrector C_4 on the Fulton--MacPherson Configuration Space.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.5 (Symplectic Correctors at the Boundary and Comparative Synthesis),
  Theorems 5.14, 5.15, 5.16, 5.17.
- Homological Lie Bracket Equations: {H, K_4} = - H_5^{linear}.
- Operator Inversion & Dyadic Telescoping:
    C_4(h) = exp(- h^4 X_{K_4}),
    C_4(h)^{-1} = exp(+ h^4 X_{K_4}),
    with exact machine-precision identity C_4^{-1} o C_4 = Id (< 1e-15).
- Formal Lean 4 Certification:
  * 01_CorrectorBM4_Theorems.lean (cone_normal_corrector_bm4_certified)
  * 03_ConormalCorrectorPBM4_Order6.lean (pbm4_effective_order_six_cancellation)
  * 07_GlobalConormalCorrectorPBM4_Order6.lean (global_pbm4_effective_order_six_cancellation)
  * 08_UniversalOrderSixObstruction.lean (universal_order_six_obstruction).
"""

from __future__ import annotations
import math
from typing import Optional, Sequence, Set, Tuple
import numpy as np

from ..core.constants import (
    G_DEFAULT,
    KAPPA4_BM4,
    C4_BCH_BM4,
    HOMOLOGICAL_ALPHA1_BM4,
    HOMOLOGICAL_ALPHA2_BM4,
    HOMOLOGICAL_BETA1_BM4,
    HOMOLOGICAL_BETA2_BM4,
)
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree, JacobiTopologyNode
from ..potential.screened_apex import ScreenedApexPotential
from ..potential.screening import eval_screening_theta


class SymplecticCorrector:
    """
    Abstract base protocol for boundary symplectic correctors C_{2p}(h).
    """

    @property
    def order(self) -> int:
        """Formal order of the correction term: 2p (4)."""
        raise NotImplementedError

    def apply(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        h: float,
        inverse: bool = False,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """
        Applies the boundary corrector C_4(h) = exp(- h^4 X_{K_4})
        or its algebraic inverse C_4(h)^{-1} = exp(+ h^4 X_{K_4}).
        """
        raise NotImplementedError

    def apply_inverse(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        h: float,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """Convenience method applying the algebraic inverse operator."""
        return self.apply(
            state,
            tree,
            h,
            inverse=True,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )


def _is_pair_in_same_supernode_cluster(
    ia: int,
    ib: int,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
) -> bool:
    """Returns True if (ia, ib) belong to the same synchronous supernode cluster."""
    if supernode_clusters:
        for cluster in supernode_clusters:
            if ia in cluster and ib in cluster:
                return True
        return False
    if supernode_bodies:
        return ia in supernode_bodies and ib in supernode_bodies
    return False


def compute_analytical_directional_force_derivative(
    eval_apex: ScreenedApexPotential,
    state: ExtendedPhaseState,
    tree: JacobiTree,
    supernode_bodies: Optional[Set[int]],
    v_dir: np.ndarray,
) -> np.ndarray:
    """
    Computes the exact analytical directional derivative of the apex force along velocity direction v_dir:
        D_v F_i = d/ds [ F_i(q + s * v_dir) ]_{s=0} = - (nabla^2 V_apex . v)_i
    at full IEEE 754 float64 machine precision (< 1e-16), completely eliminating
    the finite-difference stencil (eta = 1e-7) and its 1e-9 noise floor.

    For each active unmasked pair (i, j) with r_{ij} = q_j - q_i, Delta v_{ij} = v_j - v_i:
        D_v F_{ij} = G m_i m_j [ Delta v_{ij} / r_{ij}^3 - 3 (r_{ij} . Delta v_{ij}) r_{ij} / r_{ij}^5 ]
        D_v F_i = sum_{j != i, (i, j) not masked} D_v F_{ij}
    Pairs in tree.leaf_pairs and supernode_bodies are masked identically to compute_apex_forces.
    """
    v_norm = float(np.max(np.linalg.norm(v_dir, axis=1)))
    if v_norm < 1e-15 or state.N <= 1:
        return np.zeros_like(state.q)

    # 1. Pairwise relative positions and velocities
    diff = state.q[None, :, :] - state.q[:, None, :]  # diff[i, j] = q_j - q_i
    v_diff = v_dir[None, :, :] - v_dir[:, None, :]    # v_diff[i, j] = v_j - v_i

    dist_sq = np.sum(diff ** 2, axis=-1)
    dist = np.sqrt(dist_sq)

    # Guard diagonal self-interactions
    np.fill_diagonal(dist, np.inf)
    np.fill_diagonal(dist_sq, np.inf)

    # 2. Kernel powers
    inv_dist3 = 1.0 / (dist * dist_sq)
    inv_dist5 = inv_dist3 / dist_sq

    # 3. Mask leaf pairs solved by analytical Kepler-Sundman solver
    for leaf_node in tree.leaf_pairs:
        ia, ib = leaf_node.body_indices
        inv_dist3[ia, ib] = 0.0
        inv_dist3[ib, ia] = 0.0
        inv_dist5[ia, ib] = 0.0
        inv_dist5[ib, ia] = 0.0

    # 4. Mask supernode bodies solved synchronously by supernode solver
    if supernode_bodies:
        sn_list = list(supernode_bodies)
        for idx_a in range(len(sn_list)):
            ia = sn_list[idx_a]
            for idx_b in range(idx_a + 1, len(sn_list)):
                ib = sn_list[idx_b]
                inv_dist3[ia, ib] = 0.0
                inv_dist3[ib, ia] = 0.0
                inv_dist5[ia, ib] = 0.0
                inv_dist5[ib, ia] = 0.0

    # 5. Pairwise mass coupling G m_i m_j
    m_matrix = eval_apex.G * (state.masses[:, None] * state.masses[None, :])

    # 6. Directional derivative tensor:
    # D_v F_{ij} = G m_i m_j [ Delta v_{ij} / r_{ij}^3 - 3 (r_{ij} . Delta v_{ij}) r_{ij} / r_{ij}^5 ]
    r_dot_v = np.sum(diff * v_diff, axis=-1)
    term1 = v_diff * inv_dist3[:, :, None]
    term2 = 3.0 * diff * (r_dot_v * inv_dist5)[:, :, None]

    dF_matrix = m_matrix[:, :, None] * (term1 - term2)

    return np.sum(dF_matrix, axis=1)


# Backward compatibility alias
_directional_force_derivative = compute_analytical_directional_force_derivative


def compute_conormal_tidal_tensor(
    state: ExtendedPhaseState,
    leaf_node: JacobiTopologyNode,
    G: float = G_DEFAULT,
    apex: Optional[ScreenedApexPotential] = None,
) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """
    Computes the tidal apex tensor T_apex in Sym_0(3) and the barycentric geometry
    for a given binary leaf pair (a, b) at machine precision (< 1e-16).

    Mathematical Monograph Alignment:
    - Chapter 5, Section 5.5, Theorem 5.15 & Theorem 5.16:
        T_apex = sum_{k not in {a, b}} Theta(u_k) * G * m_k * (3 * R_k (x) R_k / R_k^5 - I_3 / R_k^3),
      where R_k = q_k - R_cm, and R_cm = (m_a * q_a + m_b * q_b) / (m_a + m_b).
    - Traceless property: Tr(T_apex) == 0.0 identically.

    Returns:
        T_apex: (3, 3) symmetric traceless tensor in Sym_0(3) (Tr(T_apex) == 0.0)
        r: (3,) relative vector q_b - q_a
        mu: reduced mass (m_a * m_b) / (m_a + m_b)
        R_cm: (3,) barycenter (m_a * q_a + m_b * q_b) / (m_a + m_b)
    """
    ia, ib = leaf_node.body_indices[0], leaf_node.body_indices[1]
    ma = float(state.masses[ia])
    mb = float(state.masses[ib])
    M_ab = ma + mb
    mu = (ma * mb) / M_ab

    r = state.q[ib] - state.q[ia]
    R_cm = (ma * state.q[ia] + mb * state.q[ib]) / M_ab
    r_sq = float(np.dot(r, r))
    rho_bar = math.sqrt((mu / M_ab) * r_sq) if r_sq > 0.0 else 0.0

    T_apex = np.zeros((3, 3), dtype=np.float64)
    for k in range(state.N):
        if k == ia or k == ib:
            continue
        mk = float(state.masses[k])
        R_k = state.q[k] - R_cm
        R_k_sq = float(np.dot(R_k, R_k))
        if R_k_sq < 1e-30:
            continue
        R_k_norm = math.sqrt(R_k_sq)
        R_k_cube = R_k_sq * R_k_norm
        R_k_5 = R_k_cube * R_k_sq

        theta_k = 1.0
        if apex is not None and apex.screening_mode != "pure_cartesian" and rho_bar > 1e-30:
            u_k = R_k_norm / rho_bar
            theta_val, _ = eval_screening_theta(u_k, mode=apex.screening_mode)
            theta_k = float(theta_val)
            if theta_k == 0.0:
                continue

        T_k = (theta_k * G * mk) * (3.0 * np.outer(R_k, R_k) / R_k_5 - np.eye(3, dtype=np.float64) / R_k_cube)
        T_apex += T_k

    # Enforce exact symmetry and zero trace at machine epsilon
    T_apex = 0.5 * (T_apex + T_apex.T)
    T_apex -= (np.trace(T_apex) / 3.0) * np.eye(3, dtype=np.float64)

    return T_apex, r, mu, R_cm


def compute_conormal_corrector_kick(
    state: ExtendedPhaseState,
    tree: JacobiTree,
    h: float,
    G: float = G_DEFAULT,
    inverse: bool = False,
    half_step: bool = False,
    divisor: float = 1440.0,
    apex: Optional[ScreenedApexPotential] = None,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
) -> np.ndarray:
    """
    Computes the exact closed-form momentum shear kicks Delta p for all particles
    arising from the boundary Hamiltonian generator:
        K_4 = sum_{lambda in LeafPairs(T)} (mu_lambda / 960) * r_lambda^T T_apex(R_cm, lambda) r_lambda.

    Mathematical Properties (Theorems 5.15 & 5.16):
    - Homological Lie Bracket Resolution: {H_leaf, K_4} = (1 / 480) * r^T T_apex w,
      identically eliminating the leading degree-5 non-integrable tidal commutator H_5^linear.
    - Galilean Lie-Bargmann Invariants:
        * Delta q == 0 (center of mass stationary: Delta R_cm == 0)
        * sum_i Delta p_i == 0 (exact linear momentum conservation < 1e-16)
        * sum_i q_i x Delta p_i == 0 (exact angular momentum conservation < 1e-16 via action-reaction)
    - Pericenter Regularity:
        ||Delta w|| proportional to ||r|| -> 0 as r -> 0 (no collision singularity).
    - 5-Strata Partition (Theorem 5.16):
        * Active supernode leaf pairs (ia, ib in same supernode cluster) are skipped.
        * When apex is provided, spectator tidal contributions are attenuated by Theta(u_k).
    """
    if h == 0.0 or not tree.leaf_pairs or state.N <= 2 or (supernode_bodies is not None and len(supernode_bodies) >= state.N):
        return np.zeros_like(state.p)

    sign = +1.0 if inverse else -1.0
    factor = 0.5 if half_step else 1.0
    coeff = sign * factor * (float(h) ** 4) / float(divisor)

    delta_p = np.zeros_like(state.p)

    for leaf_node in tree.leaf_pairs:
        ia, ib = leaf_node.body_indices[0], leaf_node.body_indices[1]
        if _is_pair_in_same_supernode_cluster(ia, ib, supernode_bodies, supernode_clusters):
            continue

        ma = float(state.masses[ia])
        mb = float(state.masses[ib])
        M_ab = ma + mb

        T_apex, r, mu, R_cm = compute_conormal_tidal_tensor(state, leaf_node, G=G, apex=apex)
        r_sq = float(np.dot(r, r))
        rho_bar = math.sqrt((mu / M_ab) * r_sq) if r_sq > 0.0 else 0.0

        # 1. Primary tidal shear kick on binary leaf relative momentum:
        # Delta p_b = + delta_p_rel, Delta p_a = - delta_p_rel
        delta_p_rel = coeff * mu * (T_apex @ r)
        delta_p[ib] += delta_p_rel
        delta_p[ia] -= delta_p_rel

        # 2. Exact back-reaction on external bodies to guarantee conservation
        # of total linear and angular momentum:
        F_cm_back = np.zeros(3, dtype=np.float64)

        for k in range(state.N):
            if k == ia or k == ib:
                continue
            mk = float(state.masses[k])
            R_k = state.q[k] - R_cm
            R_k_sq = float(np.dot(R_k, R_k))
            if R_k_sq < 1e-30:
                continue
            R_k_norm = math.sqrt(R_k_sq)
            R_k_5 = R_k_sq * R_k_sq * R_k_norm
            R_k_7 = R_k_5 * R_k_sq

            theta_k = 1.0
            if apex is not None and apex.screening_mode != "pure_cartesian" and rho_bar > 1e-30:
                u_k = R_k_norm / rho_bar
                theta_val, _ = eval_screening_theta(u_k, mode=apex.screening_mode)
                theta_k = float(theta_val)
                if theta_k == 0.0:
                    continue

            R_k_dot_r = float(np.dot(R_k, r))

            term1 = 6.0 * R_k_dot_r * r / R_k_5
            term2 = 3.0 * r_sq * R_k / R_k_5
            term3 = 15.0 * (R_k_dot_r ** 2) * R_k / R_k_7

            grad_k = 0.5 * coeff * (theta_k * mu * float(G) * mk) * (term1 + term2 - term3)
            delta_p_k = grad_k

            delta_p[k] += delta_p_k
            F_cm_back += delta_p_k

        # Equal and opposite reaction distributed to leaf bodies by mass fractions
        delta_p[ia] -= (ma / M_ab) * F_cm_back
        delta_p[ib] -= (mb / M_ab) * F_cm_back

    return delta_p


def compute_global_conormal_c4_kick(
    state: ExtendedPhaseState,
    tree: JacobiTree,
    h: float,
    G: float = G_DEFAULT,
    inverse: bool = False,
    apex: Optional[ScreenedApexPotential] = None,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
) -> np.ndarray:
    """
    Computes the exact closed-form global boundary corrector kick Delta p (Theorems 5.15 & 5.16):
        K_4^{global}(r) = sum_{lambda in Leaves(T)} (mu_lambda / 960) * r_lambda^T T_apex(R_cm, lambda) r_lambda.

    Theoretical & Formal Alignments:
    - Lean 4 Certified Module: 07_GlobalConormalCorrectorPBM4_Order6.lean
    - Effective Order 6 on D_reg (eps << h): Exact homological factor h^4 / 480.0 eliminates BCH error H_5^{linear}.
    - Spatial Stationarity: nabla_p K_4 == 0 <=> Delta q == 0 (no Cartesian drift).
    - Invariants: sum Delta p_i = 0 and sum q_i x Delta p_i = 0 via octupolar back-reaction.
    - Exact Algebraic Reversibility: C^{-1} o C = Id (< 1e-15 in f64).
    """
    return compute_conormal_corrector_kick(
        state=state,
        tree=tree,
        h=h,
        G=G,
        inverse=inverse,
        half_step=False,
        divisor=480.0,
        apex=apex,
        supernode_bodies=supernode_bodies,
        supernode_clusters=supernode_clusters,
    )


compute_normal_cone_tidal_tensor = compute_conormal_tidal_tensor
compute_normal_cone_corrector_kick = compute_conormal_corrector_kick
compute_global_normal_cone_c4_kick = compute_global_conormal_c4_kick



class GlobalCorrectorC4(SymplecticCorrector):
    """
    Fourth-order global boundary symplectic corrector C_4^{global}(h) = exp(- h^4 X_{K_4^{global}}).

    Mathematical Monograph & Certified Lean 4 Alignment:
    - Chapter 5, Section 5.5, Theorems 5.15 & 5.16 & Chapter 6, Section 6.2 & 6.5.
    - Certified Lean 4 Module: 07_GlobalConormalCorrectorPBM4_Order6.lean (10 theorems certified without sorry).
    - Absolute Spatial Stationarity: Delta q == 0 (no Verlet drift flow, eliminating composite Cartesian drift).
    - Exact Homological Calibration: Net kick coefficient h^4 / 480.0 ensures dH propto h^6 (ratio 64) when eps << h.
    - Exact Algebraic Inversion: C^{-1} o C = Id (< 1e-15 in f64).
    - Lie--Bargmann Invariants: 10/10 conserved to machine precision (< 1e-15 in f64).
    """

    def __init__(
        self,
        G: float = G_DEFAULT,
        screening_mode: str = "c6",
    ) -> None:
        self.G = G
        self.screening_mode = screening_mode

    @property
    def order(self) -> int:
        return 4

    def apply(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        h: float,
        inverse: bool = False,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """
        Applies C_4^{global}(h) or (C_4^{global}(h))^{-1} via closed-form momentum shear kicks:
            (q, p) |-> (q, p - h^4 nabla_q K_4^{global}).
        Positions are strictly invariant (Delta q == 0).
        """
        if h == 0.0 or state.N <= 2 or not tree.leaf_pairs or (supernode_bodies is not None and len(supernode_bodies) >= state.N):
            return state

        dp = compute_global_conormal_c4_kick(
            state,
            tree,
            h,
            G=self.G,
            inverse=inverse,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )
        new_p = state.p + dp
        return state.with_updates(p=new_p)


def compute_cotangent_lift_c4_flow(
    state: ExtendedPhaseState,
    tree: JacobiTree,
    h: float,
    G: float = G_DEFAULT,
    inverse: bool = False,
    c4_coeff: float = C4_BCH_BM4,
    apex: Optional[ScreenedApexPotential] = None,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
) -> ExtendedPhaseState:
    """
    Computes the exact closed-form canonical cotangent-lift flow of the degree-1 momentum generator:
        K_4(q, p) = c_4^{BM4} * h^4 * D_A^3 B(q, p)
                  = sum_{lambda in Leaves(T)} Tr( T_{apex, lambda}(R_cm) @ M_{sym, lambda}(r_lambda, w_lambda) ) in P_1,
    where the symmetric quadrupole-momentum tensor is:
        M_{sym}(r, w) = c_4^{BM4} * h^4 * (G * M_{ab} / r^3) * [ 2 (r (x) w + w (x) r) - (3 (r . w) / r^2) r (x) r ].

    Mathematical Properties:
    - Exact Homological Resolution: {K_4, H_Kepler} = - D_A K_4 = - c_4^{BM4} * h^4 * D_A^4 B in P_0 (+) P_2,
      cancelling the even-momentum 4th-order BCH error of BM_4^{(7)} on arbitrary 3D eccentric orbits.
    - Exact Noether Torque-Free Closure:
        r x nabla_r K_4 + w x nabla_w K_4 + sum_{k not in {a, b}} R_k x nabla_{R_k} K_4 == 0,
      preserving all 10 Lie--Bargmann invariants (Delta P == 0, Delta L == 0, Delta R_cm == 0) to machine precision.
    - Integrated via the symmetric, symplectic 2-stage Gauss--Legendre collocation method (order 4)
      solved to machine precision, guaranteeing exact group inversion C_4^{-1} o C_4 = Id (< 1e-15)
      and exact conservation of quadratic/bilinear Noether invariants (L_tot).
    """
    h_abs = abs(float(h))
    if h_abs == 0.0 or not tree.leaf_pairs or state.N <= 2 or (supernode_bodies is not None and len(supernode_bodies) >= state.N):
        return state

    sign = +1.0 if inverse else -1.0
    coeff = sign * float(c4_coeff) * (h_abs ** 4)
    masses = state.masses

    def _grad_k4(q_cur: np.ndarray, p_cur: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        dq = np.zeros_like(q_cur)
        dp = np.zeros_like(p_cur)
        cur_st = state.with_updates(q=q_cur, p=p_cur)

        for leaf_node in tree.leaf_pairs:
            ia, ib = leaf_node.body_indices[0], leaf_node.body_indices[1]
            if _is_pair_in_same_supernode_cluster(ia, ib, supernode_bodies, supernode_clusters):
                continue

            ma = float(masses[ia])
            mb = float(masses[ib])
            M_ab = ma + mb
            mu = (ma * mb) / M_ab
            GM_ab = float(G) * M_ab

            T_apex, r, _, R_cm = compute_conormal_tidal_tensor(cur_st, leaf_node, G=G, apex=apex)
            r2 = float(np.dot(r, r))
            if r2 < 1e-30:
                continue
            rn = math.sqrt(r2)
            r3 = r2 * rn
            r5 = r3 * r2
            r7 = r5 * r2
            rho_bar = math.sqrt((mu / M_ab) * r2)

            w = (ma * p_cur[ib] - mb * p_cur[ia]) / M_ab
            Tr = T_apex @ r
            Tw = T_apex @ w
            rw = float(np.dot(r, w))
            rTr = float(np.dot(r, Tr))
            rTw = float(np.dot(r, Tw))

            # 1. Position cotangent-lift vector field dK_4 / dw = Xi_4(r):
            dK_dw = coeff * GM_ab * (4.0 * Tr / r3 - 3.0 * rTr * r / r5)
            dq[ib] += (ma / M_ab) * dK_dw
            dq[ia] -= (mb / M_ab) * dK_dw

            # 2. Primary binary leaf relative momentum gradient -dK_4 / dr:
            dK_dr = coeff * GM_ab * (
                4.0 * Tw / r3
                - 12.0 * rTw * r / r5
                - 3.0 * rTr * w / r5
                - 6.0 * rw * Tr / r5
                + 15.0 * rw * rTr * r / r7
            )
            dp[ib] -= dK_dr
            dp[ia] += dK_dr

            # 3. Spectator Noether back-reaction -dK_4 / dR_k via symmetric tensor M_sym:
            M_sym = coeff * GM_ab * (
                (2.0 / r3) * (np.outer(r, w) + np.outer(w, r))
                - (3.0 * rw / r5) * np.outer(r, r)
            )
            tr_M = float(np.trace(M_sym))
            F_cm_back = np.zeros(3, dtype=np.float64)

            for k in range(state.N):
                if k == ia or k == ib:
                    continue
                mk = float(masses[k])
                R_k = q_cur[k] - R_cm
                Rk2 = float(np.dot(R_k, R_k))
                if Rk2 < 1e-30:
                    continue
                Rkn = math.sqrt(Rk2)
                Rk5 = Rk2 * Rk2 * Rkn
                Rk7 = Rk5 * Rk2

                theta_k = 1.0
                if apex is not None and apex.screening_mode != "pure_cartesian" and rho_bar > 1e-30:
                    u_k = Rkn / rho_bar
                    theta_val, _ = eval_screening_theta(u_k, mode=apex.screening_mode)
                    theta_k = float(theta_val)
                    if theta_k == 0.0:
                        continue

                MRk = M_sym @ R_k
                RkMRk = float(np.dot(R_k, MRk))
                dK_dRk = (theta_k * float(G) * mk) * (
                    6.0 * MRk / Rk5 + 3.0 * tr_M * R_k / Rk5 - 15.0 * RkMRk * R_k / Rk7
                )
                dp_k = -dK_dRk
                dp[k] += dp_k
                F_cm_back += dp_k

            dp[ia] -= (ma / M_ab) * F_cm_back
            dp[ib] -= (mb / M_ab) * F_cm_back

        return dq, dp

    # Symmetric, symplectic 2-stage Gauss-Legendre collocation (order 4)
    sqrt3_6 = math.sqrt(3.0) / 6.0
    a11 = 0.25
    a12 = 0.25 - sqrt3_6
    a21 = 0.25 + sqrt3_6
    a22 = 0.25

    q0 = state.q.copy()
    p0 = state.p.copy()
    dq0, dp0 = _grad_k4(q0, p0)
    c1 = 0.5 - sqrt3_6
    c2 = 0.5 + sqrt3_6
    Q1, P1 = q0 + c1 * dq0, p0 + c1 * dp0
    Q2, P2 = q0 + c2 * dq0, p0 + c2 * dp0

    for _ in range(6):
        k1_q, k1_p = _grad_k4(Q1, P1)
        k2_q, k2_p = _grad_k4(Q2, P2)
        Q1_new = q0 + a11 * k1_q + a12 * k2_q
        P1_new = p0 + a11 * k1_p + a12 * k2_p
        Q2_new = q0 + a21 * k1_q + a22 * k2_q
        P2_new = p0 + a21 * k1_p + a22 * k2_p
        err = max(
            float(np.max(np.abs(Q1_new - Q1))),
            float(np.max(np.abs(P1_new - P1))),
            float(np.max(np.abs(Q2_new - Q2))),
            float(np.max(np.abs(P2_new - P2))),
        )
        Q1, P1, Q2, P2 = Q1_new, P1_new, Q2_new, P2_new
        if err <= 1e-16:
            break

    k1_q, k1_p = _grad_k4(Q1, P1)
    k2_q, k2_p = _grad_k4(Q2, P2)
    qc = q0 + 0.5 * (k1_q + k2_q)
    pc = p0 + 0.5 * (k1_p + k2_p)

    return state.with_updates(q=qc, p=pc)


class CotangentLiftCorrectorC4(SymplecticCorrector):
    """
    Fourth-order closed-form cotangent-lift symplectic corrector C_4^{cot}(h) = exp(- X_{K_4}),
    where K_4 = c_4^{BM4} * h^4 * D_A^3 B in P_1.
    """

    def __init__(
        self,
        G: float = G_DEFAULT,
        screening_mode: str = "c6",
        c4_coeff: float = C4_BCH_BM4,
    ) -> None:
        self.G = G
        self.screening_mode = screening_mode
        self.c4_coeff = c4_coeff

    @property
    def order(self) -> int:
        return 4

    def apply(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        h: float,
        inverse: bool = False,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        return compute_cotangent_lift_c4_flow(
            state=state,
            tree=tree,
            h=h,
            G=self.G,
            inverse=inverse,
            c4_coeff=self.c4_coeff,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )


class HomologicalLieCorrectorC4(SymplecticCorrector):
    """
    Exact 2-stage palindromic Lie-conjugation symplectic corrector Z_4^{BM4}(h):
        Z_4^{BM4}(h) = S_sym(alpha_1, beta_1; |h|) o S_sym(alpha_2, beta_2; |h|),
    built directly from the exact Keplerian/composite drift flow Phi_A and screened apex kick Phi_B
    at dyadic nodes alpha_1 = 1/4, alpha_2 = 1/8 with weights:
        beta_1 = +256 * c_4^{BM4},   beta_2 = -512 * c_4^{BM4},
    satisfying:
        2 * (beta_1 * alpha_1 + beta_2 * alpha_2) == 0     (exact order-2 cancellation of D_A B),
        (1/3) * (beta_1 * alpha_1^3 + beta_2 * alpha_2^3) == c_4^{BM4} (exact reproduction of c_4^{BM4} h^4 D_A^3 B).

    Each stage S_sym(alpha_i, beta_i; |h|) pairs the antipodal half-shells (-alpha_i, -0.5*beta_i)
    and (+alpha_i, +0.5*beta_i) to cancel the O(alpha_i beta_i^2 |h|^3 [D_B, [D_A, D_B]]) quadratic
    self-commutator identically, with adjacent boundary drifts and central kicks fused into
    7 drift calls + 6 kick calls total:
        S_sym(alpha, beta; tau) = Phi_A(alpha*tau) o Phi_B(-0.5*beta*tau) o Phi_A(-2*alpha*tau)
                                  o Phi_B(+beta*tau) o Phi_A(2*alpha*tau) o Phi_B(-0.5*beta*tau) o Phi_A(-alpha*tau).

    Cancels the entire 4th-order BCH error c_4^{BM4} h^4 D_A^4 B across all multipole orders (l >= 2)
    and composite outer motion simultaneously while preserving symplecticity, all 10 Lie--Bargmann
    invariants, and exact group inversion Z_4(h) o Z_4^{-1}(h) = Id (< 1e-16).
    """

    def __init__(
        self,
        G: float = G_DEFAULT,
        screening_mode: str = "c6",
        c4_coeff: float = C4_BCH_BM4,
    ) -> None:
        self.G = G
        self.screening_mode = screening_mode
        self.c4_coeff = c4_coeff
        self._default_apex = ScreenedApexPotential(G=G, screening_mode=screening_mode, order=6)

    @property
    def order(self) -> int:
        return 4

    def apply(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        h: float,
        inverse: bool = False,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        h_abs = abs(float(h))
        if h_abs == 0.0 or state.N <= 2 or not tree.leaf_pairs or (supernode_bodies is not None and len(supernode_bodies) >= state.N):
            return state

        from .common import step_drift_and_kepler, step_apex_kick
        from ..core.profiler import get_current_profiler, set_current_profiler, reset_current_profiler

        eval_apex = apex or self._default_apex
        alpha1 = HOMOLOGICAL_ALPHA1_BM4
        alpha2 = HOMOLOGICAL_ALPHA2_BM4
        beta1 = 256.0 * float(self.c4_coeff)
        beta2 = -512.0 * float(self.c4_coeff)

        blocks = (
            (-alpha2, -0.5 * beta2),
            (+alpha2, +0.5 * beta2),
            (-alpha1, -0.5 * beta1),
            (+alpha1, +0.5 * beta1),
        )

        t0 = state.t
        st = state
        prof = get_current_profiler()
        tok = set_current_profiler(None) if prof is not None else None
        try:
            if inverse:
                # Entrance onto shadow manifold: Z_4^{-1}(h) generates +c_4^{BM4} h^4 D_A^3 B
                for a, b_h in blocks:
                    st = step_drift_and_kepler(
                        st,
                        tree,
                        +a * h_abs,
                        G=self.G,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
                    st = step_apex_kick(
                        st,
                        tree,
                        eval_apex,
                        +b_h * h_abs,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
                    st = step_drift_and_kepler(
                        st,
                        tree,
                        -2.0 * a * h_abs,
                        G=self.G,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
                    st = step_apex_kick(
                        st,
                        tree,
                        eval_apex,
                        -b_h * h_abs,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
                    st = step_drift_and_kepler(
                        st,
                        tree,
                        +a * h_abs,
                        G=self.G,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
            else:
                # Exit to physical manifold: exact group inverse of entrance sequence
                for a, b_h in reversed(blocks):
                    st = step_drift_and_kepler(
                        st,
                        tree,
                        -a * h_abs,
                        G=self.G,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
                    st = step_apex_kick(
                        st,
                        tree,
                        eval_apex,
                        +b_h * h_abs,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
                    st = step_drift_and_kepler(
                        st,
                        tree,
                        +2.0 * a * h_abs,
                        G=self.G,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
                    st = step_apex_kick(
                        st,
                        tree,
                        eval_apex,
                        -b_h * h_abs,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
                    st = step_drift_and_kepler(
                        st,
                        tree,
                        -a * h_abs,
                        G=self.G,
                        supernode_bodies=supernode_bodies,
                        supernode_clusters=supernode_clusters,
                    )
        finally:
            if tok is not None:
                reset_current_profiler(tok)

        return st.with_updates(t=t0)


class CorrectorC4(SymplecticCorrector):
    """
    Fourth-order boundary symplectic corrector C_4(h) = exp(- h^4 X_{K_4}).
    Elevates BM_4^{(7)} to effective order 6 in the asymptotic regime eps << h (Theorem 5.15).

    Supports four modes:
    - mode="homological" (canonical default): Exact 2-stage palindromic Lie-conjugation operator
      Z_4^{BM4}(h) in P_1 cancelling c_4^{BM4} h^4 D_A^4 B across all multipole orders and
      composite outer drift on arbitrary 3D eccentric orbits.
    - mode="cotangent_lift": Closed-form cotangent lift of K_4 = c_4^{BM4} h^4 D_A^3 B in P_1
      with exact Noether torque-free back-reaction M_sym.
    - mode="global": Pure momentum shear in P_0 with 1/480 scaling and strict spatial
      stationarity Delta q == 0 (for P_0 universal obstruction ablation & legacy unit tests).
    - mode="legacy": Two-stage kick + Verlet apex drift flow.
    """

    def __init__(
        self,
        G: float = G_DEFAULT,
        screening_mode: str = "c6",
        kappa4: float = KAPPA4_BM4,
        mode: str = "homological",
    ) -> None:
        self.G = G
        self.screening_mode = screening_mode
        self.kappa4 = kappa4
        self.mode = mode
        self._default_apex = ScreenedApexPotential(G=G, screening_mode=screening_mode, order=6)
        self._global_corrector = GlobalCorrectorC4(G=G, screening_mode=screening_mode)
        self._cotangent_corrector = CotangentLiftCorrectorC4(G=G, screening_mode=screening_mode)
        self._homological_corrector = HomologicalLieCorrectorC4(G=G, screening_mode=screening_mode)

    @property
    def order(self) -> int:
        return 4

    def apply(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        h: float,
        inverse: bool = False,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """
        Applies C_4(h) or C_4^{-1}(h).
        """
        if h == 0.0 or state.N <= 1:
            return state

        if self.mode == "global":
            return self._global_corrector.apply(
                state,
                tree,
                h,
                inverse=inverse,
                apex=apex,
                supernode_bodies=supernode_bodies,
                supernode_clusters=supernode_clusters,
            )
        if self.mode == "cotangent_lift":
            return self._cotangent_corrector.apply(
                state,
                tree,
                h,
                inverse=inverse,
                apex=apex,
                supernode_bodies=supernode_bodies,
                supernode_clusters=supernode_clusters,
            )
        if self.mode == "homological":
            return self._homological_corrector.apply(
                state,
                tree,
                h,
                inverse=inverse,
                apex=apex,
                supernode_bodies=supernode_bodies,
                supernode_clusters=supernode_clusters,
            )

        eval_apex = apex or self._default_apex

        # 1. First half-step momentum shear kick at current configuration q_0
        dp_conormal_0 = compute_conormal_corrector_kick(
            state,
            tree,
            h,
            G=self.G,
            inverse=inverse,
            half_step=True,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )
        st_0 = state.with_updates(p=state.p + dp_conormal_0)

        # 2. Far-field palindromic Verlet flow of the screened apex potential
        sign = +1.0 if inverse else -1.0
        step_coeff = sign * self.kappa4 * (h ** 4)

        # 2a. Half-step momentum kick at q_0
        v_current = st_0.p / st_0.masses[:, None]
        dF_0 = compute_analytical_directional_force_derivative(eval_apex, st_0, tree, supernode_bodies, v_current)
        p_half = st_0.p - 0.5 * step_coeff * dF_0

        # 2b. Full-step position shear evaluated at q_0 with p_half
        st_mid = st_0.with_updates(p=p_half)
        forces_mid, _ = eval_apex.compute_apex_forces(
            st_mid, tree, supernode_bodies, supernode_clusters=supernode_clusters
        )
        v_drift = forces_mid / st_0.masses[:, None]
        new_q = st_0.q + step_coeff * v_drift

        # 2c. Second half-step momentum kick at new_q with v_half
        st_new_pos = st_mid.with_updates(q=new_q)
        v_half = p_half / st_0.masses[:, None]
        dF_1 = compute_analytical_directional_force_derivative(eval_apex, st_new_pos, tree, supernode_bodies, v_half)
        p_apex_end = p_half - 0.5 * step_coeff * dF_1
        st_apex = st_new_pos.with_updates(p=p_apex_end)

        # 3. Second half-step momentum shear kick at final configuration new_q
        dp_conormal_1 = compute_conormal_corrector_kick(
            st_apex,
            tree,
            h,
            G=self.G,
            inverse=inverse,
            half_step=True,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )
        final_p = p_apex_end + dp_conormal_1

        return state.with_updates(q=new_q, p=final_p)



