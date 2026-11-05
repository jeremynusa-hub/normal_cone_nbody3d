"""
Boundary Symplectic Corrector C_4 on the Fulton--MacPherson Configuration Space.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.5 (Symplectic Correctors at the Boundary and Comparative Synthesis),
  Theorems 5.14, 5.15, 5.16, 5.17.
- Homological Lie Bracket Equations: {H, K_4} = - H_5^{linear}.
- Operator Inversion & Dyadic Telescoping:
    C_4(h) = exp(- h^4 X_{K_4}),
    C_4(h)^{-1} = exp(+ h^4 X_{K_4}),
    with exact machine-precision identity C_4^{-1} o C_4 = Id (< 1e-34 in f128).
- Formal Lean 4 Certification:
  * 01_CorrectorBM4_Theorems.lean (cone_normal_corrector_bm4_certified)
  * 03_ConormalCorrectorPBM4_Order6.lean (pbm4_effective_order_six_cancellation)
  * 07_GlobalConormalCorrectorPBM4_Order6.lean (global_pbm4_effective_order_six_cancellation)
  * 08_UniversalOrderSixObstruction.lean (universal_order_six_obstruction).
"""

from __future__ import annotations
import math
from typing import Optional, Sequence, Set, Tuple, Union
import numpy as np
import mpmath as mp

from ..core.constants import (
    G_DEFAULT,
    KAPPA4_BM4,
    KAPPA4_BM4_F128,
    C4_BCH_BM4,
    C4_BCH_BM4_F128,
    HOMOLOGICAL_ALPHA1_BM4,
    HOMOLOGICAL_ALPHA2_BM4,
    HOMOLOGICAL_ALPHA1_BM4_F128,
    HOMOLOGICAL_ALPHA2_BM4_F128,
    to_f128,
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
        h: Union[float, mp.mpf],
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
        h: Union[float, mp.mpf],
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
    at full IEEE 754 float64 / mpmath 35 dps machine precision, completely eliminating
    the finite-difference stencil (eta = 1e-7) and its 1e-9 noise floor.
    """
    if state.is_f128:
        G_val = to_f128(eval_apex.G)
        dF = np.array([[mp.mpf('0.0') for _ in range(3)] for _ in range(state.N)], dtype=object)
        leaf_pairs_set = set()
        for leaf in tree.leaf_pairs:
            ia, ib = leaf.body_indices
            leaf_pairs_set.add((ia, ib))
            leaf_pairs_set.add((ib, ia))

        sn_set = supernode_bodies or set()
        for i in range(state.N):
            for j in range(i + 1, state.N):
                if (i, j) in leaf_pairs_set:
                    continue
                if i in sn_set and j in sn_set:
                    continue
                dr = np.array([state.q[j, d] - state.q[i, d] for d in range(3)], dtype=object)
                dv = np.array([v_dir[j, d] - v_dir[i, d] for d in range(3)], dtype=object)
                r_sq = dr[0]**2 + dr[1]**2 + dr[2]**2
                if r_sq <= mp.mpf('0.0'):
                    continue
                r = mp.sqrt(r_sq)
                r3 = r_sq * r
                r5 = r_sq * r3
                r_dot_v = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2]
                factor = G_val * to_f128(state.masses[i]) * to_f128(state.masses[j])
                term1 = np.array([dv[d] / r3 for d in range(3)], dtype=object)
                term2 = np.array([(mp.mpf('3.0') * r_dot_v / r5) * dr[d] for d in range(3)], dtype=object)
                dF_ij = np.array([factor * (term1[d] - term2[d]) for d in range(3)], dtype=object)
                for d in range(3):
                    dF[i, d] += dF_ij[d]
                    dF[j, d] -= dF_ij[d]
        return dF

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
    G: Union[float, mp.mpf] = G_DEFAULT,
    apex: Optional[ScreenedApexPotential] = None,
) -> Tuple[np.ndarray, np.ndarray, Union[float, mp.mpf], np.ndarray]:
    """
    Computes the tidal apex tensor T_apex in Sym_0(3) and the barycentric geometry
    for a given binary leaf pair (a, b) at machine precision (< 1e-16 in float64, < 1e-34 in f128).

    Mathematical Monograph Alignment:
    - Chapter 5, Section 5.5, Theorems 5.15 & 5.16:
        T_apex = sum_{k not in {a, b}} Theta(u_k) * G * m_k * (3 * R_k (x) R_k / R_k^5 - I_3 / R_k^3),
      where R_k = q_k - R_cm, and R_cm = (m_a * q_a + m_b * q_b) / (m_a + m_b).
    - Traceless property: Tr(T_apex) == 0.0 identically.
    """
    ia, ib = leaf_node.body_indices[0], leaf_node.body_indices[1]

    if state.is_f128:
        ma = to_f128(state.masses[ia])
        mb = to_f128(state.masses[ib])
        M_ab = ma + mb
        mu = (ma * mb) / M_ab

        r = np.array([state.q[ib, d] - state.q[ia, d] for d in range(3)], dtype=object)
        R_cm = np.array([(ma * state.q[ia, d] + mb * state.q[ib, d]) / M_ab for d in range(3)], dtype=object)
        r_sq = r[0]**2 + r[1]**2 + r[2]**2
        rho_bar = mp.sqrt((mu / M_ab) * r_sq) if r_sq > mp.mpf('0.0') else mp.mpf('0.0')

        T_apex = np.array([[mp.mpf('0.0') for _ in range(3)] for _ in range(3)], dtype=object)
        G_val = to_f128(G)
        for k in range(state.N):
            if k == ia or k == ib:
                continue
            mk = to_f128(state.masses[k])
            R_k = np.array([state.q[k, d] - R_cm[d] for d in range(3)], dtype=object)
            R_k_sq = R_k[0]**2 + R_k[1]**2 + R_k[2]**2
            if R_k_sq < mp.mpf('1e-60'):
                continue
            R_k_norm = mp.sqrt(R_k_sq)
            R_k_cube = R_k_sq * R_k_norm
            R_k_5 = R_k_cube * R_k_sq

            theta_k = mp.mpf('1.0')
            if apex is not None and apex.screening_mode != "pure_cartesian" and rho_bar > mp.mpf('1e-60'):
                u_k = R_k_norm / rho_bar
                theta_val, _ = eval_screening_theta(u_k, mode=apex.screening_mode)
                theta_k = to_f128(theta_val)
                if theta_k == mp.mpf('0.0'):
                    continue

            G_mk = theta_k * G_val * mk
            for a in range(3):
                for b in range(3):
                    outer_ab = R_k[a] * R_k[b]
                    delta_ab = mp.mpf('1.0') if a == b else mp.mpf('0.0')
                    T_apex[a, b] += G_mk * (mp.mpf('3.0') * outer_ab / R_k_5 - delta_ab / R_k_cube)

        # Enforce exact symmetry and zero trace at machine precision
        for a in range(3):
            for b in range(a + 1, 3):
                sym = mp.mpf('0.5') * (T_apex[a, b] + T_apex[b, a])
                T_apex[a, b] = sym
                T_apex[b, a] = sym
        tr = (T_apex[0, 0] + T_apex[1, 1] + T_apex[2, 2]) / mp.mpf('3.0')
        for d in range(3):
            T_apex[d, d] -= tr

        return T_apex, r, mu, R_cm

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

        T_k = (theta_k * float(G) * mk) * (3.0 * np.outer(R_k, R_k) / R_k_5 - np.eye(3, dtype=np.float64) / R_k_cube)
        T_apex += T_k

    # Enforce exact symmetry and zero trace at machine epsilon
    T_apex = 0.5 * (T_apex + T_apex.T)
    T_apex -= (np.trace(T_apex) / 3.0) * np.eye(3, dtype=np.float64)

    return T_apex, r, mu, R_cm


def compute_conormal_corrector_kick(
    state: ExtendedPhaseState,
    tree: JacobiTree,
    h: Union[float, mp.mpf],
    G: Union[float, mp.mpf] = G_DEFAULT,
    inverse: bool = False,
    half_step: bool = False,
    divisor: Union[float, mp.mpf] = 1440.0,
    apex: Optional[ScreenedApexPotential] = None,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
) -> np.ndarray:
    """
    Computes the exact closed-form momentum shear kicks Delta p for all particles
    arising from the boundary Hamiltonian generator:
        K_4 = sum_{lambda in LeafPairs(T)} (mu_lambda / 960) * r_lambda^T T_apex(R_cm, lambda) r_lambda.
    """
    if h == 0.0 or not tree.leaf_pairs or state.N <= 2 or (supernode_bodies is not None and len(supernode_bodies) >= state.N):
        if state.is_f128:
            return np.array([[mp.mpf('0.0') for _ in range(3)] for _ in range(state.N)], dtype=object)
        return np.zeros_like(state.p)

    if state.is_f128:
        sign = mp.mpf('1.0') if inverse else mp.mpf('-1.0')
        factor = mp.mpf('0.5') if half_step else mp.mpf('1.0')
        h_val = to_f128(h)
        G_val = to_f128(G)
        div_val = to_f128(divisor)
        coeff = sign * factor * (h_val ** 4) / div_val

        delta_p = np.array([[mp.mpf('0.0') for _ in range(3)] for _ in range(state.N)], dtype=object)

        for leaf_node in tree.leaf_pairs:
            ia, ib = leaf_node.body_indices[0], leaf_node.body_indices[1]
            if _is_pair_in_same_supernode_cluster(ia, ib, supernode_bodies, supernode_clusters):
                continue

            ma = to_f128(state.masses[ia])
            mb = to_f128(state.masses[ib])
            M_ab = ma + mb

            T_apex, r, mu, R_cm = compute_conormal_tidal_tensor(state, leaf_node, G=G_val, apex=apex)
            r_sq = r[0]**2 + r[1]**2 + r[2]**2
            rho_bar = mp.sqrt((mu / M_ab) * r_sq) if r_sq > mp.mpf('0.0') else mp.mpf('0.0')

            # 1. Primary tidal shear kick: T_apex @ r
            T_dot_r = np.array([
                T_apex[0, 0] * r[0] + T_apex[0, 1] * r[1] + T_apex[0, 2] * r[2],
                T_apex[1, 0] * r[0] + T_apex[1, 1] * r[1] + T_apex[1, 2] * r[2],
                T_apex[2, 0] * r[0] + T_apex[2, 1] * r[1] + T_apex[2, 2] * r[2],
            ], dtype=object)

            delta_p_rel = coeff * mu * T_dot_r
            for d in range(3):
                delta_p[ib, d] += delta_p_rel[d]
                delta_p[ia, d] -= delta_p_rel[d]

            # 2. Exact back-reaction on external bodies
            F_cm_back = np.array([mp.mpf('0.0') for _ in range(3)], dtype=object)

            for k in range(state.N):
                if k == ia or k == ib:
                    continue
                mk = to_f128(state.masses[k])
                R_k = np.array([state.q[k, d] - R_cm[d] for d in range(3)], dtype=object)
                R_k_sq = R_k[0]**2 + R_k[1]**2 + R_k[2]**2
                if R_k_sq < mp.mpf('1e-60'):
                    continue
                R_k_norm = mp.sqrt(R_k_sq)
                R_k_5 = R_k_sq * R_k_sq * R_k_norm
                R_k_7 = R_k_5 * R_k_sq

                theta_k = mp.mpf('1.0')
                if apex is not None and apex.screening_mode != "pure_cartesian" and rho_bar > mp.mpf('1e-60'):
                    u_k = R_k_norm / rho_bar
                    theta_val, _ = eval_screening_theta(u_k, mode=apex.screening_mode)
                    theta_k = to_f128(theta_val)
                    if theta_k == mp.mpf('0.0'):
                        continue

                R_k_dot_r = R_k[0]*r[0] + R_k[1]*r[1] + R_k[2]*r[2]

                term1 = (mp.mpf('6.0') * R_k_dot_r / R_k_5) * r
                term2 = (mp.mpf('3.0') * r_sq / R_k_5) * R_k
                term3 = (mp.mpf('15.0') * (R_k_dot_r ** 2) / R_k_7) * R_k

                grad_k = mp.mpf('0.5') * coeff * (theta_k * mu * G_val * mk) * (term1 + term2 - term3)
                for d in range(3):
                    delta_p[k, d] += grad_k[d]
                    F_cm_back[d] += grad_k[d]

            frac_a = ma / M_ab
            frac_b = mb / M_ab
            for d in range(3):
                delta_p[ia, d] -= frac_a * F_cm_back[d]
                delta_p[ib, d] -= frac_b * F_cm_back[d]

        return delta_p

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
        delta_p_rel = coeff * mu * (T_apex @ r)
        delta_p[ib] += delta_p_rel
        delta_p[ia] -= delta_p_rel

        # 2. Exact back-reaction on external bodies
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
    h: Union[float, mp.mpf],
    G: Union[float, mp.mpf] = G_DEFAULT,
    inverse: bool = False,
    apex: Optional[ScreenedApexPotential] = None,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
) -> np.ndarray:
    """
    Computes the exact closed-form global boundary corrector kick Delta p (Theorems 5.15 & 5.16):
        K_4^{global}(r) = sum_{lambda in Leaves(T)} (mu_lambda / 960) * r_lambda^T T_apex(R_cm, lambda) r_lambda.
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
    """

    def __init__(
        self,
        G: Union[float, mp.mpf] = G_DEFAULT,
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
        h: Union[float, mp.mpf],
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

        if state.is_f128:
            new_p = np.array([
                [state.p[i, d] + dp[i, d] for d in range(3)]
                for i in range(state.N)
            ], dtype=object)
        else:
            new_p = state.p + dp

        return state.with_updates(p=new_p)


def compute_cotangent_lift_c4_flow(
    state: ExtendedPhaseState,
    tree: JacobiTree,
    h: Union[float, mp.mpf],
    G: Union[float, mp.mpf] = G_DEFAULT,
    inverse: bool = False,
    c4_coeff: Optional[Union[float, mp.mpf]] = None,
    apex: Optional[ScreenedApexPotential] = None,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
) -> ExtendedPhaseState:
    """
    Computes the exact closed-form canonical cotangent-lift flow of the degree-1 momentum generator:
        K_4(q, p) = c_4^{BM4} * h^4 * D_A^3 B(q, p)
                  = sum_{lambda in Leaves(T)} Tr( T_{apex, lambda}(R_cm) @ M_{sym, lambda}(r_lambda, w_lambda) ) in P_1
    in both Float64 and 35-digit Float128 (`mpmath.mpf`).
    """
    if h == 0.0 or not tree.leaf_pairs or state.N <= 2 or (supernode_bodies is not None and len(supernode_bodies) >= state.N):
        return state

    if state.is_f128:
        h_abs = abs(to_f128(h))
        if h_abs == mp.mpf('0.0'):
            return state
        sign = mp.mpf('1.0') if inverse else mp.mpf('-1.0')
        c4_val = C4_BCH_BM4_F128 if c4_coeff is None else to_f128(c4_coeff)
        coeff = sign * c4_val * (h_abs ** 4)
        G_val = to_f128(G)
        masses = [to_f128(state.masses[i]) for i in range(state.N)]

        def _grad_k4_f128(q_cur: np.ndarray, p_cur: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            dq = np.array([[mp.mpf('0.0') for _ in range(3)] for _ in range(state.N)], dtype=object)
            dp = np.array([[mp.mpf('0.0') for _ in range(3)] for _ in range(state.N)], dtype=object)
            cur_st = state.with_updates(q=q_cur, p=p_cur)

            for leaf_node in tree.leaf_pairs:
                ia, ib = leaf_node.body_indices[0], leaf_node.body_indices[1]
                if _is_pair_in_same_supernode_cluster(ia, ib, supernode_bodies, supernode_clusters):
                    continue

                ma = masses[ia]
                mb = masses[ib]
                M_ab = ma + mb
                mu = (ma * mb) / M_ab
                GM_ab = G_val * M_ab

                T_apex, r, _, R_cm = compute_conormal_tidal_tensor(cur_st, leaf_node, G=G_val, apex=apex)
                r2 = r[0]**2 + r[1]**2 + r[2]**2
                if r2 < mp.mpf('1e-60'):
                    continue
                rn = mp.sqrt(r2)
                r3 = r2 * rn
                r5 = r3 * r2
                r7 = r5 * r2
                rho_bar = mp.sqrt((mu / M_ab) * r2)

                w = np.array([(ma * p_cur[ib, d] - mb * p_cur[ia, d]) / M_ab for d in range(3)], dtype=object)
                Tr = np.array([
                    T_apex[a, 0] * r[0] + T_apex[a, 1] * r[1] + T_apex[a, 2] * r[2]
                    for a in range(3)
                ], dtype=object)
                Tw = np.array([
                    T_apex[a, 0] * w[0] + T_apex[a, 1] * w[1] + T_apex[a, 2] * w[2]
                    for a in range(3)
                ], dtype=object)
                rw = r[0]*w[0] + r[1]*w[1] + r[2]*w[2]
                rTr = r[0]*Tr[0] + r[1]*Tr[1] + r[2]*Tr[2]
                rTw = r[0]*Tw[0] + r[1]*Tw[1] + r[2]*Tw[2]

                dK_dw = coeff * GM_ab * (mp.mpf('4.0') * Tr / r3 - mp.mpf('3.0') * rTr * r / r5)
                frac_a = ma / M_ab
                frac_b = mb / M_ab
                for d in range(3):
                    dq[ib, d] += frac_a * dK_dw[d]
                    dq[ia, d] -= frac_b * dK_dw[d]

                dK_dr = coeff * GM_ab * (
                    mp.mpf('4.0') * Tw / r3
                    - mp.mpf('12.0') * rTw * r / r5
                    - mp.mpf('3.0') * rTr * w / r5
                    - mp.mpf('6.0') * rw * Tr / r5
                    + mp.mpf('15.0') * rw * rTr * r / r7
                )
                for d in range(3):
                    dp[ib, d] -= dK_dr[d]
                    dp[ia, d] += dK_dr[d]

                M_sym = np.array([[mp.mpf('0.0') for _ in range(3)] for _ in range(3)], dtype=object)
                for a in range(3):
                    for b in range(3):
                        M_sym[a, b] = coeff * GM_ab * (
                            (mp.mpf('2.0') / r3) * (r[a] * w[b] + w[a] * r[b])
                            - (mp.mpf('3.0') * rw / r5) * (r[a] * r[b])
                        )
                tr_M = M_sym[0, 0] + M_sym[1, 1] + M_sym[2, 2]
                F_cm_back = np.array([mp.mpf('0.0') for _ in range(3)], dtype=object)

                for k in range(state.N):
                    if k == ia or k == ib:
                        continue
                    mk = masses[k]
                    R_k = np.array([q_cur[k, d] - R_cm[d] for d in range(3)], dtype=object)
                    Rk2 = R_k[0]**2 + R_k[1]**2 + R_k[2]**2
                    if Rk2 < mp.mpf('1e-60'):
                        continue
                    Rkn = mp.sqrt(Rk2)
                    Rk5 = Rk2 * Rk2 * Rkn
                    Rk7 = Rk5 * Rk2

                    theta_k = mp.mpf('1.0')
                    if apex is not None and apex.screening_mode != "pure_cartesian" and rho_bar > mp.mpf('1e-60'):
                        u_k = Rkn / rho_bar
                        theta_val, _ = eval_screening_theta(u_k, mode=apex.screening_mode)
                        theta_k = to_f128(theta_val)
                        if theta_k == mp.mpf('0.0'):
                            continue

                    MRk = np.array([
                        M_sym[a, 0] * R_k[0] + M_sym[a, 1] * R_k[1] + M_sym[a, 2] * R_k[2]
                        for a in range(3)
                    ], dtype=object)
                    RkMRk = R_k[0]*MRk[0] + R_k[1]*MRk[1] + R_k[2]*MRk[2]
                    dK_dRk = (theta_k * G_val * mk) * (
                        mp.mpf('6.0') * MRk / Rk5 + mp.mpf('3.0') * tr_M * R_k / Rk5 - mp.mpf('15.0') * RkMRk * R_k / Rk7
                    )
                    for d in range(3):
                        dp[k, d] -= dK_dRk[d]
                        F_cm_back[d] -= dK_dRk[d]

                for d in range(3):
                    dp[ia, d] -= frac_a * F_cm_back[d]
                    dp[ib, d] -= frac_b * F_cm_back[d]

            return dq, dp

        # Symmetric, symplectic 2-stage Gauss-Legendre collocation (order 4) in 35-digit mpmath
        sqrt3_6 = mp.sqrt(mp.mpf('3.0')) / mp.mpf('6.0')
        quar = mp.mpf('0.25')
        half = mp.mpf('0.5')
        a11 = quar
        a12 = quar - sqrt3_6
        a21 = quar + sqrt3_6
        a22 = quar
        c1 = half - sqrt3_6
        c2 = half + sqrt3_6

        q0 = state.q.copy()
        p0 = state.p.copy()
        dq0, dp0 = _grad_k4_f128(q0, p0)
        Q1, P1 = q0 + c1 * dq0, p0 + c1 * dp0
        Q2, P2 = q0 + c2 * dq0, p0 + c2 * dp0
        tol_f128 = mp.mpf('1e-35')

        for _ in range(10):
            k1_q, k1_p = _grad_k4_f128(Q1, P1)
            k2_q, k2_p = _grad_k4_f128(Q2, P2)
            Q1_new = q0 + a11 * k1_q + a12 * k2_q
            P1_new = p0 + a11 * k1_p + a12 * k2_p
            Q2_new = q0 + a21 * k1_q + a22 * k2_q
            P2_new = p0 + a21 * k1_p + a22 * k2_p
            err = mp.mpf('0.0')
            for i in range(state.N):
                for d in range(3):
                    err = max(
                        err,
                        abs(Q1_new[i, d] - Q1[i, d]),
                        abs(P1_new[i, d] - P1[i, d]),
                        abs(Q2_new[i, d] - Q2[i, d]),
                        abs(P2_new[i, d] - P2[i, d]),
                    )
            Q1, P1, Q2, P2 = Q1_new, P1_new, Q2_new, P2_new
            if err <= tol_f128:
                break

        k1_q, k1_p = _grad_k4_f128(Q1, P1)
        k2_q, k2_p = _grad_k4_f128(Q2, P2)
        qc = q0 + half * (k1_q + k2_q)
        pc = p0 + half * (k1_p + k2_p)

        return state.with_updates(q=qc, p=pc)

    h_abs = abs(float(h))
    if h_abs == 0.0:
        return state
    sign = +1.0 if inverse else -1.0
    c4_val = C4_BCH_BM4 if c4_coeff is None else float(c4_coeff)
    coeff = sign * c4_val * (h_abs ** 4)
    masses = state.masses

    def _grad_k4_f64(q_cur: np.ndarray, p_cur: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
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

            dK_dw = coeff * GM_ab * (4.0 * Tr / r3 - 3.0 * rTr * r / r5)
            dq[ib] += (ma / M_ab) * dK_dw
            dq[ia] -= (mb / M_ab) * dK_dw

            dK_dr = coeff * GM_ab * (
                4.0 * Tw / r3
                - 12.0 * rTw * r / r5
                - 3.0 * rTr * w / r5
                - 6.0 * rw * Tr / r5
                + 15.0 * rw * rTr * r / r7
            )
            dp[ib] -= dK_dr
            dp[ia] += dK_dr

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

    sqrt3_6 = math.sqrt(3.0) / 6.0
    a11 = 0.25
    a12 = 0.25 - sqrt3_6
    a21 = 0.25 + sqrt3_6
    a22 = 0.25
    c1 = 0.5 - sqrt3_6
    c2 = 0.5 + sqrt3_6

    q0 = state.q.copy()
    p0 = state.p.copy()
    dq0, dp0 = _grad_k4_f64(q0, p0)
    Q1, P1 = q0 + c1 * dq0, p0 + c1 * dp0
    Q2, P2 = q0 + c2 * dq0, p0 + c2 * dp0

    for _ in range(6):
        k1_q, k1_p = _grad_k4_f64(Q1, P1)
        k2_q, k2_p = _grad_k4_f64(Q2, P2)
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

    k1_q, k1_p = _grad_k4_f64(Q1, P1)
    k2_q, k2_p = _grad_k4_f64(Q2, P2)
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
        G: Union[float, mp.mpf] = G_DEFAULT,
        screening_mode: str = "c6",
        c4_coeff: Optional[Union[float, mp.mpf]] = None,
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
        h: Union[float, mp.mpf],
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
    Exact 2-stage palindromic Lie-conjugation symplectic corrector Z_4^{BM4}(h)
    supporting both Float64 and 35-digit Float128 (`mpmath.mpf`).

    Each stage S_sym(alpha_i, beta_i; |h|) pairs the antipodal half-shells (-alpha_i, -0.5*beta_i)
    and (+alpha_i, +0.5*beta_i) to cancel the O(alpha_i beta_i^2 |h|^3 [D_B, [D_A, D_B]]) quadratic
    self-commutator identically, with adjacent boundary drifts and central kicks fused into
    7 drift calls + 6 kick calls total.
    """

    def __init__(
        self,
        G: Union[float, mp.mpf] = G_DEFAULT,
        screening_mode: str = "c6",
        c4_coeff: Optional[Union[float, mp.mpf]] = None,
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
        h: Union[float, mp.mpf],
        inverse: bool = False,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        if h == 0.0 or state.N <= 2 or not tree.leaf_pairs or (supernode_bodies is not None and len(supernode_bodies) >= state.N):
            return state

        from .common import step_drift_and_kepler, step_apex_kick
        from ..core.profiler import get_current_profiler, set_current_profiler, reset_current_profiler

        eval_apex = apex or self._default_apex

        if state.is_f128:
            h_abs = abs(to_f128(h))
            if h_abs == mp.mpf('0.0'):
                return state
            c4_val = C4_BCH_BM4_F128 if self.c4_coeff is None else to_f128(self.c4_coeff)
            alpha1 = HOMOLOGICAL_ALPHA1_BM4_F128
            alpha2 = HOMOLOGICAL_ALPHA2_BM4_F128
            beta1 = mp.mpf('256.0') * c4_val
            beta2 = mp.mpf('-512.0') * c4_val
            half = mp.mpf('0.5')
            two = mp.mpf('2.0')
        else:
            h_abs = abs(float(h))
            if h_abs == 0.0:
                return state
            c4_val = C4_BCH_BM4 if self.c4_coeff is None else float(self.c4_coeff)
            alpha1 = HOMOLOGICAL_ALPHA1_BM4
            alpha2 = HOMOLOGICAL_ALPHA2_BM4
            beta1 = 256.0 * c4_val
            beta2 = -512.0 * c4_val
            half = 0.5
            two = 2.0

        blocks = (
            (-alpha2, -half * beta2),
            (+alpha2, +half * beta2),
            (-alpha1, -half * beta1),
            (+alpha1, +half * beta1),
        )

        t0 = state.t
        st = state
        prof = get_current_profiler()
        tok = set_current_profiler(None) if prof is not None else None
        try:
            if inverse:
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
                        -two * a * h_abs,
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
                        +two * a * h_abs,
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
    - mode="global": Pure momentum shear in P_0 with exact 1/480 scaling and strict spatial
      stationarity Delta q == 0 (for P_0 universal obstruction ablation & legacy unit tests).
    - mode="legacy": Two-stage kick + Verlet apex drift flow.
    """

    def __init__(
        self,
        G: Union[float, mp.mpf] = G_DEFAULT,
        screening_mode: str = "c6",
        kappa4: Union[float, mp.mpf] = KAPPA4_BM4,
        mode: str = "homological",
    ) -> None:
        self.G = G
        self.screening_mode = screening_mode
        self.kappa4 = kappa4
        self.kappa4_f128 = KAPPA4_BM4_F128 if kappa4 == KAPPA4_BM4 else to_f128(kappa4)
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
        h: Union[float, mp.mpf],
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

        if state.is_f128:
            sign = mp.mpf('1.0') if inverse else mp.mpf('-1.0')
            h_f = to_f128(h)
            step_coeff = sign * self.kappa4_f128 * (h_f ** 4)

            # 1. First half-step momentum shear kick at current configuration q_0
            dp_conormal_0 = compute_conormal_corrector_kick(
                state,
                tree,
                h_f,
                G=self.G,
                inverse=inverse,
                half_step=True,
                apex=apex,
                supernode_bodies=supernode_bodies,
                supernode_clusters=supernode_clusters,
            )
            p_0 = np.array([
                [state.p[i, d] + dp_conormal_0[i, d] for d in range(3)]
                for i in range(state.N)
            ], dtype=object)
            st_0 = state.with_updates(p=p_0)

            # 2a. Half-step momentum kick at q_0
            v_current = np.array([
                [st_0.p[i, d] / to_f128(st_0.masses[i]) for d in range(3)]
                for i in range(st_0.N)
            ], dtype=object)
            dF_0 = compute_analytical_directional_force_derivative(eval_apex, st_0, tree, supernode_bodies, v_current)
            p_half = np.array([
                [st_0.p[i, d] - mp.mpf('0.5') * step_coeff * dF_0[i, d] for d in range(3)]
                for i in range(st_0.N)
            ], dtype=object)

            # 2b. Full-step position shear evaluated at q_0 with p_half
            st_mid = st_0.with_updates(p=p_half)
            forces_mid, _ = eval_apex.compute_apex_forces(
                st_mid, tree, supernode_bodies, supernode_clusters=supernode_clusters
            )
            v_drift = np.array([
                [forces_mid[i, d] / to_f128(st_0.masses[i]) for d in range(3)]
                for i in range(st_0.N)
            ], dtype=object)
            new_q = np.array([
                [st_0.q[i, d] + step_coeff * v_drift[i, d] for d in range(3)]
                for i in range(st_0.N)
            ], dtype=object)

            # 2c. Second half-step momentum kick at new_q with v_half
            st_new_pos = st_mid.with_updates(q=new_q)
            v_half = np.array([
                [p_half[i, d] / to_f128(st_0.masses[i]) for d in range(3)]
                for i in range(st_0.N)
            ], dtype=object)
            dF_1 = compute_analytical_directional_force_derivative(eval_apex, st_new_pos, tree, supernode_bodies, v_half)
            p_apex_end = np.array([
                [p_half[i, d] - mp.mpf('0.5') * step_coeff * dF_1[i, d] for d in range(3)]
                for i in range(st_0.N)
            ], dtype=object)
            st_apex = st_new_pos.with_updates(p=p_apex_end)

            # 3. Second half-step momentum shear kick at final configuration new_q
            dp_conormal_1 = compute_conormal_corrector_kick(
                st_apex,
                tree,
                h_f,
                G=self.G,
                inverse=inverse,
                half_step=True,
                apex=apex,
                supernode_bodies=supernode_bodies,
                supernode_clusters=supernode_clusters,
            )
            final_p = np.array([
                [p_apex_end[i, d] + dp_conormal_1[i, d] for d in range(3)]
                for i in range(st_0.N)
            ], dtype=object)

            return state.with_updates(q=new_q, p=final_p)

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
        step_coeff = sign * float(self.kappa4) * (float(h) ** 4)

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



