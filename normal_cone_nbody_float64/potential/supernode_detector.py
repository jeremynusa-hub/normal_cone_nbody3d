r"""
Autonomous Intrinsic Hill Tidal Supernode Detector and Compactness Thresholds.

Theoretical Foundations:
- Monograph Reference: Chapter 4, Section 4.3 (Intrinsic Compactness and Marchal-Sundman Energy Threshold),
  Proposition 4.5 and Equations (4.181)-(4.227).
- Dimensionless Hill Tidal Ratio:
    eta_S(q) = [bar_rho_S / ||r_P||] * (M_{ext} / M_S)^{1/3} <= eta_{crit} = 3^{-1/3} approx 0.693361.
    Strictly invariant under spatial homotheties q |-> lambda q.
- Natural Energy Threshold (Marchal-Sundman):
    delta_{cluster}^{nat}(S) = eta_{crit} * [G * M_{tot}^{5/3} * M_S^{5/6}] / (2 |E_0|).
    Entirely governed by physical first integrals (G, masses, initial energy E_0) without ad-hoc cutoffs (d_0 = 0).
- Transitive Non-Transversality:
    If Delta_{12} and Delta_{23} are clustered, transitive coalescence enforces the supernode {1, 2, 3}.
- Certified Formal Proofs:
    Lean 4 modules `01_PartitionLattice_Defs.lean` (`collision_diagonal_non_transversality`) and
    `01_UnifiedSuperNodeDef_Theorems.lean`.
"""

from __future__ import annotations
import math
from typing import Tuple, List, Set, Optional, Dict, Union
import numpy as np

from ..core.constants import (
    GRAVITATIONAL_CONSTANT,
    ETA_CRITICAL_HILL,
    ETA_ENTER_HILL,
    ETA_EXIT_HILL,
    LAMBDA_CRITICAL_ROCHE,
    LAMBDA_ENTER_ROCHE,
    LAMBDA_EXIT_ROCHE,
)
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree, JacobiTopologyNode


def compute_tidal_ratio(
    bar_rho_S: float,
    r_P: float,
    M_S: float,
    M_ext: float,
) -> float:
    """
    Computes the dimensionless Hill tidal ratio:
        eta_S = (bar_rho_S / r_P) * (M_ext / M_S)^{1/3}.

    Parameters:
        bar_rho_S: Mass-normalized hyper-radius rho_S / sqrt{M_S} in [L].
        r_P: Distance from cluster barycenter to the nearest external perturber in [L].
        M_S: Total mass of the candidate cluster.
        M_ext: Mass of the external perturber.

    Returns:
        eta_S: Dimensionless ratio (cluster is gravitationally bound if eta_S <= 3^{-1/3}).
    """
    if r_P <= 0.0 or M_S <= 0.0:
        return float('inf')
    mass_ratio = max(M_ext / M_S, 0.0)
    return float((bar_rho_S / r_P) * (mass_ratio ** (1.0 / 3.0)))


def compute_cubic_tidal_force_ratio(
    bar_rho_S: float,
    r_P: float,
    M_S: float,
    M_ext: float,
) -> float:
    """
    Computes the dimensionless cubic tidal force ratio (Section 7.2.1, Equation 7.34):
        Lambda_S = (1 / 3) * eta_S^3 = (M_ext / (3 * M_S)) * (bar_rho_S / r_P)^3.

    Parameters:
        bar_rho_S: Mass-normalized hyper-radius rho_S / sqrt{M_S} in [L].
        r_P: Distance from cluster barycenter to the nearest external perturber in [L].
        M_S: Total mass of the candidate cluster.
        M_ext: Mass of the external perturber.

    Returns:
        Lambda_S: Dimensionless cubic tidal ratio (Roche tidal disruption occurs when Lambda_S >= 0.25).
    """
    if r_P <= 0.0 or M_S <= 0.0:
        return float('inf')
    mass_ratio = max(M_ext / (3.0 * M_S), 0.0)
    return float(mass_ratio * ((bar_rho_S / r_P) ** 3.0))


def compute_natural_threshold(
    G: float,
    M_tot: float,
    M_S: float,
    E0: float,
    eta_crit: float = ETA_CRITICAL_HILL,
) -> float:
    """
    Computes the natural physical Marchal-Sundman clustering length:
        delta_{cluster}^{nat}(S) = eta_{crit} * [G * M_{tot}^{5/3} * M_S^{5/6}] / (2 |E_0|).
    """
    abs_E0 = abs(E0)
    if abs_E0 < 1e-15:
        # Zero-energy / parabolic threshold fallback (unbound systems do not form bound clusters)
        return 0.0

    numerator = G * (M_tot ** (5.0 / 3.0)) * (M_S ** (5.0 / 6.0))
    denominator = 2.0 * abs_E0
    return float(eta_crit * numerator / denominator)


def compute_internal_hierarchy_ratio(
    q: np.ndarray,
    body_indices: Tuple[int, ...] | Set[int] | List[int],
) -> float:
    """
    Computes the dimensionless internal hierarchy ratio:
        chi_S = r_max^{int} / r_min^{int}
    among all distinct pairwise separations inside the candidate cluster.

    Mathematical Alignment:
    - Monograph Reference: Chapter 6, Definition 6.1 and Theorem 6.1.
    - If chi_S >= 3.0, the subsystem contains a well-isolated binary sub-pair
      that is cleanly integrated by the hierarchical Jacobi tree (Kepler leaf).
      Collective supernode regularization is only mandated when the cluster
      is internally non-hierarchical (chi_S < 3.0), triggering the Excluded
      Middle Dilemma.
    """
    idx = list(body_indices)
    k = len(idx)
    if k < 2:
        return 1.0
    q_sub = q[idx]
    diff = q_sub[:, None, :] - q_sub[None, :, :]
    dists = np.linalg.norm(diff, axis=-1)
    i_upper = np.triu_indices(k, k=1)
    pair_dists = dists[i_upper]
    r_min = float(np.min(pair_dists))
    r_max = float(np.max(pair_dists))
    if r_min <= 1e-15:
        return 1.0  # Exact collision is non-hierarchical
    return float(r_max / r_min)


def _is_central_homothetic_collision(
    state: ExtendedPhaseState,
    G: float = GRAVITATIONAL_CONSTANT,
) -> bool:
    """
    Canonical Rule 1 (§4.3.2 & §6.4.2):
    Determines whether the full N-body system S = {0, ..., N-1} (N >= 3, without external perturbers)
    is undergoing a simultaneous radial central-configuration homothetic collapse (or equidistant
    symmetric simultaneous collision, e.g. equilateral triangle N=3, regular tetrahedron N=4, D_4 square N=4).

    Returns False for general non-homothetic N-body scattering (e.g. Burrau 3-4-5, Figure-8, Kozai-Lidov),
    ensuring those systems use the outer Jacobi atlas, DyadicScheduler k_sym, and PBM4 shadow gauge connection.
    """
    N = state.N
    if N < 3:
        return False
    m = np.asarray(state.masses, dtype=np.float64)
    q = np.asarray(state.q, dtype=np.float64)
    p = np.asarray(state.p, dtype=np.float64)
    M_tot = float(np.sum(m))
    if M_tot <= 0.0:
        return False

    q_cdm = np.sum(m[:, None] * q, axis=0) / M_tot
    p_cdm = np.sum(p, axis=0)
    x = q - q_cdm
    v = p / m[:, None] - (p_cdm / M_tot)

    r_norms = np.linalg.norm(x, axis=1)
    bar_rho_sq = float(np.sum(m * (r_norms ** 2))) / M_tot
    if bar_rho_sq <= 1e-30:
        return True

    # 1. Vanishing internal angular momentum: ||L_int||^2 <= 1e-20 * (bar_rho_S^2 + 1e-24)
    L_int = np.sum(np.cross(x, m[:, None] * v), axis=0)
    L_sq = float(np.dot(L_int, L_int))
    if L_sq > 1e-20 * (bar_rho_sq + 1e-24):
        return False

    # 2. Purely radial kinematics (at rest or radial contraction/expansion)
    for i in range(N):
        v_norm_sq = float(np.dot(v[i], v[i]))
        if r_norms[i] <= 1e-15:
            if v_norm_sq > 1e-24:
                return False
        elif v_norm_sq > 1e-24:
            cross_xv = np.cross(x[i], v[i])
            if float(np.dot(cross_xv, cross_xv)) > 1e-16 * (r_norms[i] ** 2) * v_norm_sq:
                return False

    # 3. Either exact central-configuration force alignment (a_i = -lambda * x_i, equal-rate homothety)
    #    or exact equidistant regular simplex (r_max / r_min <= 1.001 with v = 0)
    diff = q[:, None, :] - q[None, :, :]
    dists = np.linalg.norm(diff, axis=-1)
    pair_dists = dists[np.triu_indices(N, k=1)]
    r_min_p = float(np.min(pair_dists))
    r_max_p = float(np.max(pair_dists))
    if r_min_p <= 1e-15:
        return True
    if (r_max_p / r_min_p) >= 3.0:
        return False
    if (r_max_p / r_min_p) <= 1.001 and float(np.max(np.linalg.norm(v, axis=1))) <= 1e-12:
        return True

    # Evaluate internal gravitational accelerations in barycentric frame
    acc = np.zeros_like(x)
    for i in range(N):
        for j in range(N):
            if i != j:
                dr = x[j] - x[i]
                d = float(np.linalg.norm(dr))
                if d > 1e-15:
                    acc[i] += G * m[j] * dr / (d ** 3)

    lambdas = []
    for i in range(N):
        a_norm = float(np.linalg.norm(acc[i]))
        if r_norms[i] <= 1e-15:
            if a_norm > 1e-12:
                return False
            continue
        if a_norm <= 1e-30:
            return False
        cross_ax = np.cross(acc[i], x[i])
        sin_sq = float(np.dot(cross_ax, cross_ax)) / ((a_norm ** 2) * (r_norms[i] ** 2))
        if sin_sq > 1e-4:
            return False
        lam_i = -float(np.dot(acc[i], x[i])) / (r_norms[i] ** 2)
        if lam_i <= 0.0:
            return False
        lambdas.append(lam_i)

    if not lambdas:
        return False
    lam_mean = float(np.mean(lambdas))
    lam_std = float(np.std(lambdas))
    return bool(lam_mean > 0.0 and (lam_std / lam_mean) < 0.02)


def detect_supernodes(
    state: ExtendedPhaseState,
    tree: JacobiTree,
    eta_crit: float = ETA_CRITICAL_HILL,
    active_supernodes: Optional[Union[List[Tuple[int, ...]], Set[Tuple[int, ...]]]] = None,
    eta_enter: float = ETA_ENTER_HILL,
    eta_exit: float = ETA_EXIT_HILL,
    G: float = GRAVITATIONAL_CONSTANT,
    max_hierarchy_ratio: float = 3.0,
    enforce_non_hierarchical: bool = True,
    mode: str = "canonical",
    lambda_crit: float = LAMBDA_CRITICAL_ROCHE,
    lambda_enter: float = LAMBDA_ENTER_ROCHE,
    lambda_exit: float = LAMBDA_EXIT_ROCHE,
    use_energy_threshold: bool = False,
) -> List[Tuple[int, ...]]:
    """
    Detects all tight multi-body clusters (k >= 3) requiring synchronous supernode treatment,
    stabilized by a canonical spectral Schmitt trigger band [eta_enter, eta_exit] = [2^(-1/3), 2^(-1/6)]
    and Roche tidal disruption band [lambda_exit, lambda_enter] = [0.18, 0.25].

    Modes:
    - "canonical" (default, unified Chapters 4, 5 & 6 specification):
      1. Rule 1 (Full-System |S| = N >= 3): Promotes the entire system S = {0,...,N-1} if and only if
         it is undergoing an exact central-configuration radial homothetic collapse (_is_central_homothetic_collision).
      2. Rule 2 (Proper Sub-Cluster 3 <= |S| < N): Promotes S subsetneq {0,...,N-1} if and only if EITHER:
         - Channel A (Sub-Roche Tidal Disruption, §4.3 & §6.2.3 Eq. 6.135): Lambda_S >= lambda_thresh
           and the cluster is compact near the primary (bar_rho_S <= 0.1 * r_P); OR
         - Channel B (Strongly Coupled Self-Gravitating Non-Hierarchical Sub-Cluster, Prop. 4.4):
           eta_S <= eta_thresh, chi_S < 3.0, and the sub-cluster has non-negligible comparable self-gravitation
           (M_S / M_P >= 1e-6 and min_{i in S} m_i / max_{j in S} m_j >= 1e-3).
    - "hill_bound" / "roche_disruption" / "combined":
      Explicit legacy modes preserved for unit-test backward compatibility.

    Returns:
        List of disjoint body index tuples, each containing >= 3 bodies.
    """
    if mode not in ("canonical", "hill_bound", "roche_disruption", "combined"):
        raise ValueError(f"Unsupported supernode detection mode: {mode!r}")

    active_sets: Optional[List[Set[int]]] = None
    if active_supernodes is not None:
        active_sets = [set(c) for c in active_supernodes]

    candidate_clusters: List[Set[int]] = []

    # Check internal nodes representing >= 3 bodies
    for node in tree.internal_nodes:
        bodies = set(node.body_indices)
        k = len(bodies)
        if k < 3:
            continue

        # Optional internal non-hierarchy check
        if enforce_non_hierarchical:
            chi_S = compute_internal_hierarchy_ratio(state.q, bodies)
            if chi_S >= max_hierarchy_ratio:
                # Internally hierarchical subsystem: well-isolated pair integrated via Jacobi leaf
                continue

        # Nearest external perturber: Marchal-Sundman and Hill tidal evaluations apply to
        # strict sub-clusters S \subsetneq {1, ..., N} interacting with external bodies.
        ext_bodies = [b for b in range(state.N) if b not in bodies]
        if not ext_bodies:
            # Root node has no external perturber; evaluated globally below if no sub-clusters exist
            continue

        idx = list(bodies)
        m_sub = state.masses[idx]
        q_sub = state.q[idx]
        M_S = float(np.sum(m_sub))
        cdm_S = np.sum(m_sub[:, None] * q_sub, axis=0) / M_S

        # Internal hyper-radius
        rho_S = float(np.sqrt(max(0.0, np.sum(m_sub * np.sum((q_sub - cdm_S)**2, axis=1)))))
        bar_rho_S = rho_S / math.sqrt(M_S)

        # Optional global Marchal-Sundman energy threshold (only when explicitly enabled in hill_bound mode)
        if use_energy_threshold and mode == "hill_bound" and hasattr(state, "E0") and abs(state.E0) > 1e-15:
            delta_nat = compute_natural_threshold(
                G=G, M_tot=float(np.sum(state.masses)), M_S=M_S, E0=state.E0, eta_crit=eta_crit
            )
            if rho_S <= delta_nat:
                candidate_clusters.append(bodies)
                continue

        eta = 0.0
        lambda_val = 0.0
        dom_m_eb = 0.0
        dom_d_eb = float('inf')
        for eb in ext_bodies:
            d = float(np.linalg.norm(state.q[eb] - cdm_S))
            m_eb = float(state.masses[eb])
            eta_eb = compute_tidal_ratio(bar_rho_S, d, M_S, m_eb)
            if eta_eb > eta:
                eta = eta_eb
                dom_m_eb = m_eb
                dom_d_eb = d
            lam_eb = compute_cubic_tidal_force_ratio(bar_rho_S, d, M_S, m_eb)
            if lam_eb > lambda_val:
                lambda_val = lam_eb

        # Schmitt trigger hysteresis evaluation
        is_already_active = False
        if active_sets is not None:
            is_already_active = any(bodies == a or bodies.issubset(a) for a in active_sets)

        hill_match = False
        if mode in ("canonical", "hill_bound", "combined"):
            effective_thresh = eta_exit if is_already_active else (eta_enter if active_sets is not None else eta_crit)
            hill_match = (eta <= effective_thresh)
            if hill_match and mode == "canonical":
                # Rule 2 Channel B: require comparable internal masses and non-dust self-gravity
                min_m_int = float(np.min(m_sub))
                max_m_int = float(np.max(m_sub))
                comparable_internal = (max_m_int > 0.0 and (min_m_int / max_m_int) >= 1e-3)
                self_gravitating = (dom_m_eb <= 0.0 or (M_S / dom_m_eb) >= 1e-6)
                hill_match = bool(comparable_internal and self_gravitating)

        roche_match = False
        if mode in ("canonical", "roche_disruption", "combined"):
            effective_lam_thresh = (
                lambda_exit if is_already_active else (lambda_enter if active_sets is not None else lambda_crit)
            )
            roche_match = (lambda_val >= effective_lam_thresh)
            if roche_match and mode == "canonical":
                # Rule 2 Channel A: localized sub-cluster near primary (bar_rho_S <= 0.1 * r_P)
                roche_match = bool(bar_rho_S <= 0.1 * dom_d_eb)

        if hill_match or roche_match:
            candidate_clusters.append(bodies)

    # 2. Evaluation of isolated systems without external perturbers (|S| = N >= 3):
    if not candidate_clusters and state.N >= 3 and getattr(state, "E0", 0.0) <= 1e-15:
        if mode == "canonical":
            # Canonical Rule 1: Full-system promotion occurs ONLY for central-configuration homothetic collapse
            if _is_central_homothetic_collision(state, G=G):
                candidate_clusters.append(set(range(state.N)))
        elif mode != "roche_disruption":
            chi_all = compute_internal_hierarchy_ratio(state.q, tuple(range(state.N)))
            if chi_all < max_hierarchy_ratio:
                candidate_clusters.append(set(range(state.N)))

    if not candidate_clusters:
        return []

    # 3. Transitive non-transversality for symmetric multi-body collision geometries (e.g. tetrahedral N=4):
    # By Lean 4 theorem 'collision_diagonal_non_transversality', if an external body has pairwise
    # distances to all bodies of an active cluster matching its internal collision scale r_max_int,
    # it belongs to the same unstratified collision singularity and is transitively coalesced.
    augmented = []
    for cluster in candidate_clusters:
        aug = set(cluster)
        idx = list(cluster)
        q_sub = state.q[idx]
        diff = q_sub[:, None, :] - q_sub[None, :, :]
        dists = np.linalg.norm(diff, axis=-1)
        r_max_int = float(np.max(dists))

        for eb in range(state.N):
            if eb not in aug:
                d_to_bodies = np.linalg.norm(state.q[idx] - state.q[eb], axis=1)
                if np.all(d_to_bodies <= r_max_int * (1.0 + 1e-6)):
                    cand_aug = aug | {eb}
                    if not enforce_non_hierarchical or compute_internal_hierarchy_ratio(state.q, cand_aug) < max_hierarchy_ratio:
                        aug.add(eb)
        augmented.append(aug)

    # Transitive coalescence: merge any clusters sharing common bodies
    merged: List[Set[int]] = []
    for cluster in augmented:
        overlapping = [m for m in merged if not m.isdisjoint(cluster)]
        if not overlapping:
            merged.append(set(cluster))
        else:
            new_set = set(cluster)
            for m in overlapping:
                new_set.update(m)
                merged.remove(m)
            merged.append(new_set)

    # Filter to clusters with >= 3 bodies and sort tuples
    return [tuple(sorted(list(c))) for c in merged if len(c) >= 3]

