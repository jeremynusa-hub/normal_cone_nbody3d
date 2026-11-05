"""
04_SymplecticSubcycle_Interface.py
Standalone Python Interface for the Symplectic Sub-Cycle (Phase 3).

Colocalized with Chapter 6, Section 6.4 (Appendix A.25):
    "The Symplectic Sub-Cycle: Kepler--Sundman Analytic Flow and Screened Apex Kick"

Reference Implementation Interface:
Re-exports and delegates directly to the canonical reference solvers and integrators:
1. stumpff_c, Z_MIN_HYPERBOLIC: Universal Stumpff series from solvers.stumpff.
2. solve_sundman_time_equation, step_kepler_sundman, KeplerSundmanSolver:
   Unbiased incremental Gauss-Stumpff leaf flow (no non-symplectic velocity-norm
   projection) with analytical Weierstrass pericenter conditioning down to 1-e = 1e-12
   and Levi-Civita head-on bounce at d_0 = 0 (solvers.kepler_sundman).
3. step_supernode_synchronous, SupernodeSolver: Time-reversible multi-substep
   synchronous supernode flow with symmetric dyadic bisection stack and symmetric
   half-step close-pair bridge (solvers.supernode).
4. compute_intercluster_forces, ScreenedApexPotential: Screened apex kick with
   intra-cluster supernode filtering preserving inter-supernode cross-pairs for N >= 6.
5. Base Palindromic Symplectic Integrators: StrangIntegrator, Yoshida4Integrator,
   BlanesMoan4Integrator (7-stage Strang composition BM_4^{(7)}), BlanesMoan6Integrator.
"""

from __future__ import annotations
import sys
from pathlib import Path
from typing import Tuple, List, Optional, Set, Sequence, Union
import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from normal_cone_nbody_float64.solvers.stumpff import stumpff_c
from normal_cone_nbody_float64.solvers.kepler_sundman import (
    Z_MIN_HYPERBOLIC,
    solve_sundman_time_equation,
    step_kepler_sundman,
)
from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
    OptimizedKeplerSolver as KeplerSundmanSolver,
)
from normal_cone_nbody_float64.solvers.supernode import (
    step_supernode_synchronous as _canonical_step_supernode_synchronous,
)
from normal_cone_nbody_float64.solvers.optimized_supernode import (
    OptimizedSupernodeSolver as SupernodeSolver,
)
from normal_cone_nbody_float64.potential.screened_apex import (
    ScreenedApexPotential,
    _normalize_supernode_clusters,
)
from normal_cone_nbody_float64.integrators.strang import StrangIntegrator
from normal_cone_nbody_float64.integrators.yoshida4 import Yoshida4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator


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
    where V_inter excludes pairs belonging to terminal leaves and intra-supernode
    pairs (while preserving inter-supernode cross-pairs for N >= 6).
    Guarantees sum F_i == 0 and sum q_i x F_i == 0 to IEEE 754 precision (< 1e-15).
    """
    diff = q[None, :, :] - q[:, None, :]  # diff[i, j] = q[j] - q[i]
    dist_sq = np.sum(diff**2, axis=-1)
    np.fill_diagonal(dist_sq, 1.0)
    dist = np.sqrt(dist_sq)

    with np.errstate(divide="ignore", invalid="ignore"):
        inv_dist3 = np.where(dist_sq > 0.0, 1.0 / (dist_sq * dist), 0.0)
    np.fill_diagonal(inv_dist3, 0.0)

    # Mask terminal leaf pairs (handled analytically by Sundman flow)
    for i, j in leaf_pair_set:
        inv_dist3[i, j] = 0.0
        inv_dist3[j, i] = 0.0

    # Mask intra-cluster pairs of active synchronous supernodes (preserving inter-supernode pairs)
    clusters = _normalize_supernode_clusters(supernode_bodies, supernode_clusters)
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


def step_supernode_synchronous(
    q: np.ndarray,
    p: np.ndarray,
    masses: np.ndarray,
    body_indices_or_dt: Optional[Union[Tuple[int, ...], List[int], Set[int], float]] = None,
    dt: Optional[float] = None,
    G: float = 1.0,
    n_substeps: Optional[int] = None,
    drift_center_of_mass: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Delegates directly to the canonical time-reversible multi-substep synchronous
    supernode solver in solvers.supernode.
    Supports both:
    - Local cluster signature: step_supernode_synchronous(q_cluster, p_cluster, masses_cluster, dt, G=1.0)
    - Full system signature: step_supernode_synchronous(q, p, masses, body_indices, dt, G=1.0)
    """
    if isinstance(body_indices_or_dt, (tuple, list, set)):
        body_indices = tuple(int(i) for i in body_indices_or_dt)
        eff_dt = 0.0 if dt is None else float(dt)
    else:
        body_indices = tuple(range(len(masses)))
        eff_dt = float(dt if body_indices_or_dt is None else body_indices_or_dt)

    return _canonical_step_supernode_synchronous(
        q=q,
        p=p,
        masses=masses,
        body_indices=body_indices,
        dt=eff_dt,
        G=float(G),
        n_substeps=n_substeps,
        drift_center_of_mass=drift_center_of_mass,
    )


def step_strang_subcycle(
    q: np.ndarray,
    p: np.ndarray,
    masses: np.ndarray,
    leaf_pairs: List[Tuple[int, int]],
    dt: float,
    supernodes: Optional[List[Tuple[int, ...]]] = None,
    G: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Executes a complete palindromic Strang sub-cycle S_2(dt) combining:
    1. Half screened apex kick exp((dt/2) X_{V_inter}) (Phase 3.3).
    2. Analytic Kepler--Sundman flow on isolated leaf pairs (Phase 3.1),
       synchronous Sundman integration on active supernodes (Phase 3.2),
       and composite barycentric kinetic drift.
    3. Half screened apex kick exp((dt/2) X_{V_inter}) (Phase 3.3).
    """
    q_curr = np.array(q, dtype=np.float64, copy=True)
    p_curr = np.array(p, dtype=np.float64, copy=True)
    m_arr = np.asarray(masses, dtype=np.float64)

    sn_list = [tuple(int(i) for i in sn) for sn in (supernodes or [])]
    sn_bodies: Set[int] = set()
    for sn in sn_list:
        sn_bodies.update(sn)

    active_leaves: Set[Tuple[int, int]] = set()
    for i, j in leaf_pairs:
        if i not in sn_bodies and j not in sn_bodies:
            active_leaves.add((min(i, j), max(i, j)))

    # 1. Initial half apex kick (preserving inter-supernode cross-pairs)
    f_apex0 = compute_intercluster_forces(
        q_curr, m_arr, active_leaves, sn_bodies, G=G, supernode_clusters=sn_list
    )
    p_curr += 0.5 * dt * f_apex0

    # 2. Analytic leaf flows, synchronous supernodes, and free kinetic drift
    handled: Set[int] = set()
    for i, j in active_leaves:
        m_i, m_j = float(m_arr[i]), float(m_arr[j])
        M_ij = m_i + m_j
        mu_ij = (m_i * m_j) / M_ij
        R_cm = (m_i * q_curr[i] + m_j * q_curr[j]) / M_ij
        P_cm = p_curr[i] + p_curr[j]
        r0 = q_curr[i] - q_curr[j]
        w0 = (m_j * p_curr[i] - m_i * p_curr[j]) / M_ij

        r1, w1 = step_kepler_sundman(r0, w0, mu_ij, G * M_ij, dt)
        R_cm_next = R_cm + dt * (P_cm / M_ij)
        q_curr[i] = R_cm_next + (m_j / M_ij) * r1
        q_curr[j] = R_cm_next - (m_i / M_ij) * r1
        p_curr[i] = (m_i / M_ij) * P_cm + w1
        p_curr[j] = (m_j / M_ij) * P_cm - w1
        handled.add(i)
        handled.add(j)

    for sn in sn_list:
        q_curr, p_curr = _canonical_step_supernode_synchronous(
            q_curr, p_curr, m_arr, body_indices=sn, dt=dt, G=G, drift_center_of_mass=True
        )
        handled.update(sn)

    for k in range(len(m_arr)):
        if k not in handled:
            q_curr[k] += dt * (p_curr[k] / m_arr[k])

    # 3. Final half apex kick
    f_apex1 = compute_intercluster_forces(
        q_curr, m_arr, active_leaves, sn_bodies, G=G, supernode_clusters=sn_list
    )
    p_curr += 0.5 * dt * f_apex1
    return q_curr, p_curr


__all__ = [
    "Z_MIN_HYPERBOLIC",
    "stumpff_c",
    "solve_sundman_time_equation",
    "step_kepler_sundman",
    "KeplerSundmanSolver",
    "step_supernode_synchronous",
    "SupernodeSolver",
    "compute_intercluster_forces",
    "ScreenedApexPotential",
    "step_strang_subcycle",
    "StrangIntegrator",
    "Yoshida4Integrator",
    "BlanesMoan4Integrator",
    "BlanesMoan6Integrator",
]


