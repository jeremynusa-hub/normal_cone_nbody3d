"""
Direct Unsoftened Newtonian Gravitational Potential and Baseline Force Evaluator.

Theoretical Foundations:
- Monograph Reference: Chapter 1 (Unsoftened N-Body Problem on Configuration Space)
  and Section 4.1 (Complete Energy Decomposition).
- Strict Mathematical Constraints:
    - Pure unsoftened Coulomb/Newton kernel: V(q) = - sum_{i < j} G m_i m_j / ||q_i - q_j||.
    - Zero softening epsilon = 0 and zero ad-hoc distance thresholds d_0 = 0.
- Invariance Properties:
    - Translational invariance: sum_i F_i = 0 (< 1e-15).
    - Rotational invariance: sum_i q_i x F_i = 0 (< 1e-14).
    - True ground truth for total energy H = T + V without screening approximations.
"""

from __future__ import annotations
from typing import Tuple
import numpy as np

from ..core.constants import GRAVITATIONAL_CONSTANT
from ..core.state import ExtendedPhaseState


def compute_direct_newtonian_potential(
    q: np.ndarray,
    masses: np.ndarray,
    G: float = GRAVITATIONAL_CONSTANT,
) -> float:
    """
    Computes the unsoftened total Newtonian potential energy:
        V(q) = - sum_{1 <= i < j <= N} G * m_i * m_j / ||q_i - q_j||.
    """
    N = len(masses)
    if N <= 1:
        return 0.0

    diff = q[None, :, :] - q[:, None, :]  # shape (N, N, 3)
    dist_sq = np.sum(diff ** 2, axis=-1)
    np.fill_diagonal(dist_sq, 1.0)
    dist = np.sqrt(dist_sq)

    inv_dist = 1.0 / dist
    np.fill_diagonal(inv_dist, 0.0)

    m_matrix = G * (masses[:, None] * masses[None, :])
    return -0.5 * float(np.sum(m_matrix * inv_dist))


def compute_direct_newtonian_forces(
    q: np.ndarray,
    masses: np.ndarray,
    G: float = GRAVITATIONAL_CONSTANT,
) -> np.ndarray:
    """
    Computes the unsoftened Cartesian gravitational force acting on each body:
        F_i = sum_{j != i} G * m_i * m_j * (q_j - q_i) / ||q_j - q_i||^3.

    Returns:
        forces: NumPy ndarray of shape (N, 3).
    """
    N = len(masses)
    if N <= 1:
        return np.zeros((N, 3), dtype=np.float64)

    diff = q[None, :, :] - q[:, None, :]  # shape (N, N, 3), diff[i, j] = q_j - q_i
    dist_sq = np.sum(diff ** 2, axis=-1)
    np.fill_diagonal(dist_sq, 1.0)
    dist = np.sqrt(dist_sq)

    inv_dist3 = 1.0 / (dist_sq * dist)
    np.fill_diagonal(inv_dist3, 0.0)

    m_matrix = G * (masses[:, None] * masses[None, :])
    f_matrix = m_matrix * inv_dist3  # shape (N, N)

    return np.sum(f_matrix[:, :, None] * diff, axis=1)


def direct_newtonian_step(
    state: ExtendedPhaseState,
    dt: float,
    G: float = GRAVITATIONAL_CONSTANT,
) -> ExtendedPhaseState:
    """
    Executes a single Strang leapfrog (Drift-Kick-Drift) step on unsoftened direct Newtonian equations.
    Serves as an exact baseline comparison for non-singular orbits.
    """
    if dt == 0.0:
        return state

    half_dt = 0.5 * dt
    # First half drift
    v_curr = state.p / state.masses[:, None]
    q_mid = state.q + half_dt * v_curr

    # Full kick at midpoint
    forces = compute_direct_newtonian_forces(q_mid, state.masses, G=G)
    p_next = state.p + dt * forces

    # Second half drift
    v_next = p_next / state.masses[:, None]
    q_next = q_mid + half_dt * v_next

    return state.with_updates(q=q_next, p=p_next, t=state.t + dt)
