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
from typing import Tuple, Union
import numpy as np
import mpmath as mp

from ..core.constants import GRAVITATIONAL_CONSTANT, DEFAULT_DPS, to_f128
from ..core.state import ExtendedPhaseState

mp.mp.dps = DEFAULT_DPS


def _is_f128_inputs(q: np.ndarray, masses: np.ndarray, G: Union[float, mp.mpf]) -> bool:
    if isinstance(G, mp.mpf):
        return True
    if getattr(q, "dtype", None) == object or getattr(masses, "dtype", None) == object:
        return True
    return False


def compute_direct_newtonian_potential(
    q: np.ndarray,
    masses: np.ndarray,
    G: Union[float, mp.mpf] = GRAVITATIONAL_CONSTANT,
) -> Union[float, mp.mpf]:
    """
    Computes the unsoftened total Newtonian potential energy:
        V(q) = - sum_{1 <= i < j <= N} G * m_i * m_j / ||q_i - q_j||.
    """
    N = len(masses)
    is_f128 = _is_f128_inputs(q, masses, G)
    if N <= 1:
        return mp.mpf('0.0') if is_f128 else 0.0

    if is_f128:
        G_mp = to_f128(G)
        V_tot = mp.mpf('0.0')
        for i in range(N):
            m_i = to_f128(masses[i])
            for j in range(i + 1, N):
                dx = to_f128(q[j, 0]) - to_f128(q[i, 0])
                dy = to_f128(q[j, 1]) - to_f128(q[i, 1])
                dz = to_f128(q[j, 2]) - to_f128(q[i, 2])
                r_sq = dx * dx + dy * dy + dz * dz
                if r_sq == mp.mpf('0.0'):
                    return -mp.inf
                r = mp.sqrt(r_sq)
                V_tot -= G_mp * m_i * to_f128(masses[j]) / r
        return V_tot

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
    G: Union[float, mp.mpf] = GRAVITATIONAL_CONSTANT,
) -> np.ndarray:
    """
    Computes the unsoftened Cartesian gravitational force acting on each body:
        F_i = sum_{j != i} G * m_i * m_j * (q_j - q_i) / ||q_j - q_i||^3.

    Returns:
        forces: NumPy ndarray of shape (N, 3).
    """
    N = len(masses)
    is_f128 = _is_f128_inputs(q, masses, G)
    if N <= 1:
        if is_f128:
            return np.full((N, 3), mp.mpf('0.0'), dtype=object)
        return np.zeros((N, 3), dtype=np.float64)

    if is_f128:
        G_mp = to_f128(G)
        forces = np.full((N, 3), mp.mpf('0.0'), dtype=object)
        for i in range(N):
            m_i = to_f128(masses[i])
            for j in range(i + 1, N):
                dx = to_f128(q[j, 0]) - to_f128(q[i, 0])
                dy = to_f128(q[j, 1]) - to_f128(q[i, 1])
                dz = to_f128(q[j, 2]) - to_f128(q[i, 2])
                r_sq = dx * dx + dy * dy + dz * dz
                r = mp.sqrt(r_sq)
                f_mag = G_mp * m_i * to_f128(masses[j]) / (r_sq * r)
                fx, fy, fz = f_mag * dx, f_mag * dy, f_mag * dz
                forces[i, 0] += fx
                forces[i, 1] += fy
                forces[i, 2] += fz
                forces[j, 0] -= fx
                forces[j, 1] -= fy
                forces[j, 2] -= fz
        return forces

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
    dt: Union[float, mp.mpf],
    G: Union[float, mp.mpf] = GRAVITATIONAL_CONSTANT,
) -> ExtendedPhaseState:
    """
    Executes a single Strang leapfrog (Drift-Kick-Drift) step on unsoftened direct Newtonian equations.
    Serves as an exact baseline comparison for non-singular orbits.
    """
    if dt == 0.0:
        return state

    if state.is_f128:
        dt_mp = to_f128(dt)
        half_dt = mp.mpf('0.5') * dt_mp
        G_mp = to_f128(G)
        v_curr = state.p / state.masses[:, None]
        q_mid = state.q + half_dt * v_curr
        forces = compute_direct_newtonian_forces(q_mid, state.masses, G=G_mp)
        p_next = state.p + dt_mp * forces
        v_next = p_next / state.masses[:, None]
        q_next = q_mid + half_dt * v_next
        return state.with_updates(q=q_next, p=p_next, t=state.t + dt_mp)

    half_dt = 0.5 * float(dt)
    # First half drift
    v_curr = state.p / state.masses[:, None]
    q_mid = state.q + half_dt * v_curr

    # Full kick at midpoint
    forces = compute_direct_newtonian_forces(q_mid, state.masses, G=G)
    p_next = state.p + float(dt) * forces

    # Second half drift
    v_next = p_next / state.masses[:, None]
    q_next = q_mid + half_dt * v_next

    return state.with_updates(q=q_next, p=p_next, t=state.t + float(dt))

