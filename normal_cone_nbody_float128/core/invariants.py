r"""
Exact Computation and Verification of the 10 Lie--Bargmann Generators.

Mathematical Monograph Alignment:
- Chapter 1: Section 1.2 & Section 1.4 — Extended Galilean phase space P_ext and central Bargmann extension.
- Chapter 2: Section 2.3, Equations (2.12)-(2.18) — The 10-dimensional extended Galilean
  Lie--Bargmann algebra (P in R^3, L in R^3, G in R^3, H_ext in R).
- Chapter 5: Section 5.1, Equations (5.4)-(5.8) — Symplectic invariance and exact discrete
  Lie--Bargmann conservation under palindromic propagators.
- Certified Lean 4 Alignment: PaperUnifiedFoundations.lean
  (lie_bargmann_closure_dim10, BargmannLinearMomentum, BargmannAngularMomentum,
   BargmannGalileanBoost, BargmannEnergyExtended).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Union, Optional
import numpy as np
import mpmath as mp

from .state import ExtendedPhaseState
from .constants import G_DEFAULT, G_DEFAULT_F128, to_f128, to_f64


@dataclass(frozen=True)
class LieBargmannInvariants:
    """
    Summary container for the 10 canonical generators of the extended Galilean
    Lie--Bargmann group in the extended Poincaré--Sundman phase space:
    - Linear momentum P in R^3 (3 invariants)
    - Angular momentum L in R^3 (3 invariants)
    - Galilean boost generator G in R^3 (3 invariants)
    - Extended Hamiltonian energy E in R (1 invariant)

    Auxiliary fields:
    - Instantaneous Cartesian Hamiltonian H(q, p)
    - Energy residual delta_E = |H(q, p) - E0|
    - Virial-normalized relative energy residual delta_E_rel = |H(q, p) - E0| / max(|E0|, T + |V|, 1e-30)
    - Center of mass R_cm in R^3
    - Center of mass velocity V_cm in R^3
    - Coordinate time t
    """
    P: np.ndarray                            # shape (3,), float64 or object (mp.mpf)
    L: np.ndarray                            # shape (3,), float64 or object
    G: np.ndarray                            # shape (3,), float64 or object
    H: Union[float, mp.mpf]                  # Instantaneous Hamiltonian H = T + V
    E: Union[float, mp.mpf]                  # Extended energy variable E(t)
    delta_E: Union[float, mp.mpf]            # Hamiltonian energy difference |H - E0|
    R_cm: np.ndarray                         # Center of mass, shape (3,)
    V_cm: np.ndarray                         # Center of mass velocity, shape (3,)
    t: Union[float, mp.mpf]                  # Physical coordinate time
    delta_E_rel: Union[float, mp.mpf] = 0.0  # Virial-normalized relative energy difference

    @property
    def is_f128(self) -> bool:
        """Returns True if the invariant arrays and scalars are in extended precision."""
        return isinstance(self.H, mp.mpf) or (self.P.dtype == object)

    def as_vector10(self) -> np.ndarray:
        """Returns the 10 Lie--Bargmann generators as a 1D array of shape (10,)."""
        dtype = object if self.is_f128 else np.float64
        return np.concatenate([self.P, self.L, self.G, [self.H]]).astype(dtype)


def compute_linear_momentum(state: ExtendedPhaseState) -> np.ndarray:
    """
    Computes total linear momentum:
        P_tot = sum_{i=1}^N p_i in R^3.
    Conserved exactly due to spatial translation invariance (Newton's third law).
    """
    return np.sum(state.p, axis=0)


def compute_angular_momentum(state: ExtendedPhaseState) -> np.ndarray:
    """
    Computes total angular momentum:
        L_tot = sum_{i=1}^N q_i x p_i in R^3.
    Conserved exactly due to SO(3) rotational invariance.
    """
    return np.sum(np.cross(state.q, state.p), axis=0)


def compute_center_of_mass(state: ExtendedPhaseState) -> np.ndarray:
    """Computes barycentric center of mass position R_cm in R^3."""
    return state.center_of_mass


def compute_galilean_boost(state: ExtendedPhaseState) -> np.ndarray:
    """
    Computes the Galilean center-of-mass boost generator:
        G = M_tot * R_cm - t * P_tot in R^3.
    Conserved strictly along Hamiltonian trajectories (dG / dt = M_tot V_cm - P_tot = 0).
    """
    return state.total_mass * state.center_of_mass - state.t * state.total_momentum


def compute_kinetic_energy(state: ExtendedPhaseState) -> Union[float, mp.mpf]:
    """
    Computes Cartesian kinetic energy:
        T = (1/2) * sum_{i=1}^N (||p_i||^2 / m_i).
    """
    if state.is_f128:
        T_tot = mp.mpf('0.0')
        for i in range(state.N):
            p_sq = state.p[i, 0]**2 + state.p[i, 1]**2 + state.p[i, 2]**2
            T_tot += p_sq / state.masses[i]
        return mp.mpf('0.5') * T_tot
    else:
        p_sq = np.sum(state.p ** 2, axis=1)
        return float(0.5 * np.sum(p_sq / state.masses))


def compute_newtonian_potential_energy(
    state: ExtendedPhaseState,
    G: Optional[Union[float, mp.mpf]] = None,
) -> Union[float, mp.mpf]:
    """
    Computes unsoftened Newtonian pairwise gravitational potential:
        V = - sum_{1 <= i < j <= N} (G * m_i * m_j / ||q_i - q_j||).
    Strictly unsoftened: epsilon = 0.0, d_0 = 0.0.
    """
    if state.is_f128:
        G_val = G_DEFAULT_F128 if G is None else to_f128(G)
        V_tot = mp.mpf('0.0')
        N = state.N
        for i in range(N):
            for j in range(i + 1, N):
                diff = state.q[i] - state.q[j]
                r = mp.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
                if r == mp.mpf('0.0'):
                    return -mp.inf
                V_tot -= (G_val * state.masses[i] * state.masses[j]) / r
        return V_tot
    else:
        G_val = G_DEFAULT if G is None else float(G)
        N = state.N
        diff = state.q[:, np.newaxis, :] - state.q[np.newaxis, :, :]  # (N, N, 3)
        dist = np.linalg.norm(diff, axis=-1)  # (N, N)

        # Compute pairwise products m_i * m_j
        mass_matrix = np.outer(state.masses, state.masses)

        # Strict upper triangle indices (i < j)
        i_idx, j_idx = np.triu_indices(N, k=1)
        dist_pairs = dist[i_idx, j_idx]

        if np.any(dist_pairs == 0.0):
            return -float("inf")

        potential = -G_val * np.sum(mass_matrix[i_idx, j_idx] / dist_pairs)
        return float(potential)


compute_potential_energy = compute_newtonian_potential_energy


def compute_hamiltonian(
    state: ExtendedPhaseState,
    G: Optional[Union[float, mp.mpf]] = None,
) -> Union[float, mp.mpf]:
    """Computes the total classical Hamiltonian H(q, p) = T(p) + V(q)."""
    return compute_kinetic_energy(state) + compute_newtonian_potential_energy(state, G)


def compute_cartesian_hamiltonian(
    q: np.ndarray,
    p: np.ndarray,
    masses: np.ndarray,
    G: Optional[Union[float, mp.mpf]] = None,
) -> Union[float, mp.mpf]:
    """Computes total classical Hamiltonian from raw Cartesian arrays."""
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    return compute_hamiltonian(state, G=G)


def compute_invariants(
    state: ExtendedPhaseState,
    G: Optional[Union[float, mp.mpf]] = None,
) -> LieBargmannInvariants:
    """
    Evaluates all 10 Lie--Bargmann generators and auxiliary kinematics for the given state.
    For parabolic orbits (E0 = 0.0), delta_E evaluates |H(q, p) - E0| = |H(q, p)| directly
    and delta_E_rel normalizes by the virial scale max(|E0|, T + |V|, 1e-30).
    """
    p_tot = compute_linear_momentum(state)
    l_tot = compute_angular_momentum(state)
    g_boost = compute_galilean_boost(state)
    t_kin = compute_kinetic_energy(state)
    v_pot = compute_newtonian_potential_energy(state, G)
    h_energy = t_kin + v_pot
    r_cm = state.center_of_mass
    v_cm = state.center_of_mass_velocity
    delta_e = abs(h_energy - state.E0)
    if state.is_f128:
        virial_scale = max(abs(state.E0), t_kin + abs(v_pot), mp.mpf('1e-38'))
        delta_e_rel = delta_e / virial_scale if mp.isfinite(virial_scale) else mp.mpf('0.0')
    else:
        virial_scale = max(abs(state.E0), t_kin + abs(v_pot), 1e-30)
        delta_e_rel = delta_e / virial_scale if np.isfinite(virial_scale) else 0.0

    return LieBargmannInvariants(
        P=p_tot,
        L=l_tot,
        G=g_boost,
        H=h_energy,
        E=state.E,
        delta_E=delta_e,
        R_cm=r_cm,
        V_cm=v_cm,
        t=state.t,
        delta_E_rel=delta_e_rel,
    )


def compute_truncation_residual(
    current: LieBargmannInvariants,
    initial: LieBargmannInvariants,
) -> np.ndarray:
    """
    Computes the 4-component Lie--Bargmann discrete truncation residual vector Delta_trunc in R^4_+:
        Delta_trunc = [ ||Delta P||_2, ||Delta L||_2, ||Delta G||_2, |Delta H| ] in R^4_+.
    Used for machine-precision invariant audit (< 10^-14 in float64, < 10^-32 in f128).
    """
    if current.is_f128 or initial.is_f128:
        dp = current.P - initial.P
        dl = current.L - initial.L
        dg = current.G - initial.G
        delta_p = mp.sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
        delta_l = mp.sqrt(dl[0]**2 + dl[1]**2 + dl[2]**2)
        delta_g = mp.sqrt(dg[0]**2 + dg[1]**2 + dg[2]**2)
        delta_h = abs(current.H - initial.H)
        return np.array([delta_p, delta_l, delta_g, delta_h], dtype=object)
    else:
        delta_p = np.linalg.norm(current.P - initial.P)
        delta_l = np.linalg.norm(current.L - initial.L)
        delta_g = np.linalg.norm(current.G - initial.G)
        delta_h = abs(current.H - initial.H)
        return np.array([delta_p, delta_l, delta_g, delta_h], dtype=np.float64)
