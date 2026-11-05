"""
Canonical Extended Phase State on the Manifold P_ext = T*Q[N] x R_t x R_E.

Mathematical Monograph Alignment:
- Chapter 1: Section 1.2 — Extended Poincaré--Sundman phase space (P_ext, Omega_ext).
- Chapter 2: Section 2.1, Equations (2.1)-(2.4) — Zero-energy hypersurface Sigma_0 = {H - E = 0}.
- Chapter 4: Section 4.1, Equations (4.1)-(4.3) — Fulton--MacPherson configuration manifold Q[N].
- Certified Lean 4 Alignment: PaperUnifiedFoundations.lean (PhaseState, TangentVector3D).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Union, Any
import numpy as np
import mpmath as mp

from .constants import (
    to_f128,
    to_f64,
)


@dataclass(frozen=True)
class ExtendedPhaseState:
    """
    Immutable representation of an instantaneous physical state in the extended
    Poincaré--Sundman phase space:
        P_ext = T*Q[N] x R_t x R_E
    equipped with the canonical symplectic 2-form:
        Omega_ext = Omega_0 + dE ^ dt = sum_i dq_i ^ dp_i + dE ^ dt.

    Supports both IEEE-754 double precision (np.float64) and extended quadruple
    precision / multi-precision (mpmath 35 dps object arrays).

    Attributes:
        q: Cartesian positions in R^3, shape (N, 3), dtype float64 or object (mpmath.mpf).
        p: Conjugate linear momenta in R^3, shape (N, 3), dtype float64 or object.
        masses: Strictly positive point masses m_i > 0, shape (N,).
        t: Coordinate physical time t in R (float or mpmath.mpf).
        E: Extended conjugate energy variable E in R, satisfying {t, E} = -1.
        E0: Conserved initial energy baseline H(q_0, p_0) defining the zero-energy surface Sigma_0.
    """
    q: np.ndarray
    p: np.ndarray
    masses: np.ndarray
    t: Union[float, mp.mpf] = 0.0
    E: Union[float, mp.mpf] = 0.0
    E0: Union[float, mp.mpf] = 0.0
    skip_validation: bool = field(default=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        # Enforce read-only arrays for mathematical immutability
        if isinstance(self.q, np.ndarray) and self.q.flags.writeable:
            self.q.flags.writeable = False
        if isinstance(self.p, np.ndarray) and self.p.flags.writeable:
            self.p.flags.writeable = False
        if isinstance(self.masses, np.ndarray) and self.masses.flags.writeable:
            self.masses.flags.writeable = False

        if not self.skip_validation:
            self.validate()

        # Synchronize extended energy baseline if not provided
        if self.E == 0.0 and self.E0 != 0.0:
            object.__setattr__(self, "E", self.E0)

    @property
    def is_f128(self) -> bool:
        """Returns True if the state arrays store high-precision mpmath.mpf elements."""
        return self.q.dtype == object or (self.q.size > 0 and isinstance(self.q.flat[0], mp.mpf))

    def validate(self) -> None:
        """Validates dimensional consistency, data types, and positive mass bounds."""
        if not (isinstance(self.q, np.ndarray) and isinstance(self.p, np.ndarray) and isinstance(self.masses, np.ndarray)):
            raise TypeError("q, p, and masses must be numpy.ndarray instances.")

        if self.q.ndim != 2 or self.q.shape[1] != 3:
            raise ValueError(f"q must have shape (N, 3), got {self.q.shape}.")
        if self.p.ndim != 2 or self.p.shape[1] != 3:
            raise ValueError(f"p must have shape (N, 3), got {self.p.shape}.")

        n_bodies = self.q.shape[0]
        if n_bodies < 2:
            raise ValueError(f"N-body system requires at least N >= 2 bodies, got N={n_bodies}.")

        if self.p.shape[0] != n_bodies:
            raise ValueError(f"Dimension mismatch between q ({n_bodies}) and p ({self.p.shape[0]}).")

        if self.masses.ndim != 1 or self.masses.shape[0] != n_bodies:
            raise ValueError(f"masses must have shape ({n_bodies},), got {self.masses.shape}.")

        if self.is_f128:
            if any(m <= mp.mpf('0.0') for m in self.masses):
                raise ValueError("All point masses must be strictly positive (m_i > 0).")
        else:
            if np.any(self.masses <= 0.0):
                raise ValueError("All point masses must be strictly positive (m_i > 0).")

    @property
    def N(self) -> int:
        """Total number of point masses N in the system."""
        return int(self.masses.shape[0])

    @property
    def total_mass(self) -> Union[float, mp.mpf]:
        """Total gravitational mass M_tot = sum_{i=1}^N m_i."""
        tot = np.sum(self.masses)
        return tot if self.is_f128 else float(tot)

    @property
    def total_momentum(self) -> np.ndarray:
        """Total linear momentum vector P_tot = sum_{i=1}^N p_i in R^3."""
        return np.sum(self.p, axis=0)

    @property
    def center_of_mass(self) -> np.ndarray:
        """Barycentric center of mass R_cm = (1 / M_tot) * sum_{i=1}^N m_i q_i in R^3."""
        return np.sum(self.q * self.masses[:, np.newaxis], axis=0) / self.total_mass

    @property
    def center_of_mass_velocity(self) -> np.ndarray:
        """Velocity of the center of mass V_cm = P_tot / M_tot in R^3."""
        return self.total_momentum / self.total_mass

    @property
    def velocities(self) -> np.ndarray:
        """Physical velocities v_i = p_i / m_i of shape (N, 3)."""
        return self.p / self.masses[:, np.newaxis]

    def with_updates(
        self,
        q: Optional[np.ndarray] = None,
        p: Optional[np.ndarray] = None,
        masses: Optional[np.ndarray] = None,
        t: Optional[Union[float, mp.mpf]] = None,
        E: Optional[Union[float, mp.mpf]] = None,
        E0: Optional[Union[float, mp.mpf]] = None,
        skip_validation: bool = False,
    ) -> ExtendedPhaseState:
        """Functional immutable update returning a new ExtendedPhaseState instance."""
        target_dtype = self.q.dtype
        new_q = self.q if q is None else np.asarray(q, dtype=target_dtype)
        new_p = self.p if p is None else np.asarray(p, dtype=target_dtype)
        new_masses = self.masses if masses is None else np.asarray(masses, dtype=self.masses.dtype)
        new_t = self.t if t is None else (to_f128(t) if self.is_f128 else float(t))
        new_E = self.E if E is None else (to_f128(E) if self.is_f128 else float(E))
        new_E0 = self.E0 if E0 is None else (to_f128(E0) if self.is_f128 else float(E0))

        return ExtendedPhaseState(
            q=new_q,
            p=new_p,
            masses=new_masses,
            t=new_t,
            E=new_E,
            E0=new_E0,
            skip_validation=skip_validation,
        )

    def to_f128(self) -> ExtendedPhaseState:
        """Converts entire state to extended quadruple precision / mpmath object arrays."""
        if self.is_f128:
            return self
        return ExtendedPhaseState(
            q=to_f128(self.q),
            p=to_f128(self.p),
            masses=to_f128(self.masses),
            t=to_f128(self.t),
            E=to_f128(self.E),
            E0=to_f128(self.E0),
            skip_validation=self.skip_validation,
        )

    def to_f64(self) -> ExtendedPhaseState:
        """Converts state back to standard IEEE-754 float64 arrays."""
        if not self.is_f128:
            return self
        return ExtendedPhaseState(
            q=to_f64(self.q),
            p=to_f64(self.p),
            masses=to_f64(self.masses),
            t=to_f64(self.t),
            E=to_f64(self.E),
            E0=to_f64(self.E0),
            skip_validation=self.skip_validation,
        )

    def to_barycentric(self) -> ExtendedPhaseState:
        """
        Projects state onto the relative center-of-mass frame:
            q_bary = q - R_cm,  p_bary = p - (m_i / M_tot) * P_tot.
        Guarantees R_cm == 0 and P_tot == 0 to machine precision (< 1e-32 in f128).
        """
        r_cm = self.center_of_mass
        p_tot = self.total_momentum
        q_bary = self.q - r_cm
        p_bary = self.p - (self.masses[:, np.newaxis] / self.total_mass) * p_tot
        return self.with_updates(q=q_bary, p=p_bary)

    def min_pairwise_distance(self) -> Union[float, mp.mpf]:
        """Returns the minimal Euclidean distance min_{i < j} ||q_i - q_j||."""
        N = self.N
        if N < 2:
            return mp.mpf("inf") if self.is_f128 else float("inf")
        if self.is_f128:
            # High-precision Euclidean pairwise distance
            dists = []
            for i in range(N):
                for j in range(i + 1, N):
                    d_vec = self.q[i] - self.q[j]
                    d = mp.sqrt(d_vec[0]**2 + d_vec[1]**2 + d_vec[2]**2)
                    dists.append(d)
            return min(dists)
        else:
            diff = self.q[:, np.newaxis, :] - self.q[np.newaxis, :, :]
            dist_matrix = np.linalg.norm(diff, axis=-1)
            np.fill_diagonal(dist_matrix, np.inf)
            return float(np.min(dist_matrix))

    def is_collision_free(self, tol: Optional[Union[float, mp.mpf]] = None) -> bool:
        """
        Predicate: regular_iff_not_collision.
        Returns True iff min_{i < j} ||q_i - q_j|| > tol (or >= tol if tol > 0).
        By default (tol=None or tol=0), returns True iff min_{i < j} ||q_i - q_j|| > 0.
        """
        if tol is None or tol == 0.0:
            return bool(self.min_pairwise_distance() > 0.0)
        return bool(self.min_pairwise_distance() >= tol)
