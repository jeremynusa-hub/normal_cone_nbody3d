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
from typing import Optional
import numpy as np



@dataclass(frozen=True)
class ExtendedPhaseState:
    """
    Immutable representation of an instantaneous physical state in the extended
    Poincaré--Sundman phase space:
        P_ext = T*Q[N] x R_t x R_E
    equipped with the canonical symplectic 2-form:
        Omega_ext = Omega_0 + dE ^ dt = sum_i dq_i ^ dp_i + dE ^ dt.

    Attributes:
        q: Cartesian positions in R^3, shape (N, 3), dtype float64.
        p: Conjugate linear momenta in R^3, shape (N, 3), dtype float64.
        masses: Strictly positive point masses m_i > 0, shape (N,), dtype float64.
        t: Coordinate physical time t in R.
        E: Extended conjugate energy variable E in R, satisfying {t, E} = -1.
        E0: Conserved initial energy baseline H(q_0, p_0) defining the zero-energy surface Sigma_0.
    """
    q: np.ndarray
    p: np.ndarray
    masses: np.ndarray
    t: float = 0.0
    E: float = 0.0
    E0: float = 0.0
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
            object.__setattr__(self, "E", float(self.E0))

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

        if np.any(self.masses <= 0.0):
            raise ValueError("All point masses must be strictly positive (m_i > 0).")

    @property
    def N(self) -> int:
        """Total number of point masses N in the system."""
        return int(self.masses.shape[0])

    @property
    def total_mass(self) -> float:
        """Total gravitational mass M_tot = sum_{i=1}^N m_i."""
        return float(np.sum(self.masses))

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
        t: Optional[float] = None,
        E: Optional[float] = None,
        E0: Optional[float] = None,
        skip_validation: bool = False,
    ) -> ExtendedPhaseState:
        """Functional immutable update returning a new ExtendedPhaseState instance."""
        new_q = self.q if q is None else np.asarray(q, dtype=np.float64)
        new_p = self.p if p is None else np.asarray(p, dtype=np.float64)
        new_masses = self.masses if masses is None else np.asarray(masses, dtype=np.float64)
        new_t = self.t if t is None else float(t)
        new_E = self.E if E is None else float(E)
        new_E0 = self.E0 if E0 is None else float(E0)

        return ExtendedPhaseState(
            q=new_q,
            p=new_p,
            masses=new_masses,
            t=new_t,
            E=new_E,
            E0=new_E0,
            skip_validation=skip_validation,
        )

    def to_barycentric(self) -> ExtendedPhaseState:
        """
        Projects state onto the relative center-of-mass frame:
            q_bary = q - R_cm,  p_bary = p - (m_i / M_tot) * P_tot.
        Guarantees R_cm == 0 and P_tot == 0 to machine precision.
        """
        r_cm = self.center_of_mass
        p_tot = self.total_momentum
        q_bary = self.q - r_cm
        p_bary = self.p - (self.masses[:, np.newaxis] / self.total_mass) * p_tot
        return self.with_updates(q=q_bary, p=p_bary)

    def min_pairwise_distance(self) -> float:
        """Returns the minimal Euclidean distance min_{i < j} ||q_i - q_j||."""
        if self.N < 2:
            return float("inf")
        diff = self.q[:, np.newaxis, :] - self.q[np.newaxis, :, :]
        dist_matrix = np.linalg.norm(diff, axis=-1)
        # Mask diagonal
        np.fill_diagonal(dist_matrix, np.inf)
        return float(np.min(dist_matrix))

    def is_collision_free(self, tol: float = 0.0) -> bool:
        """
        Predicate: regular_iff_not_collision.
        Returns True iff min_{i < j} ||q_i - q_j|| > tol (or >= tol if tol > 0).
        By default (tol=0.0), returns True iff min_{i < j} ||q_i - q_j|| > 0.
        """
        if tol == 0.0:
            return self.min_pairwise_distance() > 0.0
        return self.min_pairwise_distance() >= tol
