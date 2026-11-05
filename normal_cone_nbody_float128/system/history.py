"""
Comprehensive Trajectory Logging and Lie--Bargmann Invariant Tracking Ledger.

Mathematical Monograph Alignment:
- Chapter 2: Section 2.3 (The 10 Lie--Bargmann Generators in Extended Phase Space).
- Chapter 5: Section 5.1 (Discrete Conservation Laws and Backward Error Metrics).
- Provides compact in-memory trajectory caching, sample subsampling, vectorized NumPy
  accessors, and metric drift estimators.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Dict, Any
import numpy as np
import mpmath as mp

from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree
from ..core.invariants import LieBargmannInvariants, compute_invariants


@dataclass(frozen=True)
class HistorySnapshot:
    """
    Immutable phase state snapshot at a single integration checkpoint:
    - t: Coordinate time (float or mp.mpf).
    - q: Configuration array, shape (N, 3).
    - p: Momentum array, shape (N, 3).
    - invariants: Exact Lie--Bargmann invariants at checkpoint.
    - dyadic_level: Refinement level k utilized for the preceding step.
    - leaf_pairs: Active Jacobi tree leaf pairs at checkpoint.
    """
    t: Any
    q: np.ndarray
    p: np.ndarray
    invariants: LieBargmannInvariants
    dyadic_level: int = 0
    leaf_pairs: tuple = ()


class HistoryTracker:
    """
    Trajectory and dynamical invariant historical tracker for N-body simulations.
    Supports both float64 and quadruple precision (mpmath 35 dps) ExtendedPhaseState.
    """

    def __init__(
        self,
        record_interval: int = 1,
        max_records: Optional[int] = None,
        store_full_trajectories: bool = True,
    ) -> None:
        if record_interval < 1:
            raise ValueError(f"record_interval must be >= 1, got {record_interval}")

        self.record_interval = int(record_interval)
        self.max_records = max_records
        self.store_full_trajectories = bool(store_full_trajectories)

        self.snapshots: List[HistorySnapshot] = []
        self._step_counter: int = 0
        self._masses: Optional[np.ndarray] = None
        self.is_f128: bool = False

    def record(
        self,
        state: ExtendedPhaseState,
        invariants: Optional[LieBargmannInvariants] = None,
        dyadic_level: int = 0,
        force: bool = False,
        leaf_pairs: tuple = (),
    ) -> bool:
        """
        Records the given state snapshot if step matches recording cadence.
        Returns True if snapshot was appended.
        """
        if self._masses is None:
            self._masses = np.array(state.masses, copy=True, dtype=object if state.is_f128 else np.float64)
        if state.is_f128:
            self.is_f128 = True

        if not force and (self._step_counter % self.record_interval != 0):
            self._step_counter += 1
            return False

        self._step_counter += 1
        inv = invariants or compute_invariants(state)

        q_rec = np.copy(state.q) if self.store_full_trajectories else np.empty((0, 3), dtype=object if self.is_f128 else np.float64)
        p_rec = np.copy(state.p) if self.store_full_trajectories else np.empty((0, 3), dtype=object if self.is_f128 else np.float64)

        snapshot = HistorySnapshot(
            t=state.t,
            q=q_rec,
            p=p_rec,
            invariants=inv,
            dyadic_level=dyadic_level,
            leaf_pairs=tuple(leaf_pairs),
        )

        if self.max_records is not None and len(self.snapshots) >= self.max_records:
            # Circular overwrite / pop oldest
            self.snapshots.pop(0)

        self.snapshots.append(snapshot)
        return True

    def __len__(self) -> int:
        return len(self.snapshots)

    @property
    def times(self) -> np.ndarray:
        """Array of recorded time stamps of shape (S,)."""
        if self.is_f128:
            return np.array([s.t for s in self.snapshots], dtype=object)
        return np.array([s.t for s in self.snapshots], dtype=np.float64)

    @property
    def positions(self) -> np.ndarray:
        """Array of recorded positions of shape (S, N, 3)."""
        if not self.snapshots or not self.store_full_trajectories:
            return np.empty((0, 0, 3), dtype=object if self.is_f128 else np.float64)
        if self.is_f128:
            return np.array([s.q for s in self.snapshots], dtype=object)
        return np.array([s.q for s in self.snapshots], dtype=np.float64)

    @property
    def momenta(self) -> np.ndarray:
        """Array of recorded momenta of shape (S, N, 3)."""
        if not self.snapshots or not self.store_full_trajectories:
            return np.empty((0, 0, 3), dtype=object if self.is_f128 else np.float64)
        if self.is_f128:
            return np.array([s.p for s in self.snapshots], dtype=object)
        return np.array([s.p for s in self.snapshots], dtype=np.float64)

    @property
    def hamiltonians(self) -> np.ndarray:
        """Array of recorded Hamiltonian energy values of shape (S,)."""
        if self.is_f128:
            return np.array([s.invariants.H for s in self.snapshots], dtype=object)
        return np.array([s.invariants.H for s in self.snapshots], dtype=np.float64)

    @property
    def linear_momenta(self) -> np.ndarray:
        """Array of total linear momentum vectors of shape (S, 3)."""
        if self.is_f128:
            return np.array([s.invariants.P for s in self.snapshots], dtype=object)
        return np.array([s.invariants.P for s in self.snapshots], dtype=np.float64)

    @property
    def angular_momenta(self) -> np.ndarray:
        """Array of total angular momentum vectors of shape (S, 3)."""
        if self.is_f128:
            return np.array([s.invariants.L for s in self.snapshots], dtype=object)
        return np.array([s.invariants.L for s in self.snapshots], dtype=np.float64)

    @property
    def boosts(self) -> np.ndarray:
        """Array of Galilean boost vectors of shape (S, 3)."""
        if self.is_f128:
            return np.array([s.invariants.G for s in self.snapshots], dtype=object)
        return np.array([s.invariants.G for s in self.snapshots], dtype=np.float64)

    @property
    def times_f64(self) -> np.ndarray:
        """Array of recorded time stamps cast to float64."""
        return np.array([float(s.t) for s in self.snapshots], dtype=np.float64)

    @property
    def positions_f64(self) -> np.ndarray:
        """Array of recorded positions cast to float64 of shape (S, N, 3)."""
        if not self.snapshots or not self.store_full_trajectories:
            return np.empty((0, 0, 3), dtype=np.float64)
        return np.array([[[float(c) for c in body] for body in s.q] for s in self.snapshots], dtype=np.float64)

    @property
    def momenta_f64(self) -> np.ndarray:
        """Array of recorded momenta cast to float64 of shape (S, N, 3)."""
        if not self.snapshots or not self.store_full_trajectories:
            return np.empty((0, 0, 3), dtype=np.float64)
        return np.array([[[float(c) for c in body] for body in s.p] for s in self.snapshots], dtype=np.float64)

    @property
    def hamiltonians_f64(self) -> np.ndarray:
        """Array of recorded Hamiltonian energy values cast to float64 of shape (S,)."""
        return np.array([float(s.invariants.H) for s in self.snapshots], dtype=np.float64)

    @property
    def linear_momenta_f64(self) -> np.ndarray:
        """Array of total linear momentum vectors cast to float64 of shape (S, 3)."""
        return np.array([[float(c) for c in s.invariants.P] for s in self.snapshots], dtype=np.float64)

    @property
    def angular_momenta_f64(self) -> np.ndarray:
        """Array of total angular momentum vectors cast to float64 of shape (S, 3)."""
        return np.array([[float(c) for c in s.invariants.L] for s in self.snapshots], dtype=np.float64)

    @property
    def boosts_f64(self) -> np.ndarray:
        """Array of Galilean boost vectors cast to float64 of shape (S, 3)."""
        return np.array([[float(c) for c in s.invariants.G] for s in self.snapshots], dtype=np.float64)

    @property
    def dyadic_levels(self) -> np.ndarray:
        """Array of recorded dyadic refinement levels k of shape (S,)."""
        return np.array([s.dyadic_level for s in self.snapshots], dtype=np.int32)

    @property
    def masses(self) -> Optional[np.ndarray]:
        """Particle masses of shape (N,) if recorded."""
        return self._masses

    def get_times(self) -> np.ndarray:
        return self.times

    def get_positions(self) -> np.ndarray:
        return self.positions

    def get_momenta(self) -> np.ndarray:
        return self.momenta

    def get_energies(self) -> np.ndarray:
        return self.hamiltonians

    def get_linear_momenta(self) -> np.ndarray:
        return self.linear_momenta

    def get_angular_momenta(self) -> np.ndarray:
        return self.angular_momenta

    def get_center_of_masses(self) -> np.ndarray:
        return self.boosts

    def get_dyadic_levels(self) -> np.ndarray:
        return self.dyadic_levels

    def max_energy_drift(self) -> Any:
        """Maximum absolute Hamiltonian variation |H(t) - H(0)|."""
        if len(self.snapshots) <= 1:
            return mp.mpf(0) if self.is_f128 else 0.0
        if self.is_f128:
            H0 = self.snapshots[0].invariants.H
            return max(abs(s.invariants.H - H0) for s in self.snapshots)
        H = self.hamiltonians
        return float(np.max(np.abs(H - H[0])))

    def max_momentum_drift(self) -> Any:
        """Maximum absolute total momentum variation ||P(t) - P(0)||."""
        if len(self.snapshots) <= 1:
            return mp.mpf(0) if self.is_f128 else 0.0
        if self.is_f128:
            P0 = self.snapshots[0].invariants.P
            diffs = [
                mp.sqrt(sum((s.invariants.P[k] - P0[k]) ** 2 for k in range(3)))
                for s in self.snapshots
            ]
            return max(diffs)
        P = self.linear_momenta
        return float(np.max(np.linalg.norm(P - P[0], axis=1)))

    def max_angular_momentum_drift(self) -> Any:
        """Maximum absolute angular momentum variation ||L(t) - L(0)||."""
        if len(self.snapshots) <= 1:
            return mp.mpf(0) if self.is_f128 else 0.0
        if self.is_f128:
            L0 = self.snapshots[0].invariants.L
            diffs = [
                mp.sqrt(sum((s.invariants.L[k] - L0[k]) ** 2 for k in range(3)))
                for s in self.snapshots
            ]
            return max(diffs)
        L = self.angular_momenta
        return float(np.max(np.linalg.norm(L - L[0], axis=1)))

    def summary_metrics(self) -> Dict[str, Any]:
        """Returns a summary dictionary of dynamical conservation errors."""
        return {
            "total_checkpoints": float(len(self.snapshots)),
            "delta_H_max": self.max_energy_drift(),
            "delta_P_max": self.max_momentum_drift(),
            "delta_L_max": self.max_angular_momentum_drift(),
        }

    def clear(self) -> None:
        """Clears all logged history."""
        self.snapshots.clear()
        self._step_counter = 0
