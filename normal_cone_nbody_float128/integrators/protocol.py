"""
Formal Protocol Definition and Data Containers for Symplectic Integrators.

Theoretical Foundations:
- Monograph Reference: Chapter 5 (Higher-Order Symplectic Integration),
  Sections 5.1–5.5.
- Formal Lean 4 Alignment: PaperUnifiedFoundations.lean
  (extended_symplectic_invariance, microscopic_time_reversibility).
- PEP 544 Runtime-Checkable Structural Subtyping for Zero-Overhead Polymorphism.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Protocol, runtime_checkable, Optional, Set, Sequence, Tuple, Union
import numpy as np
import mpmath as mp

from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree
from ..potential.screened_apex import ScreenedApexPotential


@dataclass(frozen=True)
class IntegrationStepResult:
    """
    Diagnostic container reporting step outcomes and numerical metrics:
    - state: updated ExtendedPhaseState
    - forces_evaluations: number of apex force evaluations performed in the step
    - dt: macro-step duration
    - wall_time_seconds: runtime spent in step
    """
    state: ExtendedPhaseState
    forces_evaluations: int
    dt: Union[float, mp.mpf]
    wall_time_seconds: float = 0.0


@runtime_checkable
class SymplecticIntegrator(Protocol):
    """
    Canonical structural subtyping protocol for all symplectic propagators
    acting on the extended cotangent bundle (P_ext, Omega_ext).
    """

    @property
    def name(self) -> str:
        """Human-readable identifier of the symplectic scheme."""
        ...

    @property
    def order(self) -> int:
        """Formal mathematical convergence order p."""
        ...

    @property
    def effective_order(self) -> int:
        """Effective convergence order p_eff under boundary processing."""
        ...

    @property
    def evaluations_per_step(self) -> int:
        """Nominal number of apex force evaluations per macro-step."""
        ...

    @property
    def num_force_evaluations_per_step(self) -> int:
        """Standardized alias for nominal apex force evaluations per macro-step."""
        ...

    @property
    def fsal_evaluations_per_step(self) -> int:
        """Effective number of apex force evaluations per macro-step under continuous FSAL chaining."""
        ...

    @property
    def isolated_evaluations_per_step(self) -> int:
        """Number of apex force evaluations in a single isolated macro-step (without boundary FSAL fusion)."""
        ...

    @property
    def is_palindromic(self) -> bool:
        """Whether the operator satisfies time-reversal symmetry R o Phi(-h) o R = Phi(h)^(-1)."""
        ...

    def step(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        dt: Union[float, mp.mpf],
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """
        Executes a single symplectic integration macro-step of size dt.

        Parameters:
            state: Current extended phase state (q, p, masses, t, E, E0).
            tree: Current active Jacobi binary tree topology.
            dt: Macro time step (can be positive or negative).
            apex: Optional ScreenedApexPotential evaluator (defaults to standard).
            supernode_bodies: Optional set of particle indices promoted to synchronous supernodes.
            supernode_clusters: Optional sequence of disjoint supernode body tuples S_k.

        Returns:
            Updated ExtendedPhaseState.
        """
        ...
