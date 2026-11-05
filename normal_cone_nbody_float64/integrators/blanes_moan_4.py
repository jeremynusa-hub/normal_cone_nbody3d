"""
The Optimal 7-Stage 4th-Order Palindromic Strang Composition BM_4^{(7)}.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.3 (High-Order Symmetric Compositions and Spectral Optimization),
  Definition 5.1, Equations (5.25)-(5.30).
- Formal Lean 4 Alignment: PaperUnifiedFoundations.lean &
  01_BlanesMoanShadowConfinement_Theorems.lean (blanes_moan_shadow_confinement_theorem).
- Palindromic 7-stage symmetric composition of kick-exterior Strang steps (with FSAL
  half-kick fusion b_j = (a_{j-1} + a_j)/2 across consecutive Strang sub-steps):
    S_4^{BM}(h) = S_2^*(a_1 h) o S_2^*(a_2 h) o S_2^*(a_3 h) o S_2^*(a_4 h) o S_2^*(a_3 h) o S_2^*(a_2 h) o S_2^*(a_1 h),
    where 2*(a_1 + a_2 + a_3) + a_4 = 1, 2*(a_1^3 + a_2^3 + a_3^3) + a_4^3 = 0,
    and |C_5| <= 0.01236 < (1/400) |C_5^{Yoshida}|.
"""

from __future__ import annotations
from typing import Optional, Set, Sequence, Tuple
import numpy as np

from ..core.constants import (
    G_DEFAULT,
    BM4_A1,
    BM4_A2,
    BM4_A3,
    BM4_A4,
    BM4_DRIFT_WEIGHTS,
    BM4_IMPULSE_WEIGHTS,
)
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree
from ..potential.screened_apex import ScreenedApexPotential
from .protocol import SymplecticIntegrator
from .strang import StrangIntegrator
from .common import (
    step_drift_and_kepler,
    step_apex_kick,
    set_current_integrator_stages,
    _CURRENT_INTEGRATOR_STAGES,
)

# Certified 4th-order 7-stage Blanes--Moan coefficients (imported from core.constants)
BM4_WEIGHTS: Tuple[float, ...] = BM4_DRIFT_WEIGHTS


class BlanesMoan4Integrator(SymplecticIntegrator):
    """
    Optimal 7-stage fourth-order palindromic Blanes--Moan scheme BM_4^{(7)}.

    Implements FSAL (First Same As Last) half-kick fusion, reducing force evaluations
    from 14 to 8 per isolated macro-step (7 effective under continuous chaining).
    """

    def __init__(
        self,
        G: float = G_DEFAULT,
        screening_mode: str = "c4",
    ) -> None:
        self.G = G
        self.screening_mode = screening_mode
        self._strang = StrangIntegrator(G=G, screening_mode=screening_mode, form="kick_exterior")
        self._default_apex = ScreenedApexPotential(G=G, screening_mode=screening_mode, order=4)
        self.drift_sequence = BM4_DRIFT_WEIGHTS
        self.impulse_sequence = BM4_IMPULSE_WEIGHTS

    @property
    def name(self) -> str:
        return "BlanesMoan-BM4"

    @property
    def order(self) -> int:
        return 4

    @property
    def effective_order(self) -> int:
        return 4

    @property
    def evaluations_per_step(self) -> int:
        # 8 force evaluations per isolated step, 7 under continuous FSAL chaining
        return 7

    @property
    def num_force_evaluations_per_step(self) -> int:
        return self.evaluations_per_step

    @property
    def fsal_evaluations_per_step(self) -> int:
        return 7

    @property
    def isolated_evaluations_per_step(self) -> int:
        return 8

    @property
    def is_palindromic(self) -> bool:
        return True

    def step(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        dt: float,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """
        Executes a 4th-order BM_4^{(7)} macro-step of size dt via 7 drift stages
        and 8 fused impulse (kick) evaluations.
        """
        if dt == 0.0:
            return state

        eval_apex = apex or self._default_apex
        curr_state = state

        token = set_current_integrator_stages((self.drift_sequence, self.impulse_sequence))
        try:
            # Initial half-kick K_{b_1 dt}
            curr_state = step_apex_kick(
                curr_state,
                tree,
                eval_apex,
                self.impulse_sequence[0] * dt,
                supernode_bodies,
                supernode_clusters=supernode_clusters,
            )

            # 7 drift stages interleaved with 6 fused kicks and 1 final half-kick
            for j in range(7):
                curr_state = step_drift_and_kepler(
                    curr_state,
                    tree,
                    self.drift_sequence[j] * dt,
                    G=self.G,
                    supernode_bodies=supernode_bodies,
                    supernode_clusters=supernode_clusters,
                )
                curr_state = step_apex_kick(
                    curr_state,
                    tree,
                    eval_apex,
                    self.impulse_sequence[j + 1] * dt,
                    supernode_bodies,
                    supernode_clusters=supernode_clusters,
                )

            return curr_state
        finally:
            _CURRENT_INTEGRATOR_STAGES.reset(token)
