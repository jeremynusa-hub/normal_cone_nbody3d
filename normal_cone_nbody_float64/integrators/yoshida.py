"""
The 4th-Order 3-Stage Symmetric Yoshida Composition S_4.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.2 (The Sheng--Suzuki Barrier and High-Order Symmetric Compositions),
  Equations (5.43)-(5.60).
- Formal Lean 4 Alignment: PaperUnifiedFoundations.lean
  (yoshida_order4_symmetric_composition, microscopic_time_reversibility).
- 3-stage symmetric composition of Strang steps:
    S_4(h) = S_2(w_1 h) o S_2(w_0 h) o S_2(w_1 h),
    where w_1 = 1 / (2 - 2^{1/3}) ~= 1.3512071919596576,
          w_0 = 1 - 2 w_1 = - 2^{1/3} / (2 - 2^{1/3}) ~= -1.7024143839193153.
"""

from __future__ import annotations
from typing import Optional, Set, Sequence, Tuple
import numpy as np

from ..core.constants import G_DEFAULT
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

# Exact algebraic weights for Yoshida 4th-order 3-stage composition
_CBRT_2 = 2.0 ** (1.0 / 3.0)
YOSHIDA_W1: float = 1.0 / (2.0 - _CBRT_2)
YOSHIDA_W0: float = 1.0 - 2.0 * YOSHIDA_W1
YOSHIDA_WEIGHTS: Tuple[float, float, float] = (YOSHIDA_W1, YOSHIDA_W0, YOSHIDA_W1)

# 3 drift stages for the palindromic composition
YOSHIDA_DRIFT_WEIGHTS: Tuple[float, float, float] = (YOSHIDA_W1, YOSHIDA_W0, YOSHIDA_W1)

# 4 fused impulse (kick) stages under FSAL (First Same As Last) half-kick fusion:
# The adjacent half-kicks K(w_j/2) o K(w_{j+1}/2) between consecutive Strang stages
# fuse into K((w_j + w_{j+1})/2), reducing force evaluations from 6 to 4 per isolated macro-step.
YOSHIDA_IMPULSE_WEIGHTS: Tuple[float, float, float, float] = (
    0.5 * YOSHIDA_W1,
    0.5 * (YOSHIDA_W1 + YOSHIDA_W0),
    0.5 * (YOSHIDA_W0 + YOSHIDA_W1),
    0.5 * YOSHIDA_W1,
)


class Yoshida4Integrator(SymplecticIntegrator):
    """
    Fourth-order symmetric Yoshida 3-stage palindromic composition S_4.

    Implements FSAL (First Same As Last) half-kick fusion, reducing force evaluations
    from 6 to 4 per isolated macro-step (3 effective under continuous chaining).
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
        self.drift_sequence = YOSHIDA_DRIFT_WEIGHTS
        self.impulse_sequence = YOSHIDA_IMPULSE_WEIGHTS

    @property
    def name(self) -> str:
        return "Yoshida-S4"

    @property
    def order(self) -> int:
        return 4

    @property
    def effective_order(self) -> int:
        return 4

    @property
    def evaluations_per_step(self) -> int:
        return 4

    @property
    def num_force_evaluations_per_step(self) -> int:
        return self.evaluations_per_step

    @property
    def fsal_evaluations_per_step(self) -> int:
        return 3

    @property
    def isolated_evaluations_per_step(self) -> int:
        return 4

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
        Executes a 4th-order Yoshida macro-step of size dt via 3 drift stages
        and 4 fused impulse (kick) evaluations.
        """
        if dt == 0.0:
            return state

        eval_apex = apex or self._default_apex
        curr_state = state

        token = set_current_integrator_stages((self.drift_sequence, self.impulse_sequence))
        try:
            # Initial half-kick K_{0.5 * w_1 * dt}
            curr_state = step_apex_kick(
                curr_state,
                tree,
                eval_apex,
                self.impulse_sequence[0] * dt,
                supernode_bodies,
                supernode_clusters=supernode_clusters,
            )

            # 3 drift stages interleaved with 2 fused kicks and 1 final half-kick
            for j in range(3):
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
