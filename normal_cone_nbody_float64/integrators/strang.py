"""
The 2nd-Order Palindromic Strang Splitting Propagator S_2.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.1 (The 7-Stage Palindromic Strang Propagator S_2),
  Equations (5.25)-(5.42).
- Formal Lean 4 Alignment: PaperUnifiedFoundations.lean
  (extended_symplectic_strang_order2, microscopic_time_reversibility).
- Symmetric operator splitting:
    S_2(h) = K_{h/2}^{apex} o Phi_h^{drift+leaf+super} o K_{h/2}^{apex}  (kick-exterior FSAL form)
    or Phi_{h/2} o K_h o Phi_{h/2}  (drift-exterior form).
"""

from __future__ import annotations
from typing import Optional, Set, Sequence, Tuple

from ..core.constants import G_DEFAULT
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree
from ..potential.screened_apex import ScreenedApexPotential
from .protocol import SymplecticIntegrator
from .common import (
    step_drift_and_kepler,
    step_apex_kick,
    set_current_integrator_stages,
    _CURRENT_INTEGRATOR_STAGES,
)


class StrangIntegrator(SymplecticIntegrator):
    """
    Second-order palindromic Strang splitting propagator S_2.
    """

    def __init__(
        self,
        G: float = G_DEFAULT,
        screening_mode: str = "c2",
        form: str = "kick_exterior",
    ) -> None:
        self.G = G
        self.screening_mode = screening_mode
        self.form = form
        self._default_apex = ScreenedApexPotential(G=G, screening_mode=screening_mode, order=2)

    @property
    def name(self) -> str:
        return "Strang-S2"

    @property
    def order(self) -> int:
        return 2

    @property
    def effective_order(self) -> int:
        return 2

    @property
    def evaluations_per_step(self) -> int:
        return 2

    @property
    def num_force_evaluations_per_step(self) -> int:
        return self.evaluations_per_step

    @property
    def fsal_evaluations_per_step(self) -> int:
        return 1

    @property
    def isolated_evaluations_per_step(self) -> int:
        return 2 if self.form == "kick_exterior" else 1

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
        Executes a single second-order Strang splitting macro-step of size dt.
        """
        if dt == 0.0:
            return state

        eval_apex = apex or self._default_apex

        if self.form == "kick_exterior":
            stages = ((1.0,), (0.5, 0.5))
        else:
            stages = ((0.5, 0.5), (1.0,))

        token = set_current_integrator_stages(stages)
        try:
            if self.form == "kick_exterior":
                # Kick-exterior FSAL form: K_{dt/2} o Phi_{dt} o K_{dt/2}
                s1 = step_apex_kick(
                    state, tree, eval_apex, dt * 0.5, supernode_bodies, supernode_clusters=supernode_clusters
                )
                s2 = step_drift_and_kepler(
                    s1, tree, dt, G=self.G, supernode_bodies=supernode_bodies, supernode_clusters=supernode_clusters
                )
                s3 = step_apex_kick(
                    s2, tree, eval_apex, dt * 0.5, supernode_bodies, supernode_clusters=supernode_clusters
                )
                return s3
            else:
                # Drift-exterior form: Phi_{dt/2} o K_{dt} o Phi_{dt/2}
                s1 = step_drift_and_kepler(
                    state,
                    tree,
                    dt * 0.5,
                    G=self.G,
                    supernode_bodies=supernode_bodies,
                    supernode_clusters=supernode_clusters,
                )
                s2 = step_apex_kick(
                    s1, tree, eval_apex, dt, supernode_bodies, supernode_clusters=supernode_clusters
                )
                s3 = step_drift_and_kepler(
                    s2,
                    tree,
                    dt * 0.5,
                    G=self.G,
                    supernode_bodies=supernode_bodies,
                    supernode_clusters=supernode_clusters,
                )
                return s3
        finally:
            _CURRENT_INTEGRATOR_STAGES.reset(token)
