"""
The Optimal 11-Stage 6th-Order Blanes--Moan Scheme BM_6^{(11)}.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.4 (The 11-Stage Sixth-Order Scheme BM6 and High-Precision Lie Integration),
  Equations (5.83)-(5.110), Table 1.
- Formal Lean 4 Alignment: PaperUnifiedFoundations.lean &
  BM6Order6SymmetricComposition.lean (bm6_order6_certified).
- 11-stage palindromic partitioned Runge--Kutta--Nyström composition:
    K_{b_1 h} o Phi_{a_1 h} o ... o K_{b_6/2 h} o Phi_{a_6 h} o K_{b_6/2 h} o ... o Phi_{a_1 h} o K_{b_1 h},
    strictly bounding retrograde sub-steps (|a_4| ~= 0.2463 << 1.7024).
"""

from __future__ import annotations
from typing import Optional, Set, Sequence, Tuple, Union
import numpy as np
import mpmath as mp

from ..core.constants import (
    G_DEFAULT,
    BM6_DRIFT_WEIGHTS,
    BM6_IMPULSE_WEIGHTS,
    BM6_DRIFT_WEIGHTS_F128,
    BM6_IMPULSE_WEIGHTS_F128,
    to_f128,
)
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

a1_f, a2_f, a3_f, a4_f, a5_f, a6_f = BM6_DRIFT_WEIGHTS_F128
b1_f, b2_f, b3_f, b4_f, b5_f, b6_f = BM6_IMPULSE_WEIGHTS_F128

BM6_DRIFT_SEQUENCE_F128: Tuple[mp.mpf, ...] = (
    a1_f, a2_f, a3_f, a4_f, a5_f, a6_f, a5_f, a4_f, a3_f, a2_f, a1_f
)
BM6_IMPULSE_SEQUENCE_F128: Tuple[mp.mpf, ...] = (
    b1_f, b2_f, b3_f, b4_f, b5_f,
    mp.mpf('0.5') * b6_f, mp.mpf('0.5') * b6_f,
    b5_f, b4_f, b3_f, b2_f, b1_f,
)


class BlanesMoan6Integrator(SymplecticIntegrator):
    """
    Optimal 11-stage sixth-order palindromic Blanes--Moan scheme BM_6^{(11)}.

    Mathematical Properties & Retrograde Sub-Steps (DISC-NUM-02):
    To surmount the Suzuki (1991) barrier for order p >= 3, any symmetric partitioned
    Runge-Kutta-Nyström composition must incorporate negative intermediate substeps.
    In BM_6^{(11)}, intermediate weights a_4 ~= -0.2463, b_3 ~= -0.0400, and b_5 ~= -0.0115 are retrograde.
    Microscopic time-reversibility across these negative steps is preserved bit-exactly
    via the momentum involution T : (r, w) |-> (r, -w) embedded in the universal
    Kepler-Sundman solver (kepler_sundman.py), ensuring S_6(-h) o S_6(h) = Id to
    machine precision (< 1e-15).
    """

    def __init__(
        self,
        G: float = G_DEFAULT,
        screening_mode: str = "c6",
    ) -> None:
        self.G = G
        self.screening_mode = screening_mode
        self._default_apex = ScreenedApexPotential(G=G, screening_mode=screening_mode, order=6)

        # Precompute the symmetric sequence of 11 drift weights and 12 impulse weights
        # a = (a1, a2, a3, a4, a5, a6, a5, a4, a3, a2, a1)
        a1, a2, a3, a4, a5, a6 = BM6_DRIFT_WEIGHTS
        b1, b2, b3, b4, b5, b6 = BM6_IMPULSE_WEIGHTS

        self.drift_sequence = (a1, a2, a3, a4, a5, a6, a5, a4, a3, a2, a1)
        self.impulse_sequence = (b1, b2, b3, b4, b5, 0.5 * b6, 0.5 * b6, b5, b4, b3, b2, b1)
        self.drift_sequence_f128 = BM6_DRIFT_SEQUENCE_F128
        self.impulse_sequence_f128 = BM6_IMPULSE_SEQUENCE_F128

    @property
    def name(self) -> str:
        return "BlanesMoan-BM6"

    @property
    def order(self) -> int:
        return 6

    @property
    def effective_order(self) -> int:
        return 6

    @property
    def evaluations_per_step(self) -> int:
        return 11

    @property
    def num_force_evaluations_per_step(self) -> int:
        return self.evaluations_per_step

    @property
    def fsal_evaluations_per_step(self) -> int:
        return 11

    @property
    def isolated_evaluations_per_step(self) -> int:
        return 12

    @property
    def is_palindromic(self) -> bool:
        return True

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
        Executes a 6th-order BM_6^{(11)} macro-step of size dt via 11 drift and 12 impulse stages.
        """
        if dt == 0.0:
            return state

        eval_apex = apex or self._default_apex

        if state.is_f128:
            dt_val = to_f128(dt)
            drifts = self.drift_sequence_f128
            impulses = self.impulse_sequence_f128
        else:
            dt_val = float(dt)
            drifts = self.drift_sequence
            impulses = self.impulse_sequence

        token = set_current_integrator_stages((drifts, impulses))
        try:
            curr_state = state
            # Initial kick: K_{b_1 dt}
            curr_state = step_apex_kick(
                curr_state,
                tree,
                eval_apex,
                impulses[0] * dt_val,
                supernode_bodies,
                supernode_clusters=supernode_clusters,
            )

            # Stages 1 to 5: Phi_{a_j dt} o K_{b_{j+1} dt}
            for j in range(5):
                curr_state = step_drift_and_kepler(
                    curr_state,
                    tree,
                    drifts[j] * dt_val,
                    G=self.G,
                    supernode_bodies=supernode_bodies,
                    supernode_clusters=supernode_clusters,
                )
                curr_state = step_apex_kick(
                    curr_state,
                    tree,
                    eval_apex,
                    impulses[j + 1] * dt_val,
                    supernode_bodies,
                    supernode_clusters=supernode_clusters,
                )

            # Central drift: Phi_{a_6 dt}
            curr_state = step_drift_and_kepler(
                curr_state,
                tree,
                drifts[5] * dt_val,
                G=self.G,
                supernode_bodies=supernode_bodies,
                supernode_clusters=supernode_clusters,
            )

            # Stages 7 to 11: K_{b_{j+1} dt} o Phi_{a_j dt}
            for j in range(6, 11):
                curr_state = step_apex_kick(
                    curr_state,
                    tree,
                    eval_apex,
                    impulses[j] * dt_val,
                    supernode_bodies,
                    supernode_clusters=supernode_clusters,
                )
                curr_state = step_drift_and_kepler(
                    curr_state,
                    tree,
                    drifts[j] * dt_val,
                    G=self.G,
                    supernode_bodies=supernode_bodies,
                    supernode_clusters=supernode_clusters,
                )

            # Final kick: K_{b_1 dt}
            curr_state = step_apex_kick(
                curr_state,
                tree,
                eval_apex,
                impulses[11] * dt_val,
                supernode_bodies,
                supernode_clusters=supernode_clusters,
            )

            return curr_state
        finally:
            _CURRENT_INTEGRATOR_STAGES.reset(token)
