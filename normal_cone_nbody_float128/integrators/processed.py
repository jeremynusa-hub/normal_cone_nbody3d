"""
The Processed Symplectic Integrator Architecture PBM with Dyadic Boundary Telescoping.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.5, Proposition 5.3 & Equation (5.170)
  (Boundary Operator Telescoping Identity and Comparative Synthesis Table 1).
- Operator Composition:
    Psi_{2p+2}(h) = C_{2p}(h) o S_{2p}(h) o C_{2p}(h)^{-1}.
- Global Trajectory Telescoping over M Steps:
    Psi^{(M)}(h) = C_{2p}(h) o [ prod_{m=1}^M S_{2p}(h) ] o C_{2p}(h)^{-1}.
    All internal corrector pairs telescope identically to Id (C^{-1} o C = Id),
    achieving order 2p+2 at strictly ZERO additional intermediate force evaluations!
- Formal Lean 4 Alignment: 02_CorrectorBM6_Theorems.lean
  (theorem_9_4_dyadic_telescoping_under_hysteresis).
"""

from __future__ import annotations
from typing import Optional, Set, Sequence, List, Tuple, Union
import numpy as np
import mpmath as mp

from ..core.constants import G_DEFAULT
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree
from ..potential.screened_apex import ScreenedApexPotential
from .protocol import SymplecticIntegrator
from .correctors import (
    SymplecticCorrector,
    CorrectorC4,
    GlobalCorrectorC4,
    CotangentLiftCorrectorC4,
    HomologicalLieCorrectorC4,
    compute_conormal_corrector_kick,
    compute_global_conormal_c4_kick,
    compute_cotangent_lift_c4_flow,
)
from .blanes_moan_4 import BlanesMoan4Integrator
from .blanes_moan_6 import BlanesMoan6Integrator


class ProcessedSymplecticIntegrator(SymplecticIntegrator):
    """
    Decorator implementing processed symplectic integration PBM:
    elevates base convergence order by +2 in the asymptotic hierarchical regime eps << h
    at zero marginal force evaluation cost.
    """

    def __init__(
        self,
        base_integrator: SymplecticIntegrator,
        corrector: Optional[SymplecticCorrector] = None,
    ) -> None:
        self.base_integrator = base_integrator

        if corrector is not None:
            self.corrector = corrector
        else:
            base_g = getattr(base_integrator, "G", G_DEFAULT)
            self.corrector = CorrectorC4(G=base_g, mode="homological")

    @property
    def name(self) -> str:
        return f"Processed-{self.base_integrator.name}"

    @property
    def order(self) -> int:
        return self.base_integrator.order

    @property
    def effective_order(self) -> int:
        """
        Formal effective order label: base_integrator.order + 2 (6 for PBM4 in the hierarchical regime eps << h).
        Note: By Theorems 5.14 & 5.15 of the monograph, PBM4 achieves effective order 6 on D_reg when eps << h,
        with a quadratic tidal floor O(eps^2 h^4).
        """
        return self.base_integrator.order + 2

    @property
    def evaluations_per_step(self) -> int:
        return self.base_integrator.evaluations_per_step

    @property
    def num_force_evaluations_per_step(self) -> int:
        return self.base_integrator.num_force_evaluations_per_step

    @property
    def fsal_evaluations_per_step(self) -> int:
        return self.base_integrator.fsal_evaluations_per_step

    @property
    def isolated_evaluations_per_step(self) -> int:
        return self.base_integrator.isolated_evaluations_per_step

    @property
    def is_palindromic(self) -> bool:
        return True

    @property
    def screening_mode(self) -> str:
        """Returns the screening mode adapted to the effective order of the processed integrator."""
        if hasattr(self.corrector, "screening_mode") and self.corrector.screening_mode:
            return self.corrector.screening_mode
        return "c8" if self.effective_order >= 8 else "c6"

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
        Executes an isolated single processed step of size dt:
            Psi(dt) = C(dt) o S(dt) o C^{-1}(dt).
        """
        if dt == 0.0:
            return state

        # 1. Apply algebraic inverse corrector: C^{-1}(dt) = exp(+ dt^{2p} X_K)
        s_in = self.corrector.apply_inverse(
            state,
            tree,
            dt,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )

        # 2. Propagate base symplectic step: S(dt)
        s_mid = self.base_integrator.step(
            s_in,
            tree,
            dt,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )

        # 3. Apply forward corrector: C(dt) = exp(- dt^{2p} X_K)
        s_out = self.corrector.apply(
            s_mid,
            tree,
            dt,
            inverse=False,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )

        return s_out

    def enter_shadow(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        dt: Union[float, mp.mpf],
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """Applies the inverse corrector C^{-1}(dt) to map physical state into shadow space."""
        if dt == 0.0:
            return state
        return self.corrector.apply_inverse(
            state,
            tree,
            dt,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )

    def step_shadow(
        self,
        shadow_state: ExtendedPhaseState,
        tree: JacobiTree,
        dt: Union[float, mp.mpf],
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """Advances the shadow state by one step of the base symplectic kernel S(dt)."""
        if dt == 0.0:
            return shadow_state
        return self.base_integrator.step(
            shadow_state,
            tree,
            dt,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )

    def observe_physical(
        self,
        shadow_state: ExtendedPhaseState,
        tree: JacobiTree,
        dt: Union[float, mp.mpf],
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """Non-destructively projects shadow state tilde_z_n |-> z_n = C_h(tilde_z_n) for physical readout."""
        if dt == 0.0:
            return shadow_state
        return self.corrector.apply(
            shadow_state,
            tree,
            dt,
            inverse=False,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )

    def exit_shadow(
        self,
        shadow_state: ExtendedPhaseState,
        tree: JacobiTree,
        dt: Union[float, mp.mpf],
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """Applies the forward corrector C(dt) to project shadow state back to physical space."""
        if dt == 0.0:
            return shadow_state
        return self.corrector.apply(
            shadow_state,
            tree,
            dt,
            inverse=False,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )

    def integrate(
        self,
        initial_state: ExtendedPhaseState,
        tree: JacobiTree,
        dt: Union[float, mp.mpf],
        num_steps: int,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
        record_states: bool = False,
    ) -> Tuple[ExtendedPhaseState, List[ExtendedPhaseState]]:
        """
        Integrates a trajectory of num_steps using dyadic boundary telescoping:
            Psi^{(M)}(dt) = C(dt) o [ S(dt)^M ] o C^{-1}(dt).
        Strictly zero corrector force evaluations in the interior steps!
        """
        if num_steps <= 0 or dt == 0.0:
            return initial_state, [initial_state] if record_states else []

        history: List[ExtendedPhaseState] = [initial_state] if record_states else []

        # Boundary enter: C^{-1}(dt) evaluated exactly once at t_0
        curr_state = self.enter_shadow(
            initial_state,
            tree,
            dt,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )

        # Intermediate macro-steps: pure base integrator (telescoping C o C^{-1} = Id)
        for _ in range(num_steps):
            curr_state = self.step_shadow(
                curr_state,
                tree,
                dt,
                apex=apex,
                supernode_bodies=supernode_bodies,
                supernode_clusters=supernode_clusters,
            )
            if record_states:
                # To record physical state on the true shadow manifold, apply C(dt)
                recorded = self.observe_physical(
                    curr_state,
                    tree,
                    dt,
                    apex=apex,
                    supernode_bodies=supernode_bodies,
                    supernode_clusters=supernode_clusters,
                )
                history.append(recorded)

        # Boundary exit: C(dt) evaluated exactly once at t_M
        final_state = self.exit_shadow(
            curr_state,
            tree,
            dt,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )

        if record_states and history:
            history[-1] = final_state

        return final_state, history


_UNSET_SN = object()


def apply_gauge_connection(
    state: ExtendedPhaseState,
    old_tree: JacobiTree,
    new_tree: JacobiTree,
    dt: Union[float, mp.mpf],
    corrector: Optional[SymplecticCorrector] = None,
    apex: Optional[ScreenedApexPotential] = None,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    dt_new: Optional[Union[float, mp.mpf]] = None,
    supernode_bodies_new: Any = _UNSET_SN,
    supernode_clusters_new: Any = _UNSET_SN,
) -> ExtendedPhaseState:
    """
    Applies the covariant symplectic gauge connection operator:
        T_{(A, h_A) -> (B, h_B)}^{(PBM4)} = (C_{h_B}^{(B)})^{-1} o M_{A -> B} o C_{h_A}^{(A)}
    across an atlas chart transition T_A -> T_B, dyadic step-size transition h_A -> h_B,
    and/or synchronous supernode promotion/demotion boundary on the shadow Hamiltonian manifold.

    Mathematical Monograph Alignment:
    - Chapter 5: Section 5.6, Theorems 5.17 & 5.18, Equations (5.165)-(5.195).
    - Chapter 6: Section 6.5, Theorem 6.16 & Algorithm 5 (Step 4.4).
    - Default `corrector=None` uses `CorrectorC4(mode="homological")` (Z_4^{BM4} in P_1).
    - Eliminates the interfacial shadow Hamiltonian jump Delta H = O(h^4) across chart/rung
      boundaries, elevating global trajectory conservation to O(h^6).

    Properties:
    - Identity: If both tree topology and step size are unchanged, returns state identically.
    - Exact Cech 1-cocycle closure: T_{(C,h_C)->(A,h_A)} o T_{(B,h_B)->(C,h_C)} o T_{(A,h_A)->(B,h_B)} == Id.
    - Symplecticity: Preserves the canonical Darboux 2-form and all 9 Lie--Bargmann generators (P, L, G).
    """
    if corrector is None:
        corrector = CorrectorC4(G=G_DEFAULT, mode="homological")

    h_old = dt
    h_new = dt if dt_new is None else dt_new
    sn_new = supernode_bodies if supernode_bodies_new is _UNSET_SN else supernode_bodies_new
    snc_new = supernode_clusters if supernode_clusters_new is _UNSET_SN else supernode_clusters_new

    same_chart = (
        old_tree is new_tree
        or old_tree.get_leaf_pair_indices() == new_tree.get_leaf_pair_indices()
    )
    same_step = abs(float(h_old) - float(h_new)) <= 1e-15 * max(abs(float(h_old)), abs(float(h_new)), 1.0)
    same_sn = (supernode_bodies == sn_new) and (supernode_clusters == snc_new)

    if state.N <= 1 or (float(h_old) == 0.0 and float(h_new) == 0.0) or (same_chart and same_step and same_sn):
        return state

    if isinstance(corrector, GlobalCorrectorC4) or (
        isinstance(corrector, CorrectorC4) and getattr(corrector, "mode", "homological") in ("global", "legacy")
    ):
        # Unipotent configuration-only ablation kick in P_0 (Delta q == 0, Delta p = dp_A - dp_B).
        is_global = isinstance(corrector, GlobalCorrectorC4) or getattr(corrector, "mode", "global") == "global"
        divisor = 480.0 if is_global else 1440.0
        dp_A = (
            compute_conormal_corrector_kick(
                state, old_tree, h_old, G=corrector.G, inverse=False, divisor=divisor,
                apex=apex, supernode_bodies=supernode_bodies, supernode_clusters=supernode_clusters,
            )
            if float(h_old) != 0.0
            else np.zeros_like(state.p)
        )
        dp_B = (
            compute_conormal_corrector_kick(
                state, new_tree, h_new, G=corrector.G, inverse=False, divisor=divisor,
                apex=apex, supernode_bodies=sn_new, supernode_clusters=snc_new,
            )
            if float(h_new) != 0.0
            else np.zeros_like(state.p)
        )
        if state.is_f128:
            new_p = np.array([
                [state.p[i, d] + dp_A[i, d] - dp_B[i, d] for d in range(3)]
                for i in range(state.N)
            ], dtype=object)
        else:
            new_p = state.p + dp_A - dp_B
        return state.with_updates(p=new_p)

    # 1. Forward Shadow Exit from chart A at step h_old: C_{h_A}^{(A)}(tilde_z_A) -> z_phys
    if float(h_old) != 0.0:
        s_phys = corrector.apply(
            state,
            old_tree,
            h_old,
            inverse=False,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
        )
    else:
        s_phys = state

    # 2. Backward Shadow Re-Entry into chart B at step h_new: (C_{h_B}^{(B)})^{-1}(z_phys) -> tilde_z_B
    if float(h_new) != 0.0:
        s_b = corrector.apply_inverse(
            s_phys,
            new_tree,
            h_new,
            apex=apex,
            supernode_bodies=sn_new,
            supernode_clusters=snc_new,
        )
    else:
        s_b = s_phys

    return s_b


class PBM4GaugeTransitionOperator:
    """
    Canonical covariant symplectic gauge transition operator:
        T_{(A, h_A) -> (B, h_B)}^{(PBM4)} = (C_{h_B}^{(B)})^{-1} o M_{A -> B} o C_{h_A}^{(A)}
    governing transitions between Jacobi atlas charts and dyadic step-size levels on the
    shadow manifold (Section 5.6 Theorems 5.17-5.18 & Section 6.5 Algorithm 5 Step 4.4).
    """

    def __init__(
        self,
        corrector: Optional[SymplecticCorrector] = None,
        G: Union[float, mp.mpf] = G_DEFAULT,
        screening_mode: str = "c6",
        mode: str = "homological",
    ) -> None:
        self.G = G
        self.screening_mode = screening_mode
        self.mode = mode
        self.corrector = (
            corrector
            if corrector is not None
            else CorrectorC4(G=self.G, screening_mode=screening_mode, mode=mode)
        )

    def apply(
        self,
        state: ExtendedPhaseState,
        old_tree: JacobiTree,
        new_tree: JacobiTree,
        dt: Union[float, mp.mpf],
        dt_new: Optional[Union[float, mp.mpf]] = None,
        apex: Optional[ScreenedApexPotential] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[Sequence[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """Applies T_{(A, h_A) -> (B, h_B)}^{(PBM4)} on Cartesian extended phase state."""
        return apply_gauge_connection(
            state=state,
            old_tree=old_tree,
            new_tree=new_tree,
            dt=dt,
            corrector=self.corrector,
            apex=apex,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
            dt_new=dt_new,
        )

    def verify_cech_cocycle(
        self,
        state: ExtendedPhaseState,
        tree_A: JacobiTree,
        tree_B: JacobiTree,
        tree_C: JacobiTree,
        dt_A: Union[float, mp.mpf],
        dt_B: Optional[Union[float, mp.mpf]] = None,
        dt_C: Optional[Union[float, mp.mpf]] = None,
        apex: Optional[ScreenedApexPotential] = None,
    ) -> Union[float, mp.mpf]:
        """
        Verifies the 3-chart non-commutative Cech 1-cocycle identity:
            T_{(C, h_C) -> (A, h_A)} o T_{(B, h_B) -> (C, h_C)} o T_{(A, h_A) -> (B, h_B)}(z) == z
        and returns the maximum phase-space infinity-norm residual.
        """
        hA = dt_A
        hB = hA if dt_B is None else dt_B
        hC = hA if dt_C is None else dt_C

        s_B = self.apply(state, tree_A, tree_B, dt=hA, dt_new=hB, apex=apex)
        s_C = self.apply(s_B, tree_B, tree_C, dt=hB, dt_new=hC, apex=apex)
        s_A = self.apply(s_C, tree_C, tree_A, dt=hC, dt_new=hA, apex=apex)

        if state.is_f128:
            dq = max(abs(s_A.q[i, d] - state.q[i, d]) for i in range(state.N) for d in range(3))
            dp = max(abs(s_A.p[i, d] - state.p[i, d]) for i in range(state.N) for d in range(3))
            return max(dq, dp)
        dq = float(np.max(np.abs(s_A.q - state.q)))
        dp = float(np.max(np.abs(s_A.p - state.p)))
        return max(dq, dp)


def create_pbm4_integrator(
    G: Union[float, mp.mpf] = G_DEFAULT,
    screening_mode: str = "c6",
    mode: str = "homological",
) -> ProcessedSymplecticIntegrator:
    """Factory creating PBM4 (Blanes--Moan 4 with C4 corrector in 'homological' (default), 'cotangent_lift', or 'global' mode)."""
    return ProcessedSymplecticIntegrator(
        base_integrator=BlanesMoan4Integrator(G=G),
        corrector=CorrectorC4(G=G, screening_mode=screening_mode, mode=mode),
    )



