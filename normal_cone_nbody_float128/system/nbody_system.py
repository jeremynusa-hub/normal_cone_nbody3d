r"""
Unified N-Body Simulation Facade and System Orchestration Engine.

Mathematical Monograph Alignment:
- Chapter 4: Section 4.5 (\label{sec:system_orchestration}, Discrete Skeel & Gear Synchronization).
- Chapter 5: Section 5.5 & Section 5.6 (Dyadic Boundary Telescoping \Psi = C o S^M o C^{-1}, Eq. 5.170).
- Chapter 6: Section 6.1 & Section 6.2 (Theorem 6.1 Anti-Chattering State Machine under gamma = 0.6).
- Orchestrates:
  - Jacobi coordinate charting and canonical cotangent transformations (Layer 1).
  - Analytical Kepler--Sundman universal series and synchronous supernodes (Layer 2).
  - Regularized C^{1,1} apex screening and Hill tidal detection (Layer 3).
  - High-order modular symplectic integrators S_2, S_4, BM_4, BM_6, PBM (Layer 4).
  - Dynamical atlas maintenance under gamma = 0.6 hysteresis and dyadic lattices (Layer 5).
"""

from __future__ import annotations
import warnings
from typing import Optional, Set, Tuple, List, Union
import numpy as np
import mpmath as mp

from ..core.constants import G_DEFAULT, GAMMA_HYSTERESIS, DEFAULT_MAX_DYADIC_LEVEL, to_f128
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree, build_canonical_tree
from ..core.invariants import LieBargmannInvariants, compute_invariants
from ..potential.supernode_detector import detect_supernodes
from ..potential.screened_apex import ScreenedApexPotential
from ..integrators.protocol import SymplecticIntegrator
from ..integrators.strang import StrangIntegrator
from .atlas_manager import AtlasManager
from .dyadic_scheduler import DyadicScheduler
from .history import HistoryTracker
from ..core.profiler import ComputationalProfiler, set_current_profiler, reset_current_profiler


class NBodySystem:
    """
    Unified high-level facade for unsoftened normal cone N-body gravitational simulations.

    Encapsulates:
    - Physical and extended phase state (ExtendedPhaseState).
    - Adapted hierarchical Jacobi tree topology (JacobiTree).
    - Modular symplectic numerical integrator (SymplecticIntegrator).
    - Anti-chattering atlas dynamics manager (AtlasManager).
    - Discrete dyadic multi-rate scheduler (DyadicScheduler).
    - Diagnostic history ledger (HistoryTracker).
    - Fine-grained computational cost profiler (ComputationalProfiler).
    """

    def __init__(
        self,
        state: ExtendedPhaseState,
        tree: Optional[JacobiTree] = None,
        integrator: Optional[SymplecticIntegrator] = None,
        atlas_manager: Optional[AtlasManager] = None,
        scheduler: Optional[DyadicScheduler] = None,
        history: Optional[HistoryTracker] = None,
        G: float = G_DEFAULT,
        auto_supernodes: Optional[bool] = None,
        _ablation_bypass_supernodes: bool = False,
        profiler: Optional[ComputationalProfiler] = None,
        apex_evaluator: Optional[ScreenedApexPotential] = None,
        kepler_evaluator: Optional[Any] = None,
        supernode_evaluator: Optional[Any] = None,
        desingularization_mode: str = "auto",
        enforce_leaf_dichotomy: bool = True,
        enable_dyadic_substepping: bool = True,
        max_dyadic_level: Optional[int] = None,
        in_shadow_mode: bool = True,
        use_optimized_kernels: bool = False,
        dt_macro: float = 0.01,
        supernode_mode: str = "canonical",
        record_substeps: bool = False,
        substep_record_min_level: int = 1,
        atlas_gamma: float = GAMMA_HYSTERESIS,
    ) -> None:
        self.G = float(G)
        self.state = state
        self.in_shadow_mode = bool(in_shadow_mode)
        self.use_optimized_kernels = bool(use_optimized_kernels)
        self.kepler_evaluator = kepler_evaluator
        self.supernode_evaluator = supernode_evaluator
        self.supernode_mode = str(supernode_mode)
        self.record_substeps = bool(record_substeps)
        self.substep_record_min_level = int(substep_record_min_level)

        valid_modes = {"auto", "pure_cartesian", "supernode_only", "leaf_only"}
        if desingularization_mode not in valid_modes:
            raise ValueError(f"Unknown desingularization_mode '{desingularization_mode}'. Expected one of {valid_modes}.")
        self.desingularization_mode = desingularization_mode
        self.enforce_leaf_dichotomy = bool(enforce_leaf_dichotomy)
        self.enable_dyadic_substepping = bool(enable_dyadic_substepping)
        self.max_dyadic_level = max_dyadic_level

        self.atlas_manager = atlas_manager or AtlasManager(gamma=float(atlas_gamma))
        self.tree = tree or self.atlas_manager.build_optimal_tree(state)
        self.integrator = integrator or StrangIntegrator(G=self.G)

        sched_max = max_dyadic_level if max_dyadic_level is not None else DEFAULT_MAX_DYADIC_LEVEL
        self.scheduler = scheduler or DyadicScheduler(dt_macro=float(dt_macro), max_level=sched_max, G=self.G)
        self.history = history or HistoryTracker()
        self.profiler = profiler
        if self.profiler is not None and hasattr(self.history, "profiler"):
            self.history.profiler = self.profiler

        # Architectural Invariant: Synchronous supernodes are intrinsic and mandatory
        # for unsoftened N-body normal-cone desingularization (Chapter 4).
        # Disabling supernodes is strictly confined to comparative research ablation.
        if auto_supernodes is False:
            warnings.warn(
                "Parameter 'auto_supernodes' is deprecated and should not be disabled. "
                "Disabling supernode regularization violates unsoftened normal-cone geometry. "
                "For research ablation studies, use '_ablation_bypass_supernodes=True' or desingularization_mode='pure_cartesian'.",
                category=UserWarning,
                stacklevel=2,
            )
            _ablation_bypass_supernodes = True

        if _ablation_bypass_supernodes or self.desingularization_mode in ("pure_cartesian", "leaf_only"):
            if _ablation_bypass_supernodes and self.desingularization_mode == "auto":
                self.desingularization_mode = "leaf_only"
            self.auto_supernodes = False
        else:
            self.auto_supernodes = True

        self.active_supernodes: List[Tuple[int, ...]] = []
        self._last_dt_substep: Union[float, mp.mpf] = 0.0
        self._last_supernode_bodies: Optional[Set[int]] = None
        self._last_supernode_clusters: Optional[List[Tuple[int, ...]]] = None
        self._on_shadow_manifold: bool = False
        self._shadow_state: Optional[ExtendedPhaseState] = None
        self._last_physical_state: Optional[ExtendedPhaseState] = self.state

        # Initial invariant computation and history baseline
        scr_mode = getattr(self.integrator, "screening_mode", None)
        if self.desingularization_mode == "pure_cartesian":
            scr_mode = "pure_cartesian"
        int_order = getattr(self.integrator, "effective_order", getattr(self.integrator, "order", None))
        if self.use_optimized_kernels and not getattr(state, "is_f128", False):
            from ..potential.optimized_screened_apex import OptimizedScreenedApexPotential
            from ..solvers.optimized_kepler_sundman import step_kepler_sundman_fast
            from ..solvers.jit_supernode import step_supernode_synchronous_jit
            self.apex = apex_evaluator or OptimizedScreenedApexPotential(G=self.G, screening_mode=scr_mode, order=int_order)
            if self.kepler_evaluator is None:
                self.kepler_evaluator = step_kepler_sundman_fast
            if self.supernode_evaluator is None:
                self.supernode_evaluator = step_supernode_synchronous_jit
        else:
            self.apex = apex_evaluator or ScreenedApexPotential(G=self.G, screening_mode=scr_mode, order=int_order)
        self._custom_apex_evaluator = apex_evaluator is not None or (
            self.use_optimized_kernels and not getattr(state, "is_f128", False)
        )
        inv0 = compute_invariants(self.state, G=self.G)
        self.history.record(
            self.state,
            invariants=inv0,
            dyadic_level=0,
            force=True,
            leaf_pairs=tuple(self.tree.get_leaf_pair_indices()),
        )

    @property
    def N(self) -> int:
        """Total particle count in the system."""
        return self.state.N

    @property
    def t(self) -> float:
        """Current physical simulation time."""
        return self.state.t

    def get_state(self) -> ExtendedPhaseState:
        """Returns the current physical extended phase state."""
        return self.state

    def get_physical_state(self) -> ExtendedPhaseState:
        """Non-destructive readout of the current physical state z_n = C_h(tilde_z_n)."""
        return self.state

    def observe_physical_state(self) -> ExtendedPhaseState:
        """Alias for get_physical_state() matching ProcessedSymplecticIntegrator.observe_physical."""
        return self.state

    def get_shadow_state(self) -> ExtendedPhaseState:
        """Returns the internal shadow-manifold state tilde_z_n if active, else physical state."""
        if self._on_shadow_manifold and self._shadow_state is not None:
            return self._shadow_state
        return self.state

    def enter_shadow_manifold(self, dt_substep: Optional[Union[float, mp.mpf]] = None) -> ExtendedPhaseState:
        """Explicitly maps the current physical state onto the shadow manifold."""
        from ..integrators.processed import ProcessedSymplecticIntegrator
        if not isinstance(self.integrator, ProcessedSymplecticIntegrator):
            return self.state
        if self._on_shadow_manifold and self._shadow_state is not None and self.state is self._last_physical_state:
            return self._shadow_state
        h_raw = self._last_dt_substep if dt_substep is None else dt_substep
        h_sub = abs(to_f128(h_raw) if self.state.is_f128 else float(h_raw))
        if h_sub <= 0.0:
            _, h_sub_calc, _ = self.scheduler.get_substep_params(
                t_current=self.state.t, dt_block=self.scheduler.dt_macro, state=self.state, tree=self.tree
            )
            h_sub = abs(to_f128(h_sub_calc) if self.state.is_f128 else float(h_sub_calc))
        sn_bodies = self._detect_supernodes(tree=self.tree)
        sn_clusters = list(self.active_supernodes)
        corr = self.integrator.corrector
        self._shadow_state = type(corr).apply(
            corr,
            self.state,
            self.tree,
            h_sub,
            inverse=True,
            apex=self.apex,
            supernode_bodies=sn_bodies,
            supernode_clusters=sn_clusters,
        )
        self._last_dt_substep = h_sub
        self._last_supernode_bodies = sn_bodies
        self._last_supernode_clusters = sn_clusters
        self._on_shadow_manifold = True
        self._last_physical_state = self.state
        return self._shadow_state

    def exit_shadow_manifold(self) -> ExtendedPhaseState:
        """Exits the shadow manifold and returns the physical state."""
        self._on_shadow_manifold = False
        self._shadow_state = None
        return self.state

    def get_tree(self) -> JacobiTree:
        """Returns the currently active Jacobi tree topology."""
        return self.tree

    def get_invariants(self) -> LieBargmannInvariants:
        """Computes and returns the 10 Lie--Bargmann generators at current physical state."""
        return compute_invariants(self.state, G=self.G)

    def set_integrator(self, integrator: SymplecticIntegrator) -> None:
        """Hot-swaps the symplectic integration engine without state loss."""
        self.integrator = integrator
        self._on_shadow_manifold = False
        self._shadow_state = None
        scr_mode = getattr(self.integrator, "screening_mode", None)
        int_order = getattr(self.integrator, "effective_order", getattr(self.integrator, "order", None))
        if self.use_optimized_kernels and not getattr(self.state, "is_f128", False):
            from ..potential.optimized_screened_apex import OptimizedScreenedApexPotential
            self.apex = OptimizedScreenedApexPotential(G=self.G, screening_mode=scr_mode, order=int_order)
        else:
            self.apex = ScreenedApexPotential(G=self.G, screening_mode=scr_mode, order=int_order)

    def _detect_supernodes(
        self,
        state: Optional[ExtendedPhaseState] = None,
        tree: Optional[JacobiTree] = None,
    ) -> Optional[Set[int]]:
        """Evaluates synchronous supernodes and returns set of member body indices."""
        if not self.auto_supernodes or self.desingularization_mode in ("pure_cartesian", "leaf_only"):
            self.active_supernodes = []
            return None
        st = self.state if state is None else state
        active_tree = self.tree if tree is None else tree
        enforce_non_hier = (self.desingularization_mode == "auto" and self.enforce_leaf_dichotomy)
        clusters = detect_supernodes(
            st,
            active_tree,
            active_supernodes=self.active_supernodes,
            G=self.G,
            enforce_non_hierarchical=enforce_non_hier,
            mode=self.supernode_mode,
        )
        self.active_supernodes = clusters
        if not clusters:
            return None
        sn_set: Set[int] = set()
        for c in clusters:
            sn_set.update(c)
        return sn_set

    def step(
        self,
        dt: Optional[Union[float, mp.mpf]] = None,
        in_shadow_mode: Optional[bool] = None,
        _from_evolve: bool = False,
        dt_macro: Optional[Union[float, mp.mpf]] = None,
    ) -> ExtendedPhaseState:
        """
        Executes a single macro-step of duration dt (defaulting to scheduler.dt_macro)
        following the 4-stage causal macro-boundary pipeline (Proposition 6.14).
        """
        if dt is None and dt_macro is not None:
            dt = dt_macro
        if self.state.is_f128 or isinstance(dt, mp.mpf):
            macro_dt = to_f128(self.scheduler.dt_macro if dt is None else dt)
        else:
            macro_dt = self.scheduler.dt_macro if dt is None else float(dt)
        if macro_dt == 0.0:
            return self.state
        t_start = self.state.t
        use_shadow = self.in_shadow_mode if in_shadow_mode is None else bool(in_shadow_mode)

        if self.state is not self._last_physical_state:
            self._on_shadow_manifold = False
            self._shadow_state = None

        token = None
        if self.profiler is not None and self.profiler.enabled:
            token = set_current_profiler(self.profiler)
            self.profiler.start_macro_step(step_index=len(self.history), t=self.state.t, dt_macro=macro_dt)

        apex_token = None
        if self._custom_apex_evaluator and self.apex is not None:
            from ..integrators.common import set_current_apex_evaluator
            apex_token = set_current_apex_evaluator(self.apex)

        kepler_token = None
        bypass_token = None
        if self.desingularization_mode == "pure_cartesian":
            from ..integrators.common import set_current_bypass_kepler
            bypass_token = set_current_bypass_kepler(True)
        elif self.kepler_evaluator is not None:
            from ..integrators.common import set_current_kepler_solver
            kepler_token = set_current_kepler_solver(self.kepler_evaluator)

        supernode_token = None
        if self.supernode_evaluator is not None:
            from ..integrators.common import set_current_supernode_solver
            supernode_token = set_current_supernode_solver(self.supernode_evaluator)

        from ..integrators.common import set_current_supernode_n_substeps
        sn_n_sub_token = set_current_supernode_n_substeps(1 if not self.enable_dyadic_substepping else None)

        level = 0
        try:
            # 1. Stage 1: Evaluate Jacobi tree adaptation under gamma = 0.60 / gamma_build = 0.70
            old_tree = self.tree
            if float(macro_dt) > 0.0:
                if self.profiler is not None and self.profiler.enabled:
                    with self.profiler.time_stage("atlas"):
                        new_tree, was_flipped = self.atlas_manager.evaluate_tree_adaptation(self.state, self.tree)
                        self.profiler.record_atlas(tree_flipped=was_flipped)
                else:
                    new_tree, was_flipped = self.atlas_manager.evaluate_tree_adaptation(self.state, self.tree)
                active_tree = new_tree if was_flipped else self.tree
            else:
                was_flipped = False
                active_tree = self.tree

            # 2. Stage 2: Synchronous supernode detection with Schmitt trigger hysteresis
            prev_sn_bodies = self._last_supernode_bodies
            prev_sn_clusters = self._last_supernode_clusters
            if self.profiler is not None and self.profiler.enabled:
                with self.profiler.time_stage("supernode"):
                    supernode_bodies = self._detect_supernodes(tree=active_tree)
            else:
                supernode_bodies = self._detect_supernodes(tree=active_tree)
            supernode_clusters = list(self.active_supernodes) if self.active_supernodes else None

            # 3. Stage 3: Dyadic sub-step subdivision
            prev_dt_substep = self._last_dt_substep
            if not self.enable_dyadic_substepping:
                num_substeps = 1
                dt_substep = macro_dt
                level = 0
            elif self.profiler is not None and self.profiler.enabled:
                with self.profiler.time_stage("dyadic"):
                    num_substeps, dt_substep, level = self.scheduler.get_substep_params(
                        t_current=self.state.t,
                        dt_block=macro_dt,
                        state=self.state,
                        tree=active_tree,
                    )
                    self.profiler.record_dyadic(level=level, scale_evals=self.N * (self.N - 1) // 2)
            else:
                num_substeps, dt_substep, level = self.scheduler.get_substep_params(
                    t_current=self.state.t,
                    dt_block=macro_dt,
                    state=self.state,
                    tree=active_tree,
                )
            if supernode_bodies is not None and len(supernode_bodies) >= self.N:
                num_substeps = 1
                dt_substep = macro_dt

            # 4. Stage 4: Gauge Connection & Symplectic Sub-stepping
            from ..integrators.processed import ProcessedSymplecticIntegrator, apply_gauge_connection
            is_processed = isinstance(self.integrator, ProcessedSymplecticIntegrator)
            eps_step = mp.mpf('1e-32') if self.state.is_f128 else 1e-15

            def _execute_macro_block(
                n_sub: int,
                h_sub: Union[float, mp.mpf],
                cur_level: int,
            ) -> Tuple[ExtendedPhaseState, Optional[ExtendedPhaseState], List[ExtendedPhaseState]]:
                h_abs = abs(to_f128(h_sub) if self.state.is_f128 else float(h_sub))
                substep_states: List[ExtendedPhaseState] = []
                log_subs = self.record_substeps and cur_level >= self.substep_record_min_level and n_sub > 1
                if is_processed:
                    corr = self.integrator.corrector
                    if use_shadow:
                        if not self._on_shadow_manifold or self._shadow_state is None:
                            sh_curr = type(corr).apply(
                                corr,
                                self.state,
                                active_tree,
                                h_abs,
                                inverse=True,
                                apex=self.apex,
                                supernode_bodies=supernode_bodies,
                                supernode_clusters=supernode_clusters,
                            )
                        else:
                            sh_curr = self._shadow_state
                            step_scale_changed = (
                                prev_dt_substep is not None
                                and prev_dt_substep > 0.0
                                and abs(h_abs - abs(prev_dt_substep)) > eps_step
                            )
                            supernode_changed = (prev_sn_bodies != supernode_bodies) or (prev_sn_clusters != supernode_clusters)
                            if was_flipped or step_scale_changed or supernode_changed:
                                if self.profiler is not None and self.profiler.enabled:
                                    with self.profiler.time_stage("correctors"):
                                        sh_curr = apply_gauge_connection(
                                            sh_curr,
                                            old_tree,
                                            active_tree,
                                            dt=abs(prev_dt_substep),
                                            dt_new=h_abs,
                                            corrector=corr,
                                            apex=self.apex,
                                            supernode_bodies=prev_sn_bodies,
                                            supernode_clusters=prev_sn_clusters,
                                            supernode_bodies_new=supernode_bodies,
                                            supernode_clusters_new=supernode_clusters,
                                        )
                                        self.profiler.record_corrector(applications=2, directional_evals=4)
                                else:
                                    sh_curr = apply_gauge_connection(
                                        sh_curr,
                                        old_tree,
                                        active_tree,
                                        dt=abs(prev_dt_substep),
                                        dt_new=h_abs,
                                        corrector=corr,
                                        apex=self.apex,
                                        supernode_bodies=prev_sn_bodies,
                                        supernode_clusters=prev_sn_clusters,
                                        supernode_bodies_new=supernode_bodies,
                                        supernode_clusters_new=supernode_clusters,
                                    )
                        for m_sub in range(1, n_sub + 1):
                            sh_curr = self.integrator.base_integrator.step(
                                sh_curr,
                                active_tree,
                                dt=h_sub,
                                apex=self.apex,
                                supernode_bodies=supernode_bodies,
                            )
                            if log_subs and m_sub < n_sub:
                                t_sub = t_start + (to_f128(m_sub) * dt_substep if self.state.is_f128 else float(m_sub) * float(dt_substep))
                                phys_sub = type(corr).apply(
                                    corr,
                                    sh_curr,
                                    active_tree,
                                    h_abs,
                                    inverse=False,
                                    apex=self.apex,
                                    supernode_bodies=supernode_bodies,
                                    supernode_clusters=supernode_clusters,
                                ).with_updates(t=t_sub)
                                substep_states.append(phys_sub)
                        phys_curr = type(corr).apply(
                            corr,
                            sh_curr,
                            active_tree,
                            h_abs,
                            inverse=False,
                            apex=self.apex,
                            supernode_bodies=supernode_bodies,
                            supernode_clusters=supernode_clusters,
                        )
                        return phys_curr, sh_curr, substep_states
                    else:
                        if self.profiler is not None and self.profiler.enabled:
                            with self.profiler.time_stage("correctors"):
                                st_curr = corr.apply_inverse(
                                    self.state,
                                    active_tree,
                                    h_abs,
                                    apex=self.apex,
                                    supernode_bodies=supernode_bodies,
                                )
                                self.profiler.record_corrector(applications=1, directional_evals=2)
                        else:
                            st_curr = corr.apply_inverse(
                                self.state,
                                active_tree,
                                h_abs,
                                apex=self.apex,
                                supernode_bodies=supernode_bodies,
                            )
                        for m_sub in range(1, n_sub + 1):
                            st_curr = self.integrator.base_integrator.step(
                                st_curr,
                                active_tree,
                                dt=h_sub,
                                apex=self.apex,
                                supernode_bodies=supernode_bodies,
                            )
                            if log_subs and m_sub < n_sub:
                                t_sub = t_start + (to_f128(m_sub) * dt_substep if self.state.is_f128 else float(m_sub) * float(dt_substep))
                                phys_sub = corr.apply(
                                    st_curr,
                                    active_tree,
                                    h_abs,
                                    inverse=False,
                                    apex=self.apex,
                                    supernode_bodies=supernode_bodies,
                                ).with_updates(t=t_sub)
                                substep_states.append(phys_sub)
                        if self.profiler is not None and self.profiler.enabled:
                            with self.profiler.time_stage("correctors"):
                                phys_curr = corr.apply(
                                    st_curr,
                                    active_tree,
                                    h_abs,
                                    inverse=False,
                                    apex=self.apex,
                                    supernode_bodies=supernode_bodies,
                                )
                                self.profiler.record_corrector(applications=1, directional_evals=2)
                        else:
                            phys_curr = corr.apply(
                                st_curr,
                                active_tree,
                                h_abs,
                                inverse=False,
                                apex=self.apex,
                                supernode_bodies=supernode_bodies,
                            )
                        return phys_curr, None, substep_states
                else:
                    st_curr = self.state
                    for m_sub in range(1, n_sub + 1):
                        st_curr = self.integrator.step(
                            st_curr,
                            active_tree,
                            dt=h_sub,
                            apex=self.apex,
                            supernode_bodies=supernode_bodies,
                        )
                        if log_subs and m_sub < n_sub:
                            t_sub = t_start + (to_f128(m_sub) * dt_substep if self.state.is_f128 else float(m_sub) * float(dt_substep))
                            substep_states.append(st_curr.with_updates(t=t_sub))
                    return st_curr, None, substep_states

            phys_end, sh_end, sub_states = _execute_macro_block(num_substeps, dt_substep, level)

            # Adjoint time-reversed hysteresis operator A_- under dt < 0:
            # (B_{+dt} o A_+)^{-1} = A_- o B_{-dt}
            if float(macro_dt) < 0.0 and self.N >= 3:
                rev_tree, rev_flipped = self.atlas_manager.evaluate_reverse_tree_adaptation(
                    self.state, phys_end, active_tree
                )
                if rev_flipped:
                    active_tree = rev_tree
                    was_flipped = True
                    self.active_supernodes = list(prev_sn_clusters) if prev_sn_clusters else []
                    supernode_bodies = self._detect_supernodes(tree=active_tree)
                    supernode_clusters = list(self.active_supernodes) if self.active_supernodes else None
                    if self.enable_dyadic_substepping:
                        num_substeps, dt_substep, level = self.scheduler.get_substep_params(
                            t_current=self.state.t,
                            dt_block=macro_dt,
                            state=self.state,
                            tree=active_tree,
                        )
                    if supernode_bodies is not None and len(supernode_bodies) >= self.N:
                        num_substeps = 1
                        dt_substep = macro_dt
                    phys_end, sh_end, sub_states = _execute_macro_block(num_substeps, dt_substep, level)

            # Symmetric dyadic endpoint verification: k_sym = max(k_start, k_end)
            if self.enable_dyadic_substepping and not (supernode_bodies is not None and len(supernode_bodies) >= self.N):
                k_end = self.scheduler.get_level_for_block(
                    t_current=float(t_start + macro_dt),
                    dt_block=-float(macro_dt),
                    state=phys_end,
                    tree=active_tree,
                )
                if k_end > level:
                    level = k_end
                    num_substeps = 1 << level
                    if self.state.is_f128 or isinstance(macro_dt, mp.mpf):
                        dt_substep = to_f128(macro_dt) / num_substeps
                    else:
                        dt_substep = float(macro_dt) / float(num_substeps)
                    phys_end, sh_end, sub_states = _execute_macro_block(num_substeps, dt_substep, level)

            self.tree = active_tree
            self._last_dt_substep = abs(to_f128(dt_substep) if self.state.is_f128 else float(dt_substep))
            self._last_supernode_bodies = supernode_bodies
            self._last_supernode_clusters = supernode_clusters

            t_next = t_start + macro_dt
            self.state = phys_end.with_updates(t=t_next)
            self._last_physical_state = self.state
            if is_processed and use_shadow and sh_end is not None:
                self._shadow_state = sh_end.with_updates(t=t_next)
                self._on_shadow_manifold = True
            else:
                self._shadow_state = None
                self._on_shadow_manifold = False

            # 5. History logging
            active_leaf_pairs = tuple(self.tree.get_leaf_pair_indices())
            for sub_st in sub_states:
                inv_sub = compute_invariants(sub_st, G=self.G)
                self.history.record(sub_st, invariants=inv_sub, dyadic_level=level, leaf_pairs=active_leaf_pairs)
            if self.profiler is not None and self.profiler.enabled:
                with self.profiler.time_stage("diagnostics"):
                    inv = compute_invariants(self.state, G=self.G)
                    self.history.record(self.state, invariants=inv, dyadic_level=level, leaf_pairs=active_leaf_pairs)
            else:
                inv = compute_invariants(self.state, G=self.G)
                self.history.record(self.state, invariants=inv, dyadic_level=level, leaf_pairs=active_leaf_pairs)

            return self.state
        finally:
            if sn_n_sub_token is not None:
                from ..integrators.common import reset_current_supernode_n_substeps
                reset_current_supernode_n_substeps(sn_n_sub_token)
            if bypass_token is not None:
                from ..integrators.common import _CURRENT_BYPASS_KEPLER
                _CURRENT_BYPASS_KEPLER.reset(bypass_token)
            if supernode_token is not None:
                from ..integrators.common import _CURRENT_SUPERNODE_SOLVER
                _CURRENT_SUPERNODE_SOLVER.reset(supernode_token)
            if kepler_token is not None:
                from ..integrators.common import _CURRENT_KEPLER_SOLVER
                _CURRENT_KEPLER_SOLVER.reset(kepler_token)
            if apex_token is not None:
                from ..integrators.common import reset_current_apex_evaluator
                reset_current_apex_evaluator(apex_token)
            if token is not None:
                r_min_val = self.state.min_pairwise_distance()
                self.profiler.end_macro_step(r_min=r_min_val, dyadic_level=level)
                reset_current_profiler(token)


    def evolve(self, t_span: Union[float, mp.mpf], dt: Optional[Union[float, mp.mpf]] = None) -> ExtendedPhaseState:
        """
        Evolves the N-body system forward over duration t_span from current time:
            t_final = t_initial + t_span.

        For ProcessedSymplecticIntegrator, executes global shadow telescoping:
            Psi^{(M)} = C(h) o S(h)^M o C^{-1}(h)
        with symplectic gauge connections T_{(A,h_A)->(B,h_B)}^{(PBM4)} on chart/scale flips,
        reducing corrector evaluations from 2M to 2 + 2*N_flips.
        """
        if self.state.is_f128 or isinstance(t_span, mp.mpf) or isinstance(dt, mp.mpf):
            t_span_val = to_f128(t_span)
            macro_dt = to_f128(self.scheduler.dt_macro if dt is None else dt)
            eps_time = mp.mpf('1e-32')
        else:
            t_span_val = float(t_span)
            macro_dt = self.scheduler.dt_macro if dt is None else float(dt)
            eps_time = 1e-14

        if t_span_val <= 0.0:
            return self.state

        t_target = self.state.t + t_span_val

        from ..integrators.processed import ProcessedSymplecticIntegrator
        is_processed = isinstance(self.integrator, ProcessedSymplecticIntegrator)

        if is_processed:
            # 1. Global entrance onto shadow manifold using initial substep dt_sub0
            if not self._on_shadow_manifold or self._shadow_state is None or self.state is not self._last_physical_state:
                new_tree0, was_flipped0 = self.atlas_manager.evaluate_tree_adaptation(self.state, self.tree)
                if was_flipped0:
                    self.tree = new_tree0
                sn_bodies0 = self._detect_supernodes(tree=self.tree)
                self._last_supernode_bodies = sn_bodies0
                self._last_supernode_clusters = list(self.active_supernodes) if self.active_supernodes else None
                if not self.enable_dyadic_substepping:
                    dt_sub0 = macro_dt
                else:
                    _, dt_sub0, _ = self.scheduler.get_substep_params(
                        t_current=self.state.t,
                        dt_block=macro_dt,
                        state=self.state,
                        tree=self.tree,
                    )
                self._last_dt_substep = abs(to_f128(dt_sub0) if self.state.is_f128 else float(dt_sub0))
                self._shadow_state = self.integrator.corrector.apply_inverse(
                    self.state, self.tree, self._last_dt_substep, apex=self.apex, supernode_bodies=sn_bodies0
                )
                self._on_shadow_manifold = True
                self._last_physical_state = self.state

            # 2. Interior macro-steps on shadow manifold (telescoping)
            while (t_target - self.state.t) > eps_time:
                dt_remaining = t_target - self.state.t
                dt_step = min(macro_dt, dt_remaining)
                self.step(dt=dt_step, in_shadow_mode=True, _from_evolve=True)

            # 3. Global exit from shadow manifold back to physical coordinates using exit_dt_sub
            exit_dt_sub = self._last_dt_substep if self._last_dt_substep > 0.0 else abs(macro_dt)
            exit_sn = self._last_supernode_bodies
            if self._shadow_state is not None:
                self.state = self.integrator.corrector.apply(
                    self._shadow_state, self.tree, exit_dt_sub, inverse=False, apex=self.apex, supernode_bodies=exit_sn
                ).with_updates(t=self._shadow_state.t)
                self._last_physical_state = self.state
            # Update the last history record with the true physical state
            if self.history and self.history.snapshots:
                from .history import HistorySnapshot
                inv_final = compute_invariants(self.state, G=self.G)
                last_snap = self.history.snapshots[-1]
                self.history.snapshots[-1] = HistorySnapshot(
                    t=self.state.t,
                    q=self.state.q.copy(),
                    p=self.state.p.copy(),
                    invariants=inv_final,
                    dyadic_level=last_snap.dyadic_level,
                    leaf_pairs=last_snap.leaf_pairs,
                )
        else:
            while (t_target - self.state.t) > eps_time:
                dt_remaining = t_target - self.state.t
                dt_step = min(macro_dt, dt_remaining)
                self.step(dt=dt_step)

        return self.state

    def reset(self, new_state: ExtendedPhaseState) -> None:
        """Resets the system with a new initial phase state and fresh tree."""
        self.state = new_state
        self.tree = self.atlas_manager.build_optimal_tree(new_state)
        self.atlas_manager.reset_history()
        self.history.clear()
        self._last_dt_substep = 0.0
        self._last_supernode_bodies = None
        self._last_supernode_clusters = None
        self._on_shadow_manifold = False
        self._shadow_state = None
        self._last_physical_state = self.state
        inv0 = compute_invariants(self.state, G=self.G)
        self.history.record(
            self.state,
            invariants=inv0,
            dyadic_level=0,
            force=True,
            leaf_pairs=tuple(self.tree.get_leaf_pair_indices()),
        )
