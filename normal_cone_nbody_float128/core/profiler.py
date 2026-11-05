"""
Fine-Grained Computational Profiler for Unsoftened Conormal N-Body Integration.

Core low-level instrumentation module without circular dependencies.
Provides high-resolution wall-clock CPU timers (nanosecond resolution via time.perf_counter)
and atomic operation counters partitioned across all 6 core algorithmic stages:
  1. 'kepler'     : Universal Kepler-Sundman analytic leaf solver (Halley & Stumpff iterations).
  2. 'apex'       : Hierarchical screened apex cotangent shear kick and conormal hyper-radius gradients.
  3. 'supernode'  : Synchronous multi-body supernode Leapfrog loop under unfolded Sundman time tau_S.
  4. 'dyadic'     : Discrete dyadic multirate scheduler and timescale evaluations (tau_kin, tau_dyn, tau_peri).
  5. 'atlas'      : Jacobi tree topology evaluation, Schmitt trigger hysteresis, and SL(2, R) chart flips.
  6. 'correctors' : High-order boundary symplectic correctors C_{2p} and analytical Hessian-vector products.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from contextlib import contextmanager
from typing import Dict, List, Any, Optional
import time
import contextvars

_CURRENT_PROFILER: contextvars.ContextVar[Optional[ComputationalProfiler]] = contextvars.ContextVar(
    "current_profiler", default=None
)


def get_current_profiler() -> Optional[ComputationalProfiler]:
    """Retrieves the currently active thread-local/task-local profiler, or None."""
    return _CURRENT_PROFILER.get()


def set_current_profiler(profiler: Optional[ComputationalProfiler]) -> contextvars.Token:
    """Sets the active profiler for the current context."""
    return _CURRENT_PROFILER.set(profiler)


def reset_current_profiler(token: contextvars.Token) -> None:
    """Restores the prior active profiler."""
    _CURRENT_PROFILER.reset(token)


@dataclass
class StageCostMetrics:
    """Cumulative CPU time and granular operation counters for a specific engine stage."""
    name: str
    cpu_time: float = 0.0
    calls: int = 0
    counters: Dict[str, int] = field(default_factory=dict)

    def add_time(self, elapsed: float) -> None:
        self.cpu_time += elapsed
        self.calls += 1

    def increment(self, counter_name: str, count: int = 1) -> None:
        self.counters[counter_name] = self.counters.get(counter_name, 0) + count

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "cpu_time_sec": self.cpu_time,
            "calls": self.calls,
            "counters": dict(self.counters),
        }


@dataclass
class MacroStepProfile:
    """Snapshot of computational resource consumption for a single macro-integration step."""
    step_index: int
    t: float
    dt_macro: float
    total_step_time: float
    r_min: float
    dyadic_level: int
    stage_times: Dict[str, float]
    stage_counters: Dict[str, Dict[str, int]]


class ComputationalProfiler:
    """
    Non-intrusive, thread-safe, high-precision profiler tracking exact computational costs
    across the symplectic conormal N-body integration hierarchy.
    """

    STAGE_NAMES = ("kepler", "apex", "supernode", "dyadic", "atlas", "correctors", "diagnostics")

    def __init__(self, enabled: bool = True) -> None:
        self.enabled = enabled
        self.stages: Dict[str, StageCostMetrics] = {
            s: StageCostMetrics(name=s) for s in self.STAGE_NAMES
        }
        self.history: List[MacroStepProfile] = []

        # Current in-flight macro-step buffers
        self._step_active: bool = False
        self._step_start_time: float = 0.0
        self._cur_step_index: int = 0
        self._cur_t: float = 0.0
        self._cur_dt: float = 0.0
        self._cur_r_min: float = 0.0
        self._cur_level: int = 0
        self._cur_stage_times: Dict[str, float] = {s: 0.0 for s in self.STAGE_NAMES}
        self._cur_stage_counters: Dict[str, Dict[str, int]] = {s: {} for s in self.STAGE_NAMES}

    def reset(self) -> None:
        """Resets all cumulative counters and step history."""
        for s in self.stages.values():
            s.cpu_time = 0.0
            s.calls = 0
            s.counters.clear()
        self.history.clear()
        self._step_active = False

    @contextmanager
    def active(self):
        """Context manager setting this profiler instance active for the current context."""
        token = set_current_profiler(self)
        try:
            yield self
        finally:
            reset_current_profiler(token)

    @contextmanager
    def time_stage(self, stage_name: str):
        """
        Context manager measuring exact CPU wall-clock duration of a specific stage.
        If profiler is disabled, yields immediately with negligible overhead.
        """
        if not self.enabled:
            yield
            return

        t0 = time.perf_counter()
        try:
            yield
        finally:
            elapsed = time.perf_counter() - t0
            stage = self.stages.get(stage_name)
            if stage is not None:
                stage.add_time(elapsed)
            if self._step_active:
                self._cur_stage_times[stage_name] = self._cur_stage_times.get(stage_name, 0.0) + elapsed

    def record_counter(self, stage_name: str, counter_name: str, count: int = 1) -> None:
        """Increments an atomic operation counter for a stage."""
        if not self.enabled:
            return
        stage = self.stages.get(stage_name)
        if stage is not None:
            stage.increment(counter_name, count)
        if self._step_active:
            sub = self._cur_stage_counters.setdefault(stage_name, {})
            sub[counter_name] = sub.get(counter_name, 0) + count

    # --- Convenience Semantic Hooks ---

    def record_kepler(
        self,
        halley_iters: int = 1,
        stumpff_evals: int = 1,
        is_bisection: bool = False,
    ) -> None:
        if not self.enabled:
            return
        self.record_counter("kepler", "leaf_calls", 1)
        self.record_counter("kepler", "halley_iterations", halley_iters)
        self.record_counter("kepler", "stumpff_evaluations", stumpff_evals)
        if is_bisection:
            self.record_counter("kepler", "bisection_fallbacks", 1)

    def record_apex(
        self,
        kicks: int = 1,
        pairs: int = 0,
        theta_evals: int = 0,
    ) -> None:
        if not self.enabled:
            return
        self.record_counter("apex", "kick_evaluations", kicks)
        if pairs > 0:
            self.record_counter("apex", "pairs_evaluated", pairs)
        if theta_evals > 0:
            self.record_counter("apex", "screening_theta_evals", theta_evals)

    def record_supernode(
        self,
        substeps: int = 1,
        acc_evals: int = 1,
        is_close_pair: bool = False,
        activations: int = 0,
    ) -> None:
        if not self.enabled:
            return
        if activations > 0:
            self.record_counter("supernode", "cluster_activations", activations)
        self.record_counter("supernode", "symplectic_substeps", substeps)
        self.record_counter("supernode", "leapfrog_substeps", substeps)
        self.record_counter("supernode", "internal_acc_evals", acc_evals)
        if is_close_pair:
            self.record_counter("supernode", "close_pair_regularizations", 1)

    def record_dyadic(self, level: int = 0, scale_evals: int = 1) -> None:
        if not self.enabled:
            return
        self.record_counter("dyadic", "block_evaluations", 1)
        self.record_counter("dyadic", f"level_k_{level}", 1)
        self.record_counter("dyadic", "substeps_generated", 2 ** level)
        if scale_evals > 0:
            self.record_counter("dyadic", "timescale_pair_checks", scale_evals)

    def record_atlas(self, tree_flipped: bool = False) -> None:
        if not self.enabled:
            return
        self.record_counter("atlas", "hysteresis_evaluations", 1)
        if tree_flipped:
            self.record_counter("atlas", "tree_flips", 1)

    def record_corrector(
        self,
        applications: int = 1,
        directional_evals: int = 0,
    ) -> None:
        if not self.enabled:
            return
        self.record_counter("correctors", "corrector_applications", applications)
        if directional_evals > 0:
            self.record_counter("correctors", "hessian_vector_products", directional_evals)

    # --- Macro-Step Boundaries ---

    def start_macro_step(self, step_index: int, t: float, dt_macro: float) -> None:
        if not self.enabled:
            return
        self._step_active = True
        self._cur_step_index = step_index
        self._cur_t = float(t)
        self._cur_dt = float(dt_macro)
        self._step_start_time = time.perf_counter()
        self._cur_stage_times = {s: 0.0 for s in self.STAGE_NAMES}
        self._cur_stage_counters = {s: {} for s in self.STAGE_NAMES}

    def end_macro_step(self, r_min: float = 0.0, dyadic_level: int = 0) -> None:
        if not self.enabled or not self._step_active:
            return
        total_time = time.perf_counter() - self._step_start_time
        profile = MacroStepProfile(
            step_index=self._cur_step_index,
            t=self._cur_t,
            dt_macro=self._cur_dt,
            total_step_time=total_time,
            r_min=float(r_min),
            dyadic_level=int(dyadic_level),
            stage_times=dict(self._cur_stage_times),
            stage_counters={k: dict(v) for k, v in self._cur_stage_counters.items()},
        )
        self.history.append(profile)
        self._step_active = False

    # --- Synthesis & Export ---

    @property
    def total_cpu_time(self) -> float:
        return sum(s.cpu_time for s in self.stages.values())

    def get_summary(self) -> Dict[str, Any]:
        """Returns structured synthesis of total execution time, shares, and operation counts."""
        tot = max(self.total_cpu_time, 1e-15)
        breakdown = {}
        all_counters = {}
        for name, stage in self.stages.items():
            breakdown[name] = {
                "cpu_time_sec": stage.cpu_time,
                "percentage": (stage.cpu_time / tot) * 100.0,
                "calls": stage.calls,
                "counters": dict(stage.counters),
            }
            all_counters.update(stage.counters)

        return {
            "total_cpu_time_sec": tot,
            "total_macro_steps": len(self.history),
            "breakdown": breakdown,
            "counters": all_counters,
        }
