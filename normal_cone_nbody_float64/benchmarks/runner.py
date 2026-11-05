r"""
Modular Benchmark Runner and Work-Precision Pareto Scaling Analyzer.

Theoretical Foundations:
- Monograph Reference: Chapter 7, Section 7.1 (\label{sec:unsoftened_benchmarks})
  and Section 7.2 (\label{sec:long_term_stability}); Chapter 5, Sections 5.1-5.5.
- Symplectic Work-Precision Pareto Relation:
    Work W = N_steps * S_{stages}.
    Energy Error epsilon = max_t |H(t) - H_0| / |H_0| ~ C * h^p ~ C * (W / T)^{-p}.
    On a log-log Pareto plot (W vs epsilon), an order-p integrator exhibits a slope of -p.
    Higher-order conormal schemes (BM4, BM6) achieve machine-precision floor (< 1e-14)
    with orders of magnitude less computational work than Strang S_2 or Yoshida S_4.
"""

from __future__ import annotations
import time
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Union, Tuple
import numpy as np

from ..core.constants import GRAVITATIONAL_CONSTANT
from ..core.state import ExtendedPhaseState
from ..system.nbody_system import NBodySystem
from ..integrators.protocol import SymplecticIntegrator
from ..diagnostics.metrics import compute_celestial_kpis, KPIReport
from .registry import get_benchmark, BENCHMARK_REGISTRY, get_benchmark_profile


@dataclass
class BenchmarkResult:
    """
    Structured outcome of an executed gravitational benchmark simulation.
    """
    benchmark_name: str
    integrator_name: str
    integrator_order: int
    stages_per_step: int
    dt: float
    t_span: float
    n_steps: int
    force_evaluations: int
    execution_wall_time: float
    kpi_report: KPIReport
    final_state: ExtendedPhaseState
    force_evaluations_fsal: int = 0
    fsal_stages_per_step: int = 1
    tau_per_stage: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "benchmark_name": self.benchmark_name,
            "integrator_name": self.integrator_name,
            "integrator_order": int(self.integrator_order),
            "stages_per_step": int(self.stages_per_step),
            "fsal_stages_per_step": int(self.fsal_stages_per_step),
            "dt": float(self.dt),
            "tau_per_stage": float(self.tau_per_stage),
            "t_span": float(self.t_span),
            "n_steps": int(self.n_steps),
            "force_evaluations": int(self.force_evaluations),
            "force_evaluations_fsal": int(self.force_evaluations_fsal),
            "execution_wall_time": float(self.execution_wall_time),
            "kpis": self.kpi_report.kpis.to_dict(),
            "symplecticity_passed": bool(self.kpi_report.symplecticity_passed),
        }


class BenchmarkRunner:
    """
    Automated execution engine for celestial benchmarks and cross-integrator comparisons.
    """

    def __init__(self, G: float = GRAVITATIONAL_CONSTANT, use_optimized_kernels: bool = False) -> None:
        self.G = float(G)
        self.use_optimized_kernels = bool(use_optimized_kernels)

    def run(
        self,
        benchmark: Union[str, ExtendedPhaseState],
        integrator: SymplecticIntegrator,
        t_span: float,
        dt: float,
        auto_supernodes: bool = True,
        record_interval: int = 1,
        profile: Optional[str] = None,
        use_optimized_kernels: Optional[bool] = None,
    ) -> BenchmarkResult:
        """
        Executes a single benchmark run with a specified symplectic scheme.

        Parameters:
            benchmark: Name of the registered benchmark or explicit ExtendedPhaseState.
            integrator: SymplecticIntegrator instance.
            t_span: Total simulation time.
            dt: Macro time step size.
            auto_supernodes: Enables synchronous supernode regularized clustering.
            record_interval: Cadence for checkpoint history logging.
            use_optimized_kernels: Optional override activating Section 6.6 micro-kernels.

        Returns:
            BenchmarkResult instance.
        """
        if isinstance(benchmark, str):
            bench_name = benchmark
            state0 = get_benchmark(bench_name, G=self.G)
        else:
            bench_name = "custom_state"
            state0 = benchmark

        eff_opt = self.use_optimized_kernels if use_optimized_kernels is None else bool(use_optimized_kernels)
        system = NBodySystem(
            state=state0,
            integrator=integrator,
            G=self.G,
            desingularization_mode="auto" if auto_supernodes else "leaf_only",
            use_optimized_kernels=eff_opt,
        )
        system.history.record_interval = max(1, int(record_interval))

        t_start = time.perf_counter()
        final_state = system.evolve(t_span=t_span, dt=dt)
        wall_time = time.perf_counter() - t_start

        n_steps = int(round(t_span / dt)) if dt > 0 else 0
        stages = getattr(integrator, "evaluations_per_step", 1)
        dyadic_levels = system.history.get_dyadic_levels()
        if len(dyadic_levels) > 1 and system.history.record_interval == 1:
            total_substeps = int(np.sum(2 ** np.maximum(dyadic_levels[1:], 0)))
        else:
            total_substeps = n_steps
        force_evals = total_substeps * stages
        is_palindromic = getattr(integrator, "is_palindromic", False)
        fsal_stages = (stages - 1) if (is_palindromic and stages > 1) else max(1, int(stages))
        tau_per_stage = float(dt) / float(fsal_stages) if fsal_stages > 0 else float(dt)
        force_evals_fsal = (
            (1 + total_substeps * (stages - 1))
            if (is_palindromic and stages > 1 and total_substeps > 0)
            else force_evals
        )

        kpi_profile = profile or get_benchmark_profile(bench_name)
        kpi_report = compute_celestial_kpis(
            history=system.history,
            atlas_manager=system.atlas_manager,
            profile=kpi_profile,
        )

        return BenchmarkResult(
            benchmark_name=bench_name,
            integrator_name=integrator.name,
            integrator_order=integrator.order,
            stages_per_step=stages,
            dt=float(dt),
            t_span=float(t_span),
            n_steps=n_steps,
            force_evaluations=force_evals,
            execution_wall_time=wall_time,
            kpi_report=kpi_report,
            final_state=final_state,
            force_evaluations_fsal=force_evals_fsal,
            fsal_stages_per_step=fsal_stages,
            tau_per_stage=tau_per_stage,
        )

    def run_suite(
        self,
        benchmark: Union[str, ExtendedPhaseState],
        integrators: List[SymplecticIntegrator],
        t_span: float,
        dt: float,
        auto_supernodes: bool = True,
    ) -> Dict[str, BenchmarkResult]:
        """
        Executes a comparative benchmark across a suite of integrators.
        """
        results: Dict[str, BenchmarkResult] = {}
        for integ in integrators:
            res = self.run(
                benchmark=benchmark,
                integrator=integ,
                t_span=t_span,
                dt=dt,
                auto_supernodes=auto_supernodes,
            )
            results[integ.name] = res
        return results


def fit_convergence_order(
    x_vals: Union[List[float], np.ndarray],
    error_vals: Union[List[float], np.ndarray],
    mode: str = "step_size",
    noise_floor: float = 5e-15,
) -> Tuple[float, float, float]:
    """
    Fits asymptotic power-law convergence: log(error) = a * log(x) + b.

    Parameters:
        x_vals: Step sizes dt (for mode="step_size") or work W (for mode="work").
        error_vals: Relative energy errors max |Delta H / H_0|.
        mode: "step_size" (slope +p) or "work" (slope -p).
        noise_floor: Lower bound on error to avoid numerical roundoff saturation.

    Returns:
        Tuple (observed_order, intercept, r_squared).
    """
    x = np.asarray(x_vals, dtype=np.float64)
    y = np.asarray(error_vals, dtype=np.float64)

    valid_mask = (y > noise_floor) & (x > 0.0) & np.isfinite(x) & np.isfinite(y)
    if np.sum(valid_mask) < 2:
        return (0.0, 0.0, 0.0)

    log_x = np.log(x[valid_mask])
    log_y = np.log(y[valid_mask])

    if np.ptp(log_x) <= 1e-15:
        return (0.0, 0.0, 0.0)

    p_fit = np.polyfit(log_x, log_y, 1)
    slope = float(p_fit[0])
    intercept = float(p_fit[1])

    y_pred = slope * log_x + intercept
    ss_tot = np.sum((log_y - np.mean(log_y)) ** 2)
    ss_res = np.sum((log_y - y_pred) ** 2)
    r_squared = 1.0 - (ss_res / max(ss_tot, 1e-30)) if ss_tot > 1e-30 else 1.0

    order_p = slope if mode == "step_size" else -slope
    return (float(order_p), intercept, float(r_squared))


def compute_work_precision_pareto(
    benchmark_name: str,
    integrators: List[SymplecticIntegrator],
    dt_list: List[float],
    t_span: float = 1.0,
    G: float = GRAVITATIONAL_CONSTANT,
) -> Dict[str, Dict[str, Any]]:
    """
    Computes numerical work (force evaluations) versus relative energy error
    across multiple step sizes to produce Pareto convergence curves.

    Parameters:
        benchmark_name: Registered benchmark key.
        integrators: List of SymplecticIntegrators.
        dt_list: List of time step sizes to evaluate.
        t_span: Duration of each convergence test run.
        G: Gravitational constant.

    Returns:
        Dictionary mapping integrator name to {"work": array, "error": array, "order": int, ...}.
    """
    runner = BenchmarkRunner(G=G)
    pareto_data: Dict[str, Dict[str, Any]] = {}

    for integ in integrators:
        work_vals = []
        error_vals = []
        for dt in dt_list:
            res = runner.run(
                benchmark=benchmark_name,
                integrator=integ,
                t_span=t_span,
                dt=dt,
                auto_supernodes=False,  # Isolate pure integrator convergence
                record_interval=1,
            )
            work_vals.append(res.force_evaluations)
            err = res.kpi_report.kpis.delta_E_rel_max
            error_vals.append(max(err, 1e-16))

        obs_order_dt, _, r2_dt = fit_convergence_order(dt_list, error_vals, mode="step_size")
        obs_order_work, _, r2_work = fit_convergence_order(work_vals, error_vals, mode="work")

        pareto_data[integ.name] = {
            "work": np.array(work_vals, dtype=np.int64),
            "error": np.array(error_vals, dtype=np.float64),
            "order": integ.order,
            "observed_order": obs_order_dt,
            "r_squared": r2_dt,
            "work_slope": obs_order_work,
        }

    return pareto_data
