"""
Unit Tests for Layer 7: Benchmark Runner and Pareto Work-Precision Analyzer (runner.py).

Theoretical Alignments:
- Monograph: Chapter 7 (Global Numerical Validation, Work-Precision Scaling).
- Invariants: Force evaluation stage tracking, Pareto data integrity,
  and higher-order asymptotic superiority.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.benchmarks.runner import (
    BenchmarkRunner,
    BenchmarkResult,
    compute_work_precision_pareto,
)
from normal_cone_nbody_float64.integrators.strang import StrangIntegrator
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator


def test_benchmark_runner_single_execution():
    """
    Validates executing a single benchmark run with BenchmarkRunner:
    - Accurate force evaluation count (n_steps * stages).
    - Wall time tracking > 0.
    - KPI certification report generation.
    """
    runner = BenchmarkRunner(G=1.0)
    integrator = BlanesMoan4Integrator(G=1.0)

    result = runner.run(
        benchmark="pareto_test",
        integrator=integrator,
        t_span=0.1,
        dt=0.02,
        auto_supernodes=False,
    )

    assert isinstance(result, BenchmarkResult)
    assert result.benchmark_name == "pareto_test"
    assert result.integrator_name == "BlanesMoan-BM4"
    assert result.integrator_order == 4
    assert result.n_steps == 5
    assert result.stages_per_step == 7
    assert result.force_evaluations == 5 * 7  # 35 evaluations
    assert result.execution_wall_time > 0.0
    assert result.kpi_report.symplecticity_passed is True

    d = result.to_dict()
    assert d["n_steps"] == 5
    assert d["force_evaluations"] == 35


def test_benchmark_runner_suite_comparison():
    """
    Validates comparing multiple integrators on the same benchmark.
    """
    runner = BenchmarkRunner(G=1.0)
    integrators = [StrangIntegrator(G=1.0), BlanesMoan4Integrator(G=1.0)]

    suite_results = runner.run_suite(
        benchmark="pareto_test",
        integrators=integrators,
        t_span=0.1,
        dt=0.02,
    )

    assert "Strang-S2" in suite_results
    assert "BlanesMoan-BM4" in suite_results
    assert suite_results["Strang-S2"].force_evaluations == 5 * 2
    assert suite_results["BlanesMoan-BM4"].force_evaluations == 5 * 7


def test_compute_work_precision_pareto_order_superiority():
    """
    Validates work-precision Pareto data computation across Strang and BM4.
    Confirms that at smaller step sizes, BM4 delivers significantly higher
    Hamiltonian precision than Strang.
    """
    integrators = [StrangIntegrator(G=1.0), BlanesMoan4Integrator(G=1.0)]
    dt_list = [0.04, 0.02]

    pareto = compute_work_precision_pareto(
        benchmark_name="pareto_test",
        integrators=integrators,
        dt_list=dt_list,
        t_span=0.2,
        G=1.0,
    )

    assert "Strang-S2" in pareto
    assert "BlanesMoan-BM4" in pareto

    strang_data = pareto["Strang-S2"]
    bm4_data = pareto["BlanesMoan-BM4"]

    assert len(strang_data["work"]) == 2
    assert len(bm4_data["work"]) == 2

    # Verify that at dt = 0.02, BM4 error is substantially lower than Strang error
    err_strang = strang_data["error"][-1]
    err_bm4 = bm4_data["error"][-1]
    assert err_bm4 < err_strang / 10.0, f"Expected BM4 error ({err_bm4}) < Strang ({err_strang}) / 10"


def test_benchmark_runner_multirate_force_evaluations_and_no_deprecation_warning():
    """
    Validates that BenchmarkRunner:
    1. Does not emit any deprecation UserWarning when called with auto_supernodes.
    2. Accurately counts force_evaluations = sum(2^k) * stages when dyadic substepping k >= 1 occurs.
    """
    import warnings
    from normal_cone_nbody_float64.core.state import ExtendedPhaseState

    # 3-body close encounter state that triggers dyadic substepping k >= 1 under large macro dt
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    q = np.array([
        [-0.08, 0.0, 0.0],
        [0.08, 0.0, 0.0],
        [0.0, 0.25, 0.0],
    ], dtype=np.float64)
    q -= np.sum(masses[:, None] * q, axis=0) / np.sum(masses)
    p = np.zeros_like(q)
    state0 = ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0)

    runner = BenchmarkRunner(G=1.0)
    integrator = StrangIntegrator(G=1.0)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _ = runner.run(
            benchmark=state0,
            integrator=integrator,
            t_span=0.05,
            dt=0.05,
            auto_supernodes=False,
        )
        res = runner.run(
            benchmark=state0,
            integrator=integrator,
            t_span=0.05,
            dt=0.05,
            auto_supernodes=True,
        )
    user_warnings = [w for w in caught if issubclass(w.category, UserWarning)]
    assert len(user_warnings) == 0, f"Unexpected UserWarnings emitted: {user_warnings}"

    # Verify that dyadic level k >= 1 was triggered and force_evaluations reflects sum(2^k) * stages
    assert res.n_steps == 1
    assert res.force_evaluations > res.n_steps * res.stages_per_step

