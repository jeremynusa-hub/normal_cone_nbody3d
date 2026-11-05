"""
Unit tests for Canonical Astrophysical N-Body Benchmark Registry and Runner in Float128.
Validates precision preservation, state instantiation, and cross-integrator benchmark execution.
"""

import pytest
import numpy as np
import mpmath as mp

from normal_cone_nbody_float128.benchmarks.registry import (
    get_benchmark,
    get_benchmark_profile,
    BENCHMARK_REGISTRY,
    create_eccentric_binary_benchmark,
)
from normal_cone_nbody_float128.benchmarks.runner import (
    BenchmarkRunner,
    compute_work_precision_pareto,
)
from normal_cone_nbody_float128.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float128.integrators.blanes_moan_6 import BlanesMoan6Integrator
from normal_cone_nbody_float128.integrators.processed import ProcessedSymplecticIntegrator
from normal_cone_nbody_float128.integrators.correctors import CorrectorC4
from normal_cone_nbody_float128.core.constants import G_DEFAULT_F128


def test_benchmark_registry_f128_instantiation():
    """Validates that get_benchmark instantiates both float64 and float128 states without precision loss."""
    # 1. Burrau Pythagorean
    state_f64 = get_benchmark("burrau_pythagorean", is_f128=False)
    assert not state_f64.is_f128
    assert state_f64.q.dtype == np.float64
    assert isinstance(state_f64.E0, float)

    state_f128 = get_benchmark("burrau_pythagorean", is_f128=True)
    assert state_f128.is_f128
    assert state_f128.q.dtype == object
    assert isinstance(state_f128.E0, mp.mpf)
    assert state_f128.N == 3

    # 2. Eccentric Binary
    bin_f128 = get_benchmark("eccentric_binary", is_f128=True, e=0.95, a=1.0)
    assert bin_f128.is_f128
    assert bin_f128.N == 2
    # Verify initial apoapsis distance r_max = a*(1+e) = 1.95
    r_apo = mp.mpf(bin_f128.q[1, 0]) - mp.mpf(bin_f128.q[0, 0])
    assert mp.almosteq(r_apo, mp.mpf('1.95'), rel_eps=mp.mpf('1e-30'))
    # Verify 35-digit exactness of apoapsis relative velocity v_apo^2 = G*(m0+m1)*(1-e)/(a*(1+e))
    v_rel_y = (mp.mpf(bin_f128.p[1, 1]) / mp.mpf(bin_f128.masses[1])) - (
        mp.mpf(bin_f128.p[0, 1]) / mp.mpf(bin_f128.masses[0])
    )
    v_apo_sq_exact = mp.mpf('2.0') * (mp.mpf('1') - mp.mpf('0.95')) / mp.mpf('1.95')
    assert abs(v_rel_y ** 2 - v_apo_sq_exact) < mp.mpf('1e-33'), (
        f"Expected 35-digit exact v_apo^2, got residual {abs(v_rel_y ** 2 - v_apo_sq_exact)}"
    )

    # 3. Jovian Laplace Resonance in 35-digit precision (default is_f128=True and forced equilibrium)
    lap_f128 = get_benchmark("jovian_laplace_resonance")
    assert lap_f128.is_f128
    assert lap_f128.N == 4
    p_tot_norm = mp.sqrt(sum(sum(mp.mpf(lap_f128.p[i, d]) for i in range(4)) ** 2 for d in range(3)))
    assert p_tot_norm < mp.mpf('1e-33')
    # Also verify unforced option a2 / a1 = 2^(2/3) to 35-digit precision
    lap_unforced = get_benchmark("jovian_laplace_resonance", use_forced_equilibrium=False)
    assert lap_unforced.is_f128
    a1_f128 = abs(mp.mpf(lap_unforced.q[1, 0]) - mp.mpf(lap_unforced.q[0, 0]))
    a2_f128 = abs(mp.mpf(lap_unforced.q[2, 0]) - mp.mpf(lap_unforced.q[0, 0]))
    assert abs((a2_f128 / a1_f128) ** 3 - mp.mpf('4.0')) < mp.mpf('1e-33')

    # 4. Verify that ALL registered benchmarks in Float128 return is_f128 == True by default (without passing is_f128=True)
    for b_name in BENCHMARK_REGISTRY.keys():
        st_f128 = get_benchmark(b_name)
        assert st_f128.is_f128 is True, f"Benchmark '{b_name}' did not default to is_f128=True"
        n_b = st_f128.N
        p_norm = mp.sqrt(sum(sum(mp.mpf(st_f128.p[i, d]) for i in range(n_b)) ** 2 for d in range(3)))
        cdm_norm = mp.sqrt(sum(sum(mp.mpf(st_f128.masses[i]) * mp.mpf(st_f128.q[i, d]) for i in range(n_b)) ** 2 for d in range(3)))
        assert p_norm < mp.mpf('1e-33'), f"{b_name} momentum residual {p_norm}"
        assert cdm_norm < mp.mpf('1e-33'), f"{b_name} barycenter residual {cdm_norm}"



def test_benchmark_runner_execution_f64_and_f128():
    """Validates that BenchmarkRunner runs benchmarks seamlessly under both float64 and float128."""
    bm4 = BlanesMoan4Integrator()
    runner_f64 = BenchmarkRunner(is_f128=False)
    runner_f128 = BenchmarkRunner(is_f128=True)

    # Short test run on eccentric binary
    res_f64 = runner_f64.run(
        benchmark="eccentric_binary",
        integrator=bm4,
        t_span=0.04,
        dt=0.02,
        auto_supernodes=False,
    )
    assert not res_f64.final_state.is_f128
    assert res_f64.n_steps == 2
    assert res_f64.force_evaluations == 2 * 7
    assert res_f64.kpi_report.kpis.delta_E_rel_max < 1e-4

    res_f128 = runner_f128.run(
        benchmark="eccentric_binary",
        integrator=bm4,
        t_span=0.04,
        dt=0.02,
        auto_supernodes=False,
        is_f128=True,
    )
    assert res_f128.final_state.is_f128
    assert res_f128.n_steps == 2
    assert res_f128.force_evaluations == 2 * 7
    assert res_f128.kpi_report.kpis.delta_E_rel_max < 1e-4


def test_work_precision_pareto_runner_f128():
    """Validates compute_work_precision_pareto with float128 integration."""
    bm4 = BlanesMoan4Integrator()
    pareto = compute_work_precision_pareto(
        benchmark_name="eccentric_binary",
        integrators=[bm4],
        dt_list=[0.02, 0.01],
        t_span=0.02,
        is_f128=True,
    )
    assert bm4.name in pareto
    assert len(pareto[bm4.name]["work"]) == 2
    assert len(pareto[bm4.name]["error"]) == 2
