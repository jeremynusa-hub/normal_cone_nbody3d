"""
Unit tests for ComputationalProfiler and fine-grained stage cost instrumentation.
"""

import time
import numpy as np
import pytest

from normal_cone_nbody_float64.core.profiler import (
    ComputationalProfiler,
    get_current_profiler,
    set_current_profiler,
    reset_current_profiler,
)
from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
from normal_cone_nbody_float64.benchmarks.registry import create_burrau_pythagorean_benchmark


def test_profiler_time_stage_and_counters():
    prof = ComputationalProfiler(enabled=True)
    with prof.time_stage("kepler"):
        time.sleep(0.005)
        prof.record_kepler(halley_iters=3, stumpff_evals=6, is_bisection=False)

    summary = prof.get_summary()
    assert summary["total_cpu_time_sec"] >= 0.004
    assert "kepler" in summary["breakdown"]
    assert summary["breakdown"]["kepler"]["calls"] == 1
    assert summary["breakdown"]["kepler"]["counters"]["leaf_calls"] == 1
    assert summary["breakdown"]["kepler"]["counters"]["halley_iterations"] == 3
    assert summary["breakdown"]["kepler"]["counters"]["stumpff_evaluations"] == 6

    prof.record_supernode(substeps=8, acc_evals=56, is_close_pair=True, activations=1)
    summary2 = prof.get_summary()
    sn_counters = summary2["breakdown"]["supernode"]["counters"]
    assert sn_counters["symplectic_substeps"] == 8
    assert sn_counters["leapfrog_substeps"] == 8


def test_profiler_active_context_manager():
    prof = ComputationalProfiler(enabled=True)
    assert get_current_profiler() is None

    with prof.active():
        assert get_current_profiler() is prof
        prof.record_apex(kicks=2, pairs=3)

    assert get_current_profiler() is None
    summary = prof.get_summary()
    assert summary["breakdown"]["apex"]["counters"]["kick_evaluations"] == 2
    assert summary["breakdown"]["apex"]["counters"]["pairs_evaluated"] == 3


def test_nbody_system_profiler_integration():
    state = create_burrau_pythagorean_benchmark()
    prof = ComputationalProfiler(enabled=True)
    system = NBodySystem(
        state=state,
        integrator=BlanesMoan6Integrator(),
        profiler=prof,
    )

    # Step forward 3 macro-steps
    for _ in range(3):
        system.step(dt=0.01)

    assert len(prof.history) == 3
    summary = prof.get_summary()
    assert summary["total_macro_steps"] == 3
    assert summary["total_cpu_time_sec"] > 0.0

    # Ensure core conormal stages recorded work
    breakdown = summary["breakdown"]
    assert breakdown["dyadic"]["calls"] >= 3
    assert breakdown["atlas"]["calls"] >= 3
    assert breakdown["kepler"]["calls"] >= 1
    assert breakdown["apex"]["calls"] >= 1
    assert "kick_evaluations" in breakdown["apex"]["counters"]

    # Check step history properties
    p0 = prof.history[0]
    assert p0.step_index == 1
    assert p0.dt_macro == 0.01
    assert p0.r_min > 0.0
    assert "apex" in p0.stage_times
    assert "kepler" in p0.stage_times


def test_render_computational_cost_dashboard(tmp_path):
    from normal_cone_nbody_float64.diagnostics.visualizer import PublicationVisualizer
    state = create_burrau_pythagorean_benchmark()
    prof = ComputationalProfiler(enabled=True)
    system = NBodySystem(
        state=state,
        integrator=BlanesMoan6Integrator(),
        profiler=prof,
    )
    for _ in range(4):
        system.step(dt=0.01)

    viz = PublicationVisualizer(dpi=100)
    out_file = str(tmp_path / "test_cost_dash.png")
    fig = viz.render_computational_cost_dashboard(
        profiler=prof,
        history=system.history,
        title="Test Cost Dashboard",
        integrator_name="BM6",
        case_name="Burrau",
        save_path=out_file,
    )
    assert fig is not None
    import os
    assert os.path.exists(out_file)
    viz.close(fig)

