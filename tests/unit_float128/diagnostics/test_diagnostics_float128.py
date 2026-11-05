"""
Unit Tests for Extended Precision (Float128 / Mpmath 35 dps) Diagnostics & KPIs.

Theoretical Alignments:
- Monograph: Chapter 5 Section 5.3 & Chapter 6 Global Symplectic Metrology.
- Validates compute_celestial_kpis, Brouwer diffusion envelope, and CollisionRadar
  under high-precision mpmath representations.
"""

import numpy as np
import pytest
import mpmath as mp

from normal_cone_nbody_float128.core.state import ExtendedPhaseState
from normal_cone_nbody_float128.core.invariants import compute_invariants
from normal_cone_nbody_float128.system.history import HistoryTracker
from normal_cone_nbody_float128.diagnostics.metrics import (
    compute_celestial_kpis,
    compute_brouwer_diffusion_series,
    compute_inertia_radius_and_dimensionless_invariants,
    STANDARD_KPI_PROFILES,
)
from normal_cone_nbody_float128.diagnostics.collision_radar import (
    CollisionRadar,
    predict_pair_approach,
)


def test_celestial_kpis_extended_precision_profile():
    """
    Validates compute_celestial_kpis with 'extended_precision' profile on float128 history.
    """
    tracker = HistoryTracker(record_interval=1)
    masses = np.array([mp.mpf("1.0"), mp.mpf("1.0"), mp.mpf("1.0")], dtype=object)

    # 3 bodies in symmetric Lagrange configuration with exact energy
    # E0 = -1.5, P = 0, L = constant
    for step in range(6):
        t = mp.mpf(step) * mp.mpf("0.1")
        # Symmetric equilateral triangle
        angle = mp.mpf(step) * mp.mpf("0.05")
        cos_a = mp.cos(angle)
        sin_a = mp.sin(angle)

        q = np.array([
            [cos_a, sin_a, mp.mpf("0.0")],
            [-mp.mpf("0.5") * cos_a - mp.sqrt(3)/2 * sin_a, -mp.mpf("0.5") * sin_a + mp.sqrt(3)/2 * cos_a, mp.mpf("0.0")],
            [-mp.mpf("0.5") * cos_a + mp.sqrt(3)/2 * sin_a, -mp.mpf("0.5") * sin_a - mp.sqrt(3)/2 * cos_a, mp.mpf("0.0")],
        ], dtype=object)

        # Velocities orthogonal to positions
        p = np.array([
            [-mp.mpf("0.5") * sin_a, mp.mpf("0.5") * cos_a, mp.mpf("0.0")],
            [mp.mpf("0.25") * sin_a - mp.mpf("0.25") * mp.sqrt(3) * cos_a, -mp.mpf("0.25") * cos_a - mp.mpf("0.25") * mp.sqrt(3) * sin_a, mp.mpf("0.0")],
            [mp.mpf("0.25") * sin_a + mp.mpf("0.25") * mp.sqrt(3) * cos_a, -mp.mpf("0.25") * cos_a + mp.mpf("0.25") * mp.sqrt(3) * sin_a, mp.mpf("0.0")],
        ], dtype=object)

        state = ExtendedPhaseState(q=q, p=p, masses=masses, t=t)
        tracker.record(state)

    report = compute_celestial_kpis(tracker, profile="extended_precision")

    assert report.symplecticity_passed is True
    assert report.energy_conserved is True
    assert report.momentum_conserved is True
    assert report.angular_momentum_conserved is True
    assert report.galilean_boost_conserved is True

    # Check precision levels
    assert report.kpis.delta_E_rel_max < 1e-25
    assert report.kpis.delta_P_max < 1e-28
    assert report.kpis.delta_L_max < 1e-28


def test_brouwer_diffusion_and_inertia_radius_float128():
    """
    Validates Brouwer diffusion calculation scaling at u_mach ~ 1.93e-34 for float128.
    """
    tracker = HistoryTracker(record_interval=1)
    masses = np.array([mp.mpf("1.0"), mp.mpf("1.0")], dtype=object)

    for step in range(5):
        q = np.array([[mp.mpf("1.0"), mp.mpf("0.0"), mp.mpf("0.0")],
                      [mp.mpf("-1.0"), mp.mpf("0.0"), mp.mpf("0.0")]], dtype=object)
        p = np.array([[mp.mpf("0.0"), mp.mpf("0.5"), mp.mpf("0.0")],
                      [mp.mpf("0.0"), mp.mpf("-0.5"), mp.mpf("0.0")]], dtype=object)
        state = ExtendedPhaseState(q=q, p=p, masses=masses, t=mp.mpf(step) * mp.mpf("0.1"))
        tracker.record(state)

    sigma_series, max_brouwer_bound, ratio_p = compute_brouwer_diffusion_series(tracker)
    assert len(sigma_series) == 5
    # Brouwer bound in float128 should be near 1e-32
    assert max_brouwer_bound < 1e-25

    R_I_series, eps_L, eps_G, max_R_I, max_eps_L, max_eps_G = (
        compute_inertia_radius_and_dimensionless_invariants(tracker)
    )
    assert len(R_I_series) == 5
    assert abs(max_R_I - 1.0) < 1e-12
    assert max_eps_L < 1e-25
    assert max_eps_G < 1e-25


def test_collision_radar_with_float128_state():
    """
    Validates that CollisionRadar handles ExtendedPhaseState(is_f128=True) seamlessly.
    """
    radar = CollisionRadar()
    masses = np.array([mp.mpf("1.0"), mp.mpf("2.0"), mp.mpf("0.5")], dtype=object)
    q = np.array([
        [mp.mpf("0.1"), mp.mpf("0.0"), mp.mpf("0.0")],
        [mp.mpf("-0.1"), mp.mpf("0.0"), mp.mpf("0.0")],
        [mp.mpf("5.0"), mp.mpf("5.0"), mp.mpf("0.0")],
    ], dtype=object)
    p = np.array([
        [mp.mpf("-0.1"), mp.mpf("0.0"), mp.mpf("0.0")],  # Approaching body 1
        [mp.mpf("0.1"), mp.mpf("0.0"), mp.mpf("0.0")],
        [mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")],
    ], dtype=object)
    state = ExtendedPhaseState(q=q, p=p, masses=masses, t=mp.mpf("0.0"))

    report = radar.scan(state)

    assert report.closest_pair == (0, 1)
    assert abs(report.min_separation - 0.2) < 1e-12
    assert 0.0 <= report.hazard_score <= 1.0
    assert len(report.pair_predictions) == 3
    # Pair (0, 1) is approaching
    pred_01 = next(p for p in report.pair_predictions if p.body_i == 0 and p.body_j == 1)
    assert pred_01.is_approaching is True
    assert pred_01.radial_velocity < 0.0


def test_radar_alias_and_visualizer_float128(tmp_path):
    """
    Validates that float128 diagnostics.radar imports cleanly and that
    PublicationVisualizer renders 10-panel, 6-panel, and 8-panel dashboards
    on a HistoryTracker(is_f128=True) without TypeError (D-146-04, D-146-05).
    """
    import inspect
    from normal_cone_nbody_float128.diagnostics.radar import (
        CollisionRadar as RadarAlias,
        CollisionWarning,
        predict_closest_approach,
        compute_pairwise_free_fall_times,
    )
    from normal_cone_nbody_float128.diagnostics import visualizer as vis_mod
    from normal_cone_nbody_float128.diagnostics.visualizer import PublicationVisualizer

    assert RadarAlias is CollisionRadar
    assert "processed_bm6" not in inspect.getsource(vis_mod)

    tracker = HistoryTracker(record_interval=1)
    masses = np.array([mp.mpf("1.0"), mp.mpf("1.0"), mp.mpf("0.5")], dtype=object)
    for step in range(4):
        t = mp.mpf(step) * mp.mpf("0.05")
        q = np.array([
            [mp.mpf("-0.5") + mp.mpf("0.01") * step, mp.mpf("0.0"), mp.mpf("0.0")],
            [mp.mpf("0.5") - mp.mpf("0.01") * step, mp.mpf("0.0"), mp.mpf("0.0")],
            [mp.mpf("0.0"), mp.mpf("4.0"), mp.mpf("0.0")],
        ], dtype=object)
        p = np.array([
            [mp.mpf("0.0"), mp.mpf("-0.6"), mp.mpf("0.0")],
            [mp.mpf("0.0"), mp.mpf("0.6"), mp.mpf("0.0")],
            [mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")],
        ], dtype=object)
        st = ExtendedPhaseState(q=q, p=p, masses=masses, t=t)
        tracker.record(st)
        if step == 0:
            ff = compute_pairwise_free_fall_times(st)
            assert (0, 1) in ff

    vis = PublicationVisualizer(dpi=72)
    fig10 = vis.render_10panel_dashboard(tracker, title="Float128 10-Panel Test")
    vis.close(fig10)

    fig6 = vis.render_multi_integrator_invariants_dashboard({"blanes_moan_4": tracker}, title="Float128 6-Panel Test")
    vis.close(fig6)

    fig8 = vis.render_refined_dashboard(tracker, title="Float128 8-Panel Test")
    vis.close(fig8)


def test_laskar_naff_and_forward_backward_reversibility_float128():
    """
    Validates Laskar's NAFF frequency drift, 35-digit Forward-Backward time-reversibility
    (< 1e-32), and Laskar's tau_per_stage metric in float128.
    """
    from normal_cone_nbody_float128.diagnostics.metrics import (
        extract_laskar_fundamental_frequency,
        compute_laskar_frequency_drift,
        compute_forward_backward_reversibility,
    )
    from normal_cone_nbody_float128.benchmarks.runner import BenchmarkRunner
    from normal_cone_nbody_float128.integrators.processed import create_pbm4_integrator

    # 1. Exact analytical signal frequency extraction
    times = np.linspace(0.0, 20.0, 1025)
    omega0 = 3.141592653589793
    sig = np.exp(1j * omega0 * times)
    assert abs(extract_laskar_fundamental_frequency(times, sig) - omega0) < 1e-12

    # 2. 35-digit Forward-Backward time-reversibility (< 1e-32) under PBM4
    masses = np.array([mp.mpf("1.0"), mp.mpf("1.0"), mp.mpf("0.2")], dtype=object)
    q0 = np.array([
        [mp.mpf("-0.5"), mp.mpf("0.0"), mp.mpf("0.0")],
        [mp.mpf("0.5"), mp.mpf("0.0"), mp.mpf("0.0")],
        [mp.mpf("0.0"), mp.mpf("5.0"), mp.mpf("0.0")],
    ], dtype=object)
    p0 = np.array([
        [mp.mpf("0.0"), mp.mpf("-0.6"), mp.mpf("0.0")],
        [mp.mpf("0.0"), mp.mpf("0.6"), mp.mpf("0.0")],
        [mp.mpf("-0.08"), mp.mpf("0.0"), mp.mpf("0.0")],
    ], dtype=object)
    s0 = ExtendedPhaseState(q=q0, p=p0, masses=masses, t=mp.mpf("0.0"))
    pbm4 = create_pbm4_integrator(G=mp.mpf("1.0"))

    rev = compute_forward_backward_reversibility(
        state0=s0,
        integrator=pbm4,
        dt=mp.mpf("0.02"),
        n_steps=6,
        G=mp.mpf("1.0"),
    )
    assert rev["delta_q_norm"] < 1e-32
    assert rev["delta_p_norm"] < 1e-32
    assert rev["delta_E_rel"] < 1e-32
    assert rev["dyadic_palindrome_matched"] is True

    # 3. BenchmarkRunner tau_per_stage in float128
    runner = BenchmarkRunner(G=mp.mpf("1.0"), is_f128=True)
    res = runner.run(benchmark=s0, integrator=pbm4, t_span=mp.mpf("0.04"), dt=mp.mpf("0.02"), auto_supernodes=False, is_f128=True)
    assert res.fsal_stages_per_step == 6
    assert res.tau_per_stage == pytest.approx(0.02 / 6.0, rel=1e-14)

    # 4. Osculating elements and Laplace resonant angle on forced equilibrium in float128
    from normal_cone_nbody_float128.benchmarks.registry import get_benchmark
    from normal_cone_nbody_float128.diagnostics.metrics import (
        compute_osculating_elements,
        compute_laplace_resonant_angle,
    )

    lap_s0 = get_benchmark("jovian_laplace_resonance")
    assert lap_s0.is_f128 is True
    el_io = compute_osculating_elements(lap_s0, primary_idx=0, body_idx=1)
    assert abs(el_io["e"] - 0.0041) < 1e-4
    lap_ang_rad = compute_laplace_resonant_angle(lap_s0)
    assert abs(np.degrees(lap_ang_rad) - 180.0) < 1e-9


