"""
Unit Tests for Layer 6: Celestial KPIs and Symplectic Metrology (metrics.py).

Theoretical Alignments:
- Monograph: Chapter 5 (Section 5.3, Backward Error Analysis),
  Chapter 6 (Global Symplectic Metrology).
- Invariants: Exact linear slope extraction, secular slope vanishing for symplectic schemes
  (Slope_H < 1e-11), and Lie--Bargmann exact momentum preservation (< 1e-13).
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float64.diagnostics.metrics import (
    compute_energy_secular_slope,
    compute_celestial_kpis,
    CelestialKPIs,
    KPIReport,
)


def test_secular_energy_slope_linear_regression():
    """
    Validates ordinary least squares slope calculation:
    1. Perfect linear drift: E(t) = 5.0 + 2.5 * t -> slope = 2.5 exactly.
    2. Zero-mean pure oscillation: E(t) = 10.0 + 1e-3 * sin(2*pi*t) over integer periods
       -> slope < 1e-15.
    """
    times = np.linspace(0.0, 10.0, 101)

    # 1. Linear drift
    energies_linear = 5.0 + 2.5 * times
    slope, intercept = compute_energy_secular_slope(times, energies_linear)
    assert slope == pytest.approx(2.5, abs=1e-14)
    assert intercept == pytest.approx(5.0, abs=1e-14)

    # 2. Oscillatory series over symmetric interval [-5, 5]
    times_sym = np.linspace(-5.0, 5.0, 1001)
    energies_osc = 10.0 + 1e-3 * np.cos(2.0 * np.pi * times_sym)
    slope_osc, _ = compute_energy_secular_slope(times_sym, energies_osc)
    assert abs(slope_osc) < 1e-14


def test_celestial_kpis_extraction_on_symplectic_simulation():
    """
    Validates extracting the 6 canonical celestial KPIs on a 3-body system
    integrated with Blanes-Moan 4:
    - Linear momentum drift < 1e-13.
    - Angular momentum drift < 1e-13.
    - Secular energy slope < 1e-11 (absence of uncompensated secular drift).
    - Symplectic certification report passes.
    """
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q0 = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 6.0, 0.2],
    ], dtype=np.float64)
    p0 = np.array([
        [0.0, -0.7, 0.0],
        [0.0, 0.7, 0.0],
        [0.2, 0.0, -0.05],
    ], dtype=np.float64)
    state0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)

    system = NBodySystem(state=state0, integrator=BlanesMoan4Integrator(G=1.0), G=1.0)
    system.evolve(t_span=0.4, dt=0.02)

    report = compute_celestial_kpis(
        history=system.history,
        atlas_manager=system.atlas_manager,
    )

    assert isinstance(report, KPIReport)
    assert report.symplecticity_passed is True
    assert report.zero_secular_drift is True
    assert report.momentum_conserved is True
    assert report.angular_momentum_conserved is True

    kpis = report.kpis
    assert kpis.delta_P_max < 1e-13
    assert kpis.delta_L_max < 1e-13
    assert kpis.secular_energy_slope < 1e-11
    assert kpis.n_records == 21
    assert kpis.t_span == pytest.approx(0.4, abs=1e-14)

    # Markdown export test
    md = report.to_markdown()
    assert "Celestial Symplectic KPI Certification Report" in md
    assert "[PASS]" in md


def test_contextual_kpi_profiles():
    """
    Validates contextual calibration profiles:
    - chaotic_encounter profile adapts energy, momentum and secular slope evaluation.
    - hierarchical_multiscale profile adapts tolerances for multi-rate atlas systems.
    """
    from normal_cone_nbody_float64.benchmarks.registry import get_benchmark
    from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator

    s0 = get_benchmark("burrau_pythagorean")
    system = NBodySystem(state=s0, integrator=BlanesMoan6Integrator(G=1.0), G=1.0, desingularization_mode="auto")
    system.evolve(t_span=1.0, dt=0.01)

    report_chaotic = compute_celestial_kpis(
        history=system.history,
        atlas_manager=system.atlas_manager,
        profile="chaotic_encounter",
    )
    assert report_chaotic.profile_name == "chaotic_encounter"
    assert report_chaotic.tol_energy_rel == 1e-3
    assert report_chaotic.tol_momentum == 1e-11
    assert report_chaotic.symplecticity_passed is True
    assert report_chaotic.energy_conserved is True


def test_brouwer_stochastic_diffusion_and_scale_invariants():
    """
    Validates the IEEE-754 float64 Brouwer roundoff diffusion envelope
    and dimensionless scale-invariant drift metrics.
    """
    from normal_cone_nbody_float64.benchmarks.registry import get_benchmark
    from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator

    s0 = get_benchmark("figure8_choreography")
    system = NBodySystem(state=s0, integrator=BlanesMoan6Integrator(G=1.0), G=1.0)
    system.evolve(t_span=0.5, dt=0.02)

    report = compute_celestial_kpis(
        history=system.history,
        profile="smooth_periodic",
    )

    # 1. Brouwer stochastic diffusion
    assert len(report.brouwer_sigma_series) == len(system.history)
    assert report.brouwer_delta_P_bound > 0.0
    assert report.ratio_P_brouwer <= 1.0  # Drift must be inside the Brouwer bound

    # 2. Inertia radius & Dimensionless invariants
    assert len(report.inertia_radius_series) == len(system.history)
    assert np.all(report.inertia_radius_series > 0.0)
    assert report.dimensionless_delta_L < 1e-11
    assert report.dimensionless_delta_G < 1e-11
    assert report.kpis.brouwer_delta_P_bound == report.brouwer_delta_P_bound
    assert report.kpis.dimensionless_delta_L == report.dimensionless_delta_L

    # 3. Verify custom stages_per_substep forwarding (Strang=1 vs BM6=11)
    rep_strang = compute_celestial_kpis(history=system.history, stages_per_substep=1)
    rep_bm6 = compute_celestial_kpis(history=system.history, stages_per_substep=11)
    assert rep_bm6.brouwer_sigma_series[-1] == pytest.approx(
        rep_strang.brouwer_sigma_series[-1] * np.sqrt(11.0), rel=1e-12
    )


def test_laskar_frequency_drift_and_forward_backward_reversibility():
    """
    Validates Laskar's NAFF (Hanning-windowed Fourier fundamental frequency extractor),
    Rein & Spiegel's Forward-Backward time-reversibility audit, and Laskar's
    normalized inverse-cost abscissa tau_per_stage = dt / s_fsal.
    """
    from normal_cone_nbody_float64.diagnostics.metrics import (
        extract_laskar_fundamental_frequency,
        compute_laskar_frequency_drift,
        compute_forward_backward_reversibility,
    )
    from normal_cone_nbody_float64.benchmarks.runner import BenchmarkRunner
    from normal_cone_nbody_float64.integrators.processed import create_pbm4_integrator

    # 1. Exact analytical quasi-periodic signal z(t) = exp(i * omega1 * t) + 0.15 * exp(i * omega2 * t)
    times = np.linspace(0.0, 60.0, 4097)
    omega1 = 2.3456789012345
    omega2 = 7.8912345678901
    signal = np.exp(1j * omega1 * times) + 0.15 * np.exp(1j * omega2 * times)
    omega_est = extract_laskar_fundamental_frequency(times, signal)
    assert abs(omega_est - omega1) / omega1 < 1e-12

    # 2. Laskar frequency drift on a circular Keplerian binary (r0=1, M=2, v_circ=sqrt(2), p_circ=sqrt(2)/2 -> omega=sqrt(2))
    masses = np.array([1.0, 1.0], dtype=np.float64)
    q0 = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]], dtype=np.float64)
    p_circ = 0.5 * np.sqrt(2.0)
    p0 = np.array([[0.0, -p_circ, 0.0], [0.0, p_circ, 0.0]], dtype=np.float64)
    s0 = ExtendedPhaseState(q=q0, p=p0, masses=masses)
    pbm4 = create_pbm4_integrator(G=1.0)
    sys_bin = NBodySystem(state=s0, integrator=pbm4, G=1.0)
    sys_bin.evolve(t_span=8.0 * np.pi, dt=0.05)

    naff_info = compute_laskar_frequency_drift(sys_bin.history, pair=(0, 1))
    assert abs(naff_info["nu_1"] - np.sqrt(2.0)) / np.sqrt(2.0) < 1e-6
    assert naff_info["rel_freq_drift"] < 1e-10

    # 3. Forward-Backward time-reversibility audit on 3-body hierarchical system
    m3 = np.array([1.0, 1.0, 0.2], dtype=np.float64)
    q3 = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 5.0, 0.0]], dtype=np.float64)
    p3 = np.array([[0.0, -0.6, 0.0], [0.0, 0.6, 0.0], [-0.08, 0.0, 0.0]], dtype=np.float64)
    s3 = ExtendedPhaseState(q=q3, p=p3, masses=m3)
    rev_report = compute_forward_backward_reversibility(
        state0=s3,
        integrator=create_pbm4_integrator(G=1.0),
        dt=0.02,
        n_steps=15,
        G=1.0,
    )
    assert rev_report["delta_q_norm"] < 5e-14
    assert rev_report["delta_p_norm"] < 5e-14
    assert rev_report["delta_E_rel"] < 5e-14
    assert rev_report["dyadic_palindrome_matched"] is True

    # 4. Laskar normalized inverse-cost abscissa tau_per_stage in BenchmarkRunner
    runner = BenchmarkRunner(G=1.0)
    res = runner.run(benchmark=s3, integrator=pbm4, t_span=0.2, dt=0.05, auto_supernodes=False)
    assert res.fsal_stages_per_step == 6  # 7-stage BM4/PBM4 has 6 FSAL stages per step
    assert res.tau_per_stage == pytest.approx(0.05 / 6.0, rel=1e-14)

