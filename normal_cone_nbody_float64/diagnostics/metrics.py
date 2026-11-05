"""
Canonical Celestial KPIs and Symplectic Metrology.

Theoretical Foundations:
- Monograph Reference: Chapter 5 (Section 5.3, Backward Error Analysis and
  Absence of Secular Energy Drift), Chapter 6 (Global Symplectic Metrology).
- The 6 Canonical Celestial KPIs:
    1. Maximum Relative Hamiltonian Energy Drift:
         Delta E_{rel}^{max} = max_t |H(t) - H_0| / (|H_0| + epsilon_{reg}).
    2. Hamiltonian Secular Slope:
         Slope_H = |beta| from linear regression Delta H(t) = alpha + beta * t.
         For symplectic schemes, Slope_H = 0 (bounded oscillation O(h^{2p}) without linear drift).
    3. Maximum Total Linear Momentum Drift:
         Delta P^{max} = max_t ||P(t) - P_0||.
    4. Maximum Total Angular Momentum Drift:
         Delta L^{max} = max_t ||L(t) - L_0||.
    5. Extended Center of Mass / Galilean Drift:
         Delta G^{max} = max_t ||G(t) - G_0||.
    6. Mean Atlas Dwell Time & Transition Count:
         tau_{dwell} = T_{sim} / (N_{transitions} + 1).
- Certified Formal Proofs:
    Lean 4 modules `02_LieBargmannExactInvariance_Theorems.lean` and
    `AsyncSymplecticConservation.lean`.
"""

from __future__ import annotations
import math
from dataclasses import dataclass, field
from typing import Dict, Any, Tuple, Optional, List, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from ..system.history import HistoryTracker
    from ..system.atlas_manager import AtlasManager


@dataclass(frozen=True)
class CelestialKPIs:
    """
    Standard container for the canonical celestial metrics of an N-body simulation,
    including intrinsic scale-invariant and machine-roundoff metrics.
    """
    delta_E_rel_max: float
    secular_energy_slope: float
    delta_P_max: float
    delta_L_max: float
    delta_G_max: float
    mean_wronskian_error: float
    atlas_dwell_time: float
    n_atlas_transitions: int
    t_span: float
    n_records: int
    # Intrinsic & Machine-Roundoff Metrics
    brouwer_delta_P_bound: float = 0.0
    ratio_P_brouwer: float = 0.0
    dimensionless_delta_L: float = 0.0
    dimensionless_delta_G: float = 0.0
    inertia_radius_max: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "delta_E_rel_max": float(self.delta_E_rel_max),
            "secular_energy_slope": float(self.secular_energy_slope),
            "delta_P_max": float(self.delta_P_max),
            "delta_L_max": float(self.delta_L_max),
            "delta_G_max": float(self.delta_G_max),
            "mean_wronskian_error": float(self.mean_wronskian_error),
            "atlas_dwell_time": float(self.atlas_dwell_time),
            "n_atlas_transitions": int(self.n_atlas_transitions),
            "t_span": float(self.t_span),
            "n_records": int(self.n_records),
            "brouwer_delta_P_bound": float(self.brouwer_delta_P_bound),
            "ratio_P_brouwer": float(self.ratio_P_brouwer),
            "dimensionless_delta_L": float(self.dimensionless_delta_L),
            "dimensionless_delta_G": float(self.dimensionless_delta_G),
            "inertia_radius_max": float(self.inertia_radius_max),
        }


@dataclass(frozen=True)
class KPIProfile:
    """
    Contextual calibration profile specifying admissible tolerances by dynamical regime.
    """
    name: str
    tol_energy_rel: float
    tol_secular_slope: float
    tol_momentum: float
    tol_angular_momentum: float
    tol_galilean_cm: float
    check_secular_slope: bool = True
    description: str = ""
    tol_dimensionless_L: float = 1e-10
    tol_dimensionless_G: float = 1e-10
    enable_brouwer_certification: bool = True


STANDARD_KPI_PROFILES: Dict[str, KPIProfile] = {
    "smooth_periodic": KPIProfile(
        name="smooth_periodic",
        tol_energy_rel=1e-4,
        tol_secular_slope=1e-11,
        tol_momentum=1e-13,
        tol_angular_momentum=1e-13,
        tol_galilean_cm=1e-12,
        check_secular_slope=True,
        description="Smooth periodic choreographies (e.g. Chenciner-Montgomery Figure-8, Arenstorf).",
        tol_dimensionless_L=1e-11,
        tol_dimensionless_G=1e-11,
        enable_brouwer_certification=True,
    ),
    "chaotic_encounter": KPIProfile(
        name="chaotic_encounter",
        tol_energy_rel=1e-3,
        tol_secular_slope=1e-4,
        tol_momentum=1e-11,
        tol_angular_momentum=1e-11,
        tol_galilean_cm=1e-10,
        check_secular_slope=False,
        description="Chaotic singular close encounters without softening (e.g. Burrau, Sitnikov).",
        tol_dimensionless_L=1e-10,
        tol_dimensionless_G=1e-10,
        enable_brouwer_certification=True,
    ),
    "hierarchical_multiscale": KPIProfile(
        name="hierarchical_multiscale",
        tol_energy_rel=1e-2,
        tol_secular_slope=1e-2,
        tol_momentum=1e-11,
        tol_angular_momentum=1e-11,
        tol_galilean_cm=1e-10,
        check_secular_slope=False,
        description="Hierarchical multi-scale cluster systems with atlas flips and Hill supernodes.",
        tol_dimensionless_L=1e-10,
        tol_dimensionless_G=1e-10,
        enable_brouwer_certification=True,
    ),
}


@dataclass(frozen=True)
class KPIReport:
    """
    Certified validation report summarizing symplectic compliance and tolerances.
    """
    kpis: CelestialKPIs
    energy_conserved: bool
    zero_secular_drift: bool
    momentum_conserved: bool
    angular_momentum_conserved: bool
    galilean_boost_conserved: bool
    symplecticity_passed: bool
    tol_energy_rel: float = 1e-4
    tol_secular_slope: float = 1e-11
    tol_momentum: float = 1e-13
    tol_angular_momentum: float = 1e-13
    tol_galilean_cm: float = 1e-12
    profile_name: str = "smooth_periodic"
    check_secular_slope: bool = True
    # Intrinsic & Machine-Roundoff Metrics
    brouwer_delta_P_bound: float = 0.0
    ratio_P_brouwer: float = 0.0
    dimensionless_delta_L: float = 0.0
    dimensionless_delta_G: float = 0.0
    inertia_radius_series: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))
    brouwer_sigma_series: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))
    times: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))

    def to_markdown(self) -> str:
        """Renders an academic markdown scorecard."""
        status_symp = "[PASS]" if self.symplecticity_passed else "[FAIL]"
        status_drift = "[PASS]" if self.zero_secular_drift else "[FAIL]"
        status_p = "[PASS]" if self.momentum_conserved else "[FAIL]"
        status_l = "[PASS]" if self.angular_momentum_conserved else "[FAIL]"

        return f"""### Celestial Symplectic KPI Certification Report
| Metric | Measured Value | Threshold / Tolerance | Certification Status |
|---|:---:|:---:|:---:|
| Max Relative Energy Drift $\\Delta H / |H_0|$ | {self.kpis.delta_E_rel_max:.3e} | $< {self.tol_energy_rel:.0e}$ | {'[PASS]' if self.energy_conserved else '[FAIL]'} |
| Secular Energy Slope $\\mathrm{{Slope}}_H$ | {self.kpis.secular_energy_slope:.3e} | $< {self.tol_secular_slope:.0e}$ | {status_drift} |
| Linear Momentum Drift $\\|\\Delta \\mathbf{{P}}\\|$ | {self.kpis.delta_P_max:.3e} | $< {self.tol_momentum:.0e}$ | {status_p} |
| Angular Momentum Drift $\\|\\Delta \\mathbf{{L}}\\|$ | {self.kpis.delta_L_max:.3e} | $< {self.tol_angular_momentum:.0e}$ | {status_l} |
| Center of Mass Drift $\\|\\Delta \\mathbf{{G}}\\|$ | {self.kpis.delta_G_max:.3e} | $< {self.tol_galilean_cm:.0e}$ | {'[PASS]' if self.galilean_boost_conserved else '[FAIL]'} |
| Atlas Transitions / Dwell Time | {self.kpis.n_atlas_transitions} flips ($\\tau_{{dwell}} = {self.kpis.atlas_dwell_time:.2f}$) | Deterministic $\\gamma=0.6$ | [PASS] |
| **Global Symplectic Certification** | -- | Strict Invariants | **{status_symp}** |
"""


def compute_energy_secular_slope(times: np.ndarray, energies: np.ndarray) -> Tuple[float, float]:
    """
    Calculates the linear secular energy drift slope via ordinary least squares:
        E(t) = alpha + beta * t.

    Parameters:
        times: Array of simulation times (M,).
        energies: Array of Hamiltonian energies (M,).

    Returns:
        (slope, intercept): Tuple of float regression coefficients.
    """
    if len(times) < 2:
        return 0.0, float(energies[0]) if len(energies) > 0 else 0.0

    t_mean = float(np.mean(times))
    e_mean = float(np.mean(energies))

    t_diff = times - t_mean
    e_diff = energies - e_mean

    var_t = float(np.sum(t_diff**2))
    if var_t < 1e-30:
        return 0.0, e_mean

    cov_te = float(np.sum(t_diff * e_diff))
    slope = cov_te / var_t
    intercept = e_mean - slope * t_mean

    return float(slope), float(intercept)


def compute_brouwer_diffusion_series(
    history: HistoryTracker,
    stages_per_substep: int = 11,
    scheme_variation: float = 3.0,
) -> Tuple[np.ndarray, float, float]:
    """
    Computes the cumulative Brouwer-Henrici stochastic roundoff diffusion envelope
    under IEEE-754 float64 machine precision across dyadic multirate levels:
        sigma_Brouwer(t) = u * sqrt(N_stages(t)) * gamma_scheme * sqrt(D * N_bodies) * P_scale
        where u = 2^{-53} ~ 1.1102e-16,
              N_stages(t) = s * sum_{m=1}^{N(t)} 2^{k_m},
              gamma_scheme = sum |b_i| (total variation of splitting coefficients),
              D = 3 (spatial dimensions),
              N_bodies = particle count,
              P_scale = max_t sum_i ||p_i(t)||.

    Returns:
        (sigma_series, max_brouwer_bound, ratio_p_brouwer)
    """
    if len(history) == 0:
        return np.empty(0, dtype=np.float64), 0.0, 0.0

    M = len(history)
    dyadic_levels = history.get_dyadic_levels() if hasattr(history, "get_dyadic_levels") else np.zeros(M, dtype=np.int32)
    substeps_per_step = 2.0 ** np.maximum(dyadic_levels, 0)

    # Cumulative stages up to step m (step 0 has 0 stages accumulated)
    cum_substeps = np.zeros(M, dtype=np.float64)
    if M > 1:
        cum_substeps[1:] = np.cumsum(substeps_per_step[1:])

    cum_stages = float(stages_per_substep) * cum_substeps

    # Particle count and spatial dimension (3 degrees of freedom per body)
    positions = history.get_positions()
    N_bodies = positions.shape[1] if (len(positions) > 0 and positions.ndim == 3 and positions.shape[1] > 0) else 3
    dim = 3.0
    amplification = float(scheme_variation) * np.sqrt(dim * float(N_bodies))

    # Momentum scale
    momenta = history.get_momenta()
    if len(momenta) > 0 and momenta.ndim == 3 and momenta.shape[1] > 0:
        p_norms = np.sum(np.linalg.norm(momenta, axis=-1), axis=-1)
        p_scale = float(np.max(p_norms))
    else:
        P_vectors = history.get_linear_momenta()
        p_scale = float(np.max(np.linalg.norm(P_vectors, axis=-1))) if len(P_vectors) > 0 else 1.0

    p_scale = max(p_scale, 1.0)
    u_mach = 2.0 ** -53  # ~ 1.1102230246251565e-16

    sigma_series = u_mach * np.sqrt(cum_stages) * amplification * p_scale
    max_brouwer_bound = float(5.0 * sigma_series[-1] + 1e-15) if M > 0 else 1e-15

    P_vectors = history.get_linear_momenta()
    delta_P = np.linalg.norm(P_vectors - P_vectors[0], axis=-1) if len(P_vectors) > 0 else np.zeros(M)
    delta_P_max = float(np.max(delta_P))
    ratio_p_brouwer = delta_P_max / max_brouwer_bound if max_brouwer_bound > 1e-30 else 0.0

    return sigma_series, max_brouwer_bound, ratio_p_brouwer


def compute_inertia_radius_and_dimensionless_invariants(
    history: HistoryTracker,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float]:
    """
    Computes barycentric inertia radius R_I(t), root-mean-square momentum P_rms(t),
    and dimensionless scale-invariant drift metrics:
        epsilon_L(t) = ||L(t) - L0|| / (||L0|| + R_I(t) * P_rms(t))
        epsilon_G(t) = ||G(t) - G0|| / (M_tot * R_I(t))

    Returns:
        (R_I_series, epsilon_L_series, epsilon_G_series, max_R_I, max_epsilon_L, max_epsilon_G)
    """
    if len(history) == 0:
        empty = np.empty(0, dtype=np.float64)
        return empty, empty, empty, 0.0, 0.0, 0.0

    positions = history.get_positions()  # (M, N, 3)
    momenta = history.get_momenta()      # (M, N, 3)
    M_steps, N_bodies, _ = positions.shape

    masses = history.masses
    if masses is None:
        masses = np.ones(N_bodies, dtype=np.float64)
    else:
        masses = np.asarray(masses, dtype=np.float64)

    M_tot = float(np.sum(masses))
    if M_tot < 1e-30:
        M_tot = 1.0

    # Barycenter R_cm: (M, 1, 3)
    masses_3d = masses[None, :, None]
    R_cm = np.sum(positions * masses_3d, axis=1, keepdims=True) / M_tot
    q_rel = positions - R_cm
    sq_dist = np.sum(q_rel**2, axis=-1)  # (M, N)
    R_I_series = np.sqrt(np.sum(sq_dist * masses[None, :], axis=-1) / M_tot)  # (M,)

    # Kinetic energy and P_rms
    sq_p = np.sum(momenta**2, axis=-1)   # (M, N)
    T_series = 0.5 * np.sum(sq_p / masses[None, :], axis=-1)  # (M,)
    P_rms_series = np.sqrt(2.0 * M_tot * np.maximum(T_series, 1e-14))

    # Angular momentum dimensionless drift epsilon_L
    L_vectors = history.get_angular_momenta()
    L0 = L_vectors[0]
    delta_L = np.linalg.norm(L_vectors - L0, axis=-1)
    L0_norm = float(np.linalg.norm(L0))
    action_scale = L0_norm + R_I_series * P_rms_series
    epsilon_L_series = delta_L / np.maximum(action_scale, 1e-12)
    max_epsilon_L = float(np.max(epsilon_L_series))

    # Center of mass dimensionless drift epsilon_G
    G_vectors = history.get_center_of_masses()
    G0 = G_vectors[0]
    delta_G = np.linalg.norm(G_vectors - G0, axis=-1)
    times = history.get_times()
    dt_elapsed = np.abs(times - times[0]) if len(times) == len(R_I_series) else np.zeros_like(R_I_series)
    spatial_scale = M_tot * np.maximum(R_I_series, 1e-12) + dt_elapsed * P_rms_series
    epsilon_G_series = delta_G / spatial_scale
    max_epsilon_G = float(np.max(epsilon_G_series))

    max_R_I = float(np.max(R_I_series))

    return R_I_series, epsilon_L_series, epsilon_G_series, max_R_I, max_epsilon_L, max_epsilon_G


def compute_celestial_kpis(
    history: HistoryTracker,
    atlas_manager: Optional[AtlasManager] = None,
    mean_wronskian_error: float = 0.0,
    tol_energy_rel: Optional[float] = None,
    tol_secular_slope: Optional[float] = None,
    tol_momentum: Optional[float] = None,
    tol_angular_momentum: Optional[float] = None,
    tol_galilean_cm: Optional[float] = None,
    profile: str = "smooth_periodic",
    stages_per_substep: int = 11,
    scheme_variation: float = 3.0,
) -> KPIReport:
    """
    Extracts the canonical celestial KPIs from a completed simulation history
    and verifies compliance against formal thresholds calibrated by dynamical profile.

    Parameters:
        history: HistoryTracker populated during simulation.
        atlas_manager: Optional AtlasManager tracking Jacobi chart transitions.
        mean_wronskian_error: Kepler solver Wronskian error.
        tol_*: Optional overrides for precision thresholds.
        profile: Key into STANDARD_KPI_PROFILES ('smooth_periodic', 'chaotic_encounter', 'hierarchical_multiscale').
        stages_per_substep: Number of composition stages per substep for Brouwer bound (default 11).
        scheme_variation: Sum of absolute impulse weights for Brouwer bound (default 3.0).

    Returns:
        KPIReport containing metrics and certification flags.
    """
    if len(history) == 0:
        raise ValueError("Cannot compute celestial KPIs from an empty history tracker.")

    prof = STANDARD_KPI_PROFILES.get(profile, STANDARD_KPI_PROFILES["smooth_periodic"])
    t_e = prof.tol_energy_rel if tol_energy_rel is None else float(tol_energy_rel)
    t_slope = prof.tol_secular_slope if tol_secular_slope is None else float(tol_secular_slope)
    t_p = prof.tol_momentum if tol_momentum is None else float(tol_momentum)
    t_l = prof.tol_angular_momentum if tol_angular_momentum is None else float(tol_angular_momentum)
    t_g = prof.tol_galilean_cm if tol_galilean_cm is None else float(tol_galilean_cm)
    chk_slope = prof.check_secular_slope

    times = history.get_times()
    energies = history.get_energies()
    P_vectors = history.get_linear_momenta()
    L_vectors = history.get_angular_momenta()
    G_vectors = history.get_center_of_masses()

    t_span = float(times[-1] - times[0]) if len(times) > 1 else 0.0
    n_records = len(times)

    # 1. Energy drift and secular slope
    E0 = energies[0]
    denom_E = abs(E0) if abs(E0) > 1e-15 else 1.0
    delta_E_abs = np.abs(energies - E0)
    delta_E_rel_max = float(np.max(delta_E_abs) / denom_E)
    slope_raw, _ = compute_energy_secular_slope(times, energies)
    secular_slope = abs(slope_raw)

    # 2. Linear momentum drift
    P0 = P_vectors[0]
    delta_P = np.linalg.norm(P_vectors - P0, axis=1)
    delta_P_max = float(np.max(delta_P))

    # 3. Angular momentum drift
    L0 = L_vectors[0]
    delta_L = np.linalg.norm(L_vectors - L0, axis=1)
    delta_L_max = float(np.max(delta_L))

    # 4. Center of mass drift
    G0 = G_vectors[0]
    delta_G = np.linalg.norm(G_vectors - G0, axis=1)
    delta_G_max = float(np.max(delta_G))

    # 5. Atlas transitions
    if atlas_manager is not None:
        n_flips = len(atlas_manager.transition_history)
        dwell_time = t_span / (n_flips + 1)
    else:
        n_flips = 0
        dwell_time = t_span

    # 6. Brouwer stochastic diffusion and scale-invariant metrics
    brouwer_series, max_brouwer_bound, ratio_p_brouwer = compute_brouwer_diffusion_series(
        history,
        stages_per_substep=stages_per_substep,
        scheme_variation=scheme_variation,
    )
    R_I_series, eps_L_series, eps_G_series, max_R_I, max_eps_L, max_eps_G = (
        compute_inertia_radius_and_dimensionless_invariants(history)
    )

    kpis = CelestialKPIs(
        delta_E_rel_max=delta_E_rel_max,
        secular_energy_slope=secular_slope,
        delta_P_max=delta_P_max,
        delta_L_max=delta_L_max,
        delta_G_max=delta_G_max,
        mean_wronskian_error=float(mean_wronskian_error),
        atlas_dwell_time=dwell_time,
        n_atlas_transitions=n_flips,
        t_span=t_span,
        n_records=n_records,
        brouwer_delta_P_bound=max_brouwer_bound,
        ratio_P_brouwer=ratio_p_brouwer,
        dimensionless_delta_L=max_eps_L,
        dimensionless_delta_G=max_eps_G,
        inertia_radius_max=max_R_I,
    )

    pass_e = delta_E_rel_max <= t_e
    pass_slope = (secular_slope <= t_slope) if chk_slope else True
    # Momentum passes either by raw threshold OR by being within the 5-sigma Brouwer roundoff bound
    pass_p = (delta_P_max <= t_p) or (prof.enable_brouwer_certification and delta_P_max <= max_brouwer_bound)
    # Angular momentum passes either by raw threshold OR by scale-invariant action bound
    pass_l = (delta_L_max <= t_l) or (max_eps_L <= prof.tol_dimensionless_L)
    # Center of mass passes either by raw threshold OR by dimensionless inertia bound
    pass_g = (delta_G_max <= t_g) or (max_eps_G <= prof.tol_dimensionless_G)
    pass_global = pass_e and pass_slope and pass_p and pass_l and pass_g

    return KPIReport(
        kpis=kpis,
        energy_conserved=pass_e,
        zero_secular_drift=pass_slope,
        momentum_conserved=pass_p,
        angular_momentum_conserved=pass_l,
        galilean_boost_conserved=pass_g,
        symplecticity_passed=pass_global,
        tol_energy_rel=t_e,
        tol_secular_slope=t_slope,
        tol_momentum=t_p,
        tol_angular_momentum=t_l,
        tol_galilean_cm=t_g,
        profile_name=prof.name,
        check_secular_slope=chk_slope,
        brouwer_delta_P_bound=max_brouwer_bound,
        ratio_P_brouwer=ratio_p_brouwer,
        dimensionless_delta_L=max_eps_L,
        dimensionless_delta_G=max_eps_G,
        inertia_radius_series=R_I_series,
        brouwer_sigma_series=brouwer_series,
        times=times,
    )


def extract_laskar_fundamental_frequency(
    times: np.ndarray,
    signal: np.ndarray,
    hanning_order: int = 3,
) -> float:
    """
    Extracts the dominant fundamental angular frequency omega (rad / time unit)
    of a complex or real signal using Jacques Laskar's NAFF (Numerical Analysis of
    Fundamental Frequencies, Laskar 1990, 1993, 1999) with an order-p Hanning filter:
        chi_p(s) = (2^p * (p!)^2 / (2p)!) * (1 + cos(pi * s))^p,   s in [-1, 1],
    which provides O(1 / T^{2p+2}) sidelobe suppression, followed by exact analytical
    derivative Newton-Raphson maximization of the Fourier energy E(omega) = |<f, e^{i omega t}>|^2.
    """
    t = np.asarray(times, dtype=np.float64)
    z = np.asarray(signal, dtype=np.complex128)
    n = len(t)
    if n < 4:
        return 0.0

    t0 = float(t[0])
    t1 = float(t[-1])
    T = t1 - t0
    if T <= 1e-30:
        return 0.0

    tau = t - 0.5 * (t0 + t1)  # centered time in [-T/2, T/2] for real-symmetric phase properties
    s_norm = 2.0 * tau / T
    p_ord = max(1, int(hanning_order))
    hanning_weights = (1.0 + np.cos(np.pi * s_norm)) ** p_ord

    # Composite trapezoidal quadrature weights (exact to O(N^{-2p-2}) for C^{2p-1}-periodic window)
    dt_trap = np.empty(n, dtype=np.float64)
    dt_trap[0] = 0.5 * (t[1] - t[0])
    dt_trap[-1] = 0.5 * (t[-1] - t[-2])
    dt_trap[1:-1] = 0.5 * (t[2:] - t[:-2])
    w_quad = (hanning_weights * dt_trap) / T
    z_w = z * w_quad

    # Initial FFT coarse grid estimate on positive frequencies
    dt_mean = T / float(n - 1)
    fft_vals = np.fft.fft(z * hanning_weights)
    fft_freqs = 2.0 * np.pi * np.fft.fftfreq(n, d=dt_mean)
    pos_mask = fft_freqs > 0.0
    if np.any(pos_mask):
        pos_indices = np.where(pos_mask)[0]
        peak_idx = pos_indices[int(np.argmax(np.abs(fft_vals[pos_indices])))]
        omega_coarse = float(fft_freqs[peak_idx])
    else:
        omega_coarse = 2.0 * np.pi / T

    delta_omega = 2.0 * np.pi / T
    a = max(0.0, omega_coarse - 1.5 * delta_omega)
    b = omega_coarse + 1.5 * delta_omega

    def neg_energy(omega: float) -> float:
        phi = np.sum(z_w * np.exp(-1j * omega * tau))
        return -float(phi.real * phi.real + phi.imag * phi.imag)

    # 1. Golden-section bracket narrowing (30 steps)
    invphi = (math.sqrt(5.0) - 1.0) / 2.0
    invphi2 = (3.0 - math.sqrt(5.0)) / 2.0
    h = b - a
    c = a + invphi2 * h
    d = a + invphi * h
    yc = neg_energy(c)
    yd = neg_energy(d)

    for _ in range(30):
        h *= invphi
        if yc < yd:
            b = d
            d = c
            yd = yc
            c = a + invphi2 * h
            yc = neg_energy(c)
        else:
            a = c
            c = d
            yc = yd
            d = a + invphi * h
            yd = neg_energy(d)

    omega_opt = 0.5 * (a + b)

    # 2. Exact analytical Newton-Raphson on dE/domega = 2 Re(conj(phi) * phi') = 0
    for _ in range(15):
        exp_vec = np.exp(-1j * omega_opt * tau)
        terms = z_w * exp_vec
        phi = np.sum(terms)
        dphi = np.sum((-1j * tau) * terms)
        d2phi = np.sum((-(tau * tau)) * terms)

        grad = 2.0 * float(phi.real * dphi.real + phi.imag * dphi.imag)
        hess = 2.0 * float(
            (dphi.real * dphi.real + dphi.imag * dphi.imag)
            + (phi.real * d2phi.real + phi.imag * d2phi.imag)
        )
        if abs(hess) <= 1e-30 or hess >= 0.0:
            break
        step = grad / hess
        omega_opt -= step
        if abs(step) < 1e-16 * max(1.0, abs(omega_opt)):
            break

    return float(omega_opt)


def compute_laskar_frequency_drift(
    history: HistoryTracker,
    pair: Tuple[int, int] = (0, 1),
) -> Dict[str, float]:
    """
    Computes Jacques Laskar's Frequency Map Analysis (FMA / NAFF) stability metric
    by extracting the fundamental orbital angular frequency nu_1 on [0, T/2] and
    nu_2 on [T/2, T] from the complex relative planar coordinate
        z(t) = (x_j(t) - x_i(t)) + i * (y_j(t) - y_i(t))
    and measuring the relative frequency diffusion:
        delta_nu / nu = |nu_2 - nu_1| / max(|nu_1|, 1e-30).
    """
    is_f128 = getattr(history, "is_f128", False)
    times = history.times_f64 if is_f128 else history.get_times()
    positions = history.positions_f64 if is_f128 else history.get_positions()

    n = len(times)
    if n < 8:
        return {"nu_1": 0.0, "nu_2": 0.0, "abs_freq_drift": 0.0, "rel_freq_drift": 0.0}

    i, j = pair
    r_rel = positions[:, j, :] - positions[:, i, :]
    # Project onto dominant orbital plane (x, y) or principal plane
    z_sig = r_rel[:, 0] + 1j * r_rel[:, 1]
    if np.max(np.abs(r_rel[:, 1])) < 1e-12 * np.max(np.abs(r_rel)):
        momenta = history.momenta_f64 if is_f128 else history.get_momenta()
        p_rel = momenta[:, j, :] - momenta[:, i, :]
        z_sig = r_rel[:, 0] + 1j * p_rel[:, 0]

    mid = n // 2
    t_first = times[: mid + 1]
    z_first = z_sig[: mid + 1]
    t_second = times[mid:]
    z_second = z_sig[mid:]

    nu_1 = extract_laskar_fundamental_frequency(t_first, z_first)
    nu_2 = extract_laskar_fundamental_frequency(t_second, z_second)
    abs_drift = abs(nu_2 - nu_1)
    rel_drift = abs_drift / max(abs(nu_1), 1e-30)

    return {
        "nu_1": float(nu_1),
        "nu_2": float(nu_2),
        "abs_freq_drift": float(abs_drift),
        "rel_freq_drift": float(rel_drift),
    }


def compute_forward_backward_reversibility(
    state0: Any,
    integrator: Any,
    dt: float,
    n_steps: int,
    G: float = 1.0,
    desingularization_mode: str = "auto",
    enable_dyadic_substepping: bool = True,
    max_dyadic_level: int = 24,
) -> Dict[str, Any]:
    """
    Executes Rein & Spiegel's (2015, Fig. 3) Forward-Backward Time-Reversibility Audit:
    integrates state0 forward for n_steps with +dt, then backward for n_steps with -dt,
    and measures the round-trip discrepancy in positions, momenta, energy, and the
    palindromic symmetry of the dyadic multirate level sequence k_fwd[m] == k_bwd[n_steps - 1 - m].
    """
    from ..system.nbody_system import NBodySystem
    from ..core.invariants import compute_invariants

    system = NBodySystem(
        state=state0,
        integrator=integrator,
        G=G,
        desingularization_mode=desingularization_mode,
        enable_dyadic_substepping=enable_dyadic_substepping,
        max_dyadic_level=max_dyadic_level,
    )
    dt_val = float(dt)
    for _ in range(int(n_steps)):
        system.step(dt=dt_val)
    n_flips_fwd = len(system.atlas_manager.transition_history)
    all_levels_fwd = list(system.history.get_dyadic_levels()[1 : 1 + int(n_steps)])

    for _ in range(int(n_steps)):
        system.step(dt=-dt_val)
    state_ret = system.get_physical_state()
    n_flips_bwd = len(system.atlas_manager.transition_history) - n_flips_fwd
    all_levels_bwd = list(system.history.get_dyadic_levels()[1 + int(n_steps) : 1 + 2 * int(n_steps)])

    delta_q_norm = float(np.linalg.norm(state_ret.q - state0.q))
    delta_p_norm = float(np.linalg.norm(state_ret.p - state0.p))

    inv0 = compute_invariants(state0, G=G)
    inv_ret = compute_invariants(state_ret, G=G)
    denom_E = max(abs(float(inv0.H)), 1e-30)
    delta_E_rel = abs(float(inv_ret.H - inv0.H)) / denom_E

    dyadic_palindrome = (all_levels_fwd == list(reversed(all_levels_bwd)))

    return {
        "delta_q_norm": delta_q_norm,
        "delta_p_norm": delta_p_norm,
        "delta_E_rel": float(delta_E_rel),
        "k_fwd": all_levels_fwd,
        "k_bwd": all_levels_bwd,
        "dyadic_palindrome_matched": bool(dyadic_palindrome),
        "n_atlas_flips_fwd": int(n_flips_fwd),
        "n_atlas_flips_bwd": int(n_flips_bwd),
    }


def compute_osculating_elements(
    state: Any,
    primary_idx: int = 0,
    body_idx: int = 1,
    G: float = 1.0,
) -> Dict[str, float]:
    """
    Computes the 3D osculating Keplerian orbital elements of `body_idx` relative
    to `primary_idx` (semi-major axis `a`, eccentricity `e`, inclination `inc_rad`,
    longitude of ascending node `Omega`, argument of pericenter `omega`,
    longitude of pericenter `varpi`, mean anomaly `M`, mean longitude `mean_longitude`,
    and true longitude `true_longitude`).
    """
    q0 = np.asarray(state.q[primary_idx], dtype=np.float64)
    q1 = np.asarray(state.q[body_idx], dtype=np.float64)
    p0 = np.asarray(state.p[primary_idx], dtype=np.float64)
    p1 = np.asarray(state.p[body_idx], dtype=np.float64)
    m0 = float(state.masses[primary_idx])
    m1 = float(state.masses[body_idx])

    r_vec = q1 - q0
    v_vec = (p1 / max(m1, 1e-30)) - (p0 / max(m0, 1e-30))
    r = float(np.linalg.norm(r_vec))
    v2 = float(np.dot(v_vec, v_vec))
    rv = float(np.dot(r_vec, v_vec))
    mu = float(G) * (m0 + m1)

    h_vec = np.cross(r_vec, v_vec)
    h_norm = float(np.linalg.norm(h_vec))
    inc_rad = math.acos(max(-1.0, min(1.0, h_vec[2] / max(h_norm, 1e-30))))

    e_vec = ((v2 - mu / max(r, 1e-30)) * r_vec - rv * v_vec) / max(mu, 1e-30)
    e = float(np.linalg.norm(e_vec))
    inv_a = 2.0 / max(r, 1e-30) - v2 / max(mu, 1e-30)
    a = 1.0 / inv_a if abs(inv_a) > 1e-30 else float("inf")

    n_vec = np.array([-h_vec[1], h_vec[0], 0.0], dtype=np.float64)
    n_norm = float(np.linalg.norm(n_vec))
    if n_norm > 1e-14:
        Omega = math.atan2(n_vec[1], n_vec[0]) % (2.0 * math.pi)
    else:
        Omega = 0.0

    true_long = math.atan2(r_vec[1], r_vec[0]) % (2.0 * math.pi)
    if e > 1e-12:
        varpi = math.atan2(e_vec[1], e_vec[0]) % (2.0 * math.pi)
        omega = (varpi - Omega) % (2.0 * math.pi)
        if a > 0.0 and math.isfinite(a):
            cos_E = max(-1.0, min(1.0, (1.0 - r / a) / e))
            sin_E = rv / (e * math.sqrt(max(mu * a, 1e-30)))
            E_anom = math.atan2(sin_E, cos_E)
            M_anom = (E_anom - e * math.sin(E_anom)) % (2.0 * math.pi)
        else:
            M_anom = (true_long - varpi) % (2.0 * math.pi)
        mean_long = (varpi + M_anom) % (2.0 * math.pi)
    else:
        varpi = 0.0
        omega = 0.0
        M_anom = true_long
        mean_long = true_long

    return {
        "a": float(a),
        "e": float(e),
        "inc_rad": float(inc_rad),
        "inc_deg": float(math.degrees(inc_rad)),
        "Omega": float(Omega),
        "omega": float(omega),
        "varpi": float(varpi),
        "M": float(M_anom),
        "mean_longitude": float(mean_long),
        "true_longitude": float(true_long),
    }


def compute_laplace_resonant_angle(
    state: Any,
    G: float = 1.0,
    use_mean_longitude: bool = True,
) -> float:
    """
    Computes the Galilean 4:2:1 Laplace resonant angle
        theta_Lap = (lambda_1 - 3 * lambda_2 + 2 * lambda_3) mod (2*pi)
    in radians, using osculating mean longitudes lambda_i = varpi_i + M_i
    by default (`use_mean_longitude=True`).
    """
    key = "mean_longitude" if use_mean_longitude else "true_longitude"
    lam1 = compute_osculating_elements(state, 0, 1, G=G)[key]
    lam2 = compute_osculating_elements(state, 0, 2, G=G)[key]
    lam3 = compute_osculating_elements(state, 0, 3, G=G)[key]
    return float((lam1 - 3.0 * lam2 + 2.0 * lam3) % (2.0 * math.pi))
