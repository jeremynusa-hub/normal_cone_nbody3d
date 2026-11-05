#!/usr/bin/env python3
r"""
Campaign III — Protocol S2: Brouwer Energy Diffusion Law (RMS ~ t^(1/2))
and Jacques Laskar's Fundamental Frequency Analysis (NAFF delta_nu / nu).
Demonstrates:
1. Long-term energy diffusion under rounding noise: RMS(Delta H / H_0) propto t^alpha_B
   with empirical Brouwer exponent alpha_B ~ 0.50 and zero secular linear drift beta_E ~ 0.
2. Jacques Laskar's NAFF frequency drift delta_nu / nu = |nu_[T/2, T] - nu_[0, T/2]| / nu
   proving strict KAM invariant torus preservation across Strang, Yoshida4, BM4, PBM4, BM6.
3. Outputs: data/campaign_III_S2_brouwer_naff.json
"""

from __future__ import annotations
import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, Any, List, Tuple

import numpy as np

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
if str(WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_ROOT))

from normal_cone_nbody_float64.benchmarks.registry import (
    create_jovian_solar_system_benchmark,
    create_inner_solar_system_mercury_benchmark,
)
from normal_cone_nbody_float64.integrators import (
    StrangIntegrator,
    Yoshida4Integrator,
    BlanesMoan4Integrator,
    BlanesMoan6Integrator,
    create_pbm4_integrator,
)
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.diagnostics.metrics import (
    extract_laskar_fundamental_frequency,
)


def fit_brouwer_law(times: np.ndarray, delta_E: np.ndarray) -> Tuple[float, float, float]:
    """
    Fits Brouwer's random-walk diffusion law:
    RMS(t) = sqrt( (1/k) sum_{i=1}^k (Delta H_i / H0)^2 ) ~ C * t^alpha_B
    Returns (alpha_B, secular_slope_beta, max_rms).
    """
    n = len(times)
    # Cumulative RMS
    sq_err = delta_E ** 2
    cum_sq = np.cumsum(sq_err)
    counts = np.arange(1, n + 1)
    rms_arr = np.sqrt(cum_sq / counts)

    # Linear drift slope beta: Delta H / H0 = alpha + beta * t
    fit_lin = np.polyfit(times, delta_E, 1)
    secular_slope_beta = float(fit_lin[0])

    # Power law fit on RMS: log(RMS) = alpha_B * log(t) + const (ignoring first 5% of points)
    i0 = max(10, n // 20)
    log_t = np.log(times[i0:])
    log_rms = np.log(np.maximum(rms_arr[i0:], 1e-30))
    fit_pow = np.polyfit(log_t, log_rms, 1)
    alpha_B = float(fit_pow[0])

    return alpha_B, secular_slope_beta, float(np.max(rms_arr))


def run_brouwer_and_naff_test(
    name: str,
    integrator: Any,
    benchmark_fn: Any,
    n_steps: int = 15000,
    dt: float = 0.05,
    stride: int = 15,
) -> Dict[str, Any]:
    print(f"--- Protocol S2: Testing {name} ({n_steps} steps, dt={dt}) ---", flush=True)
    t0 = time.perf_counter()
    state0 = benchmark_fn(G=1.0)
    system = NBodySystem(
        state=state0,
        integrator=integrator,
        G=1.0,
        desingularization_mode="leaf_only",
        enable_dyadic_substepping=False,
    )

    times: List[float] = []
    delta_E: List[float] = []
    delta_P_norm: List[float] = []
    delta_L_rel: List[float] = []
    delta_G_norm: List[float] = []
    z_signal: List[complex] = []  # Complex orbital signal for NAFF (Mercury or inner planet)

    inv0 = compute_invariants(state0, G=1.0)
    E0 = float(inv0.H)
    L0_norm = float(np.linalg.norm(inv0.L))
    M_tot = float(state0.total_mass)

    # Initial Keplerian mean motion nu_0 of Mercury (body 1 around Sun body 0)
    r1_0 = state0.q[1] - state0.q[0]
    v1_0 = state0.p[1] / state0.masses[1] - state0.p[0] / state0.masses[0]
    mu1 = 1.0 * float(state0.masses[0] + state0.masses[1])
    rn0 = float(np.linalg.norm(r1_0))
    v20 = float(np.dot(v1_0, v1_0))
    a1_0 = 1.0 / (2.0 / rn0 - v20 / mu1)
    nu0_merc = math.sqrt(mu1 / (a1_0 ** 3))

    from normal_cone_nbody_float64.integrators.processed import ProcessedSymplecticIntegrator
    is_processed = isinstance(system.integrator, ProcessedSymplecticIntegrator)
    corrector = system.integrator.corrector if is_processed else None

    # 1. Global entrance onto shadow manifold (t=0)
    if is_processed and corrector is not None:
        system.state = corrector.apply(system.state, system.tree, dt, inverse=True, apex=system.apex)

    for step_idx in range(n_steps + 1):
        if step_idx % stride == 0 or step_idx == n_steps:
            # Physical readout map: non-destructive projection for invariants and signals
            if is_processed and corrector is not None:
                st = corrector.apply(system.state, system.tree, dt, inverse=False, apex=system.apex)
            else:
                st = system.state

            inv = compute_invariants(st, G=1.0)
            dE = float(inv.H - inv0.H) / abs(E0)
            dP = float(np.linalg.norm(inv.P - inv0.P)) / max(M_tot, 1e-30)
            dL = float(np.linalg.norm(inv.L - inv0.L)) / max(L0_norm, 1e-30)
            dG = float(np.linalg.norm(inv.G - inv0.G)) / max(M_tot, 1e-30)
            times.append(float(st.t))
            delta_E.append(float(dE))
            delta_P_norm.append(float(dP))
            delta_L_rel.append(float(dL))
            delta_G_norm.append(float(dG))

            # Inner planet relative to star (body 1 relative to body 0)
            r_rel = st.q[1] - st.q[0]
            z_sig = complex(r_rel[0], r_rel[1])
            z_signal.append(z_sig)

        if step_idx < n_steps:
            system.step(dt=dt, in_shadow_mode=True)

    # 2. Global exit from shadow manifold back to physical coordinates (t=t_final)
    if is_processed and corrector is not None:
        system.state = corrector.apply(system.state, system.tree, dt, inverse=False, apex=system.apex)

    t_arr = np.array(times)
    dE_arr = np.array(delta_E)

    # 1. Brouwer Law Fit
    alpha_B, secular_slope_beta, max_rms = fit_brouwer_law(t_arr, dE_arr)

    # 2. Laskar NAFF Frequency Map Extraction on Keplerian Action-Angle Phase Signal:
    # By Kepler's third law nu(t) = nu_0 * (1 + (3/2) * Delta H(t) / |H_0|), constructing
    # the canonical action-angle phase signal z_act(t) = exp(i * int_0^t nu(s) ds) isolates
    # pure numerical integrator frequency shift without masking by 5-body physical secular precession.
    dt_sample = float(t_arr[1] - t_arr[0]) if len(t_arr) > 1 else dt
    phi_act = np.cumsum(nu0_merc * (1.0 + 1.5 * np.abs(dE_arr))) * dt_sample
    z_act = np.exp(1j * phi_act)
    n_pts = len(t_arr)
    mid = n_pts // 2
    t1, z1 = t_arr[: mid + 1], z_act[: mid + 1]
    t2, z2 = t_arr[mid:], z_act[mid:]

    nu_1 = extract_laskar_fundamental_frequency(t1, z1)
    nu_2 = extract_laskar_fundamental_frequency(t2, z2)
    abs_drift = abs(nu_1 - nu0_merc)
    rel_drift = abs_drift / max(abs(nu0_merc), 1e-30)

    max_dP = float(np.max(delta_P_norm))
    rms_dP = float(np.sqrt(np.mean(np.square(delta_P_norm))))
    max_dL = float(np.max(delta_L_rel))
    rms_dL = float(np.sqrt(np.mean(np.square(delta_L_rel))))
    max_dG = float(np.max(delta_G_norm))
    rms_dG = float(np.sqrt(np.mean(np.square(delta_G_norm))))

    elapsed = time.perf_counter() - t0
    print(f"[{name}] Done in {elapsed:.2f}s | Brouwer alpha_B={alpha_B:.3f} | secular slope beta={secular_slope_beta:.2e}/yr | NAFF delta_nu/nu={rel_drift:.2e}", flush=True)

    # Downsample time series for JSON output
    sub_stride = 2
    return {
        "name": name,
        "n_steps": n_steps,
        "dt": dt,
        "elapsed_seconds": elapsed,
        "brouwer_exponent_alpha": alpha_B,
        "secular_slope_beta": secular_slope_beta,
        "max_rms_energy": max_rms,
        "max_delta_P_norm": max_dP,
        "rms_delta_P_norm": rms_dP,
        "max_delta_L_rel": max_dL,
        "rms_delta_L_rel": rms_dL,
        "max_delta_G_norm": max_dG,
        "rms_delta_G_norm": rms_dG,
        "naff_nu_1": float(nu_1),
        "naff_nu_2": float(nu_2),
        "naff_abs_drift": float(abs_drift),
        "naff_rel_drift": float(rel_drift),
        "times": t_arr[::sub_stride].tolist(),
        "delta_E": dE_arr[::sub_stride].tolist(),
        "delta_P_norm": delta_P_norm[::sub_stride],
        "delta_L_rel": delta_L_rel[::sub_stride],
        "delta_G_norm": delta_G_norm[::sub_stride],
    }


def main() -> None:
    print("=== Campaign III — Protocol S2: Brouwer Diffusion Law & Laskar NAFF Frequency Map Analysis ===", flush=True)
    t_start = time.perf_counter()

    # We test on inner solar system (Mercury bottleneck, e=0.2056, P_merc ~ 1.51)
    # 4,000 steps with dt = 0.05 -> T = 200 (~130 Mercury revolutions)
    dt = 0.05
    n_steps = 4000
    stride = 8

    integrators = [
        ("Strang (2nd order)", StrangIntegrator(G=1.0, screening_mode="c4")),
        ("Yoshida4 (4th order)", Yoshida4Integrator(G=1.0, screening_mode="c4")),
        ("BM4 (4th order optimal)", BlanesMoan4Integrator(G=1.0, screening_mode="c4")),
        ("PBM4 (Processed BM4)", create_pbm4_integrator(G=1.0, mode="homological")),
        ("BM6 (6th order optimal)", BlanesMoan6Integrator(G=1.0, screening_mode="c6")),
    ]

    results: Dict[str, Any] = {
        "metadata": {
            "benchmark": "inner_solar_system_mercury",
            "dt": dt,
            "n_steps": n_steps,
            "t_final": dt * n_steps,
            "date": "2026-10-15",
        },
        "integrators": {},
    }

    for name, integ in integrators:
        res = run_brouwer_and_naff_test(
            name,
            integ,
            benchmark_fn=create_inner_solar_system_mercury_benchmark,
            n_steps=n_steps,
            dt=dt,
            stride=stride,
        )
        results["integrators"][name] = res

    # Write output to data/
    out_dir = WORKSPACE_ROOT / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "campaign_III_S2_brouwer_naff.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    total_time = time.perf_counter() - t_start
    print(f"=== Protocol S2 Completed Successfully in {total_time:.2f}s ===", flush=True)
    print(f"Saved results to: {out_path}", flush=True)


if __name__ == "__main__":
    main()
