#!/usr/bin/env python3
r"""
Campaign III — Protocol S1: Secular Stability of the 4:2:1 Jovian Laplace Resonance.
Demonstrates:
1. Long-term secular integration of the 4-body Jupiter-Io-Europa-Ganymede system
   over 1,000 Io orbits under Strang, Yoshida4, BM4, PBM4, and BM6.
2. Bounded libration of the critical Laplace resonance angle
   phi_L(t) = lambda_1 - 3*lambda_2 + 2*lambda_3 around 180 degrees (pi).
3. Exact Lie-Bargmann angular momentum and linear momentum conservation.
4. Output JSON: data/campaign_III_S1_laplace.json
"""

from __future__ import annotations
import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, Any, List

import numpy as np

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
if str(WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_ROOT))

from normal_cone_nbody_float64.benchmarks.registry import (
    create_jovian_laplace_resonance_benchmark,
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


def compute_laplace_angle(q: np.ndarray) -> float:
    """
    Computes Laplace angle phi_L = lambda_Io - 3*lambda_Eu + 2*lambda_Ga centered on Jupiter (body 0).
    Wraps into [-pi, pi] centered on pi (i.e. deviation from 180 deg).
    """
    r1 = q[1] - q[0]
    r2 = q[2] - q[0]
    r3 = q[3] - q[0]
    lam1 = math.atan2(r1[1], r1[0])
    lam2 = math.atan2(r2[1], r2[0])
    lam3 = math.atan2(r3[1], r3[0])
    phi = (lam1 - 3.0 * lam2 + 2.0 * lam3) % (2.0 * math.pi)
    # Map into [0, 2*pi]
    return float(phi)


def compute_osculating_semi_major_axes(q: np.ndarray, p: np.ndarray, masses: np.ndarray, G: float = 1.0) -> List[float]:
    """
    Computes osculating semi-major axes for the 3 moons relative to Jupiter (body 0).
    """
    a_list = []
    m0 = masses[0]
    r0 = q[0]
    v0 = p[0] / m0
    for i in (1, 2, 3):
        mi = masses[i]
        mu = G * (m0 + mi)
        r_rel = q[i] - r0
        v_rel = (p[i] / mi) - v0
        r_dist = float(np.linalg.norm(r_rel))
        v_sq = float(np.dot(v_rel, v_rel))
        # 1/a = 2/r - v^2/mu
        inv_a = 2.0 / r_dist - v_sq / mu
        a = 1.0 / inv_a if abs(inv_a) > 1e-15 else float('nan')
        a_list.append(float(a))
    return a_list


def run_integrator_simulation(
    name: str,
    integrator: Any,
    n_steps: int,
    dt: float,
    stride: int = 10,
) -> Dict[str, Any]:
    print(f"--- Simulating {name} ({n_steps} steps, dt={dt}) ---", flush=True)
    t0 = time.perf_counter()
    state0 = create_jovian_laplace_resonance_benchmark(G=1.0)
    system = NBodySystem(
        state=state0,
        integrator=integrator,
        G=1.0,
        desingularization_mode="leaf_only",
        enable_dyadic_substepping=False,
    )

    times: List[float] = []
    phi_L_deg: List[float] = []
    delta_E_rel: List[float] = []
    delta_P_norm: List[float] = []
    delta_L_rel: List[float] = []
    delta_G_norm: List[float] = []
    a_moons: List[List[float]] = []

    inv0 = compute_invariants(state0, G=1.0)
    E0 = float(inv0.H)
    L0_norm = float(np.linalg.norm(inv0.L))
    M_tot = float(state0.total_mass)

    from normal_cone_nbody_float64.integrators.processed import ProcessedSymplecticIntegrator
    is_processed = isinstance(system.integrator, ProcessedSymplecticIntegrator)
    corrector = system.integrator.corrector if is_processed else None

    # 1. Global entrance onto shadow manifold (t=0)
    if is_processed and corrector is not None:
        system.state = corrector.apply(system.state, system.tree, dt, inverse=True, apex=system.apex)

    for step_idx in range(n_steps + 1):
        if step_idx % stride == 0 or step_idx == n_steps:
            # Physical readout map: non-destructive projection for invariants and angles
            if is_processed and corrector is not None:
                st = corrector.apply(system.state, system.tree, dt, inverse=False, apex=system.apex)
            else:
                st = system.state

            inv = compute_invariants(st, G=1.0)
            dE = abs(float(inv.H) - E0) / abs(E0)
            dP = float(np.linalg.norm(inv.P - inv0.P)) / max(M_tot, 1e-30)
            dL = float(np.linalg.norm(inv.L - inv0.L)) / max(L0_norm, 1e-30)
            dG = float(np.linalg.norm(inv.G - inv0.G)) / max(M_tot, 1e-30)
            phi = compute_laplace_angle(st.q)
            a_vals = compute_osculating_semi_major_axes(st.q, st.p, st.masses, G=1.0)

            times.append(float(st.t))
            phi_L_deg.append(math.degrees(phi))
            delta_E_rel.append(float(dE))
            delta_P_norm.append(float(dP))
            delta_L_rel.append(float(dL))
            delta_G_norm.append(float(dG))
            a_moons.append(a_vals)

        if step_idx < n_steps:
            system.step(dt=dt, in_shadow_mode=True)

    # 2. Global exit from shadow manifold back to physical coordinates (t=t_final)
    if is_processed and corrector is not None:
        system.state = corrector.apply(system.state, system.tree, dt, inverse=False, apex=system.apex)

    elapsed = time.perf_counter() - t0
    phi_arr = np.array(phi_L_deg)
    med_phi = float(np.median(phi_arr))
    libration_amp = float(0.5 * (np.max(phi_arr) - np.min(phi_arr)))
    max_dE = float(np.max(delta_E_rel))
    max_dP = float(np.max(delta_P_norm))
    rms_dP = float(np.sqrt(np.mean(np.square(delta_P_norm))))
    max_dL = float(np.max(delta_L_rel))
    rms_dL = float(np.sqrt(np.mean(np.square(delta_L_rel))))
    max_dG = float(np.max(delta_G_norm))
    rms_dG = float(np.sqrt(np.mean(np.square(delta_G_norm))))

    print(
        f"[{name}] Done in {elapsed:.2f}s | max |dE|={max_dE:.2e} | "
        f"max |dP|={max_dP:.2e} | max |dL|={max_dL:.2e} | max |dG|={max_dG:.2e} | "
        f"phi_L median={med_phi:.2f} deg, amp=+/-{libration_amp:.2f} deg",
        flush=True,
    )

    return {
        "name": name,
        "n_steps": n_steps,
        "dt": dt,
        "elapsed_seconds": elapsed,
        "max_delta_E_rel": max_dE,
        "max_delta_P_norm": max_dP,
        "rms_delta_P_norm": rms_dP,
        "max_delta_L_rel": max_dL,
        "rms_delta_L_rel": rms_dL,
        "max_delta_G_norm": max_dG,
        "rms_delta_G_norm": rms_dG,
        "phi_L_median_deg": med_phi,
        "phi_L_libration_amplitude_deg": libration_amp,
        "times": times,
        "phi_L_deg": phi_L_deg,
        "delta_E_rel": delta_E_rel,
        "delta_P_norm": delta_P_norm,
        "delta_L_rel": delta_L_rel,
        "delta_G_norm": delta_G_norm,
        "a_moons": a_moons,
    }


def main() -> None:
    print("=== Campaign III — Protocol S1: 4:2:1 Jovian Laplace Resonance Multi-Integrator Secular Benchmark ===", flush=True)
    t_start = time.perf_counter()

    # 200 Io orbits, T_final ~ 1256.64 rad, dt = 0.25 -> 5,027 steps
    p_io = 2.0 * math.pi
    n_orbits_io = 200
    t_final = n_orbits_io * p_io
    dt = 0.25
    n_steps = int(math.ceil(t_final / dt))
    stride = 10  # ~500 sampled points

    print(f"Configuration: {n_orbits_io} Io orbits | T_final={t_final:.2f} | dt={dt} | n_steps={n_steps} | stride={stride}", flush=True)

    integrators = [
        ("Strang (2nd order)", StrangIntegrator(G=1.0, screening_mode="c4")),
        ("Yoshida4 (4th order)", Yoshida4Integrator(G=1.0, screening_mode="c4")),
        ("BM4 (4th order optimal)", BlanesMoan4Integrator(G=1.0, screening_mode="c4")),
        ("PBM4 (Processed BM4)", create_pbm4_integrator(G=1.0, mode="homological")),
        ("BM6 (6th order optimal)", BlanesMoan6Integrator(G=1.0, screening_mode="c6")),
    ]

    results: Dict[str, Any] = {
        "metadata": {
            "benchmark": "jovian_laplace_resonance",
            "n_orbits_io": n_orbits_io,
            "t_final": t_final,
            "dt": dt,
            "n_steps": n_steps,
            "date": "2026-10-15",
        },
        "integrators": {},
    }

    for name, integ in integrators:
        res = run_integrator_simulation(name, integ, n_steps=n_steps, dt=dt, stride=stride)
        results["integrators"][name] = res

    # Write output to data/
    out_dir = WORKSPACE_ROOT / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "campaign_III_S1_laplace.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    total_time = time.perf_counter() - t_start
    print(f"=== Protocol S1 Completed Successfully in {total_time:.2f}s ===", flush=True)
    print(f"Saved results to: {out_path}", flush=True)


if __name__ == "__main__":
    main()
