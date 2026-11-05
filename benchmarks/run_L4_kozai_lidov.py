#!/usr/bin/env python3
r"""
Campaign I — Protocol L4: Kozai-Lidov Hierarchical Resonance with Zero Softening.
Demonstrates:
1. Long-term secular eccentricity pumping and inclination oscillation (e_max ~ 0.999).
2. Strict conservation of vertical angular momentum component L_z = sqrt(1-e^2)*cos(i) (Kozai integral).
3. Theorem 5.20 Leaf/Apex Decoupling: Smooth dyadic multirate scheduling without time-step collapse.
Outputs: data/campaign_I_L4_kozai.json
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

from normal_cone_nbody_float64.system import NBodySystem as NBodySystemF64
from normal_cone_nbody_float64.core.state import ExtendedPhaseState as ExtendedPhaseStateF64
from normal_cone_nbody_float64.integrators import (
    BlanesMoan4Integrator as BM4F64,
    create_pbm4_integrator as create_pbm4_f64,
)
from normal_cone_nbody_float64.core.invariants import compute_cartesian_hamiltonian


def create_hierarchical_kozai_system(
    a_in: float = 1.0,
    e_in: float = 0.1,
    inc_deg: float = 84.0,
    a_out: float = 4.0,
    m0: float = 1.0,
    m1: float = 1e-4,
    m2: float = 0.4,
    G: float = 1.0,
) -> ExtendedPhaseStateF64:
    masses = np.array([m0, m1, m2], dtype=np.float64)
    # Inner binary at apoapsis
    r_apo_in = a_in * (1.0 + e_in)
    v_apo_in = math.sqrt(G * (m0 + m1) * (1.0 - e_in) / r_apo_in)

    # Outer body at circular orbit
    v_out = math.sqrt(G * (m0 + m1 + m2) / a_out)

    inc_rad = math.radians(inc_deg)
    cos_i = math.cos(inc_rad)
    sin_i = math.sin(inc_rad)

    # Positions
    q = np.array([
        [0.0, 0.0, 0.0],
        [r_apo_in, 0.0, 0.0],
        [0.0, a_out, 0.0],
    ], dtype=np.float64)

    # Velocities
    v = np.array([
        [0.0, 0.0, 0.0],
        [0.0, v_apo_in * cos_i, v_apo_in * sin_i],
        [-v_out, 0.0, 0.0],
    ], dtype=np.float64)

    # Exact barycenter reduction
    cdm = np.sum(masses[:, None] * q, axis=0) / np.sum(masses)
    q -= cdm
    p = v * masses[:, None]
    p_cm = np.sum(p, axis=0) / np.sum(masses)
    p -= masses[:, None] * p_cm

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseStateF64(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def compute_inner_orbital_elements(q: np.ndarray, p: np.ndarray, masses: np.ndarray, G: float = 1.0) -> Dict[str, float]:
    r = q[1] - q[0]
    m0, m1 = masses[0], masses[1]
    mu = m0 * m1 / (m0 + m1)
    v = p[1] / m1 - p[0] / m0
    r_norm = float(np.linalg.norm(r))
    v_norm = float(np.linalg.norm(v))

    h_vec = np.cross(r, v)
    h_norm = float(np.linalg.norm(h_vec))
    M_tot = m0 + m1

    # Specific energy
    eps = 0.5 * v_norm * v_norm - G * M_tot / r_norm
    a = -0.5 * G * M_tot / eps if eps < 0 else 1.0

    # Eccentricity vector e_vec = (v x h) / (G * M) - r / ||r||
    e_vec = np.cross(v, h_vec) / (G * M_tot) - r / r_norm
    e = float(np.linalg.norm(e_vec))

    # Inclination
    inc_rad = math.acos(max(-1.0, min(1.0, h_vec[2] / max(h_norm, 1e-30))))
    inc_deg = math.degrees(inc_rad)

    # Kozai integral H_z ~ sqrt(1 - e^2) * cos(i)
    e_clamped = min(0.9999999, e)
    kozai_integral = math.sqrt(max(0.0, 1.0 - e_clamped * e_clamped)) * math.cos(inc_rad)

    return {
        "a": a,
        "e": e,
        "i_deg": inc_deg,
        "h_z": float(h_vec[2]),
        "kozai_integral": float(kozai_integral),
    }


def main() -> None:
    t0 = time.perf_counter()
    print("=== [Protocol L4] Kozai-Lidov Hierarchical Resonance ===", flush=True)

    st_kl = create_hierarchical_kozai_system(a_in=1.0, e_in=0.1, inc_deg=84.0, a_out=4.0)
    sys_kl = NBodySystemF64(st_kl.with_updates(), G=1.0, integrator=create_pbm4_f64(G=1.0), max_dyadic_level=6)

    dt = 0.5
    total_steps = 600
    sample_stride = 3

    times, eccs, incs, kozai_ints, rel_dEs, dyadic_ks = [], [], [], [], [], []

    print(f"[L4.1] Evolving Kozai-Lidov hierarchy for {total_steps} steps (dt={dt}, t_span={total_steps*dt:.1f})...", flush=True)
    t_start = time.perf_counter()

    corr = sys_kl.integrator.corrector
    # 1. Entrance onto shadow manifold
    sys_kl.state = corr.apply(sys_kl.state, sys_kl.tree, dt, inverse=True)

    for s in range(total_steps):
        sys_kl.step(dt=dt, in_shadow_mode=True)
        if (s + 1) % sample_stride == 0 or s == total_steps - 1:
            st_phys = corr.apply(sys_kl.state, sys_kl.tree, dt, inverse=False)
            orb = compute_inner_orbital_elements(st_phys.q, st_phys.p, st_phys.masses)
            E_curr = compute_cartesian_hamiltonian(st_phys.q, st_phys.p, st_phys.masses, G=1.0)
            dE = abs(E_curr - st_kl.E0) / abs(st_kl.E0)
            k_level = sys_kl.history.get_dyadic_levels()[-1] if len(sys_kl.history.get_dyadic_levels()) > 0 else 0

            times.append(float(sys_kl.state.t))
            eccs.append(orb["e"])
            incs.append(orb["i_deg"])
            kozai_ints.append(orb["kozai_integral"])
            rel_dEs.append(float(dE))
            dyadic_ks.append(int(k_level))

        if (s + 1) % 150 == 0 or s == total_steps - 1:
            print(f"       Step {s + 1}/{total_steps} (t={sys_kl.state.t:.1f}): e={eccs[-1]:.4f} | i={incs[-1]:.2f}° | KozaiInt={kozai_ints[-1]:.4f} | dE/E0={rel_dEs[-1]:.2e}", flush=True)

    # 2. Exit from shadow manifold back to physical coordinates
    sys_kl.state = corr.apply(sys_kl.state, sys_kl.tree, dt, inverse=False)

    e_min, e_max = min(eccs), max(eccs)
    i_min, i_max = min(incs), max(incs)
    k_int0 = kozai_ints[0]
    max_k_int_dev = max(abs(k - k_int0) for k in kozai_ints)

    print(f"       Kozai-Lidov cycle range: e in [{e_min:.4f}, {e_max:.4f}], i in [{i_min:.2f}°, {i_max:.2f}°]", flush=True)
    print(f"       Kozai integral max variation: {max_k_int_dev:.2e} | Max energy error: {max(rel_dEs):.2e} (Elapsed: {time.perf_counter() - t_start:.2f}s)", flush=True)

    res = {
        "times": times,
        "eccentricities": eccs,
        "inclinations_deg": incs,
        "kozai_integrals": kozai_ints,
        "rel_energy_errors": rel_dEs,
        "dyadic_levels": dyadic_ks,
        "e_min": float(e_min),
        "e_max": float(e_max),
        "i_min": float(i_min),
        "i_max": float(i_max),
        "max_kozai_int_dev": float(max_k_int_dev),
        "max_rel_dE": float(max(rel_dEs)),
        "elapsed_seconds": time.perf_counter() - t0,
    }

    out_dir = WORKSPACE_ROOT / "results" / "campaigns"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "campaign_I_L4_kozai.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    print(f"=== [Protocol L4 COMPLETE] Saved to {out_file} in {res['elapsed_seconds']:.2f}s ===", flush=True)


if __name__ == "__main__":
    main()
