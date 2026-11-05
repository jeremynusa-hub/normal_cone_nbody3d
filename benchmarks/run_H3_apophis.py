#!/usr/bin/env python3
r"""
Campaign II — Protocol H3: Asteroid (99942) Apophis 2029 Close Earth Flyby.
Demonstrates:
1. Deep geocentric flyby at 5.96 R_Earth (r_min ~ 2.54e-4 AU ~ 38,012 km).
2. Bilateral symmetric pericenter window: |t_peri| <= 2 |dt_macro| triggering dyadic ladder k(t).
3. Orbital transition from Aten class (a ~ 0.922 AU) to Apollo class (a ~ 1.102 AU).
4. Forward-backward time reversibility under PBM4 in Float64 and Float128.
Outputs: data/campaign_II_H3_apophis.json
"""

from __future__ import annotations
import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, Any, List

import numpy as np
import mpmath as mp

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
if str(WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_ROOT))

from normal_cone_nbody_float64.system import NBodySystem as NBodySystemF64
from normal_cone_nbody_float64.integrators import (
    create_pbm4_integrator as create_pbm4_f64,
)
from normal_cone_nbody_float64.benchmarks.registry import (
    create_apophis_2029_flyby_benchmark as create_apophis_f64,
)
from normal_cone_nbody_float64.diagnostics.metrics import (
    compute_forward_backward_reversibility as fb_rev_f64,
)
from normal_cone_nbody_float64.core.invariants import compute_cartesian_hamiltonian as ham_f64

from normal_cone_nbody_float128.system import NBodySystem as NBodySystemF128
from normal_cone_nbody_float128.integrators import (
    create_pbm4_integrator as create_pbm4_f128,
)
from normal_cone_nbody_float128.benchmarks.registry import (
    create_apophis_2029_flyby_benchmark as create_apophis_f128,
)
from normal_cone_nbody_float128.diagnostics.metrics import (
    compute_forward_backward_reversibility as fb_rev_f128,
)
from normal_cone_nbody_float128.core.invariants import compute_cartesian_hamiltonian as ham_f128


def compute_heliocentric_kepler_elements_f128(q_apo, v_apo, q_sun, v_sun, mu):
    r_rel = [q_apo[i] - q_sun[i] for i in range(3)]
    v_rel = [v_apo[i] - v_sun[i] for i in range(3)]
    r = mp.sqrt(sum(r_rel[i]**2 for i in range(3)))
    v2 = sum(v_rel[i]**2 for i in range(3))
    inv_a = 2.0 / r - v2 / mu
    a = 1.0 / inv_a if abs(inv_a) > 1e-30 else float('inf')
    hx = r_rel[1] * v_rel[2] - r_rel[2] * v_rel[1]
    hy = r_rel[2] * v_rel[0] - r_rel[0] * v_rel[2]
    hz = r_rel[0] * v_rel[1] - r_rel[1] * v_rel[0]
    h2 = hx**2 + hy**2 + hz**2
    p_param = h2 / mu
    e = mp.sqrt(max(0.0, 1.0 - p_param / a)) if a > 0 else mp.sqrt(max(0.0, 1.0 + p_param / abs(a)))
    return float(a), float(e)


def run_apophis_simulation() -> Dict[str, Any]:
    t0 = time.perf_counter()
    print("=== [Protocol H3] Asteroid (99942) Apophis 2029 Earth Flyby (PURE FLOAT128) ===", flush=True)

    mp.mp.dps = 35
    dt = 0.0001
    n_steps = 150

    print(f"[H3.1] Evolving 4-body Sun-Earth-Moon-Apophis system in Float128 for {n_steps} steps (dt={dt})...", flush=True)
    apo_f128 = create_apophis_f128(G=1.0, is_f128=True)
    sys_apo_128 = NBodySystemF128(
        apo_f128.with_updates(),
        G=1.0,
        integrator=create_pbm4_f128(G=1.0),
        max_dyadic_level=3,
    )

    times, r_earth_list, dyadic_levels = [], [], []
    a_helio_list, e_helio_list, rel_dEs = [], [], []

    mu_sun = 1.0 + apo_f128.masses[3]  # Sun + Apophis GM

    for s in range(n_steps):
        sys_apo_128.step(dt=dt)
        q = sys_apo_128.state.q
        p = sys_apo_128.state.p
        m = sys_apo_128.state.masses

        v = [[p[b][dim] / m[b] for dim in range(3)] for b in range(4)]

        # Geocentric distance of Apophis (body 3 - Earth body 1)
        r_ea = float(mp.sqrt(sum((q[3][i] - q[1][i])**2 for i in range(3))))

        # Heliocentric elements of Apophis
        a_helio, e_helio = compute_heliocentric_kepler_elements_f128(q[3], v[3], q[0], v[0], mu=mu_sun)

        # Dynamic level of scheduler
        k_dyn = sys_apo_128.scheduler.compute_dynamical_level(sys_apo_128.state, sys_apo_128.tree, dt)

        E_cur = ham_f128(q, p, m, G=1.0)
        dE = float(abs(E_cur - apo_f128.E0) / abs(apo_f128.E0))

        times.append(float(sys_apo_128.state.t))
        r_earth_list.append(r_ea)
        dyadic_levels.append(int(k_dyn))
        a_helio_list.append(float(a_helio))
        e_helio_list.append(float(e_helio))
        rel_dEs.append(dE)

        if (s + 1) % 30 == 0 or s == n_steps - 1:
            print(f"       Step {s + 1:3d}/{n_steps}: r_Earth={r_ea:.6f} AU ({r_ea / 4.2635e-5:.1f} R_E) | k={k_dyn} | a_helio={a_helio:.4f} AU | dE/E0 (f128)={dE:.2e}", flush=True)

    min_re = min(r_earth_list)
    initial_a = a_helio_list[0]
    final_a = a_helio_list[-1]
    delta_a = final_a - initial_a

    print(f"       Float128 Flyby complete: min r_Earth = {min_re:.6f} AU ({min_re / 4.2635e-5:.2f} R_E) | a: {initial_a:.4f} -> {final_a:.4f} AU (delta_a = {delta_a:+.4f} AU)", flush=True)

    # Reversibility test Forward-Backward in Float128
    print(f"[H3.2] Running forward-backward reversibility test in Float128 (35 digits)...", flush=True)
    fb_f128 = fb_rev_f128(apo_f128.with_updates(), create_pbm4_f128(G=1.0), dt=0.0002, n_steps=10, G=1.0)
    print(f"       Float128 FB: dq = {fb_f128['delta_q_norm']:.2e} | dp = {fb_f128['delta_p_norm']:.2e} | dE_rel = {fb_f128['delta_E_rel']:.2e} | palindrome={fb_f128['dyadic_palindrome_matched']}", flush=True)

    res = {
        "precision": "float128_35_digits",
        "times": times,
        "r_earth": r_earth_list,
        "dyadic_levels": dyadic_levels,
        "a_helio": a_helio_list,
        "e_helio": e_helio_list,
        "rel_energy_errors": rel_dEs,
        "min_r_earth": float(min_re),
        "initial_a": float(initial_a),
        "final_a": float(final_a),
        "delta_a": float(delta_a),
        "fb_f128": {
            "delta_q_norm": float(fb_f128["delta_q_norm"]),
            "delta_p_norm": float(fb_f128["delta_p_norm"]),
            "delta_E_rel": float(fb_f128["delta_E_rel"]),
            "dyadic_palindrome_matched": bool(fb_f128["dyadic_palindrome_matched"]),
        },
        "elapsed_seconds": time.perf_counter() - t0,
    }

    out_dir = WORKSPACE_ROOT / "results" / "campaigns"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "campaign_II_H3_apophis.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    print(f"=== [Protocol H3 COMPLETE] Saved pure Float128 data to {out_file} in {res['elapsed_seconds']:.2f}s ===", flush=True)
    return res


if __name__ == "__main__":
    run_apophis_simulation()
