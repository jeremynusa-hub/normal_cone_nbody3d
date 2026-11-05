#!/usr/bin/env python3
r"""
Campaign II — Protocol H2: Comet D/1993 F2 (Shoemaker-Levy 9) 1992 Tidal Disruption.
Demonstrates:
1. Passage inside Jupiter's fluid Roche limit (r_min ~ 1.33 R_J ~ 6.35e-4 AU).
2. Hill tidal field detection: eta_S(t) exceeding threshold eta_crit = 0.25.
3. Automatic supernode regularized evolution and tidal fragment elongation (string of pearls).
4. Energy conservation with zero softening (eps = 0, d0 = 0).
Outputs: data/campaign_II_H2_shoemaker_levy.json
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
    create_shoemaker_levy_9_1992_benchmark as create_sl9_f64,
)
from normal_cone_nbody_float64.potential.supernode_detector import (
    detect_supernodes as detect_sn_f64,
    compute_tidal_ratio,
)
from normal_cone_nbody_float64.core.invariants import compute_cartesian_hamiltonian as ham_f64

from normal_cone_nbody_float128.system import NBodySystem as NBodySystemF128
from normal_cone_nbody_float128.integrators import (
    create_pbm4_integrator as create_pbm4_f128,
)
from normal_cone_nbody_float128.benchmarks.registry import (
    create_shoemaker_levy_9_1992_benchmark as create_sl9_f128,
)
from normal_cone_nbody_float128.core.invariants import compute_cartesian_hamiltonian as ham_f128


from normal_cone_nbody_float128.potential.supernode_detector import (
    detect_supernodes as detect_sn_f128,
    compute_tidal_ratio as compute_tidal_ratio_f128,
)


def run_shoemaker_levy_simulation() -> Dict[str, Any]:
    t0 = time.perf_counter()
    print("=== [Protocol H2] Shoemaker-Levy 9 1992 Tidal Disruption (PURE FLOAT128) ===", flush=True)

    mp.mp.dps = 35
    dt = 0.00015
    n_steps = 135

    print(f"[H2.1] Evolving 5-body Sun-Jupiter-3-Fragment system in Float128 for {n_steps} steps (dt={dt})...", flush=True)
    sl9_f128 = create_sl9_f128(G=1.0, is_f128=True)
    sys_sl9_128 = NBodySystemF128(
        sl9_f128.with_updates(),
        G=1.0,
        integrator=create_pbm4_f128(G=1.0),
        max_dyadic_level=3,
    )

    times, r_jup_list, eta_s_list = [], [], []
    fragment_spans, rel_dEs, sn_active = [], [], []

    for s in range(n_steps):
        sys_sl9_128.step(dt=dt)
        q = sys_sl9_128.state.q
        p = sys_sl9_128.state.p
        m = sys_sl9_128.state.masses

        # Jovicentric distance of comet center (body 3 - body 1)
        r_j = float(mp.sqrt(sum((q[3][i] - q[1][i])**2 for i in range(3))))

        # Inter-fragment elongation span: max distance between fragments 2, 3, 4
        d23 = float(mp.sqrt(sum((q[2][i] - q[3][i])**2 for i in range(3))))
        d24 = float(mp.sqrt(sum((q[2][i] - q[4][i])**2 for i in range(3))))
        d34 = float(mp.sqrt(sum((q[3][i] - q[4][i])**2 for i in range(3))))
        span = max(d23, d24, d34)

        # Tidal ratio on cluster {2, 3, 4} from Jupiter (body 1)
        cluster = [2, 3, 4]
        M_S = sum(m[i] for i in cluster)
        q_cm = [sum(m[i] * q[i][dim] for i in cluster) / M_S for dim in range(3)]
        bar_rho_S = float(mp.sqrt(sum(m[i] * sum((q[i][dim] - q_cm[dim])**2 for dim in range(3)) for i in cluster) / M_S))
        r_P = float(mp.sqrt(sum((q_cm[dim] - q[1][dim])**2 for dim in range(3))))
        M_ext = float(m[1])
        eta_val = float(compute_tidal_ratio_f128(bar_rho_S, r_P, float(M_S), M_ext))

        sns = detect_sn_f128(sys_sl9_128.state, sys_sl9_128.tree)
        is_sn = len(sns) > 0

        E_cur = ham_f128(q, p, m, G=1.0)
        dE = float(abs(E_cur - sl9_f128.E0) / abs(sl9_f128.E0))

        times.append(float(sys_sl9_128.state.t))
        r_jup_list.append(r_j)
        eta_s_list.append(eta_val)
        fragment_spans.append(span)
        rel_dEs.append(dE)
        sn_active.append(is_sn)

        if (s + 1) % 30 == 0 or s == n_steps - 1:
            print(f"       Step {s + 1:3d}/{n_steps}: r_jup={r_j:.5f} AU | eta_S={eta_val:.3f} | span={span:.2e} | sn={is_sn} | dE/E0 (f128)={dE:.2e}", flush=True)

    min_rj = min(r_jup_list)
    max_eta = max(eta_s_list)
    initial_span = fragment_spans[0]
    final_span = fragment_spans[-1]
    elongation_ratio = final_span / max(initial_span, 1e-30)

    print(f"       Float128 SL9 complete: min r_jup = {min_rj:.5f} AU ({min_rj / 4.77e-4:.2f} R_J), max eta_S = {max_eta:.3f}, elongation = {elongation_ratio:.2f}x, max dE/E0 = {max(rel_dEs):.2e}", flush=True)

    res = {
        "precision": "float128_35_digits",
        "times": times,
        "r_jup": r_jup_list,
        "eta_s": eta_s_list,
        "fragment_spans": fragment_spans,
        "rel_energy_errors": rel_dEs,
        "supernode_active": sn_active,
        "min_r_jup": float(min_rj),
        "max_eta_s": float(max_eta),
        "elongation_ratio": float(elongation_ratio),
        "elapsed_seconds": time.perf_counter() - t0,
    }

    out_dir = WORKSPACE_ROOT / "results" / "campaigns"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "campaign_II_H2_shoemaker_levy.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    print(f"=== [Protocol H2 COMPLETE] Saved pure Float128 data to {out_file} in {res['elapsed_seconds']:.2f}s ===", flush=True)
    return res


if __name__ == "__main__":
    run_shoemaker_levy_simulation()
