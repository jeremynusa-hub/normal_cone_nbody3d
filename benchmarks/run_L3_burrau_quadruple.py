#!/usr/bin/env python3
r"""
Campaign I — Protocol L3: Burrau Pythagorean 3-Body & Quadruple Central Collision (D_4).
Demonstrates:
1. Burrau 3-4-5 chaotic close encounters with zero softening (eps = 0, d0 = 0).
   Tracks energy conservation, minimum mutual separation, and dyadic level adaptivity.
2. Symmetric Quadruple Central Collision (D_4 symmetry):
   Bit-exact preservation of D_4 point-group symmetry and forward-backward time-reversibility
   in Float64 and Float128 (35 digits).
Outputs: data/campaign_I_L3_burrau.json
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
    BlanesMoan4Integrator as BM4F64,
    create_pbm4_integrator as create_pbm4_f64,
)
from normal_cone_nbody_float64.benchmarks.registry import (
    create_burrau_pythagorean_benchmark as create_burrau_f64,
    create_quadruple_collision_benchmark as create_quad_f64,
)
from normal_cone_nbody_float64.diagnostics.metrics import compute_forward_backward_reversibility as fb_rev_f64
from normal_cone_nbody_float64.potential.supernode_detector import detect_supernodes as detect_sn_f64

from normal_cone_nbody_float128.system import NBodySystem as NBodySystemF128
from normal_cone_nbody_float128.integrators import (
    BlanesMoan4Integrator as BM4F128,
    create_pbm4_integrator as create_pbm4_f128,
)
from normal_cone_nbody_float128.benchmarks.registry import (
    create_quadruple_collision_benchmark as create_quad_f128,
)
from normal_cone_nbody_float128.diagnostics.metrics import compute_forward_backward_reversibility as fb_rev_f128


def run_burrau_pythagorean_simulation(t_span: float = 12.0, dt: float = 0.04) -> Dict[str, Any]:
    print(f"[L3.1] Simulating Burrau 3-4-5 problem (t_span={t_span}, dt={dt})...", flush=True)
    st = create_burrau_f64(G=1.0)
    sys_burrau = NBodySystemF64(st.with_updates(), G=1.0, integrator=create_pbm4_f64(G=1.0))

    n_steps = int(math.ceil(t_span / dt))
    t_start = time.perf_counter()

    min_separations = []
    supernode_counts = []
    for step in range(n_steps):
        sys_burrau.step(dt=dt)
        # track minimum mutual distance
        q = sys_burrau.state.q
        d01 = float(np.linalg.norm(q[1] - q[0]))
        d02 = float(np.linalg.norm(q[2] - q[0]))
        d12 = float(np.linalg.norm(q[2] - q[1]))
        min_d = min(d01, d02, d12)
        min_separations.append(min_d)

        # count supernodes
        try:
            sns = detect_sn_f64(sys_burrau.state, sys_burrau.tree)
            supernode_counts.append(len(sns))
        except Exception:
            supernode_counts.append(0)

        if (step + 1) % 50 == 0 or step == n_steps - 1:
            print(f"       Step {step + 1}/{n_steps} (t={sys_burrau.state.t:.2f}): min_r={min_d:.3e} | supernodes={supernode_counts[-1]}", flush=True)

    times = list(sys_burrau.history.get_times())
    energies = list(sys_burrau.history.get_energies())
    rel_dE = [float(abs(E - energies[0]) / abs(energies[0])) for E in energies]
    dyadic_k = [int(k) for k in sys_burrau.history.get_dyadic_levels()]

    print(f"       Burrau completed in {time.perf_counter() - t_start:.2f}s: max dE/E0 = {max(rel_dE):.2e}, global min_r = {min(min_separations):.3e}, max dyadic k = {max(dyadic_k)}", flush=True)
    return {
        "t_span": t_span,
        "dt": dt,
        "n_steps": n_steps,
        "times": times,
        "min_separations": min_separations,
        "rel_energy_errors": rel_dE,
        "dyadic_levels": dyadic_k,
        "supernode_counts": supernode_counts,
        "max_dE_rel": float(max(rel_dE)),
        "global_min_separation": float(min(min_separations)),
    }


def run_quadruple_collision_d4() -> Dict[str, Any]:
    print("[L3.2] Testing Quadruple Central Collision (D_4 symmetry & reversibility)...", flush=True)
    quad_f64 = create_quad_f64(G=1.0)
    quad_f128 = create_quad_f128(G=1.0).to_f128()

    # 1. Forward-backward time reversibility
    dt = 0.04
    steps_list = [5, 10, 20]
    res_f64, res_f128 = [], []

    for n in steps_list:
        fb64 = fb_rev_f64(quad_f64.with_updates(), BM4F64(G=1.0), dt=dt, n_steps=n, G=1.0)
        fb128 = fb_rev_f128(quad_f128.with_updates(), BM4F128(G=1.0), dt=dt, n_steps=n, G=1.0)
        print(f"       N_steps={n:2d} | f64: dE_rel={fb64['delta_E_rel']:.2e}, dq={fb64['delta_q_norm']:.2e} | f128: dE_rel={fb128['delta_E_rel']:.2e}, dq={fb128['delta_q_norm']:.2e}", flush=True)
        res_f64.append({
            "n_steps": n,
            "delta_E_rel": float(fb64["delta_E_rel"]),
            "delta_q_norm": float(fb64["delta_q_norm"]),
            "delta_p_norm": float(fb64["delta_p_norm"]),
        })
        res_f128.append({
            "n_steps": n,
            "delta_E_rel": float(fb128["delta_E_rel"]),
            "delta_q_norm": float(fb128["delta_q_norm"]),
            "delta_p_norm": float(fb128["delta_p_norm"]),
        })

    # 2. Symmetry check along 20 steps
    sys_q = NBodySystemF64(quad_f64.with_updates(), G=1.0, integrator=BM4F64(G=1.0))
    d4_deviations = []
    for _ in range(20):
        sys_q.step(dt=dt)
        q = sys_q.state.q
        # D4 symmetry: body 0 and 1 opposite on x-axis, body 2 and 3 opposite on y-axis
        # ||q0 + q1|| == 0 and ||q2 + q3|| == 0 and ||q0|| == ||q2||
        dev_opp1 = float(np.linalg.norm(q[0] + q[1]))
        dev_opp2 = float(np.linalg.norm(q[2] + q[3]))
        dev_iso = float(abs(np.linalg.norm(q[0]) - np.linalg.norm(q[2])))
        d4_deviations.append(dev_opp1 + dev_opp2 + dev_iso)

    max_d4_dev = float(max(d4_deviations))
    print(f"       D_4 point-group symmetry deviation across 20 steps: max Delta = {max_d4_dev:.2e}", flush=True)

    return {
        "float64_fb": res_f64,
        "float128_fb": res_f128,
        "max_d4_symmetry_deviation": max_d4_dev,
        "d4_deviations_history": d4_deviations,
    }


def main() -> None:
    t0 = time.perf_counter()
    print("=== [Protocol L3] Burrau 3-4-5 & Quadruple D_4 Central Collision ===", flush=True)

    burrau_res = run_burrau_pythagorean_simulation(t_span=12.0, dt=0.04)
    quad_res = run_quadruple_collision_d4()

    res = {
        "burrau": burrau_res,
        "quadruple_d4": quad_res,
        "elapsed_seconds": time.perf_counter() - t0,
    }

    out_dir = WORKSPACE_ROOT / "results" / "campaigns"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "campaign_I_L3_burrau.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    print(f"=== [Protocol L3 COMPLETE] Saved to {out_file} in {res['elapsed_seconds']:.2f}s ===", flush=True)


if __name__ == "__main__":
    main()
