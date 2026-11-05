#!/usr/bin/env python3
r"""
Campaign I — Protocol L2: McGehee Collision Regularization and C_4 Kick Linear Scaling (r -> 0).
Demonstrates:
1. Exact McGehee regularized traversal of two-body head-on collision (r -> 0, eps = 0, d0 = 0).
2. Theorem 5.15: Linear cancellation of the boundary symplectic corrector kick ||Delta p_{C4}|| = O(r).
3. Forward-backward time-reversibility round-trip across collision passage in Float64 and Float128.
Outputs: data/campaign_I_L2_mcgehee.json
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

from normal_cone_nbody_float64.system import NBodySystem as NBodySystemF64, AtlasManager as AtlasManagerF64
from normal_cone_nbody_float64.core.state import ExtendedPhaseState as ExtendedPhaseStateF64
from normal_cone_nbody_float64.integrators import (
    BlanesMoan4Integrator as BM4F64,
    create_pbm4_integrator as create_pbm4_f64,
)
from normal_cone_nbody_float64.benchmarks.registry import create_radial_collision_benchmark as create_rad_f64
from normal_cone_nbody_float64.diagnostics.metrics import compute_forward_backward_reversibility as fb_rev_f64
from normal_cone_nbody_float64.integrators.correctors import compute_global_normal_cone_c4_kick as c4_kick_f64

from normal_cone_nbody_float128.system import NBodySystem as NBodySystemF128
from normal_cone_nbody_float128.integrators import (
    BlanesMoan4Integrator as BM4F128,
    create_pbm4_integrator as create_pbm4_f128,
)
from normal_cone_nbody_float128.benchmarks.registry import create_radial_collision_benchmark as create_rad_f128
from normal_cone_nbody_float128.diagnostics.metrics import compute_forward_backward_reversibility as fb_rev_f128


def run_mcgehee_collision_traversal() -> Dict[str, Any]:
    print("[L2.1] Simulating head-on radial collision traversal (r0 = 2.0 -> collision)...", flush=True)
    rad_st = create_rad_f64(r0=2.0, m0=1.0, m1=1.0, G=1.0)
    # Analytical collision time t_coll = pi * 2^{1.5} / 4 = pi * sqrt(2) / 2 ~ 2.221441469
    t_coll = math.pi * (2.0 ** 1.5) / 4.0
    dt = 0.02
    n_steps = int(math.ceil(t_coll / dt)) + 10  # traverse through and past collision

    sys_pbm4 = NBodySystemF64(rad_st.with_updates(), G=1.0, integrator=create_pbm4_f64(G=1.0))
    for s in range(n_steps):
        sys_pbm4.step(dt=dt)

    times = list(sys_pbm4.history.get_times())
    pos = sys_pbm4.history.get_positions()
    r_sep = [float(np.linalg.norm(pos[step, 1] - pos[step, 0])) for step in range(len(times))]
    energies = list(sys_pbm4.history.get_energies())
    dE_rel = [float(abs(E - energies[0]) / abs(energies[0])) for E in energies]
    min_r = min(r_sep)

    print(f"       Collision traversal completed: {len(times)} steps, min separation r_min = {min_r:.2e}, max dE/E0 = {max(dE_rel):.2e}", flush=True)
    return {
        "t_coll_analytical": float(t_coll),
        "times": times,
        "separations": r_sep,
        "rel_energy_errors": dE_rel,
        "min_separation": float(min_r),
    }


def run_c4_kick_scaling_sweep() -> Dict[str, Any]:
    print("[L2.2] Sweeping separation r in [1e-5, 2.0] for C_4 corrector kick ||Delta p_{C4}||...", flush=True)
    r_grid = np.logspace(-5, 0.3, 25)
    kick_norms = []
    atlas = AtlasManagerF64()

    for idx, r_val in enumerate(r_grid):
        # 3-body configuration: binary of mass (1, 1) separated by r, plus 3rd apex body of mass 0.5 at distance R = 5.0
        masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
        q = np.array([
            [-0.5 * r_val, 0.0, 0.0],
            [0.5 * r_val, 0.0, 0.0],
            [0.0, 5.0, 0.0],
        ], dtype=np.float64)
        p = np.zeros_like(q)
        st = ExtendedPhaseStateF64(q=q, p=p, masses=masses)
        tree = atlas.build_optimal_tree(st)
        dp_c4 = c4_kick_f64(st, tree, h=0.05, G=1.0)
        norm_dp = float(np.linalg.norm(dp_c4))
        kick_norms.append(norm_dp)

    # Compute empirical power-law exponent alpha: log(||dp||) ~ alpha * log(r) + const
    # fit in the linear regime r in [1e-5, 1e-1]
    mask = (r_grid >= 1e-5) & (r_grid <= 1e-1)
    log_r = np.log10(r_grid[mask])
    log_dp = np.log10(np.array(kick_norms)[mask] + 1e-35)
    coeffs = np.polyfit(log_r, log_dp, 1)
    alpha = float(coeffs[0])
    print(f"       C_4 kick scaling exponent alpha = {alpha:.4f} (Theorem 5.15 expects alpha = 1.0000)", flush=True)

    return {
        "separations": [float(r) for r in r_grid],
        "c4_kick_norms": kick_norms,
        "empirical_power_law_alpha": alpha,
    }


def run_forward_backward_collision() -> Dict[str, Any]:
    print("[L2.3] Testing Forward-Backward Time-Reversibility across collision passage...", flush=True)
    rad_f64 = create_rad_f64(r0=2.0, m0=1.0, m1=1.0, G=1.0)
    rad_f128 = create_rad_f128(r0=2.0, m0=1.0, m1=1.0, G=1.0).to_f128()

    steps_list = [10, 20, 40]
    dt = 0.02
    res_f64, res_f128 = [], []

    for n in steps_list:
        fb64 = fb_rev_f64(rad_f64.with_updates(), create_pbm4_f64(G=1.0), dt=dt, n_steps=n, G=1.0)
        fb128 = fb_rev_f128(rad_f128.with_updates(), create_pbm4_f128(G=1.0), dt=dt, n_steps=n, G=1.0)
        print(f"       N_steps={n:2d} | f64: dE_rel={fb64['delta_E_rel']:.2e}, dq={fb64['delta_q_norm']:.2e} | f128: dE_rel={fb128['delta_E_rel']:.2e}, dq={fb128['delta_q_norm']:.2e}", flush=True)
        res_f64.append({
            "n_steps": n,
            "delta_E_rel": float(fb64["delta_E_rel"]),
            "delta_q_norm": float(fb64["delta_q_norm"]),
            "delta_p_norm": float(fb64["delta_p_norm"]),
            "dyadic_palindrome_matched": bool(fb64["dyadic_palindrome_matched"]),
        })
        res_f128.append({
            "n_steps": n,
            "delta_E_rel": float(fb128["delta_E_rel"]),
            "delta_q_norm": float(fb128["delta_q_norm"]),
            "delta_p_norm": float(fb128["delta_p_norm"]),
            "dyadic_palindrome_matched": bool(fb128["dyadic_palindrome_matched"]),
        })

    return {
        "float64_fb": res_f64,
        "float128_fb": res_f128,
    }


def main() -> None:
    t0 = time.perf_counter()
    print("=== [Protocol L2] McGehee Regularization & C_4 Kick Scaling ===", flush=True)

    traversal = run_mcgehee_collision_traversal()
    c4_sweep = run_c4_kick_scaling_sweep()
    fb_reversibility = run_forward_backward_collision()

    res = {
        "traversal": traversal,
        "c4_kick_scaling": c4_sweep,
        "forward_backward": fb_reversibility,
        "elapsed_seconds": time.perf_counter() - t0,
    }

    out_dir = WORKSPACE_ROOT / "results" / "campaigns"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "campaign_I_L2_mcgehee.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    print(f"=== [Protocol L2 COMPLETE] Saved to {out_file} in {res['elapsed_seconds']:.2f}s ===", flush=True)


if __name__ == "__main__":
    main()
