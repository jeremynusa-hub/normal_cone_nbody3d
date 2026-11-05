#!/usr/bin/env python3
r"""
Campaign IV — Dual Precision Convergence Sweep (Float64 vs Float128 35 Digits)
and Experimental Validation of Theorem 5.15 (PBM4 Quadratic Floor O(epsilon^2 h^4) Crossover).

According to Farrés & Laskar (2013, Appendix D, Celest. Mech. Dyn. Astron.):
1. Dual Precision Convergence on Figure-8 Choreography:
   - Measures log10(Delta H / H0) and log10(||Delta q||_2) vs log10(tau / s) where tau/s = T / N_eval = h / s.
   - Compares Float64 (saturating at 1e-16 with Brouwer (tau/s)^(-1/2) walk)
     vs Float128 (pure slopes p = 2, 4, 6 down to 1e-32).
2. Theorem 5.15 Graphical Proof (PBM4 Order 6 -> 4 Transition):
   - Demonstrates the pre-asymptotic order 6 slope (C6 * eps * h^6) for h >> h_cross
   - Demonstrates the crossover knee at h_cross proportional to eps
   - Demonstrates the asymptotic order 4 slope governed by the quadratic floor C4 * eps^2 * h^4 for h << h_cross.

Outputs:
  data/campaign_IV_dual_convergence.json
"""

from __future__ import annotations
import json
import math
import time
from pathlib import Path
from typing import Dict, Any, List

import numpy as np
import mpmath as mp

# Workspace setup
WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = WORKSPACE_ROOT / "data"

import sys
if str(WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_ROOT))

# Import Float64 infrastructure
from normal_cone_nbody_float64.benchmarks.registry import create_figure8_choreography_benchmark as b64_fig8
from normal_cone_nbody_float64.system.nbody_system import NBodySystem as Sys64
from normal_cone_nbody_float64.core.invariants import compute_invariants as inv64
from normal_cone_nbody_float64.integrators.strang import StrangIntegrator as Strang64
from normal_cone_nbody_float64.integrators.yoshida import Yoshida4Integrator as Y4_64
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator as BM4_64
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator as BM6_64
from normal_cone_nbody_float64.integrators.processed import (
    ProcessedSymplecticIntegrator as PBM4_64,
    create_pbm4_integrator as create_pbm4_64,
)

# Import Float128 infrastructure
from normal_cone_nbody_float128.benchmarks.registry import create_figure8_choreography_benchmark as b128_fig8
from normal_cone_nbody_float128.system.nbody_system import NBodySystem as Sys128
from normal_cone_nbody_float128.core.invariants import compute_invariants as inv128
from normal_cone_nbody_float128.core.state import ExtendedPhaseState as State128
from normal_cone_nbody_float128.core.jacobi import JacobiTree, JacobiTopologyNode
from normal_cone_nbody_float128.integrators.strang import StrangIntegrator as Strang128
from normal_cone_nbody_float128.integrators.yoshida import Yoshida4Integrator as Y4_128
from normal_cone_nbody_float128.integrators.blanes_moan_4 import BlanesMoan4Integrator as BM4_128
from normal_cone_nbody_float128.integrators.blanes_moan_6 import BlanesMoan6Integrator as BM6_128
from normal_cone_nbody_float128.integrators.correctors import (
    GlobalCorrectorC4,
    HomologicalLieCorrectorC4,
)
from normal_cone_nbody_float128.integrators.processed import (
    ProcessedSymplecticIntegrator as PBM4_128,
    create_pbm4_integrator as create_pbm4_128,
)


def make_hierarchical_system_f128(eps_val: float) -> tuple[State128, JacobiTree, Any]:
    """
    Constructs a bound hierarchical 3-body system parameterized by perturbation strength eps:
    Inner circular binary (m0=1, m1=1, r=1, v_rel = sqrt(G*M01/r)), Outer perturber at distance R(eps)
    with mass m2 = 0.5 such that the tidal ratio epsilon = (m2/M01) * (r/R)^3 matches eps_val.
    """
    m0 = mp.mpf('1.0')
    m1 = mp.mpf('1.0')
    G = mp.mpf('1.0')
    r_bin = mp.mpf('1.0')
    M01 = m0 + m1
    mu01 = (m0 * m1) / M01
    v_rel_in = mp.sqrt(G * M01 / r_bin)

    m2 = mp.mpf('0.5')
    R_out = r_bin * (m2 / (M01 * mp.mpf(str(eps_val)))) ** (mp.mpf('1.0') / mp.mpf('3.0'))
    M_tot = M01 + m2
    mu_tot = (M01 * m2) / M_tot
    v_rel_out = mp.sqrt(G * M_tot / R_out)

    c45 = mp.cos(mp.pi / 4)
    s45 = mp.sin(mp.pi / 4)
    half = mp.mpf('0.5')

    q = np.array([
        [-half * r_bin * c45, -half * r_bin * s45, mp.mpf('0.0')],
        [half * r_bin * c45, half * r_bin * s45, mp.mpf('0.0')],
        [mp.mpf('0.0'), R_out, mp.mpf('0.0')],
    ], dtype=object)

    v = np.array([
        [half * v_rel_in * s45, -half * v_rel_in * c45, mp.mpf('0.0')],
        [-half * v_rel_in * s45, half * v_rel_in * c45, mp.mpf('0.0')],
        [-v_rel_out, mp.mpf('0.0'), mp.mpf('0.0')],
    ], dtype=object)

    masses = np.array([m0, m1, m2], dtype=object)
    for d in range(3):
        v_cm = (m0 * v[0, d] + m1 * v[1, d] + m2 * v[2, d]) / M_tot
        q_cm = (m0 * q[0, d] + m1 * q[1, d] + m2 * q[2, d]) / M_tot
        for i in range(3):
            v[i, d] -= v_cm
            q[i, d] -= q_cm

    p = np.array([[masses[i] * v[i, d] for d in range(3)] for i in range(3)], dtype=object)
    state = State128(q=q, p=p, masses=masses)

    l0 = JacobiTopologyNode(id=0, is_leaf=True, body_index=0, M=m0, body_indices=(0,))
    l1 = JacobiTopologyNode(id=1, is_leaf=True, body_index=1, M=m1, body_indices=(1,))
    l2 = JacobiTopologyNode(id=2, is_leaf=True, body_index=2, M=m2, body_indices=(2,))
    p01 = JacobiTopologyNode(id=10, is_leaf=False, left=l0, right=l1, body_indices=(0, 1),
                             M=M01, mu=mu01, c_left=m0 / M01, c_right=m1 / M01)
    root = JacobiTopologyNode(id=11, is_leaf=False, left=p01, right=l2, body_indices=(0, 1, 2),
                              M=M_tot, mu=mu_tot, c_left=M01 / M_tot, c_right=m2 / M_tot)
    tree = JacobiTree(root=root, masses=masses, node_map={0: l0, 1: l1, 2: l2, 10: p01, 11: root},
                      internal_nodes=(p01, root), leaf_pairs=(p01,), composite_nodes=(root,))
    return state, tree, R_out


def run_figure8_sweep() -> Dict[str, Any]:
    print("=== Running Part 1: Figure-8 Dual Precision Convergence Sweep ===", flush=True)
    T_target = 6.0  # Duration for convergence measurement

    # Stages s (force evaluations per macro step)
    stages_map = {
        "Strang (2nd order)": 1,
        "Yoshida4 (4th order)": 3,
        "BM4 (4th order optimal)": 6,
        "PBM4 (Processed BM4)": 6,
        "BM6 (6th order optimal)": 10,
    }

    # 1. Compute reference solution in Float64 with BM6 at very small step dt=0.0005
    print("Computing high-precision reference trajectory (BM6, dt=0.0005)...", flush=True)
    st_ref = b64_fig8()
    sys_ref = Sys64(state=st_ref, integrator=BM6_64(), G=1.0, desingularization_mode="leaf_only", enable_dyadic_substepping=False)
    n_ref = int(round(T_target / 0.0005))
    for _ in range(n_ref):
        sys_ref.step(dt=0.0005)
    q_ref_f64 = np.copy(sys_ref.state.q)
    H_ref_f64 = inv64(st_ref, G=1.0).H

    # 2. Float64 Sweep (broad range of steps to capture both order slope and Brouwer floor)
    f64_steps = [0.16, 0.08, 0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625, 0.0003125]
    f64_results: Dict[str, Any] = {}

    integrators_f64 = {
        "Strang (2nd order)": lambda: Strang64(),
        "Yoshida4 (4th order)": lambda: Y4_64(),
        "BM4 (4th order optimal)": lambda: BM4_64(),
        "PBM4 (Processed BM4)": lambda: create_pbm4_64(mode="homological"),
        "BM6 (6th order optimal)": lambda: BM6_64(),
    }

    for name, factory in integrators_f64.items():
        s = stages_map[name]
        tau_over_s_list = []
        dE_list = []
        dq_list = []
        print(f"  [Float64] Simulating {name} (s={s})...", flush=True)

        for h in f64_steps:
            st = b64_fig8()
            integ = factory()
            is_proc = isinstance(integ, PBM4_64)
            corr = integ.corrector if is_proc else None
            sys = Sys64(state=st, integrator=integ, G=1.0, desingularization_mode="leaf_only", enable_dyadic_substepping=False)
            n_steps = int(round(T_target / h))

            # Proper canonical telescoping
            if is_proc and corr is not None:
                sys.state = corr.apply(sys.state, sys.tree, h, inverse=True, apex=sys.apex)
            for _ in range(n_steps):
                sys.step(dt=h, in_shadow_mode=True)
            if is_proc and corr is not None:
                sys.state = corr.apply(sys.state, sys.tree, h, inverse=False, apex=sys.apex)

            H1 = inv64(sys.state, G=1.0).H
            dE = abs(H1 - H_ref_f64) / abs(H_ref_f64)
            dq = float(np.linalg.norm(sys.state.q - q_ref_f64))
            tau_over_s = h / s

            tau_over_s_list.append(tau_over_s)
            dE_list.append(max(float(dE), 1e-16))
            dq_list.append(max(float(dq), 1e-16))

        f64_results[name] = {
            "stages_s": s,
            "h_steps": f64_steps,
            "tau_over_s": tau_over_s_list,
            "delta_E_rel": dE_list,
            "delta_q_norm": dq_list,
        }

    # 3. Float128 Sweep (mpmath 35 dps, demonstrates asymptotic slopes down to 1e-32)
    print("=== [Float128 35 Digits] Simulating Asymptotic Slopes ===", flush=True)
    mp.mp.dps = 35
    f128_steps = [mp.mpf('0.16'), mp.mpf('0.08'), mp.mpf('0.04'), mp.mpf('0.02'), mp.mpf('0.01')]
    f128_results: Dict[str, Any] = {}

    integrators_f128 = {
        "Strang (2nd order)": lambda: Strang128(),
        "Yoshida4 (4th order)": lambda: Y4_128(),
        "BM4 (4th order optimal)": lambda: BM4_128(),
        "PBM4 (Processed BM4)": lambda: create_pbm4_128(G=mp.mpf('1.0'), mode="homological"),
        "BM6 (6th order optimal)": lambda: BM6_128(),
    }

    # Reference H0 in Float128
    st_init_128 = b128_fig8()
    H0_128 = inv128(st_init_128, G=mp.mpf('1.0')).H

    for name, factory in integrators_f128.items():
        s = stages_map[name]
        tau_over_s_list = []
        dE_list = []
        print(f"  [Float128] Simulating {name} (s={s})...", flush=True)

        for h in f128_steps:
            st = b128_fig8()
            integ = factory()
            is_proc = isinstance(integ, PBM4_128)
            corr = integ.corrector if is_proc else None
            sys = Sys128(state=st, integrator=integ, G=mp.mpf('1.0'), desingularization_mode="leaf_only", enable_dyadic_substepping=False)
            n_steps = int(round(T_target / float(h)))

            # Proper canonical telescoping
            if is_proc and corr is not None:
                sys.state = corr.apply(sys.state, sys.tree, h, inverse=True, apex=sys.apex)
            for _ in range(n_steps):
                sys.step(dt=h, in_shadow_mode=True)
            if is_proc and corr is not None:
                sys.state = corr.apply(sys.state, sys.tree, h, inverse=False, apex=sys.apex)

            H1 = inv128(sys.state, G=mp.mpf('1.0')).H
            dE = abs(H1 - H0_128) / abs(H0_128)
            tau_over_s = float(h) / s

            tau_over_s_list.append(tau_over_s)
            dE_list.append(float(dE))

        f128_results[name] = {
            "stages_s": s,
            "h_steps": [float(h) for h in f128_steps],
            "tau_over_s": tau_over_s_list,
            "delta_E_rel": dE_list,
        }

    return {"float64": f64_results, "float128": f128_results}


def run_theorem_5_15_crossover_sweep() -> Dict[str, Any]:
    print("=== Running Part 2: Theorem 5.15 Quadratic Floor O(eps^2 h^4) Crossover Sweep ===", flush=True)
    mp.mp.dps = 35
    eps_values = [0.01, 0.001, 0.0001]
    h_steps = [mp.mpf('0.40'), mp.mpf('0.20'), mp.mpf('0.10'), mp.mpf('0.05'), mp.mpf('0.025'), mp.mpf('0.0125')]
    T_tot = mp.mpf('2.0')

    bm4 = BM4_128(G=mp.mpf('1.0'))
    corrector = HomologicalLieCorrectorC4(G=mp.mpf('1.0'))

    crossover_data: Dict[str, Any] = {}

    for eps in eps_values:
        eps_label = f"eps_{eps}"
        print(f"  Simulating Hierarchical System for epsilon = {eps}...", flush=True)
        st0, tree, R_out = make_hierarchical_system_f128(eps)
        H0 = inv128(st0, G=mp.mpf('1.0')).H

        dE_pbm4: List[float] = []
        dE_bm4: List[float] = []

        for h in h_steps:
            n_steps = int(round(T_tot / float(h)))
            stride = max(1, n_steps // 20)
            checkpoints = set(range(1, n_steps + 1, stride)) | {n_steps}

            st_sh = corrector.apply(State128(q=np.copy(st0.q), p=np.copy(st0.p), masses=st0.masses), tree, h, inverse=True)
            st_bm = State128(q=np.copy(st0.q), p=np.copy(st0.p), masses=st0.masses)

            max_pbm = mp.mpf('0.0')
            max_bm = mp.mpf('0.0')

            for s_idx in range(1, n_steps + 1):
                st_sh = bm4.step(st_sh, tree, dt=h)
                st_bm = bm4.step(st_bm, tree, dt=h)
                if s_idx in checkpoints:
                    st_phys = corrector.apply(st_sh, tree, h, inverse=False)
                    H_pbm = inv128(st_phys, G=mp.mpf('1.0')).H
                    H_bm = inv128(st_bm, G=mp.mpf('1.0')).H
                    max_pbm = max(max_pbm, abs(H_pbm - H0) / abs(H0))
                    max_bm = max(max_bm, abs(H_bm - H0) / abs(H0))

            dE_pbm4.append(float(max_pbm))
            dE_bm4.append(float(max_bm))

        # Local finite-difference slopes p_local = log(dE[i]/dE[i+1]) / log(h[i]/h[i+1])
        local_slopes_pbm4 = []
        for i in range(len(dE_pbm4) - 1):
            ratio_E = dE_pbm4[i] / dE_pbm4[i + 1]
            ratio_h = float(h_steps[i]) / float(h_steps[i + 1])
            slope = math.log(ratio_E) / math.log(ratio_h)
            local_slopes_pbm4.append(float(slope))

        crossover_data[eps_label] = {
            "epsilon": eps,
            "R_out": float(R_out),
            "h_steps": [float(h) for h in h_steps],
            "dE_pbm4": dE_pbm4,
            "dE_bm4": dE_bm4,
            "local_slopes_pbm4": local_slopes_pbm4,
        }

    return crossover_data


def main() -> None:
    t_start = time.perf_counter()
    print("=================================================================", flush=True)
    print("Dual Precision Convergence & Theorem 5.15 Execution", flush=True)
    print("=================================================================", flush=True)

    out_file = DATA_DIR / "campaign_IV_dual_convergence.json"
    crossover_only = "--crossover-only" in sys.argv

    if crossover_only and out_file.exists():
        with open(out_file, "r", encoding="utf-8") as f:
            existing_payload = json.load(f)
        fig8_data = existing_payload.get("figure8_convergence", {})
    else:
        fig8_data = run_figure8_sweep()

    crossover_data = run_theorem_5_15_crossover_sweep()

    full_payload = {
        "metadata": {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "benchmark_fig8": "figure8_choreography_T6",
            "benchmark_hierarchical": "hierarchical_3body_crossover",
            "elapsed_seconds": round(time.perf_counter() - t_start, 2),
        },
        "figure8_convergence": fig8_data,
        "theorem_5_15_crossover": crossover_data,
    }

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(full_payload, f, indent=2)

    print(f"\n=== Successfully saved convergence data to: {out_file} ===", flush=True)
    print(f"Total execution time: {time.perf_counter() - t_start:.2f} seconds.", flush=True)


if __name__ == "__main__":
    main()

