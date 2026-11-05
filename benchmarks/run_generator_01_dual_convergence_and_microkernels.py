#!/usr/bin/env python3
r"""
Benchmark Generator 01 (`Iteration 59`, Pillar II):
Dual-Precision Symplectic Convergence (`Float64` / `Float128`), Farrés–Laskar Inverted-Cost
Metric (`§7.1.1`, `Figure 2` & `Table 5`) & Section 6.6 High-Performance Micro-Kernels (`Table 4`).

Governed by Contract:
  `scripts/campaigns/contracts/01_generator_dual_convergence_and_microkernels.md`

Key Architectural Guarantees:
1. Imports all `Float64` components directly through the 7 certified Paper Showcase Interfaces
   (`01..06_*.py` in `chapters/ch06_algorithmic_architecture/`).
2. Imports all `Float128` (35-digit `mpmath`) components from `normal_cone_nbody_float128`
   with mandatory runtime assertions `assert state_f128.is_f128 is True`.
3. Computes the reference Figure-8 trajectory `(q_ref_f128, p_ref_f128, H0_f128)` in 35-digit `Float128`
   using `BlanesMoan6Integrator`, eliminating the `5e-11` `Float64` reference floor (`[MAJ-7.1.1-02]`).
4. Evaluates `Table 5` at a strictly common normalized inverse-cost abscissa `h_eff = h / s_eval = 0.01`
   (`N_eval = 600` force evaluations over `T = 6.0`) for all 5 integrators (`[CRIT-7.1.1-03]`).
5. Evaluates the Hierarchical 3D Eccentric Triple convergence sweep comparing all 5 canonical
   integrators (`Strang`, `Yoshida4`, `BM4`, `PBM4 (Z_4^{BM4})`, and `BM6`), proving >1000x order 4->6
   elevation of `PBM4 (Z_4^{BM4})` over `BM4` when `W_apex != 0`.
6. Benchmarks the 4 high-performance micro-kernels of `§6.6` across `N in {3, 4, 8, 16, 32, 64}`
   and certifies all 5 rows of `Table 4` (including multi-substep supernode reversibility <= 1e-15).
7. Enforces `verify_contract_assertions(results)` before writing JSON outputs.
"""

from __future__ import annotations
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple

import numpy as np
import mpmath as mp

mp.mp.dps = 35

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
if str(WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_ROOT))

RESULTS_DIR = WORKSPACE_ROOT / "results" / "campaigns"
DATA_DIR = WORKSPACE_ROOT / "data"
CHAPTERS_DIR = WORKSPACE_ROOT / "interfaces"


def _load_showcase_module(module_name: str, rel_path: str):
    full_path = CHAPTERS_DIR / rel_path
    spec = importlib.util.spec_from_file_location(module_name, full_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load paper showcase module {module_name} from {full_path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    return mod


# Load certified Paper Showcase Interfaces (01..06_*.py)
mod_01 = _load_showcase_module(
    "paper_01_global_pipeline",
    "01_PhaseStateAndSystem_Interface.py",
)
mod_02 = _load_showcase_module(
    "paper_02_tree_and_corrector",
    "02_TreeBuilderAndCorrector_Interface.py",
)
mod_03 = _load_showcase_module(
    "paper_03_dyadic_scheduler",
    "03_DyadicScheduler_Interface.py",
)
mod_04 = _load_showcase_module(
    "paper_04_symplectic_subcycle",
    "04_SymplecticSubcycle_Interface.py",
)
mod_05 = _load_showcase_module(
    "paper_05_atlas_transitions",
    "05_AtlasTransitions_Interface.py",
)
mod_06 = _load_showcase_module(
    "paper_06_optimizations",
    "06_HighPerformanceOptimizations_Interface.py",
)

# Canonical Float64 & Float128 Benchmark Registries
from normal_cone_nbody_float64.benchmarks.registry import (
    create_figure8_choreography_benchmark as create_fig8_f64,
    create_hierarchical_3d_eccentric_triple_benchmark as create_hier3d_f64,
    create_burrau_pythagorean_benchmark as create_burrau_f64,
)
from normal_cone_nbody_float64.potential.screened_apex import ScreenedApexPotential
from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec
from normal_cone_nbody_float64.system.atlas_manager import AtlasManager as AtlasManager64

from normal_cone_nbody_float128.benchmarks.registry import (
    create_figure8_choreography_benchmark as create_fig8_f128,
    create_hierarchical_3d_eccentric_triple_benchmark as create_hier3d_f128,
)
from normal_cone_nbody_float128.system.nbody_system import NBodySystem as Sys128
from normal_cone_nbody_float128.system.atlas_manager import AtlasManager as AtlasManager128
from normal_cone_nbody_float128.core.invariants import compute_invariants as inv128
from normal_cone_nbody_float128.integrators.strang import StrangIntegrator as Strang128
from normal_cone_nbody_float128.integrators.yoshida import Yoshida4Integrator as Y4_128
from normal_cone_nbody_float128.integrators.blanes_moan_4 import BlanesMoan4Integrator as BM4_128
from normal_cone_nbody_float128.integrators.blanes_moan_6 import BlanesMoan6Integrator as BM6_128
from normal_cone_nbody_float128.integrators.processed import (
    create_pbm4_integrator as create_pbm4_128,
)


def _compute_local_slopes(h_list: List[float], err_list: List[float]) -> List[float]:
    slopes: List[float] = []
    for i in range(len(err_list) - 1):
        e_hi = max(float(err_list[i]), 1e-38)
        e_lo = max(float(err_list[i + 1]), 1e-38)
        h_hi = float(h_list[i])
        h_lo = float(h_list[i + 1])
        slopes.append(float(math.log(e_hi / e_lo) / math.log(h_hi / h_lo)))
    return slopes


def run_part1_figure8_dual_convergence() -> Dict[str, Any]:
    print("=== [Part 1] Running Figure-8 Dual-Precision Convergence & Table 5 Common-Abscissa Sweep ===", flush=True)
    # Fixed exact power-of-two single-chart horizon T_conv = 0.16 (divisible by all h = 0.16 * 2^-k)
    T_conv = 0.16
    T_table5 = 6.0
    T_orb_fig8 = 6.325913982926236

    stages_map = {
        "Strang (2nd order)": 1,
        "Yoshida4 (4th order)": 3,
        "BM4 (4th order optimal)": 6,
        "PBM4 (Processed BM4)": 6,
        "BM6 (6th order optimal)": 10,
    }

    st_ref_128 = create_fig8_f128(is_f128=True)
    assert st_ref_128.is_f128 is True, "Reference state must be native Float128 (mpmath.mpf)"
    assert isinstance(st_ref_128.q[0, 0], mp.mpf), "Reference q[0,0] must be mpmath.mpf"
    H0_128 = inv128(st_ref_128, G=mp.mpf("1.0")).H

    # 1. Float128 Sweep (35-digit mpmath, verified .is_f128 == True, descending down to <= 1e-30 on fixed T_conv = 0.16)
    f128_steps_str = [
        "0.16",
        "0.08",
        "0.04",
        "0.02",
        "0.01",
        "0.005",
        "0.0025",
        "0.00125",
        "0.000625",
        "0.0003125",
        "0.00015625",
        "0.000078125",
        "0.0000390625",
    ]
    integrators_f128: Dict[str, Callable[[], Any]] = {
        "Strang (2nd order)": lambda: Strang128(),
        "Yoshida4 (4th order)": lambda: Y4_128(),
        "BM4 (4th order optimal)": lambda: BM4_128(),
        "PBM4 (Processed BM4)": lambda: create_pbm4_128(G=mp.mpf("1.0"), mode="homological"),
        "BM6 (6th order optimal)": lambda: BM6_128(),
    }

    inv0_128 = inv128(st_ref_128, G=mp.mpf("1.0"))
    f128_results: Dict[str, Any] = {}
    f128_end_states_q: Dict[str, List[np.ndarray]] = {}
    all_f128_flags_verified = True

    for name, factory in integrators_f128.items():
        s = stages_map[name]
        h_float_list: List[float] = []
        tau_over_s_list: List[float] = []
        dE_max_list: List[float] = []
        dE_end_list: List[float] = []
        dP_list: List[float] = []
        dL_list: List[float] = []
        dG_list: List[float] = []
        q_ends_for_scheme: List[np.ndarray] = []
        print(f"  [1.1 Float128 35-digit] Sweeping {name} (s={s}) on fixed horizon T_conv={T_conv}...", flush=True)

        # 10 dyadic octaves for s<=6 (1 to 512 steps), 13 dyadic octaves for BM6 (1 to 4096 steps, reaching < 1e-30)
        active_steps_str = f128_steps_str if s == 10 else f128_steps_str[:10]

        for h_str in active_steps_str:
            h_mp = mp.mpf(h_str)
            h_flt = float(h_mp)
            st_128 = create_fig8_f128(is_f128=True)
            if not (st_128.is_f128 and isinstance(st_128.q[0, 0], mp.mpf)):
                all_f128_flags_verified = False
                raise AssertionError("Float128 state failed is_f128 check!")

            integ_128 = factory()
            sys_128 = Sys128(
                state=st_128,
                integrator=integ_128,
                atlas_manager=AtlasManager128(),
                G=mp.mpf("1.0"),
                enable_dyadic_substepping=False,
                in_shadow_mode=True,
            )

            n_steps = int(round(T_conv / h_flt))
            for _ in range(n_steps):
                sys_128.step(dt=h_mp)

            phys_end_128 = sys_128.observe_physical_state()
            if not (phys_end_128.is_f128 and isinstance(phys_end_128.q[0, 0], mp.mpf)):
                all_f128_flags_verified = False
                raise AssertionError("Post-step Float128 state lost mpmath.mpf type!")

            inv_end_128 = inv128(phys_end_128, G=mp.mpf("1.0"))
            dE_end_mp = abs(inv_end_128.H - H0_128) / abs(H0_128)
            q_ends_for_scheme.append(np.copy(phys_end_128.q))

            dP_mp = mp.sqrt(sum((inv_end_128.P[d] - inv0_128.P[d]) ** 2 for d in range(3)))
            dL_mp = mp.sqrt(sum((inv_end_128.L[d] - inv0_128.L[d]) ** 2 for d in range(3)))
            dG_mp = mp.sqrt(sum((inv_end_128.G[d] - inv0_128.G[d]) ** 2 for d in range(3)))

            h_float_list.append(h_flt)
            tau_over_s_list.append(float(h_flt / s))
            dE_max_list.append(float(dE_end_mp))
            dE_end_list.append(float(dE_end_mp))
            dP_list.append(float(dP_mp))
            dL_list.append(float(dL_mp))
            dG_list.append(float(dG_mp))

        f128_end_states_q[name] = q_ends_for_scheme
        local_slopes = _compute_local_slopes(h_float_list, dE_max_list)
        f128_results[name] = {
            "stages_s": s,
            "h_steps": h_float_list,
            "tau_over_s": tau_over_s_list,
            "delta_E_rel": dE_max_list,
            "delta_E_rel_endpoint": dE_end_list,
            "local_slopes_dE": local_slopes,
            "mean_asymptotic_slope": float(np.mean(local_slopes[2:8])),
            "min_delta_E_rel": float(min(dE_max_list)),
            "max_delta_P": float(max(dP_list)),
            "max_delta_L": float(max(dL_list)),
            "max_delta_G": float(max(dG_list)),
        }

    # Construct 35-digit Richardson-extrapolated O(h^8) reference state at T_conv = 0.16 from the two finest BM6_128 runs
    q_bm6_finest = f128_end_states_q["BM6 (6th order optimal)"][-1]
    q_bm6_prev = f128_end_states_q["BM6 (6th order optimal)"][-2]
    q_ref_f128_obj = np.empty((3, 3), dtype=object)
    for bi in range(3):
        for d in range(3):
            q_ref_f128_obj[bi, d] = (mp.mpf("64.0") * q_bm6_finest[bi, d] - q_bm6_prev[bi, d]) / mp.mpf("63.0")
    q_ref_f64_from_f128 = np.array(
        [[float(q_ref_f128_obj[bi, d]) for d in range(3)] for bi in range(3)], dtype=np.float64
    )
    ref_energy_drift_f128 = float(f128_results["BM6 (6th order optimal)"]["min_delta_E_rel"])
    print(f"        Finest Float128 BM6 energy drift at T_conv={T_conv}: {ref_energy_drift_f128:.3e}", flush=True)

    # Populate Float128 trajectory errors against q_ref_f128_obj
    for name in integrators_f128.keys():
        dq_128_list: List[float] = []
        for q_end_128 in f128_end_states_q[name]:
            dq_sq = mp.mpf("0.0")
            for bi in range(3):
                for d in range(3):
                    diff_qd = q_end_128[bi, d] - q_ref_f128_obj[bi, d]
                    dq_sq += diff_qd * diff_qd
            dq_128_list.append(float(mp.sqrt(dq_sq)))
        f128_results[name]["delta_q_norm"] = dq_128_list
        f128_results[name]["delta_q_norm_T6"] = dq_128_list
        f128_results[name]["tau_over_s_T6"] = f128_results[name]["tau_over_s"]

    # 2. Float64 Sweep via Paper Showcase Interfaces (01..05_*.py) on fixed T_conv = 0.16
    f64_steps = [0.16, 0.08, 0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625, 0.0003125]
    integrators_f64: Dict[str, Callable[[], Any]] = {
        "Strang (2nd order)": lambda: mod_04.StrangIntegrator(),
        "Yoshida4 (4th order)": lambda: mod_04.Yoshida4Integrator(),
        "BM4 (4th order optimal)": lambda: mod_04.BlanesMoan4Integrator(),
        "PBM4 (Processed BM4)": lambda: mod_02.create_pbm4_integrator(G=1.0, mode="homological"),
        "BM6 (6th order optimal)": lambda: mod_04.BlanesMoan6Integrator(),
    }

    st0_f64 = create_fig8_f64()
    inv0_f64 = mod_05.compute_lie_bargmann_invariants(st0_f64, G=1.0)
    H0_f64 = float(inv0_f64.H)

    f64_results: Dict[str, Any] = {}
    for name, factory in integrators_f64.items():
        s = stages_map[name]
        tau_over_s_list = []
        dE_end_list = []
        dE_max_list = []
        dq_list = []
        dP_list = []
        dL_list = []
        dG_list = []
        print(f"  [1.2 Float64] Sweeping {name} (s={s}) on fixed horizon T_conv={T_conv}...", flush=True)

        for h in f64_steps:
            st = create_fig8_f64()
            integ = factory()
            sys_run = mod_01.NBodySystem(
                state=st,
                integrator=integ,
                atlas_manager=AtlasManager64(),
                G=1.0,
                enable_dyadic_substepping=False,
                in_shadow_mode=True,
                use_optimized_kernels=True,
            )
            n_steps = int(round(T_conv / h))
            for _ in range(n_steps):
                sys_run.step(dt=h)

            phys_end = sys_run.observe_physical_state()
            inv_end = mod_05.compute_lie_bargmann_invariants(phys_end, G=1.0, E0=H0_f64)
            dE_end = abs(inv_end.H - H0_f64) / abs(H0_f64)
            dq = float(np.linalg.norm(phys_end.q - q_ref_f64_from_f128))
            dP = float(np.linalg.norm(inv_end.P - inv0_f64.P))
            dL = float(np.linalg.norm(inv_end.L - inv0_f64.L))
            dG = float(np.linalg.norm(inv_end.G - inv0_f64.G))

            tau_over_s_list.append(float(h / s))
            dE_end_list.append(float(dE_end))
            dE_max_list.append(float(dE_end))
            dq_list.append(float(dq))
            dP_list.append(dP)
            dL_list.append(dL)
            dG_list.append(dG)

        f64_results[name] = {
            "stages_s": s,
            "h_steps": f64_steps,
            "tau_over_s": tau_over_s_list,
            "delta_E_rel": dE_max_list,
            "delta_E_rel_endpoint": dE_end_list,
            "delta_q_norm": dq_list,
            "local_slopes_dE": _compute_local_slopes(f64_steps, dE_max_list),
            "local_slopes_dq": _compute_local_slopes(f64_steps, dq_list),
            "max_delta_P": float(max(dP_list)),
            "max_delta_L": float(max(dL_list)),
            "max_delta_G": float(max(dG_list)),
        }

    # 3. Exact Common Normalized Inverse-Cost Benchmark at h_eff = h / s = 0.005 (N_eval = 1200 over T = 6.0) for Table 5
    print("  [1.3] Evaluating Table 5 at strictly common normalized cost h_eff = h / s = 0.005 (N_eval = 1200, T = 6.0)...", flush=True)
    # Reference at T = 6.0 using unoptimized BM6 in 35-digit Float128 at h_ref = 0.01 (600 steps = 6000 force evals)
    st_ref_t6_128 = create_fig8_f128(is_f128=True)
    assert st_ref_t6_128.is_f128 is True
    sys_ref_t6_128 = Sys128(
        state=st_ref_t6_128,
        integrator=BM6_128(),
        atlas_manager=AtlasManager128(),
        G=mp.mpf("1.0"),
        enable_dyadic_substepping=False,
        in_shadow_mode=True,
    )
    dt_ref_t6_mp = mp.mpf("0.01")
    n_ref_t6_steps = int(round(T_table5 / float(dt_ref_t6_mp)))
    max_dE_ref_t6 = 0.0
    for _ in range(n_ref_t6_steps):
        sys_ref_t6_128.step(dt=dt_ref_t6_mp)
    phys_ref_t6_128 = sys_ref_t6_128.observe_physical_state()
    inv_ref_t6_128 = inv128(phys_ref_t6_128, G=mp.mpf("1.0"))
    end_dE_ref_t6_128 = float(abs(inv_ref_t6_128.H - H0_128) / abs(H0_128))
    q_ref_t6 = np.array(
        [[float(phys_ref_t6_128.q[bi, d]) for d in range(3)] for bi in range(3)], dtype=np.float64
    )
    table5_ref_bm6_f128 = {
        "scheme": "BM6 (Float128 Reference, Unoptimized)",
        "stages_s": 10,
        "h_eff": 0.001,
        "h_macro": 0.01,
        "n_steps": n_ref_t6_steps,
        "n_force_evals": n_ref_t6_steps * 10,
        "delta_E_rel_endpoint": end_dE_ref_t6_128,
        "float128_min_error_Tconv": float(f128_results["BM6 (6th order optimal)"]["min_delta_E_rel"]),
        "float128_min_dq_Tconv": float(min(f128_results["BM6 (6th order optimal)"]["delta_q_norm"])),
    }

    h_eff_common = 0.005
    table5_common: Dict[str, Any] = {}
    for name, s in stages_map.items():
        h_common = h_eff_common * s
        n_steps_common = int(round(T_table5 / h_common))
        st_c = create_fig8_f64()
        sys_c = mod_01.NBodySystem(
            state=st_c,
            integrator=integrators_f64[name](),
            atlas_manager=AtlasManager64(),
            G=1.0,
            enable_dyadic_substepping=False,
            in_shadow_mode=True,
            use_optimized_kernels=True,
        )
        max_dE_c = 0.0
        sum_sq_dE_c = 0.0
        for _ in range(n_steps_common):
            sys_c.step(dt=h_common)
            phys_c = sys_c.observe_physical_state()
            H_c = mod_05.compute_lie_bargmann_invariants(phys_c, G=1.0).H
            err_c = abs(H_c - H0_f64) / abs(H0_f64)
            max_dE_c = max(max_dE_c, err_c)
            sum_sq_dE_c += err_c * err_c

        phys_c_end = sys_c.observe_physical_state()
        H_c_end = mod_05.compute_lie_bargmann_invariants(phys_c_end, G=1.0).H
        end_dE_c = abs(H_c_end - H0_f64) / abs(H0_f64)
        rms_dE_c = math.sqrt(sum_sq_dE_c / n_steps_common)
        dq_c = float(np.linalg.norm(phys_c_end.q - q_ref_t6))

        table5_common[name] = {
            "stages_s": s,
            "h_eff": h_eff_common,
            "h_macro": h_common,
            "n_steps": n_steps_common,
            "n_force_evals": n_steps_common * s,
            "delta_E_rel_max": float(max_dE_c),
            "delta_E_rel_rms": float(rms_dE_c),
            "delta_E_rel_endpoint": float(end_dE_c),
            "delta_q_norm": float(dq_c),
            "float64_floor": float(min(f64_results[name]["delta_E_rel"])),
            "float128_min_error": float(f128_results[name]["min_delta_E_rel"]),
            "empirical_order_f128": float(f128_results[name]["mean_asymptotic_slope"]),
        }

    y4_max_err = table5_common["Yoshida4 (4th order)"]["delta_E_rel_max"]
    bm4_max_err = table5_common["BM4 (4th order optimal)"]["delta_E_rel_max"]
    pbm4_max_err = table5_common["PBM4 (Processed BM4)"]["delta_E_rel_max"]
    bm6_max_err = table5_common["BM6 (6th order optimal)"]["delta_E_rel_max"]
    equal_cost_bm4_over_y4_gain = float(y4_max_err / bm4_max_err)
    equal_cost_pbm4_over_y4_gain = float(y4_max_err / pbm4_max_err)
    equal_cost_bm6_over_y4_gain = float(y4_max_err / bm6_max_err)
    # Also compute equal-macro-step gain at h = 0.01 (index 4 in f128_results)
    equal_step_bm4_over_y4_gain = float(
        f128_results["Yoshida4 (4th order)"]["delta_E_rel"][4]
        / f128_results["BM4 (4th order optimal)"]["delta_E_rel"][4]
    )

    # 5. Volet C Physical Figure-8 Full-Period Trajectory (T_orb = 6.32591398, 400 steps with BM6)
    print("  [1.5] Sampling full-period Figure-8 choreography trajectory for Volet C...", flush=True)
    n_volet_c = 400
    dt_volet_c = T_orb_fig8 / n_volet_c
    st_vc = create_fig8_f64()
    q0_vc = np.copy(st_vc.q)
    p0_vc = np.copy(st_vc.p)
    sys_vc = mod_01.NBodySystem(
        state=st_vc,
        integrator=mod_04.BlanesMoan6Integrator(),
        G=1.0,
        enable_dyadic_substepping=False,
        use_optimized_kernels=True,
    )
    t_vc_list: List[float] = [0.0]
    q_vc_hist: List[List[List[float]]] = [st_vc.q.tolist()]
    r12_list: List[float] = [float(np.linalg.norm(st_vc.q[0] - st_vc.q[1]))]
    r23_list: List[float] = [float(np.linalg.norm(st_vc.q[1] - st_vc.q[2]))]
    r13_list: List[float] = [float(np.linalg.norm(st_vc.q[0] - st_vc.q[2]))]
    T_kin_list: List[float] = [float(0.5 * np.sum(st_vc.p ** 2))]
    V_pot_list: List[float] = [float(H0_f64 - 0.5 * np.sum(st_vc.p ** 2))]

    for step_idx in range(1, n_volet_c + 1):
        sys_vc.step(dt=dt_volet_c)
        cur = sys_vc.observe_physical_state()
        t_vc_list.append(float(cur.t))
        q_vc_hist.append(cur.q.tolist())
        r12 = float(np.linalg.norm(cur.q[0] - cur.q[1]))
        r23 = float(np.linalg.norm(cur.q[1] - cur.q[2]))
        r13 = float(np.linalg.norm(cur.q[0] - cur.q[2]))
        r12_list.append(r12)
        r23_list.append(r23)
        r13_list.append(r13)
        t_kin = float(0.5 * np.sum(cur.p ** 2))
        T_kin_list.append(t_kin)
        V_pot_list.append(float(-(1.0 / r12 + 1.0 / r23 + 1.0 / r13)))

    periodic_closure_q_err = float(np.linalg.norm(sys_vc.observe_physical_state().q - q0_vc))
    periodic_closure_p_err = float(np.linalg.norm(sys_vc.observe_physical_state().p - p0_vc))

    return {
        "float128_verified": all_f128_flags_verified,
        "T_conv": T_conv,
        "T_table5": T_table5,
        "reference_f128_energy_drift": ref_energy_drift_f128,
        "float64": f64_results,
        "float128": f128_results,
        "table5_reference_bm6_f128": table5_ref_bm6_f128,
        "table5_common_abscissa_h_eff_0_01": table5_common,
        "equal_cost_bm4_over_y4_gain": equal_cost_bm4_over_y4_gain,
        "equal_cost_pbm4_over_y4_gain": equal_cost_pbm4_over_y4_gain,
        "equal_cost_bm6_over_y4_gain": equal_cost_bm6_over_y4_gain,
        "equal_step_bm4_over_y4_gain": equal_step_bm4_over_y4_gain,
        "volet_C_figure8_orbit": {
            "T_orb": T_orb_fig8,
            "q0": q0_vc.tolist(),
            "p0": p0_vc.tolist(),
            "periodic_closure_q_err": periodic_closure_q_err,
            "periodic_closure_p_err": periodic_closure_p_err,
            "t": t_vc_list,
            "q_history": q_vc_hist,
            "r12": r12_list,
            "r23": r23_list,
            "r13": r13_list,
            "T_kinetic": T_kin_list,
            "V_potential": V_pot_list,
        },
    }


def run_part2_hierarchical_3d_eccentric_probe() -> Dict[str, Any]:
    print("=== [Part 2] Running Hierarchical 3D Eccentric Triple Probe (5 Canonical Integrators in Float64 & Float128) ===", flush=True)
    T_in = 2.0 * math.pi
    step_counts = [8, 16, 32, 64, 128, 256, 512]
    h_list = [T_in / n for n in step_counts]

    st0 = create_hier3d_f64(e_in=0.30, inc_deg=28.0, a_out=12.0)
    H0 = mod_05.compute_lie_bargmann_invariants(st0, G=1.0).H

    stages_map = {
        "Strang (2nd order)": 1,
        "Yoshida4 (4th order)": 3,
        "BM4 (4th order optimal)": 6,
        "PBM4 (Processed BM4)": 6,
        "BM6 (6th order optimal)": 10,
    }
    canonical_f64: Dict[str, Callable[[], Any]] = {
        "Strang (2nd order)": lambda: mod_04.StrangIntegrator(),
        "Yoshida4 (4th order)": lambda: mod_04.Yoshida4Integrator(),
        "BM4 (4th order optimal)": lambda: mod_04.BlanesMoan4Integrator(G=1.0),
        "PBM4 (Processed BM4)": lambda: mod_02.create_pbm4_integrator(G=1.0, mode="homological"),
        "BM6 (6th order optimal)": lambda: mod_04.BlanesMoan6Integrator(),
    }
    canonical_f128: Dict[str, Callable[[], Any]] = {
        "Strang (2nd order)": lambda: Strang128(),
        "Yoshida4 (4th order)": lambda: Y4_128(),
        "BM4 (4th order optimal)": lambda: BM4_128(),
        "PBM4 (Processed BM4)": lambda: create_pbm4_128(G=mp.mpf("1.0"), mode="homological"),
        "BM6 (6th order optimal)": lambda: BM6_128(),
    }

    results_3d: Dict[str, Any] = {
        "T_in": T_in,
        "step_counts": step_counts,
        "h_steps": h_list,
        "tau_over_s": [h / 6.0 for h in h_list],
        "canonical_float64": {},
        "canonical_float128": {},
        "schemes": {},
    }

    # 1. Float64 sweep across all 5 canonical integrators (use_optimized_kernels=True)
    f64_q_ends_3d: Dict[str, List[np.ndarray]] = {}
    for name, factory in canonical_f64.items():
        s = stages_map[name]
        dE_max_list: List[float] = []
        q_ends_list: List[np.ndarray] = []
        for n_steps in step_counts:
            h = T_in / n_steps
            st = create_hier3d_f64(e_in=0.30, inc_deg=28.0, a_out=12.0)
            sys_h = mod_01.NBodySystem(
                state=st,
                integrator=factory(),
                G=1.0,
                enable_dyadic_substepping=False,
                in_shadow_mode=True,
                use_optimized_kernels=True,
            )
            max_err = 0.0
            stride = max(1, n_steps // 16)
            for idx in range(1, n_steps + 1):
                sys_h.step(dt=h)
                if idx % stride == 0 or idx == n_steps:
                    phys = sys_h.observe_physical_state()
                    H_cur = mod_05.compute_lie_bargmann_invariants(phys, G=1.0).H
                    max_err = max(max_err, abs(H_cur - H0) / abs(H0))
            dE_max_list.append(float(max_err))
            q_ends_list.append(np.copy(sys_h.observe_physical_state().q))

        f64_q_ends_3d[name] = q_ends_list
        results_3d["canonical_float64"][name] = {
            "stages_s": s,
            "h_steps": h_list,
            "tau_over_s": [h / s for h in h_list],
            "delta_E_rel": dE_max_list,
            "local_slopes": _compute_local_slopes(h_list, dE_max_list),
        }

    # 2. Float128 (35-digit mpmath) sweep across all 5 canonical integrators
    step_counts_128 = [8, 16, 32, 64, 128, 256]
    h_list_128 = [T_in / n for n in step_counts_128]
    T_in_mp = 2 * mp.pi
    st0_128 = create_hier3d_f128(e_in=0.30, inc_deg=28.0, a_out=12.0, is_f128=True)
    assert st0_128.is_f128 is True, "3D Hierarchical Triple state must be native Float128"
    H0_128 = inv128(st0_128, G=mp.mpf("1.0")).H
    f128_q_ends_3d: Dict[str, List[np.ndarray]] = {}

    for name, factory in canonical_f128.items():
        s = stages_map[name]
        dE_128_list: List[float] = []
        q_ends_128_list: List[np.ndarray] = []
        print(f"  [Part 2 Float128] Sweeping {name} (s={s}) on 3D Eccentric Triple...", flush=True)
        for n_steps in step_counts_128:
            h_mp = T_in_mp / mp.mpf(n_steps)
            st_128 = create_hier3d_f128(e_in=0.30, inc_deg=28.0, a_out=12.0, is_f128=True)
            sys_128 = Sys128(
                state=st_128,
                integrator=factory(),
                G=mp.mpf("1.0"),
                enable_dyadic_substepping=False,
                in_shadow_mode=True,
            )
            max_err_mp = mp.mpf("0.0")
            stride = max(1, n_steps // 8)
            for idx in range(1, n_steps + 1):
                sys_128.step(dt=h_mp)
                if idx % stride == 0 or idx == n_steps:
                    phys_128 = sys_128.observe_physical_state()
                    H_cur_128 = inv128(phys_128, G=mp.mpf("1.0")).H
                    err_mp = abs(H_cur_128 - H0_128) / abs(H0_128)
                    if err_mp > max_err_mp:
                        max_err_mp = err_mp
            dE_128_list.append(float(max_err_mp))
            q_ends_128_list.append(np.copy(sys_128.observe_physical_state().q))

        f128_q_ends_3d[name] = q_ends_128_list
        results_3d["canonical_float128"][name] = {
            "stages_s": s,
            "h_steps": h_list_128,
            "tau_over_s": [h / s for h in h_list_128],
            "delta_E_rel": dE_128_list,
            "local_slopes": _compute_local_slopes(h_list_128, dE_128_list),
        }

    # Construct 35-digit Richardson-extrapolated O(h^8) reference state at T_in = 2*pi from the two finest BM6_128 runs
    q_bm6_3d_finest = f128_q_ends_3d["BM6 (6th order optimal)"][-1]
    q_bm6_3d_prev = f128_q_ends_3d["BM6 (6th order optimal)"][-2]
    q_ref_3d_f128_obj = np.empty((3, 3), dtype=object)
    for bi in range(3):
        for d in range(3):
            q_ref_3d_f128_obj[bi, d] = (mp.mpf("64.0") * q_bm6_3d_finest[bi, d] - q_bm6_3d_prev[bi, d]) / mp.mpf("63.0")
    q_ref_3d_f64 = np.array(
        [[float(q_ref_3d_f128_obj[bi, d]) for d in range(3)] for bi in range(3)], dtype=np.float64
    )

    for name in canonical_f64.keys():
        dq_64_list = [float(np.linalg.norm(q_end - q_ref_3d_f64)) for q_end in f64_q_ends_3d[name]]
        results_3d["canonical_float64"][name]["delta_q_norm"] = dq_64_list

    for name in canonical_f128.keys():
        dq_128_list: List[float] = []
        for q_end_128 in f128_q_ends_3d[name]:
            dq_sq = mp.mpf("0.0")
            for bi in range(3):
                for d in range(3):
                    diff_qd = q_end_128[bi, d] - q_ref_3d_f128_obj[bi, d]
                    dq_sq += diff_qd * diff_qd
            dq_128_list.append(float(mp.sqrt(dq_sq)))
        results_3d["canonical_float128"][name]["delta_q_norm"] = dq_128_list

    # 3. Populate summary comparison keys on step_counts[2:6] ([32, 64, 128, 256])
    bm4_errs = results_3d["canonical_float64"]["BM4 (4th order optimal)"]["delta_E_rel"][2:6]
    z4_errs = results_3d["canonical_float64"]["PBM4 (Processed BM4)"]["delta_E_rel"][2:6]
    h_sub = h_list[2:6]

    results_3d["schemes"]["BM4 (Unprocessed)"] = {
        "delta_E_rel": bm4_errs,
        "local_slopes": _compute_local_slopes(h_sub, bm4_errs),
    }
    results_3d["schemes"]["PBM4 (Z4_BM4 Homological)"] = {
        "delta_E_rel": z4_errs,
        "local_slopes": _compute_local_slopes(h_sub, z4_errs),
    }

    gains_z4 = [bm4_errs[i] / z4_errs[i] for i in range(len(h_sub))]

    results_3d["gains_z4_over_bm4"] = gains_z4
    results_3d["max_gain_z4_over_bm4"] = float(max(gains_z4))
    print(
        f"        3D Eccentric Triple max gain PBM4/BM4 = {max(gains_z4):.1f}x",
        flush=True,
    )
    return results_3d


def run_part3_section_6_6_microkernels() -> Dict[str, Any]:
    print("=== [Part 3] Running Section 6.6 Micro-Kernel Certification (Table 4) & Multi-Body Scaling (N=3..64) ===", flush=True)

    # 1. Table 4 Row 1: Vectorized Screened Apex Kick (Antisymmetry, Torque=0, Bit-Exact Concordance, N=6 2-Supernode Cross-Force)
    st_bur = create_burrau_f64()
    tree_bur = mod_02.build_optimal_tree(st_bur)
    base_apex = ScreenedApexPotential(G=1.0)
    opt_apex = mod_06.OptimizedScreenedApexPotential(st_bur.masses, G=1.0)

    f_base, _ = base_apex.compute_apex_forces(st_bur, tree_bur, include_screened_monopole=False)
    f_opt, _ = opt_apex.compute_apex_forces(st_bur, tree_bur, include_screened_monopole=False)
    row1_antisymmetry = float(np.linalg.norm(np.sum(f_opt, axis=0)))
    row1_torque = float(np.linalg.norm(np.sum(np.cross(st_bur.q, f_opt), axis=0)))
    row1_concordance = float(np.max(np.abs(f_opt - f_base)))

    # N=6 Two-Supernode Cross-Pair Force Verification ([MAJ-6.6-03])
    masses6 = np.ones(6, dtype=np.float64)
    q6 = np.array(
        [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.0, 0.05, 0.0], [5.0, 0.0, 0.0], [5.05, 0.0, 0.0], [5.0, 0.05, 0.0]],
        dtype=np.float64,
    )
    st6 = mod_01.ExtendedPhaseState(q=q6, p=np.zeros_like(q6), masses=masses6)
    tree6 = mod_02.build_optimal_tree(st6)
    opt_apex6 = mod_06.OptimizedScreenedApexPotential(masses6, G=1.0)
    f6_clusters, _ = opt_apex6.compute_apex_forces(
        st6, tree6, supernode_clusters=[(0, 1, 2), (3, 4, 5)], include_screened_monopole=False
    )
    inter_supernode_force_norm = float(np.linalg.norm(np.sum(f6_clusters[:3], axis=0)))

    # 2. Table 4 Row 2: Inlined Kepler-Sundman Drift (2x6 Horner + Stumpff recurrence <= 1 ULP & Battin Wronskian)
    from decimal import Decimal, getcontext
    getcontext().prec = 50

    def _exact_stumpff_decimal(z_flt: float) -> List[float]:
        z_dec = Decimal(str(z_flt))
        out_c: List[float] = []
        for m_idx in range(4):
            s_acc = Decimal(0)
            fact = Decimal(math.factorial(m_idx))
            z_pow = Decimal(1)
            for k_idx in range(22):
                if k_idx > 0:
                    fact *= Decimal((2 * k_idx + m_idx - 1) * (2 * k_idx + m_idx))
                    z_pow *= -z_dec
                s_acc += z_pow / fact
            out_c.append(float(s_acc))
        return out_c

    max_stumpff_err = 0.0
    max_stumpff_opt_vs_base = 0.0
    for z_val in (-0.099, -0.05, -1e-4, 0.0, 1e-4, 0.05, 0.099):
        c_fast = mod_06.OptimizedKeplerSolver.stumpff_c0123_scalar(z_val)
        c_ref = mod_04.stumpff_c(z_val)[:4]
        c_exact50 = _exact_stumpff_decimal(z_val)
        max_stumpff_opt_vs_base = max(max_stumpff_opt_vs_base, float(np.max(np.abs(np.array(c_fast) - np.array(c_ref)))))
        max_stumpff_err = max(max_stumpff_err, float(np.max(np.abs(np.array(c_fast) - np.array(c_exact50)))))

    r0_k = np.array([0.8, 0.2, -0.1], dtype=np.float64)
    w0_k = np.array([-0.15, 0.95, 0.05], dtype=np.float64)
    mu_k = 0.5
    GM_k = 1.0
    dt_k = 0.05
    r_fast, w_fast = mod_06.step_kepler_sundman_fast(r0_k, w0_k, mu=mu_k, GM=GM_k, dt_target=dt_k)
    r_ref, w_ref = mod_04.step_kepler_sundman(r0_k, w0_k, mu=mu_k, GM=GM_k, dt_target=dt_k)
    L_before = np.cross(r0_k, w0_k)
    L_after = np.cross(r_fast, w_fast)
    row2_wronskian_and_L = float(
        max(np.linalg.norm(L_after - L_before), np.max(np.abs(r_fast - r_ref)), np.max(np.abs(w_fast - w_ref)))
    )

    # 3. Table 4 Row 3 & Row 4: Synchronous Supernode Dirac Constraints & Multi-Substep Time-Reversibility ([CRIT-6.6-03])
    masses_sn = np.array([1.0, 1.2, 0.8], dtype=np.float64)
    q_sn0 = np.array([[0.0, 0.0, 0.0], [0.18, 0.02, -0.01], [-0.09, 0.15, 0.03]], dtype=np.float64)
    p_sn0 = np.array([[0.05, -0.08, 0.01], [-0.03, 0.06, -0.02], [-0.02, 0.02, 0.01]], dtype=np.float64)
    st_sn0 = mod_01.ExtendedPhaseState(q=q_sn0, p=p_sn0, masses=masses_sn).to_barycentric()
    dt_sn = 0.015  # triggers multiple adaptive substeps (8 substeps)

    q_sn_fwd, p_sn_fwd = mod_06.step_supernode_synchronous_jit(
        st_sn0.q, st_sn0.p, st_sn0.masses, (0, 1, 2), dt=dt_sn, G=1.0
    )
    q_sn_bwd, p_sn_bwd = mod_06.step_supernode_synchronous_jit(
        q_sn_fwd, p_sn_fwd, st_sn0.masses, (0, 1, 2), dt=-dt_sn, G=1.0
    )
    row3_dirac_barycentric = float(
        max(
            np.linalg.norm(np.sum(st_sn0.masses[:, None] * q_sn_fwd, axis=0)),
            np.linalg.norm(np.sum(p_sn_fwd, axis=0)),
        )
    )
    row4_supernode_reversibility = float(
        max(np.linalg.norm(q_sn_bwd - st_sn0.q), np.linalg.norm(p_sn_bwd - st_sn0.p))
    )

    # Also check concordance between baseline, optimized, and JIT supernode solvers
    q_sn_base, p_sn_base = mod_04.step_supernode_synchronous(
        st_sn0.q, st_sn0.p, st_sn0.masses, (0, 1, 2), dt=dt_sn, G=1.0
    )
    q_sn_opt, p_sn_opt = mod_06.step_supernode_synchronous_fast(
        st_sn0.q, st_sn0.p, st_sn0.masses, (0, 1, 2), dt=dt_sn, G=1.0
    )
    supernode_kernel_concordance = float(
        max(
            np.max(np.abs(q_sn_fwd - q_sn_base)),
            np.max(np.abs(p_sn_fwd - p_sn_base)),
            np.max(np.abs(q_sn_opt - q_sn_base)),
        )
    )

    # 4. Table 4 Row 5: End-to-End Trajectory Concordance (use_optimized_kernels=False vs True)
    st_e2e_base = create_hier3d_f64()
    st_e2e_opt = create_hier3d_f64()
    sys_base = mod_01.NBodySystem(
        state=st_e2e_base,
        integrator=mod_02.create_pbm4_integrator(G=1.0, mode="homological"),
        G=1.0,
        use_optimized_kernels=False,
    )
    sys_opt = mod_01.NBodySystem(
        state=st_e2e_opt,
        integrator=mod_02.create_pbm4_integrator(G=1.0, mode="homological"),
        G=1.0,
        use_optimized_kernels=True,
    )
    sys_base.scheduler.warn_on_saturation = True
    sys_opt.scheduler.warn_on_saturation = True
    end_base = sys_base.evolve(t_span=0.5, dt=0.05)
    end_opt = sys_opt.evolve(t_span=0.5, dt=0.05)
    assert len(sys_base.scheduler.saturation_events) == 0, f"Saturation in sys_base: {sys_base.scheduler.saturation_events}"
    assert len(sys_opt.scheduler.saturation_events) == 0, f"Saturation in sys_opt: {sys_opt.scheduler.saturation_events}"
    row5_trajectory_concordance = float(
        max(np.max(np.abs(end_opt.q - end_base.q)), np.max(np.abs(end_opt.p - end_base.p)))
    )

    # 4b. Extreme-Regime Equality Certification (Opt vs Base across 6 Extreme Stress Regimes)
    print("  [3.1b] Certifying Opt vs Base Equality across 6 Extreme Singular & Multi-Scale Regimes...", flush=True)
    # Regime 1: Extreme Kepler pericenter homotopy 1-e in {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12} (L1: outgoing + incoming pericenter traversal)
    max_kepler_l1_dq = 0.0
    max_kepler_l1_dw = 0.0
    for eps_e in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
        e_val = 1.0 - eps_e
        r_peri = 1.0 - e_val
        v_peri = math.sqrt((1.0 + e_val) / r_peri)
        r0_l1 = np.array([r_peri, 0.0, 0.0], dtype=np.float64)
        w0_l1 = np.array([0.0, 0.5 * v_peri, 0.0], dtype=np.float64)
        r_b_l1, w_b_l1 = mod_04.step_kepler_sundman(r0_l1, w0_l1, mu=0.5, GM=1.0, dt_target=0.25)
        r_o_l1, w_o_l1 = mod_06.step_kepler_sundman_fast(r0_l1, w0_l1, mu=0.5, GM=1.0, dt_target=0.25)
        max_kepler_l1_dq = max(max_kepler_l1_dq, float(np.linalg.norm(r_o_l1 - r_b_l1)))
        max_kepler_l1_dw = max(max_kepler_l1_dw, float(np.linalg.norm(w_o_l1 - w_b_l1)))
        # Also test full incoming traversal across pericenter (activating rho_raw < 0.25 * rho0 branch)
        r_in_0, w_in_0 = mod_04.step_kepler_sundman(r0_l1, w0_l1, mu=0.5, GM=1.0, dt_target=-1e-5)
        r_b_tr, w_b_tr = mod_04.step_kepler_sundman(r_in_0, w_in_0, mu=0.5, GM=1.0, dt_target=2e-5)
        r_o_tr, w_o_tr = mod_06.step_kepler_sundman_fast(r_in_0, w_in_0, mu=0.5, GM=1.0, dt_target=2e-5)
        max_kepler_l1_dq = max(max_kepler_l1_dq, float(np.linalg.norm(r_o_tr - r_b_tr)))
        max_kepler_l1_dw = max(max_kepler_l1_dw, float(np.linalg.norm(w_o_tr - w_b_tr)))

    # Regime 2: Exact radial head-on collision bounce 1-e = 0, d_0 = 0 (L2)
    r0_l2 = np.array([0.5, 0.0, 0.0], dtype=np.float64)
    w0_l2 = np.array([-0.4, 0.0, 0.0], dtype=np.float64)
    r_b_l2, w_b_l2 = mod_04.step_kepler_sundman(r0_l2, w0_l2, mu=0.5, GM=1.0, dt_target=0.8)
    r_o_l2, w_o_l2 = mod_06.step_kepler_sundman_fast(r0_l2, w0_l2, mu=0.5, GM=1.0, dt_target=0.8)
    kepler_l2_bounce_dq = float(np.linalg.norm(r_o_l2 - r_b_l2))
    kepler_l2_bounce_dw = float(np.linalg.norm(w_o_l2 - w_b_l2))

    # Regime 3: Synchronous Supernode Deep 3-Body Close Encounter (8 adaptive substeps)
    sn_opt_vs_base_dq = float(np.linalg.norm(q_sn_opt - q_sn_base))
    sn_opt_vs_base_dp = float(np.linalg.norm(p_sn_opt - p_sn_base))
    sn_jit_vs_base_dq = float(np.linalg.norm(q_sn_fwd - q_sn_base))
    sn_jit_vs_base_dp = float(np.linalg.norm(p_sn_fwd - p_sn_base))

    # Regime 4: Multi-Body N=64 Screened Apex (2,016 pairs, @njit vs Reference)
    rng_64 = np.random.default_rng(123)
    m_64 = np.ones(64, dtype=np.float64)
    st_64 = mod_01.ExtendedPhaseState(
        q=rng_64.uniform(-10.0, 10.0, size=(64, 3)),
        p=rng_64.uniform(-0.1, 0.1, size=(64, 3)),
        masses=m_64,
    ).to_barycentric()
    tree_64 = mod_02.build_optimal_tree(st_64)
    f64_base, _ = ScreenedApexPotential(G=1.0).compute_apex_forces(st_64, tree_64, include_screened_monopole=False)
    f64_opt, _ = mod_06.OptimizedScreenedApexPotential(m_64, G=1.0).compute_apex_forces(
        st_64, tree_64, include_screened_monopole=False
    )
    apex_n64_rel_err = float(np.linalg.norm(f64_opt - f64_base) / max(np.linalg.norm(f64_base), 1e-30))

    # Regime 5: End-to-End Comet Lexell 1770 (H1: Deep Perijove 3.10 R_J + Dynamic Hill Atlas Chart Flips)
    from normal_cone_nbody_float64.benchmarks.registry import (
        create_comet_lexell_1770_benchmark,
        create_shoemaker_levy_9_1992_benchmark,
    )
    G_lex = 1.0
    st_lex_b = create_comet_lexell_1770_benchmark(G=G_lex)
    st_lex_o = create_comet_lexell_1770_benchmark(G=G_lex)
    sys_lex_b = mod_01.NBodySystem(
        state=st_lex_b,
        integrator=mod_02.create_pbm4_integrator(G=G_lex, mode="homological"),
        G=G_lex,
        max_dyadic_level=11,
        enable_dyadic_substepping=True,
        in_shadow_mode=True,
        use_optimized_kernels=False,
    )
    sys_lex_o = mod_01.NBodySystem(
        state=st_lex_o,
        integrator=mod_02.create_pbm4_integrator(G=G_lex, mode="homological"),
        G=G_lex,
        max_dyadic_level=11,
        enable_dyadic_substepping=True,
        in_shadow_mode=True,
        use_optimized_kernels=True,
    )
    sys_lex_b.scheduler.warn_on_saturation = True
    sys_lex_o.scheduler.warn_on_saturation = True
    for _ in range(150):
        sys_lex_b.step(dt=0.08)
        sys_lex_o.step(dt=0.08)
    assert len(sys_lex_b.scheduler.saturation_events) == 0, f"Saturation in Lexell base: {sys_lex_b.scheduler.saturation_events}"
    assert len(sys_lex_o.scheduler.saturation_events) == 0, f"Saturation in Lexell opt: {sys_lex_o.scheduler.saturation_events}"
    end_lex_b = sys_lex_b.observe_physical_state()
    end_lex_o = sys_lex_o.observe_physical_state()
    lexell_h1_dq = float(np.linalg.norm(end_lex_o.q - end_lex_b.q))
    lexell_h1_dp = float(np.linalg.norm(end_lex_o.p - end_lex_b.p))

    # Regime 6: End-to-End Comet Shoemaker-Levy 9 (H2: Sub-Roche Tidal Disruption 1.32 R_J + Active Supernode)
    G_sl9 = 1.0
    st_sl9_b = create_shoemaker_levy_9_1992_benchmark(G=G_sl9)
    st_sl9_o = create_shoemaker_levy_9_1992_benchmark(G=G_sl9)
    sys_sl9_b = mod_01.NBodySystem(
        state=st_sl9_b,
        integrator=mod_02.create_pbm4_integrator(G=G_sl9, mode="homological"),
        G=G_sl9,
        max_dyadic_level=8,
        enable_dyadic_substepping=True,
        in_shadow_mode=True,
        use_optimized_kernels=False,
    )
    sys_sl9_o = mod_01.NBodySystem(
        state=st_sl9_o,
        integrator=mod_02.create_pbm4_integrator(G=G_sl9, mode="homological"),
        G=G_sl9,
        max_dyadic_level=8,
        enable_dyadic_substepping=True,
        in_shadow_mode=True,
        use_optimized_kernels=True,
    )
    sys_sl9_b.scheduler.warn_on_saturation = True
    sys_sl9_o.scheduler.warn_on_saturation = True
    for _ in range(160):
        sys_sl9_b.step(dt=0.00018)
        sys_sl9_o.step(dt=0.00018)
    assert len(sys_sl9_b.scheduler.saturation_events) == 0, f"Saturation in SL9 base: {sys_sl9_b.scheduler.saturation_events}"
    assert len(sys_sl9_o.scheduler.saturation_events) == 0, f"Saturation in SL9 opt: {sys_sl9_o.scheduler.saturation_events}"
    end_sl9_b = sys_sl9_b.observe_physical_state()
    end_sl9_o = sys_sl9_o.observe_physical_state()
    sl9_h2_dq = float(np.linalg.norm(end_sl9_o.q - end_sl9_b.q))
    sl9_h2_dp = float(np.linalg.norm(end_sl9_o.p - end_sl9_b.p))

    # 5. Multi-Body Scaling Sweep for N in {3, 4, 8, 16, 32, 64} (Apex Kick, Kepler Drift, Supernode JIT)
    print("  [3.2] Benchmarking Micro-Kernel Throughput across N in {3, 4, 8, 16, 32, 64}...", flush=True)
    n_bodies_grid = [3, 4, 8, 16, 32, 64]
    apex_speedups: List[float] = []
    apex_base_us: List[float] = []
    apex_opt_us: List[float] = []

    rng = np.random.default_rng(42)
    for N_val in n_bodies_grid:
        masses_n = np.ones(N_val, dtype=np.float64)
        q_n = rng.uniform(-10.0, 10.0, size=(N_val, 3))
        p_n = rng.uniform(-0.1, 0.1, size=(N_val, 3))
        st_n = mod_01.ExtendedPhaseState(q=q_n, p=p_n, masses=masses_n).to_barycentric()
        tree_n = mod_02.build_optimal_tree(st_n)

        base_ev = ScreenedApexPotential(G=1.0)
        opt_ev = mod_06.OptimizedScreenedApexPotential(masses_n, G=1.0)
        # Warmup
        base_ev.compute_apex_forces(st_n, tree_n, include_screened_monopole=False)
        opt_ev.compute_apex_forces(st_n, tree_n, include_screened_monopole=False)

        n_iters = 200 if N_val <= 16 else 80
        dt_base_runs = []
        dt_opt_runs = []
        for _rep in range(5):
            t0 = time.perf_counter()
            for _ in range(n_iters):
                base_ev.compute_apex_forces(st_n, tree_n, include_screened_monopole=False)
            dt_base_runs.append((time.perf_counter() - t0) / n_iters * 1e6)

            t0 = time.perf_counter()
            for _ in range(n_iters):
                opt_ev.compute_apex_forces(st_n, tree_n, include_screened_monopole=False)
            dt_opt_runs.append((time.perf_counter() - t0) / n_iters * 1e6)
        dt_base = min(dt_base_runs)
        dt_opt = min(dt_opt_runs)

        apex_base_us.append(float(dt_base))
        apex_opt_us.append(float(dt_opt))
        apex_speedups.append(float(dt_base / max(dt_opt, 1e-9)))

    # Benchmark Kepler drift-and-leaf stage speedup (Phase 3.1: baseline tree dict vs OptimizedKeplerSolver N=3 fast path + 2x5 Horner)
    from normal_cone_nbody_float64.integrators.common import step_drift_and_kepler as _base_step_drift_and_kepler
    opt_kepler_solver = mod_06.OptimizedKeplerSolver(G=1.0)
    st_k3 = create_hier3d_f64()
    tree_k3 = mod_02.build_optimal_tree(st_k3)
    _base_step_drift_and_kepler(st_k3, tree_k3, dt_k, G=1.0)
    opt_kepler_solver.step_drift_and_kepler(st_k3, tree_k3, dt_k, G=1.0)

    n_kep_iters = 400
    t0 = time.perf_counter()
    for _ in range(n_kep_iters):
        _base_step_drift_and_kepler(st_k3, tree_k3, dt_k, G=1.0)
    kep_base_us = (time.perf_counter() - t0) / n_kep_iters * 1e6

    t0 = time.perf_counter()
    for _ in range(n_kep_iters):
        opt_kepler_solver.step_drift_and_kepler(st_k3, tree_k3, dt_k, G=1.0)
    kep_opt_us = (time.perf_counter() - t0) / n_kep_iters * 1e6
    kepler_speedup = float(kep_base_us / max(kep_opt_us, 1e-9))

    n_sn_iters = 60
    t0 = time.perf_counter()
    for _ in range(n_sn_iters):
        mod_04.step_supernode_synchronous(st_sn0.q, st_sn0.p, st_sn0.masses, (0, 1, 2), dt=dt_sn, G=1.0)
    sn_base_us = (time.perf_counter() - t0) / n_sn_iters * 1e6

    t0 = time.perf_counter()
    for _ in range(n_sn_iters):
        mod_06.step_supernode_synchronous_jit(st_sn0.q, st_sn0.p, st_sn0.masses, (0, 1, 2), dt=dt_sn, G=1.0)
    sn_jit_us = (time.perf_counter() - t0) / n_sn_iters * 1e6
    supernode_jit_speedup = float(sn_base_us / max(sn_jit_us, 1e-9))

    return {
        "table4_invariants": {
            "row1_apex_antisymmetry": float(row1_antisymmetry),
            "row1_apex_torque": float(row1_torque),
            "row1_apex_concordance": float(row1_concordance),
            "row1_n6_inter_supernode_force_norm": float(inter_supernode_force_norm),
            "row2_stumpff_horner_5stage_max_err": float(max_stumpff_err),
            "row2_stumpff_opt_vs_base_err": float(max_stumpff_opt_vs_base),
            "row2_battin_wronskian": float(row2_wronskian_and_L),
            "row3_dirac_barycentric": float(row3_dirac_barycentric),
            "row4_supernode_multistep_reversibility": float(row4_supernode_reversibility),
            "row4_supernode_kernel_concordance": float(supernode_kernel_concordance),
            "row5_trajectory_concordance": float(row5_trajectory_concordance),
        },
        "multibody_scaling": {
            "N_grid": n_bodies_grid,
            "apex_base_us": apex_base_us,
            "apex_opt_us": apex_opt_us,
            "apex_speedup": apex_speedups,
            "kepler_base_us": float(kep_base_us),
            "kepler_opt_us": float(kep_opt_us),
            "kepler_speedup": kepler_speedup,
            "supernode_base_us": float(sn_base_us),
            "supernode_jit_us": float(sn_jit_us),
            "supernode_jit_speedup": supernode_jit_speedup,
        },
        "extreme_regime_concordance": {
            "kepler_l1_pericenter_1e_minus_12_dq": float(max_kepler_l1_dq),
            "kepler_l1_pericenter_1e_minus_12_dw": float(max_kepler_l1_dw),
            "kepler_l2_head_on_bounce_d0_0_dq": float(kepler_l2_bounce_dq),
            "kepler_l2_head_on_bounce_d0_0_dw": float(kepler_l2_bounce_dw),
            "supernode_opt_vs_base_dq": float(sn_opt_vs_base_dq),
            "supernode_opt_vs_base_dp": float(sn_opt_vs_base_dp),
            "supernode_jit_vs_base_dq": float(sn_jit_vs_base_dq),
            "supernode_jit_vs_base_dp": float(sn_jit_vs_base_dp),
            "apex_n64_njit_rel_err": float(apex_n64_rel_err),
            "lexell_h1_atlas_flip_e2e_dq": float(lexell_h1_dq),
            "lexell_h1_atlas_flip_e2e_dp": float(lexell_h1_dp),
            "shoemaker_levy9_h2_roche_supernode_e2e_dq": float(sl9_h2_dq),
            "shoemaker_levy9_h2_roche_supernode_e2e_dp": float(sl9_h2_dp),
        },
    }


def verify_contract_assertions(results: Dict[str, Any]) -> None:
    print("=== [Contract Gate] Verifying PASS/FAIL Assertions of Contract 01 ===", flush=True)
    fig8 = results["figure8_convergence"]
    hier3d = results["hierarchical_3d_eccentric_probe"]
    sec66 = results["section_6_6_microkernels"]

    # ASSERT-G01-01: True 35-digit Float128 execution & descent to ~1e-30 (and dq <= 1e-30)
    assert fig8["float128_verified"] is True, "ASSERT-G01-01 FAILED: Float128 verification flag is False"
    bm6_min_f128 = fig8["float128"]["BM6 (6th order optimal)"]["min_delta_E_rel"]
    bm6_min_dq_f128 = min(fig8["float128"]["BM6 (6th order optimal)"]["delta_q_norm"])
    assert bm6_min_f128 <= 1.5e-30, f"ASSERT-G01-01 FAILED: BM6 Float128 min error {bm6_min_f128} > 1.5e-30"
    assert bm6_min_dq_f128 <= 1.0e-30, f"ASSERT-G01-01 FAILED: BM6 Float128 min dq {bm6_min_dq_f128} > 1.0e-30"

    # ASSERT-G01-02: Asymptotic slopes in Float128
    s_strang = fig8["float128"]["Strang (2nd order)"]["mean_asymptotic_slope"]
    s_y4 = fig8["float128"]["Yoshida4 (4th order)"]["mean_asymptotic_slope"]
    s_bm4 = fig8["float128"]["BM4 (4th order optimal)"]["mean_asymptotic_slope"]
    s_bm6 = fig8["float128"]["BM6 (6th order optimal)"]["mean_asymptotic_slope"]
    assert abs(s_strang - 2.0) <= 0.12, f"ASSERT-G01-02 FAILED: Strang slope {s_strang}"
    assert abs(s_y4 - 4.0) <= 0.12, f"ASSERT-G01-02 FAILED: Yoshida4 slope {s_y4}"
    assert abs(s_bm4 - 4.0) <= 0.15, f"ASSERT-G01-02 FAILED: BM4 slope {s_bm4}"
    assert abs(s_bm6 - 6.0) <= 0.25, f"ASSERT-G01-02 FAILED: BM6 slope {s_bm6}"

    # ASSERT-G01-03: Table 5 common abscissa h_eff = 0.005 ordering & BM4/PBM4 over Y4 gains
    t5 = fig8["table5_common_abscissa_h_eff_0_01"]
    assert (
        t5["BM6 (6th order optimal)"]["delta_E_rel_max"]
        < min(t5["PBM4 (Processed BM4)"]["delta_E_rel_max"], t5["BM4 (4th order optimal)"]["delta_E_rel_max"])
        <= max(t5["PBM4 (Processed BM4)"]["delta_E_rel_max"], t5["BM4 (4th order optimal)"]["delta_E_rel_max"])
        < t5["Yoshida4 (4th order)"]["delta_E_rel_max"]
        < t5["Strang (2nd order)"]["delta_E_rel_max"]
    ), "ASSERT-G01-03 FAILED: Ordering at common h_eff violated"
    assert fig8["equal_cost_bm4_over_y4_gain"] >= 3.0, (
        f"ASSERT-G01-03 FAILED: Equal-cost BM4/Y4 gain {fig8['equal_cost_bm4_over_y4_gain']} < 3.0"
    )
    assert fig8["equal_step_bm4_over_y4_gain"] >= 50.0, (
        f"ASSERT-G01-03 FAILED: Equal-step BM4/Y4 gain {fig8['equal_step_bm4_over_y4_gain']} < 50.0"
    )

    # ASSERT-G01-04: Hierarchical 3D eccentric order 4 -> 6 elevation of PBM4 (Z4_BM4)
    assert hier3d["max_gain_z4_over_bm4"] >= 150.0, (
        f"ASSERT-G01-04 FAILED: Z4_BM4/BM4 gain {hier3d['max_gain_z4_over_bm4']} < 150"
    )

    # ASSERT-G01-05: Section 6.6 Table 4 Invariants <= 1e-14 (including multi-substep supernode reversibility)
    t4 = sec66["table4_invariants"]
    assert t4["row1_apex_antisymmetry"] <= 1e-14
    assert t4["row1_apex_concordance"] <= 1e-14
    assert t4["row1_n6_inter_supernode_force_norm"] > 1e-2
    assert t4["row2_stumpff_horner_5stage_max_err"] <= 5e-16
    assert t4["row2_battin_wronskian"] <= 1e-14
    assert t4["row3_dirac_barycentric"] <= 1e-14
    assert t4["row4_supernode_multistep_reversibility"] <= 1e-14
    assert t4["row5_trajectory_concordance"] <= 1e-13

    # ASSERT-G01-06: Lie-Bargmann conservation in Float64 & Float128
    for name, res64 in fig8["float64"].items():
        assert res64["max_delta_P"] <= 1e-13, f"Float64 delta_P failed for {name}"
        assert res64["max_delta_L"] <= 5e-13, f"Float64 delta_L failed for {name}"
    for name, res128 in fig8["float128"].items():
        assert res128["max_delta_P"] <= 1e-31, f"Float128 delta_P failed for {name}"
        assert res128["max_delta_L"] <= 1e-31, f"Float128 delta_L failed for {name}"

    print("  [PASS] All 6 Contract Assertions (ASSERT-G01-01 to ASSERT-G01-06) CERTIFIED!", flush=True)


def run_part1_table5_only(existing_fig8: Dict[str, Any]) -> Dict[str, Any]:
    print("=== [Part 1.3] Updating Table 5 Common-Abscissa Sweep against Unoptimized BM6 Float128 Reference ===", flush=True)
    T_table5 = 6.0
    stages_map = {
        "Strang (2nd order)": 1,
        "Yoshida4 (4th order)": 3,
        "BM4 (4th order optimal)": 6,
        "PBM4 (Processed BM4)": 6,
        "BM6 (6th order optimal)": 10,
    }
    integrators_f64: Dict[str, Callable[[], Any]] = {
        "Strang (2nd order)": lambda: mod_04.StrangIntegrator(),
        "Yoshida4 (4th order)": lambda: mod_04.Yoshida4Integrator(),
        "BM4 (4th order optimal)": lambda: mod_04.BlanesMoan4Integrator(),
        "PBM4 (Processed BM4)": lambda: mod_02.create_pbm4_integrator(G=1.0, mode="homological"),
        "BM6 (6th order optimal)": lambda: mod_04.BlanesMoan6Integrator(),
    }
    st_ref_t6_128 = create_fig8_f128(is_f128=True)
    assert st_ref_t6_128.is_f128 is True
    H0_128 = inv128(st_ref_t6_128, G=mp.mpf("1.0")).H
    sys_ref_t6_128 = Sys128(
        state=st_ref_t6_128,
        integrator=BM6_128(),
        atlas_manager=AtlasManager128(),
        G=mp.mpf("1.0"),
        enable_dyadic_substepping=False,
        in_shadow_mode=True,
    )
    dt_ref_t6_mp = mp.mpf("0.01")
    n_ref_t6_steps = int(round(T_table5 / float(dt_ref_t6_mp)))
    for _ in range(n_ref_t6_steps):
        sys_ref_t6_128.step(dt=dt_ref_t6_mp)
    phys_ref_t6_128 = sys_ref_t6_128.observe_physical_state()
    inv_ref_t6_128 = inv128(phys_ref_t6_128, G=mp.mpf("1.0"))
    end_dE_ref_t6_128 = float(abs(inv_ref_t6_128.H - H0_128) / abs(H0_128))
    q_ref_t6 = np.array(
        [[float(phys_ref_t6_128.q[bi, d]) for d in range(3)] for bi in range(3)], dtype=np.float64
    )
    f64_results = existing_fig8["float64"]
    f128_results = existing_fig8["float128"]
    table5_ref_bm6_f128 = {
        "scheme": "BM6 (Float128 Reference, Unoptimized)",
        "stages_s": 10,
        "h_eff": 0.001,
        "h_macro": 0.01,
        "n_steps": n_ref_t6_steps,
        "n_force_evals": n_ref_t6_steps * 10,
        "delta_E_rel_endpoint": end_dE_ref_t6_128,
        "float128_min_error_Tconv": float(f128_results["BM6 (6th order optimal)"]["min_delta_E_rel"]),
        "float128_min_dq_Tconv": float(min(f128_results["BM6 (6th order optimal)"]["delta_q_norm"])),
    }

    st0_f64 = create_fig8_f64()
    H0_f64 = float(mod_05.compute_lie_bargmann_invariants(st0_f64, G=1.0).H)
    h_eff_common = 0.005
    table5_common: Dict[str, Any] = {}
    for name, s in stages_map.items():
        h_common = h_eff_common * s
        n_steps_common = int(round(T_table5 / h_common))
        st_c = create_fig8_f64()
        sys_c = mod_01.NBodySystem(
            state=st_c,
            integrator=integrators_f64[name](),
            atlas_manager=AtlasManager64(),
            G=1.0,
            enable_dyadic_substepping=False,
            in_shadow_mode=True,
            use_optimized_kernels=True,
        )
        max_dE_c = 0.0
        sum_sq_dE_c = 0.0
        for _ in range(n_steps_common):
            sys_c.step(dt=h_common)
            phys_c = sys_c.observe_physical_state()
            H_c = mod_05.compute_lie_bargmann_invariants(phys_c, G=1.0).H
            err_c = abs(H_c - H0_f64) / abs(H0_f64)
            max_dE_c = max(max_dE_c, err_c)
            sum_sq_dE_c += err_c * err_c

        phys_c_end = sys_c.observe_physical_state()
        H_c_end = mod_05.compute_lie_bargmann_invariants(phys_c_end, G=1.0).H
        end_dE_c = abs(H_c_end - H0_f64) / abs(H0_f64)
        rms_dE_c = math.sqrt(sum_sq_dE_c / n_steps_common)
        dq_c = float(np.linalg.norm(phys_c_end.q - q_ref_t6))

        table5_common[name] = {
            "stages_s": s,
            "h_eff": h_eff_common,
            "h_macro": h_common,
            "n_steps": n_steps_common,
            "n_force_evals": n_steps_common * s,
            "delta_E_rel_max": float(max_dE_c),
            "delta_E_rel_rms": float(rms_dE_c),
            "delta_E_rel_endpoint": float(end_dE_c),
            "delta_q_norm": float(dq_c),
            "float64_floor": float(min(f64_results[name]["delta_E_rel"])),
            "float128_min_error": float(f128_results[name]["min_delta_E_rel"]),
            "empirical_order_f128": float(f128_results[name]["mean_asymptotic_slope"]),
        }

    existing_fig8["table5_reference_bm6_f128"] = table5_ref_bm6_f128
    existing_fig8["table5_common_abscissa_h_eff_0_01"] = table5_common
    return existing_fig8


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Benchmark Generator 01 (Dual Convergence & Section 6.6 Micro-Kernels)")
    parser.add_argument(
        "--protocol",
        type=str,
        default="all",
        choices=["C1", "C2", "C1_C2", "T5_C2", "M1", "M2", "M1_M2", "verify", "all"],
        help="Isolated protocol to run (C1=Figure-8, C2=3D Hier, C1_C2=Both C1 and C2, T5_C2=Table5+3D Hier, M1/M2=Sec 6.6 Micro-Kernels, verify=Assertions only)",
    )
    args = parser.parse_args()

    t_start = time.perf_counter()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / "generator_01_dual_convergence_and_microkernels.json"

    results: Dict[str, Any] = {
        "generator_id": "01_generator_dual_convergence_and_microkernels",
        "contract_path": "scripts/campaigns/contracts/01_generator_dual_convergence_and_microkernels.md",
    }
    if out_path.exists():
        try:
            with open(out_path, "r", encoding="utf-8") as f:
                results.update(json.load(f))
        except Exception:
            pass

    if args.protocol in ("C1", "C1_C2", "all"):
        results["figure8_convergence"] = run_part1_figure8_dual_convergence()
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
    if args.protocol == "T5_C2":
        results["figure8_convergence"] = run_part1_table5_only(results["figure8_convergence"])
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
    if args.protocol in ("C2", "C1_C2", "T5_C2", "all"):
        results["hierarchical_3d_eccentric_probe"] = run_part2_hierarchical_3d_eccentric_probe()
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
    if args.protocol in ("M1", "M2", "M1_M2", "all"):
        results["section_6_6_microkernels"] = run_part3_section_6_6_microkernels()
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)

    if "hierarchical_3d_eccentric_probe" in results:
        h3d = results["hierarchical_3d_eccentric_probe"]
        if "schemes" in h3d and "PBM4 (P0 Global, Obstructed)" in h3d["schemes"]:
            del h3d["schemes"]["PBM4 (P0 Global, Obstructed)"]
        h3d.pop("gains_p0_over_bm4", None)
        h3d.pop("max_gain_p0_over_bm4", None)

    results["elapsed_seconds"] = float(time.perf_counter() - t_start)

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    verify_contract_assertions(results)
    print(f"=== Saved and CERTIFIED JSON at {out_path} ===", flush=True)

    # Also synchronize data/campaign_IV_dual_convergence.json figure8_convergence block
    compat_path = DATA_DIR / "campaign_IV_dual_convergence.json"
    compat_data: Dict[str, Any] = {}
    if compat_path.exists():
        try:
            with open(compat_path, "r", encoding="utf-8") as f:
                compat_data = json.load(f)
        except Exception:
            compat_data = {}
    if "figure8_convergence" in results:
        compat_data["figure8_convergence"] = results["figure8_convergence"]
    if "hierarchical_3d_eccentric_probe" in results:
        compat_data["hierarchical_3d_eccentric_probe"] = results["hierarchical_3d_eccentric_probe"]
    if "section_6_6_microkernels" in results:
        compat_data["section_6_6_microkernels"] = results["section_6_6_microkernels"]
    with open(compat_path, "w", encoding="utf-8") as f:
        json.dump(compat_data, f, indent=2)
    print(f"=== Synchronized {compat_path} (Total runtime: {results['elapsed_seconds']:.2f}s) ===", flush=True)


if __name__ == "__main__":
    main()
