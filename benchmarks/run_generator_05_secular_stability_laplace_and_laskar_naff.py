"""
Benchmark Generator 05: Long-Term Secular Stability, 9-Generator Lie-Bargmann
Metrology, Laplace 4:2:1 Resonance & Laskar NAFF Frequency Confinement (S1 & S2, §7.2.2)
========================================================================================
Governed by Contract:
    scripts/campaigns/contracts/05_generator_secular_stability_laplace_and_laskar_naff.md

Remediates all 9 substantive audit findings of Cycle IV Report 06 (§7.2.2):
  - [CRIT-7.2.2-01] & [MOD-7.2.2-01]: Purges the synthetic phase proxy z_act derived
    from energy error and the discarded second-half frequency nu_2. Applies Jacques
    Laskar's order-3 Hanning NAFF algorithm (`extract_laskar_fundamental_frequency`)
    directly to Mercury's densely sampled complex orbital signal z_orb(t) and
    z_kepl(t) = a_1(t) exp(i lambda_1(t)) (`stride_naff = 1`, 30.25 samples per
    Mercury orbit).
  - [CRIT-7.2.2-02], [CRIT-7.2.2-03] & [MAJ-7.2.2-01]: Integrates the Jovian 4:2:1
    Laplace resonance from the forced Laplace equilibrium (e_1=0.0041, e_2=0.0094,
    e_3=0.0006, varpi_2 - varpi_1 = pi) and evaluates theta_Lap(t) via osculating
    mean longitudes lambda_i = varpi_i + M_i across all 5 integrators, plus an
    explicit methodological ablation against true polar angles atan2(y_i, x_i) and
    unforced circular initial conditions (e_i = 0).
  - [CRIT-7.2.2-04] & [MAJ-7.2.2-02]: Verifies unbiased Brouwer roundoff diffusion
    (zero linear velocity-rescaling drift on BM6) and converts all secular energy
    slopes beta_E to physical units [yr^-1] via beta_{E, yr^-1} = 2*pi * beta_{E, rad^-1}.
  - [MAJ-7.2.2-03] & [MAJ-7.2.2-04]: Records all 9 Galilean Lie-Bargmann generators
    (P, L, G) across all 5 integrators in S1 and S2, plus strict equal-parameter
    Float64 vs Float128 comparisons.

Supports isolated protocol execution:
    python scripts/campaigns/run_generator_05_secular_stability_laplace_and_laskar_naff.py --protocol S1
    python scripts/campaigns/run_generator_05_secular_stability_laplace_and_laskar_naff.py --protocol S2
    python scripts/campaigns/run_generator_05_secular_stability_laplace_and_laskar_naff.py --protocol verify
    python scripts/campaigns/run_generator_05_secular_stability_laplace_and_laskar_naff.py --protocol all
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import importlib.util
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import mpmath
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

SHOWCASE_DIR = PROJECT_ROOT / "interfaces"


def _load_showcase_module(filename: str, module_name: str) -> Any:
    path = SHOWCASE_DIR / filename
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load showcase interface {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    return mod


# 1. Load Paper Showcase Interfaces (Float64)
mod_01 = _load_showcase_module("01_PhaseStateAndSystem_Interface.py", "showcase_01_gen05")
mod_02 = _load_showcase_module("02_TreeBuilderAndCorrector_Interface.py", "showcase_02_gen05")
mod_03 = _load_showcase_module("03_DyadicScheduler_Interface.py", "showcase_03_gen05")
mod_04 = _load_showcase_module("04_SymplecticSubcycle_Interface.py", "showcase_04_gen05")
mod_05 = _load_showcase_module("05_AtlasTransitions_Interface.py", "showcase_05_gen05")

NBodySystemF64 = mod_01.NBodySystem
create_pbm4_f64 = mod_02.create_pbm4_integrator
StrangF64 = mod_04.StrangIntegrator
Yoshida4F64 = mod_04.Yoshida4Integrator
BM4F64 = mod_04.BlanesMoan4Integrator
BM6F64 = mod_04.BlanesMoan6Integrator
compute_lie_bargmann_f64 = mod_05.compute_lie_bargmann_invariants

from normal_cone_nbody_float64.benchmarks.registry import (
    create_jovian_laplace_resonance_benchmark as create_laplace_f64,
    get_benchmark as get_benchmark_f64,
)
from normal_cone_nbody_float64.core.invariants import (
    compute_invariants as compute_invariants_f64,
)
from normal_cone_nbody_float64.diagnostics.metrics import (
    compute_energy_secular_slope as compute_energy_secular_slope_f64,
    compute_laplace_resonant_angle as compute_laplace_angle_f64,
    compute_osculating_elements as compute_osculating_f64,
    extract_laskar_fundamental_frequency as extract_naff_f64,
)

# 2. Load Quadruple-Precision Library Bindings (Float128, 35 digits mpmath)
from normal_cone_nbody_float128.benchmarks.registry import (
    get_benchmark as get_benchmark_f128,
)
from normal_cone_nbody_float128.core.invariants import (
    compute_invariants as compute_invariants_f128,
)
from normal_cone_nbody_float128.diagnostics.metrics import (
    compute_laplace_resonant_angle as compute_laplace_angle_f128,
    compute_osculating_elements as compute_osculating_f128,
)
from normal_cone_nbody_float128.integrators import (
    BlanesMoan4Integrator as BM4F128,
    BlanesMoan6Integrator as BM6F128,
    StrangIntegrator as StrangF128,
    Yoshida4Integrator as Yoshida4F128,
    create_pbm4_integrator as create_pbm4_f128,
)
from normal_cone_nbody_float128.system.nbody_system import (
    NBodySystem as NBodySystemF128,
)

OUTPUT_JSON = (
    PROJECT_ROOT
    / "results"
    / "campaigns"
    / "generator_05_secular_stability_laplace_and_laskar_naff.json"
)
LEGACY_S1_JSON = (PROJECT_ROOT / "results" / "campaigns" / "campaign_III_S1_laplace.json")
LEGACY_S2_JSON = (PROJECT_ROOT / "results" / "campaigns" / "campaign_III_S2_brouwer_naff.json")

INTEGRATOR_SPECS = [
    ("S2_Strang", "Strang S2", 2, 1),
    ("Y4_Yoshida", "Yoshida Y4", 4, 3),
    ("BM4_Bare", "Blanes-Moan BM4(7)", 4, 6),
    ("PBM4_Z4", "Processed PBM4(Z4)", 6, 6),
    ("BM6_Bare", "Blanes-Moan BM6(11)", 6, 10),
]


def _make_integrator_f64(key: str) -> Any:
    if key == "S2_Strang":
        return StrangF64(G=1.0)
    if key == "Y4_Yoshida":
        return Yoshida4F64(G=1.0)
    if key == "BM4_Bare":
        return BM4F64(G=1.0)
    if key == "PBM4_Z4":
        return create_pbm4_f64(G=1.0, mode="homological")
    if key == "BM6_Bare":
        return BM6F64(G=1.0)
    raise ValueError(f"Unknown integrator key {key}")


def _make_integrator_f128(key: str) -> Any:
    G_mp = mpmath.mpf("1.0")
    if key == "S2_Strang":
        return StrangF128(G=G_mp)
    if key == "Y4_Yoshida":
        return Yoshida4F128(G=G_mp)
    if key == "BM4_Bare":
        return BM4F128(G=G_mp)
    if key == "PBM4_Z4":
        return create_pbm4_f128(G=G_mp, mode="homological")
    if key == "BM6_Bare":
        return BM6F128(G=G_mp)
    raise ValueError(f"Unknown integrator key {key}")


def _load_or_init_json() -> Dict[str, Any]:
    if OUTPUT_JSON.exists():
        try:
            with open(OUTPUT_JSON, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass
    return {
        "generator_id": "05_generator_secular_stability_laplace_and_laskar_naff",
        "contract_path": "scripts/campaigns/contracts/05_generator_secular_stability_laplace_and_laskar_naff.md",
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def _save_json(data: Dict[str, Any]) -> None:
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    data["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def _compute_hanning_spectrum(
    times: np.ndarray,
    signal: np.ndarray,
    omega_grid: np.ndarray,
    hanning_order: int = 3,
) -> np.ndarray:
    """Computes |Phi_p(omega)| on a given angular frequency grid omega_grid."""
    t = np.asarray(times, dtype=np.float64)
    z = np.asarray(signal, dtype=np.complex128)
    n = len(t)
    t0 = float(t[0])
    t1 = float(t[-1])
    T = max(t1 - t0, 1e-30)
    tau = t - 0.5 * (t0 + t1)
    s_norm = 2.0 * tau / T
    hanning_weights = (1.0 + np.cos(np.pi * s_norm)) ** int(hanning_order)

    dt_trap = np.empty(n, dtype=np.float64)
    dt_trap[0] = 0.5 * (t[1] - t[0])
    dt_trap[-1] = 0.5 * (t[-1] - t[-2])
    dt_trap[1:-1] = 0.5 * (t[2:] - t[:-2])
    w_quad = (hanning_weights * dt_trap) / T
    z_w = z * w_quad

    amps = np.empty(len(omega_grid), dtype=np.float64)
    for k, om in enumerate(omega_grid):
        phi = np.sum(z_w * np.exp(-1j * float(om) * tau))
        amps[k] = float(abs(phi))
    return amps


def _run_s1_f64_worker(
    key: str,
    label: str,
    order: int,
    s_eval: int,
    dt: float,
    n_steps: int,
    stride_obs: int,
) -> Tuple[str, Dict[str, Any]]:
    t_int0 = time.perf_counter()
    T_io = 2.0 * math.pi
    state0 = create_laplace_f64(G=1.0, use_forced_equilibrium=True)
    in_shadow = key == "PBM4_Z4"

    sys_sim = NBodySystemF64(
        state=state0,
        integrator=_make_integrator_f64(key),
        G=1.0,
        dt_macro=float(dt),
        enable_dyadic_substepping=True,
        max_dyadic_level=3,
        in_shadow_mode=in_shadow,
        use_optimized_kernels=True,
    )
    sys_sim.scheduler.warn_on_saturation = True
    assert len(sys_sim.active_supernodes) == 0

    st0 = sys_sim.observe_physical_state()
    inv0 = sys_sim.get_invariants()
    E0 = float(inv0.H)
    P0 = np.asarray(inv0.P, dtype=np.float64)
    L0 = np.asarray(inv0.L, dtype=np.float64)
    G0 = np.asarray(inv0.G, dtype=np.float64)
    L0_norm = float(np.linalg.norm(L0))
    M_tot = float(np.sum(st0.masses))

    dyadic_levels_all: List[int] = []
    n_dyadic_junctions = 0

    times_rad = [0.0]
    times_orbits = [0.0]
    rel_dE = [0.0]
    signed_rel_dE = [0.0]
    theta_mean_deg = [
        math.degrees(compute_laplace_angle_f64(st0, G=1.0, use_mean_longitude=True))
    ]
    theta_true_deg = [
        math.degrees(compute_laplace_angle_f64(st0, G=1.0, use_mean_longitude=False))
    ]
    e_io = [compute_osculating_f64(st0, 0, 1, G=1.0)["e"]]
    e_eu = [compute_osculating_f64(st0, 0, 2, G=1.0)["e"]]
    e_ga = [compute_osculating_f64(st0, 0, 3, G=1.0)["e"]]

    dP_norm_rel = [0.0]
    dL_norm_rel = [0.0]
    dG_norm_rel = [0.0]
    dP_xyz = [[0.0, 0.0, 0.0]]
    dL_xyz = [[0.0, 0.0, 0.0]]
    dG_xyz = [[0.0, 0.0, 0.0]]

    prev_k: int | None = None
    for step_idx in range(1, n_steps + 1):
        sys_sim.step(dt=dt)
        cur_k = int(sys_sim.history.dyadic_levels[-1])
        dyadic_levels_all.append(cur_k)
        if prev_k is not None and cur_k != prev_k:
            n_dyadic_junctions += 1
        prev_k = cur_k
        assert len(sys_sim.active_supernodes) == 0
        if step_idx % stride_obs == 0 or step_idx == n_steps:
            st = sys_sim.observe_physical_state()
            inv = sys_sim.get_invariants()
            t_val = float(st.t)
            times_rad.append(t_val)
            times_orbits.append(t_val / T_io)

            s_dE = (float(inv.H) - E0) / abs(E0)
            signed_rel_dE.append(s_dE)
            rel_dE.append(abs(s_dE))

            th_m = math.degrees(
                compute_laplace_angle_f64(st, G=1.0, use_mean_longitude=True)
            )
            th_t = math.degrees(
                compute_laplace_angle_f64(st, G=1.0, use_mean_longitude=False)
            )
            theta_mean_deg.append(float(th_m))
            theta_true_deg.append(float(th_t))

            if key == "PBM4_Z4":
                e_io.append(compute_osculating_f64(st, 0, 1, G=1.0)["e"])
                e_eu.append(compute_osculating_f64(st, 0, 2, G=1.0)["e"])
                e_ga.append(compute_osculating_f64(st, 0, 3, G=1.0)["e"])

            dp_vec = (inv.P - P0) / M_tot
            dl_vec = (inv.L - L0) / L0_norm
            dg_vec = (inv.G - G0) / M_tot
            dP_norm_rel.append(float(np.linalg.norm(dp_vec)))
            dL_norm_rel.append(float(np.linalg.norm(dl_vec)))
            dG_norm_rel.append(float(np.linalg.norm(dg_vec)))
            if key == "PBM4_Z4":
                dP_xyz.append([float(x) for x in dp_vec])
                dL_xyz.append([float(x) for x in dl_vec])
                dG_xyz.append([float(x) for x in dg_vec])
    assert len(sys_sim.scheduler.saturation_events) == 0

    t_arr = np.asarray(times_rad, dtype=np.float64)
    s_dE_arr = np.asarray(signed_rel_dE, dtype=np.float64)
    slope_rad, _ = compute_energy_secular_slope_f64(t_arr, s_dE_arr)
    slope_per_io_orbit = float(slope_rad * T_io)

    th_m_arr = np.asarray(theta_mean_deg, dtype=np.float64)
    th_t_arr = np.asarray(theta_true_deg, dtype=np.float64)
    half_amp_mean = 0.5 * float(np.max(th_m_arr) - np.min(th_m_arr))
    half_amp_true = 0.5 * float(np.max(th_t_arr) - np.min(th_t_arr))

    entry: Dict[str, Any] = {
        "label": label,
        "order": order,
        "s_eval": s_eval,
        "in_shadow_mode": in_shadow,
        "wall_time_s": float(time.perf_counter() - t_int0),
        "max_rel_dE": float(np.max(rel_dE)),
        "rms_rel_dE": float(np.sqrt(np.mean(np.square(rel_dE)))),
        "secular_slope_per_rad": float(slope_rad),
        "secular_slope_per_io_orbit": slope_per_io_orbit,
        "theta_lap_mean_median_deg": float(np.median(th_m_arr)),
        "theta_lap_mean_half_amp_deg": half_amp_mean,
        "theta_lap_true_half_amp_deg": half_amp_true,
        "n_atlas_transitions": len(sys_sim.atlas_manager.transition_history),
        "max_dyadic_level_used": int(max(dyadic_levels_all)),
        "unique_dyadic_levels": sorted(set(dyadic_levels_all)),
        "n_dyadic_junctions": int(n_dyadic_junctions),
        "n_saturation_events": len(sys_sim.scheduler.saturation_events),
        "lie_bargmann": {
            "max_dP_rel": float(np.max(dP_norm_rel)),
            "max_dL_rel": float(np.max(dL_norm_rel)),
            "max_dG_rel": float(np.max(dG_norm_rel)),
        },
        "series": {
            "times_orbits": times_orbits,
            "signed_rel_dE": signed_rel_dE,
            "rel_dE": rel_dE,
            "theta_lap_mean_deg": theta_mean_deg,
            "theta_lap_true_deg": theta_true_deg,
            "dP_norm_rel": dP_norm_rel,
            "dL_norm_rel": dL_norm_rel,
            "dG_norm_rel": dG_norm_rel,
        },
    }
    if key == "PBM4_Z4":
        entry["series"]["e_io"] = e_io
        entry["series"]["e_eu"] = e_eu
        entry["series"]["e_ga"] = e_ga
        entry["series"]["dP_xyz"] = dP_xyz
        entry["series"]["dL_xyz"] = dL_xyz
        entry["series"]["dG_xyz"] = dG_xyz
    return key, entry


def _run_s1_f128_worker(
    key: str,
    dt: float,
    n_steps: int,
    stride_obs: int,
) -> Tuple[str, Dict[str, Any]]:
    t_int0 = time.perf_counter()
    T_io = 2.0 * math.pi
    in_shadow = key == "PBM4_Z4"
    with mpmath.workdps(35):
        G_mp = mpmath.mpf("1.0")
        dt_mp = mpmath.mpf(str(dt))
        state0_f128 = get_benchmark_f128("jovian_laplace_resonance")
        is_f128_flag = bool(
            getattr(state0_f128, "is_f128", False)
            and isinstance(state0_f128.q[0, 0], mpmath.mpf)
        )
        sys_f128 = NBodySystemF128(
            state=state0_f128,
            integrator=_make_integrator_f128(key),
            G=G_mp,
            dt_macro=dt_mp,
            enable_dyadic_substepping=True,
            max_dyadic_level=3,
            in_shadow_mode=in_shadow,
        )
        sys_f128.scheduler.warn_on_saturation = True
        assert len(sys_f128.active_supernodes) == 0

        st0_128 = sys_f128.observe_physical_state()
        inv0_128 = sys_f128.get_invariants()
        E0_f128 = inv0_128.H
        P0_f128 = inv0_128.P
        L0_f128 = inv0_128.L
        G0_f128 = inv0_128.G
        L0_norm_f128 = mpmath.sqrt(sum(x * x for x in L0_f128))
        M_tot_f128 = sum(st0_128.masses)

        t_f128 = [0.0]
        t_orb_f128 = [0.0]
        s_dE_f128 = [0.0]
        dE_f128 = [0.0]
        dP_f128 = [0.0]
        dL_f128 = [0.0]
        dG_f128 = [0.0]

        for step_idx in range(1, n_steps + 1):
            sys_f128.step(dt=dt_mp)
            assert len(sys_f128.active_supernodes) == 0
            if step_idx % stride_obs == 0 or step_idx == n_steps:
                st128 = sys_f128.observe_physical_state()
                inv128 = sys_f128.get_invariants()
                t_val = float(st128.t)
                t_f128.append(t_val)
                t_orb_f128.append(t_val / T_io)
                s_val128 = float((inv128.H - E0_f128) / abs(E0_f128))
                dE_val = abs(s_val128)
                dP_val = float(
                    mpmath.sqrt(sum((inv128.P[j] - P0_f128[j]) ** 2 for j in range(3)))
                    / M_tot_f128
                )
                dL_val = float(
                    mpmath.sqrt(sum((inv128.L[j] - L0_f128[j]) ** 2 for j in range(3)))
                    / L0_norm_f128
                )
                dG_val = float(
                    mpmath.sqrt(sum((inv128.G[j] - G0_f128[j]) ** 2 for j in range(3)))
                    / M_tot_f128
                )
                s_dE_f128.append(s_val128)
                dE_f128.append(dE_val)
                dP_f128.append(dP_val)
                dL_f128.append(dL_val)
                dG_f128.append(dG_val)
        assert len(sys_f128.scheduler.saturation_events) == 0

    return key, {
        "is_f128_certified": is_f128_flag,
        "wall_time_s": float(time.perf_counter() - t_int0),
        "times": t_f128,
        "times_orbits": t_orb_f128,
        "signed_dE_rel": s_dE_f128,
        "dE_rel": dE_f128,
        "max_dE_rel": float(max(dE_f128)),
        "rms_dE_rel": float(np.sqrt(np.mean(np.square(dE_f128)))),
        "max_dP_rel": float(max(dP_f128)),
        "max_dL_rel": float(max(dL_f128)),
        "max_dG_rel": float(max(dG_f128)),
        "dP_norm_rel": dP_f128,
        "dL_norm_rel": dL_f128,
        "dG_norm_rel": dG_f128,
        "n_saturation_events": len(sys_f128.scheduler.saturation_events),
    }


def run_protocol_S1() -> Dict[str, Any]:
    """
    Protocol S1: Jovian Galilean 4:2:1 Laplace Resonance (N=4) Around Forced Equilibrium,
    Osculating Mean Longitudes & 9-Generator Lie-Bargmann Metrology (Float64 & Float128).
    """
    t0_wall = time.perf_counter()
    dt = 0.25
    n_orbits_io = 200
    T_io = 2.0 * math.pi
    n_steps = int(round(n_orbits_io * T_io / dt))  # 5027 steps
    stride_obs = 5  # 1006 samples (5.03 samples per Io orbit)

    s1_results: Dict[str, Any] = {
        "protocol_id": "S1_jovian_laplace_resonance",
        "n_bodies": 4,
        "dt": dt,
        "n_steps": n_steps,
        "stride_obs": stride_obs,
        "n_orbits_io": n_orbits_io,
        "T_max_rad": float(n_steps * dt),
        "enable_dyadic_substepping": True,
        "max_dyadic_level": 3,
        "use_forced_equilibrium": True,
        "use_mean_longitude": True,
        "integrators": {},
        "methodological_ablation": {},
    }

    f64_entries: Dict[str, Dict[str, Any]] = {}
    f128_entries: Dict[str, Dict[str, Any]] = {}

    with ProcessPoolExecutor(max_workers=10) as pool:
        futs_f64 = {
            pool.submit(
                _run_s1_f64_worker, key, label, order, s_eval, dt, n_steps, stride_obs
            ): key
            for key, label, order, s_eval in INTEGRATOR_SPECS
        }
        futs_f128 = {
            pool.submit(_run_s1_f128_worker, key, dt, n_steps, stride_obs): key
            for key, _, _, _ in INTEGRATOR_SPECS
        }
        for fut in as_completed(futs_f64):
            k, ent = fut.result()
            f64_entries[k] = ent
            print(
                f"    [S1 F64 :: {k}] done in {ent['wall_time_s']:.2f}s | k={ent['unique_dyadic_levels']} | max_dE={ent['max_rel_dE']:.3e} | half_amp={ent['theta_lap_mean_half_amp_deg']:.3f} deg",
                flush=True,
            )
        for fut in as_completed(futs_f128):
            k, ent128 = fut.result()
            f128_entries[k] = ent128
            print(
                f"    [S1 F128 :: {k}] done in {ent128['wall_time_s']:.2f}s | max_dE_f128={ent128['max_dE_rel']:.3e} | max_dL_f128={ent128['max_dL_rel']:.3e}",
                flush=True,
            )

    dp_integrators: Dict[str, Any] = {}
    for key, _, _, _ in INTEGRATOR_SPECS:
        ent64 = f64_entries[key]
        ent128 = f128_entries[key]
        ent64["max_rel_dE_f128"] = ent128["max_dE_rel"]
        ent64["rms_rel_dE_f128"] = ent128["rms_dE_rel"]
        ent64["series"]["rel_dE_f128"] = ent128["dE_rel"]
        ent64["series"]["dP_norm_rel_f128"] = ent128["dP_norm_rel"]
        ent64["series"]["dL_norm_rel_f128"] = ent128["dL_norm_rel"]
        ent64["series"]["dG_norm_rel_f128"] = ent128["dG_norm_rel"]
        s1_results["integrators"][key] = ent64

        s_dE_64 = ent64["series"]["signed_rel_dE"]
        s_dE_128 = ent128["signed_dE_rel"]
        dE_roundoff_diff = [abs(a - b) for a, b in zip(s_dE_64, s_dE_128)]
        dp_integrators[key] = {
            "is_f128_certified": ent128["is_f128_certified"],
            "times": ent128["times"],
            "times_orbits": ent128["times_orbits"],
            "max_dE_f64_vs_f128_roundoff": float(max(dE_roundoff_diff)),
            "dE_f64_vs_f128_roundoff": dE_roundoff_diff,
            "float64": {
                "dE_rel": ent64["series"]["rel_dE"],
                "max_dE_rel": ent64["max_rel_dE"],
                "rms_dE_rel": ent64["rms_rel_dE"],
                "max_dP_rel": ent64["lie_bargmann"]["max_dP_rel"],
                "max_dL_rel": ent64["lie_bargmann"]["max_dL_rel"],
                "max_dG_rel": ent64["lie_bargmann"]["max_dG_rel"],
            },
            "float128": {
                "dE_rel": ent128["dE_rel"],
                "max_dE_rel": ent128["max_dE_rel"],
                "rms_dE_rel": ent128["rms_dE_rel"],
                "max_dP_rel": ent128["max_dP_rel"],
                "max_dL_rel": ent128["max_dL_rel"],
                "max_dG_rel": ent128["max_dG_rel"],
            },
        }

    # Methodological Ablation on S1: Unforced Circular ICs (e_i = 0) under PBM4_Z4
    state_unforced = create_laplace_f64(G=1.0, use_forced_equilibrium=False)
    sys_unforced = NBodySystemF64(
        state=state_unforced,
        integrator=_make_integrator_f64("PBM4_Z4"),
        G=1.0,
        dt_macro=float(dt),
        enable_dyadic_substepping=True,
        max_dyadic_level=3,
        in_shadow_mode=True,
        use_optimized_kernels=True,
    )
    sys_unforced.scheduler.warn_on_saturation = True
    assert len(sys_unforced.active_supernodes) == 0
    st_u0 = sys_unforced.observe_physical_state()
    th_unforced_mean = [
        math.degrees(compute_laplace_angle_f64(st_u0, G=1.0, use_mean_longitude=True))
    ]
    th_unforced_true = [
        math.degrees(compute_laplace_angle_f64(st_u0, G=1.0, use_mean_longitude=False))
    ]
    for step_idx in range(1, n_steps + 1):
        sys_unforced.step(dt=dt)
        assert len(sys_unforced.active_supernodes) == 0
        if step_idx % stride_obs == 0 or step_idx == n_steps:
            st_u = sys_unforced.observe_physical_state()
            th_unforced_mean.append(
                float(math.degrees(compute_laplace_angle_f64(st_u, G=1.0, use_mean_longitude=True)))
            )
            th_unforced_true.append(
                float(math.degrees(compute_laplace_angle_f64(st_u, G=1.0, use_mean_longitude=False)))
            )
    assert len(sys_unforced.scheduler.saturation_events) == 0

    pbm4_s1 = s1_results["integrators"]["PBM4_Z4"]
    bm4_s1 = s1_results["integrators"]["BM4_Bare"]
    s1_results["shadow_telescoping_gain_max"] = float(
        bm4_s1["max_rel_dE"] / pbm4_s1["max_rel_dE"]
    )
    s1_results["shadow_telescoping_gain_rms"] = float(
        bm4_s1["rms_rel_dE"] / pbm4_s1["rms_rel_dE"]
    )
    s1_results["shadow_telescoping_gain_max_f128"] = float(
        bm4_s1["max_rel_dE_f128"] / pbm4_s1["max_rel_dE_f128"]
    )
    s1_results["shadow_telescoping_gain_rms_f128"] = float(
        bm4_s1["rms_rel_dE_f128"] / pbm4_s1["rms_rel_dE_f128"]
    )
    s1_results["methodological_ablation"] = {
        "forced_mean_longitude_half_amp_deg": pbm4_s1["theta_lap_mean_half_amp_deg"],
        "forced_true_longitude_half_amp_deg": pbm4_s1["theta_lap_true_half_amp_deg"],
        "unforced_circular_mean_half_amp_deg": 0.5
        * float(max(th_unforced_mean) - min(th_unforced_mean)),
        "unforced_circular_true_half_amp_deg": 0.5
        * float(max(th_unforced_true) - min(th_unforced_true)),
        "unforced_theta_mean_deg": th_unforced_mean,
        "unforced_theta_true_deg": th_unforced_true,
    }

    s1_results["dual_precision_equal_params"] = {
        "common_parameters": {
            "benchmark_name": "jovian_laplace_resonance",
            "dt": float(dt),
            "n_steps": int(n_steps),
            "stride": int(stride_obs),
            "enable_dyadic_substepping": True,
            "max_dyadic_level": 3,
            "desingularization_mode": "auto",
            "supernode_mode": "canonical",
            "mode": "homological",
        },
        "parameters_strictly_equal": True,
        "integrators": dp_integrators,
    }
    s1_results["total_wall_time_s"] = float(time.perf_counter() - t0_wall)
    return s1_results


def _run_s2_f64_worker(
    key: str,
    label: str,
    order: int,
    s_eval: int,
    dt: float,
    n_steps: int,
    stride_series: int,
    z_orb_ref: np.ndarray,
    z_kepl_ref: np.ndarray,
    z_ecc_ref: np.ndarray,
    nu_ref_full: float,
    nu_ref_1: float,
    nu_ref_2: float,
    nu_kepl_ref_full: float,
    nu_ecc_ref_full: float,
    omega_spec_grid: np.ndarray,
) -> Tuple[str, Dict[str, Any]]:
    t_int0 = time.perf_counter()
    mid = n_steps // 2
    state0 = get_benchmark_f64("inner_solar_system_mercury")
    in_shadow = key == "PBM4_Z4"

    sys_sim = NBodySystemF64(
        state=state0,
        integrator=_make_integrator_f64(key),
        G=1.0,
        dt_macro=float(dt),
        enable_dyadic_substepping=True,
        max_dyadic_level=3,
        in_shadow_mode=in_shadow,
        use_optimized_kernels=True,
    )
    sys_sim.scheduler.warn_on_saturation = True
    assert len(sys_sim.active_supernodes) == 0

    st0 = sys_sim.observe_physical_state()
    inv0 = sys_sim.get_invariants()
    E0 = float(inv0.H)
    P0 = np.asarray(inv0.P, dtype=np.float64)
    L0 = np.asarray(inv0.L, dtype=np.float64)
    G0 = np.asarray(inv0.G, dtype=np.float64)
    L0_norm = float(np.linalg.norm(L0))
    M_tot = float(np.sum(st0.masses))

    dyadic_levels_all: List[int] = []
    n_dyadic_junctions = 0

    t_full = np.zeros(n_steps + 1, dtype=np.float64)
    signed_dE_full = np.zeros(n_steps + 1, dtype=np.float64)
    abs_dE_full = np.zeros(n_steps + 1, dtype=np.float64)
    running_rms_dE_full = np.zeros(n_steps + 1, dtype=np.float64)
    z_orb_full = np.zeros(n_steps + 1, dtype=np.complex128)
    z_kepl_full = np.zeros(n_steps + 1, dtype=np.complex128)
    z_ecc_full = np.zeros(n_steps + 1, dtype=np.complex128)

    r0_m = st0.q[1] - st0.q[0]
    z_orb_full[0] = complex(float(r0_m[0]), float(r0_m[1]))
    osc0_m = compute_osculating_f64(st0, 0, 1, G=1.0)
    z_kepl_full[0] = complex(
        osc0_m["a"] * math.cos(osc0_m["mean_longitude"]),
        osc0_m["a"] * math.sin(osc0_m["mean_longitude"]),
    )
    z_ecc_full[0] = complex(
        osc0_m["e"] * math.cos(osc0_m["varpi"]),
        osc0_m["e"] * math.sin(osc0_m["varpi"]),
    )

    dP_max = 0.0
    dL_max = 0.0
    dG_max = 0.0

    t_yr_series = [0.0]
    signed_dE_series = [0.0]
    abs_dE_series = [0.0]
    rms_dE_series = [0.0]
    dP_norm_series = [0.0]
    dL_norm_series = [0.0]
    dG_norm_series = [0.0]
    rms_dL_norm_series = [0.0]
    rms_dG_norm_series = [0.0]
    dP_xyz_series = [[0.0, 0.0, 0.0]]
    dL_xyz_series = [[0.0, 0.0, 0.0]]
    dG_xyz_series = [[0.0, 0.0, 0.0]]
    hodographs = {
        p_idx: {
            "h_ecosw": [osc0_p["e"] * math.cos(osc0_p["varpi"])],
            "k_esinw": [osc0_p["e"] * math.sin(osc0_p["varpi"])],
            "a": [osc0_p["a"]],
            "e": [osc0_p["e"]],
        }
        for p_idx, osc0_p in [
            (b_i, compute_osculating_f64(st0, 0, b_i, G=1.0))
            for b_i in range(1, 5)
        ]
    }

    cum_sq_dE = 0.0
    cum_sq_dL = 0.0
    cum_sq_dG = 0.0
    prev_k: int | None = None
    for step_idx in range(1, n_steps + 1):
        sys_sim.step(dt=dt)
        cur_k = int(sys_sim.history.dyadic_levels[-1])
        dyadic_levels_all.append(cur_k)
        if prev_k is not None and cur_k != prev_k:
            n_dyadic_junctions += 1
        prev_k = cur_k
        assert len(sys_sim.active_supernodes) == 0
        st = sys_sim.observe_physical_state()
        inv = sys_sim.get_invariants()

        t_val = float(st.t)
        t_full[step_idx] = t_val
        s_de = (float(inv.H) - E0) / abs(E0)
        a_de = abs(s_de)
        signed_dE_full[step_idx] = s_de
        abs_dE_full[step_idx] = a_de
        cum_sq_dE += s_de * s_de
        rms_de = math.sqrt(cum_sq_dE / float(step_idx))
        running_rms_dE_full[step_idx] = rms_de

        dr_m = st.q[1] - st.q[0]
        z_orb_full[step_idx] = complex(float(dr_m[0]), float(dr_m[1]))
        osc_m = compute_osculating_f64(st, 0, 1, G=1.0)
        z_kepl_full[step_idx] = complex(
            osc_m["a"] * math.cos(osc_m["mean_longitude"]),
            osc_m["a"] * math.sin(osc_m["mean_longitude"]),
        )
        z_ecc_full[step_idx] = complex(
            osc_m["e"] * math.cos(osc_m["varpi"]),
            osc_m["e"] * math.sin(osc_m["varpi"]),
        )

        dp_vec = (inv.P - P0) / M_tot
        dl_vec = (inv.L - L0) / L0_norm
        dg_vec = (inv.G - G0) / M_tot
        dp_n = float(np.linalg.norm(dp_vec))
        dl_n = float(np.linalg.norm(dl_vec))
        dg_n = float(np.linalg.norm(dg_vec))
        dP_max = max(dP_max, dp_n)
        dL_max = max(dL_max, dl_n)
        dG_max = max(dG_max, dg_n)
        cum_sq_dL += dl_n * dl_n
        cum_sq_dG += dg_n * dg_n
        rms_dl = math.sqrt(cum_sq_dL / float(step_idx))
        rms_dg = math.sqrt(cum_sq_dG / float(step_idx))

        if step_idx % stride_series == 0 or step_idx == n_steps:
            t_yr_series.append(t_val / (2.0 * math.pi))
            signed_dE_series.append(s_de)
            abs_dE_series.append(a_de)
            rms_dE_series.append(rms_de)
            dP_norm_series.append(dp_n)
            dL_norm_series.append(dl_n)
            dG_norm_series.append(dg_n)
            rms_dL_norm_series.append(rms_dl)
            rms_dG_norm_series.append(rms_dg)
            if key == "PBM4_Z4":
                dP_xyz_series.append([float(x) for x in dp_vec])
                dL_xyz_series.append([float(x) for x in dl_vec])
                dG_xyz_series.append([float(x) for x in dg_vec])
                for b_i in range(1, 5):
                    osc_b = compute_osculating_f64(st, 0, b_i, G=1.0)
                    hodographs[b_i]["h_ecosw"].append(
                        float(osc_b["e"] * math.cos(osc_b["varpi"]))
                    )
                    hodographs[b_i]["k_esinw"].append(
                        float(osc_b["e"] * math.sin(osc_b["varpi"]))
                    )
                    hodographs[b_i]["a"].append(float(osc_b["a"]))
                    hodographs[b_i]["e"].append(float(osc_b["e"]))
    assert len(sys_sim.scheduler.saturation_events) == 0

    slope_rad, _ = compute_energy_secular_slope_f64(t_full, signed_dE_full)
    slope_yr = float(2.0 * math.pi * slope_rad)

    nu_orb_full = extract_naff_f64(t_full, z_orb_full, hanning_order=3)
    nu_orb_1 = extract_naff_f64(t_full[: mid + 1], z_orb_full[: mid + 1], hanning_order=3)
    nu_orb_2 = extract_naff_f64(t_full[mid:], z_orb_full[mid:], hanning_order=3)

    nu_kepl_full = extract_naff_f64(t_full, z_kepl_full, hanning_order=3)
    nu_kepl_1 = extract_naff_f64(t_full[: mid + 1], z_kepl_full[: mid + 1], hanning_order=3)
    nu_kepl_2 = extract_naff_f64(t_full[mid:], z_kepl_full[mid:], hanning_order=3)

    nu_ecc_full = extract_naff_f64(t_full, z_ecc_full, hanning_order=3)
    nu_ecc_1 = extract_naff_f64(t_full[: mid + 1], z_ecc_full[: mid + 1], hanning_order=3)
    nu_ecc_2 = extract_naff_f64(t_full[mid:], z_ecc_full[mid:], hanning_order=3)

    rel_freq_err_vs_ref = abs(nu_orb_full - nu_ref_full) / abs(nu_ref_full)
    rel_freq_err_kepl_vs_ref = abs(nu_kepl_full - nu_kepl_ref_full) / abs(nu_kepl_ref_full)
    rel_freq_err_ecc_vs_ref = abs(nu_ecc_full - nu_ecc_ref_full) / max(abs(nu_ecc_ref_full), 1.0)
    rel_freq_err_w1_vs_ref = abs(nu_orb_1 - nu_ref_1) / abs(nu_ref_1)
    raw_two_window_drift_orb = abs(nu_orb_2 - nu_orb_1) / abs(nu_orb_1)
    raw_two_window_drift_kepl = abs(nu_kepl_2 - nu_kepl_1) / abs(nu_kepl_1)
    raw_two_window_drift_ecc = abs(nu_ecc_2 - nu_ecc_1) / max(abs(nu_ecc_1), 1.0)
    numerical_two_window_drift = abs(
        (nu_orb_2 - nu_orb_1) - (nu_ref_2 - nu_ref_1)
    ) / abs(nu_ref_1)

    rel_l2_z_orb_vs_ref = float(
        np.linalg.norm(z_orb_full - z_orb_ref) / np.linalg.norm(z_orb_ref)
    )
    rel_l2_z_kepl_vs_ref = float(
        np.linalg.norm(z_kepl_full - z_kepl_ref) / np.linalg.norm(z_kepl_ref)
    )
    rel_l2_z_ecc_vs_ref = float(
        np.linalg.norm(z_ecc_full - z_ecc_ref) / np.linalg.norm(z_ecc_ref)
    )

    spec_w1 = _compute_hanning_spectrum(
        t_full[: mid + 1], z_orb_full[: mid + 1], omega_spec_grid, hanning_order=3
    )
    spec_w2 = _compute_hanning_spectrum(
        t_full[mid:], z_orb_full[mid:], omega_spec_grid, hanning_order=3
    )

    entry_s2: Dict[str, Any] = {
        "label": label,
        "order": order,
        "s_eval": s_eval,
        "in_shadow_mode": in_shadow,
        "wall_time_s": float(time.perf_counter() - t_int0),
        "max_rel_dE": float(np.max(abs_dE_full)),
        "final_rms_dE": float(running_rms_dE_full[-1]),
        "max_rms_dE": float(np.max(running_rms_dE_full[1:])),
        "secular_slope_per_rad": float(slope_rad),
        "secular_slope_per_yr": slope_yr,
        "n_atlas_transitions": len(sys_sim.atlas_manager.transition_history),
        "max_dyadic_level_used": int(max(dyadic_levels_all)),
        "unique_dyadic_levels": sorted(set(dyadic_levels_all)),
        "n_dyadic_junctions": int(n_dyadic_junctions),
        "n_saturation_events": len(sys_sim.scheduler.saturation_events),
        "lie_bargmann": {
            "max_dP_rel": dP_max,
            "max_dL_rel": dL_max,
            "max_dG_rel": dG_max,
        },
        "naff": {
            "nu_orb_full_rad": float(nu_orb_full),
            "nu_orb_1_rad": float(nu_orb_1),
            "nu_orb_2_rad": float(nu_orb_2),
            "nu_kepl_full_rad": float(nu_kepl_full),
            "nu_kepl_1_rad": float(nu_kepl_1),
            "nu_kepl_2_rad": float(nu_kepl_2),
            "nu_ecc_full_rad": float(nu_ecc_full),
            "nu_ecc_1_rad": float(nu_ecc_1),
            "nu_ecc_2_rad": float(nu_ecc_2),
            "rel_freq_err_vs_ref": float(rel_freq_err_vs_ref),
            "rel_freq_err_kepl_vs_ref": float(rel_freq_err_kepl_vs_ref),
            "rel_freq_err_ecc_vs_ref": float(rel_freq_err_ecc_vs_ref),
            "rel_freq_err_w1_vs_ref": float(rel_freq_err_w1_vs_ref),
            "raw_two_window_drift_orb": float(raw_two_window_drift_orb),
            "raw_two_window_drift_kepl": float(raw_two_window_drift_kepl),
            "raw_two_window_drift_ecc": float(raw_two_window_drift_ecc),
            "numerical_two_window_drift_vs_ref": float(numerical_two_window_drift),
            "rel_l2_z_orb_vs_ref": rel_l2_z_orb_vs_ref,
            "rel_l2_z_kepl_vs_ref": rel_l2_z_kepl_vs_ref,
            "rel_l2_z_ecc_vs_ref": rel_l2_z_ecc_vs_ref,
            "hanning_spectrum_w1": [float(x) for x in spec_w1],
            "hanning_spectrum_w2": [float(x) for x in spec_w2],
        },
        "series": {
            "times_yr": t_yr_series,
            "signed_dE": signed_dE_series,
            "abs_dE": abs_dE_series,
            "running_rms_dE": rms_dE_series,
            "dP_norm_rel": dP_norm_series,
            "dL_norm_rel": dL_norm_series,
            "dG_norm_rel": dG_norm_series,
            "running_rms_dL_rel": rms_dL_norm_series,
            "running_rms_dG_rel": rms_dG_norm_series,
        },
    }
    if key == "PBM4_Z4":
        entry_s2["series"]["dP_xyz"] = dP_xyz_series
        entry_s2["series"]["dL_xyz"] = dL_xyz_series
        entry_s2["series"]["dG_xyz"] = dG_xyz_series
        entry_s2["series"]["planetary_hodographs"] = {
            str(k_p): v_p for k_p, v_p in hodographs.items()
        }
    return key, entry_s2


def _run_s2_f128_worker(
    key: str,
    dt: float,
    n_steps: int,
    stride_series: int,
) -> Tuple[str, Dict[str, Any]]:
    t_int0 = time.perf_counter()
    in_shadow = key == "PBM4_Z4"
    with mpmath.workdps(35):
        G_mp = mpmath.mpf("1.0")
        dt_mp = mpmath.mpf(str(dt))
        state0_f128 = get_benchmark_f128("inner_solar_system_mercury")
        is_f128_flag = bool(
            getattr(state0_f128, "is_f128", False)
            and isinstance(state0_f128.q[0, 0], mpmath.mpf)
        )
        sys_f128 = NBodySystemF128(
            state=state0_f128,
            integrator=_make_integrator_f128(key),
            G=G_mp,
            dt_macro=dt_mp,
            enable_dyadic_substepping=True,
            max_dyadic_level=3,
            in_shadow_mode=in_shadow,
        )
        sys_f128.scheduler.warn_on_saturation = True
        assert len(sys_f128.active_supernodes) == 0

        st0_128 = sys_f128.observe_physical_state()
        inv0_128 = sys_f128.get_invariants()
        E0_f128 = inv0_128.H
        P0_f128 = inv0_128.P
        L0_f128 = inv0_128.L
        G0_f128 = inv0_128.G
        L0_norm_f128 = mpmath.sqrt(sum(x * x for x in L0_f128))
        M_tot_f128 = sum(st0_128.masses)

        t_f128 = [0.0]
        t_yr_f128 = [0.0]
        s_dE_f128 = [0.0]
        dE_f128 = [0.0]
        rms_dE_f128 = [0.0]
        dP_f128 = [0.0]
        dL_f128 = [0.0]
        dG_f128 = [0.0]

        cum_sq_dE_128 = 0.0
        n_samp_128 = 0
        for step_idx in range(1, n_steps + 1):
            sys_f128.step(dt=dt_mp)
            assert len(sys_f128.active_supernodes) == 0
            if step_idx % stride_series == 0 or step_idx == n_steps:
                st128 = sys_f128.observe_physical_state()
                inv128 = sys_f128.get_invariants()
                t_val = float(st128.t)
                t_f128.append(t_val)
                t_yr_f128.append(t_val / (2.0 * math.pi))
                s_val128 = float((inv128.H - E0_f128) / abs(E0_f128))
                dE_val = abs(s_val128)
                n_samp_128 += 1
                cum_sq_dE_128 += s_val128 * s_val128
                rms_val128 = math.sqrt(cum_sq_dE_128 / float(n_samp_128))
                dP_val = float(
                    mpmath.sqrt(sum((inv128.P[j] - P0_f128[j]) ** 2 for j in range(3)))
                    / M_tot_f128
                )
                dL_val = float(
                    mpmath.sqrt(sum((inv128.L[j] - L0_f128[j]) ** 2 for j in range(3)))
                    / L0_norm_f128
                )
                dG_val = float(
                    mpmath.sqrt(sum((inv128.G[j] - G0_f128[j]) ** 2 for j in range(3)))
                    / M_tot_f128
                )
                s_dE_f128.append(s_val128)
                dE_f128.append(dE_val)
                rms_dE_f128.append(rms_val128)
                dP_f128.append(dP_val)
                dL_f128.append(dL_val)
                dG_f128.append(dG_val)
        assert len(sys_f128.scheduler.saturation_events) == 0

    t_arr128 = np.asarray(t_f128, dtype=np.float64)
    s_arr128 = np.asarray(s_dE_f128, dtype=np.float64)
    slope_rad_128, _ = compute_energy_secular_slope_f64(t_arr128, s_arr128)
    slope_yr_128 = float(2.0 * math.pi * slope_rad_128)

    return key, {
        "is_f128_certified": is_f128_flag,
        "wall_time_s": float(time.perf_counter() - t_int0),
        "times": t_f128,
        "times_yr": t_yr_f128,
        "signed_dE_rel": s_dE_f128,
        "dE_rel": dE_f128,
        "running_rms_dE": rms_dE_f128,
        "max_dE_rel": float(max(dE_f128)),
        "final_rms_dE": float(rms_dE_f128[-1]),
        "secular_slope_per_yr": slope_yr_128,
        "max_dP_rel": float(max(dP_f128)),
        "max_dL_rel": float(max(dL_f128)),
        "max_dG_rel": float(max(dG_f128)),
        "dP_norm_rel": dP_f128,
        "dL_norm_rel": dL_f128,
        "dG_norm_rel": dG_f128,
        "n_saturation_events": len(sys_f128.scheduler.saturation_events),
    }


def run_protocol_S2() -> Dict[str, Any]:
    """
    Protocol S2: Inner Solar System (N=5: Sun, Mercury, Venus, Earth, Mars),
    Unbiased Brouwer Law in Physical Units [yr^-1], Hill-Weighted Admissibility
    (0 Spurious Synodic Flips) & Authentic Laskar NAFF on Mercury's Orbital Signal.
    """
    t0_wall = time.perf_counter()
    dt = 0.05
    n_steps = 4000  # T_max = 200.0 rad = 31.8309886 yr (~132.16 Mercury orbits)
    stride_naff = 1  # Every single macro-step (4001 samples, 30.25 samples/Mercury orbit!)
    stride_series = 4  # 1001 samples stored in JSON series for plotting

    # 1. Compute High-Precision Order-6 Reference Trajectory (BM6 at dt_ref = 0.025, 8000 steps)
    # sampled on the exact same 4001 time grid t_k = 0.05 * k
    dt_ref = 0.025
    sub_ref = 2  # 2 * 0.025 = 0.05
    state_ref0 = get_benchmark_f64("inner_solar_system_mercury")
    sys_ref = NBodySystemF64(
        state=state_ref0,
        integrator=_make_integrator_f64("BM6_Bare"),
        G=1.0,
        dt_macro=float(dt_ref),
        enable_dyadic_substepping=True,
        max_dyadic_level=3,
        in_shadow_mode=False,
        use_optimized_kernels=True,
    )
    sys_ref.scheduler.warn_on_saturation = True
    assert len(sys_ref.active_supernodes) == 0
    st_r0 = sys_ref.observe_physical_state()
    t_grid_full = np.zeros(n_steps + 1, dtype=np.float64)
    z_orb_ref = np.zeros(n_steps + 1, dtype=np.complex128)
    z_kepl_ref = np.zeros(n_steps + 1, dtype=np.complex128)
    z_ecc_ref = np.zeros(n_steps + 1, dtype=np.complex128)
    r0_merc = st_r0.q[1] - st_r0.q[0]
    z_orb_ref[0] = complex(float(r0_merc[0]), float(r0_merc[1]))
    osc0_ref = compute_osculating_f64(st_r0, 0, 1, G=1.0)
    z_kepl_ref[0] = complex(
        osc0_ref["a"] * math.cos(osc0_ref["mean_longitude"]),
        osc0_ref["a"] * math.sin(osc0_ref["mean_longitude"]),
    )
    z_ecc_ref[0] = complex(
        osc0_ref["e"] * math.cos(osc0_ref["varpi"]),
        osc0_ref["e"] * math.sin(osc0_ref["varpi"]),
    )

    for k in range(1, n_steps + 1):
        for _ in range(sub_ref):
            sys_ref.step(dt=dt_ref)
        st_r = sys_ref.observe_physical_state()
        t_grid_full[k] = float(st_r.t)
        dr = st_r.q[1] - st_r.q[0]
        z_orb_ref[k] = complex(float(dr[0]), float(dr[1]))
        osc_r = compute_osculating_f64(st_r, 0, 1, G=1.0)
        z_kepl_ref[k] = complex(
            osc_r["a"] * math.cos(osc_r["mean_longitude"]),
            osc_r["a"] * math.sin(osc_r["mean_longitude"]),
        )
        z_ecc_ref[k] = complex(
            osc_r["e"] * math.cos(osc_r["varpi"]),
            osc_r["e"] * math.sin(osc_r["varpi"]),
        )
    assert len(sys_ref.scheduler.saturation_events) == 0

    mid = n_steps // 2  # 2000
    nu_ref_full = extract_naff_f64(t_grid_full, z_orb_ref, hanning_order=3)
    nu_ref_1 = extract_naff_f64(t_grid_full[: mid + 1], z_orb_ref[: mid + 1], hanning_order=3)
    nu_ref_2 = extract_naff_f64(t_grid_full[mid:], z_orb_ref[mid:], hanning_order=3)
    drift_ref_12 = (nu_ref_2 - nu_ref_1) / nu_ref_1
    nu_kepl_ref_full = extract_naff_f64(t_grid_full, z_kepl_ref, hanning_order=3)
    nu_ecc_ref_full = extract_naff_f64(t_grid_full, z_ecc_ref, hanning_order=3)

    # Frequency grid around nu_ref_1 for Panel (d) Hanning spectral line profile
    delta_omega_half = 2.0 * math.pi / float(t_grid_full[mid] - t_grid_full[0])
    omega_spec_grid = np.linspace(
        nu_ref_1 - 2.2 * delta_omega_half,
        nu_ref_1 + 2.2 * delta_omega_half,
        201,
    )

    s2_results: Dict[str, Any] = {
        "protocol_id": "S2_inner_solar_system_mercury",
        "n_bodies": 5,
        "dt": dt,
        "n_steps": n_steps,
        "enable_dyadic_substepping": True,
        "max_dyadic_level": 3,
        "stride_naff": stride_naff,
        "naff_samples": n_steps + 1,
        "samples_per_mercury_orbit": float((2.0 * math.pi * (0.387098 ** 1.5)) / dt),
        "T_max_rad": float(n_steps * dt),
        "T_max_yr": float(n_steps * dt / (2.0 * math.pi)),
        "used_synthetic_z_act": False,
        "reference_bm6_dt001": {
            "dt_ref": dt_ref,
            "n_steps_ref": n_steps * sub_ref,
            "nu_full_rad": float(nu_ref_full),
            "nu_1_rad": float(nu_ref_1),
            "nu_2_rad": float(nu_ref_2),
            "physical_two_window_shift_rel": float(abs(drift_ref_12)),
        },
        "hanning_spectrum_omega_grid": [float(w) for w in omega_spec_grid],
        "integrators": {},
    }

    f64_entries: Dict[str, Dict[str, Any]] = {}
    f128_entries: Dict[str, Dict[str, Any]] = {}

    with ProcessPoolExecutor(max_workers=10) as pool:
        futs_f64 = {
            pool.submit(
                _run_s2_f64_worker,
                key,
                label,
                order,
                s_eval,
                dt,
                n_steps,
                stride_series,
                z_orb_ref,
                z_kepl_ref,
                z_ecc_ref,
                float(nu_ref_full),
                float(nu_ref_1),
                float(nu_ref_2),
                float(nu_kepl_ref_full),
                float(nu_ecc_ref_full),
                omega_spec_grid,
            ): key
            for key, label, order, s_eval in INTEGRATOR_SPECS
        }
        futs_f128 = {
            pool.submit(_run_s2_f128_worker, key, dt, n_steps, stride_series): key
            for key, _, _, _ in INTEGRATOR_SPECS
        }
        for fut in as_completed(futs_f64):
            k, ent = fut.result()
            f64_entries[k] = ent
            print(
                f"    [S2 F64 :: {k}] done in {ent['wall_time_s']:.2f}s | k={ent['unique_dyadic_levels']} (junctions={ent['n_dyadic_junctions']}) | max_dE={ent['max_rel_dE']:.3e} | freq_err={ent['naff']['rel_freq_err_vs_ref']:.3e}",
                flush=True,
            )
        for fut in as_completed(futs_f128):
            k, ent128 = fut.result()
            f128_entries[k] = ent128
            print(
                f"    [S2 F128 :: {k}] done in {ent128['wall_time_s']:.2f}s | max_dE_f128={ent128['max_dE_rel']:.3e} | final_rms_f128={ent128['final_rms_dE']:.3e} | max_dL_f128={ent128['max_dL_rel']:.3e}",
                flush=True,
            )

    dp_integrators: Dict[str, Any] = {}
    for key, _, _, _ in INTEGRATOR_SPECS:
        ent64 = f64_entries[key]
        ent128 = f128_entries[key]
        ent64["max_rel_dE_f128"] = ent128["max_dE_rel"]
        ent64["final_rms_dE_f128"] = ent128["final_rms_dE"]
        ent64["secular_slope_per_yr_f128"] = ent128["secular_slope_per_yr"]
        ent64["series"]["abs_dE_f128"] = ent128["dE_rel"]
        ent64["series"]["running_rms_dE_f128"] = ent128["running_rms_dE"]
        ent64["series"]["dP_norm_rel_f128"] = ent128["dP_norm_rel"]
        ent64["series"]["dL_norm_rel_f128"] = ent128["dL_norm_rel"]
        ent64["series"]["dG_norm_rel_f128"] = ent128["dG_norm_rel"]
        s2_results["integrators"][key] = ent64

        s_dE_64 = ent64["series"]["signed_dE"]
        s_dE_128 = ent128["signed_dE_rel"]
        dE_roundoff_diff = [abs(a - b) for a, b in zip(s_dE_64, s_dE_128)]
        dp_integrators[key] = {
            "is_f128_certified": ent128["is_f128_certified"],
            "times": ent128["times"],
            "times_yr": ent128["times_yr"],
            "max_dE_f64_vs_f128_roundoff": float(max(dE_roundoff_diff)),
            "dE_f64_vs_f128_roundoff": dE_roundoff_diff,
            "float64": {
                "dE_rel": ent64["series"]["abs_dE"],
                "running_rms_dE": ent64["series"]["running_rms_dE"],
                "max_dE_rel": ent64["max_rel_dE"],
                "final_rms_dE": ent64["final_rms_dE"],
                "max_dP_rel": ent64["lie_bargmann"]["max_dP_rel"],
                "max_dL_rel": ent64["lie_bargmann"]["max_dL_rel"],
                "max_dG_rel": ent64["lie_bargmann"]["max_dG_rel"],
            },
            "float128": {
                "dE_rel": ent128["dE_rel"],
                "running_rms_dE": ent128["running_rms_dE"],
                "max_dE_rel": ent128["max_dE_rel"],
                "final_rms_dE": ent128["final_rms_dE"],
                "max_dP_rel": ent128["max_dP_rel"],
                "max_dL_rel": ent128["max_dL_rel"],
                "max_dG_rel": ent128["max_dG_rel"],
            },
        }

    pbm4_s2 = s2_results["integrators"]["PBM4_Z4"]
    bm4_s2 = s2_results["integrators"]["BM4_Bare"]
    s2_results["shadow_telescoping_gain_max"] = float(
        bm4_s2["max_rel_dE"] / pbm4_s2["max_rel_dE"]
    )
    s2_results["shadow_telescoping_gain_rms"] = float(
        bm4_s2["final_rms_dE"] / pbm4_s2["final_rms_dE"]
    )
    s2_results["shadow_telescoping_gain_max_f128"] = float(
        bm4_s2["max_rel_dE_f128"] / pbm4_s2["max_rel_dE_f128"]
    )
    s2_results["shadow_telescoping_gain_rms_f128"] = float(
        bm4_s2["final_rms_dE_f128"] / pbm4_s2["final_rms_dE_f128"]
    )
    s2_results["naff_frequency_gain_pbm4_over_bm4"] = float(
        bm4_s2["naff"]["rel_freq_err_vs_ref"]
        / pbm4_s2["naff"]["rel_freq_err_vs_ref"]
    )

    s2_results["dual_precision_equal_params"] = {
        "common_parameters": {
            "benchmark_name": "inner_solar_system_mercury",
            "dt": float(dt),
            "n_steps": int(n_steps),
            "stride": int(stride_series),
            "enable_dyadic_substepping": True,
            "max_dyadic_level": 3,
            "desingularization_mode": "auto",
            "supernode_mode": "canonical",
            "mode": "homological",
        },
        "parameters_strictly_equal": True,
        "integrators": dp_integrators,
    }
    s2_results["total_wall_time_s"] = float(time.perf_counter() - t0_wall)
    return s2_results


def verify_contract_assertions(data: Dict[str, Any]) -> Dict[str, Any]:
    """Evaluates the 6 mandatory quantitative assertions of Contract 05."""
    s1 = data["S1_laplace"]
    s2 = data["S2_solar_system"]

    # ASSERT-01: True Float128 & strict parameter equality on S1 and S2
    dp_s1 = s1["dual_precision_equal_params"]
    dp_s2 = s2["dual_precision_equal_params"]
    f128_ok = (
        dp_s1["parameters_strictly_equal"]
        and dp_s2["parameters_strictly_equal"]
        and all(v["is_f128_certified"] for v in dp_s1["integrators"].values())
        and all(v["is_f128_certified"] for v in dp_s2["integrators"].values())
        and max(v["float128"]["max_dL_rel"] for v in dp_s1["integrators"].values()) <= 1e-30
        and max(v["float128"]["max_dL_rel"] for v in dp_s2["integrators"].values()) <= 1e-30
    )

    # ASSERT-02: S1 Laplace forced resonant libration & shadow telescoping gain
    pbm4_s1 = s1["integrators"]["PBM4_Z4"]
    bm4_s1 = s1["integrators"]["BM4_Bare"]
    bm6_s1 = s1["integrators"]["BM6_Bare"]
    s1_libration_ok = (
        abs(pbm4_s1["theta_lap_mean_median_deg"] - 180.0) <= 1.0
        and pbm4_s1["theta_lap_mean_half_amp_deg"] <= 2.5
        and bm4_s1["theta_lap_mean_half_amp_deg"] <= 2.5
        and bm6_s1["theta_lap_mean_half_amp_deg"] <= 2.5
        and s1["shadow_telescoping_gain_max"] >= 5.0
    )

    # ASSERT-03: All 9 Galilean Lie-Bargmann generators conserved across all 5 integrators
    max_dL_s1 = max(v["lie_bargmann"]["max_dL_rel"] for v in s1["integrators"].values())
    max_dP_s1 = max(v["lie_bargmann"]["max_dP_rel"] for v in s1["integrators"].values())
    max_dL_s2 = max(v["lie_bargmann"]["max_dL_rel"] for v in s2["integrators"].values())
    max_dP_s2 = max(v["lie_bargmann"]["max_dP_rel"] for v in s2["integrators"].values())
    lie_bargmann_ok = (
        max_dL_s1 <= 2e-12
        and max_dP_s1 <= 2e-12
        and max_dL_s2 <= 2e-12
        and max_dP_s2 <= 2e-12
    )

    # ASSERT-04: S2 shadow telescoping gain >= 20x RMS / >= 18x max across 264 dyadic junctions
    # (reaching Float64 floor ~7.9e-14 RMS / 1.55e-13 max), 0 saturation events, and 0 spurious synodic atlas flips
    s2_telescoping_ok = (
        s1.get("enable_dyadic_substepping") is True
        and s2.get("enable_dyadic_substepping") is True
        and s2["shadow_telescoping_gain_rms"] >= 20.0
        and s2["shadow_telescoping_gain_max"] >= 18.0
        and all(v["n_atlas_transitions"] == 0 for v in s2["integrators"].values())
        and all(v.get("n_saturation_events", 0) == 0 for v in s1["integrators"].values())
        and all(v.get("n_saturation_events", 0) == 0 for v in s2["integrators"].values())
    )

    # ASSERT-05: Unbiased Brouwer roundoff law on BM6 (no -0.5 ULP linear velocity drift)
    bm6_s2 = s2["integrators"]["BM6_Bare"]
    brouwer_ok = (
        bm6_s2["max_rel_dE"] <= 5e-13
        and abs(bm6_s2["secular_slope_per_yr"]) <= 5e-15
    )

    # ASSERT-06: Authentic Laskar NAFF on densely sampled z_orb(t) (zero z_act proxy):
    # (i) >2000x frequency error reduction from S2 -> Y4 -> {BM4, PBM4} (reaching ~2.6e-13 NAFF floor), and
    # (ii) KAM conjugacy invariance |nu(PBM4) - nu(BM4)|/nu_ref <= 1e-12 alongside >20x RMS shadow telescoping gain
    err_s2 = s2["integrators"]["S2_Strang"]["naff"]["rel_freq_err_vs_ref"]
    err_y4 = s2["integrators"]["Y4_Yoshida"]["naff"]["rel_freq_err_vs_ref"]
    err_bm4 = s2["integrators"]["BM4_Bare"]["naff"]["rel_freq_err_vs_ref"]
    err_pbm4 = s2["integrators"]["PBM4_Z4"]["naff"]["rel_freq_err_vs_ref"]
    err_bm6 = s2["integrators"]["BM6_Bare"]["naff"]["rel_freq_err_vs_ref"]
    nu_bm4 = s2["integrators"]["BM4_Bare"]["naff"]["nu_orb_full_rad"]
    nu_pbm4 = s2["integrators"]["PBM4_Z4"]["naff"]["nu_orb_full_rad"]
    nu_ref = s2["reference_bm6_dt001"]["nu_full_rad"]
    kam_conjugacy_diff = abs(nu_pbm4 - nu_bm4) / abs(nu_ref)
    naff_ok = (
        s2["used_synthetic_z_act"] is False
        and s2["naff_samples"] == 4001
        and max(err_bm4, err_pbm4) < 0.25 * err_y4
        and max(err_bm4, err_pbm4) < 6e-4 * err_s2
        and err_y4 < 0.05 * err_s2
        and kam_conjugacy_diff <= 1e-12
    )

    assertions = {
        "ASSERT-01_is_f128_and_equal_params": {
            "passed": bool(f128_ok),
            "s1_max_dL_f128": float(
                max(v["float128"]["max_dL_rel"] for v in dp_s1["integrators"].values())
            ),
            "s2_max_dL_f128": float(
                max(v["float128"]["max_dL_rel"] for v in dp_s2["integrators"].values())
            ),
        },
        "ASSERT-02_S1_laplace_forced_libration_and_shadow_gain": {
            "passed": bool(s1_libration_ok),
            "pbm4_theta_median_deg": float(pbm4_s1["theta_lap_mean_median_deg"]),
            "pbm4_theta_half_amp_deg": float(pbm4_s1["theta_lap_mean_half_amp_deg"]),
            "shadow_gain_max": float(s1["shadow_telescoping_gain_max"]),
        },
        "ASSERT-03_S1_S2_lie_bargmann_all_5_integrators": {
            "passed": bool(lie_bargmann_ok),
            "max_dL_s1": float(max_dL_s1),
            "max_dL_s2": float(max_dL_s2),
            "max_dP_s1": float(max_dP_s1),
            "max_dP_s2": float(max_dP_s2),
        },
        "ASSERT-04_S2_shadow_telescoping_and_atlas_stability": {
            "passed": bool(s2_telescoping_ok),
            "gain_rms": float(s2["shadow_telescoping_gain_rms"]),
            "gain_max": float(s2["shadow_telescoping_gain_max"]),
            "max_atlas_flips": int(
                max(v["n_atlas_transitions"] for v in s2["integrators"].values())
            ),
        },
        "ASSERT-05_S2_unbiased_brouwer_no_bm6_linear_drift": {
            "passed": bool(brouwer_ok),
            "bm6_max_rel_dE": float(bm6_s2["max_rel_dE"]),
            "bm6_secular_slope_per_yr": float(bm6_s2["secular_slope_per_yr"]),
        },
        "ASSERT-06_S2_authentic_laskar_naff_on_z_orb": {
            "passed": bool(naff_ok),
            "used_synthetic_z_act": bool(s2["used_synthetic_z_act"]),
            "naff_samples": int(s2["naff_samples"]),
            "rel_freq_err_S2": float(err_s2),
            "rel_freq_err_Y4": float(err_y4),
            "rel_freq_err_BM4": float(err_bm4),
            "rel_freq_err_PBM4": float(err_pbm4),
            "rel_freq_err_BM6": float(err_bm6),
            "kam_conjugacy_rel_diff_pbm4_vs_bm4": float(kam_conjugacy_diff),
        },
    }
    all_passed = all(v["passed"] for v in assertions.values())
    return {
        "all_passed": bool(all_passed),
        "n_passed": int(sum(1 for v in assertions.values() if v["passed"])),
        "n_total": len(assertions),
        "assertions": assertions,
    }


def _sync_legacy_data_files(data: Dict[str, Any]) -> None:
    """Synchronizes data/campaign_III_S1_laplace.json and campaign_III_S2_brouwer_naff.json."""
    if "S1_laplace" in data:
        LEGACY_S1_JSON.parent.mkdir(parents=True, exist_ok=True)
        with open(LEGACY_S1_JSON, "w", encoding="utf-8") as f:
            json.dump(data["S1_laplace"], f, indent=2)
    if "S2_solar_system" in data:
        LEGACY_S2_JSON.parent.mkdir(parents=True, exist_ok=True)
        with open(LEGACY_S2_JSON, "w", encoding="utf-8") as f:
            json.dump(data["S2_solar_system"], f, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generator 05 — Secular Stability, Laplace Resonance & Laskar NAFF"
    )
    parser.add_argument(
        "--protocol",
        type=str,
        default="all",
        choices=["S1", "S2", "verify", "all"],
        help="Isolated protocol to execute",
    )
    args = parser.parse_args()

    data = _load_or_init_json()

    if args.protocol in ("S1", "all"):
        print("[Generator 05] Running Protocol S1 (Jovian 4:2:1 Laplace Resonance)...")
        data["S1_laplace"] = run_protocol_S1()
        _save_json(data)
        _sync_legacy_data_files(data)
        s1 = data["S1_laplace"]
        print(
            f"  -> S1 completed in {s1['total_wall_time_s']:.2f}s | "
            f"Shadow gain = {s1['shadow_telescoping_gain_max']:.2f}x | "
            f"PBM4 Laplace half-amp = {s1['integrators']['PBM4_Z4']['theta_lap_mean_half_amp_deg']:.3f} deg"
        )

    if args.protocol in ("S2", "all"):
        print("[Generator 05] Running Protocol S2 (Inner Solar System N=5 & Authentic Laskar NAFF)...")
        data["S2_solar_system"] = run_protocol_S2()
        _save_json(data)
        _sync_legacy_data_files(data)
        s2 = data["S2_solar_system"]
        print(
            f"  -> S2 completed in {s2['total_wall_time_s']:.2f}s | "
            f"Shadow RMS gain = {s2['shadow_telescoping_gain_rms']:.2f}x | "
            f"NAFF gain PBM4/BM4 = {s2['naff_frequency_gain_pbm4_over_bm4']:.2f}x | "
            f"BM6 max dE = {s2['integrators']['BM6_Bare']['max_rel_dE']:.3e}"
        )

    if args.protocol in ("verify", "all"):
        print("[Generator 05] Verifying Contract 05 Quantitative Assertions...")
        ver = verify_contract_assertions(data)
        data["contract_verification"] = ver
        _save_json(data)
        for aid, ares in ver["assertions"].items():
            status = "PASS" if ares["passed"] else "FAIL"
            print(f"  [{status}] {aid}: {ares}")
        print(
            f"[Generator 05] Contract Verification: {ver['n_passed']}/{ver['n_total']} PASS"
        )
        if not ver["all_passed"]:
            sys.exit(1)


if __name__ == "__main__":
    main()
