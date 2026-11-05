#!/usr/bin/env python3
"""
Cycle VI - Iteration 75: Head-to-Head Benchmark B4
Definitive Reversibility of the Chaotic Pythagorean 3-4-5 Problem Beyond the Float64 Lyapunov Horizon
References:
  - Burrau (1913), Astron. Nachr. 195, 113 (Pythagorean 3-4-5 three-body problem)
  - Szebehely & Peters (1967), AJ 72, 876 (Complete solution & final ejection at t ~ 60)
  - Boekholt & Portegies Zwart (2015), Comput. Astrophys. Cosmol. 2:2 (Brutus arbitrary-precision reversibility)
  - Rein & Spiegel (2015), MNRAS 446, 1424 (REBOUND IAS15)
  - Rein et al. (2019), MNRAS 485, 5490 (REBOUND MERCURIUS)
  - Lu, Hernandez & Rein (2024), MNRAS 533, 3708 (REBOUND TRACE)

Design Principles:
  1. Direct head-to-head execution of REBOUND 5.2.2 (IAS15, TRACE BS, TRACE IAS15, MERCURIUS,
     WHFast c=17, WHFast c=0) against our Unsoftened Normal-Cone Symplectic Splitting Solver
     (Normal-Cone BM6 & Normal-Cone PBM4) in both IEEE-754 Float64 and 35-digit Float128.
  2. Dynamical Jacobi chart flips (AtlasManager + non-commutative covariant shadow gauge
     connection T_{A->B}^{(PBM4)}) always active on all Normal-Cone solvers.
  3. Explicit separation of Category 1 (Adaptive / Multi-Scale Solvers) and
     Category 2 (Fixed-Step Solvers k=0 / Fixed dt).
  4. Zero arbitrary floor/clamping on any physical observable or error metric.
"""

from __future__ import annotations
import json
import math
import multiprocessing as mp_proc
import os
import sys
import time
from typing import Any, Dict, List, Tuple

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import numpy as np
import mpmath as mp
try:
    import rebound
except ImportError:
    rebound = None  # type: ignore

from normal_cone_nbody_float64.core.state import (
    ExtendedPhaseState as ExtendedPhaseStateF64,
)
from normal_cone_nbody_float64.core.invariants import (
    compute_invariants as compute_invariants_f64,
)
from normal_cone_nbody_float64.system.nbody_system import (
    NBodySystem as NBodySystemF64,
)
from normal_cone_nbody_float64.integrators.blanes_moan_4 import (
    BlanesMoan4Integrator as BlanesMoan4IntegratorF64,
)
from normal_cone_nbody_float64.integrators.blanes_moan_6 import (
    BlanesMoan6Integrator as BlanesMoan6IntegratorF64,
)
from normal_cone_nbody_float64.integrators.processed import (
    ProcessedSymplecticIntegrator as ProcessedSymplecticIntegratorF64,
    CorrectorC4 as CorrectorC4F64,
)

from normal_cone_nbody_float128.core.constants import (
    to_f128,
)
from normal_cone_nbody_float128.core.state import (
    ExtendedPhaseState as ExtendedPhaseStateF128,
)
from normal_cone_nbody_float128.core.invariants import (
    compute_invariants as compute_invariants_f128,
    compute_truncation_residual as compute_truncation_residual_f128,
)
from normal_cone_nbody_float128.system.nbody_system import (
    NBodySystem as NBodySystemF128,
)
from normal_cone_nbody_float128.integrators.blanes_moan_4 import (
    BlanesMoan4Integrator as BlanesMoan4IntegratorF128,
)
from normal_cone_nbody_float128.integrators.blanes_moan_6 import (
    BlanesMoan6Integrator as BlanesMoan6IntegratorF128,
)
from normal_cone_nbody_float128.integrators.processed import (
    ProcessedSymplecticIntegrator as ProcessedSymplecticIntegratorF128,
    CorrectorC4 as CorrectorC4F128,
)


def _to_f128_array(arr: Any) -> np.ndarray:
    return to_f128(np.asarray(arr, dtype=object))


# ============================================================================
# 1. Exact Initial Conditions & Common Observables (Burrau 1913, 3-4-5)
# ============================================================================

def build_b4_initial_conditions_f64(G: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns exact barycentric initial conditions (masses, q0, v0, p0) in Float64
    for Burrau's (1913) / Szebehely & Peters' (1967) Pythagorean 3-4-5 problem:
      m = (3, 4, 5), q1 = (1, 3, 0), q2 = (-2, -1, 0), q3 = (1, -1, 0), v_i = 0.
    """
    masses = np.array([3.0, 4.0, 5.0], dtype=np.float64)
    q0 = np.array([
        [1.0, 3.0, 0.0],
        [-2.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
    ], dtype=np.float64)
    v0 = np.zeros((3, 3), dtype=np.float64)
    M_tot = float(np.sum(masses))
    q0 -= np.sum(masses[:, None] * q0, axis=0) / M_tot
    p0 = masses[:, None] * v0
    return masses, q0, v0, p0


def build_b4_initial_state_f128() -> ExtendedPhaseStateF128:
    """
    Returns exact 35-digit Float128 barycentric ExtendedPhaseState for Burrau's 3-4-5 problem.
    Center of mass is exactly R_cm = (-1/12, 0, 0).
    """
    mp.mp.dps = 35
    shift_x = mp.mpf(1) / mp.mpf(12)
    masses = _to_f128_array([mp.mpf(3), mp.mpf(4), mp.mpf(5)])
    q0 = _to_f128_array([
        [mp.mpf(1) + shift_x, mp.mpf(3), mp.mpf(0)],
        [mp.mpf(-2) + shift_x, mp.mpf(-1), mp.mpf(0)],
        [mp.mpf(1) + shift_x, mp.mpf(-1), mp.mpf(0)],
    ])
    p0 = _to_f128_array([[mp.mpf(0)] * 3 for _ in range(3)])
    st = ExtendedPhaseStateF128(q=q0, p=p0, masses=masses, t=to_f128(0))
    return st


def compute_common_observables_f64(
    masses: np.ndarray,
    q: np.ndarray,
    v: np.ndarray,
    t: float = 0.0,
    G: float = 1.0,
) -> Dict[str, Any]:
    """Computes Hamiltonian H, Lie-Bargmann Noether invariants (P, L, G_boost), and min pairwise distance."""
    p = masses[:, None] * v
    st = ExtendedPhaseStateF64(q=q.copy(), p=p.copy(), masses=masses.copy(), t=float(t))
    inv = compute_invariants_f64(st, G=float(G))
    d01 = float(np.linalg.norm(q[0] - q[1]))
    d02 = float(np.linalg.norm(q[0] - q[2]))
    d12 = float(np.linalg.norm(q[1] - q[2]))
    r_min = min(d01, d02, d12)
    if r_min == d01:
        closest_pair = [0, 1]
    elif r_min == d02:
        closest_pair = [0, 2]
    else:
        closest_pair = [1, 2]
    return {
        "H": float(inv.H),
        "P": np.asarray(inv.P, dtype=np.float64).copy(),
        "L": np.asarray(inv.L, dtype=np.float64).copy(),
        "G_boost": np.asarray(inv.G, dtype=np.float64).copy(),
        "r_min": float(r_min),
        "closest_pair": closest_pair,
    }


# ============================================================================
# 2. REBOUND 5.2.2 Runner (IAS15, TRACE BS, TRACE IAS15, MERCURIUS, WHFast)
# ============================================================================

def _configure_rebound_integrator(
    sim: rebound.Simulation,
    solver_key: str,
    step_by_step_monitoring: bool = False,
) -> str:
    if solver_key == "rebound_ias15":
        sim.integrator = "ias15"
        sim.integrator.epsilon = 1e-9
        sim.integrator.adaptive_mode = "PRS23"
        return "adaptive_multiscale"
    elif solver_key == "rebound_trace_bs":
        sim.integrator = "trace"
        sim.integrator.r_crit_hill = 3.0
        sim.integrator.peri_crit_eta = 1.0
        sim.integrator.peri_mode = "FULL_BS"
        return "adaptive_multiscale"
    elif solver_key == "rebound_trace_ias15":
        sim.integrator = "trace"
        sim.integrator.r_crit_hill = 3.0
        sim.integrator.peri_crit_eta = 1.0
        sim.integrator.peri_mode = "FULL_IAS15"
        return "adaptive_multiscale"
    elif solver_key == "rebound_mercurius":
        sim.integrator = "mercurius"
        sim.integrator.r_crit_hill = 3.0
        sim.integrator.safe_mode = 1 if step_by_step_monitoring else 0
        return "adaptive_multiscale"
    elif solver_key == "rebound_whfast_c17":
        sim.integrator = "whfast"
        sim.integrator.corrector = 17
        sim.integrator.safe_mode = 0
        sim.integrator.keep_unsynchronized = 1
        return "fixed_step"
    elif solver_key == "rebound_whfast_c0":
        sim.integrator = "whfast"
        sim.integrator.corrector = 0
        sim.integrator.safe_mode = 0
        sim.integrator.keep_unsynchronized = 1
        return "fixed_step"
    else:
        raise ValueError(f"Unknown REBOUND solver_key: {solver_key}")


def run_rebound_forward_b4a(
    solver_key: str,
    t_end: float = 70.0,
    dt: float = 0.04,
    G: float = 1.0,
) -> Dict[str, Any]:
    masses, q0, v0, _ = build_b4_initial_conditions_f64(G=G)
    obs0 = compute_common_observables_f64(masses, q0, v0, t=0.0, G=G)
    H0 = obs0["H"]
    n_steps = int(round(t_end / dt))

    # Pure timing run
    sim_pure = rebound.Simulation()
    sim_pure.G = float(G)
    for i in range(3):
        sim_pure.add(
            m=float(masses[i]),
            x=float(q0[i, 0]), y=float(q0[i, 1]), z=float(q0[i, 2]),
            vx=float(v0[i, 0]), vy=float(v0[i, 1]), vz=float(v0[i, 2]),
        )
    category = _configure_rebound_integrator(sim_pure, solver_key)
    sim_pure.dt = float(dt)
    t0_pure = time.perf_counter()
    if solver_key == "rebound_ias15":
        sim_pure.integrate(float(t_end), exact_finish_time=1)
    else:
        sim_pure.steps(n_steps)
        sim_pure.synchronize()
    wall_time_pure_s = float(time.perf_counter() - t0_pure)

    # Monitored trajectory run
    sim = rebound.Simulation()
    sim.G = float(G)
    for i in range(3):
        sim.add(
            m=float(masses[i]),
            x=float(q0[i, 0]), y=float(q0[i, 1]), z=float(q0[i, 2]),
            vx=float(v0[i, 0]), vy=float(v0[i, 1]), vz=float(v0[i, 2]),
        )
    _configure_rebound_integrator(sim, solver_key, step_by_step_monitoring=True)
    sim.dt = float(dt)

    series_t: List[float] = [0.0]
    series_rel_dE: List[float] = [0.0]
    series_r_min: List[float] = [float(obs0["r_min"])]
    series_dP: List[float] = [0.0]
    series_dL: List[float] = [0.0]
    series_dG: List[float] = [0.0]
    max_rel_dE = 0.0
    min_r_overall = float(obs0["r_min"])

    for k in range(1, n_steps + 1):
        t_target = float(k * dt)
        if solver_key == "rebound_ias15":
            sim.integrate(t_target, exact_finish_time=1)
        else:
            sim.steps(1)
            sim.synchronize()
        q_cur = np.array([[p.x, p.y, p.z] for p in sim.particles], dtype=np.float64)
        v_cur = np.array([[p.vx, p.vy, p.vz] for p in sim.particles], dtype=np.float64)
        obs = compute_common_observables_f64(masses, q_cur, v_cur, t=float(sim.t), G=G)
        rel_dE = abs(obs["H"] - H0) / abs(H0)
        if rel_dE > max_rel_dE:
            max_rel_dE = rel_dE
        if obs["r_min"] < min_r_overall:
            min_r_overall = obs["r_min"]
        series_t.append(float(sim.t))
        series_rel_dE.append(float(rel_dE))
        series_r_min.append(float(obs["r_min"]))
        series_dP.append(float(np.linalg.norm(obs["P"] - obs0["P"])))
        series_dL.append(float(np.linalg.norm(obs["L"] - obs0["L"])))
        series_dG.append(float(np.linalg.norm(obs["G_boost"] - obs0["G_boost"])))

    q_end = np.array([[p.x, p.y, p.z] for p in sim.particles], dtype=np.float64)
    v_end = np.array([[p.vx, p.vy, p.vz] for p in sim.particles], dtype=np.float64)
    obs_end = compute_common_observables_f64(masses, q_end, v_end, t=float(sim.t), G=G)
    end_rel_dE = abs(obs_end["H"] - H0) / abs(H0)
    dP_norm = float(np.linalg.norm(obs_end["P"] - obs0["P"]))
    dL_norm = float(np.linalg.norm(obs_end["L"] - obs0["L"]))
    dG_norm = float(np.linalg.norm(obs_end["G_boost"] - obs0["G_boost"]))

    return {
        "solver": solver_key,
        "category": category,
        "precision": "Float64",
        "t_end": float(t_end),
        "dt": float(dt),
        "n_steps": int(n_steps),
        "max_rel_dE": float(max_rel_dE),
        "end_rel_dE": float(end_rel_dE),
        "min_r": float(min_r_overall),
        "dP_norm": float(dP_norm),
        "dL_norm": float(dL_norm),
        "dG_norm": float(dG_norm),
        "wall_time_pure_s": float(wall_time_pure_s),
        "ms_per_time_unit": float(wall_time_pure_s * 1000.0 / max(1e-9, float(t_end))),
        "series_t": series_t,
        "series_rel_dE": series_rel_dE,
        "series_r_min": series_r_min,
        "series_dP": series_dP,
        "series_dL": series_dL,
        "series_dG": series_dG,
    }


def run_rebound_reversibility_b4b(
    solver_key: str,
    t_turnaround: float,
    dt: float = 0.04,
    G: float = 1.0,
) -> Dict[str, Any]:
    masses, q0, v0, p0 = build_b4_initial_conditions_f64(G=G)
    obs0 = compute_common_observables_f64(masses, q0, v0, t=0.0, G=G)
    H0 = obs0["H"]
    n_steps = int(round(t_turnaround / dt))

    sim = rebound.Simulation()
    sim.G = float(G)
    for i in range(3):
        sim.add(
            m=float(masses[i]),
            x=float(q0[i, 0]), y=float(q0[i, 1]), z=float(q0[i, 2]),
            vx=float(v0[i, 0]), vy=float(v0[i, 1]), vz=float(v0[i, 2]),
        )
    category = _configure_rebound_integrator(sim, solver_key)
    sim.dt = float(dt)

    t0_wall = time.perf_counter()
    if solver_key == "rebound_ias15":
        sim.integrate(float(t_turnaround), exact_finish_time=1)
    else:
        sim.steps(n_steps)
        sim.synchronize()

    q_turn = np.array([[p.x, p.y, p.z] for p in sim.particles], dtype=np.float64)
    v_turn = np.array([[p.vx, p.vy, p.vz] for p in sim.particles], dtype=np.float64)

    sim_rev = rebound.Simulation()
    sim_rev.G = float(G)
    for i in range(3):
        sim_rev.add(
            m=float(masses[i]),
            x=float(q_turn[i, 0]), y=float(q_turn[i, 1]), z=float(q_turn[i, 2]),
            vx=float(-v_turn[i, 0]), vy=float(-v_turn[i, 1]), vz=float(-v_turn[i, 2]),
        )
    _configure_rebound_integrator(sim_rev, solver_key)
    sim_rev.dt = float(dt)
    if solver_key == "rebound_ias15":
        sim_rev.integrate(float(t_turnaround), exact_finish_time=1)
    else:
        sim_rev.steps(n_steps)
        sim_rev.synchronize()
    wall_time_s = float(time.perf_counter() - t0_wall)

    q_rev = np.array([[p.x, p.y, p.z] for p in sim_rev.particles], dtype=np.float64)
    v_rev = np.array([[-p.vx, -p.vy, -p.vz] for p in sim_rev.particles], dtype=np.float64)
    p_rev = masses[:, None] * v_rev
    obs_rev = compute_common_observables_f64(masses, q_rev, v_rev, t=0.0, G=G)

    rev_dq_norm = float(np.linalg.norm(q_rev - q0))
    rev_dp_norm = float(np.linalg.norm(p_rev - p0))
    rev_rel_dE = float(abs(obs_rev["H"] - H0) / abs(H0))

    return {
        "solver": solver_key,
        "category": category,
        "precision": "Float64",
        "t_turnaround": float(t_turnaround),
        "rev_dq_norm": rev_dq_norm,
        "rev_dp_norm": rev_dp_norm,
        "rev_rel_dE": rev_rel_dE,
        "wall_time_s": wall_time_s,
    }


# ============================================================================
# 3. Normal-Cone Symplectic Splitting Solver Runner (Float64 & Float128)
# ============================================================================

def _build_normal_cone_f64(solver_key: str, G: float = 1.0) -> Tuple[Any, bool, bool, str]:
    if solver_key == "normal_cone_bm6_f64_dyadic":
        return BlanesMoan6IntegratorF64(G=G), True, False, "adaptive_multiscale"
    elif solver_key == "normal_cone_pbm4_f64_dyadic":
        return (
            ProcessedSymplecticIntegratorF64(
                BlanesMoan4IntegratorF64(G=G),
                CorrectorC4F64(G=G, mode="homological"),
            ),
            True,
            True,
            "adaptive_multiscale",
        )
    elif solver_key == "normal_cone_bm6_f64_k0":
        return BlanesMoan6IntegratorF64(G=G), False, False, "fixed_step"
    elif solver_key == "normal_cone_pbm4_f64_k0":
        return (
            ProcessedSymplecticIntegratorF64(
                BlanesMoan4IntegratorF64(G=G),
                CorrectorC4F64(G=G, mode="homological"),
            ),
            False,
            True,
            "fixed_step",
        )
    else:
        raise ValueError(f"Unknown Float64 Normal-Cone solver_key: {solver_key}")


def _build_normal_cone_f128(solver_key: str, G: float = 1.0) -> Tuple[Any, bool, bool, str]:
    mp.mp.dps = 35
    if solver_key == "normal_cone_bm6_f128_dyadic":
        return BlanesMoan6IntegratorF128(G=G), True, False, "adaptive_multiscale"
    elif solver_key == "normal_cone_pbm4_f128_dyadic":
        return (
            ProcessedSymplecticIntegratorF128(
                BlanesMoan4IntegratorF128(G=G),
                CorrectorC4F128(G=G, mode="homological"),
            ),
            True,
            True,
            "adaptive_multiscale",
        )
    else:
        raise ValueError(f"Unknown Float128 Normal-Cone solver_key: {solver_key}")


def run_normal_cone_forward_b4a_f64(
    solver_key: str,
    t_end: float = 70.0,
    dt: float = 0.04,
    max_dyadic_level: int = 16,
    G: float = 1.0,
    atlas_gamma: float = 1.0,
) -> Dict[str, Any]:
    masses, q0, v0, p0 = build_b4_initial_conditions_f64(G=G)
    obs0 = compute_common_observables_f64(masses, q0, v0, t=0.0, G=G)
    H0 = obs0["H"]
    n_steps = int(round(t_end / dt))

    integ, enable_dyad, in_sh, category = _build_normal_cone_f64(solver_key, G=G)
    st0 = ExtendedPhaseStateF64(q=q0.copy(), p=p0.copy(), masses=masses.copy(), t=0.0)
    sys_nb = NBodySystemF64(
        state=st0,
        integrator=integ,
        dt_macro=float(dt),
        max_dyadic_level=int(max_dyadic_level),
        in_shadow_mode=in_sh,
        enable_dyadic_substepping=enable_dyad,
        use_optimized_kernels=True,
        atlas_gamma=float(atlas_gamma),
    )
    if enable_dyad:
        sys_nb.scheduler.warn_on_saturation = True

    series_t: List[float] = [0.0]
    series_rel_dE: List[float] = [0.0]
    series_r_min: List[float] = [float(obs0["r_min"])]
    series_dyadic_k: List[int] = [0]
    series_leaf_pair: List[List[int]] = [
        list(sys_nb.tree.get_leaf_pair_indices()[0]) if sys_nb.tree.get_leaf_pair_indices() else [1, 2]
    ]
    series_dP: List[float] = [0.0]
    series_dL: List[float] = [0.0]
    series_dG: List[float] = [0.0]

    t0_wall = time.perf_counter()
    for _ in range(n_steps):
        sys_nb.step(float(dt))
        st_phys = sys_nb.observe_physical_state()
        q_cur = np.asarray(st_phys.q, dtype=np.float64)
        v_cur = np.asarray(st_phys.p, dtype=np.float64) / masses[:, None]
        obs = compute_common_observables_f64(masses, q_cur, v_cur, t=float(st_phys.t), G=G)
        rel_dE = abs(obs["H"] - H0) / abs(H0)
        series_t.append(float(st_phys.t))
        series_rel_dE.append(float(rel_dE))
        series_r_min.append(float(obs["r_min"]))
        snap_last = sys_nb.history.snapshots[-1]
        series_dyadic_k.append(int(snap_last.dyadic_level))
        lp = sys_nb.tree.get_leaf_pair_indices()
        series_leaf_pair.append(list(lp[0]) if lp else [1, 2])
        series_dP.append(float(np.linalg.norm(obs["P"] - obs0["P"])))
        series_dL.append(float(np.linalg.norm(obs["L"] - obs0["L"])))
        series_dG.append(float(np.linalg.norm(obs["G_boost"] - obs0["G_boost"])))
    wall_time_pure_s = float(time.perf_counter() - t0_wall)

    flips_info = [
        {
            "t": float(tr.t),
            "old_leaf_pairs": [list(p) for p in tr.old_leaf_pairs],
            "new_leaf_pairs": [list(p) for p in tr.new_leaf_pairs],
        }
        for tr in sys_nb.atlas_manager.transition_history
    ]

    return {
        "solver": solver_key,
        "category": category,
        "precision": "Float64",
        "t_end": float(t_end),
        "dt": float(dt),
        "n_steps": int(n_steps),
        "max_rel_dE": float(max(series_rel_dE)),
        "end_rel_dE": float(series_rel_dE[-1]),
        "min_r": float(min(series_r_min)),
        "max_dyadic_level": int(max(series_dyadic_k)),
        "saturation_events": list(sys_nb.scheduler.saturation_events),
        "n_atlas_flips": int(len(flips_info)),
        "atlas_flips": flips_info,
        "dP_norm": float(series_dP[-1]),
        "dL_norm": float(series_dL[-1]),
        "dG_norm": float(series_dG[-1]),
        "max_dP_norm": float(max(series_dP)),
        "max_dL_norm": float(max(series_dL)),
        "max_dG_norm": float(max(series_dG)),
        "wall_time_pure_s": float(wall_time_pure_s),
        "ms_per_time_unit": float(wall_time_pure_s * 1000.0 / max(1e-9, float(t_end))),
        "series_t": series_t,
        "series_rel_dE": series_rel_dE,
        "series_r_min": series_r_min,
        "series_dyadic_k": series_dyadic_k,
        "series_leaf_pair": series_leaf_pair,
        "series_dP": series_dP,
        "series_dL": series_dL,
        "series_dG": series_dG,
    }


def run_normal_cone_forward_b4a_f128(
    solver_key: str,
    t_end: float = 12.0,
    dt: float = 0.04,
    max_dyadic_level: int = 16,
    G: float = 1.0,
    atlas_gamma: float = 1.0,
) -> Dict[str, Any]:
    with mp.workdps(35):
        st0 = build_b4_initial_state_f128()
        inv0 = compute_invariants_f128(st0, G=to_f128(G))
        H0 = inv0.H
        n_steps = int(round(t_end / dt))
        dt_f128 = to_f128(str(dt))

        integ, enable_dyad, in_sh, category = _build_normal_cone_f128(solver_key, G=G)
        sys_nb = NBodySystemF128(
            state=st0.with_updates(),
            integrator=integ,
            G=to_f128(G),
            dt_macro=dt_f128,
            max_dyadic_level=int(max_dyadic_level),
            in_shadow_mode=in_sh,
            enable_dyadic_substepping=enable_dyad,
            atlas_gamma=float(atlas_gamma),
        )
        if enable_dyad:
            sys_nb.scheduler.warn_on_saturation = True

        q0_f64 = np.asarray(st0.q, dtype=np.float64)
        r0 = min(
            float(np.linalg.norm(q0_f64[0] - q0_f64[1])),
            float(np.linalg.norm(q0_f64[0] - q0_f64[2])),
            float(np.linalg.norm(q0_f64[1] - q0_f64[2])),
        )
        series_t: List[float] = [0.0]
        series_rel_dE: List[float] = [0.0]
        series_r_min: List[float] = [r0]
        series_dyadic_k: List[int] = [0]
        series_leaf_pair: List[List[int]] = [
            list(sys_nb.tree.get_leaf_pair_indices()[0]) if sys_nb.tree.get_leaf_pair_indices() else [1, 2]
        ]
        series_dP: List[float] = [0.0]
        series_dL: List[float] = [0.0]
        series_dG: List[float] = [0.0]

        t0_wall = time.perf_counter()
        for _ in range(n_steps):
            sys_nb.step(dt=dt_f128)
            st_phys = sys_nb.observe_physical_state()
            inv_cur = compute_invariants_f128(st_phys, G=to_f128(G))
            rel_dE = float(abs(inv_cur.H - H0) / abs(H0))
            res_vec = compute_truncation_residual_f128(inv_cur, inv0)
            q_f64 = np.asarray(st_phys.q, dtype=np.float64)
            r_min = min(
                float(np.linalg.norm(q_f64[0] - q_f64[1])),
                float(np.linalg.norm(q_f64[0] - q_f64[2])),
                float(np.linalg.norm(q_f64[1] - q_f64[2])),
            )
            series_t.append(float(st_phys.t))
            series_rel_dE.append(rel_dE)
            series_r_min.append(r_min)
            snap_last = sys_nb.history.snapshots[-1]
            series_dyadic_k.append(int(snap_last.dyadic_level))
            lp = sys_nb.tree.get_leaf_pair_indices()
            series_leaf_pair.append(list(lp[0]) if lp else [1, 2])
            series_dP.append(float(res_vec[0]))
            series_dL.append(float(res_vec[1]))
            series_dG.append(float(res_vec[2]))
        wall_time_pure_s = float(time.perf_counter() - t0_wall)

        flips_info = [
            {
                "t": float(tr.t),
                "old_leaf_pairs": [list(p) for p in tr.old_leaf_pairs],
                "new_leaf_pairs": [list(p) for p in tr.new_leaf_pairs],
            }
            for tr in sys_nb.atlas_manager.transition_history
        ]

        return {
            "solver": solver_key,
            "category": category,
            "precision": "Float128",
            "t_end": float(t_end),
            "dt": float(dt),
            "n_steps": int(n_steps),
            "max_rel_dE": float(max(series_rel_dE)),
            "end_rel_dE": float(series_rel_dE[-1]),
            "min_r": float(min(series_r_min)),
            "max_dyadic_level": int(max(series_dyadic_k)),
            "saturation_events": list(sys_nb.scheduler.saturation_events),
            "n_atlas_flips": int(len(flips_info)),
            "atlas_flips": flips_info,
            "dP_norm": float(series_dP[-1]),
            "dL_norm": float(series_dL[-1]),
            "dG_norm": float(series_dG[-1]),
            "max_dP_norm": float(max(series_dP)),
            "max_dL_norm": float(max(series_dL)),
            "max_dG_norm": float(max(series_dG)),
            "wall_time_pure_s": float(wall_time_pure_s),
            "ms_per_time_unit": float(wall_time_pure_s * 1000.0 / max(1e-9, float(t_end))),
            "series_t": series_t,
            "series_rel_dE": series_rel_dE,
            "series_r_min": series_r_min,
            "series_dyadic_k": series_dyadic_k,
            "series_leaf_pair": series_leaf_pair,
            "series_dP": series_dP,
            "series_dL": series_dL,
            "series_dG": series_dG,
        }


def run_normal_cone_reversibility_b4b_f64(
    solver_key: str,
    t_turnaround: float,
    dt: float = 0.04,
    max_dyadic_level: int = 16,
    G: float = 1.0,
    atlas_gamma: float = 1.0,
) -> Dict[str, Any]:
    masses, q0, v0, p0 = build_b4_initial_conditions_f64(G=G)
    obs0 = compute_common_observables_f64(masses, q0, v0, t=0.0, G=G)
    H0 = obs0["H"]
    n_steps = int(round(t_turnaround / dt))

    integ, enable_dyad, in_sh, category = _build_normal_cone_f64(solver_key, G=G)
    st0 = ExtendedPhaseStateF64(q=q0.copy(), p=p0.copy(), masses=masses.copy(), t=0.0)
    sys_nb = NBodySystemF64(
        state=st0,
        integrator=integ,
        dt_macro=float(dt),
        max_dyadic_level=int(max_dyadic_level),
        in_shadow_mode=in_sh,
        enable_dyadic_substepping=enable_dyad,
        use_optimized_kernels=True,
        atlas_gamma=float(atlas_gamma),
    )
    if enable_dyad:
        sys_nb.scheduler.warn_on_saturation = True

    t0_wall = time.perf_counter()
    for _ in range(n_steps):
        sys_nb.step(float(dt))
    n_flips_fwd = len(sys_nb.atlas_manager.transition_history)
    k_fwd = [int(s.dyadic_level) for s in sys_nb.history.snapshots[1 : 1 + n_steps]]

    for _ in range(n_steps):
        sys_nb.step(-float(dt))
    wall_time_s = float(time.perf_counter() - t0_wall)

    st_ret = sys_nb.observe_physical_state()
    n_flips_bwd = len(sys_nb.atlas_manager.transition_history) - n_flips_fwd
    k_bwd = [int(s.dyadic_level) for s in sys_nb.history.snapshots[1 + n_steps : 1 + 2 * n_steps]]
    max_k_all = max(k_fwd + k_bwd) if (k_fwd or k_bwd) else 0

    q_rev = np.asarray(st_ret.q, dtype=np.float64)
    p_rev = np.asarray(st_ret.p, dtype=np.float64)
    v_rev = p_rev / masses[:, None]
    obs_rev = compute_common_observables_f64(masses, q_rev, v_rev, t=0.0, G=G)

    rev_dq_norm = float(np.linalg.norm(q_rev - q0))
    rev_dp_norm = float(np.linalg.norm(p_rev - p0))
    rev_rel_dE = float(abs(obs_rev["H"] - H0) / abs(H0))
    palindrome_ok = (k_fwd == list(reversed(k_bwd)))

    return {
        "solver": solver_key,
        "category": category,
        "precision": "Float64",
        "t_turnaround": float(t_turnaround),
        "rev_dq_norm": rev_dq_norm,
        "rev_dp_norm": rev_dp_norm,
        "rev_rel_dE": rev_rel_dE,
        "max_dyadic_level": int(max_k_all),
        "saturation_events": list(sys_nb.scheduler.saturation_events),
        "n_flips_fwd": int(n_flips_fwd),
        "n_flips_bwd": int(n_flips_bwd),
        "dyadic_palindrome_matched": bool(palindrome_ok),
        "wall_time_s": wall_time_s,
    }


def run_normal_cone_reversibility_b4b_f128(
    solver_key: str,
    t_turnaround: float,
    dt: float = 0.04,
    max_dyadic_level: int = 16,
    G: float = 1.0,
    atlas_gamma: float = 1.0,
) -> Dict[str, Any]:
    with mp.workdps(35):
        st0 = build_b4_initial_state_f128()
        inv0 = compute_invariants_f128(st0, G=to_f128(G))
        H0 = inv0.H
        n_steps = int(round(t_turnaround / dt))
        dt_f128 = to_f128(str(dt))

        integ, enable_dyad, in_sh, category = _build_normal_cone_f128(solver_key, G=G)
        sys_nb = NBodySystemF128(
            state=st0.with_updates(),
            integrator=integ,
            G=to_f128(G),
            dt_macro=dt_f128,
            max_dyadic_level=int(max_dyadic_level),
            in_shadow_mode=in_sh,
            enable_dyadic_substepping=enable_dyad,
            atlas_gamma=float(atlas_gamma),
        )
        if enable_dyad:
            sys_nb.scheduler.warn_on_saturation = True

        t0_wall = time.perf_counter()
        for _ in range(n_steps):
            sys_nb.step(dt=dt_f128)
        n_flips_fwd = len(sys_nb.atlas_manager.transition_history)
        k_fwd = [int(s.dyadic_level) for s in sys_nb.history.snapshots[1 : 1 + n_steps]]

        for _ in range(n_steps):
            sys_nb.step(dt=-dt_f128)
        wall_time_s = float(time.perf_counter() - t0_wall)

        st_ret = sys_nb.observe_physical_state()
        n_flips_bwd = len(sys_nb.atlas_manager.transition_history) - n_flips_fwd
        k_bwd = [int(s.dyadic_level) for s in sys_nb.history.snapshots[1 + n_steps : 1 + 2 * n_steps]]
        max_k_all = max(k_fwd + k_bwd) if (k_fwd or k_bwd) else 0

        dq_sq = to_f128(0)
        dp_sq = to_f128(0)
        for i in range(3):
            for d in range(3):
                diff_q = st_ret.q[i, d] - st0.q[i, d]
                diff_p = st_ret.p[i, d] - st0.p[i, d]
                dq_sq += diff_q * diff_q
                dp_sq += diff_p * diff_p
        rev_dq_norm = float(mp.sqrt(dq_sq))
        rev_dp_norm = float(mp.sqrt(dp_sq))

        inv_ret = compute_invariants_f128(st_ret, G=to_f128(G))
        rev_rel_dE = float(abs(inv_ret.H - H0) / abs(H0))
        res_ret = compute_truncation_residual_f128(inv_ret, inv0)
        palindrome_ok = (k_fwd == list(reversed(k_bwd)))

        return {
            "solver": solver_key,
            "category": category,
            "precision": "Float128",
            "t_turnaround": float(t_turnaround),
            "rev_dq_norm": rev_dq_norm,
            "rev_dp_norm": rev_dp_norm,
            "rev_rel_dE": rev_rel_dE,
            "rev_dP_norm": float(res_ret[0]),
            "rev_dL_norm": float(res_ret[1]),
            "rev_dG_norm": float(res_ret[2]),
            "max_dyadic_level": int(max_k_all),
            "saturation_events": list(sys_nb.scheduler.saturation_events),
            "n_flips_fwd": int(n_flips_fwd),
            "n_flips_bwd": int(n_flips_bwd),
            "dyadic_palindrome_matched": bool(palindrome_ok),
            "wall_time_s": wall_time_s,
        }


def _worker_task(task: Tuple[str, Dict[str, Any]]) -> Tuple[str, Any]:
    kind, kwargs = task
    if kind == "rebound_b4a":
        return ("b4a", run_rebound_forward_b4a(**kwargs))
    elif kind == "nc_f64_b4a":
        return ("b4a", run_normal_cone_forward_b4a_f64(**kwargs))
    elif kind == "nc_f128_b4a":
        return ("b4a", run_normal_cone_forward_b4a_f128(**kwargs))
    elif kind == "rebound_b4b":
        return ("b4b", run_rebound_reversibility_b4b(**kwargs))
    elif kind == "nc_f64_b4b":
        return ("b4b", run_normal_cone_reversibility_b4b_f64(**kwargs))
    elif kind == "nc_f128_b4b":
        return ("b4b", run_normal_cone_reversibility_b4b_f128(**kwargs))
    else:
        raise ValueError(f"Unknown task kind: {kind}")


def main(protocol: str = "all", max_k: int = 8) -> None:
    print(f"=== Cycle VI - Iteration 75: Benchmark B4 [protocol={protocol}, max_k={max_k}] ===", flush=True)
    dt_base = 0.04

    rebound_solvers = [
        "rebound_ias15",
        "rebound_trace_bs",
        "rebound_trace_ias15",
        "rebound_mercurius",
        "rebound_whfast_c17",
        "rebound_whfast_c0",
    ]
    nc_f64_solvers = [
        "normal_cone_bm6_f64_dyadic",
        "normal_cone_pbm4_f64_dyadic",
        "normal_cone_bm6_f64_k0",
        "normal_cone_pbm4_f64_k0",
    ]
    nc_f128_solvers = [
        "normal_cone_bm6_f128_dyadic",
        "normal_cone_pbm4_f128_dyadic",
    ]

    horizons_f128 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0]
    horizons_f64 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0]

    tasks: List[Tuple[str, Dict[str, Any]]] = []

    if protocol in ("all", "nc_all", "nc_f128", "b4a_f128"):
        for sk in nc_f128_solvers:
            tasks.append(("nc_f128_b4a", {"solver_key": sk, "t_end": 70.0, "dt": dt_base, "max_dyadic_level": max_k}))
    if protocol in ("all", "nc_all", "nc_f64", "b4a_f64"):
        for sk in nc_f64_solvers:
            tasks.append(("nc_f64_b4a", {"solver_key": sk, "t_end": 70.0, "dt": dt_base, "max_dyadic_level": max_k}))
    if protocol in ("all", "rebound", "b4a_f64") and rebound is not None:
        for sk in rebound_solvers:
            tasks.append(("rebound_b4a", {"solver_key": sk, "t_end": 70.0, "dt": dt_base}))

    if protocol in ("all", "nc_all", "nc_f128", "b4b_f128"):
        for T in horizons_f128:
            for sk in nc_f128_solvers:
                tasks.append(("nc_f128_b4b", {"solver_key": sk, "t_turnaround": T, "dt": dt_base, "max_dyadic_level": max_k}))
    if protocol in ("all", "nc_all", "nc_f64", "b4b_f64"):
        for T in horizons_f64:
            for sk in nc_f64_solvers:
                tasks.append(("nc_f64_b4b", {"solver_key": sk, "t_turnaround": T, "dt": dt_base, "max_dyadic_level": max_k}))
    if protocol in ("all", "rebound", "b4b_f64") and rebound is not None:
        for T in horizons_f64:
            for sk in rebound_solvers:
                tasks.append(("rebound_b4b", {"solver_key": sk, "t_turnaround": T, "dt": dt_base}))

    n_workers = min(os.cpu_count() or 4, 16)
    print(f"Dispatching {len(tasks)} benchmark tasks across {n_workers} parallel workers...", flush=True)
    t0_all = time.perf_counter()

    out_dir = os.path.join(_REPO_ROOT, "results", "campaigns")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "rebound_b4_pythagorean_chaos.json")

    b4a_results: Dict[str, Any] = {}
    b4b_by_T: Dict[float, Dict[str, Any]] = {T: {} for T in horizons_f64}
    prev_wall = 0.0
    if os.path.exists(out_path):
        try:
            with open(out_path, "r", encoding="utf-8") as f:
                existing = json.load(f)
            b4a_results.update(existing.get("protocol_B4a_forward_chaos", {}))
            for row in existing.get("protocol_B4b_reversibility_sweep", []):
                T_ex = float(row["t_turnaround"])
                b4b_by_T.setdefault(T_ex, {}).update(row.get("solvers", {}))
            prev_wall = float(existing.get("total_wall_time_s", 0.0))
        except Exception:
            pass

    def _save_checkpoint() -> float:
        rows_b4b = [
            {"t_turnaround": float(T), "solvers": b4b_by_T[T]}
            for T in sorted(b4b_by_T.keys())
        ]
        el = prev_wall + float(time.perf_counter() - t0_all)
        output_data = {
            "benchmark": "B4_pythagorean_3_4_5_chaos_and_reversibility",
            "references": [
                "Burrau (1913), Astron. Nachr. 195, 113",
                "Szebehely & Peters (1967), AJ 72, 876",
                "Boekholt & Portegies Zwart (2015), Comput. Astrophys. Cosmol. 2:2 (Brutus)",
                "Rein & Spiegel (2015), MNRAS 446, 1424 (REBOUND IAS15)",
                "Rein et al. (2019), MNRAS 485, 5490 (REBOUND MERCURIUS)",
                "Lu, Hernandez & Rein (2024), MNRAS 533, 3708 (REBOUND TRACE)",
            ],
            "parameters": {
                "masses": [3.0, 4.0, 5.0],
                "q0_raw": [[1.0, 3.0, 0.0], [-2.0, -1.0, 0.0], [1.0, -1.0, 0.0]],
                "v0_raw": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                "E0_exact": -769.0 / 60.0,
                "dt_base": dt_base,
                "max_dyadic_level": max_k,
            },
            "protocol_B4a_forward_chaos": b4a_results,
            "protocol_B4b_reversibility_sweep": rows_b4b,
            "total_wall_time_s": el,
        }
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(output_data, f, indent=2)
        return el

    with mp_proc.Pool(processes=n_workers) as pool:
        for kind_tag, res in pool.imap_unordered(_worker_task, tasks):
            if kind_tag == "b4a":
                sk = res["solver"]
                b4a_results[sk] = res
                print(
                    f"  [B4-a DONE] {sk:28s} | T={res['t_end']:5.1f} | max_dE={res['max_rel_dE']:.3e} | "
                    f"end_dE={res['end_rel_dE']:.3e} | dL={res['dL_norm']:.3e} | wall={res['wall_time_pure_s']:.2f}s",
                    flush=True,
                )
            elif kind_tag == "b4b":
                sk = res["solver"]
                T = float(res["t_turnaround"])
                b4b_by_T.setdefault(T, {})[sk] = res
                print(
                    f"  [B4-b DONE] T={T:5.1f} | {sk:28s} | ||dq_FB||_2={res['rev_dq_norm']:.3e} | "
                    f"dE_FB={res['rev_rel_dE']:.3e} | wall={res['wall_time_s']:.2f}s",
                    flush=True,
                )
            _save_checkpoint()

    el_final = _save_checkpoint()
    print(f"Saved certified B4 JSON to {out_path} in {el_final:.2f}s", flush=True)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=str, default="all")
    parser.add_argument("--max-k", type=int, default=16)
    args = parser.parse_args()
    main(protocol=args.protocol, max_k=args.max_k)
