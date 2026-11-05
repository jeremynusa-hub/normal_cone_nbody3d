#!/usr/bin/env python3
"""
Cycle VI - Iteration 73: Head-to-Head Benchmark B1
Perturbed Star-Grazing / Ultra-High-Eccentricity Pericenter Passage (1 - e1 in [1e-1, 1e-9] and e1 = 1.0)
References:
  - Lu, Hernandez & Rein (2024), MNRAS 533, 3708 (TRACE, Section 4.2)
  - Rauch & Holman (1999), AJ 117, 1087
  - Hernandez & Dehnen (2017), MNRAS 468, 2614

Design Principles:
  1. Strictly equal floating-point precision: IEEE-754 Float64 on both sides.
  2. Simple, modular, mirrored architecture between REBOUND 5.2.2 and Normal-Cone Symplectic Splitting Solver (PBM4 & BM6).
  3. Zero arbitrary floor/clamping on any physical observable or error metric.
  4. Comparative wall-clock timing recorded in every benchmark regime (both pure integration
     timing and per-step monitored timing).
"""

from __future__ import annotations
import json
import math
import os
import sys
import time
from typing import Dict, Any, List, Tuple

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import numpy as np
import rebound

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
from normal_cone_nbody_float64.integrators.processed import (
    ProcessedSymplecticIntegrator,
    CorrectorC4,
)


# ============================================================================
# 1. Shared Initial Conditions & Common Physical Observables (Float64)
# ============================================================================

def build_b1_initial_conditions(
    one_minus_e: float,
    a_grazer: float = 1.0,
    a_perturber: float = 5.2,
    e_perturber: float = 0.05,
    inc_perturber: float = 0.02,
    f_perturber: float = 0.7,
    swap_planet_order: bool = False,
    G: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Constructs bit-identical barycentric initial conditions (masses, q0, v0, p0) in Float64.
    The star-grazing planet starts at apocenter (f = pi) so both solvers initialize in the
    smooth orbital regime before plunging through the pericenter singularity.
    When one_minus_e == 0.0 (e = 1.0), initializes the exact rectilinear collision orbit at
    apocenter distance r_apo = 2 * a_grazer with zero transverse velocity relative to the star.
    """
    sim = rebound.Simulation()
    sim.G = float(G)
    sim.add(m=1.0)  # Index 0: Central Star

    def _add_grazer(s: rebound.Simulation) -> None:
        if one_minus_e > 0.0:
            s.add(m=1e-3, a=float(a_grazer), e=float(1.0 - one_minus_e), f=math.pi)
        else:
            # Exact radial orbit e = 1.0 at apocenter r = 2 * a_grazer (v_rel = 0 at turnaround)
            p_star = s.particles[0]
            s.add(
                m=1e-3,
                x=float(p_star.x - 2.0 * a_grazer),
                y=float(p_star.y),
                z=float(p_star.z),
                vx=float(p_star.vx),
                vy=float(p_star.vy),
                vz=float(p_star.vz),
            )

    def _add_perturber(s: rebound.Simulation) -> None:
        s.add(
            m=1e-3,
            a=float(a_perturber),
            e=float(e_perturber),
            inc=float(inc_perturber),
            f=float(f_perturber),
        )

    if not swap_planet_order:
        _add_grazer(sim)
        _add_perturber(sim)
    else:
        _add_perturber(sim)
        _add_grazer(sim)

    sim.move_to_com()

    masses = np.array([p.m for p in sim.particles], dtype=np.float64)
    q0 = np.array([[p.x, p.y, p.z] for p in sim.particles], dtype=np.float64)
    v0 = np.array([[p.vx, p.vy, p.vz] for p in sim.particles], dtype=np.float64)
    p0 = masses[:, None] * v0
    return masses, q0, v0, p0


def compute_common_observables(
    masses: np.ndarray,
    q: np.ndarray,
    v: np.ndarray,
    t: float,
    G: float = 1.0,
) -> Dict[str, Any]:
    """
    Identical Float64 physical observable evaluator shared by REBOUND and Conormal Solver:
      - Total Newtonian energy H = T + V (unsoftened, epsilon = 0)
      - Total linear momentum P = sum_i m_i v_i
      - Total angular momentum L = sum_i q_i x (m_i v_i)
      - Galilean boost generator G_boost = M_tot R_cm - t P
    """
    p = masses[:, None] * v
    T = float(0.5 * np.sum(masses * np.sum(v * v, axis=1)))
    V = 0.0
    n = len(masses)
    for i in range(n):
        for j in range(i + 1, n):
            rij = float(np.linalg.norm(q[i] - q[j]))
            V -= float(G * masses[i] * masses[j]) / rij
    H = T + V
    P = np.sum(p, axis=0)
    L = np.sum(np.cross(q, p), axis=0)
    M_tot = float(np.sum(masses))
    R_cm = np.sum(masses[:, None] * q, axis=0) / M_tot
    G_boost = M_tot * R_cm - float(t) * P
    return {
        "H": float(H),
        "P": P,
        "L": L,
        "G_boost": G_boost,
    }


# ============================================================================
# 2. Mirrored Solver Runners: REBOUND 5.2.2 vs Normal-Cone Solver (Float64)
# ============================================================================

SOLVER_METADATA: Dict[str, Dict[str, Any]] = {
    # --- Category 1: Adaptive / Multi-Scale Solvers ---
    "normal_cone_bm6_dyadic": {
        "category": "adaptive_multiscale",
        "order_label": "6th-ord. sympl. (Chart Flip + Palindromic Dyadic)",
    },
    "normal_cone_pbm4_dyadic": {
        "category": "adaptive_multiscale",
        "order_label": "4th/6th-ord. Gauge (Chart Flip + Palindromic Dyadic)",
    },
    "rebound_ias15": {
        "category": "adaptive_multiscale",
        "order_label": "15th-ord. Gauss-Radau (Adaptive PRS23)",
    },
    "rebound_trace_bs": {
        "category": "adaptive_multiscale",
        "order_label": "2nd-ord. DHC + 16th-ord. BS (Adaptive Switching)",
    },
    "rebound_trace_ias15": {
        "category": "adaptive_multiscale",
        "order_label": "2nd-ord. DHC + 15th-ord. IAS15 (Adaptive Switching)",
    },
    "rebound_mercurius": {
        "category": "adaptive_multiscale",
        "order_label": "2nd-ord. DHC + 15th-ord. IAS15 (Adaptive Hybrid)",
    },
    # --- Category 2: Fixed-Step Solvers ---
    "normal_cone_bm6_k0": {
        "category": "fixed_step",
        "order_label": "6th-ord. sympl. (Chart Flip, Fixed Step k=0)",
    },
    "normal_cone_pbm4_k0": {
        "category": "fixed_step",
        "order_label": "4th/6th-ord. Gauge (Chart Flip, Fixed Step k=0)",
    },
    "rebound_whfast_c17": {
        "category": "fixed_step",
        "order_label": "2nd-ord. + 17th-ord. Corrector (Static Jacobi, Fixed Step)",
    },
    "rebound_whfast_c0": {
        "category": "fixed_step",
        "order_label": "2nd-ord. Wisdom-Holman (Static Jacobi, Fixed Step)",
    },
}


def _configure_rebound_integrator(sim: rebound.Simulation, solver_key: str) -> None:
    """Applies verified REBOUND 5.2.2 integrator parameters after setting sim.integrator."""
    if solver_key == "rebound_whfast_c0":
        sim.integrator = "whfast"
        sim.integrator.corrector = 0
        sim.integrator.safe_mode = 0
        sim.integrator.keep_unsynchronized = 1
    elif solver_key == "rebound_whfast_c17":
        sim.integrator = "whfast"
        sim.integrator.corrector = 17
        sim.integrator.safe_mode = 0
        sim.integrator.keep_unsynchronized = 1
    elif solver_key == "rebound_mercurius":
        sim.integrator = "mercurius"
        sim.integrator.r_crit_hill = 3.0
        sim.integrator.safe_mode = 0
    elif solver_key == "rebound_trace_bs":
        sim.integrator = "trace"
        sim.integrator.r_crit_hill = 3.0
        sim.integrator.peri_crit_eta = 1.0
        sim.integrator.peri_mode = "FULL_BS"
    elif solver_key == "rebound_trace_ias15":
        sim.integrator = "trace"
        sim.integrator.r_crit_hill = 3.0
        sim.integrator.peri_crit_eta = 1.0
        sim.integrator.peri_mode = "FULL_IAS15"
    elif solver_key == "rebound_ias15":
        sim.integrator = "ias15"
        sim.integrator.epsilon = 1e-9
        sim.integrator.adaptive_mode = "PRS23"
    else:
        raise ValueError(f"Unknown REBOUND solver_key: {solver_key}")


def run_rebound_solver(
    masses: np.ndarray,
    q0: np.ndarray,
    v0: np.ndarray,
    dt: float,
    n_steps: int,
    solver_key: str,
    n_orbits: int,
    test_reversibility: bool = True,
    record_series: bool = False,
    G: float = 1.0,
) -> Dict[str, Any]:
    """
    Mirrored runner for REBOUND 5.2.2:
      1. Measures pure uninterrupted forward integration wall-clock time (wall_time_pure_s).
      2. Measures step-by-step monitored trajectory (max_rel_dE, end_rel_dE, Noether invariants,
         and round-trip time-reversibility error ||dq_FB||_2).
    """
    # --- Part A: Pure integration timing & exact uninterrupted round-trip reversibility ---
    sim_pure = rebound.Simulation()
    sim_pure.G = float(G)
    for i in range(len(masses)):
        sim_pure.add(
            m=float(masses[i]),
            x=float(q0[i, 0]), y=float(q0[i, 1]), z=float(q0[i, 2]),
            vx=float(v0[i, 0]), vy=float(v0[i, 1]), vz=float(v0[i, 2]),
        )
    _configure_rebound_integrator(sim_pure, solver_key)
    if solver_key in ("rebound_whfast_c0", "rebound_whfast_c17"):
        sim_pure.integrator.keep_unsynchronized = 0
    sim_pure.dt = float(dt)

    t_start_pure = time.perf_counter()
    if solver_key == "rebound_ias15":
        sim_pure.integrate(float(n_steps * dt), exact_finish_time=1)
    else:
        sim_pure.steps(int(n_steps))
        sim_pure.synchronize()
    wall_time_pure_s = float(time.perf_counter() - t_start_pure)

    rev_dq_norm = float("nan")
    rev_dp_norm = float("nan")
    if test_reversibility:
        q_turn = np.array([[p.x, p.y, p.z] for p in sim_pure.particles], dtype=np.float64)
        v_turn = np.array([[p.vx, p.vy, p.vz] for p in sim_pure.particles], dtype=np.float64)
        sim_bwd = rebound.Simulation()
        sim_bwd.G = float(G)
        for i in range(len(masses)):
            sim_bwd.add(
                m=float(masses[i]),
                x=float(q_turn[i, 0]), y=float(q_turn[i, 1]), z=float(q_turn[i, 2]),
                vx=float(-v_turn[i, 0]), vy=float(-v_turn[i, 1]), vz=float(-v_turn[i, 2]),
            )
        _configure_rebound_integrator(sim_bwd, solver_key)
        if solver_key in ("rebound_whfast_c0", "rebound_whfast_c17"):
            sim_bwd.integrator.keep_unsynchronized = 0
        sim_bwd.dt = float(dt)
        if solver_key == "rebound_ias15":
            sim_bwd.integrate(float(n_steps * dt), exact_finish_time=1)
        else:
            sim_bwd.steps(int(n_steps))
            sim_bwd.synchronize()
        q_rev = np.array([[p.x, p.y, p.z] for p in sim_bwd.particles], dtype=np.float64)
        v_rev = -np.array([[p.vx, p.vy, p.vz] for p in sim_bwd.particles], dtype=np.float64)
        p_rev = masses[:, None] * v_rev
        p0 = masses[:, None] * v0
        rev_dq_norm = float(np.linalg.norm(q_rev - q0))
        rev_dp_norm = float(np.linalg.norm(p_rev - p0))

    # --- Part B: Step-monitored integration ---
    sim = rebound.Simulation()
    sim.G = float(G)
    for i in range(len(masses)):
        sim.add(
            m=float(masses[i]),
            x=float(q0[i, 0]), y=float(q0[i, 1]), z=float(q0[i, 2]),
            vx=float(v0[i, 0]), vy=float(v0[i, 1]), vz=float(v0[i, 2]),
        )
    obs0 = compute_common_observables(masses, q0, v0, t=0.0, G=G)
    H0 = obs0["H"]
    _configure_rebound_integrator(sim, solver_key)
    sim.dt = float(dt)

    max_rel_dE = 0.0
    series_t: List[float] = []
    series_rel_dE: List[float] = []

    t_start_mon = time.perf_counter()
    for k in range(1, n_steps + 1):
        if solver_key == "rebound_ias15":
            sim.integrate(float(k * dt), exact_finish_time=1)
        else:
            sim.steps(1)
            sim.synchronize()
        rel_dE = abs(float(sim.energy()) - H0) / abs(H0)
        if rel_dE > max_rel_dE:
            max_rel_dE = rel_dE
        if record_series:
            series_t.append(float(sim.t))
            series_rel_dE.append(float(rel_dE))
    wall_time_monitored_s = float(time.perf_counter() - t_start_mon)

    q_cur = np.array([[p.x, p.y, p.z] for p in sim.particles], dtype=np.float64)
    v_cur = np.array([[p.vx, p.vy, p.vz] for p in sim.particles], dtype=np.float64)
    obs = compute_common_observables(masses, q_cur, v_cur, t=float(sim.t), G=G)
    end_rel_dE = float(abs(obs["H"] - H0) / abs(H0))
    dP_norm = float(np.linalg.norm(obs["P"] - obs0["P"]))
    dL_norm = float(np.linalg.norm(obs["L"] - obs0["L"]))
    dG_norm = float(np.linalg.norm(obs["G_boost"] - obs0["G_boost"]))

    res: Dict[str, Any] = {
        "solver": solver_key,
        "category": SOLVER_METADATA[solver_key]["category"],
        "order_label": SOLVER_METADATA[solver_key]["order_label"],
        "max_rel_dE": float(max_rel_dE),
        "end_rel_dE": float(end_rel_dE),
        "dP_norm": dP_norm,
        "dL_norm": dL_norm,
        "dG_norm": dG_norm,
        "rev_dq_norm": rev_dq_norm,
        "rev_dp_norm": rev_dp_norm,
        "wall_time_pure_s": wall_time_pure_s,
        "wall_time_monitored_s": wall_time_monitored_s,
        "ms_per_orbit_pure": float(wall_time_pure_s * 1000.0 / max(1, n_orbits)),
        "us_per_step_pure": float(wall_time_pure_s * 1e6 / max(1, n_steps)),
    }
    if record_series:
        res["series_t"] = series_t
        res["series_rel_dE"] = series_rel_dE
    return res


def _build_normal_cone_integrator(solver_key: str, G: float = 1.0) -> Tuple[Any, bool]:
    """
    Constructs the Normal-Cone integrator (PBM4 or BM6) and returns (integrator, enable_dyadic_substepping).
    All variants ALWAYS operate with dynamical Jacobi chart flips (AtlasManager + shadow gauge connection).
    """
    if solver_key == "normal_cone_pbm4_dyadic":
        return (
            ProcessedSymplecticIntegrator(
                BlanesMoan4Integrator(G=G),
                CorrectorC4(G=G, mode="homological"),
            ),
            True,
        )
    elif solver_key == "normal_cone_bm6_dyadic":
        return BlanesMoan6Integrator(G=G), True
    elif solver_key == "normal_cone_pbm4_k0":
        return (
            ProcessedSymplecticIntegrator(
                BlanesMoan4Integrator(G=G),
                CorrectorC4(G=G, mode="homological"),
            ),
            False,
        )
    elif solver_key == "normal_cone_bm6_k0":
        return BlanesMoan6Integrator(G=G), False
    else:
        raise ValueError(f"Unknown Normal-Cone solver_key: {solver_key}")


def run_normal_cone_solver(
    masses: np.ndarray,
    q0: np.ndarray,
    v0: np.ndarray,
    dt: float,
    n_steps: int,
    solver_key: str,
    n_orbits: int,
    test_reversibility: bool = True,
    record_series: bool = False,
    G: float = 1.0,
) -> Dict[str, Any]:
    """
    Mirrored runner for our Unsoftened Normal-Cone Symplectic Splitting Solver (PBM4 & BM6 in Float64):
      - Always runs with dynamical Jacobi chart flips (AtlasManager) + shadow gauge connection.
      - Runs either in Adaptive/Multi-Scale mode (enable_dyadic_substepping=True, max_dyadic_level=10)
        or in Fixed-Step mode (enable_dyadic_substepping=False, k=0).
      - Measures integration wall-clock time, step-by-step energy/Noether invariants,
        active Jacobi leaf pair, dyadic levels, and exact round-trip reversibility ||dq_FB||_2.
    """
    p0 = masses[:, None] * v0
    obs0 = compute_common_observables(masses, q0, v0, t=0.0, G=G)
    H0 = obs0["H"]

    integ, enable_dyad = _build_normal_cone_integrator(solver_key, G=G)
    in_sh = "pbm4" in solver_key
    st0 = ExtendedPhaseState(q=q0.copy(), p=p0.copy(), masses=masses.copy(), t=0.0)
    sys_nb = NBodySystem(
        state=st0,
        integrator=integ,
        dt_macro=float(dt),
        max_dyadic_level=10,
        in_shadow_mode=in_sh,
        enable_dyadic_substepping=enable_dyad,
        use_optimized_kernels=True,
    )
    if enable_dyad:
        sys_nb.scheduler.warn_on_saturation = True

    t_start_pure = time.perf_counter()
    for _ in range(n_steps):
        sys_nb.step(float(dt))
    wall_time_pure_s = float(time.perf_counter() - t_start_pure)
    wall_time_monitored_s = wall_time_pure_s

    snaps = sys_nb.history.snapshots[1:]
    series_t = [float(s.t) for s in snaps]
    series_rel_dE = [float(abs(s.invariants.H - H0) / abs(H0)) for s in snaps]
    max_rel_dE = float(max(series_rel_dE)) if series_rel_dE else 0.0
    end_rel_dE = float(series_rel_dE[-1]) if series_rel_dE else 0.0

    st_end = sys_nb.observe_physical_state()
    q_end = np.asarray(st_end.q, dtype=np.float64)
    v_end = np.asarray(st_end.p, dtype=np.float64) / masses[:, None]
    obs_end = compute_common_observables(masses, q_end, v_end, t=float(st_end.t), G=G)

    dP_norm = float(np.linalg.norm(obs_end["P"] - obs0["P"]))
    dL_norm = float(np.linalg.norm(obs_end["L"] - obs0["L"]))
    dG_norm = float(np.linalg.norm(obs_end["G_boost"] - obs0["G_boost"]))
    dyad_arr = np.asarray(sys_nb.history.dyadic_levels, dtype=int)
    max_k = int(np.max(dyad_arr)) if dyad_arr.size > 0 else 0
    total_substeps = int(np.sum(1 << dyad_arr[1:])) if dyad_arr.size > 1 else n_steps
    active_leaf = [list(pair) for pair in sys_nb.tree.get_leaf_pair_indices()]
    sat_events_fwd = list(sys_nb.scheduler.saturation_events)
    if enable_dyad:
        assert len(sat_events_fwd) == 0, f"Dyadic saturation on {solver_key}: {sat_events_fwd}"

    rev_dq_norm = float("nan")
    rev_dp_norm = float("nan")
    if test_reversibility:
        if in_sh and sys_nb._shadow_state is not None:
            sys_nb._shadow_state = ExtendedPhaseState(
                q=sys_nb._shadow_state.q.copy(),
                p=-sys_nb._shadow_state.p.copy(),
                masses=masses.copy(),
                t=sys_nb._shadow_state.t,
            )
            sys_nb.state = ExtendedPhaseState(
                q=sys_nb._shadow_state.q.copy(),
                p=sys_nb._shadow_state.p.copy(),
                masses=masses.copy(),
                t=sys_nb._shadow_state.t,
            )
        else:
            sys_nb.state = ExtendedPhaseState(
                q=sys_nb.state.q.copy(),
                p=-sys_nb.state.p.copy(),
                masses=masses.copy(),
                t=sys_nb.state.t,
            )
        for _ in range(n_steps):
            sys_nb.step(float(dt))
        if in_sh and sys_nb._shadow_state is not None:
            sys_nb._shadow_state = ExtendedPhaseState(
                q=sys_nb._shadow_state.q.copy(),
                p=-sys_nb._shadow_state.p.copy(),
                masses=masses.copy(),
                t=sys_nb._shadow_state.t,
            )
            sys_nb.state = ExtendedPhaseState(
                q=sys_nb._shadow_state.q.copy(),
                p=sys_nb._shadow_state.p.copy(),
                masses=masses.copy(),
                t=sys_nb._shadow_state.t,
            )
        else:
            sys_nb.state = ExtendedPhaseState(
                q=sys_nb.state.q.copy(),
                p=-sys_nb.state.p.copy(),
                masses=masses.copy(),
                t=sys_nb.state.t,
            )
        st_rev = sys_nb.observe_physical_state()
        q_rev = np.asarray(st_rev.q, dtype=np.float64)
        p_rev = np.asarray(st_rev.p, dtype=np.float64)
        rev_dq_norm = float(np.linalg.norm(q_rev - q0))
        rev_dp_norm = float(np.linalg.norm(p_rev - p0))
        if enable_dyad:
            assert len(sys_nb.scheduler.saturation_events) == 0, (
                f"Dyadic saturation on reverse pass of {solver_key}: {sys_nb.scheduler.saturation_events}"
            )

    res: Dict[str, Any] = {
        "solver": solver_key,
        "category": SOLVER_METADATA[solver_key]["category"],
        "order_label": SOLVER_METADATA[solver_key]["order_label"],
        "max_rel_dE": float(max_rel_dE),
        "end_rel_dE": float(end_rel_dE),
        "dP_norm": dP_norm,
        "dL_norm": dL_norm,
        "dG_norm": dG_norm,
        "rev_dq_norm": rev_dq_norm,
        "rev_dp_norm": rev_dp_norm,
        "max_dyadic_level": max_k,
        "total_substeps": total_substeps,
        "saturation_events": list(sys_nb.scheduler.saturation_events),
        "active_leaf": active_leaf,
        "wall_time_pure_s": wall_time_pure_s,
        "wall_time_monitored_s": wall_time_monitored_s,
        "ms_per_orbit_pure": float(wall_time_pure_s * 1000.0 / max(1, n_orbits)),
        "us_per_step_pure": float(wall_time_pure_s * 1e6 / max(1, n_steps)),
    }
    if record_series:
        res["series_t"] = series_t
        res["series_rel_dE"] = series_rel_dE
    return res


def warmup_normal_cone_jit() -> None:
    """Executes a 2-step warmup on Normal-Cone solvers so Numba JIT compilation is excluded from timers."""
    masses, q0, v0, p0 = build_b1_initial_conditions(1e-2)
    dt = 0.02
    for key in ["normal_cone_pbm4_dyadic", "normal_cone_bm6_dyadic", "normal_cone_pbm4_k0", "normal_cone_bm6_k0"]:
        integ, enable_dyad = _build_normal_cone_integrator(key, G=1.0)
        st = ExtendedPhaseState(q=q0.copy(), p=p0.copy(), masses=masses.copy(), t=0.0)
        sys_nb = NBodySystem(
            state=st,
            integrator=integ,
            dt_macro=dt,
            max_dyadic_level=10,
            in_shadow_mode=("pbm4" in key),
            enable_dyadic_substepping=enable_dyad,
            use_optimized_kernels=True,
        )
        sys_nb.step(dt)
        sys_nb.step(-dt)


def main() -> None:
    print("Warming up Normal-Cone Numba JIT kernels...")
    warmup_normal_cone_jit()

    ecc_grid = [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-9, 0.0]
    T_orb = 2.0 * math.pi / math.sqrt(1.0 + 1e-3)
    steps_per_orbit = 200
    dt = T_orb / float(steps_per_orbit)
    n_orbits_sweep = 5
    n_steps_sweep = n_orbits_sweep * steps_per_orbit

    normal_cone_adaptive_keys = [
        "normal_cone_bm6_dyadic",
        "normal_cone_pbm4_dyadic",
    ]
    normal_cone_fixed_keys = [
        "normal_cone_bm6_k0",
        "normal_cone_pbm4_k0",
    ]
    normal_cone_keys = normal_cone_adaptive_keys + normal_cone_fixed_keys

    rebound_adaptive_keys = [
        "rebound_ias15",
        "rebound_trace_bs",
        "rebound_trace_ias15",
        "rebound_mercurius",
    ]
    rebound_fixed_keys = [
        "rebound_whfast_c17",
        "rebound_whfast_c0",
    ]
    rebound_keys = rebound_adaptive_keys + rebound_fixed_keys

    # 1. Sweep B1-a: Canonical Inside-Out Ordering [Star, Grazer(a=1.0), Perturber(a=5.2)]
    print(f"\n=== Running B1-a Sweep ({n_orbits_sweep} orbits, dt = T_1/{steps_per_orbit}) ===")
    b1a_results: List[Dict[str, Any]] = []
    for one_minus_e in ecc_grid:
        masses, q0, v0, _ = build_b1_initial_conditions(one_minus_e, swap_planet_order=False)
        entry: Dict[str, Any] = {
            "one_minus_e": float(one_minus_e),
            "e1": float(1.0 - one_minus_e),
            "r_peri_au": float(one_minus_e),
            "solvers": {},
        }
        print(f"  [B1-a] 1 - e1 = {one_minus_e:.0e}:")
        for skey in normal_cone_keys:
            r = run_normal_cone_solver(masses, q0, v0, dt, n_steps_sweep, skey, n_orbits_sweep, test_reversibility=True)
            entry["solvers"][skey] = r
            print(
                f"    [{r['category']:19s}] {skey:22s}: max|dE/E0|={r['max_rel_dE']:.3e} | rev|dq|={r['rev_dq_norm']:.2e} "
                f"| k_max={r['max_dyadic_level']} | leaf={r['active_leaf']} | time={r['wall_time_pure_s']*1e3:7.2f} ms ({r['ms_per_orbit_pure']:.2f} ms/orb)"
            )
        for skey in rebound_keys:
            r = run_rebound_solver(masses, q0, v0, dt, n_steps_sweep, skey, n_orbits_sweep, test_reversibility=True)
            entry["solvers"][skey] = r
            print(
                f"    [{r['category']:19s}] {skey:22s}: max|dE/E0|={r['max_rel_dE']:.3e} | rev|dq|={r['rev_dq_norm']:.2e} "
                f"| time={r['wall_time_pure_s']*1e3:7.2f} ms ({r['ms_per_orbit_pure']:.2f} ms/orb)"
            )
        b1a_results.append(entry)

    # 2. Sweep B1-b: Permuted Ordering [Star, Perturber(a=5.2), Grazer(a=1.0)]
    print(f"\n=== Running B1-b Topological Ordering Stress Test ({n_orbits_sweep} orbits) ===")
    b1b_results: List[Dict[str, Any]] = []
    for one_minus_e in ecc_grid:
        masses, q0, v0, _ = build_b1_initial_conditions(one_minus_e, swap_planet_order=True)
        entry = {
            "one_minus_e": float(one_minus_e),
            "e2": float(1.0 - one_minus_e),
            "r_peri_au": float(one_minus_e),
            "solvers": {},
        }
        print(f"  [B1-b] 1 - e2 = {one_minus_e:.0e}:")
        for skey in normal_cone_keys:
            r = run_normal_cone_solver(masses, q0, v0, dt, n_steps_sweep, skey, n_orbits_sweep, test_reversibility=True)
            entry["solvers"][skey] = r
            print(
                f"    [{r['category']:19s}] {skey:22s}: max|dE/E0|={r['max_rel_dE']:.3e} | rev|dq|={r['rev_dq_norm']:.2e} "
                f"| k_max={r['max_dyadic_level']} | leaf={r['active_leaf']} | time={r['wall_time_pure_s']*1e3:7.2f} ms"
            )
        for skey in rebound_keys:
            r = run_rebound_solver(masses, q0, v0, dt, n_steps_sweep, skey, n_orbits_sweep, test_reversibility=True)
            entry["solvers"][skey] = r
            print(
                f"    [{r['category']:19s}] {skey:22s}: max|dE/E0|={r['max_rel_dE']:.3e} | rev|dq|={r['rev_dq_norm']:.2e} "
                f"| time={r['wall_time_pure_s']*1e3:7.2f} ms"
            )
        b1b_results.append(entry)

    # 3. 10-Orbit Time Series & Round-Trip Reversibility at Extreme Star-Grazing 1 - e1 = 1e-6
    n_orbits_ts = 10
    n_steps_ts = n_orbits_ts * steps_per_orbit
    one_minus_e_ts = 1e-6
    print(f"\n=== Running {n_orbits_ts}-Orbit Time Series at 1 - e1 = {one_minus_e_ts:.0e} ===")
    masses_ts, q0_ts, v0_ts, _ = build_b1_initial_conditions(one_minus_e_ts, swap_planet_order=False)
    ts_solvers: Dict[str, Any] = {}
    for skey in normal_cone_keys:
        ts_solvers[skey] = run_normal_cone_solver(
            masses_ts, q0_ts, v0_ts, dt, n_steps_ts, skey, n_orbits_ts, test_reversibility=True, record_series=True
        )
    for skey in rebound_keys:
        ts_solvers[skey] = run_rebound_solver(
            masses_ts, q0_ts, v0_ts, dt, n_steps_ts, skey, n_orbits_ts, test_reversibility=True, record_series=True
        )

    payload = {
        "benchmark": "B1_Perturbed_Star_Grazing_Pericenter_Lu2024_RauchHolman1999",
        "rebound_version": rebound.__version__,
        "precision": "Float64 (IEEE-754 double precision on both sides)",
        "solver_categories": {
            "adaptive_multiscale": normal_cone_adaptive_keys + rebound_adaptive_keys,
            "fixed_step": normal_cone_fixed_keys + rebound_fixed_keys,
        },
        "solver_metadata": SOLVER_METADATA,
        "parameters": {
            "G": 1.0,
            "m_star": 1.0,
            "m_grazer": 1e-3,
            "a_grazer": 1.0,
            "m_perturber": 1e-3,
            "a_perturber": 5.2,
            "e_perturber": 0.05,
            "inc_perturber": 0.02,
            "f_perturber": 0.7,
            "T_orb_grazer": T_orb,
            "steps_per_orbit": steps_per_orbit,
            "dt_macro": dt,
            "n_orbits_sweep": n_orbits_sweep,
            "n_orbits_timeseries": n_orbits_ts,
            "one_minus_e_timeseries": one_minus_e_ts,
        },
        "b1a_canonical_sweep": b1a_results,
        "b1b_permuted_sweep": b1b_results,
        "timeseries_10orbits_1e_6": ts_solvers,
    }

    os.makedirs("results/campaigns", exist_ok=True)
    out_path = "results/campaigns/rebound_b1_star_grazing_pericenter.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"\nSaved certified JSON to {out_path}")


if __name__ == "__main__":
    main()
