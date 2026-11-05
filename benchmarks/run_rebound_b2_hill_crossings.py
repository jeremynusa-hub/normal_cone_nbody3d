#!/usr/bin/env python3
"""
Cycle VI - Iteration 74: Head-to-Head Benchmark B2
Repeated Hill-Sphere Crossings, Shadow Gauge Connection & Hybrid Time-Reversibility
References:
  - Rein et al. (2019), MNRAS 485, 5490 (MERCURIUS)
  - Wisdom (2017), MNRAS 464, 2350 (Symplectic correctors & Hill switching boundary)
  - Hernandez & Dehnen (2023), MNRAS 522, 4639 (Hybrid time-reversibility breakdown)
  - Lu, Hernandez & Rein (2024), MNRAS 533, 3708 (TRACE)

Design Principles:
  1. Strictly equal floating-point precision: IEEE-754 Float64 on both sides.
  2. Simple, modular, mirrored architecture between REBOUND 5.2.2 and
     Normal-Cone Symplectic Splitting Solver (Normal-Cone PBM4 & Normal-Cone BM6).
  3. Zero arbitrary floor/clamping on any physical observable or error metric.
  4. Comparative wall-clock timing recorded in every benchmark regime.
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
from normal_cone_nbody_float64.core.invariants import compute_invariants
from normal_cone_nbody_float64.system.nbody_system import NBodySystem
from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator
from normal_cone_nbody_float64.integrators.blanes_moan_6 import BlanesMoan6Integrator
from normal_cone_nbody_float64.integrators.processed import (
    ProcessedSymplecticIntegrator,
    CorrectorC4,
    apply_gauge_connection,
)


# ============================================================================
# 1. Shared Initial Conditions & Common Physical Observables (Float64)
# ============================================================================

def build_b2a_initial_conditions(
    rp_target: float,
    m_intruder: float = 1e-6,
    r0_jup: float = 0.70,
    G: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Protocol B2-a: Deep Planetocentric Hill-Sphere Crossing & Jacobi Chart Flip
    (T_helio = [(0, 2)] -> T_Jovi = [(0, 2), (1, 3)] -> T_helio = [(0, 2)]).
      - Body 0: Sun (m0 = 1.0)
      - Body 1: Jupiter (m1 = 9.5458e-4, a1 = 5.20336 AU, R_Hill,J = 0.355121 AU)
      - Body 2: Earth (m2 = 3.0035e-6, a2 = 1.00000 AU)
      - Body 3: Planetary/satellite intruder (m3 = 1e-6) starting outside Jupiter's
        Hill switching sphere at r0_J = 0.70 AU (~1.97 R_Hill,J) and plunging through
        perijove rp_target in [0.10, 1e-4] AU before exiting back beyond 0.70 AU.
    Returns bit-identical barycentric (masses, q0, v0, p0, r_hill_jup) in Float64.
    """
    m_sun = 1.0
    m_jup = 9.5458e-4
    m_earth = 3.0035e-6
    masses = np.array([m_sun, m_jup, m_earth, float(m_intruder)], dtype=np.float64)

    a_jup = 5.20336
    a_earth = 1.00000
    r_hill_jup = a_jup * ((m_jup / (3.0 * (m_sun + m_jup))) ** (1.0 / 3.0))

    v_jup = math.sqrt(float(G) * (m_sun + m_jup) / a_jup)
    v_earth = math.sqrt(float(G) * (m_sun + m_earth) / a_earth)
    v_inf_j = 0.25 * v_jup
    mu_j = float(G) * (m_jup + float(m_intruder))
    v0_j = math.sqrt(v_inf_j * v_inf_j + 2.0 * mu_j / float(r0_jup))
    vp_j = math.sqrt(v_inf_j * v_inf_j + 2.0 * mu_j / float(rp_target))

    # Coriolis/solar-tidal baseline offset from r0_J = 0.70 AU plus impact parameter
    sin_base = -0.0258152
    sin_imp = (float(rp_target) * vp_j) / (float(r0_jup) * v0_j)
    sin_gamma = max(-0.99, min(0.99, sin_base + sin_imp))
    cos_gamma = math.sqrt(1.0 - sin_gamma * sin_gamma)

    q = np.array([
        [0.0, 0.0, 0.0],
        [a_jup, 0.0, 0.0],
        [0.0, a_earth, 0.0],
        [a_jup + float(r0_jup), 0.0, 1.0e-3],
    ], dtype=np.float64)
    v = np.array([
        [0.0, 0.0, 0.0],
        [0.0, v_jup, 0.0],
        [-v_earth, 0.0, 0.0],
        [-v0_j * cos_gamma, v_jup + v0_j * sin_gamma, 0.0],
    ], dtype=np.float64)

    M_tot = float(np.sum(masses))
    q -= np.sum(masses[:, None] * q, axis=0) / M_tot
    v -= np.sum(masses[:, None] * v, axis=0) / M_tot
    p = masses[:, None] * v
    return masses, q, v, p, float(r_hill_jup)


def build_b2b_initial_conditions(
    G: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Protocol B2-b: Repeated Multi-Encounter Two-Giant-Planet Hill-Sphere Crossings
    (Rein et al. 2019 MERCURIUS §4 & Hernandez & Dehnen 2023):
      - Body 0: Central Star (m0 = 1.0)
      - Body 1: Inner Giant Planet (m1 = 1e-3, a1 = 1.00 AU, e1 = 0.05, inc1 = 0.01 rad, f1 = 0)
      - Body 2: Outer Giant Planet (m2 = 1e-3, a2 = 1.22 AU, e2 = 0.05, inc2 = 0.01 rad, pomega2 = 0.3, f2 = pi)
    Returns bit-identical barycentric (masses, q0, v0, p0, r_hill_mutual) in Float64.
    """
    sim = rebound.Simulation()
    sim.G = float(G)
    sim.add(m=1.0)
    sim.add(m=1e-3, a=1.00, e=0.05, inc=0.01, f=0.0)
    sim.add(m=1e-3, a=1.22, e=0.05, inc=0.01, pomega=0.3, f=math.pi)
    sim.move_to_com()

    masses = np.array([p.m for p in sim.particles], dtype=np.float64)
    q0 = np.array([[p.x, p.y, p.z] for p in sim.particles], dtype=np.float64)
    v0 = np.array([[p.vx, p.vy, p.vz] for p in sim.particles], dtype=np.float64)
    p0 = masses[:, None] * v0
    r_hill_mut = 0.5 * (1.00 + 1.22) * (((1e-3 + 1e-3) / 3.0) ** (1.0 / 3.0))
    return masses, q0, v0, p0, float(r_hill_mut)


def compute_common_observables(
    masses: np.ndarray,
    q: np.ndarray,
    v: np.ndarray,
    t: float,
    G: float = 1.0,
) -> Dict[str, Any]:
    """
    Identical Float64 physical observable evaluator shared by REBOUND and Normal-Cone Solver:
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

def _configure_rebound_integrator(
    sim: rebound.Simulation,
    solver_key: str,
    step_by_step_monitoring: bool = False,
) -> None:
    """
    Applies verified REBOUND 5.2.2 integrator parameters after setting sim.integrator.
    Note: For MERCURIUS, safe_mode=1 is mandatory during step-by-step monitoring because
    reb_integrator_mercurius_synchronize does not set recalculate_coordinates_this_timestep=1
    when safe_mode=0.
    """
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
        sim.integrator.safe_mode = 1 if step_by_step_monitoring else 0
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
    n_orbits: float,
    pair_indices: Tuple[int, int] = (1, 2),
    r_hill_ref: float = 1.0,
    test_reversibility: bool = True,
    record_series: bool = False,
    G: float = 1.0,
) -> Dict[str, Any]:
    """
    Mirrored runner for REBOUND 5.2.2:
      1. Measures pure integration wall-clock time (wall_time_pure_s) without per-step Python overhead.
      2. Measures step-by-step monitored trajectory (max_rel_dE, end_rel_dE, min pair separation,
         Hill crossing count, Noether invariants, and round-trip time-reversibility ||dq_FB||_2).
    """
    # --- Pass 1: Pure integration timing & exact forward-backward reversibility ---
    sim_pure = rebound.Simulation()
    sim_pure.G = float(G)
    for i in range(len(masses)):
        sim_pure.add(
            m=float(masses[i]),
            x=float(q0[i, 0]), y=float(q0[i, 1]), z=float(q0[i, 2]),
            vx=float(v0[i, 0]), vy=float(v0[i, 1]), vz=float(v0[i, 2]),
        )
    _configure_rebound_integrator(sim_pure, solver_key, step_by_step_monitoring=False)
    sim_pure.dt = float(dt)
    t_end = float(n_steps * dt)

    t_start_pure = time.perf_counter()
    if solver_key == "rebound_ias15":
        sim_pure.integrate(t_end, exact_finish_time=1)
    else:
        sim_pure.steps(n_steps)
        sim_pure.synchronize()
    wall_time_pure_s = float(time.perf_counter() - t_start_pure)

    rev_dq_norm = float("nan")
    rev_dp_norm = float("nan")
    if test_reversibility:
        if solver_key in ("rebound_whfast_c0", "rebound_whfast_c17"):
            sim_pure.integrator.keep_unsynchronized = 0
            sim_pure.synchronize()
        q_turn = np.array([[p.x, p.y, p.z] for p in sim_pure.particles], dtype=np.float64)
        v_turn = np.array([[p.vx, p.vy, p.vz] for p in sim_pure.particles], dtype=np.float64)

        sim_rev = rebound.Simulation()
        sim_rev.G = float(G)
        for i in range(len(masses)):
            sim_rev.add(
                m=float(masses[i]),
                x=float(q_turn[i, 0]), y=float(q_turn[i, 1]), z=float(q_turn[i, 2]),
                vx=float(-v_turn[i, 0]), vy=float(-v_turn[i, 1]), vz=float(-v_turn[i, 2]),
            )
        _configure_rebound_integrator(sim_rev, solver_key, step_by_step_monitoring=False)
        sim_rev.dt = float(dt)
        if solver_key == "rebound_ias15":
            sim_rev.integrate(t_end, exact_finish_time=1)
        else:
            sim_rev.steps(n_steps)
            sim_rev.synchronize()
        q_rev = np.array([[p.x, p.y, p.z] for p in sim_rev.particles], dtype=np.float64)
        v_rev = np.array([[-p.vx, -p.vy, -p.vz] for p in sim_rev.particles], dtype=np.float64)
        p_rev = masses[:, None] * v_rev
        p0 = masses[:, None] * v0
        rev_dq_norm = float(np.linalg.norm(q_rev - q0))
        rev_dp_norm = float(np.linalg.norm(p_rev - p0))

    # --- Pass 2: Step-by-step monitored trajectory ---
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
    _configure_rebound_integrator(sim, solver_key, step_by_step_monitoring=True)
    sim.dt = float(dt)

    ia, ib = pair_indices
    max_rel_dE = 0.0
    min_pair_dist = float(np.linalg.norm(q0[ia] - q0[ib]))
    in_enc = (min_pair_dist < 3.0 * r_hill_ref)
    n_hill_crossings = 1 if in_enc else 0

    series_t: List[float] = []
    series_rel_dE: List[float] = []
    series_r_over_rhill: List[float] = []

    t_start_mon = time.perf_counter()
    for k in range(1, n_steps + 1):
        if solver_key == "rebound_ias15":
            sim.integrate(float(k * dt), exact_finish_time=1)
        else:
            sim.steps(1)
            if solver_key in ("rebound_whfast_c0", "rebound_whfast_c17"):
                sim.synchronize()
        rel_dE = abs(float(sim.energy()) - H0) / abs(H0)
        if rel_dE > max_rel_dE:
            max_rel_dE = rel_dE

        pa, pb = sim.particles[ia], sim.particles[ib]
        r_ab = math.sqrt((pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2)
        if r_ab < min_pair_dist:
            min_pair_dist = r_ab
        if r_ab < 3.0 * r_hill_ref and not in_enc:
            in_enc = True
            n_hill_crossings += 1
        elif r_ab >= 3.0 * r_hill_ref and in_enc:
            in_enc = False

        if record_series:
            series_t.append(float(sim.t))
            series_rel_dE.append(float(rel_dE))
            series_r_over_rhill.append(float(r_ab / r_hill_ref))
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
        "max_rel_dE": float(max_rel_dE),
        "end_rel_dE": float(end_rel_dE),
        "min_pair_dist_au": float(min_pair_dist),
        "min_pair_over_rhill": float(min_pair_dist / r_hill_ref),
        "n_hill_crossings": int(n_hill_crossings),
        "dP_norm": dP_norm,
        "dL_norm": dL_norm,
        "dG_norm": dG_norm,
        "rev_dq_norm": rev_dq_norm,
        "rev_dp_norm": rev_dp_norm,
        "wall_time_pure_s": wall_time_pure_s,
        "wall_time_monitored_s": wall_time_monitored_s,
        "ms_per_orbit_pure": float(wall_time_pure_s * 1000.0 / max(1e-9, float(n_orbits))),
        "us_per_step_pure": float(wall_time_pure_s * 1e6 / max(1, n_steps)),
    }
    if record_series:
        res["series_t"] = series_t
        res["series_rel_dE"] = series_rel_dE
        res["series_r_over_rhill"] = series_r_over_rhill
    return res


def _build_normal_cone_integrator(solver_key: str, G: float = 1.0) -> Tuple[Any, bool, bool]:
    """Returns (integrator_instance, enable_dyadic_substepping, in_shadow_mode) for Normal-Cone PBM4 and BM6."""
    if solver_key == "normal_cone_pbm4_k0":
        return (
            ProcessedSymplecticIntegrator(
                BlanesMoan4Integrator(G=G),
                CorrectorC4(G=G, mode="homological"),
            ),
            False,
            True,
        )
    elif solver_key == "normal_cone_bm6_k0":
        return BlanesMoan6Integrator(G=G), False, False
    elif solver_key == "normal_cone_pbm4_dyadic":
        return (
            ProcessedSymplecticIntegrator(
                BlanesMoan4Integrator(G=G),
                CorrectorC4(G=G, mode="homological"),
            ),
            True,
            True,
        )
    elif solver_key == "normal_cone_bm6_dyadic":
        return BlanesMoan6Integrator(G=G), True, False
    else:
        raise ValueError(f"Unknown Normal-Cone solver_key: {solver_key}")


def run_normal_cone_solver(
    masses: np.ndarray,
    q0: np.ndarray,
    v0: np.ndarray,
    dt: float,
    n_steps: int,
    solver_key: str,
    n_orbits: float,
    pair_indices: Tuple[int, int] = (1, 2),
    r_hill_ref: float = 1.0,
    test_reversibility: bool = True,
    record_series: bool = False,
    record_substeps: bool = False,
    G: float = 1.0,
) -> Dict[str, Any]:
    """
    Mirrored runner for our Unsoftened Normal-Cone Symplectic Splitting Solver (PBM4 & BM6 in Float64):
      1. Measures pure integration wall-clock time (wall_time_pure_s) with optimized kernels.
      2. Measures trajectory metrics (max_rel_dE, end_rel_dE, Noether invariants, chart flips,
         max dyadic level, and round-trip time-reversibility ||dq_FB||_2 via canonical velocity reversal).
    """
    p0 = masses[:, None] * v0
    obs0 = compute_common_observables(masses, q0, v0, t=0.0, G=G)
    H0 = obs0["H"]
    ia, ib = pair_indices

    integ, enable_dyad, in_sh = _build_normal_cone_integrator(solver_key, G=G)
    st0 = ExtendedPhaseState(q=q0.copy(), p=p0.copy(), masses=masses.copy(), t=0.0)
    sys_nb = NBodySystem(
        state=st0,
        integrator=integ,
        dt_macro=float(dt),
        max_dyadic_level=16,
        in_shadow_mode=in_sh,
        enable_dyadic_substepping=enable_dyad,
        use_optimized_kernels=True,
        record_substeps=record_substeps,
        substep_record_min_level=4,
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
    series_r = [float(np.linalg.norm(s.q[ia] - s.q[ib])) for s in snaps]
    series_r_over_rhill = [float(r / r_hill_ref) for r in series_r]
    series_dyadic_k = [int(s.dyadic_level) for s in snaps]

    max_rel_dE = float(max(series_rel_dE)) if series_rel_dE else 0.0
    end_rel_dE = float(series_rel_dE[-1]) if series_rel_dE else 0.0
    min_pair_dist = float(min(series_r)) if series_r else float(np.linalg.norm(q0[ia] - q0[ib]))

    in_enc = False
    n_hill_crossings = 0
    for r_val in series_r:
        if r_val < 3.0 * r_hill_ref and not in_enc:
            in_enc = True
            n_hill_crossings += 1
        elif r_val >= 3.0 * r_hill_ref and in_enc:
            in_enc = False

    st_end = sys_nb.observe_physical_state()
    q_end = np.asarray(st_end.q, dtype=np.float64)
    v_end = np.asarray(st_end.p, dtype=np.float64) / masses[:, None]
    obs_end = compute_common_observables(masses, q_end, v_end, t=float(st_end.t), G=G)

    dP_norm = float(np.linalg.norm(obs_end["P"] - obs0["P"]))
    dL_norm = float(np.linalg.norm(obs_end["L"] - obs0["L"]))
    dG_norm = float(np.linalg.norm(obs_end["G_boost"] - obs0["G_boost"]))
    max_k = int(max(series_dyadic_k)) if series_dyadic_k else 0
    n_flips = len(sys_nb.atlas_manager.transition_history)
    flips_info = [
        {
            "t": float(tr.t),
            "old_leaf_pairs": [list(p) for p in tr.old_leaf_pairs],
            "new_leaf_pairs": [list(p) for p in tr.new_leaf_pairs],
        }
        for tr in sys_nb.atlas_manager.transition_history
    ]
    sat_events_fwd = list(sys_nb.scheduler.saturation_events)
    if enable_dyad:
        assert len(sat_events_fwd) == 0, f"Dyadic saturation on {solver_key}: {sat_events_fwd}"

    rev_dq_norm = float("nan")
    rev_dp_norm = float("nan")
    if test_reversibility:
        sys_nb.record_substeps = False
        if len(masses) == 4:
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
        else:
            for _ in range(n_steps):
                sys_nb.step(-float(dt))
        st_back = sys_nb.observe_physical_state()
        q_rev = np.asarray(st_back.q, dtype=np.float64)
        p_rev = np.asarray(st_back.p, dtype=np.float64)
        rev_dq_norm = float(np.linalg.norm(q_rev - q0))
        rev_dp_norm = float(np.linalg.norm(p_rev - p0))
        if enable_dyad:
            assert len(sys_nb.scheduler.saturation_events) == 0, (
                f"Dyadic saturation on reverse pass of {solver_key}: {sys_nb.scheduler.saturation_events}"
            )

    res: Dict[str, Any] = {
        "solver": solver_key,
        "max_rel_dE": float(max_rel_dE),
        "end_rel_dE": float(end_rel_dE),
        "min_pair_dist_au": float(min_pair_dist),
        "min_pair_over_rhill": float(min_pair_dist / r_hill_ref),
        "n_hill_crossings": int(n_hill_crossings),
        "n_atlas_flips": int(n_flips),
        "atlas_flips": flips_info,
        "dP_norm": dP_norm,
        "dL_norm": dL_norm,
        "dG_norm": dG_norm,
        "rev_dq_norm": rev_dq_norm,
        "rev_dp_norm": rev_dp_norm,
        "max_dyadic_level": max_k,
        "saturation_events": list(sys_nb.scheduler.saturation_events),
        "wall_time_pure_s": wall_time_pure_s,
        "wall_time_monitored_s": wall_time_monitored_s,
        "ms_per_orbit_pure": float(wall_time_pure_s * 1000.0 / max(1e-9, float(n_orbits))),
        "us_per_step_pure": float(wall_time_pure_s * 1e6 / max(1, n_steps)),
    }
    if record_series:
        res["series_t"] = series_t
        res["series_rel_dE"] = series_rel_dE
        res["series_r_over_rhill"] = series_r_over_rhill
        res["series_dyadic_k"] = series_dyadic_k
    return res


def _warmup_jit() -> None:
    """Warms up Numba JIT kernels before benchmark timing."""
    masses, q0, v0, _, r_hill = build_b2b_initial_conditions()
    dt = 2.0 * math.pi / 100.0
    for key in ["normal_cone_pbm4_dyadic", "normal_cone_bm6_dyadic"]:
        run_normal_cone_solver(
            masses, q0, v0, dt, n_steps=5, solver_key=key,
            n_orbits=0.05, pair_indices=(1, 2), r_hill_ref=r_hill, test_reversibility=False,
        )


def compute_gauge_connection_ablation_b2a(
    rp_targets: List[float],
    dt: float,
    n_steps: int,
    rows_b2a: List[Dict[str, Any]],
    G: float = 1.0,
) -> List[Dict[str, Any]]:
    """
    Directly measures the intruder orbital energy discontinuity at the Hill-sphere
    chart flip T_A = [(0, 2)] -> T_B = [(0, 2), (1, 3)] with our non-commutative
    covariant shadow gauge connection T_{A->B}^{(PBM4)} vs naive shadow transfer (no gauge).
    """
    def _intruder_energy(state: ExtendedPhaseState) -> float:
        m = state.masses
        q = state.q
        v = state.p / m[:, None]
        e3 = 0.5 * float(m[3]) * float(np.dot(v[3], v[3]))
        for j in range(3):
            r3j = float(np.linalg.norm(q[3] - q[j]))
            e3 -= float(G * m[3] * m[j]) / r3j
        return e3

    ablation_rows: List[Dict[str, Any]] = []
    for idx, rp_t in enumerate(rp_targets):
        masses, q0, v0, p0, r_hill_jup = build_b2a_initial_conditions(rp_t, m_intruder=1e-6, G=G)
        integ, enable_dyad, in_sh = _build_normal_cone_integrator("normal_cone_pbm4_dyadic", G=G)
        st0 = ExtendedPhaseState(q=q0.copy(), p=p0.copy(), masses=masses.copy(), t=0.0)
        sys_nb = NBodySystem(
            state=st0,
            integrator=integ,
            dt_macro=float(dt),
            max_dyadic_level=16,
            in_shadow_mode=in_sh,
            enable_dyadic_substepping=enable_dyad,
            use_optimized_kernels=True,
        )
        flip_step = 0
        flip_t = 0.0
        jump_with_gauge = 0.0
        jump_naive = 0.0
        for step_idx in range(1, n_steps + 1):
            n_prev = len(sys_nb.atlas_manager.transition_history)
            sh_before = sys_nb.get_shadow_state()
            tree_before = sys_nb.tree
            h_before = float(sys_nb._last_dt_substep)
            sys_nb.step(float(dt))
            if len(sys_nb.atlas_manager.transition_history) > n_prev:
                flip_step = step_idx
                flip_t = float(sys_nb.atlas_manager.transition_history[-1].t)
                tree_after = sys_nb.tree
                phys_before = integ.observe_physical(sh_before, tree_before, h_before, apex=sys_nb.apex)
                sh_gauge = apply_gauge_connection(
                    sh_before, tree_before, tree_after, dt=h_before, corrector=integ.corrector, apex=sys_nb.apex
                )
                phys_with_gauge = integ.observe_physical(sh_gauge, tree_after, h_before, apex=sys_nb.apex)
                phys_naive = integ.observe_physical(sh_before, tree_after, h_before, apex=sys_nb.apex)
                e3_ref = _intruder_energy(phys_before)
                e3_gauge = _intruder_energy(phys_with_gauge)
                e3_naive = _intruder_energy(phys_naive)
                jump_with_gauge = abs(e3_gauge - e3_ref) / abs(e3_ref)
                jump_naive = abs(e3_naive - e3_ref) / abs(e3_ref)
                break
        ablation_rows.append({
            "rp_target_au": float(rp_t),
            "r_min_over_rhill": float(rows_b2a[idx]["r_min_over_rhill"]),
            "flip_step": int(flip_step),
            "flip_t": float(flip_t),
            "jump_rel_dE3_with_gauge": float(jump_with_gauge),
            "jump_rel_dE3_naive_no_gauge": float(jump_naive),
        })
    return ablation_rows


# ============================================================================
# 3. Full Benchmark B2 Execution (B2-a + B2-b + Time Series)
# ============================================================================

SOLVER_REGISTRY: Dict[str, Dict[str, Any]] = {
    # Category 1: Adaptive / Multi-Scale Solvers
    "normal_cone_bm6_dyadic": {
        "category": "adaptive_multiscale",
        "label": "Normal-Cone BM6 (6th-ord., Chart Flip + Dyadic)",
        "order": "6th-order Symplectic + Chart Flip + Dyadic",
    },
    "normal_cone_pbm4_dyadic": {
        "category": "adaptive_multiscale",
        "label": "Normal-Cone PBM4 (4th/6th-ord. Gauge, Chart Flip + Dyadic)",
        "order": "4th-ord. Kernel + 6th-ord. C4 Gauge + Chart Flip + Dyadic",
    },
    "rebound_ias15": {
        "category": "adaptive_multiscale",
        "label": "REBOUND IAS15 (15th-ord. Gauss-Radau, Adaptive)",
        "order": "15th-order Non-Symplectic Adaptive",
    },
    "rebound_trace_bs": {
        "category": "adaptive_multiscale",
        "label": "REBOUND TRACE (2nd DHC + 16th BS, Adaptive)",
        "order": "2nd-ord. DHC + Adaptive Bulirsch-Stoer",
    },
    "rebound_trace_ias15": {
        "category": "adaptive_multiscale",
        "label": "REBOUND TRACE (2nd DHC + 15th IAS15, Adaptive)",
        "order": "2nd-ord. DHC + Adaptive IAS15",
    },
    "rebound_mercurius": {
        "category": "adaptive_multiscale",
        "label": "REBOUND MERCURIUS (2nd DHC + 15th IAS15, Adaptive)",
        "order": "2nd-ord. DHC + Smooth IAS15 Switch",
    },
    # Category 2: Fixed-Step Solvers (k=0 / Fixed dt)
    "normal_cone_bm6_k0": {
        "category": "fixed_step",
        "label": "Normal-Cone BM6 (6th-ord., Chart Flip, Fixed Step k=0)",
        "order": "6th-order Symplectic + Chart Flip (Fixed dt, k=0)",
    },
    "normal_cone_pbm4_k0": {
        "category": "fixed_step",
        "label": "Normal-Cone PBM4 (4th/6th-ord. Gauge, Chart Flip, Fixed Step k=0)",
        "order": "4th-ord. Kernel + 6th-ord. C4 Gauge + Chart Flip (Fixed dt, k=0)",
    },
    "rebound_whfast_c17": {
        "category": "fixed_step",
        "label": "REBOUND WHFast (2nd+17th corr., c=17, Fixed Step)",
        "order": "2nd-ord. WH + 17th-ord. Corrector (Fixed dt)",
    },
    "rebound_whfast_c0": {
        "category": "fixed_step",
        "label": "REBOUND WHFast (2nd-ord., c=0, Fixed Step)",
        "order": "2nd-order Wisdom-Holman (Fixed dt)",
    },
}

def main() -> None:
    print("Warming up Numba JIT kernels (Float64)...")
    _warmup_jit()

    rebound_solvers = [
        "rebound_whfast_c0",
        "rebound_whfast_c17",
        "rebound_mercurius",
        "rebound_trace_bs",
        "rebound_trace_ias15",
        "rebound_ias15",
    ]
    normal_cone_solvers = [
        "normal_cone_pbm4_k0",
        "normal_cone_bm6_k0",
        "normal_cone_pbm4_dyadic",
        "normal_cone_bm6_dyadic",
    ]

    # ------------------------------------------------------------------------
    # Protocol B2-a: Deep Planetocentric Hill-Sphere Crossing & Chart Flip
    # ------------------------------------------------------------------------
    rp_targets = [0.10, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001]
    dt_b2a = 0.08
    n_steps_b2a = 150
    n_orbits_b2a = float(n_steps_b2a * dt_b2a) / (2.0 * math.pi)  # in Earth-orbit units (~1.91 yr)

    results_b2a: List[Dict[str, Any]] = []
    print("\n=== Running Protocol B2-a (Deep Planetocentric Hill-Sphere Crossing & Chart Flip) ===")
    for rp_t in rp_targets:
        masses, q0, v0, p0, r_hill_jup = build_b2a_initial_conditions(rp_t, m_intruder=1e-6)
        row: Dict[str, Any] = {
            "rp_target_au": float(rp_t),
            "r_hill_jup_au": float(r_hill_jup),
            "dt": float(dt_b2a),
            "n_steps": int(n_steps_b2a),
            "solvers": {},
        }
        for s_key in rebound_solvers:
            res = run_rebound_solver(
                masses, q0, v0, dt_b2a, n_steps_b2a, s_key,
                n_orbits=n_orbits_b2a, pair_indices=(1, 3), r_hill_ref=r_hill_jup,
                test_reversibility=True, record_series=False,
            )
            row["solvers"][s_key] = res
        for s_key in normal_cone_solvers:
            res = run_normal_cone_solver(
                masses, q0, v0, dt_b2a, n_steps_b2a, s_key,
                n_orbits=n_orbits_b2a, pair_indices=(1, 3), r_hill_ref=r_hill_jup,
                test_reversibility=True, record_series=False,
                record_substeps=(s_key == "normal_cone_bm6_dyadic"),
            )
            row["solvers"][s_key] = res

        # True minimum Jovicentric distance from dyadic BM6
        r_min_actual = row["solvers"]["normal_cone_bm6_dyadic"]["min_pair_dist_au"]
        row["r_min_actual_au"] = float(r_min_actual)
        row["r_min_over_rhill"] = float(r_min_actual / r_hill_jup)
        results_b2a.append(row)

        print(
            f"  rp_target={rp_t:<7.4f} (r_min/R_H={row['r_min_over_rhill']:.2e}) | "
            f"WHFast(c17)={row['solvers']['rebound_whfast_c17']['max_rel_dE']:.2e} | "
            f"MERCURIUS={row['solvers']['rebound_mercurius']['max_rel_dE']:.2e} | "
            f"TRACE={row['solvers']['rebound_trace_bs']['max_rel_dE']:.2e} | "
            f"NC-PBM4(dyad)={row['solvers']['normal_cone_pbm4_dyadic']['max_rel_dE']:.2e} | "
            f"NC-BM6(dyad)={row['solvers']['normal_cone_bm6_dyadic']['max_rel_dE']:.2e} | "
            f"NC-BM6(k=0)={row['solvers']['normal_cone_bm6_k0']['max_rel_dE']:.2e}"
        )

    # ------------------------------------------------------------------------
    # Protocol B2-b: Repeated Multi-Encounter Two-Giant-Planet Hill Crossings
    # ------------------------------------------------------------------------
    orbit_horizons = [1, 2, 3, 5, 10, 13, 15, 20]
    steps_per_orbit = 100
    dt_b2b = (2.0 * math.pi) / float(steps_per_orbit)
    masses_b2b, q0_b2b, v0_b2b, p0_b2b, r_hill_mut = build_b2b_initial_conditions()

    results_b2b: List[Dict[str, Any]] = []
    time_series_20orb: Dict[str, Any] = {}
    print("\n=== Running Protocol B2-b (Repeated Multi-Encounter Two-Giant-Planet Hill Crossings) ===")
    for n_orb in orbit_horizons:
        n_steps_b2b = int(n_orb * steps_per_orbit)
        is_last = (n_orb == 20)
        row_b: Dict[str, Any] = {
            "n_orbits": int(n_orb),
            "n_steps": int(n_steps_b2b),
            "dt": float(dt_b2b),
            "r_hill_mutual_au": float(r_hill_mut),
            "solvers": {},
        }
        for s_key in rebound_solvers:
            res = run_rebound_solver(
                masses_b2b, q0_b2b, v0_b2b, dt_b2b, n_steps_b2b, s_key,
                n_orbits=float(n_orb), pair_indices=(1, 2), r_hill_ref=r_hill_mut,
                test_reversibility=True, record_series=is_last,
            )
            if is_last:
                time_series_20orb[s_key] = {
                    "series_t": res.pop("series_t"),
                    "series_rel_dE": res.pop("series_rel_dE"),
                    "series_r_over_rhill": res.pop("series_r_over_rhill"),
                }
            row_b["solvers"][s_key] = res

        for s_key in normal_cone_solvers:
            res = run_normal_cone_solver(
                masses_b2b, q0_b2b, v0_b2b, dt_b2b, n_steps_b2b, s_key,
                n_orbits=float(n_orb), pair_indices=(1, 2), r_hill_ref=r_hill_mut,
                test_reversibility=True, record_series=is_last,
                record_substeps=False,
            )
            if is_last:
                time_series_20orb[s_key] = {
                    "series_t": res.pop("series_t"),
                    "series_rel_dE": res.pop("series_rel_dE"),
                    "series_r_over_rhill": res.pop("series_r_over_rhill"),
                    "series_dyadic_k": res.pop("series_dyadic_k"),
                }
            row_b["solvers"][s_key] = res

        row_b["min_r12_over_rhill"] = float(row_b["solvers"]["rebound_ias15"]["min_pair_over_rhill"])
        row_b["n_hill_crossings"] = int(row_b["solvers"]["rebound_ias15"]["n_hill_crossings"])
        results_b2b.append(row_b)

        print(
            f"  N_orb={n_orb:2d} (cross={row_b['n_hill_crossings']}, r_min={row_b['min_r12_over_rhill']:.2f} R_H) | "
            f"dE: WHFast17={row_b['solvers']['rebound_whfast_c17']['max_rel_dE']:.2e}, "
            f"MERC={row_b['solvers']['rebound_mercurius']['max_rel_dE']:.2e}, "
            f"TRACE={row_b['solvers']['rebound_trace_bs']['max_rel_dE']:.2e}, "
            f"NC-BM6={row_b['solvers']['normal_cone_bm6_dyadic']['max_rel_dE']:.2e} | "
            f"rev ||dq||: MERC={row_b['solvers']['rebound_mercurius']['rev_dq_norm']:.2e}, "
            f"TRACE={row_b['solvers']['rebound_trace_bs']['rev_dq_norm']:.2e}, "
            f"NC-BM6={row_b['solvers']['normal_cone_bm6_dyadic']['rev_dq_norm']:.2e}"
        )

    gauge_ablation = compute_gauge_connection_ablation_b2a(
        rp_targets, dt_b2a, n_steps_b2a, results_b2a
    )

    payload = {
        "benchmark_id": "B2",
        "title": "Repeated Hill-Sphere Crossings, Shadow Gauge Connection & Hybrid Time-Reversibility",
        "precision": "IEEE-754 Float64",
        "rebound_version": rebound.__version__,
        "solver_registry": SOLVER_REGISTRY,
        "protocol_B2a_deep_hill_crossing": results_b2a,
        "protocol_B2b_repeated_hill_crossings": results_b2b,
        "time_series_B2b_20orbits": time_series_20orb,
        "gauge_connection_ablation_B2a": gauge_ablation,
    }

    out_paths = [
        "rebound_b2_hill_crossings.json",
        os.path.join("results", "campaigns", "rebound_b2_hill_crossings.json"),
    ]
    for p in out_paths:
        d = os.path.dirname(p)
        if d:
            os.makedirs(d, exist_ok=True)
        with open(p, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
    print("\nSaved B2 results to:", out_paths[0])


if __name__ == "__main__":
    main()
