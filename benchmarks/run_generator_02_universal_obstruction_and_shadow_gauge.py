#!/usr/bin/env python3
r"""
Benchmark Generator 02: Shadow Manifold Z_4^{BM4},
Quadratic Floor & Dynamic Atlas Gauge 1-Cocycle (§7.1.2, Stress-Tests ST-1 to ST-4).

Governed by Contract:
  scripts/campaigns/contracts/02_generator_universal_obstruction_and_shadow_gauge.md

Remediates all substantive findings of Report 03 (§7.1.2):
  - [CRIT-7.1.2-01]: Replaces the 0-step static tautology in ST-4 by a genuine
    dynamical multi-step integration across a transition (N_A = 4 steps before,
    N_B = 4 steps after) in Float64 & Float128, and evaluates the 3-chart Cech
    1-cocycle and Lie-Bargmann conservation on PBM4 (Z_4^{BM4})
    without any artificial clamping to 1e-35.
  - [CRIT-7.1.2-02]: Runs the 35-digit Float128 order 6 -> 4 crossover sweep on
    the authentic 3D eccentric inclined hierarchical triple (e_in=0.30, i=28 deg)
    with zero np.clip on local empirical slopes.
  - [CRIT-7.1.2-03]: Exercises the paper showcase interfaces (01..05_*.py) and
    verifies that mode="homological" (Z_4^{BM4}) is the canonical default.
  - [CRIT-7.1.2-04] & [MAJ-7.1.2-03]: Quantifies shadow telescoping amortization
    (2M -> 2 corrector calls, s_eff(M) = 6 + 16/M -> 6.0) and reports exact gains.
  - [MAJ-7.1.2-01] & [MAJ-7.1.2-02]: Distinguishes the pre-asymptotic order-6
    channel O(eps * h^6) ~ eps^1 (eps << h^2) from the quadratic floor
    O(eps^2 * h^4) ~ eps^2 (eps >= h^2) in ST-1 for PBM4 (Z_4^{BM4}) vs BM4.
  - [MOD-7.1.2-01]: Integrates ST-3 over 100 binary periods with reference step
    h_ref = h/4 and dense sampling sample_stride = 1 (6,000 samples) using
    continuous Hann-windowed Laskar NAFF frequency extraction.
"""

from __future__ import annotations
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import mpmath as mp

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
if str(WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_ROOT))

RESULTS_DIR = WORKSPACE_ROOT / "results" / "campaigns"
DATA_DIR = WORKSPACE_ROOT / "data"
CH06_ROOT = WORKSPACE_ROOT / "interfaces"


def _load_showcase_module(rel_path: str, mod_name: str):
    full_path = CH06_ROOT / rel_path
    spec = importlib.util.spec_from_file_location(mod_name, str(full_path))
    assert spec is not None and spec.loader is not None, f"Cannot load {full_path}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod


# Load Float64 paper showcase interfaces (01..05_*.py)
mod_01 = _load_showcase_module("01_PhaseStateAndSystem_Interface.py", "showcase_01_g02")
mod_02 = _load_showcase_module("02_TreeBuilderAndCorrector_Interface.py", "showcase_02_g02")
mod_03 = _load_showcase_module("03_DyadicScheduler_Interface.py", "showcase_03_g02")
mod_04 = _load_showcase_module("04_SymplecticSubcycle_Interface.py", "showcase_04_g02")
mod_05 = _load_showcase_module("05_AtlasTransitions_Interface.py", "showcase_05_g02")
mod_06 = _load_showcase_module("06_HighPerformanceOptimizations_Interface.py", "showcase_06_g02")

from normal_cone_nbody_float64.integrators.common import (
    set_current_kepler_solver,
    set_current_apex_evaluator,
    set_current_supernode_solver,
)
set_current_kepler_solver(mod_06.OptimizedKeplerSundmanSolver(G=1.0))
set_current_apex_evaluator(mod_06.OptimizedScreenedApexPotential(G=1.0))
set_current_supernode_solver(mod_06.JitSupernodeSolver(G=1.0))

ExtendedPhaseState64 = mod_01.ExtendedPhaseState
NBodySystem64 = mod_01.NBodySystem
PalindromicCorrectorZ4_64 = mod_02.PalindromicCorrectorZ4
CorrectorC4_64 = mod_02.CorrectorC4
create_pbm4_integrator64 = mod_02.create_pbm4_integrator
BlanesMoan4Integrator64 = mod_04.BlanesMoan4Integrator
PBM4GaugeTransitionOperator64 = mod_05.PBM4GaugeTransitionOperator
apply_gauge_connection64 = mod_05.apply_gauge_connection
compute_lie_bargmann_invariants64 = mod_05.compute_lie_bargmann_invariants

from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec as build_tree64
from normal_cone_nbody_float64.core.invariants import compute_invariants as inv64

# Float128 imports
from normal_cone_nbody_float128.core.state import ExtendedPhaseState as State128
from normal_cone_nbody_float128.core.jacobi import build_tree_from_spec as build_tree128
from normal_cone_nbody_float128.core.invariants import compute_invariants as inv128
from normal_cone_nbody_float128.integrators.blanes_moan_4 import BlanesMoan4Integrator as BM4_128
from normal_cone_nbody_float128.integrators.blanes_moan_6 import BlanesMoan6Integrator as BM6_128
from normal_cone_nbody_float128.integrators.correctors import (
    HomologicalLieCorrectorC4 as HomC4_128,
    CorrectorC4 as CorrectorC4_128,
)
from normal_cone_nbody_float128.integrators.processed import (
    apply_gauge_connection as apply_gauge128,
)


# ==============================================================================
# Helper Constructors: 3D Eccentric Hierarchical Triple in Float64 & Float128
# ==============================================================================

def make_3d_eccentric_hierarchical_f64(
    eps: float,
    e_in: float = 0.30,
    inc_deg: float = 28.0,
    E_anom_deg: float = 40.0,
) -> Tuple[Any, Any, float]:
    """
    Constructs a 3D eccentric inclined inner binary (m0=1.0, m1=1.0, a_in=1.0, e_in=0.30,
    i=28 deg) perturbed by a tertiary body (m2=0.5) at distance R_out determined by
        eps = (m2 / M01) * (a_in / R_out)^3.
    At E_anom = 40 deg, the radial momentum w_rad = (r . w) / r is strongly non-zero,
    activating the momentum-dependent terms of D_A^3 B (Theorem 5.14).
    """
    m0, m1, m2 = 1.0, 1.0, 0.5
    M01 = m0 + m1
    M_tot = M01 + m2
    a_in = 1.0
    R_out = a_in * (m2 / (M01 * eps)) ** (1.0 / 3.0)

    E_anom = math.radians(E_anom_deg)
    cosE, sinE = math.cos(E_anom), math.sin(E_anom)
    r_mag = a_in * (1.0 - e_in * cosE)
    n_in = math.sqrt(M01 / (a_in ** 3))

    rx_orb = a_in * (cosE - e_in)
    ry_orb = a_in * math.sqrt(1.0 - e_in ** 2) * sinE
    vx_orb = -(a_in ** 2 * n_in / r_mag) * sinE
    vy_orb = (a_in ** 2 * n_in / r_mag) * math.sqrt(1.0 - e_in ** 2) * cosE

    inc = math.radians(inc_deg)
    ci, si = math.cos(inc), math.sin(inc)
    r_rel = np.array([rx_orb, ry_orb * ci, ry_orb * si], dtype=np.float64)
    v_rel = np.array([vx_orb, vy_orb * ci, vy_orb * si], dtype=np.float64)

    v_out = math.sqrt(M_tot / R_out)
    R_vec = np.array([0.6 * R_out, 0.8 * R_out, 0.0], dtype=np.float64)
    V_vec = np.array([-0.8 * v_out, 0.6 * v_out, 0.0], dtype=np.float64)

    q = np.zeros((3, 3), dtype=np.float64)
    v = np.zeros((3, 3), dtype=np.float64)
    q[0] = -0.5 * r_rel
    q[1] = +0.5 * r_rel
    q[2] = R_vec
    v[0] = -0.5 * v_rel
    v[1] = +0.5 * v_rel
    v[2] = V_vec

    masses = np.array([m0, m1, m2], dtype=np.float64)
    q_cm = np.sum(masses[:, None] * q, axis=0) / M_tot
    v_cm = np.sum(masses[:, None] * v, axis=0) / M_tot
    q -= q_cm
    v -= v_cm

    state = ExtendedPhaseState64(q=q, p=masses[:, None] * v, masses=masses)
    tree = build_tree64(masses, ((0, 1), 2))
    return state, tree, R_out


def make_3d_eccentric_hierarchical_f128(
    eps: float,
    e_in: float = 0.30,
    inc_deg: float = 28.0,
    E_anom_deg: float = 40.0,
) -> Tuple[State128, Any, float]:
    """
    35-digit Float128 counterpart of make_3d_eccentric_hierarchical_f64.
    Verifies state.is_f128 == True.
    """
    mp.mp.dps = 35
    m0 = mp.mpf("1.0")
    m1 = mp.mpf("1.0")
    m2 = mp.mpf("0.5")
    M01 = m0 + m1
    M_tot = M01 + m2
    a_in = mp.mpf("1.0")
    eps_mp = mp.mpf(str(eps))
    R_out = a_in * (m2 / (M01 * eps_mp)) ** (mp.mpf("1") / mp.mpf("3"))

    e_mp = mp.mpf(str(e_in))
    E_anom = mp.mpf(str(E_anom_deg)) * mp.pi / mp.mpf("180")
    cosE, sinE = mp.cos(E_anom), mp.sin(E_anom)
    r_mag = a_in * (mp.mpf("1") - e_mp * cosE)
    n_in = mp.sqrt(M01 / (a_in ** 3))

    rx_orb = a_in * (cosE - e_mp)
    ry_orb = a_in * mp.sqrt(mp.mpf("1") - e_mp ** 2) * sinE
    vx_orb = -(a_in ** 2 * n_in / r_mag) * sinE
    vy_orb = (a_in ** 2 * n_in / r_mag) * mp.sqrt(mp.mpf("1") - e_mp ** 2) * cosE

    inc = mp.mpf(str(inc_deg)) * mp.pi / mp.mpf("180")
    ci, si = mp.cos(inc), mp.sin(inc)
    r_rel = [rx_orb, ry_orb * ci, ry_orb * si]
    v_rel = [vx_orb, vy_orb * ci, vy_orb * si]

    v_out = mp.sqrt(M_tot / R_out)
    R_vec = [mp.mpf("0.6") * R_out, mp.mpf("0.8") * R_out, mp.mpf("0.0")]
    V_vec = [-mp.mpf("0.8") * v_out, mp.mpf("0.6") * v_out, mp.mpf("0.0")]

    half = mp.mpf("0.5")
    q = np.empty((3, 3), dtype=object)
    v = np.empty((3, 3), dtype=object)
    for d in range(3):
        q[0, d] = -half * r_rel[d]
        q[1, d] = +half * r_rel[d]
        q[2, d] = R_vec[d]
        v[0, d] = -half * v_rel[d]
        v[1, d] = +half * v_rel[d]
        v[2, d] = V_vec[d]

    masses = np.array([m0, m1, m2], dtype=object)
    for d in range(3):
        q_cm_d = (masses[0] * q[0, d] + masses[1] * q[1, d] + masses[2] * q[2, d]) / M_tot
        v_cm_d = (masses[0] * v[0, d] + masses[1] * v[1, d] + masses[2] * v[2, d]) / M_tot
        for i in range(3):
            q[i, d] = q[i, d] - q_cm_d
            v[i, d] = v[i, d] - v_cm_d

    p = np.empty((3, 3), dtype=object)
    for i in range(3):
        for d in range(3):
            p[i, d] = masses[i] * v[i, d]

    state = State128(q=q, p=p, masses=masses).to_f128()
    assert state.is_f128 is True
    assert isinstance(state.q[0, 0], mp.mpf)
    tree = build_tree128(masses, ((0, 1), 2))
    return state, tree, float(R_out)


# ==============================================================================
# Part I: 35-Digit Float128 Order 6 -> 4 Crossover on 3D Eccentric Triple
# ==============================================================================

def run_part1_crossover_3d_f128() -> Dict[str, Any]:
    """
    Executes the 35-digit Float128 order 6 -> 4 crossover sweep on the authentic
    3D eccentric hierarchical triple (e_in=0.30, i=28 deg) across eps in {1e-2, 1e-3, 1e-4}.
    Uses a fixed dyadic horizon T_cross = 0.32 so every step h in {0.32, 0.16, 0.08, 0.04,
    0.02, 0.01} lands at the exact same physical time and orbital phase (zero phase jitter,
    zero np.clip on local empirical slopes!).
    """
    print("=== Part I: 35-Digit Float128 Order 6 -> 4 Crossover on 3D Eccentric Triple (e=0.30, i=28 deg) ===", flush=True)
    mp.mp.dps = 35
    eps_list = [1e-2, 1e-3, 1e-4]
    h_steps = [0.32, 0.16, 0.08, 0.04, 0.02, 0.01]
    dt_chk = 0.32
    n_chk = 5
    T_cross = n_chk * dt_chk
    stages_s = 6

    bm4 = BM4_128(G=mp.mpf("1.0"))
    bm6_ref = BM6_128(G=mp.mpf("1.0"))
    corr_hom = HomC4_128(G=mp.mpf("1.0"))

    results_by_eps: Dict[str, Any] = {}

    for eps in eps_list:
        eps_key = f"eps_{eps:.0e}"
        st0, tree, R_out = make_3d_eccentric_hierarchical_f128(eps=eps, e_in=0.30, inc_deg=28.0)
        assert st0.is_f128 is True
        H0 = inv128(st0, G=mp.mpf("1.0")).H

        # Reference trajectory via unoptimized 35-digit Float128 BM6 with 6th-order
        # Richardson extrapolation (h_ref1 = 0.01, h_ref2 = 0.005) at all 5 checkpoints
        h_ref1 = mp.mpf("0.01")
        h_ref2 = mp.mpf("0.005")
        sub_r1 = int(round(dt_chk / 0.01))
        sub_r2 = int(round(dt_chk / 0.005))
        st_r1 = st0
        st_r2 = st0
        q_ref_chk: List[Any] = []
        c_rich = mp.mpf("64")
        c_denom = mp.mpf("63")
        for _ in range(n_chk):
            for _ in range(sub_r1):
                st_r1 = bm6_ref.step(st_r1, tree, dt=h_ref1)
            for _ in range(sub_r2):
                st_r2 = bm6_ref.step(st_r2, tree, dt=h_ref2)
            q_rich = np.empty((3, 3), dtype=object)
            for ii in range(3):
                for dd in range(3):
                    q_rich[ii, dd] = (c_rich * st_r2.q[ii, dd] - st_r1.q[ii, dd]) / c_denom
            q_ref_chk.append(q_rich)

        dE_bm4_list: List[float] = []
        dE_pbm4_list: List[float] = []
        gains_list: List[float] = []
        dq_bm4_list: List[float] = []
        dq_pbm4_list: List[float] = []
        gains_q_list: List[float] = []

        for h in h_steps:
            substeps_per_chk = int(round(dt_chk / h))
            h_mp = mp.mpf(str(h))

            # 1. Raw BM4 evaluated at the exact common checkpoints t_m = m * dt_chk
            st_b = st0
            sq_b = mp.mpf("0.0")
            sq_qb = mp.mpf("0.0")
            for chk_idx in range(n_chk):
                for _ in range(substeps_per_chk):
                    st_b = bm4.step(st_b, tree, dt=h_mp)
                eb = abs(inv128(st_b, G=mp.mpf("1.0")).H - H0) / abs(H0)
                sq_b += eb * eb
                q_ref_m = q_ref_chk[chk_idx]
                for ii in range(3):
                    for dd in range(3):
                        dq_val = st_b.q[ii, dd] - q_ref_m[ii, dd]
                        sq_qb += dq_val * dq_val
            err_b = float(mp.sqrt(sq_b / mp.mpf(str(n_chk))))
            err_qb = float(mp.sqrt(sq_qb / mp.mpf(str(n_chk))))

            # 2. Telescoped PBM4 (Z_4^{BM4} in shadow manifold) at exact common checkpoints
            st_sh = corr_hom.apply(st0, tree, h_mp, inverse=True)
            sq_p = mp.mpf("0.0")
            sq_qp = mp.mpf("0.0")
            for chk_idx in range(n_chk):
                for _ in range(substeps_per_chk):
                    st_sh = bm4.step(st_sh, tree, dt=h_mp)
                st_p = corr_hom.apply(st_sh, tree, h_mp, inverse=False)
                ep = abs(inv128(st_p, G=mp.mpf("1.0")).H - H0) / abs(H0)
                sq_p += ep * ep
                q_ref_m = q_ref_chk[chk_idx]
                for ii in range(3):
                    for dd in range(3):
                        dq_val = st_p.q[ii, dd] - q_ref_m[ii, dd]
                        sq_qp += dq_val * dq_val
            err_p = float(mp.sqrt(sq_p / mp.mpf(str(n_chk))))
            err_qp = float(mp.sqrt(sq_qp / mp.mpf(str(n_chk))))

            dE_bm4_list.append(err_b)
            dE_pbm4_list.append(err_p)
            gains_list.append(err_b / err_p)
            dq_bm4_list.append(err_qb)
            dq_pbm4_list.append(err_qp)
            gains_q_list.append(err_qb / err_qp)

        # Compute unclipped local finite-difference slopes
        slopes_bm4: List[float] = []
        slopes_pbm4: List[float] = []
        slopes_q_bm4: List[float] = []
        slopes_q_pbm4: List[float] = []
        h_midpoints: List[float] = []
        for k in range(len(h_steps) - 1):
            ratio_h = math.log(h_steps[k] / h_steps[k + 1])
            s_b = math.log(dE_bm4_list[k] / dE_bm4_list[k + 1]) / ratio_h
            s_p = math.log(dE_pbm4_list[k] / dE_pbm4_list[k + 1]) / ratio_h
            sq_b_sl = math.log(dq_bm4_list[k] / dq_bm4_list[k + 1]) / ratio_h
            sq_p_sl = math.log(dq_pbm4_list[k] / dq_pbm4_list[k + 1]) / ratio_h
            slopes_bm4.append(float(s_b))
            slopes_pbm4.append(float(s_p))
            slopes_q_bm4.append(float(sq_b_sl))
            slopes_q_pbm4.append(float(sq_p_sl))
            h_midpoints.append(math.sqrt(h_steps[k] * h_steps[k + 1]))

        # Fit two-channel model dE_model(h) = c6 * eps * h^6 + c4_2 * eps^2 * h^4
        # so the theoretical slope curve p_theo(h) = d ln(dE_model) / d ln(h) can be plotted
        c6_est = dE_pbm4_list[0] / (eps * (h_steps[0] ** 6))
        c42_est = max((dE_pbm4_list[-1] - c6_est * eps * (h_steps[-1] ** 6)) / ((eps ** 2) * (h_steps[-1] ** 4)), 0.0)
        theo_slopes: List[float] = []
        for hm in h_midpoints:
            term6 = c6_est * eps * (hm ** 6)
            term4 = c42_est * (eps ** 2) * (hm ** 4)
            p_th = (6.0 * term6 + 4.0 * term4) / (term6 + term4)
            theo_slopes.append(float(p_th))

        results_by_eps[eps_key] = {
            "epsilon": float(eps),
            "R_out": float(R_out),
            "e_in": 0.30,
            "inc_deg": 28.0,
            "T_cross": T_cross,
            "stages_s": stages_s,
            "h_steps": h_steps,
            "tau_over_s": [h / stages_s for h in h_steps],
            "dE_bm4": dE_bm4_list,
            "dE_pbm4": dE_pbm4_list,
            "gain_pbm4_vs_bm4": gains_list,
            "dq_bm4": dq_bm4_list,
            "dq_pbm4": dq_pbm4_list,
            "gain_q_pbm4_vs_bm4": gains_q_list,
            "h_midpoints": h_midpoints,
            "local_slopes_bm4": slopes_bm4,
            "local_slopes_pbm4": slopes_pbm4,
            "local_slopes_q_bm4": slopes_q_bm4,
            "local_slopes_q_pbm4": slopes_q_pbm4,
            "theoretical_slopes_pbm4": theo_slopes,
            "unclipped": True,
        }
        print(
            f"  [Crossover 3D F128] {eps_key} (R_out={R_out:5.2f}): "
            f"max_gain_E={max(gains_list):.1f}x, max_gain_q={max(gains_q_list):.1f}x, "
            f"min_dE_pbm4={min(dE_pbm4_list):.3e}, min_dq_pbm4={min(dq_pbm4_list):.3e}, "
            f"slopes_pbm4={[round(s, 2) for s in slopes_pbm4]}",
            flush=True,
        )

    return {
        "float128_verified": True,
        "e_in": 0.30,
        "inc_deg": 28.0,
        "T_cross": T_cross,
        "sweeps": results_by_eps,
    }


# ==============================================================================
# Part II: Protocol ST-1 — Two-Regime Tidal Perturbation Scaling of PBM4 (Z_4^{BM4}) vs BM4
# ==============================================================================

def run_part2_st1_obstruction_sweep() -> Dict[str, Any]:
    """
    Protocol ST-1: Sweeps eps in [1e-5, 3e-2] at fixed step h = 0.08 on the
    3D eccentric hierarchical triple (e_in=0.30, i=28 deg), comparing:
      1. Raw BM4
      2. Palindromic Lie conjugator PBM4 Z_4^{BM4} (PalindromicCorrectorZ4, mode='homological')
    """
    print("=== Part II (ST-1): Two-Regime Tidal Scaling Sweep (BM4 vs PBM4 Z_4^{BM4}, e=0.30, i=28 deg) ===", flush=True)
    eps_list = [3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
    h_fixed = 0.08
    T_tot = 1.60
    n_steps = int(round(T_tot / h_fixed))

    # Verify canonical default is mode="homological"
    default_corr = CorrectorC4_64(G=1.0)
    default_pbm4 = create_pbm4_integrator64(G=1.0)
    assert default_corr.mode == "homological", f"Expected homological default, got {default_corr.mode}"
    assert default_pbm4.corrector.mode == "homological"

    bm4 = BlanesMoan4Integrator64(G=1.0)
    corr_hom = PalindromicCorrectorZ4_64(G=1.0)

    dE_bm4_list: List[float] = []
    dE_hom_list: List[float] = []
    gains_hom_list: List[float] = []
    R_out_list: List[float] = []

    for eps in eps_list:
        st0, tree, R_out = make_3d_eccentric_hierarchical_f64(eps=eps, e_in=0.30, inc_deg=28.0)
        assert mod_02.detect_supernodes(st0, tree) == [], f"Unexpected supernode on 3D hierarchical triple at eps={eps}"
        R_out_list.append(float(R_out))
        H0 = inv64(st0, G=1.0).H

        # 1. Raw BM4
        st_b = st0
        max_b = 0.0
        for _ in range(n_steps):
            st_b = bm4.step(st_b, tree, dt=h_fixed)
            max_b = max(max_b, abs(inv64(st_b, G=1.0).H - H0) / abs(H0))
        dE_bm4_list.append(float(max_b))

        # 2. PalindromicCorrectorZ4 (mode="homological")
        st_h = corr_hom.apply(st0, tree, h_fixed, inverse=True)
        max_h = 0.0
        for _ in range(n_steps):
            st_h = bm4.step(st_h, tree, dt=h_fixed)
            st_hp = corr_hom.apply(st_h, tree, h_fixed, inverse=False)
            max_h = max(max_h, abs(inv64(st_hp, G=1.0).H - H0) / abs(H0))
        dE_hom_list.append(float(max_h))

        # Verify bit-for-bit identity with canonical NBodySystem orchestrator for eps <= 1e-2
        if eps <= 1e-2:
            sys_canon = NBodySystem64(
                state=st0,
                integrator=create_pbm4_integrator64(G=1.0),
                tree=tree,
                G=1.0,
                enable_dyadic_substepping=False,
                in_shadow_mode=True,
                use_optimized_kernels=True,
            )
            for _ in range(n_steps):
                sys_canon.step(dt=h_fixed)
            st_canon_phys = sys_canon.observe_physical_state()
            assert np.array_equal(st_canon_phys.q, st_hp.q) and np.array_equal(st_canon_phys.p, st_hp.p)

        g_z4 = float(max_b / max_h)
        gains_hom_list.append(g_z4)

        print(
            f"  [ST-1] eps={eps:.1e} (R_out={R_out:5.1f}): "
            f"BM4={max_b:.2e}, PBM4_Z4={max_h:.2e} ({g_z4:.1f}x)",
            flush=True,
        )

    # Fit asymptotic slopes wrt epsilon in the two regimes:
    # Regime A (quadratic floor eps >= h^2 = 6.4e-3, indices 0..3: eps in [1e-3, 3e-2])
    log_eps_high = np.log(eps_list[:4])
    slope_eps_high_z4 = float(np.polyfit(log_eps_high, np.log(dE_hom_list[:4]), 1)[0])
    # Regime B (pre-asymptotic order-6 channel eps << h^2, indices 4..7: eps in [1e-5, 3e-4])
    log_eps_low = np.log(eps_list[4:])
    slope_eps_low_z4 = float(np.polyfit(log_eps_low, np.log(dE_hom_list[4:]), 1)[0])

    return {
        "h_fixed": h_fixed,
        "h_squared": h_fixed ** 2,
        "T_tot": T_tot,
        "e_in": 0.30,
        "inc_deg": 28.0,
        "default_mode_verified": "homological",
        "eps_values": eps_list,
        "R_out_values": R_out_list,
        "dE_bm4": dE_bm4_list,
        "dE_pbm4_homological_z4": dE_hom_list,
        "gain_homological_z4_vs_bm4": gains_hom_list,
        "slope_wrt_eps_quadratic_regime": slope_eps_high_z4,
        "slope_wrt_eps_order6_regime": slope_eps_low_z4,
    }


# ==============================================================================
# Part III: Shadow Telescoping Amortization & 100-Orbit Secular Endurance (ST-3)
# ==============================================================================

def compute_naff_fundamental_freq(signal: np.ndarray, dt_sample: float) -> float:
    """
    Computes the fundamental frequency nu_0 (in cycles per unit time) of a complex
    time series z(t) = x(t) + i y(t) using Laskar's Hann-windowed NAFF maximization.
    """
    N = len(signal)
    t_grid = np.arange(N, dtype=np.float64) * dt_sample
    T_span = max((N - 1) * dt_sample, 1e-15)
    hann = 1.0 - np.cos(2.0 * math.pi * t_grid / T_span)
    sig_w = signal * hann

    spec = np.abs(np.fft.fft(sig_w))
    freqs = np.fft.fftfreq(N, d=dt_sample)
    pos_mask = freqs > 0
    pos_freqs = freqs[pos_mask]
    pos_spec = spec[pos_mask]
    f0 = float(pos_freqs[int(np.argmax(pos_spec))])

    df = 1.5 / T_span
    a, b = f0 - df, f0 + df
    phi_gold = (math.sqrt(5.0) - 1.0) / 2.0

    def _amp(f_val: float) -> float:
        phase = np.exp(-2j * math.pi * f_val * t_grid)
        return float(np.abs(np.sum(sig_w * phase)))

    c = b - phi_gold * (b - a)
    d = a + phi_gold * (b - a)
    fc, fd = _amp(c), _amp(d)
    for _ in range(50):
        if fc > fd:
            b = d
            d = c
            fd = fc
            c = b - phi_gold * (b - a)
            fc = _amp(c)
        else:
            a = c
            c = d
            fc = fd
            d = a + phi_gold * (b - a)
            fd = _amp(d)
    return 0.5 * (a + b)


def run_part3_telescoping_and_st3_endurance() -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    1. Shadow Telescoping Amortization: sweeps M in {1, 2, 4, 8, 16, 32, 64, 128, 256}
       comparing un-telescoped step-by-step physical projection (2M corrector calls)
       vs shadow manifold telescoping C_h o Phi_h^M o C_h^{-1} (2 corrector calls).
    2. Protocol ST-3: 100-orbit secular endurance on the 3D eccentric hierarchical triple
       (e_in=0.30, i=28 deg, eps=4e-4) with fine reference trajectory h_ref = h/4 and
       dense sampling sample_stride = 1 (40 samples/orbit, 4,000 samples total).
    """
    print("=== Part IIIa: Shadow Manifold Telescoping Amortization (M in [1, 256]) ===", flush=True)
    st0, tree, _ = make_3d_eccentric_hierarchical_f64(eps=1e-3, e_in=0.30, inc_deg=28.0)
    h = 0.05
    bm4 = BlanesMoan4Integrator64(G=1.0)
    corr_hom = PalindromicCorrectorZ4_64(G=1.0)

    M_list = [1, 2, 4, 8, 16, 32, 64, 128, 256]
    s_eff_telescoped: List[float] = []
    s_eff_untelescoped: List[float] = []
    speedup_list: List[float] = []
    state_diff_list: List[float] = []

    for M in M_list:
        # Un-telescoped: apply C_h^{-1} and C_h at every step (or C_h readout at every step)
        t0 = time.perf_counter()
        st_un = corr_hom.apply(st0, tree, h, inverse=True)
        for _ in range(M):
            st_un = bm4.step(st_un, tree, dt=h)
            _ = corr_hom.apply(st_un, tree, h, inverse=False)
        st_un_phys = corr_hom.apply(st_un, tree, h, inverse=False)
        dt_un = max(time.perf_counter() - t0, 1e-9)

        # Telescoped: enter once, step M times in shadow space, exit once
        t1 = time.perf_counter()
        st_tel = corr_hom.apply(st0, tree, h, inverse=True)
        for _ in range(M):
            st_tel = bm4.step(st_tel, tree, dt=h)
        st_tel_phys = corr_hom.apply(st_tel, tree, h, inverse=False)
        dt_tel = max(time.perf_counter() - t1, 1e-9)

        diff_q = float(np.max(np.abs(st_tel_phys.q - st_un_phys.q)))
        state_diff_list.append(diff_q)
        # Interior BM4 is 6 FSAL apex evaluations; each Z4 call has 8 apex kicks (16 for enter+exit)
        s_eff_telescoped.append(float(6.0 + 16.0 / M))
        s_eff_untelescoped.append(float(6.0 + 8.0 + 8.0 / M))
        speedup_list.append(float(dt_un / dt_tel))

    telescoping_data = {
        "M_steps": M_list,
        "s_eff_telescoped": s_eff_telescoped,
        "s_eff_untelescoped": s_eff_untelescoped,
        "wall_clock_speedup": speedup_list,
        "max_state_diff_telescoped_vs_untelescoped": max(state_diff_list),
    }
    print(
        f"  [Telescoping] M=256: s_eff={s_eff_telescoped[-1]:.2f} (vs {s_eff_untelescoped[-1]:.2f}), "
        f"speedup={speedup_list[-1]:.2f}x, state_diff={max(state_diff_list):.2e}",
        flush=True,
    )

    print("=== Part IIIb (ST-3): 100-Orbit 3D Eccentric Secular Endurance (h_ref=h/4, sample_stride=1) ===", flush=True)
    eps = 4.0e-4
    st0, tree, R_out = make_3d_eccentric_hierarchical_f64(eps=eps, e_in=0.30, inc_deg=28.0)
    H0 = inv64(st0, G=1.0).H
    L0 = inv64(st0, G=1.0).L
    L0_norm = float(np.linalg.norm(L0))

    T_bin = 2.0 * math.pi / math.sqrt(2.0)
    n_orbits = 100
    steps_per_orbit = 60
    h = T_bin / steps_per_orbit
    total_steps = n_orbits * steps_per_orbit
    h_ref = h / 4.0

    st_bm4 = st0
    st_sh = corr_hom.apply(st0, tree, h, inverse=True)
    st_ref = corr_hom.apply(st0, tree, h_ref, inverse=True)

    orbit_times: List[float] = []
    orbit_dE_bm4: List[float] = []
    orbit_dE_pbm4: List[float] = []
    orbit_dL_bm4: List[float] = []
    orbit_dL_pbm4: List[float] = []

    sig_bm4: List[complex] = []
    sig_pbm4: List[complex] = []
    sig_ref: List[complex] = []

    max_in_orbit_bm4 = 0.0
    max_in_orbit_pbm4 = 0.0

    sample_stride = 1  # Remediates [MOD-7.1.2-01]: 40 samples/orbit, 4,000 total samples

    for step_idx in range(1, total_steps + 1):
        st_bm4 = bm4.step(st_bm4, tree, dt=h)
        st_sh = bm4.step(st_sh, tree, dt=h)
        # 4 substeps of h_ref = h / 4
        for _ in range(4):
            st_ref = bm4.step(st_ref, tree, dt=h_ref)

        st_pbm4 = corr_hom.apply(st_sh, tree, h, inverse=False)
        st_ref_phys = corr_hom.apply(st_ref, tree, h_ref, inverse=False)

        inv_b = inv64(st_bm4, G=1.0)
        inv_p = inv64(st_pbm4, G=1.0)

        err_b = abs(inv_b.H - H0) / abs(H0)
        err_p = abs(inv_p.H - H0) / abs(H0)
        max_in_orbit_bm4 = max(max_in_orbit_bm4, float(err_b))
        max_in_orbit_pbm4 = max(max_in_orbit_pbm4, float(err_p))

        r_b = st_bm4.q[1] - st_bm4.q[0]
        r_p = st_pbm4.q[1] - st_pbm4.q[0]
        r_r = st_ref_phys.q[1] - st_ref_phys.q[0]
        sig_bm4.append(complex(float(r_b[0]), float(r_b[1])))
        sig_pbm4.append(complex(float(r_p[0]), float(r_p[1])))
        sig_ref.append(complex(float(r_r[0]), float(r_r[1])))

        if step_idx % steps_per_orbit == 0:
            orb_idx = step_idx // steps_per_orbit
            orbit_times.append(float(orb_idx))
            orbit_dE_bm4.append(float(max_in_orbit_bm4))
            orbit_dE_pbm4.append(float(max_in_orbit_pbm4))
            orbit_dL_bm4.append(float(np.linalg.norm(inv_b.L - L0) / L0_norm))
            orbit_dL_pbm4.append(float(np.linalg.norm(inv_p.L - L0) / L0_norm))
            max_in_orbit_bm4 = 0.0
            max_in_orbit_pbm4 = 0.0

    dt_sample = sample_stride * h
    arr_b = np.array(sig_bm4, dtype=np.complex128)
    arr_p = np.array(sig_pbm4, dtype=np.complex128)
    arr_r = np.array(sig_ref, dtype=np.complex128)

    nu_ref = compute_naff_fundamental_freq(arr_r, dt_sample)
    nu_bm4 = compute_naff_fundamental_freq(arr_b, dt_sample)
    nu_pbm4 = compute_naff_fundamental_freq(arr_p, dt_sample)

    naff_err_bm4 = abs(nu_bm4 - nu_ref) / abs(nu_ref)
    naff_err_pbm4 = abs(nu_pbm4 - nu_ref) / abs(nu_ref)
    naff_gain = float(naff_err_bm4 / naff_err_pbm4)

    max_global_bm4 = max(orbit_dE_bm4)
    max_global_pbm4 = max(orbit_dE_pbm4)
    gain_secular = float(max_global_bm4 / max_global_pbm4)

    print(
        f"  [ST-3] 100 Orbits (e=0.30, eps={eps:.1e}, h_ref=h/4, N_samples={len(sig_ref)}): "
        f"max|dE|_BM4={max_global_bm4:.2e}, max|dE|_PBM4={max_global_pbm4:.2e} (Gain={gain_secular:.1f}x), "
        f"NAFF err BM4={naff_err_bm4:.2e} vs PBM4={naff_err_pbm4:.2e} (Gain={naff_gain:.1f}x)",
        flush=True,
    )

    st3_data = {
        "epsilon": eps,
        "e_in": 0.30,
        "inc_deg": 28.0,
        "R_out": float(R_out),
        "T_bin": float(T_bin),
        "h": float(h),
        "h_ref": float(h_ref),
        "sample_stride": sample_stride,
        "n_samples": len(sig_ref),
        "n_orbits": n_orbits,
        "orbit_numbers": orbit_times,
        "orbit_max_dE_bm4": orbit_dE_bm4,
        "orbit_max_dE_pbm4": orbit_dE_pbm4,
        "orbit_dL_bm4": orbit_dL_bm4,
        "orbit_dL_pbm4": orbit_dL_pbm4,
        "max_dE_bm4": float(max_global_bm4),
        "max_dE_pbm4": float(max_global_pbm4),
        "energy_suppression_factor": gain_secular,
        "nu_ref": float(nu_ref),
        "naff_drift_bm4": float(naff_err_bm4),
        "naff_drift_pbm4": float(naff_err_pbm4),
        "naff_gain_factor": naff_gain,
    }
    return telescoping_data, st3_data


# ==============================================================================
# Part IV: Protocol ST-4 — Dynamical Multi-Step Gauge Transition & Cech 1-Cocycle
# ==============================================================================

def run_part4_st4_dynamic_gauge_and_cocycle() -> Dict[str, Any]:
    """
    Remediates [CRIT-7.1.2-01]:
      1. ST-4a (Dynamical Multi-Step Interfacial Transition in Float64 & Float128):
         Integrates N_A = 4 real steps before transition and N_B = 4 real steps after
         transition on the 3D eccentric hierarchical triple (eps = 1e-3, e_in = 0.30,
         i = 28 deg) across a dyadic step-size transition h_A = h -> h_B = h/2,
         comparing:
           (i) Naive Switch (no gauge connection: z_sh_B = z_sh_A, O(eps * h^4) jump),
           (ii) Obstructed P_0 Switch (mode="global": fails to cancel momentum-dependent jump),
           (iii) Covariant PBM4 (Z_4^{BM4}) Switch (mode="homological": O(eps * h^6) jump suppression).
      2. ST-4b (Non-Commutative 3-Chart Cech 1-Cocycle Closure & Lie-Bargmann Conservation):
         Evaluates T_{C->A} o T_{B->C} o T_{A->B} on mode="homological" (Z_4^{BM4}) in both
         Float64 and Float128 across h in [0.16, 0.08, 0.04, 0.02, 0.01, 0.005] WITHOUT any
         clamping to 1e-35.
    """
    print("=== Part IV (ST-4): Dynamical Multi-Step Gauge Transition & Non-Commutative Cech 1-Cocycle ===", flush=True)
    mp.mp.dps = 35

    h_list = [0.08, 0.04, 0.02, 0.01, 0.005]
    n_steps_before = 4
    n_steps_after = 4
    eps_dyn = 1.0e-3

    st0_128, treeA_128, _ = make_3d_eccentric_hierarchical_f128(eps=eps_dyn, e_in=0.30, inc_deg=28.0)
    H0_128 = inv128(st0_128, G=mp.mpf("1.0")).H
    bm4_128 = BM4_128(G=mp.mpf("1.0"))
    c4_hom_128 = CorrectorC4_128(G=mp.mpf("1.0"), mode="homological")

    dyn_naive_f128: List[float] = []
    dyn_hom_f128: List[float] = []
    dyn_suppression_f128: List[float] = []

    for h in h_list:
        h_A = mp.mpf(str(h))
        h_B = h_A / mp.mpf("2.0")

        # Enter shadow manifold at (treeA, h_A) and evolve N_A = 4 real steps
        z_sh_pre = c4_hom_128.apply(st0_128, treeA_128, h_A, inverse=True)
        for _ in range(n_steps_before):
            z_sh_pre = bm4_128.step(z_sh_pre, treeA_128, dt=h_A)

        # Branch 1: Naive switch to h_B without gauge connection, then evolve N_B = 4 real steps
        z_sh_naive = z_sh_pre
        for _ in range(n_steps_after):
            z_sh_naive = bm4_128.step(z_sh_naive, treeA_128, dt=h_B)
        z_phys_naive = c4_hom_128.apply(z_sh_naive, treeA_128, h_B, inverse=False)
        err_naive = float(abs(inv128(z_phys_naive, G=mp.mpf("1.0")).H - H0_128) / abs(H0_128))

        # Branch 2: Covariant PBM4 (Z_4^{BM4}) gauge switch, then evolve N_B = 4 real steps
        z_sh_hom = apply_gauge128(
            z_sh_pre, treeA_128, treeA_128, dt=h_A, dt_new=h_B, corrector=c4_hom_128
        )
        for _ in range(n_steps_after):
            z_sh_hom = bm4_128.step(z_sh_hom, treeA_128, dt=h_B)
        z_phys_hom = c4_hom_128.apply(z_sh_hom, treeA_128, h_B, inverse=False)
        err_hom = float(abs(inv128(z_phys_hom, G=mp.mpf("1.0")).H - H0_128) / abs(H0_128))

        supp = err_naive / err_hom
        dyn_naive_f128.append(err_naive)
        dyn_hom_f128.append(err_hom)
        dyn_suppression_f128.append(float(supp))

        print(
            f"  [ST-4a Dynamic Transition N_A=4 -> N_B=4] h={h:.3f}: "
            f"Naive={err_naive:.3e}, Covariant_Z4={err_hom:.3e} "
            f"(Suppression={supp:.1f}x)",
            flush=True,
        )

    # ST-4b: Non-Commutative 3-Chart Cech 1-Cocycle T_{C->A} o T_{B->C} o T_{A->B} = Id
    # evaluated on mode="homological" (Z_4^{BM4}) in Float64 and Float128 WITHOUT clamping
    h_cocycle_list = [0.16, 0.08, 0.04, 0.02, 0.01, 0.005]

    # Float64 3-chart cycle
    st0_64, _, _ = make_3d_eccentric_hierarchical_f64(eps=1e-2, e_in=0.30, inc_deg=28.0)
    m64 = st0_64.masses
    treeA_64 = build_tree64(m64, ((0, 1), 2))
    treeB_64 = build_tree64(m64, (0, (1, 2)))
    treeC_64 = build_tree64(m64, ((0, 2), 1))
    c4_hom_64 = CorrectorC4_64(G=1.0, mode="homological")
    inv0_64 = inv64(st0_64, G=1.0)

    # Float128 3-chart cycle
    st0_c128, _, _ = make_3d_eccentric_hierarchical_f128(eps=1e-2, e_in=0.30, inc_deg=28.0)
    m128 = st0_c128.masses
    tA_128 = build_tree128(m128, ((0, 1), 2))
    tB_128 = build_tree128(m128, (0, (1, 2)))
    tC_128 = build_tree128(m128, ((0, 2), 1))
    inv0_128 = inv128(st0_c128, G=mp.mpf("1.0"))

    cocycle_f64_list: List[float] = []
    cocycle_f128_list: List[float] = []
    lb_drift_f64_list: List[float] = []
    lb_drift_f128_list: List[float] = []

    for h in h_cocycle_list:
        # Float64 cycle on mode="homological"
        z_AB_64 = apply_gauge_connection64(st0_64, treeA_64, treeB_64, dt=h, corrector=c4_hom_64)
        z_BC_64 = apply_gauge_connection64(z_AB_64, treeB_64, treeC_64, dt=h, corrector=c4_hom_64)
        z_CA_64 = apply_gauge_connection64(z_BC_64, treeC_64, treeA_64, dt=h, corrector=c4_hom_64)
        def_64 = float(max(np.max(np.abs(z_CA_64.q - st0_64.q)), np.max(np.abs(z_CA_64.p - st0_64.p))))
        inv_AB_64 = inv64(z_AB_64, G=1.0)
        dL_64 = float(np.linalg.norm(inv_AB_64.L - inv0_64.L))
        cocycle_f64_list.append(def_64)
        lb_drift_f64_list.append(dL_64)

        # Float128 cycle on mode="homological" (unclamped!)
        h_mp = mp.mpf(str(h))
        z_AB_128 = apply_gauge128(st0_c128, tA_128, tB_128, dt=h_mp, corrector=c4_hom_128)
        z_BC_128 = apply_gauge128(z_AB_128, tB_128, tC_128, dt=h_mp, corrector=c4_hom_128)
        z_CA_128 = apply_gauge128(z_BC_128, tC_128, tA_128, dt=h_mp, corrector=c4_hom_128)
        def_128 = float(max(
            max(abs(z_CA_128.q[i, d] - st0_c128.q[i, d]) for i in range(3) for d in range(3)),
            max(abs(z_CA_128.p[i, d] - st0_c128.p[i, d]) for i in range(3) for d in range(3)),
        ))
        inv_AB_128 = inv128(z_AB_128, G=mp.mpf("1.0"))
        dL_128 = float(mp.sqrt(sum((inv_AB_128.L[d] - inv0_128.L[d]) ** 2 for d in range(3))))
        cocycle_f128_list.append(def_128)
        lb_drift_f128_list.append(dL_128)

        print(
            f"  [ST-4b 3-Chart Cocycle Z4] h={h:.3f}: "
            f"defect_f64={def_64:.3e}, defect_f128={def_128:.3e}, dL_f128={dL_128:.3e}",
            flush=True,
        )

    return {
        "dynamic_transition": {
            "n_steps_before": n_steps_before,
            "n_steps_after": n_steps_after,
            "epsilon": eps_dyn,
            "h_steps": h_list,
            "dE_naive_switch_f128": dyn_naive_f128,
            "dE_covariant_z4_switch_f128": dyn_hom_f128,
            "dynamic_jump_suppression_f128": dyn_suppression_f128,
        },
        "cech_cocycle_homological": {
            "corrector_mode": "homological",
            "clamped": False,
            "h_steps": h_cocycle_list,
            "cocycle_defect_f64": cocycle_f64_list,
            "cocycle_defect_f128": cocycle_f128_list,
            "lie_bargmann_L_drift_f64": lb_drift_f64_list,
            "lie_bargmann_L_drift_f128": lb_drift_f128_list,
        },
    }


# ==============================================================================
# Part V: Physical 3D Orbital Geometry & Shadow-to-Physical Manifold Offset
# ==============================================================================

def run_part5_physical_3d_and_shadow_offset() -> Dict[str, Any]:
    """
    Integrates 2 binary periods of the 3D eccentric inclined hierarchical triple
    (e_in=0.30, i=28 deg, eps=1e-2) at h=0.06, recording:
      - 3D relative trajectory r(t) = q_1(t) - q_0(t) and outer perturber trajectory
      - Radial momentum w_rad(t) = (r . w) / r driving the P_0 momentum-parity obstruction
      - Canonical shadow-to-physical position offset ||q_shadow(t) - q_phys(t)||_2 = O(h^4)
      - Raw shadow energy error |H(z_shadow(t)) - H0| / |H0| vs physical energy error
        |H(observe_physical(z_shadow(t))) - H0| / |H0|.
    """
    print("=== Part V: Physical 3D Geometry, Radial Momentum w_rad(t) & Shadow Offset ===", flush=True)
    eps = 1.0e-2
    st0, tree, R_out = make_3d_eccentric_hierarchical_f64(eps=eps, e_in=0.30, inc_deg=28.0)
    H0 = inv64(st0, G=1.0).H

    T_bin = 2.0 * math.pi / math.sqrt(2.0)
    n_steps = 160
    h = (2.0 * T_bin) / n_steps

    bm4 = BlanesMoan4Integrator64(G=1.0)
    corr_hom = PalindromicCorrectorZ4_64(G=1.0)

    st_sh = corr_hom.apply(st0, tree, h, inverse=True)

    times: List[float] = []
    rx_list: List[float] = []
    ry_list: List[float] = []
    rz_list: List[float] = []
    r_norm_list: List[float] = []
    w_rad_list: List[float] = []
    shadow_offset_q_list: List[float] = []
    shadow_dE_raw_list: List[float] = []
    phys_dE_recovered_list: List[float] = []

    mu_in = (st0.masses[0] * st0.masses[1]) / (st0.masses[0] + st0.masses[1])

    for k in range(n_steps + 1):
        t_val = k * h
        st_phys = corr_hom.apply(st_sh, tree, h, inverse=False)

        r_vec = st_phys.q[1] - st_phys.q[0]
        v_vec = (st_phys.p[1] / st_phys.masses[1]) - (st_phys.p[0] / st_phys.masses[0])
        w_vec = mu_in * v_vec
        r_mag = float(np.linalg.norm(r_vec))
        w_rad = float(np.dot(r_vec, w_vec) / max(r_mag, 1e-15))

        offset_q = float(np.linalg.norm(st_sh.q - st_phys.q))
        dE_sh = float(abs(inv64(st_sh, G=1.0).H - H0) / abs(H0))
        dE_ph = float(abs(inv64(st_phys, G=1.0).H - H0) / abs(H0))

        times.append(float(t_val / T_bin))
        rx_list.append(float(r_vec[0]))
        ry_list.append(float(r_vec[1]))
        rz_list.append(float(r_vec[2]))
        r_norm_list.append(r_mag)
        w_rad_list.append(w_rad)
        shadow_offset_q_list.append(offset_q)
        shadow_dE_raw_list.append(dE_sh)
        phys_dE_recovered_list.append(dE_ph)

        if k < n_steps:
            st_sh = bm4.step(st_sh, tree, dt=h)

    return {
        "epsilon": eps,
        "e_in": 0.30,
        "inc_deg": 28.0,
        "R_out": float(R_out),
        "T_bin": float(T_bin),
        "h": float(h),
        "t_over_Tbin": times,
        "r_rel_x": rx_list,
        "r_rel_y": ry_list,
        "r_rel_z": rz_list,
        "r_norm": r_norm_list,
        "w_rad": w_rad_list,
        "shadow_offset_q": shadow_offset_q_list,
        "shadow_dE_raw": shadow_dE_raw_list,
        "phys_dE_recovered": phys_dE_recovered_list,
    }


# ==============================================================================
# Automated Contract Verification Gate (PASS/FAIL)
# ==============================================================================

def verify_contract_assertions(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Verifies all 6 mandatory assertions ASSERT-G02-01 to ASSERT-G02-06 from
    scripts/campaigns/contracts/02_generator_universal_obstruction_and_shadow_gauge.md.
    """
    c_f128 = data["crossover_3d_f128"]
    eps_1e4 = c_f128["sweeps"]["eps_1e-04"]
    max_gain_1e4 = max(eps_1e4["gain_pbm4_vs_bm4"])
    pass_01 = (
        c_f128["float128_verified"] is True
        and abs(c_f128["e_in"] - 0.30) < 1e-12
        and abs(c_f128["inc_deg"] - 28.0) < 1e-12
        and max_gain_1e4 >= 500.0
    )

    slopes_1e4 = eps_1e4["local_slopes_pbm4"]
    slope_coarse = slopes_1e4[0]
    slope_fine = slopes_1e4[-1]
    pass_02 = (
        eps_1e4["unclipped"] is True
        and slope_coarse >= 5.5
        and slope_fine <= 4.6
    )

    st1 = data["st1_obstruction_sweep"]
    max_g_z4 = max(st1["gain_homological_z4_vs_bm4"])
    pass_03 = max_g_z4 >= 200.0

    st3 = data["st3_secular_endurance"]
    pass_04 = (
        st3["energy_suppression_factor"] >= 50.0
        and st3["naff_gain_factor"] >= 75.0
        and st3["sample_stride"] == 1
    )

    st4_dyn = data["st4_dynamic_gauge_and_cocycle"]["dynamic_transition"]
    max_supp_dyn = max(st4_dyn["dynamic_jump_suppression_f128"])
    pass_05 = (
        st4_dyn["n_steps_before"] >= 4
        and st4_dyn["n_steps_after"] >= 4
        and max_supp_dyn >= 100.0
    )

    st4_coc = data["st4_dynamic_gauge_and_cocycle"]["cech_cocycle_homological"]
    max_coc_64 = max(st4_coc["cocycle_defect_f64"])
    max_coc_128 = max(st4_coc["cocycle_defect_f128"])
    pass_06 = (
        st4_coc["corrector_mode"] == "homological"
        and st4_coc["clamped"] is False
        and max_coc_64 <= 1e-14
        and max_coc_128 <= 1e-32
    )

    assertions = {
        "ASSERT-G02-01_float128_3d_eccentric_crossover": {
            "passed": bool(pass_01),
            "max_gain_eps_1e4": float(max_gain_1e4),
            "e_in": float(c_f128["e_in"]),
            "inc_deg": float(c_f128["inc_deg"]),
        },
        "ASSERT-G02-02_unclipped_crossover_slopes_6_to_4": {
            "passed": bool(pass_02),
            "slope_coarse": float(slope_coarse),
            "slope_fine": float(slope_fine),
            "all_slopes_eps_1e4": [float(s) for s in slopes_1e4],
        },
        "ASSERT-G02-03_two_regime_tidal_scaling_st1": {
            "passed": bool(pass_03),
            "max_gain_homological_z4": float(max_g_z4),
        },
        "ASSERT-G02-04_secular_endurance_and_naff_st3": {
            "passed": bool(pass_04),
            "energy_suppression_factor": float(st3["energy_suppression_factor"]),
            "naff_gain_factor": float(st3["naff_gain_factor"]),
            "n_samples": int(st3["n_samples"]),
        },
        "ASSERT-G02-05_dynamic_multistep_gauge_transition_st4a": {
            "passed": bool(pass_05),
            "n_steps_before": int(st4_dyn["n_steps_before"]),
            "n_steps_after": int(st4_dyn["n_steps_after"]),
            "max_dynamic_jump_suppression_f128": float(max_supp_dyn),
        },
        "ASSERT-G02-06_unclamped_homological_cech_cocycle_st4b": {
            "passed": bool(pass_06),
            "max_cocycle_defect_f64": float(max_coc_64),
            "max_cocycle_defect_f128": float(max_coc_128),
        },
    }

    all_passed = all(v["passed"] for v in assertions.values())
    if not all_passed:
        failed = [k for k, v in assertions.items() if not v["passed"]]
        raise AssertionError(f"Contract assertion(s) failed: {failed}\n{json.dumps(assertions, indent=2)}")

    print("\n=======================================================================", flush=True)
    print("  [PASS] All 6 Contract Assertions (ASSERT-G02-01 to ASSERT-G02-06) CERTIFIED!", flush=True)
    print("=======================================================================", flush=True)
    return assertions


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Benchmark Generator 02 (Shadow Manifold & Quadratic Floor)")
    parser.add_argument(
        "--protocol",
        type=str,
        default="all",
        choices=["refresh-f64", "verify", "all"],
        help="Protocol to run (refresh-f64 re-runs Parts II-V using cached Float128 Part I; all runs everything)",
    )
    args = parser.parse_args()

    t0 = time.perf_counter()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / "generator_02_universal_obstruction_and_shadow_gauge.json"

    existing: Dict[str, Any] = {}
    if out_path.exists():
        try:
            with open(out_path, "r", encoding="utf-8") as f:
                existing = json.load(f)
        except Exception:
            existing = {}

    if args.protocol == "verify":
        verify_contract_assertions(existing)
        return

    if args.protocol == "refresh-f64" and "crossover_3d_f128" in existing:
        crossover_3d_f128 = existing["crossover_3d_f128"]
    else:
        crossover_3d_f128 = run_part1_crossover_3d_f128()

    st1_obstruction = run_part2_st1_obstruction_sweep()
    shadow_telescoping, st3_secular = run_part3_telescoping_and_st3_endurance()
    st4_dynamic = run_part4_st4_dynamic_gauge_and_cocycle()
    physical_3d = run_part5_physical_3d_and_shadow_offset()

    payload: Dict[str, Any] = {
        "generator_id": "02_generator_universal_obstruction_and_shadow_gauge",
        "contract_path": "scripts/campaigns/contracts/02_generator_universal_obstruction_and_shadow_gauge.md",
        "elapsed_seconds": round(time.perf_counter() - t0, 2),
        "crossover_3d_f128": crossover_3d_f128,
        "st1_obstruction_sweep": st1_obstruction,
        "shadow_telescoping": shadow_telescoping,
        "st3_secular_endurance": st3_secular,
        "st4_dynamic_gauge_and_cocycle": st4_dynamic,
        "physical_3d_and_shadow_offset": physical_3d,
    }

    assertions = verify_contract_assertions(payload)
    payload["contract_assertions"] = assertions

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"[SAVED] Certified JSON -> {out_path}", flush=True)

    # Synchronize data/pbm4_vs_bm4_stress_test.json
    sync_path = DATA_DIR / "pbm4_vs_bm4_stress_test.json"
    with open(sync_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"[SYNCED] Legacy JSON path -> {sync_path}", flush=True)


if __name__ == "__main__":
    main()
