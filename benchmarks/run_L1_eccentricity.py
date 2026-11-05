#!/usr/bin/env python3
r"""
Campaign I — Protocol L1: Extreme Keplerian Eccentricity Sweep (e -> 1^-).
Demonstrates Theorem 5.20 (Leaf/Apex Decoupling: k(t) == 0 and machine-floor accuracy
in Float64 and Float128 vs (1-e)^(-5/2) Cartesian symplectic divergence).
Outputs: data/campaign_I_L1_eccentricity.json
"""

from __future__ import annotations
import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, Any

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
from normal_cone_nbody_float64.benchmarks.registry import create_eccentric_binary_benchmark as create_ecc_f64
from normal_cone_nbody_float64.core.invariants import compute_cartesian_hamiltonian as ham_f64

from normal_cone_nbody_float128.system import NBodySystem as NBodySystemF128
from normal_cone_nbody_float128.integrators import (
    BlanesMoan4Integrator as BM4F128,
    create_pbm4_integrator as create_pbm4_f128,
)
from normal_cone_nbody_float128.benchmarks.registry import create_eccentric_binary_benchmark as create_ecc_f128
from normal_cone_nbody_float128.core.invariants import compute_cartesian_hamiltonian as ham_f128


def _cartesian_verlet_eccentric_binary_error(one_minus_e: float, n_steps: int = 200) -> float:
    e = 1.0 - float(one_minus_e)
    a, m0, m1, G = 1.0, 1.0, 1.0, 1.0
    st = create_ecc_f64(e=e, a=a, m0=m0, m1=m1, G=G)
    q = st.q.copy()
    p = st.p.copy()
    masses = st.masses.copy()
    E0 = float(st.E0)
    T_orb = 2.0 * math.pi * (a ** 1.5) / math.sqrt(G * (m0 + m1))
    dt = T_orb / float(n_steps)

    def _force(q_curr: np.ndarray) -> np.ndarray:
        r_vec = q_curr[1] - q_curr[0]
        r_norm = float(np.linalg.norm(r_vec))
        if r_norm < 1e-15:
            return np.zeros_like(q_curr)
        f01 = G * m0 * m1 * r_vec / (r_norm ** 3)
        return np.array([f01, -f01], dtype=np.float64)

    max_rel_err = 0.0
    f = _force(q)
    for _ in range(n_steps):
        p = p + 0.5 * dt * f
        q = q + dt * (p / masses[:, None])
        f = _force(q)
        p = p + 0.5 * dt * f
        E_curr = ham_f64(q, p, masses, G=G)
        rel_err = abs(E_curr - E0) / abs(E0)
        if rel_err > max_rel_err:
            max_rel_err = rel_err
    return float(max_rel_err)


def main() -> None:
    t0 = time.perf_counter()
    print("=== [Protocol L1] Extreme Eccentricity Sweep (e -> 1^-) ===", flush=True)
    grid = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
    n_points = len(grid)

    cart_errs, f64_bm4_errs, f64_pbm4_errs = [], [], []
    f128_bm4_errs, f128_pbm4_errs = [], []
    max_k_f64, max_k_f128 = [], []

    T_orb = 2.0 * math.pi / math.sqrt(2.0)
    dt_f64 = T_orb / 40.0
    dt_f128 = T_orb / 20.0

    for idx, ome in enumerate(grid, 1):
        t_pt0 = time.perf_counter()
        e_val = 1.0 - ome
        print(f"[{idx}/{n_points}] (1 - e) = {ome:.0e} (e = {e_val:.5f})...", end=" ", flush=True)

        # 1. Cartesian Verlet
        c_err = _cartesian_verlet_eccentric_binary_error(ome, n_steps=200)
        cart_errs.append(c_err)

        # 2. Float64 BM4 & PBM4
        st_f64 = create_ecc_f64(e=e_val, a=1.0, m0=1.0, m1=1.0, G=1.0)
        sys_b64 = NBodySystemF64(st_f64.with_updates(), G=1.0, integrator=BM4F64(G=1.0))
        sys_b64.evolve(t_span=T_orb, dt=dt_f64)
        k_f64 = sys_b64.history.get_dyadic_levels()
        max_k_f64.append(int(np.max(k_f64)) if len(k_f64) > 0 else 0)
        E_f64 = sys_b64.history.get_energies()
        err_b64 = float(np.max(np.abs(E_f64 - E_f64[0])) / abs(E_f64[0])) if len(E_f64) > 0 else 0.0
        f64_bm4_errs.append(max(err_b64, 1e-16))

        sys_p64 = NBodySystemF64(st_f64.with_updates(), G=1.0, integrator=create_pbm4_f64(G=1.0))
        sys_p64.evolve(t_span=T_orb, dt=dt_f64)
        E_p64 = sys_p64.history.get_energies()
        err_p64 = float(np.max(np.abs(E_p64 - E_p64[0])) / abs(E_p64[0])) if len(E_p64) > 0 else 0.0
        f64_pbm4_errs.append(max(err_p64, 1e-16))

        # 3. Float128 (35-digit mpmath) BM4 & PBM4 (5 steps)
        mp.mp.dps = 35
        e_mpf = mp.mpf(1) - mp.mpf(str(ome))
        st_f128 = create_ecc_f128(e=e_mpf, a=1.0, m0=1.0, m1=1.0, G=1.0, is_f128=True)
        sys_b128 = NBodySystemF128(st_f128.with_updates(), G=1.0, integrator=BM4F128(G=1.0))
        for _ in range(5):
            sys_b128.step(dt=dt_f128)
        k_128 = sys_b128.history.get_dyadic_levels()
        max_k_f128.append(int(np.max(k_128)) if len(k_128) > 0 else 0)
        E_b128 = ham_f128(sys_b128.state.q, sys_b128.state.p, sys_b128.state.masses, G=1.0)
        err_b128 = float(abs(E_b128 - st_f128.E0) / abs(st_f128.E0))
        f128_bm4_errs.append(max(err_b128, 1e-35))

        sys_p128 = NBodySystemF128(st_f128.with_updates(), G=1.0, integrator=create_pbm4_f128(G=1.0))
        for _ in range(5):
            sys_p128.step(dt=dt_f128)
        E_p128 = ham_f128(sys_p128.state.q, sys_p128.state.p, sys_p128.state.masses, G=1.0)
        err_p128 = float(abs(E_p128 - st_f128.E0) / abs(st_f128.E0))
        f128_pbm4_errs.append(max(err_p128, 1e-35))

        dt_pt = time.perf_counter() - t_pt0
        print(f"CartVerlet={c_err:.2e} | f64_PBM4={err_p64:.2e} | f128_PBM4={err_p128:.2e} | k_max={max_k_f64[-1]} ({dt_pt:.2f}s)", flush=True)

    res = {
        "one_minus_e": grid,
        "cartesian_verlet_rel_dE": cart_errs,
        "float64_bm4_rel_dE": f64_bm4_errs,
        "float64_pbm4_rel_dE": f64_pbm4_errs,
        "float128_bm4_rel_dE": f128_bm4_errs,
        "float128_pbm4_rel_dE": f128_pbm4_errs,
        "max_dyadic_level_f64": max_k_f64,
        "max_dyadic_level_f128": max_k_f128,
        "elapsed_seconds": time.perf_counter() - t0,
    }

    out_dir = WORKSPACE_ROOT / "results" / "campaigns"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "campaign_I_L1_eccentricity.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    print(f"=== [Protocol L1 COMPLETE] Saved to {out_file} in {res['elapsed_seconds']:.2f}s ===", flush=True)


if __name__ == "__main__":
    main()
