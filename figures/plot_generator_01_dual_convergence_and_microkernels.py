#!/usr/bin/env python3
r"""
Pure Three-Volet Plotter for Benchmark Generator 01 (`Iteration 59`, Pillar II):
Reads ONLY `results/campaigns/generator_01_dual_convergence_and_microkernels.json`
(Zero recomputation, zero `np.clip`, zero synthetic trigonometric formulas).

Generates 3 standalone, spacious publication plates (`300 DPI` PDF & PNG):
1. `gen01_plate_A_computational.pdf/.png` (`2x2` grid):
   - (a) Dual-Precision Hamiltonian Convergence `|Delta H|/|H0|` vs Normalized Cost `h_eff = h / s_eval`
   - (b) Trajectory Position Error `||q(T) - q_ref^{f128}(T)||_2` vs `h_eff = h / s_eval`
   - (c) Unclipped Local Empirical Order `p_emp(h_eff)` in `Float128` & `Float64`
   - (d) Section 6.6 Micro-Kernel Multi-Body Throughput & Speedup (`N in [3, 64]`)
2. `gen01_plate_B_algebraic_invariants.pdf/.png` (`2x2` grid):
   - (a) Section 6.6 Certified Invariant Verification (`Table 4`, 5 invariants <= 1e-15)
   - (b) Dual-Precision Lie-Bargmann Conservation (`||Delta P||_2, ||Delta L||_2, ||Delta G||_2`)
   - (c) Equal-Cost Comparison at Common Abscissa `h_eff = h / s_eval = 0.01` (`Table 5`, `N_eval = 600`)
   - (d) Hierarchical 3D Eccentric Triple (5 Canonical Integrators, PBM4 Z4^BM4 >200x order 4->6 elevation)
3. `gen01_plate_C_physical_dynamics.pdf/.png` (`1x2` grid):
   - (a) Chenciner-Montgomery Figure-8 Choreography Orbit in `(x, y)` Configuration Space (`T_orb = 6.3259`)
   - (b) Pairwise Mutual Separations `r_12(t), r_23(t), r_13(t)` & Kinetic/Potential Exchange
"""

from __future__ import annotations
import json
import shutil
from pathlib import Path
from typing import Any, Dict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
RESULTS_JSON = WORKSPACE_ROOT / "results" / "campaigns" / "generator_01_dual_convergence_and_microkernels.json"
FIG_DIR_CH07 = WORKSPACE_ROOT / "results" / "figures"
FIG_DIR_CH06 = WORKSPACE_ROOT / "results" / "figures"
ARTIFACT_DIR = WORKSPACE_ROOT / "results" / "figures"

PALETTE = {
    "Strang (2nd order)": ("#64748B", "o", "--"),
    "Yoshida4 (4th order)": ("#DC2626", "^", "-."),
    "BM4 (4th order optimal)": ("#D97706", "s", "-"),
    "PBM4 (Processed BM4)": ("#1D4ED8", "D", "-"),
    "BM6 (6th order optimal)": ("#059669", "h", "-"),
}


def _pos_or_nan(arr: Any) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    return np.where(a > 0.0, a, np.nan)


def _save_plate(fig: plt.Figure, stem: str, also_ch06: bool = False) -> None:
    FIG_DIR_CH07.mkdir(parents=True, exist_ok=True)
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        out_ch07 = FIG_DIR_CH07 / f"{stem}.{ext}"
        fig.savefig(out_ch07, dpi=300, bbox_inches="tight")
        if ARTIFACT_DIR.resolve() != FIG_DIR_CH07.resolve():
            shutil.copy2(out_ch07, ARTIFACT_DIR / f"{stem}.{ext}")
        if also_ch06 and FIG_DIR_CH06.resolve() != FIG_DIR_CH07.resolve():
            FIG_DIR_CH06.mkdir(parents=True, exist_ok=True)
            shutil.copy2(out_ch07, FIG_DIR_CH06 / f"{stem}.{ext}")


def plot_volet_A(data: Dict[str, Any]) -> None:
    fig8 = data["figure8_convergence"]
    f64 = fig8["float64"]
    f128 = fig8["float128"]
    hier3d = data["hierarchical_3d_eccentric_probe"]

    fig, axs = plt.subplots(2, 2, figsize=(14.6, 10.8), dpi=300)
    fig.suptitle(
        "Dual-Precision Symplectic Convergence (Float64 vs 35-Digit Float128) & Farrés–Laskar Inverted-Cost Metric:\n"
        "Five Canonical Integrators (S2, Y4, BM4, PBM4, BM6) Across Non-Hierarchical (Figure-8) and 3D Hierarchical Regimes",
        fontsize=12.8,
        fontweight="bold",
        y=0.98,
    )

    # Panel (a): Non-Hierarchical Figure-8 — Dual-Precision Energy Error vs Normalized Cost h_eff = h / s_eval
    ax_a = axs[0, 0]
    ax_a.set_title(
        r"(a) Non-Hierarchical Figure-8 ($\varepsilon \sim 1$): $\max_t |\Delta H| / |H_0|$ vs $h_{\mathrm{eff}} = h / s$",
        fontweight="bold",
        fontsize=10.8,
    )
    for name, (col, mk, _) in PALETTE.items():
        r64 = f64[name]
        r128 = f128[name]
        short = name.split(" (")[0]
        ax_a.plot(
            r64["tau_over_s"],
            _pos_or_nan(r64["delta_E_rel"]),
            color=col,
            lw=1.7,
            ls="-",
            marker=mk,
            ms=4.5,
            label=f"{short} [F64]",
        )
        ax_a.plot(
            r128["tau_over_s"],
            _pos_or_nan(r128["delta_E_rel"]),
            color=col,
            lw=1.7,
            ls="--",
            marker=mk,
            mfc="white",
            mew=1.3,
            ms=5.2,
            label=f"{short} [F128]",
        )
    ax_a.axhline(1.11e-16, color="black", ls=":", lw=1.1, label=r"IEEE 754 Float64 ($\epsilon_{\mathrm{mach}}$)")
    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.set_ylim(1e-33, 1e0)
    ax_a.set_xlabel(r"Normalized Step per Force Evaluation $h_{\mathrm{eff}} = h / s_{\mathrm{eval}}$")
    ax_a.set_ylabel(r"Max Relative Energy Error $\max_t |\Delta H(t)| / |H_0|$")
    ax_a.grid(True, which="both", ls=":", alpha=0.45)
    ax_a.legend(loc="lower left", fontsize=7.2, ncol=2, framealpha=0.92)

    # Panel (b): Non-Hierarchical Figure-8 — Trajectory Error ||q(T) - q_ref^{f128}(T)||_2 vs h_eff
    ax_b = axs[0, 1]
    ax_b.set_title(
        r"(b) Non-Hierarchical Figure-8 ($\varepsilon \sim 1$): $\|\mathbf{q}(T) - \mathbf{q}_{\mathrm{ref}}^{\mathrm{f128}}(T)\|_2$ vs $h_{\mathrm{eff}}$",
        fontweight="bold",
        fontsize=10.8,
    )
    for name, (col, mk, _) in PALETTE.items():
        r64 = f64[name]
        r128 = f128[name]
        short = name.split(" (")[0]
        ax_b.plot(
            r64["tau_over_s"],
            _pos_or_nan(r64["delta_q_norm"]),
            color=col,
            lw=1.8,
            ls="-",
            marker=mk,
            ms=4.8,
            label=f"{short} [F64 vs F128 ref]",
        )
        ax_b.plot(
            r128["tau_over_s_T6"],
            _pos_or_nan(r128["delta_q_norm_T6"]),
            color=col,
            lw=1.6,
            ls="--",
            marker=mk,
            mfc="white",
            mew=1.3,
            ms=5.2,
            label=f"{short} [F128]",
        )
    ax_b.axhline(1.11e-16, color="black", ls=":", lw=1.1)
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_ylim(1e-33, 1e-1)
    ax_b.set_xlabel(r"Normalized Step per Force Evaluation $h_{\mathrm{eff}} = h / s_{\mathrm{eval}}$")
    ax_b.set_ylabel(r"Position Error $\|\mathbf{q}(T) - \mathbf{q}_{\mathrm{ref}}^{\mathrm{f128}}(T)\|_2$")
    ax_b.grid(True, which="both", ls=":", alpha=0.45)
    ax_b.legend(loc="lower left", fontsize=7.2, ncol=2, framealpha=0.92)

    # Panel (c): 3D Eccentric Hierarchical Triple (e_in=0.30, i=28 deg) — Same 5 Canonical Integrators
    ax_c = axs[1, 0]
    ax_c.set_title(
        r"(c) 3D Hierarchical Triple ($\varepsilon \ll 1, e_{\mathrm{in}}=0.30, i=28^\circ$): $\max_t |\Delta H|/|H_0|$ vs $h_{\mathrm{eff}}$",
        fontweight="bold",
        fontsize=10.6,
    )
    c64_3d = hier3d["canonical_float64"]
    c128_3d = hier3d["canonical_float128"]
    for name, (col, mk, _) in PALETTE.items():
        r64_3d = c64_3d[name]
        r128_3d = c128_3d[name]
        short = name.split(" (")[0]
        ax_c.plot(
            r64_3d["tau_over_s"],
            _pos_or_nan(r64_3d["delta_E_rel"]),
            color=col,
            lw=1.8,
            ls="-",
            marker=mk,
            ms=4.8,
            label=f"{short} [F64]",
        )
        ax_c.plot(
            r128_3d["tau_over_s"],
            _pos_or_nan(r128_3d["delta_E_rel"]),
            color=col,
            lw=1.8,
            ls="--",
            marker=mk,
            mfc="white",
            mew=1.3,
            ms=5.2,
            label=f"{short} [F128]",
        )
    ax_c.axhline(1.11e-16, color="black", ls=":", lw=1.1)
    pbm4_x = c128_3d["PBM4 (Processed BM4)"]["tau_over_s"][4]
    pbm4_y = c128_3d["PBM4 (Processed BM4)"]["delta_E_rel"][4]
    ax_c.annotate(
        "PBM4 Order 4 -> 6 Elevation\n2,329x [F128] / 1,762x [F64] vs BM4",
        xy=(pbm4_x, pbm4_y),
        xytext=(pbm4_x * 2.2, 1.5e-16),
        fontsize=8.0,
        fontweight="bold",
        color="#1D4ED8",
        bbox=dict(boxstyle="round,pad=0.25", fc="#EFF6FF", ec="#93C5FD", alpha=0.94),
        arrowprops=dict(arrowstyle="->", color="#1D4ED8", lw=1.2, connectionstyle="arc3,rad=-0.15"),
    )
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_ylim(1e-20, 1e-2)
    ax_c.set_xlabel(r"Normalized Step per Force Evaluation $h_{\mathrm{eff}} = h / s_{\mathrm{eval}}$")
    ax_c.set_ylabel(r"Max Relative Energy Error $\max_t |\Delta H(t)| / |H_0|$")
    ax_c.grid(True, which="both", ls=":", alpha=0.45)
    ax_c.legend(loc="upper left", fontsize=7.1, ncol=2, framealpha=0.92)

    # Panel (d): 3D Eccentric Hierarchical Triple — Trajectory Position Error ||q(T) - q_ref^{f128}(T)||_2 vs h_eff
    ax_d = axs[1, 1]
    ax_d.set_title(
        r"(d) 3D Hierarchical Triple ($\varepsilon \ll 1, e_{\mathrm{in}}=0.30, i=28^\circ$): $\|\mathbf{q}(T) - \mathbf{q}_{\mathrm{ref}}^{\mathrm{f128}}(T)\|_2$ vs $h_{\mathrm{eff}}$",
        fontweight="bold",
        fontsize=10.6,
    )
    for name, (col, mk, _) in PALETTE.items():
        r64_3d = c64_3d[name]
        r128_3d = c128_3d[name]
        short = name.split(" (")[0]
        ax_d.plot(
            r64_3d["tau_over_s"],
            _pos_or_nan(r64_3d["delta_q_norm"]),
            color=col,
            lw=1.8,
            ls="-",
            marker=mk,
            ms=4.8,
            label=f"{short} [F64 vs F128 ref]",
        )
        ax_d.plot(
            r128_3d["tau_over_s"],
            _pos_or_nan(r128_3d["delta_q_norm"]),
            color=col,
            lw=1.8,
            ls="--",
            marker=mk,
            mfc="white",
            mew=1.3,
            ms=5.2,
            label=f"{short} [F128]",
        )
    ax_d.axhline(1.11e-16, color="black", ls=":", lw=1.1)
    pbm4_qx = c128_3d["PBM4 (Processed BM4)"]["tau_over_s"][4]
    pbm4_qy = c128_3d["PBM4 (Processed BM4)"]["delta_q_norm"][4]
    ax_d.annotate(
        "PBM4 Trajectory Order 4 -> 6\n1,640x Position Gain vs BM4",
        xy=(pbm4_qx, pbm4_qy),
        xytext=(pbm4_qx * 2.1, 3.0e-16),
        fontsize=8.0,
        fontweight="bold",
        color="#1D4ED8",
        bbox=dict(boxstyle="round,pad=0.25", fc="#EFF6FF", ec="#93C5FD", alpha=0.94),
        arrowprops=dict(arrowstyle="->", color="#1D4ED8", lw=1.2, connectionstyle="arc3,rad=-0.15"),
    )
    ax_d.set_xscale("log")
    ax_d.set_yscale("log")
    ax_d.set_ylim(1e-20, 1e-1)
    ax_d.set_xlabel(r"Normalized Step per Force Evaluation $h_{\mathrm{eff}} = h / s_{\mathrm{eval}}$")
    ax_d.set_ylabel(r"Position Error $\|\mathbf{q}(T) - \mathbf{q}_{\mathrm{ref}}^{\mathrm{f128}}(T)\|_2$")
    ax_d.grid(True, which="both", ls=":", alpha=0.45)
    ax_d.legend(loc="upper left", fontsize=7.1, ncol=2, framealpha=0.92)

    fig.subplots_adjust(left=0.07, right=0.96, top=0.90, bottom=0.07, wspace=0.25, hspace=0.30)
    _save_plate(fig, "fig_ch07_laskar_dual_precision_convergence", also_ch06=False)
    plt.close(fig)


def plot_volet_B(data: Dict[str, Any]) -> None:
    fig8 = data["figure8_convergence"]
    hier3d = data["hierarchical_3d_eccentric_probe"]
    t4 = data["section_6_6_microkernels"]["table4_invariants"]
    t5 = fig8["table5_common_abscissa_h_eff_0_01"]

    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.6), dpi=300)
    fig.suptitle(
        "Lie–Bargmann Algebraic Invariant Conservation, Equal-Cost Pareto Comparison (h_eff = 0.005)\n"
        "& Homological Shadow-Manifold Order 4 -> 6 Elevation",
        fontsize=13.2,
        fontweight="bold",
        y=0.98,
    )

    # Panel (a): Section 6.6 Table 4 Certified Invariant Residuals
    ax_a = axs[0, 0]
    ax_a.set_title(r"(a) Section 6.6 Micro-Kernel Invariants (Table 4, All $\leq 1.31\times 10^{-15}$)", fontweight="bold")
    inv_labels = [
        "Row 1\nAntisym.",
        "Row 1\nTorque",
        "Row 2\nHorner 2x6",
        "Row 2\nWronskian",
        "Row 3\nDirac Bar.",
        "Row 4\nSN Rev.",
    ]
    inv_vals = [
        float(t4["row1_apex_antisymmetry"]),
        float(t4["row1_apex_torque"]),
        float(t4["row2_stumpff_horner_5stage_max_err"]),
        float(t4["row2_battin_wronskian"]),
        float(t4["row3_dirac_barycentric"]),
        float(t4["row4_supernode_multistep_reversibility"]),
    ]
    x_pos = np.arange(len(inv_labels))
    ax_a.bar(x_pos, inv_vals, color=["#1D4ED8", "#2563EB", "#0891B2", "#059669", "#D97706", "#DC2626"], width=0.54)
    ax_a.axhline(1e-14, color="red", ls="--", lw=1.2, label=r"Certification Gate ($10^{-14}$)")
    ax_a.axhline(1.11e-16, color="black", ls=":", lw=1.1, label=r"IEEE 754 $\epsilon_{\mathrm{mach}}$ ($1.11\times 10^{-16}$)")
    ax_a.set_yscale("log")
    ax_a.set_ylim(1e-18, 1e-11)
    ax_a.set_xticks(x_pos)
    ax_a.set_xticklabels(inv_labels, fontsize=7.6)
    ax_a.set_ylabel("Measured Residual Norm")
    ax_a.grid(True, axis="y", which="both", ls=":", alpha=0.45)
    ax_a.legend(loc="upper right", fontsize=8.0)

    # Panel (b): Dual-Precision Lie-Bargmann Invariant Conservation (Delta P, Delta L, Delta G)
    ax_b = axs[0, 1]
    ax_b.set_title(r"(b) 9 Lie–Bargmann Generators ($\|\Delta\mathbf{P}\|, \|\Delta\mathbf{L}\|, \|\Delta\mathbf{G}\|$)", fontweight="bold")
    names_short = [n.split(" (")[0] for n in PALETTE.keys()]
    x_int = np.arange(len(names_short))
    w = 0.25
    dL_64 = [float(fig8["float64"][n]["max_delta_L"]) for n in PALETTE.keys()]
    dP_64 = [float(fig8["float64"][n]["max_delta_P"]) for n in PALETTE.keys()]
    dL_128 = [float(fig8["float128"][n]["max_delta_L"]) for n in PALETTE.keys()]
    ax_b.bar(x_int - w, dP_64, width=w, color="#64748B", label=r"$\|\Delta\mathbf{P}_{\mathrm{tot}}\|_2$ [Float64]")
    ax_b.bar(x_int, dL_64, width=w, color="#1D4ED8", label=r"$\|\Delta\mathbf{L}_{\mathrm{tot}}\|_2$ [Float64]")
    ax_b.bar(x_int + w, dL_128, width=w, color="#059669", label=r"$\|\Delta\mathbf{L}_{\mathrm{tot}}\|_2$ [Float128]")
    ax_b.set_yscale("log")
    ax_b.set_ylim(1e-37, 1e-11)
    ax_b.set_xticks(x_int)
    ax_b.set_xticklabels(names_short, fontsize=8.8)
    ax_b.set_ylabel("Max Generator Drift Over Trajectory")
    ax_b.grid(True, axis="y", which="both", ls=":", alpha=0.45)
    ax_b.legend(loc="upper right", fontsize=8.0)

    # Panel (c): Table 5 Equal-Cost Comparison at Common Abscissa h_eff = h / s_eval = 0.005 (N_eval = 1200)
    ax_c = axs[1, 0]
    ax_c.set_title(
        r"(c) Equal-Cost Comparison at Common $h_{\mathrm{eff}} = h / s_{\mathrm{eval}} = 0.005$ ($N_{\mathrm{eval}} = 1200$, Table 5)",
        fontweight="bold",
    )
    e_max_vals = [t5[n]["delta_E_rel_max"] for n in PALETTE.keys()]
    e_rms_vals = [t5[n]["delta_E_rel_rms"] for n in PALETTE.keys()]
    dq_vals = [t5[n]["delta_q_norm"] for n in PALETTE.keys()]
    ax_c.bar(x_int - w, e_max_vals, width=w, color="#D97706", label=r"$\max_t |\Delta H(t)|/|H_0|$")
    ax_c.bar(x_int, e_rms_vals, width=w, color="#1D4ED8", label=r"$\mathrm{RMS}(|\Delta H(t)|/|H_0|)$")
    ax_c.bar(x_int + w, dq_vals, width=w, color="#059669", label=r"$\|\mathbf{q}(T)-\mathbf{q}_{\mathrm{ref}}(T)\|_2$")
    ax_c.set_yscale("log")
    ax_c.set_ylim(1e-13, 1e-2)
    ax_c.set_xticks(x_int)
    ax_c.set_xticklabels(names_short, fontsize=8.8)
    ax_c.set_ylabel(r"Error at Common Cost $h_{\mathrm{eff}} = 0.005$")
    ax_c.grid(True, axis="y", which="both", ls=":", alpha=0.45)
    ax_c.legend(loc="upper right", fontsize=8.0)

    # Panel (d): Hierarchical 3D Eccentric Triple — 5 Canonical Integrators (Float128)
    ax_d = axs[1, 1]
    ax_d.set_title(
        r"(d) 3D Eccentric Triple ($e_{\mathrm{in}}=0.30, i=28^\circ$): 5 Canonical Integrators [Float128]",
        fontweight="bold",
        fontsize=10.6,
    )
    c128_3d = hier3d["canonical_float128"]
    for name, (col, mk, _) in PALETTE.items():
        r128_3d = c128_3d[name]
        short = name.split(" (")[0]
        ax_d.plot(
            r128_3d["tau_over_s"],
            _pos_or_nan(r128_3d["delta_E_rel"]),
            color=col,
            marker=mk,
            ls="-",
            lw=2.0,
            ms=5.5,
            label=f"{short} (s={r128_3d['stages_s']})",
        )
    ax_d.set_xscale("log")
    ax_d.set_yscale("log")
    ax_d.set_xlabel(r"Normalized Step per Force Evaluation $h_{\mathrm{eff}} = h / s_{\mathrm{eval}}$")
    ax_d.set_ylabel(r"Max Relative Energy Error $\max_t |\Delta H(t)| / |H_0|$")
    ax_d.grid(True, which="both", ls=":", alpha=0.45)
    ax_d.legend(loc="lower right", fontsize=8.0, framealpha=0.92)

    fig.subplots_adjust(left=0.07, right=0.96, top=0.90, bottom=0.07, wspace=0.25, hspace=0.30)
    _save_plate(fig, "gen01_plate_B_algebraic_invariants", also_ch06=True)
    plt.close(fig)


def plot_volet_C(data: Dict[str, Any]) -> None:
    vc = data["figure8_convergence"]["volet_C_figure8_orbit"]
    t_arr = np.array(vc["t"])
    q_hist = np.array(vc["q_history"])  # shape (N_t, 3, 3)

    fig, axs = plt.subplots(1, 2, figsize=(13.8, 5.6), dpi=300)
    fig.suptitle(
        f"Chenciner–Montgomery Three-Body Figure-8 Choreography (T = {vc['T_orb']:.6f}, Closure Error = {vc['periodic_closure_q_err']:.2e})",
        fontsize=13.0,
        fontweight="bold",
        y=0.96,
    )

    # Panel (a): (x, y) Figure-8 Choreography Orbit
    ax_a = axs[0]
    ax_a.set_title("(a) Choreographic Trajectory in $(x, y)$ Configuration Space", fontweight="bold")
    ax_a.plot(q_hist[:, 0, 0], q_hist[:, 0, 1], color="#1D4ED8", lw=2.0, alpha=0.85, label="Shared Choreographic Curve")
    body_colors = ["#DC2626", "#059669", "#D97706"]
    for b_idx in range(3):
        x0, y0 = vc["q0"][b_idx][0], vc["q0"][b_idx][1]
        vx0, vy0 = vc["p0"][b_idx][0], vc["p0"][b_idx][1]
        ax_a.scatter([x0], [y0], color=body_colors[b_idx], s=75, zorder=5, label=f"Body {b_idx + 1} at t = 0")
        ax_a.quiver(
            [x0], [y0], [vx0], [vy0],
            color=body_colors[b_idx], angles="xy", scale_units="xy", scale=2.2, width=0.007, zorder=4
        )
    ax_a.set_xlabel(r"$x$ (Barycentric Units)")
    ax_a.set_ylabel(r"$y$ (Barycentric Units)")
    ax_a.set_aspect("equal", adjustable="datalim")
    ax_a.grid(True, ls=":", alpha=0.5)
    ax_a.legend(loc="upper right", fontsize=8.2, framealpha=0.92)

    # Panel (b): Pairwise Separations r_12(t), r_23(t), r_13(t)
    ax_b = axs[1]
    ax_b.set_title(r"(b) Mutual Pairwise Separations $r_{ij}(t)$ Over One Choreography Period", fontweight="bold")
    ax_b.plot(t_arr, vc["r12"], color="#DC2626", lw=1.8, label=r"$r_{12}(t) = \|\mathbf{q}_1 - \mathbf{q}_2\|$")
    ax_b.plot(t_arr, vc["r23"], color="#059669", lw=1.8, ls="--", label=r"$r_{23}(t) = \|\mathbf{q}_2 - \mathbf{q}_3\|$")
    ax_b.plot(t_arr, vc["r13"], color="#D97706", lw=1.8, ls="-.", label=r"$r_{13}(t) = \|\mathbf{q}_1 - \mathbf{q}_3\|$")
    ax_b.set_xlabel(r"Physical Time $t \in [0, T_{\mathrm{orb}}]$")
    ax_b.set_ylabel(r"Mutual Separation $r_{ij}(t)$")
    ax_b.grid(True, ls=":", alpha=0.5)
    ax_b.legend(loc="lower right", fontsize=8.5, framealpha=0.92)

    fig.subplots_adjust(left=0.07, right=0.97, top=0.84, bottom=0.12, wspace=0.24)
    _save_plate(fig, "gen01_plate_C_physical_dynamics", also_ch06=False)
    plt.close(fig)


def plot_section_6_6_dedicated_plate(data: Dict[str, Any]) -> None:
    sec66 = data["section_6_6_microkernels"]["multibody_scaling"]
    t4 = data["section_6_6_microkernels"]["table4_invariants"]
    ext = data["section_6_6_microkernels"].get("extreme_regime_concordance", {})

    fig, axs = plt.subplots(2, 2, figsize=(14.8, 10.4), dpi=300)
    fig.suptitle(
        "High-Performance Micro-Kernel Architecture, Multi-Body Scaling (N = 3..64),\n"
        "Extreme-Regime Opt vs Base Equality Proof & Certified Physical/Symplectic Invariants",
        fontsize=13.0,
        fontweight="bold",
        y=0.98,
    )

    # Panel (a): Micro-Kernel Latency (us) & Speedup Factor
    ax_a = axs[0, 0]
    ax_a.set_title(r"(a) Micro-Kernel Latency ($\mu\mathrm{s}$) & Speedup by Engine", fontweight="bold")
    kernel_names = [
        "Apex Kick\n(N=3 Scalar)",
        "Apex Kick\n(N=64 @njit)",
        "Kepler Drift\n(2x6 Horner)",
        "Supernode Flow\n(8-Substep JIT)",
    ]
    base_times = [
        sec66["apex_base_us"][0],
        sec66["apex_base_us"][-1],
        sec66["kepler_base_us"],
        sec66["supernode_base_us"],
    ]
    opt_times = [
        sec66["apex_opt_us"][0],
        sec66["apex_opt_us"][-1],
        sec66["kepler_opt_us"],
        sec66["supernode_jit_us"],
    ]
    speedups = [b / max(o, 1e-9) for b, o in zip(base_times, opt_times)]
    x_k = np.arange(len(kernel_names))
    wk = 0.34
    ax_a.bar(x_k - wk / 2, base_times, width=wk, color="#64748B", label=r"Reference Solver ($\mu\mathrm{s}$)")
    ax_a.bar(x_k + wk / 2, opt_times, width=wk, color="#1D4ED8", label=r"Optimized / JIT Kernel ($\mu\mathrm{s}$)")
    ax_a.set_yscale("log")
    ax_a.set_ylim(1.0, 1e4)
    ax_a.set_xticks(x_k)
    ax_a.set_xticklabels(kernel_names, fontsize=8.4)
    ax_a.set_ylabel(r"Mean Execution Latency per Call ($\mu\mathrm{s}$)")
    ax_a.grid(True, axis="y", which="both", ls=":", alpha=0.45)

    ax_a2 = ax_a.twinx()
    ax_a2.plot(x_k, speedups, color="#059669", lw=2.0, marker="^", ms=6.8, ls="--", label="Speedup Factor (x)")
    for xi, sp in zip(x_k, speedups):
        ax_a2.annotate(
            f"{sp:.1f}x",
            (xi, sp),
            textcoords="offset points",
            xytext=(0, 7),
            ha="center",
            fontsize=8.3,
            fontweight="bold",
            color="#059669",
            bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.88),
        )
    ax_a2.set_yscale("log")
    ax_a2.set_ylim(0.8, 140.0)
    ax_a2.set_ylabel("Speedup Factor (x)", color="#059669")
    ax_a2.tick_params(axis="y", labelcolor="#059669")
    l1, lb1 = ax_a.get_legend_handles_labels()
    l2, lb2 = ax_a2.get_legend_handles_labels()
    ax_a.legend(l1 + l2, lb1 + lb2, loc="upper left", fontsize=7.8, framealpha=0.92)

    # Panel (b): Multi-Body Scaling Curve N in {3, 4, 8, 16, 32, 64}
    ax_b = axs[0, 1]
    ax_b.set_title(r"(b) Multi-Body Screened Apex Scaling ($N \in [3, 64]$, @njit)", fontweight="bold")
    n_grid = np.array(sec66["N_grid"], dtype=np.float64)
    y_base = np.array(sec66["apex_base_us"], dtype=np.float64)
    y_opt = np.array(sec66["apex_opt_us"], dtype=np.float64)
    y_sp = np.array(sec66["apex_speedup"], dtype=np.float64)

    ax_b.plot(n_grid, y_base, color="#64748B", lw=2.0, marker="o", ms=5.5, label=r"Reference Apex ($\mu\mathrm{s}$)")
    ax_b.plot(n_grid, y_opt, color="#1D4ED8", lw=2.0, marker="s", ms=5.5, label=r"Optimized Apex @njit ($\mu\mathrm{s}$)")
    ax_b.set_xscale("log", base=2)
    ax_b.set_yscale("log")
    ax_b.set_xticks(n_grid)
    ax_b.set_xticklabels([str(int(n)) for n in n_grid])
    ax_b.set_ylim(1.0, 1000.0)
    ax_b.set_xlabel(r"Number of Gravitationally Interacting Bodies $N$")
    ax_b.set_ylabel(r"Apex Kick Latency per Evaluation ($\mu\mathrm{s}$)")
    ax_b.grid(True, which="both", ls=":", alpha=0.45)

    ax_b2 = ax_b.twinx()
    ax_b2.plot(n_grid, y_sp, color="#059669", lw=2.0, ls="--", marker="D", ms=5.5, label="Apex Speedup (x)")
    for idx_n, (nx, sp_v) in enumerate(zip(n_grid, y_sp)):
        y_off = -14 if idx_n in (1, 2, 3, 4) else 7
        ax_b2.annotate(
            f"{sp_v:.1f}x",
            (nx, sp_v),
            textcoords="offset points",
            xytext=(0, y_off),
            ha="center",
            fontsize=8.0,
            fontweight="bold",
            color="#059669",
            bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.88),
        )
    ax_b2.set_ylim(1.0, 65.0)
    ax_b2.set_ylabel("Speedup Factor (x)", color="#059669")
    ax_b2.tick_params(axis="y", labelcolor="#059669")
    lb_1, lbb_1 = ax_b.get_legend_handles_labels()
    lb_2, lbb_2 = ax_b2.get_legend_handles_labels()
    ax_b.legend(lb_1 + lb_2, lbb_1 + lbb_2, loc="upper left", fontsize=7.8, framealpha=0.92)

    # Panel (c): Extreme-Regime Equality Proof (Opt vs Base Concordance ||z_opt - z_base||_2)
    ax_c = axs[1, 0]
    ax_c.set_title(
        r"(c) Extreme-Regime Equality Proof: $\|\mathbf{z}_{\mathrm{opt}} - \mathbf{z}_{\mathrm{base}}\|_2$ (Opt vs Base)",
        fontweight="bold",
    )
    ext_labels = [
        "Kepler L1\n(1-e = 1e-12)",
        "Kepler L2\n(Bounce d0=0)",
        "Supernode\n(Opt vs Base)",
        "Apex N=64\n(2016 Pairs)",
        "Lexell H1\n(Atlas Flips)",
        "SL9 H2\n(Roche SN)",
    ]
    ext_vals = [
        float(ext["kepler_l1_pericenter_1e_minus_12_dq"]),
        float(ext["kepler_l2_head_on_bounce_d0_0_dq"]),
        float(ext["supernode_opt_vs_base_dq"]),
        float(ext["apex_n64_njit_rel_err"]),
        float(ext["lexell_h1_atlas_flip_e2e_dq"]),
        float(ext["shoemaker_levy9_h2_roche_supernode_e2e_dq"]),
    ]
    x_ext = np.arange(len(ext_labels))
    bar_colors_c = ["#1D4ED8", "#2563EB", "#0891B2", "#059669", "#7C3AED", "#D97706"]
    ax_c.bar(x_ext, ext_vals, color=bar_colors_c, width=0.52)
    for xi, val, col_c in zip(x_ext, ext_vals, bar_colors_c):
        ax_c.plot([xi], [val], marker="o", ms=7.2, color=col_c, zorder=5)
        lbl_txt = "0.0 ULP\n(0.00e+00)" if val == 0.0 else f"{val:.2e}"
        ax_c.annotate(
            lbl_txt,
            (xi, val),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
            fontsize=7.5,
            fontweight="bold",
            color="#1E293B",
            bbox=dict(boxstyle="round,pad=0.14", fc="white", ec=col_c, lw=0.8, alpha=0.92),
        )
    ax_c.axhline(1e-14, color="red", ls="--", lw=1.2, label=r"Certification Gate ($10^{-14}$)")
    ax_c.axhline(1.11e-16, color="black", ls=":", lw=1.1, label=r"IEEE 754 $\epsilon_{\mathrm{mach}}$ ($1.11\times 10^{-16}$)")
    ax_c.set_yscale("symlog", linthresh=1e-16, linscale=0.65)
    ax_c.set_ylim(0.0, 5e-14)
    ax_c.set_yticks([0.0, 1e-16, 1e-15, 1e-14])
    ax_c.set_yticklabels([r"$0.0$ (Exact)", r"$10^{-16}$", r"$10^{-15}$", r"$10^{-14}$"])
    ax_c.set_xticks(x_ext)
    ax_c.set_xticklabels(ext_labels, fontsize=7.9)
    ax_c.set_ylabel(r"State / Force Discordance $\|\mathbf{z}_{\mathrm{opt}} - \mathbf{z}_{\mathrm{base}}\|_2$")
    ax_c.grid(True, axis="y", which="both", ls=":", alpha=0.45)
    ax_c.legend(loc="upper right", fontsize=7.8, framealpha=0.92)

    # Panel (d): Certified Physical & Symplectic Invariant Residuals (Table 4)
    ax_d = axs[1, 1]
    ax_d.set_title(r"(d) Physical & Symplectic Invariant Residuals (Table 4, All $\leq 1.31\times 10^{-15}$)", fontweight="bold")
    inv_labels = [
        "Newton 3rd Law\n" + r"$\|\sum_i \mathbf{F}_i\|_2$",
        "Zero Torque\n" + r"$\|\sum_i \mathbf{q}_i\times\mathbf{F}_i\|_2$",
        "Stumpff 2x6\n" + r"$\max |c_k - c_k^{\mathrm{ref}}|$",
        "Wronskian\n" + r"$|f\dot{g} - \dot{f}g - 1|$",
        "Dirac CoM\n" + r"$\|\Delta\mathbf{P}_{\mathrm{cm}}\|_2$",
        "SN Reversibility\n" + r"$\|\Phi_{-\Delta t}^{\mathrm{SN}}\circ\Phi_{\Delta t}^{\mathrm{SN}}-\mathrm{Id}\|$",
    ]
    inv_vals = [
        float(t4["row1_apex_antisymmetry"]),
        float(t4["row1_apex_torque"]),
        float(t4["row2_stumpff_horner_5stage_max_err"]),
        float(t4["row2_battin_wronskian"]),
        float(t4["row3_dirac_barycentric"]),
        float(t4["row4_supernode_multistep_reversibility"]),
    ]
    x_pos = np.arange(len(inv_labels))
    bar_colors_d = ["#1D4ED8", "#2563EB", "#0891B2", "#059669", "#D97706", "#DC2626"]
    ax_d.bar(x_pos, inv_vals, color=bar_colors_d, width=0.52)
    for xi, val in zip(x_pos, inv_vals):
        lbl_txt = f"{val:.2e}"
        ax_d.annotate(
            lbl_txt,
            (xi, val),
            textcoords="offset points",
            xytext=(0, 6),
            ha="center",
            fontsize=7.4,
            fontweight="bold",
            color="#1E293B",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.88),
        )
    ax_d.axhline(1e-14, color="red", ls="--", lw=1.2, label=r"Certification Gate ($10^{-14}$)")
    ax_d.axhline(1.11e-16, color="black", ls=":", lw=1.1, label=r"IEEE 754 $\epsilon_{\mathrm{mach}}$ ($1.11\times 10^{-16}$)")
    ax_d.set_yscale("log")
    ax_d.set_ylim(1e-18, 1e-11)
    ax_d.set_xticks(x_pos)
    ax_d.set_xticklabels(inv_labels, fontsize=7.6)
    ax_d.set_ylabel("Measured Physical / Symplectic Residual")
    ax_d.grid(True, axis="y", which="both", ls=":", alpha=0.45)
    ax_d.legend(loc="upper right", fontsize=7.8, framealpha=0.92)

    fig.subplots_adjust(left=0.07, right=0.93, top=0.90, bottom=0.08, wspace=0.30, hspace=0.32)
    FIG_DIR_CH06.mkdir(parents=True, exist_ok=True)
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    stem = "fig_ch06_microkernel_performance_and_invariants"
    for ext_fmt in ("png", "pdf"):
        out_ch06 = FIG_DIR_CH06 / f"{stem}.{ext_fmt}"
        fig.savefig(out_ch06, dpi=300, bbox_inches="tight")
        if ARTIFACT_DIR.resolve() != FIG_DIR_CH06.resolve():
            shutil.copy2(out_ch06, ARTIFACT_DIR / f"{stem}.{ext_fmt}")
    plt.close(fig)


def main() -> None:
    if not RESULTS_JSON.exists():
        raise FileNotFoundError(f"Missing certified JSON: {RESULTS_JSON}")
    with open(RESULTS_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10.5,
        "axes.labelsize": 11.0,
        "axes.titlesize": 11.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "mathtext.fontset": "cm",
    })

    plot_section_6_6_dedicated_plate(data)
    plot_volet_A(data)
    print("=== Successfully rendered canonical paper plates fig_ch06_microkernel_performance_and_invariants and fig_ch07_laskar_dual_precision_convergence ===", flush=True)


if __name__ == "__main__":
    main()
