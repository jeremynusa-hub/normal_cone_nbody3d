#!/usr/bin/env python3
"""
Cycle VI - Iteration 73: Publication-Grade 4-Panel Comparative Plate for Benchmark B1
Perturbed Star-Grazing / Ultra-High-Eccentricity Pericenter Passage (1 - e1 in [1e-1, 1e-9] and e1 = 1.0)
Normal-Cone Symplectic Splitting Solver (BM6 & PBM4, Float64) vs REBOUND 5.2.2 (WHFast, MERCURIUS, TRACE, IAS15, Float64).

Strict Scientific Integrity:
  - Zero floor/clamping on any data series.
  - Comparative wall-clock compute times included in legends, annotations, and Panel (d).
"""

from __future__ import annotations
import json
import os
from typing import Dict, Any, List

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


SOLVER_STYLES: Dict[str, Dict[str, Any]] = {
    # --- Category 1: Adaptive / Multi-Scale Solvers ---
    "normal_cone_bm6_dyadic": {
        "category": "adaptive_multiscale",
        "label": "Normal-Cone BM6 (6th-ord., Chart Flip + Dyadic)",
        "color": "#0072B2",
        "marker": "o",
        "mfc": "#0072B2",
        "ls": "-",
        "lw": 2.4,
        "ms": 6.5,
        "zorder": 10,
    },
    "normal_cone_pbm4_dyadic": {
        "category": "adaptive_multiscale",
        "label": "Normal-Cone PBM4 (4th/6th-ord. Gauge, Chart Flip + Dyadic)",
        "color": "#009E73",
        "marker": "D",
        "mfc": "#009E73",
        "ls": "-",
        "lw": 2.2,
        "ms": 5.8,
        "zorder": 9,
    },
    "rebound_ias15": {
        "category": "adaptive_multiscale",
        "label": "REBOUND IAS15 (15th-ord. Gauss-Radau, Adaptive)",
        "color": "#CC79A7",
        "marker": "P",
        "mfc": "#CC79A7",
        "ls": "-",
        "lw": 2.0,
        "ms": 6.2,
        "zorder": 8,
    },
    "rebound_trace_bs": {
        "category": "adaptive_multiscale",
        "label": "REBOUND TRACE (2nd DHC + 16th BS, Adaptive)",
        "color": "#D55E00",
        "marker": "X",
        "mfc": "#D55E00",
        "ls": "--",
        "lw": 1.9,
        "ms": 6.0,
        "zorder": 7,
    },
    "rebound_trace_ias15": {
        "category": "adaptive_multiscale",
        "label": "REBOUND TRACE (2nd DHC + 15th IAS15, Adaptive)",
        "color": "#8C510A",
        "marker": "d",
        "mfc": "#8C510A",
        "ls": ":",
        "lw": 1.8,
        "ms": 5.5,
        "zorder": 6,
    },
    "rebound_mercurius": {
        "category": "adaptive_multiscale",
        "label": "REBOUND MERCURIUS (2nd DHC + 15th IAS15, Adaptive)",
        "color": "#B2182B",
        "marker": "*",
        "mfc": "#B2182B",
        "ls": ":",
        "lw": 1.8,
        "ms": 7.5,
        "zorder": 5,
    },
    # --- Category 2: Fixed-Step Solvers ---
    "normal_cone_bm6_k0": {
        "category": "fixed_step",
        "label": "Normal-Cone BM6 (6th-ord., Chart Flip, Fixed Step $k=0$)",
        "color": "#56B4E9",
        "marker": "s",
        "mfc": "white",
        "ls": "--",
        "lw": 2.0,
        "ms": 5.8,
        "zorder": 9,
    },
    "normal_cone_pbm4_k0": {
        "category": "fixed_step",
        "label": "Normal-Cone PBM4 (4th/6th-ord. Gauge, Chart Flip, Fixed Step $k=0$)",
        "color": "#66C2A5",
        "marker": "^",
        "mfc": "white",
        "ls": "--",
        "lw": 1.9,
        "ms": 5.8,
        "zorder": 8,
    },
    "rebound_whfast_c17": {
        "category": "fixed_step",
        "label": "REBOUND WHFast (2nd+17th corr., $c=17$, Fixed Step)",
        "color": "#E69F00",
        "marker": "v",
        "mfc": "white",
        "ls": "-.",
        "lw": 1.9,
        "ms": 6.0,
        "zorder": 7,
    },
    "rebound_whfast_c0": {
        "category": "fixed_step",
        "label": "REBOUND WHFast (2nd-ord., $c=0$, Fixed Step)",
        "color": "#7F7F7F",
        "marker": "<",
        "mfc": "white",
        "ls": "--",
        "lw": 1.7,
        "ms": 5.8,
        "zorder": 6,
    },
}

ADAPTIVE_SOLVER_KEYS: List[str] = [
    "normal_cone_bm6_dyadic",
    "normal_cone_pbm4_dyadic",
    "rebound_ias15",
    "rebound_trace_bs",
    "rebound_trace_ias15",
    "rebound_mercurius",
]

FIXED_STEP_SOLVER_KEYS: List[str] = [
    "normal_cone_bm6_k0",
    "normal_cone_pbm4_k0",
    "rebound_whfast_c17",
    "rebound_whfast_c0",
]


def _add_categorized_header(ax: plt.Axes, title_text: str) -> None:
    """Adds an invisible line handle acting as a section header inside the legend."""
    ax.plot([], [], linestyle="none", marker="", label=title_text)


def main() -> None:
    json_path = "results/campaigns/rebound_b1_star_grazing_pericenter.json"
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    b1a = data["b1a_canonical_sweep"]
    b1b = data["b1b_permuted_sweep"]
    ts_data = data["timeseries_10orbits_1e_6"]
    T_orb = float(data["parameters"]["T_orb_grazer"])

    # Discrete x positions: index 0..6 corresponding to 1 - e1 = [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-9, 0.0 (e=1)]
    x_idx = np.arange(len(b1a), dtype=float)
    x_labels = [
        "$10^{-1}$",
        "$10^{-2}$",
        "$10^{-3}$",
        "$10^{-4}$",
        "$10^{-6}$",
        "$10^{-9}$",
        "$0$\n($e_1=1$)",
    ]

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.titlesize": 10.8,
        "axes.labelsize": 10.2,
        "legend.fontsize": 7.4,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
    })

    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(2, 2, figsize=(17.2, 13.0), dpi=300)
    fig.subplots_adjust(
        left=0.065, right=0.98, top=0.820, bottom=0.060, wspace=0.19, hspace=0.27
    )
    fig.suptitle(
        "Head-to-Head Benchmark #1 (B1 — Float64): Perturbed Star-Grazing Pericenter ($1 - e_1 \\in [10^{-1}, 10^{-9}]$ & $e_1 = 1.0$, Lu et al. 2024 TRACE §4.2)\n"
        "Category 1 (Adaptive / Multi-Scale: Normal-Cone BM6 & PBM4 [Chart Flip + Dyadic] vs IAS15, TRACE, MERCURIUS)  |  Category 2 (Fixed-Step: $k=0$ vs WHFast)",
        fontsize=10.5,
        fontweight="bold",
        y=0.978,
    )

    # Top shared figure legend (3 columns, exact visual match across all 4 panels)
    top_handles = [Line2D([], [], linestyle="none", label="[Cat. 1] Adaptive / Multi-Scale Solvers (filled markers):")]
    for skey in ADAPTIVE_SOLVER_KEYS:
        st = SOLVER_STYLES[skey]
        mean_ms = float(np.mean([float(entry["solvers"][skey]["ms_per_orbit_pure"]) for entry in b1a]))
        end_dE_e1 = float(b1a[-1]["solvers"][skey]["max_rel_dE"])
        top_handles.append(
            Line2D(
                [0], [0],
                color=st["color"], ls=st["ls"], lw=st["lw"],
                marker=st["marker"], ms=st["ms"], mfc=st["mfc"], mec=st["color"], mew=1.2,
                label=f"{st['label']} [{mean_ms:.2f} ms/orb | $e_1=1$: {end_dE_e1:.1e}]",
            )
        )
    top_handles.append(Line2D([], [], linestyle="none", label="[Cat. 2] Fixed-Step Solvers $k=0$ / Fixed $\\Delta t$ (open markers):"))
    for skey in FIXED_STEP_SOLVER_KEYS:
        st = SOLVER_STYLES[skey]
        mean_ms = float(np.mean([float(entry["solvers"][skey]["ms_per_orbit_pure"]) for entry in b1a]))
        end_dE_e1 = float(b1a[-1]["solvers"][skey]["max_rel_dE"])
        top_handles.append(
            Line2D(
                [0], [0],
                color=st["color"], ls=st["ls"], lw=st["lw"],
                marker=st["marker"], ms=st["ms"], mfc=st["mfc"], mec=st["color"], mew=1.4,
                label=f"{st['label']} [{mean_ms:.2f} ms/orb | $e_1=1$: {end_dE_e1:.1e}]",
            )
        )
    leg_top = fig.legend(
        handles=top_handles,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.944),
        ncol=3,
        frameon=True,
        framealpha=0.96,
        edgecolor="#cccccc",
        fontsize=6.9,
        columnspacing=0.75,
        handletextpad=0.35,
    )
    for txt in leg_top.get_texts():
        if txt.get_text().startswith("[Cat."):
            txt.set_weight("bold")
            txt.set_color("#111111")

    # ------------------------------------------------------------------------
    # Panel (a): B1-a Canonical Ordering [Star, Grazer(a=1.0), Perturber(a=5.2)]
    # ------------------------------------------------------------------------
    ax_a = axes[0, 0]
    for skey in ADAPTIVE_SOLVER_KEYS + FIXED_STEP_SOLVER_KEYS:
        st = SOLVER_STYLES[skey]
        y_vals = [float(entry["solvers"][skey]["max_rel_dE"]) for entry in b1a]
        ax_a.plot(
            x_idx,
            y_vals,
            color=st["color"],
            marker=st["marker"],
            markerfacecolor=st["mfc"],
            markeredgecolor=st["color"],
            markeredgewidth=1.4 if st["mfc"] == "white" else 1.1,
            linestyle=st["ls"],
            linewidth=st["lw"],
            markersize=st["ms"],
            zorder=st["zorder"],
        )

    ax_a.axvline(5.5, color="#555555", linestyle="--", linewidth=1.1, alpha=0.75)
    ax_a.text(
        5.55,
        2e-12,
        "Rectilinear\nCollision\n$e_1 = 1.0$",
        fontsize=7.8,
        color="#333333",
        va="bottom",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#F8F9FA", edgecolor="#CCCCCC", alpha=0.9),
    )
    bm6_e1_dE = float(b1a[-1]["solvers"]["normal_cone_bm6_dyadic"]["max_rel_dE"])
    pbm4_e1_dE = float(b1a[-1]["solvers"]["normal_cone_pbm4_dyadic"]["max_rel_dE"])
    ias15_e1_dE = float(b1a[-1]["solvers"]["rebound_ias15"]["max_rel_dE"])
    gain_ias15_e1 = ias15_e1_dE / max(bm6_e1_dE, 1e-30)
    ax_a.annotate(
        f"Flagship Domain Crossover ($1-e_1 \\leq 10^{{-6}} \\to e_1=1.0$):\n"
        f"IAS15 & TRACE diverge ($|\\Delta H/H_0| = {ias15_e1_dE:.1f}$ at $e_1=1.0$),\n"
        f"Normal-Cone BM6 (${bm6_e1_dE:.1e}$) & PBM4 (${pbm4_e1_dE:.1e}$) gain ${gain_ias15_e1:.1e}\\times$!",
        xy=(6.0, bm6_e1_dE),
        xytext=(1.75, 1.5e-15),
        arrowprops=dict(arrowstyle="->", color="#0072B2", lw=1.2),
        fontsize=7.2,
        fontweight="bold",
        color="#004B76",
        bbox=dict(boxstyle="round,pad=0.20", facecolor="#EEF7FC", edgecolor="#0072B2", alpha=0.95),
    )
    ax_a.set_yscale("log")
    ax_a.set_ylim(1e-16, 3e3)
    ax_a.set_xticks(x_idx)
    ax_a.set_xticklabels(x_labels)
    ax_a.set_xlabel("Pericenter Factor $1 - e_1 = r_{\\mathrm{peri}} / a_1$ (Decreasing Pericenter Distance $\\rightarrow$)")
    ax_a.set_ylabel("Peak Relative Energy Error $\\max_t |\\Delta H(t) / H_0|$")
    ax_a.set_title("(a) B1-a Canonical Hierarchy $[\\mathrm{Star}, \\mathrm{Grazer}(1\\,\\mathrm{AU}), \\mathrm{Giant}(5.2\\,\\mathrm{AU})]$ (5 Orbits)")
    ax_a.grid(True, which="both", linestyle=":", alpha=0.45)

    # ------------------------------------------------------------------------
    # Panel (b): B1-b Permuted Ordering [Star, Perturber(a=5.2), Grazer(a=1.0)]
    # ------------------------------------------------------------------------
    ax_b = axes[0, 1]
    for skey in ADAPTIVE_SOLVER_KEYS + FIXED_STEP_SOLVER_KEYS:
        st = SOLVER_STYLES[skey]
        y_vals = [float(entry["solvers"][skey]["max_rel_dE"]) for entry in b1b]
        ax_b.plot(
            x_idx,
            y_vals,
            color=st["color"],
            marker=st["marker"],
            markerfacecolor=st["mfc"],
            markeredgecolor=st["color"],
            markeredgewidth=1.4 if st["mfc"] == "white" else 1.1,
            linestyle=st["ls"],
            linewidth=st["lw"],
            markersize=st["ms"],
            zorder=st["zorder"],
        )

    ax_b.axvline(5.5, color="#555555", linestyle="--", linewidth=1.1, alpha=0.75)
    wh_b1b_peak = max(float(row["solvers"]["rebound_whfast_c17"]["max_rel_dE"]) for row in b1b)
    bm6_b1b_e1 = float(b1b[-1]["solvers"]["normal_cone_bm6_dyadic"]["max_rel_dE"])
    pbm4_b1b_4 = float(b1b[4]["solvers"]["normal_cone_pbm4_dyadic"]["max_rel_dE"])
    ax_b.annotate(
        f"Fixed-Step WHFast Static Jacobi Breakdown ($10^1\\text{{–}}{wh_b1b_peak:.1e}$)\nwhen Star-Grazer is at Particle Index 2",
        xy=(3.0, float(b1b[3]["solvers"]["rebound_whfast_c17"]["max_rel_dE"])),
        xytext=(1.35, 3.0e1),
        arrowprops=dict(arrowstyle="->", color="#E69F00", lw=1.2),
        fontsize=7.8,
        fontweight="bold",
        color="#8C2D04",
        bbox=dict(boxstyle="round,pad=0.24", facecolor="#FFF5EB", edgecolor="#E69F00", alpha=0.95),
    )
    ax_b.annotate(
        f"Normal-Cone AtlasManager Auto-Flips to Leaf $[(0,2)]$\n($|\\Delta H/H_0| = {bm6_b1b_e1:.1e}$ for BM6 at $e_2 = 1.0$, $0$ saturation)",
        xy=(4.0, pbm4_b1b_4),
        xytext=(2.20, 2.0e-14),
        arrowprops=dict(arrowstyle="->", color="#0072B2", lw=1.2),
        fontsize=7.4,
        fontweight="bold",
        color="#004B76",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="#EEF7FC", edgecolor="#0072B2", alpha=0.95),
    )
    ax_b.set_yscale("log")
    ax_b.set_ylim(1e-16, 3e6)
    ax_b.set_xticks(x_idx)
    ax_b.set_xticklabels(x_labels)
    ax_b.set_xlabel("Pericenter Factor $1 - e_2 = r_{\\mathrm{peri}} / a_2$ (Decreasing Pericenter Distance $\\rightarrow$)")
    ax_b.set_ylabel("Peak Relative Energy Error $\\max_t |\\Delta H(t) / H_0|$")
    ax_b.set_title("(b) B1-b Topological Index Permutation $[\\mathrm{Star}, \\mathrm{Giant}(5.2\\,\\mathrm{AU}), \\mathrm{Grazer}(1\\,\\mathrm{AU})]$")
    ax_b.grid(True, which="both", linestyle=":", alpha=0.45)

    # ------------------------------------------------------------------------
    # Panel (c): 10-Orbit Time Series & Round-Trip Reversibility at 1 - e1 = 1e-6
    # ------------------------------------------------------------------------
    ax_c = axes[1, 0]
    for skey in ADAPTIVE_SOLVER_KEYS + FIXED_STEP_SOLVER_KEYS:
        st = SOLVER_STYLES[skey]
        r_ts = ts_data[skey]
        t_orb_arr = np.asarray(r_ts["series_t"], dtype=float) / T_orb
        y_arr = np.asarray(r_ts["series_rel_dE"], dtype=float)
        stride_c = max(1, len(t_orb_arr) // 18)
        ax_c.plot(
            t_orb_arr,
            y_arr,
            color=st["color"],
            linestyle=st["ls"],
            linewidth=st["lw"] * 0.9,
            marker=st["marker"],
            markersize=st["ms"] * 0.85,
            markerfacecolor=st["mfc"],
            markeredgecolor=st["color"],
            markeredgewidth=1.3 if st["mfc"] == "white" else 1.0,
            markevery=stride_c,
            alpha=0.92,
            zorder=st["zorder"],
        )

    rev_bm6 = float(ts_data["normal_cone_bm6_dyadic"]["rev_dq_norm"])
    rev_pbm4 = float(ts_data["normal_cone_pbm4_dyadic"]["rev_dq_norm"])
    rev_ias = float(ts_data["rebound_ias15"]["rev_dq_norm"])
    rev_tr_bs = float(ts_data["rebound_trace_bs"]["rev_dq_norm"])
    rev_merc = float(ts_data["rebound_mercurius"]["rev_dq_norm"])
    gain_rev_ias = rev_ias / max(rev_pbm4, 1e-30)
    ax_c.text(
        0.97, 0.06,
        f"10-Orbit Round-Trip Reversibility $\\|\\Delta\\mathbf{{q}}_{{\\mathrm{{FB}}}}\\|_2$ at $1-e_1 = 10^{{-6}}$:\n"
        f"• Normal-Cone PBM4 (Dyadic): ${rev_pbm4:.2e}$ AU  |  BM6 (Dyadic): ${rev_bm6:.2e}$ AU\n"
        f"• REBOUND IAS15: ${rev_ias:.2e}$ AU (${gain_rev_ias:,.0f}\\times$ worse than PBM4/BM6!)\n"
        f"• REBOUND TRACE (BS): ${rev_tr_bs:.2e}$ AU  |  MERCURIUS: ${rev_merc:.2e}$ AU",
        transform=ax_c.transAxes,
        fontsize=7.4,
        va="bottom",
        ha="right",
        color="#003366",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#F4F9FF", edgecolor="#0072B2", alpha=0.95),
    )

    ax_c.set_yscale("log")
    ax_c.set_xlim(0.0, 10.0)
    ax_c.set_ylim(1e-17, 1e2)
    ax_c.set_xlabel("Elapsed Time $t / T_1$ (Number of Star-Grazing Pericenter Passages at $r_{\\mathrm{peri}} = 10^{-6}\\,\\mathrm{AU}$)")
    ax_c.set_ylabel("Instantaneous Relative Energy Error $|\\Delta H(t) / H_0|$")
    ax_c.set_title("(c) 10-Orbit Evolution & Round-Trip Reversibility $\\|\\Delta\\mathbf{q}_{\\mathrm{FB}}\\|_2$ at $1 - e_1 = 10^{-6}$")
    ax_c.grid(True, which="both", linestyle=":", alpha=0.45)

    # ------------------------------------------------------------------------
    # Panel (d): Comparative Compute Time (ms / orbit) Across All Eccentricities
    # ------------------------------------------------------------------------
    ax_d = axes[1, 1]
    for skey in ADAPTIVE_SOLVER_KEYS + FIXED_STEP_SOLVER_KEYS:
        st = SOLVER_STYLES[skey]
        ms_vals = [float(entry["solvers"][skey]["ms_per_orbit_pure"]) for entry in b1a]
        ax_d.plot(
            x_idx,
            ms_vals,
            color=st["color"],
            marker=st["marker"],
            markerfacecolor=st["mfc"],
            markeredgecolor=st["color"],
            markeredgewidth=1.4 if st["mfc"] == "white" else 1.1,
            linestyle=st["ls"],
            linewidth=st["lw"],
            markersize=st["ms"],
            zorder=st["zorder"],
        )

    ax_d.axvline(5.5, color="#555555", linestyle="--", linewidth=1.1, alpha=0.75)
    bm6_ms_list = [float(row["solvers"]["normal_cone_bm6_dyadic"]["ms_per_orbit_pure"]) for row in b1a]
    pbm4_ms_list = [float(row["solvers"]["normal_cone_pbm4_dyadic"]["ms_per_orbit_pure"]) for row in b1a]
    max_k_b1a = max(int(row["solvers"]["normal_cone_bm6_dyadic"]["max_dyadic_level"]) for row in b1a)
    merc_ms_list = [float(row["solvers"]["rebound_mercurius"]["ms_per_orbit_pure"]) for row in b1a]
    merc_peak_idx = int(np.argmax(merc_ms_list))
    merc_peak_ms = merc_ms_list[merc_peak_idx]
    merc_slowdown = merc_peak_ms / max(merc_ms_list[0], 1e-9)
    ax_d.annotate(
        f"Normal-Cone BM6 & PBM4 (Numba, Chart Flip + Dyadic):\n"
        f"flat $\\sim {min(bm6_ms_list):.0f}\\text{{–}}{max(pbm4_ms_list):.0f}\\,\\mathrm{{ms/orb}}$ ($k_{{\\max}}\\leq {max_k_b1a}$, 0 saturation down to $e_1 = 1.0$)",
        xy=(1.0, bm6_ms_list[1]),
        xytext=(0.15, 2.5e-2),
        arrowprops=dict(arrowstyle="->", color="#0072B2", lw=1.2),
        fontsize=7.5,
        color="#004B76",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="#EEF7FC", edgecolor="#0072B2", alpha=0.95),
    )
    ax_d.annotate(
        f"MERCURIUS (C99 Adaptive): slows down ${merc_slowdown:.0f}\\times$\n"
        f"($\\sim {merc_peak_ms:.0f}\\,\\mathrm{{ms/orb}}$ at $10^{{-9}}$, while $|\\Delta H/H_0| = \\mathcal{{O}}(1)$!)",
        xy=(float(merc_peak_idx), merc_peak_ms),
        xytext=(1.85, 3.0e-3),
        arrowprops=dict(arrowstyle="->", color="#B2182B", lw=1.2),
        fontsize=7.5,
        color="#7F0000",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="#FFF0F0", edgecolor="#B2182B", alpha=0.95),
    )
    ax_d.set_yscale("log")
    ax_d.set_ylim(1.0e-3, 4.0e3)
    ax_d.set_xticks(x_idx)
    ax_d.set_xticklabels(x_labels)
    ax_d.set_xlabel("Pericenter Factor $1 - e_1 = r_{\\mathrm{peri}} / a_1$ (Decreasing Pericenter Distance $\\rightarrow$)")
    ax_d.set_ylabel("Pure Integration Wall-Clock Time per Orbit [ms / orbit]")
    ax_d.set_title("(d) Comparative Compute Time per Orbit vs Pericenter Factor ($N_{\\mathrm{steps}}=200$/orbit)")
    ax_d.grid(True, which="both", linestyle=":", alpha=0.45)

    out_dir = os.path.join("results", "figures")
    os.makedirs(out_dir, exist_ok=True)
    png_path = os.path.join(out_dir, "fig_rebound_b1_star_grazing_head_to_head.png")
    pdf_path = os.path.join(out_dir, "fig_rebound_b1_star_grazing_head_to_head.pdf")
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Saved {png_path} and {pdf_path}")


if __name__ == "__main__":
    main()
