#!/usr/bin/env python3
"""
Generates the 4-panel (2x2) publication plate for Benchmark #3 (B4):
Definitive Reversibility of the Chaotic Pythagorean 3-4-5 Problem Beyond the Float64 Lyapunov Horizon
(Burrau 1913; Szebehely & Peters 1967; Boekholt & Portegies Zwart 2015 Brutus; Rein & Spiegel 2015 IAS15;
Rein et al. 2019 MERCURIUS; Lu et al. 2024 TRACE), comparing REBOUND 5.2.2 against our
Unsoftened Normal-Cone Symplectic Splitting Solver (Normal-Cone BM6 & Normal-Cone PBM4)
in IEEE-754 Float64 and 35-digit Float128.

Strict compliance rules:
  - Zero arbitrary clamping or floors (no np.maximum(..., 1e-17)).
  - Canonical Normal-Cone vs Conormal terminology (conormal_vs_normal_cone_foundations.md).
  - Explicit separation of Category 1 (Adaptive / Multi-Scale Solvers, filled markers)
    and Category 2 (Fixed-Step Solvers k=0 / Fixed dt, open markers).
  - Unobstructed legends on ALL 4 panels placed in dedicated bottom bands via extended ylim.
"""

from __future__ import annotations
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def _positive_xy(x_list, y_list):
    xs, ys = [], []
    for x, y in zip(x_list, y_list):
        if y is not None and np.isfinite(y) and y > 0.0:
            xs.append(float(x))
            ys.append(float(y))
    return np.array(xs, dtype=np.float64), np.array(ys, dtype=np.float64)


def main() -> None:
    json_candidates = [
        "rebound_b4_pythagorean_chaos.json",
        os.path.join("results", "campaigns", "rebound_b4_pythagorean_chaos.json"),
    ]
    json_path = None
    for c in json_candidates:
        if os.path.exists(c):
            json_path = c
            break
    if json_path is None:
        raise FileNotFoundError("Could not find rebound_b4_pythagorean_chaos.json")

    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    b4a = data["protocol_B4a_forward_chaos"]
    rows_b4b = data["protocol_B4b_reversibility_sweep"]

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.titlesize": 10.2,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 7.1,
        "figure.titlesize": 11.5,
    })

    fig, axes = plt.subplots(2, 2, figsize=(17.2, 13.0), dpi=300)
    fig.subplots_adjust(
        left=0.065, right=0.965, top=0.825, bottom=0.055, wspace=0.20, hspace=0.26
    )

    fig.suptitle(
        "Head-to-Head Benchmark #3 (B4 — Float64 & Float128): Pythagorean 3–4–5 Chaotic Scattering & Time-Reversibility\n"
        "Category 1 (Adaptive / Multi-Scale: Normal-Cone BM6 & PBM4 [Chart Flip + Dyadic, Float64 & Float128] vs REBOUND IAS15, TRACE, MERCURIUS)  |  Category 2 (Fixed-Step: k=0 vs WHFast)",
        fontsize=10.1,
        fontweight="bold",
        y=0.978,
    )

    style_map = {
        # Category 1: Adaptive / Multi-Scale Solvers in Float128 & Float64 (filled markers)
        "normal_cone_bm6_f128_dyadic": {
            "category": "adaptive_multiscale",
            "color": "#00264d", "marker": "*", "ls": "-", "lw": 2.5, "ms": 7.8, "mfc": "#00264d", "mew": 1.0,
            "name": "Normal-Cone BM6 Float128 (6th, Chart Flip + Dyadic)",
        },
        "normal_cone_pbm4_f128_dyadic": {
            "category": "adaptive_multiscale",
            "color": "#006644", "marker": "P", "ls": "-", "lw": 2.2, "ms": 6.6, "mfc": "#006644", "mew": 1.0,
            "name": "Normal-Cone PBM4 Float128 (4th/6th Gauge, Chart Flip + Dyadic)",
        },
        "normal_cone_bm6_f64_dyadic": {
            "category": "adaptive_multiscale",
            "color": "#0055b3", "marker": "^", "ls": "-", "lw": 2.4, "ms": 6.8, "mfc": "#0055b3", "mew": 1.0,
            "name": "Normal-Cone BM6 Float64 (6th, Chart Flip + Dyadic)",
        },
        "normal_cone_pbm4_f64_dyadic": {
            "category": "adaptive_multiscale",
            "color": "#009e73", "marker": "p", "ls": "-", "lw": 2.1, "ms": 6.4, "mfc": "#009e73", "mew": 1.0,
            "name": "Normal-Cone PBM4 Float64 (4th/6th Gauge, Chart Flip + Dyadic)",
        },
        "rebound_ias15": {
            "category": "adaptive_multiscale",
            "color": "#222222", "marker": "o", "ls": ":", "lw": 2.3, "ms": 5.8, "mfc": "#222222", "mew": 1.0,
            "name": "REBOUND IAS15 Float64 (15th Gauss-Radau, Adaptive)",
        },
        "rebound_trace_bs": {
            "category": "adaptive_multiscale",
            "color": "#8c564b", "marker": "D", "ls": "-", "lw": 1.7, "ms": 5.0, "mfc": "#8c564b", "mew": 1.0,
            "name": "REBOUND TRACE Float64 (2nd DHC + 16th BS, Adaptive)",
        },
        "rebound_trace_ias15": {
            "category": "adaptive_multiscale",
            "color": "#e377c2", "marker": "h", "ls": "-.", "lw": 1.6, "ms": 5.0, "mfc": "#e377c2", "mew": 1.0,
            "name": "REBOUND TRACE Float64 (2nd DHC + 15th IAS15, Adaptive)",
        },
        "rebound_mercurius": {
            "category": "adaptive_multiscale",
            "color": "#9467bd", "marker": "s", "ls": "-", "lw": 1.8, "ms": 5.2, "mfc": "#9467bd", "mew": 1.0,
            "name": "REBOUND MERCURIUS Float64 (2nd DHC + 15th IAS15, Adaptive)",
        },
        # Category 2: Fixed-Step Solvers in Float64 (open markers)
        "normal_cone_bm6_f64_k0": {
            "category": "fixed_step",
            "color": "#1f77b4", "marker": "o", "ls": "--", "lw": 1.8, "ms": 5.6, "mfc": "white", "mew": 1.5,
            "name": "Normal-Cone BM6 Float64 (6th, Chart Flip, Fixed k=0)",
        },
        "normal_cone_pbm4_f64_k0": {
            "category": "fixed_step",
            "color": "#2ca02c", "marker": "^", "ls": "--", "lw": 1.8, "ms": 5.6, "mfc": "white", "mew": 1.5,
            "name": "Normal-Cone PBM4 Float64 (4th/6th Gauge, Chart Flip, Fixed k=0)",
        },
        "rebound_whfast_c17": {
            "category": "fixed_step",
            "color": "#ff7f0e", "marker": "v", "ls": "--", "lw": 1.6, "ms": 5.4, "mfc": "white", "mew": 1.4,
            "name": "REBOUND WHFast Float64 (2nd+17th corr., c=17, Fixed)",
        },
        "rebound_whfast_c0": {
            "category": "fixed_step",
            "color": "#d62728", "marker": "x", "ls": "--", "lw": 1.6, "ms": 5.6, "mfc": "#d62728", "mew": 1.5,
            "name": "REBOUND WHFast Float64 (2nd-ord., c=0, Fixed)",
        },
    }

    cat1_order = [
        "normal_cone_bm6_f128_dyadic",
        "normal_cone_pbm4_f128_dyadic",
        "normal_cone_bm6_f64_dyadic",
        "normal_cone_pbm4_f64_dyadic",
        "rebound_ias15",
        "rebound_trace_bs",
        "rebound_trace_ias15",
        "rebound_mercurius",
    ]
    cat2_order = [
        "normal_cone_bm6_f64_k0",
        "normal_cone_pbm4_f64_k0",
        "rebound_whfast_c17",
        "rebound_whfast_c0",
    ]

    # Top shared figure legend (3 columns, strictly Float64 vs Float64)
    top_handles = [Line2D([], [], linestyle="none", label="[Cat. 1] Adaptive / Multi-Scale (filled markers):")]
    for sk in cat1_order:
        st = style_map[sk]
        ms_str = f" [{b4a[sk]['ms_per_time_unit']:.1f} ms/t]" if sk in b4a else ""
        top_handles.append(
            Line2D([0], [0], color=st["color"], ls=st["ls"], lw=st["lw"],
                   marker=st["marker"], ms=st["ms"], mfc=st["mfc"], mew=st["mew"],
                   label=f"{st['name']}{ms_str}")
        )
    top_handles.append(Line2D([], [], linestyle="none", label="[Cat. 2] Fixed-Step k=0 / dt=0.04 (open markers):"))
    for sk in cat2_order:
        st = style_map[sk]
        ms_str = f" [{b4a[sk]['ms_per_time_unit']:.1f} ms/t]" if sk in b4a else ""
        top_handles.append(
            Line2D([0], [0], color=st["color"], ls=st["ls"], lw=st["lw"],
                   marker=st["marker"], ms=st["ms"], mfc=st["mfc"], mew=st["mew"],
                   label=f"{st['name']}{ms_str}")
        )
    fig.legend(
        handles=top_handles,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.942),
        ncol=3,
        frameon=True,
        framealpha=0.95,
        edgecolor="#cccccc",
        fontsize=7.6,
    )

    # ========================================================================
    # Panel (a): B4-a Relative Energy Error |dH(t)|/|H0| over Chaotic Scattering
    # ========================================================================
    ax_a = axes[0, 0]
    for sk in cat1_order + cat2_order:
        if sk not in b4a:
            continue
        res = b4a[sk]
        st = style_map[sk]
        ts, ys = _positive_xy(res["series_t"], res["series_rel_dE"])
        if len(ts) == 0:
            continue
        stride = max(1, len(ts) // 18)
        ax_a.plot(
            ts, ys,
            color=st["color"], ls=st["ls"], lw=st["lw"],
            marker=st["marker"], ms=st["ms"], mfc=st["mfc"], mew=st["mew"],
            markevery=stride, alpha=0.92,
        )

    ax_a.axvline(59.42, color="#333333", ls="--", lw=1.1, alpha=0.7)
    ax_a.text(58.0, 3e-14, "Final Triple Encounter\n& Ejection ($t \\approx 59.4$)", fontsize=8.0, ha="right", va="center",
              color="#222222", bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#cccccc", alpha=0.92))
    ax_a.set_yscale("log")
    ax_a.set_xlim(0.0, 70.0)
    ax_a.set_ylim(1e-34, 1e4)
    ax_a.set_xlabel("Time $t$ (N-body units, $G=1, m=(3,4,5)$)")
    ax_a.set_ylabel("Relative Energy Error $|\\Delta H(t)| / |H_0|$")
    ax_a.set_title("(a) Protocol B4-a: Energy Conservation Under Democratic 3–4–5 Chaotic Scattering")
    ax_a.grid(True, which="both", ls=":", alpha=0.45)

    # ========================================================================
    # Panel (b): B4-b Round-Trip Time-Reversibility ||dq_FB(2T)||_2 vs Horizon T
    #            (Clean single y-axis [1e-34, 1e2], Float64 & Float128 solvers)
    # ========================================================================
    ax_b = axes[0, 1]
    b4b_solvers = [
        "normal_cone_bm6_f128_dyadic",
        "normal_cone_pbm4_f128_dyadic",
        "normal_cone_bm6_f64_dyadic",
        "normal_cone_pbm4_f64_dyadic",
        "normal_cone_bm6_f64_k0",
        "normal_cone_pbm4_f64_k0",
        "rebound_ias15",
        "rebound_trace_bs",
        "rebound_trace_ias15",
    ]
    for sk in b4b_solvers:
        st = style_map[sk]
        T_vals, dq_vals = [], []
        for r in rows_b4b:
            if sk in r["solvers"]:
                val = r["solvers"][sk]["rev_dq_norm"]
                if val is not None and np.isfinite(val) and val > 0.0:
                    T_vals.append(float(r["t_turnaround"]))
                    dq_vals.append(float(val))
        if len(T_vals) > 0:
            ax_b.plot(
                T_vals, dq_vals,
                color=st["color"], ls=st["ls"], lw=st["lw"],
                marker=st["marker"], ms=st["ms"], mfc=st["mfc"], mew=st["mew"],
                alpha=0.94,
            )

    bm6_120_f128 = float(rows_b4b[-1]["solvers"]["normal_cone_bm6_f128_dyadic"]["rev_dq_norm"])
    pbm4_120_f128 = float(rows_b4b[-1]["solvers"]["normal_cone_pbm4_f128_dyadic"]["rev_dq_norm"])
    ias_120 = float(rows_b4b[-1]["solvers"]["rebound_ias15"]["rev_dq_norm"])
    ax_b.axhline(1.0, color="#b22222", ls="--", lw=1.2, alpha=0.8)
    ax_b.text(119.0, 2.2, f"$\\mathcal{{O}}(1)$ Macroscopic Reversibility Loss Threshold (IAS15 = ${ias_120:.2f}$ at $T=120$)",
              fontsize=7.6, color="#b22222", ha="right", va="bottom",
              bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#e6b0aa", alpha=0.92))
    ax_b.text(
        48.0, 1e-33,
        f"Normal-Cone BM6 & PBM4 Float128 (Dyadic, 0 sat):\nExact reversibility beyond Float64 Lyapunov horizon\n(${bm6_120_f128:.2e}$ & ${pbm4_120_f128:.2e}$ at $T=120$!)",
        fontsize=7.2, fontweight="bold", color="#00264d", ha="left", va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", fc="#eef6ff", ec="#00264d", alpha=0.94),
    )
    ax_b.set_yscale("log")
    ax_b.set_xlim(0.0, 123.0)
    ax_b.set_ylim(1e-34, 1e2)
    ax_b.set_xlabel("Turnaround Horizon $T$ in Round-Trip Test ($0 \\to T \\to 0$)")
    ax_b.set_ylabel("Round-Trip Position Error $\\|\\mathbf{q}_{\\mathrm{ret}}(0) - \\mathbf{q}(0)\\|_2$")
    ax_b.set_title("(b) Protocol B4-b: Round-Trip Time-Reversibility Error (Float64 & Float128)")
    ax_b.grid(True, which="both", ls=":", alpha=0.45)

    # ========================================================================
    # Panel (c): B4-a Minimum Inter-Body Separation r_min(t), Dyadic Level k(t) & Active Jacobi Charts
    # ========================================================================
    ax_c = axes[1, 0]
    ref_nc = b4a.get("normal_cone_bm6_f64_dyadic") or b4a.get("normal_cone_pbm4_f64_dyadic")
    if ref_nc is not None:
        ts_c = np.array(ref_nc["series_t"], dtype=np.float64)
        rmin_c = np.array(ref_nc["series_r_min"], dtype=np.float64)
        k_c = np.array(ref_nc["series_dyadic_k"], dtype=np.float64)
        lps = ref_nc.get("series_leaf_pair", [])

        ax_c.plot(ts_c, rmin_c, color="#0055b3", lw=1.6, label="Normal-Cone BM6: Min Separation $r_{\\min}(t) = \\min_{i<j} r_{ij}(t)$")
        if "rebound_ias15" in b4a:
            ts_ias = np.array(b4a["rebound_ias15"]["series_t"], dtype=np.float64)
            rm_ias = np.array(b4a["rebound_ias15"]["series_r_min"], dtype=np.float64)
            ax_c.plot(ts_ias, rm_ias, color="#222222", ls=":", lw=1.3, alpha=0.85, label="REBOUND IAS15: Min Separation $r_{\\min}(t)$")

        # Shade active democratic Jacobi chart leaf pair along top strip
        chart_colors = {"(0, 1)": "#ffcccc", "(0, 2)": "#cce5ff", "(1, 2)": "#d5f5e3"}
        if len(lps) == len(ts_c):
            cur_pair = str(tuple(sorted(lps[0])))
            t_start_seg = ts_c[0]
            for idx in range(1, len(ts_c)):
                p_str = str(tuple(sorted(lps[idx])))
                if p_str != cur_pair or idx == len(ts_c) - 1:
                    t_end_seg = ts_c[idx]
                    ax_c.axvspan(t_start_seg, t_end_seg, ymin=0.88, ymax=0.98,
                                 color=chart_colors.get(cur_pair, "#eeeeee"), alpha=0.85)
                    cur_pair = p_str
                    t_start_seg = t_end_seg

        ax_c2 = ax_c.twinx()
        ax_c2.step(ts_c, k_c, where="post", color="#d95f02", lw=1.2, alpha=0.65, label="Dyadic Substepping Level $k_{\\mathrm{sym}}(t)$")
        ax_c2.set_ylabel("Dyadic Level $k_{\\mathrm{sym}}(t)$", color="#d95f02")
        ax_c2.set_ylim(-2.0, max(12.0, float(np.max(k_c)) + 2.0))
        ax_c2.tick_params(axis="y", labelcolor="#d95f02")

    ax_c.set_yscale("log")
    ax_c.set_xlim(0.0, 70.0)
    ax_c.set_ylim(1e-6, 1e2)
    ax_c.set_xlabel("Time $t$ (N-body units)")
    ax_c.set_ylabel("Minimum Pairwise Separation $r_{\\min}(t)$")
    n_flips_bm6 = ref_nc["n_atlas_flips"] if ref_nc else 0
    ax_c.set_title(f"(c) Democratic Jacobi Chart Flips ({n_flips_bm6} Transitions), Close Encounters $r_{{\\min}}(t)$ & Dyadic $k(t)$")
    ax_c.grid(True, which="both", ls=":", alpha=0.45)

    from matplotlib.patches import Patch
    handles_c = [
        Line2D([0], [0], color="#0055b3", lw=1.6, label="Normal-Cone BM6: $r_{\\min}(t)$"),
        Line2D([0], [0], color="#222222", ls=":", lw=1.3, label="REBOUND IAS15: $r_{\\min}(t)$"),
        Line2D([0], [0], color="#d95f02", lw=1.2, label="Active Dyadic Level $k_{\\mathrm{sym}}(t)$ (right axis)"),
        Patch(facecolor="#d5f5e3", edgecolor="#2ca02c", label="Jacobi Chart $\\mathcal{T}_{12} = [(1, 2)]$ (Bodies $m=4,5$)"),
        Patch(facecolor="#cce5ff", edgecolor="#1f77b4", label="Jacobi Chart $\\mathcal{T}_{02} = [(0, 2)]$ (Bodies $m=3,5$)"),
        Patch(facecolor="#ffcccc", edgecolor="#d62728", label="Jacobi Chart $\\mathcal{T}_{01} = [(0, 1)]$ (Bodies $m=3,4$)"),
    ]
    ax_c.legend(handles=handles_c, loc="lower left", ncol=2, frameon=True, framealpha=0.93, fontsize=7.0)

    # ========================================================================
    # Panel (d): Dedicated Reversibility Gain Factor Over REBOUND IAS15
    #            G_solver(T) = ||dq_FB(IAS15)||_2 / ||dq_FB(Solver)||_2
    # ========================================================================
    ax_d = axes[1, 1]
    gain_solvers = [
        "normal_cone_bm6_f128_dyadic",
        "normal_cone_pbm4_f128_dyadic",
        "normal_cone_bm6_f64_dyadic",
        "normal_cone_pbm4_f64_dyadic",
        "normal_cone_bm6_f64_k0",
        "normal_cone_pbm4_f64_k0",
        "rebound_trace_bs",
        "rebound_trace_ias15",
    ]
    ax_d.axhspan(1.0, 1e26, color="#e8f8f5", alpha=0.55)
    ax_d.axhline(1.0, color="#222222", ls=":", lw=1.8, alpha=0.9)
    ax_d.text(
        3.5, 2.5e24,
        "Green Zone ($> 10^0$): Solver MORE Reversible than REBOUND IAS15",
        fontsize=7.8, fontweight="bold", color="#006644", ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.22", fc="#f2fbf9", ec="#006644", lw=0.8, alpha=0.94),
    )

    for sk in gain_solvers:
        st = style_map[sk]
        T_g, G_g = [], []
        for r in rows_b4b:
            sols = r["solvers"]
            if "rebound_ias15" in sols and sk in sols:
                dq_ias = sols["rebound_ias15"]["rev_dq_norm"]
                dq_sol = sols[sk]["rev_dq_norm"]
                if dq_ias is not None and dq_sol is not None and dq_ias > 0.0 and dq_sol > 0.0:
                    T_g.append(float(r["t_turnaround"]))
                    G_g.append(float(dq_ias / dq_sol))
        if len(T_g) > 0:
            ax_d.plot(
                T_g, G_g,
                color=st["color"], ls=st["ls"], lw=st["lw"],
                marker=st["marker"], ms=st["ms"], mfc=st["mfc"], mew=st["mew"],
                alpha=0.94,
            )
            if sk == "normal_cone_bm6_f128_dyadic":
                for t_target, dy_fac in [(80.0, 3.0), (120.0, 3.0)]:
                    matches = [g for t_val, g in zip(T_g, G_g) if abs(t_val - t_target) < 1e-6]
                    if matches:
                        g_val = matches[0]
                        ax_d.text(
                            t_target, g_val * dy_fac,
                            f"BM6 F128: $\\mathbf{{{g_val:.1e}\\times}}$\n($T={int(t_target)}$)",
                            color="#00264d", fontsize=7.1, ha="center", va="bottom",
                            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#00264d", lw=0.8, alpha=0.94),
                        )
            elif sk == "normal_cone_pbm4_f128_dyadic":
                for t_target, dy_fac in [(100.0, 0.04)]:
                    matches = [g for t_val, g in zip(T_g, G_g) if abs(t_val - t_target) < 1e-6]
                    if matches:
                        g_val = matches[0]
                        ax_d.text(
                            t_target, g_val * dy_fac,
                            f"PBM4 F128: $\\mathbf{{{g_val:.1e}\\times}}$\n($T={int(t_target)}$)",
                            color="#006644", fontsize=7.1, ha="center", va="top",
                            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#006644", lw=0.8, alpha=0.94),
                        )
            elif sk == "normal_cone_pbm4_f64_dyadic":
                for t_target, dy_fac in [(40.0, 0.04), (120.0, 25.0)]:
                    matches = [g for t_val, g in zip(T_g, G_g) if abs(t_val - t_target) < 1e-6]
                    if matches:
                        g_val = matches[0]
                        lbl_g = f"{g_val:.1f}\\times" if g_val >= 0.1 else f"{g_val:.1e}\\times"
                        ax_d.text(
                            t_target, g_val * dy_fac,
                            f"PBM4 F64: $\\mathbf{{{lbl_g}}}$\n($T={int(t_target)}$)",
                            color="#009e73", fontsize=7.0, ha="center", va="bottom" if dy_fac > 1 else "top",
                            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#009e73", lw=0.8, alpha=0.94),
                        )

    ax_d.set_yscale("log")
    ax_d.set_xlim(0.0, 123.0)
    ax_d.set_ylim(1e-7, 1e26)
    ax_d.set_xlabel("Turnaround Horizon $T$ in Round-Trip Test ($0 \\to T \\to 0$)")
    ax_d.set_ylabel("Reversibility Gain vs IAS15 $\\|\\Delta\\mathbf{q}_{\\mathrm{FB}}^{\\mathrm{IAS15}}\\|_2 / \\|\\Delta\\mathbf{q}_{\\mathrm{FB}}^{\\mathrm{Solver}}\\|_2$")
    ax_d.set_title("(d) Protocol B4-b: Reversibility Gain Factor Relative to REBOUND IAS15 (Float64 & Float128)")
    ax_d.grid(True, which="both", ls=":", alpha=0.45)

    out_dirs = [
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "figures"),
        # (root duplicate removed)
    ]
    for d in out_dirs:
        if os.path.isdir(d):
            png_p = os.path.join(d, "fig_rebound_b4_pythagorean_chaos_head_to_head.png")
            pdf_p = os.path.join(d, "fig_rebound_b4_pythagorean_chaos_head_to_head.pdf")
            fig.savefig(png_p, dpi=300)
            fig.savefig(pdf_p)
            print(f"Saved {png_p} and {pdf_p}")
    plt.close(fig)


if __name__ == "__main__":
    main()
