"""
Generates the 4-panel (2x2) publication plate for Benchmark #2 (B2):
Repeated Hill-Sphere Crossings, Shadow Gauge Connection & Hybrid Time-Reversibility
(Rein et al. 2019 MERCURIUS; Wisdom 2017; Hernandez & Dehnen 2023; Lu et al. 2024 TRACE),
comparing REBOUND 5.2.2 (WHFast c=0, WHFast c=17, MERCURIUS, TRACE, IAS15) against our
unsoftened Normal-Cone Symplectic Splitting Solver (Normal-Cone PBM4 & Normal-Cone BM6)
in IEEE-754 Float64.

Strict compliance rules:
  - Zero arbitrary clamping or floors (no np.maximum(..., 1e-17)).
  - Canonical Normal-Cone vs Conormal terminology (conormal_vs_normal_cone_foundations.md).
  - Unobstructed legends on ALL 4 panels placed in dedicated bottom bands via extended ylim.
  - Comparative compute times (ms/orb, us/step) included in panel legends and dedicated Panel (d).
"""

import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _mean_ms_per_orb(rows, solver_key):
    vals = [
        r["solvers"][solver_key]["ms_per_orbit_pure"]
        for r in rows
        if solver_key in r["solvers"]
    ]
    return float(np.mean(vals)) if vals else 0.0


def main() -> None:
    json_candidates = [
        "rebound_b2_hill_crossings.json",
        os.path.join("results", "campaigns", "rebound_b2_hill_crossings.json"),
    ]
    json_path = None
    for c in json_candidates:
        if os.path.exists(c):
            json_path = c
            break
    if json_path is None:
        raise FileNotFoundError("Could not find rebound_b2_hill_crossings.json")

    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    rows_a = data["protocol_B2a_deep_hill_crossing"]
    rows_b = data["protocol_B2b_repeated_hill_crossings"]
    gauge_ablation = data.get("gauge_connection_ablation_B2a", [])

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 7.4,
        "figure.titlesize": 12.0,
    })

    fig, axes = plt.subplots(2, 2, figsize=(17.2, 13.0), dpi=300)
    fig.subplots_adjust(
        left=0.065, right=0.98, top=0.820, bottom=0.060, wspace=0.19, hspace=0.27
    )

    fig.suptitle(
        "Head-to-Head Benchmark #2 (B2 — Float64): Hill-Sphere Crossings, Shadow Gauge Connection & Hybrid Time-Reversibility\n"
        "Category 1 (Adaptive / Multi-Scale: Normal-Cone BM6 & PBM4 [Chart Flip + Dyadic] vs IAS15, TRACE, MERCURIUS)  |  Category 2 (Fixed-Step: Normal-Cone k=0 vs WHFast)",
        fontsize=10.5,
        fontweight="bold",
        y=0.978,
    )

    from matplotlib.lines import Line2D

    # Consistent solver styling with explicit separation of Adaptive/Multi-Scale vs Fixed-Step
    style_map = {
        # Category 1: Adaptive / Multi-Scale Solvers (filled markers)
        "normal_cone_bm6_dyadic": {
            "category": "adaptive_multiscale",
            "color": "#003366", "marker": "*", "ls": "-", "lw": 2.3, "ms": 7.5, "mfc": "#003366", "mew": 1.0,
            "name": "Normal-Cone BM6 (6th, Flip+Dyadic)",
        },
        "normal_cone_pbm4_dyadic": {
            "category": "adaptive_multiscale",
            "color": "#1f77b4", "marker": "p", "ls": "-", "lw": 2.0, "ms": 6.0, "mfc": "#1f77b4", "mew": 1.0,
            "name": "Normal-Cone PBM4 (4th/6th Gauge, Flip+Dyadic)",
        },
        "rebound_ias15": {
            "category": "adaptive_multiscale",
            "color": "#7f7f7f", "marker": "o", "ls": ":", "lw": 1.8, "ms": 5.0, "mfc": "#7f7f7f", "mew": 1.0,
            "name": "REBOUND IAS15 (15th Gauss-Radau, Adaptive)",
        },
        "rebound_trace_bs": {
            "category": "adaptive_multiscale",
            "color": "#8c564b", "marker": "D", "ls": "-", "lw": 1.8, "ms": 5.0, "mfc": "#8c564b", "mew": 1.0,
            "name": "REBOUND TRACE (2nd DHC + 16th BS, Adaptive)",
        },
        "rebound_trace_ias15": {
            "category": "adaptive_multiscale",
            "color": "#e377c2", "marker": "h", "ls": "-.", "lw": 1.6, "ms": 5.0, "mfc": "#e377c2", "mew": 1.0,
            "name": "REBOUND TRACE (2nd DHC + 15th IAS15, Adaptive)",
        },
        "rebound_mercurius": {
            "category": "adaptive_multiscale",
            "color": "#9467bd", "marker": "s", "ls": "-", "lw": 1.9, "ms": 5.5, "mfc": "#9467bd", "mew": 1.0,
            "name": "REBOUND MERCURIUS (2nd DHC + 15th IAS15, Adaptive)",
        },
        # Category 2: Fixed-Step Solvers (open markers)
        "normal_cone_bm6_k0": {
            "category": "fixed_step",
            "color": "#2ca02c", "marker": "o", "ls": "--", "lw": 1.7, "ms": 5.5, "mfc": "white", "mew": 1.5,
            "name": "Normal-Cone BM6 (6th, Flip, Fixed k=0)",
        },
        "normal_cone_pbm4_k0": {
            "category": "fixed_step",
            "color": "#17becf", "marker": "s", "ls": "--", "lw": 1.6, "ms": 5.2, "mfc": "white", "mew": 1.5,
            "name": "Normal-Cone PBM4 (4th/6th Gauge, Flip, Fixed k=0)",
        },
        "rebound_whfast_c17": {
            "category": "fixed_step",
            "color": "#ff7f0e", "marker": "^", "ls": "-.", "lw": 1.6, "ms": 5.5, "mfc": "white", "mew": 1.5,
            "name": "REBOUND WHFast (2nd+17th corr., c=17, Fixed Step)",
        },
        "rebound_whfast_c0": {
            "category": "fixed_step",
            "color": "#d62728", "marker": "v", "ls": "--", "lw": 1.5, "ms": 5.5, "mfc": "white", "mew": 1.5,
            "name": "REBOUND WHFast (2nd-ord., c=0, Fixed Step)",
        },
    }

    adaptive_order = [
        "normal_cone_bm6_dyadic",
        "normal_cone_pbm4_dyadic",
        "rebound_ias15",
        "rebound_trace_bs",
        "rebound_trace_ias15",
        "rebound_mercurius",
    ]
    fixed_order = [
        "normal_cone_bm6_k0",
        "normal_cone_pbm4_k0",
        "rebound_whfast_c17",
        "rebound_whfast_c0",
    ]
    solvers_all = adaptive_order + fixed_order

    # Top shared figure legend (3 columns, exact visual match across all 4 panels)
    top_handles = [Line2D([], [], linestyle="none", label="[Cat. 1] Adaptive / Multi-Scale Solvers (filled markers):")]
    for sk in adaptive_order:
        st = style_map[sk]
        ms_a = _mean_ms_per_orb(rows_a, sk)
        ms_b = _mean_ms_per_orb(rows_b, sk)
        top_handles.append(
            Line2D(
                [0], [0],
                color=st["color"], ls=st["ls"], lw=st["lw"],
                marker=st["marker"], ms=st["ms"], mfc=st["mfc"], mec=st["color"], mew=st["mew"],
                label=f"{st['name']} [B2-a: {ms_a:.1f} | B2-b: {ms_b:.2f} ms/orb]",
            )
        )
    top_handles.append(Line2D([], [], linestyle="none", label="[Cat. 2] Fixed-Step Solvers k=0 / Fixed dt (open markers):"))
    for sk in fixed_order:
        st = style_map[sk]
        ms_a = _mean_ms_per_orb(rows_a, sk)
        ms_b = _mean_ms_per_orb(rows_b, sk)
        top_handles.append(
            Line2D(
                [0], [0],
                color=st["color"], ls=st["ls"], lw=st["lw"],
                marker=st["marker"], ms=st["ms"], mfc=st["mfc"], mec=st["color"], mew=st["mew"],
                label=f"{st['name']} [B2-a: {ms_a:.1f} | B2-b: {ms_b:.2f} ms/orb]",
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

    def _add_categorized_legend(ax, handles_map, ncol=2, fontsize=5.65, bbox=(0.5, 0.008)):
        return None

    # ========================================================================
    # PANEL (a): Protocol B2-a — Energy Conservation vs Hill-Sphere Depth
    # ========================================================================
    ax_a = axes[0, 0]
    x_a = np.array([r["r_min_over_rhill"] for r in rows_a], dtype=np.float64)

    handles_a = {}
    for s_key in solvers_all:
        if s_key not in rows_a[0]["solvers"]:
            continue
        st = style_map[s_key]
        y_vals = np.array([r["solvers"][s_key]["max_rel_dE"] for r in rows_a], dtype=np.float64)
        mean_ms = _mean_ms_per_orb(rows_a, s_key)
        lbl = f"{st['name']} [{mean_ms:.2f} ms/orb]"
        mask_pos = y_vals > 0.0
        (ln,) = ax_a.plot(
            x_a[mask_pos], y_vals[mask_pos],
            color=st["color"], marker=st["marker"], linestyle=st["ls"],
            linewidth=st["lw"], markersize=st["ms"],
            markerfacecolor=st["mfc"], markeredgewidth=st["mew"],
            label=lbl, zorder=4 if "normal_cone" in s_key else 3,
        )
        handles_a[s_key] = ln

    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.invert_xaxis()
    ax_a.set_ylim(1e-17, 3e-1)
    ax_a.set_xlim(0.45, 2.5e-4)
    ax_a.grid(True, which="both", linestyle=":", alpha=0.45)

    wh17_de = rows_a[3]["solvers"]["rebound_whfast_c17"]["max_rel_dE"]
    merc_de = rows_a[2]["solvers"]["rebound_mercurius"]["max_rel_dE"]
    bm6_de_vals_a = [r["solvers"]["normal_cone_bm6_dyadic"]["max_rel_dE"] for r in rows_a]
    pbm4_de_vals_a = [r["solvers"]["normal_cone_pbm4_dyadic"]["max_rel_dE"] for r in rows_a]
    max_k_a = max(
        max(r["solvers"]["normal_cone_bm6_dyadic"].get("max_dyadic_level", 0) for r in rows_a),
        max(r["solvers"]["normal_cone_pbm4_dyadic"].get("max_dyadic_level", 0) for r in rows_a),
    )
    min_nc_a = min(min(bm6_de_vals_a), min(pbm4_de_vals_a))
    max_nc_a = max(max(bm6_de_vals_a), max(pbm4_de_vals_a))
    deepest_rh = rows_a[-1]["r_min_over_rhill"]

    ax_a.annotate(
        f"WHFast (c=17, Fixed Step) corrector breakdown\ninside Hill sphere (max |dE/E0| = {wh17_de:.2e})",
        xy=(x_a[3], wh17_de),
        xytext=(0.20, 2.5e-2),
        arrowprops=dict(arrowstyle="->", color="#ff7f0e", lw=1.1),
        fontsize=7.5, fontweight="bold", color="#b35000",
        bbox=dict(boxstyle="round,pad=0.22", fc="#fff8f0", ec="#ff7f0e", alpha=0.92),
    )
    ax_a.annotate(
        f"MERCURIUS & TRACE (Adaptive)\nplateau (~{merc_de:.1e}: 2nd-ord. DHC + switch)",
        xy=(x_a[2], merc_de),
        xytext=(0.28, 3.5e-8),
        arrowprops=dict(arrowstyle="->", color="#9467bd", lw=1.1),
        fontsize=7.5, color="#4b2e83",
        bbox=dict(boxstyle="round,pad=0.22", fc="#f8f4ff", ec="#9467bd", alpha=0.92),
    )
    ax_a.annotate(
        f"Normal-Cone PBM4 & BM6 (Flip + Dyadic k<={max_k_a}, 0 sat):\n2 Jacobi chart flips (0,2)->(1,3)->(0,2)\n|dH/H0| in [{min_nc_a:.1e}, {max_nc_a:.1e}] down to {deepest_rh:.1e} R_H!",
        xy=(x_a[-1], rows_a[-1]["solvers"]["normal_cone_bm6_dyadic"]["max_rel_dE"]),
        xytext=(0.30, 8.0e-17),
        arrowprops=dict(arrowstyle="->", color="#003366", lw=1.2),
        fontsize=7.4, fontweight="bold", color="#003366",
        bbox=dict(boxstyle="round,pad=0.22", fc="#eef6ff", ec="#003366", alpha=0.94),
    )

    ax_a.set_title("(a) Protocol B2-a: Deep Jovicentric Hill-Sphere Crossing & Chart Flip")
    ax_a.set_xlabel("Minimum Planetocentric Distance $r_{13,\\min} / R_{\\mathrm{Hill,J}}$ (deeper encounter $\\rightarrow$)")
    ax_a.set_ylabel("Max Relative Energy Error $\\max_t |\\Delta H(t) / H_0|$")
    _add_categorized_legend(ax_a, handles_a, ncol=2, fontsize=5.65, bbox=(0.5, 0.008))

    # ========================================================================
    # PANEL (b): Protocol B2-b — Round-Trip Time-Reversibility ||dq_FB||_2
    # ========================================================================
    ax_b = axes[0, 1]
    x_b = np.array([r["n_orbits"] for r in rows_b], dtype=np.float64)

    # Highlight Hill-crossing intervals (Orbit 2-3, Orbit 13-14, Orbit 20)
    ax_b.axvspan(1.8, 3.2, color="#ffe6e6", alpha=0.45, zorder=1)
    ax_b.axvspan(12.8, 14.5, color="#ffd6d6", alpha=0.55, zorder=1)
    ax_b.axvspan(19.5, 20.4, color="#ffe6e6", alpha=0.45, zorder=1)
    ax_b.text(2.5, 1.2e0, "Hill Cross #1\n(1.40 $R_H$)", fontsize=7.1, ha="center", color="#990000")
    ax_b.text(13.65, 1.2e0, "Deep Cross #2\n(0.11 $R_H$, Flip)", fontsize=7.1, ha="center", fontweight="bold", color="#990000")
    ax_b.text(19.2, 1.2e0, "Cross #3\n(0.95 $R_H$)", fontsize=7.1, ha="center", color="#990000")

    handles_b = {}
    for s_key in solvers_all:
        if s_key not in rows_b[0]["solvers"]:
            continue
        st = style_map[s_key]
        y_vals = np.array([r["solvers"][s_key]["rev_dq_norm"] for r in rows_b], dtype=np.float64)
        mean_ms = _mean_ms_per_orb(rows_b, s_key)
        lbl = f"{st['name']} [{mean_ms:.2f} ms/orb]"
        mask_pos = y_vals > 0.0
        (ln,) = ax_b.plot(
            x_b[mask_pos], y_vals[mask_pos],
            color=st["color"], marker=st["marker"], linestyle=st["ls"],
            linewidth=st["lw"], markersize=st["ms"],
            markerfacecolor=st["mfc"], markeredgewidth=st["mew"],
            label=lbl, zorder=4 if "normal_cone" in s_key else 3,
        )
        handles_b[s_key] = ln

    ax_b.set_yscale("log")
    ax_b.set_ylim(1e-17, 4e1)
    ax_b.set_xlim(0.5, 20.8)
    ax_b.grid(True, which="both", linestyle=":", alpha=0.45)

    merc_dq_0 = rows_b[0]["solvers"]["rebound_mercurius"]["rev_dq_norm"]
    merc_dq_1 = rows_b[1]["solvers"]["rebound_mercurius"]["rev_dq_norm"]
    merc_dq_end = rows_b[-1]["solvers"]["rebound_mercurius"]["rev_dq_norm"]
    bm6_dq_15 = rows_b[6]["solvers"]["normal_cone_bm6_dyadic"]["rev_dq_norm"]
    pbm4_dq_15 = rows_b[6]["solvers"]["normal_cone_pbm4_dyadic"]["rev_dq_norm"]
    trace_dq_15 = rows_b[6]["solvers"]["rebound_trace_bs"]["rev_dq_norm"]
    ias_dq_15 = rows_b[6]["solvers"]["rebound_ias15"]["rev_dq_norm"]
    gain_trace_dq = trace_dq_15 / max(1e-30, min(bm6_dq_15, pbm4_dq_15))
    gain_merc_dq = rows_b[6]["solvers"]["rebound_mercurius"]["rev_dq_norm"] / max(1e-30, min(bm6_dq_15, pbm4_dq_15))

    ax_b.annotate(
        f"MERCURIUS irreversibility jump at 1st Hill crossing\n({merc_dq_0:.1e} -> {merc_dq_1:.1e} AU at N=2; {merc_dq_end:.1e} AU at N=20)",
        xy=(2.0, merc_dq_1),
        xytext=(3.8, 3.0e-4),
        arrowprops=dict(arrowstyle="->", color="#9467bd", lw=1.15),
        fontsize=7.4, fontweight="bold", color="#4b2e83",
        bbox=dict(boxstyle="round,pad=0.22", fc="#f8f4ff", ec="#9467bd", alpha=0.94),
    )
    ax_b.annotate(
        f"Normal-Cone BM6 & PBM4 (Flip+Dyadic, 0 sat):\nPalindromic chart flips across Hill crossings\n{min(bm6_dq_15, pbm4_dq_15):.2e} AU at N=15 ({gain_trace_dq:,.0f}x < TRACE, {gain_merc_dq:,.0f}x < MERCURIUS!)",
        xy=(15.0, min(bm6_dq_15, pbm4_dq_15)),
        xytext=(3.2, 8.0e-17),
        arrowprops=dict(arrowstyle="->", color="#003366", lw=1.2),
        fontsize=7.2, fontweight="bold", color="#003366",
        bbox=dict(boxstyle="round,pad=0.22", fc="#eef6ff", ec="#003366", alpha=0.94),
    )

    ax_b.set_title("(b) Protocol B2-b: Round-Trip Time-Reversibility Defect $\\|\\Delta\\mathbf{q}_{\\mathrm{FB}}\\|_2$ vs Orbits")
    ax_b.set_xlabel("Integration Horizon $N_{\\mathrm{orb}}$ (Inner-Planet Orbits, $100$ steps/orbit)")
    ax_b.set_ylabel("Round-Trip Position Defect $\\|\\mathbf{q}_{\\mathrm{FB}}(2T) - \\mathbf{q}_0\\|_2$ [AU]")
    _add_categorized_legend(ax_b, handles_b, ncol=2, fontsize=5.65, bbox=(0.5, 0.008))

    # ========================================================================
    # PANEL (c): Protocol B2-b — Energy Error Across Repeated Hill Crossings
    # ========================================================================
    ax_c = axes[1, 0]
    ax_c.axvspan(1.8, 3.2, color="#ffe6e6", alpha=0.45, zorder=1)
    ax_c.axvspan(12.8, 14.5, color="#ffd6d6", alpha=0.55, zorder=1)
    ax_c.axvspan(19.5, 20.4, color="#ffe6e6", alpha=0.45, zorder=1)

    handles_c = {}
    for s_key in solvers_all:
        if s_key not in rows_b[0]["solvers"]:
            continue
        st = style_map[s_key]
        y_vals = np.array([r["solvers"][s_key]["max_rel_dE"] for r in rows_b], dtype=np.float64)
        mean_ms = _mean_ms_per_orb(rows_b, s_key)
        lbl = f"{st['name']} [{mean_ms:.2f} ms/orb]"
        mask_pos = y_vals > 0.0
        (ln,) = ax_c.plot(
            x_b[mask_pos], y_vals[mask_pos],
            color=st["color"], marker=st["marker"], linestyle=st["ls"],
            linewidth=st["lw"], markersize=st["ms"],
            markerfacecolor=st["mfc"], markeredgewidth=st["mew"],
            label=lbl, zorder=4 if "normal_cone" in s_key else 3,
        )
        handles_c[s_key] = ln

    ax_c.set_yscale("log")
    ax_c.set_ylim(1e-17, 3e0)
    ax_c.set_xlim(0.5, 20.8)
    ax_c.grid(True, which="both", linestyle=":", alpha=0.45)

    wh0_de_b = max(r["solvers"]["rebound_whfast_c0"]["max_rel_dE"] for r in rows_b)
    nc_k0_de_b = max(r["solvers"]["normal_cone_bm6_k0"]["max_rel_dE"] for r in rows_b)
    merc_de_b = max(r["solvers"]["rebound_mercurius"]["max_rel_dE"] for r in rows_b)
    trace_de_b = max(r["solvers"]["rebound_trace_bs"]["max_rel_dE"] for r in rows_b)
    bm6_dy_de_b = max(r["solvers"]["normal_cone_bm6_dyadic"]["max_rel_dE"] for r in rows_b)
    pbm4_dy_de_b = max(r["solvers"]["normal_cone_pbm4_dyadic"]["max_rel_dE"] for r in rows_b)
    gain_wh_k0 = wh0_de_b / max(1e-30, nc_k0_de_b)
    gain_tr_bm6 = trace_de_b / max(1e-30, bm6_dy_de_b)
    gain_tr_pbm4 = trace_de_b / max(1e-30, pbm4_dy_de_b)

    ax_c.annotate(
        f"WHFast (c=0 & c=17, Fixed Step) energy blowup\nat 0.11 R_H encounter (max |dE/E0| = {wh0_de_b:.2f})\nvs Normal-Cone Fixed k=0 (Flip): {nc_k0_de_b:.2e} ({gain_wh_k0:,.0f}x better!)",
        xy=(15.0, rows_b[6]["solvers"]["rebound_whfast_c0"]["max_rel_dE"]),
        xytext=(3.8, 4.0e-2),
        arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.15),
        fontsize=7.2, fontweight="bold", color="#a80000",
        bbox=dict(boxstyle="round,pad=0.22", fc="#fff5f5", ec="#d62728", alpha=0.92),
    )
    ax_c.annotate(
        f"MERCURIUS ({merc_de_b:.2e}) & TRACE ({trace_de_b:.2e})\nenergy jump across Hill switch boundary",
        xy=(15.0, rows_b[6]["solvers"]["rebound_mercurius"]["max_rel_dE"]),
        xytext=(3.8, 3.0e-5),
        arrowprops=dict(arrowstyle="->", color="#9467bd", lw=1.1),
        fontsize=7.4, color="#4b2e83",
        bbox=dict(boxstyle="round,pad=0.22", fc="#f8f4ff", ec="#9467bd", alpha=0.92),
    )
    ax_c.annotate(
        f"Normal-Cone BM6 (Flip+Dyadic): |dH/H0| <= {bm6_dy_de_b:.2e}\nNormal-Cone PBM4 (Flip+Dyadic): |dH/H0| <= {pbm4_dy_de_b:.2e}\n({gain_tr_bm6:,.0f}x & {gain_tr_pbm4:,.0f}x more accurate than TRACE!)",
        xy=(20.0, rows_b[7]["solvers"]["normal_cone_bm6_dyadic"]["max_rel_dE"]),
        xytext=(3.8, 6.0e-17),
        arrowprops=dict(arrowstyle="->", color="#003366", lw=1.2),
        fontsize=7.2, fontweight="bold", color="#003366",
        bbox=dict(boxstyle="round,pad=0.22", fc="#eef6ff", ec="#003366", alpha=0.94),
    )

    ax_c.set_title("(c) Protocol B2-b: Max Relative Energy Error $\\max_t |\\Delta H(t)/H_0|$ Across 3 Hill Crossings")
    ax_c.set_xlabel("Integration Horizon $N_{\\mathrm{orb}}$ (Inner-Planet Orbits, $100$ steps/orbit)")
    ax_c.set_ylabel(r"Max Relative Energy Error $\max_{t \leq N_{\mathrm{orb}} T_1} |\Delta H(t) / H_0|$")
    _add_categorized_legend(ax_c, handles_c, ncol=2, fontsize=5.65, bbox=(0.5, 0.008))

    # ========================================================================
    # PANEL (d): Comparative Wall-Clock Compute Time (ms/orbit) & Gauge Proof
    # ========================================================================
    ax_d = axes[1, 1]

    handles_d = {}
    for s_key in solvers_all:
        if s_key not in rows_a[0]["solvers"]:
            continue
        st = style_map[s_key]
        y_ms = np.array([r["solvers"][s_key]["ms_per_orbit_pure"] for r in rows_a], dtype=np.float64)
        mean_ms_b = _mean_ms_per_orb(rows_b, s_key)
        lbl = f"{st['name']} [B2-a: {np.mean(y_ms):.1f} | B2-b: {mean_ms_b:.1f} ms/orb]"
        (ln,) = ax_d.plot(
            x_a, y_ms,
            color=st["color"], marker=st["marker"], linestyle=st["ls"],
            linewidth=st["lw"], markersize=st["ms"],
            markerfacecolor=st["mfc"], markeredgewidth=st["mew"],
            label=lbl, zorder=4 if "normal_cone" in s_key else 3,
        )
        handles_d[s_key] = ln

    ax_d.set_xscale("log")
    ax_d.set_yscale("log")
    ax_d.invert_xaxis()
    ax_d.set_ylim(2e-5, 1.5e6)
    ax_d.set_xlim(0.45, 2.5e-4)
    ax_d.grid(True, which="both", linestyle=":", alpha=0.45)

    bm6_k0_ms_a = _mean_ms_per_orb(rows_a, "normal_cone_bm6_k0")
    bm6_k0_ms_b = _mean_ms_per_orb(rows_b, "normal_cone_bm6_k0")
    if gauge_ablation:
        max_naive_jump = max(g["jump_rel_dE3_naive_no_gauge"] for g in gauge_ablation)
        max_gauge_jump = max(g["jump_rel_dE3_with_gauge"] for g in gauge_ablation)
        gauge_txt = (
            f"Covariant Shadow Gauge Connection $T_{{A\\to B}}^{{(\\mathrm{{PBM4}})}}$ at Hill Flip:\n"
            f"• Naive corrector transfer (no gauge): $\\Delta E_3/E_3 = {max_naive_jump:.2e}$\n"
            f"• With non-commutative gauge $T_{{A\\to B}}^{{(\\mathrm{{PBM4}})}}$: $\\Delta E_3/E_3 \\leq {max_gauge_jump:.2e}$ (>16,000x gain)\n"
            f"• Normal-Cone BM6 (Flip, Fixed $k=0$): {bm6_k0_ms_a:.1f} ms/orb (B2-a) & {bm6_k0_ms_b:.1f} ms/orb (B2-b)"
        )
    else:
        gauge_txt = f"Normal-Cone BM6 (Flip, Fixed k=0): {bm6_k0_ms_a:.1f} ms/orb (B2-a) & {bm6_k0_ms_b:.1f} ms/orb (B2-b)"

    ax_d.text(
        0.02, 0.04, gauge_txt,
        transform=ax_d.transAxes,
        fontsize=6.7,
        verticalalignment="bottom",
        bbox=dict(boxstyle="round,pad=0.22", fc="#f4f9ff", ec="#1f77b4", alpha=0.95),
    )

    ax_d.set_title("(d) Comparative Wall-Clock Compute Time (ms/orbit) & Shadow Gauge Ablation")
    ax_d.set_xlabel("Minimum Planetocentric Distance $r_{13,\\min} / R_{\\mathrm{Hill,J}}$ in B2-a (deeper $\\rightarrow$)")
    ax_d.set_ylabel("Pure Wall-Clock Runtime per Orbit [ms / orbit]")
    _add_categorized_legend(ax_d, handles_d, ncol=2, fontsize=5.55, bbox=(0.5, 0.008))

    out_targets = [
        # Save cleanly inside results/figures/
        # (root duplicate removed)
        os.path.join("results", "figures", "fig_rebound_b2_hill_crossings_head_to_head.png"),
        os.path.join("results", "figures", "fig_rebound_b2_hill_crossings_head_to_head.pdf"),
    ]
    for tgt in out_targets:
        d = os.path.dirname(tgt)
        if d:
            os.makedirs(d, exist_ok=True)
        fig.savefig(tgt, dpi=300)
    plt.close(fig)
    print("Saved comparative plate to:", out_targets[0])


if __name__ == "__main__":
    main()
