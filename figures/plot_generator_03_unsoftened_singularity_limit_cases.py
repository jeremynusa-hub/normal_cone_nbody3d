#!/usr/bin/env python3
r"""
Pure Three-Volet Plotter for Benchmark Generator 03 (Unsoftened Singularity Traversals & Limit Cases L1-L4, §7.1.3).

Reads ONLY:
    results/campaigns/generator_03_unsoftened_singularity_limit_cases.json

Generates 3 standalone Laskar-standard plates (PDF + PNG) in:
    - results/figures/
    - Artifact directory (for visual inspection)
    - Legacy dashboard path fig_ch07_unsoftened_limit_cases_dashboard.png
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
from matplotlib.ticker import LogLocator, NullFormatter

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = WORKSPACE_ROOT / "results" / "campaigns" / "generator_03_unsoftened_singularity_limit_cases.json"
CH07_FIG_DIR = WORKSPACE_ROOT / "results" / "figures"
ARTIFACT_DIR = WORKSPACE_ROOT / "results" / "figures"


def _configure_matplotlib() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.2,
        "figure.titlesize": 12.5,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": "--",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.22,
    })


def _pos_or_nan(arr: Any) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    return np.where(a > 0.0, a, np.nan)


def _clean_log_axes(ax: plt.Axes, clean_x: bool = True, clean_y: bool = True) -> None:
    if clean_x:
        ax.xaxis.set_minor_formatter(NullFormatter())
    if clean_y:
        ax.yaxis.set_minor_formatter(NullFormatter())


def plot_volet_A(data: Dict[str, Any]) -> None:
    l1 = data["L1_eccentricity_homotopy"]
    l3 = data["L3_burrau_and_quadruple_d4"]
    l4 = data["L4_kozai_lidov"]

    fig, axes = plt.subplots(2, 2, figsize=(14.6, 10.4), dpi=300)
    fig.subplots_adjust(wspace=0.28, hspace=0.35)
    fig.suptitle(
        "Computational Conditioning, Dyadic Adaptivity & Dual-Precision Singularity Scaling",
        fontweight="bold",
        y=0.98,
    )

    # Panel (a): L1 Pericenter Conditioning Across 1 - e in [1e-12, 1e-1] & e = 1 (Pure N=2, 1e-12 on Left)
    ax = axes[0, 0]
    ome = np.array(l1["one_minus_e"], dtype=float)
    verlet = _pos_or_nan(l1["cartesian_verlet_rel_dE"])
    legacy_uncond = _pos_or_nan(l1["legacy_unconditioned_f64_rel_dE"])
    leaf64 = _pos_or_nan(l1["float64_pbm4_rel_dE"])
    leaf128 = _pos_or_nan(l1["float128_pbm4_rel_dE"])

    ax.loglog(ome, verlet, "o--", color="#DC2626", lw=1.7, ms=5, label=r"Unregularized Cartesian Verlet ($\propto (1-e)^{-5/2}$)")
    ax.loglog(ome, legacy_uncond, "d:", color="#7C3AED", lw=1.6, ms=4.8, label=r"Unconditioned Stumpff Subtraction ($\propto (1-e)^{-2}\epsilon_{\mathrm{mach}}$)")
    ax.loglog(ome, leaf64, "^-", color="#1D4ED8", lw=2.0, ms=5.5, label=r"Conditioned Sundman–Stumpff Leaf (Float64, $e=1$: " + f"{l1['rectilinear_e1_f64_rel_dE']:.1e})")
    ax.loglog(ome, leaf128, "h-", color="#059669", lw=2.0, ms=5.5, label=r"Conditioned Sundman–Stumpff Leaf (Float128, $e=1$: " + f"{l1['rectilinear_e1_f128_rel_dE']:.1e})")

    ax.set_xlim(5e-13, 2e-1)
    ax.set_ylim(1e-38, 1e6)
    ax.set_xlabel(r"Pericenter Proximity Parameter $1 - e = r_{\min} / a$")
    ax.set_ylabel(r"Full-Orbit Relative Energy Error $\max_{t\in[0,T_{\mathrm{orb}}]} |\Delta E|/|E_0|$")
    ax.set_title(r"(a) L1: Analytic Pericenter Conditioning ($N=2, 1.0\,T_{\mathrm{orb}}$)")
    ax.legend(loc="center right", bbox_to_anchor=(0.98, 0.36), framealpha=0.94, fontsize=7.6)
    _clean_log_axes(ax, True, True)

    # Panel (b): Symmetric Dyadic Level k_sym(t) & Leaf/Apex Decoupling
    ax = axes[0, 1]
    t_b = np.array(l3["burrau"]["times"], dtype=float)
    k_b = np.array(l3["burrau"]["dyadic_levels"], dtype=float)[-len(t_b):]
    t_kl = np.array(l4["times"], dtype=float) / 25.0  # normalized to [0, 12] abscissa for joint comparison
    k_kl = np.array(l4["dyadic_levels"], dtype=float)[-len(t_kl):]

    ax.step(t_b, k_b, where="post", color="#1D4ED8", lw=1.7, label=r"L3a Burrau 3-4-5 Chaotic Scattering ($k_{\max} = " + f"{int(np.max(k_b))}$)")
    ax.step(t_kl, k_kl, where="post", color="#D97706", lw=1.5, alpha=0.85, label=r"L4 Kozai–Lidov ($t/25$, $e_{\max}=0.9989$, $k_{\max} = " + f"{int(np.max(k_kl))}$)")
    ax.axhline(0.0, color="#059669", lw=2.4, ls="-", label=r"L1 Kepler Homotopy ($1-e \in [10^{-1}, 10^{-12}]$): $k_{\mathrm{sym}}(t) \equiv 0$")

    ax.set_xlim(0.0, 12.0)
    ax.set_ylim(-0.6, max(float(np.max(k_b)), float(np.max(k_kl))) + 1.5)
    ax.set_xlabel(r"Physical Time $t$ (Burrau) / Rescaled Time $t/25$ (Kozai–Lidov)")
    ax.set_ylabel(r"Symmetric Dyadic Refinement Level $k_{\mathrm{sym}}(t)$")
    ax.set_title(r"(b) Leaf/Apex Decoupling ($k_{\mathrm{sym}}\equiv 0$ on L1) & Multirate Adaptivity")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=8.0)

    # Panel (c): L3b D_4 Quadruple Collision Reversibility Across Singularity (t_max = 2.50 > t_coll = 1.1353)
    ax = axes[1, 0]
    qd = l3["quadruple_d4"]
    t_fb = np.array([x["t_max"] for x in qd["float64_fb"]], dtype=float)
    fb_dE64 = _pos_or_nan([x["delta_E_rel"] for x in qd["float64_fb"]])
    fb_dq64 = _pos_or_nan([x["delta_q_norm"] for x in qd["float64_fb"]])
    fb_dE128 = _pos_or_nan([x["delta_E_rel"] for x in qd["float128_fb"]])
    fb_dq128 = _pos_or_nan([x["delta_q_norm"] for x in qd["float128_fb"]])
    t_coll_d4 = qd["t_coll_analytical"]

    ax.semilogy(t_fb, fb_dE64, "o-", color="#1D4ED8", lw=1.8, ms=6, label=r"Float64 Round-Trip $|\Delta E_{\mathrm{FB}}|/|E_0|$")
    ax.semilogy(t_fb, fb_dq64, "s--", color="#0891B2", lw=1.6, ms=5.5, label=r"Float64 Round-Trip $\|\Delta\mathbf{q}_{\mathrm{FB}}\|_2$")
    ax.semilogy(t_fb, fb_dE128, "^-", color="#059669", lw=1.8, ms=6, label=r"Float128 Round-Trip $|\Delta E_{\mathrm{FB}}|/|E_0|$")
    ax.semilogy(t_fb, fb_dq128, "h--", color="#D97706", lw=1.6, ms=5.5, label=r"Float128 Round-Trip $\|\Delta\mathbf{q}_{\mathrm{FB}}\|_2$")
    ax.axvline(t_coll_d4, color="#DC2626", lw=1.8, ls=":", label=rf"Exact $D_4$ Collision $t_{{\mathrm{{coll}}}}^{{(D_4)}} = {t_coll_d4:.4f}$")
    ax.axvspan(t_coll_d4, 2.60, color="#FEE2E2", alpha=0.35, label=r"Post-Collision Rebound Regime ($R \to 0^+$ crossed)")

    ax.set_xlim(0.50, 2.60)
    ax.set_ylim(1e-38, 1e-10)
    ax.set_xlabel(r"Forward Integration Horizon $T_{\max} = N \Delta t$ ($\Delta t = 0.025$)")
    ax.set_ylabel(r"Round-Trip Reversibility Error ($\Phi_{-\Delta t}^N \circ \Phi_{\Delta t}^N$)")
    ax.set_title(r"(c) L3b: Reversibility Across $D_4$ Quadruple Collision ($T_{\max}=2.50 > t_{\mathrm{coll}}^{(D_4)}$)")
    ax.legend(loc="center left", framealpha=0.92, fontsize=7.6)
    _clean_log_axes(ax, False, True)

    # Panel (d): Dual-Precision (Float64 & Float128) Energy Evolution on L3a Burrau & L4 Kozai-Lidov
    ax = axes[1, 1]
    dE_b64 = _pos_or_nan(l3["burrau"]["rel_energy_errors"])
    dE_b128 = _pos_or_nan(l3["burrau"]["rel_energy_errors_f128"])
    t_kl64 = np.array(l4["times"], dtype=float) / 25.0
    dE_kl64 = _pos_or_nan(l4["rel_energy_errors"])
    t_kl128 = np.array(l4["times_f128"], dtype=float) / 25.0
    dE_kl128 = _pos_or_nan(l4["rel_energy_errors_f128"])

    ax.semilogy(t_b, dE_b64, "-", color="#1D4ED8", lw=1.5, label=rf"L3a Burrau Float64 ($\max = {np.nanmax(dE_b64):.2e}$)")
    ax.semilogy(t_b, dE_b128, "--", color="#0891B2", lw=1.5, label=rf"L3a Burrau Float128 ($\max = {np.nanmax(dE_b128):.2e}$)")
    ax.semilogy(t_kl64, dE_kl64, "-", color="#D97706", lw=1.6, label=rf"L4 Kozai–Lidov Float64 ($t/25$, $\max = {np.nanmax(dE_kl64):.2e}$)")
    ax.semilogy(t_kl128, dE_kl128, "--", color="#059669", lw=1.6, label=rf"L4 Kozai–Lidov Float128 ($t/25$, $\max = {np.nanmax(dE_kl128):.2e}$)")

    ax.set_xlim(0.0, 12.0)
    ax.set_ylim(1e-20, 1e-6)
    ax.set_xlabel(r"Physical Time $t$ (Burrau) / Rescaled Time $t/25$ (Kozai–Lidov)")
    ax.set_ylabel(r"Relative Hamiltonian Error $|\Delta H(t)| / |H_0|$")
    ax.set_title(r"(d) L3a & L4: Certified Dual-Precision (Float64 & Float128) Energy Boundedness")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=7.8)
    _clean_log_axes(ax, False, True)

    for out_dir in (CH07_FIG_DIR, ARTIFACT_DIR):
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "gen03_plate_A_computational.pdf")
        fig.savefig(out_dir / "gen03_plate_A_computational.png", dpi=300)
    plt.close(fig)


def plot_volet_B(data: Dict[str, Any]) -> None:
    l2 = data["L2_mcgehee_and_corrector"]
    l3 = data["L3_burrau_and_quadruple_d4"]
    l4 = data["L4_kozai_lidov"]

    fig, axes = plt.subplots(2, 2, figsize=(14.8, 10.8), dpi=300)
    fig.subplots_adjust(top=0.91, bottom=0.10, left=0.08, right=0.96, wspace=0.28, hspace=0.36)
    fig.suptitle(
        r"Conormal Bounce Invariants, $\mathrm{PBM}_4$ Free-Fall Regularity & Noether vs. Delaunay Actions",
        fontweight="bold",
        y=0.97,
    )

    # Panel (a): L2.1 Energy Conservation Across d_0 = 0 Conormal Bounce (N=2 Sundman & N=3 Lagrange)
    ax = axes[0, 0]
    tr_n2 = l2["traversal"]
    tr_lag = l2["lagrange_n3_mcgehee"]
    # Normalize time by full cycle 2 * t_coll so both N=2 and N=3 span [0, 1] and cross d_0 = 0 at 0.50
    t_n2 = np.array(tr_n2["times"][1:], dtype=float) / (2.0 * tr_n2["t_coll_analytical"])
    dE_n2_64 = _pos_or_nan(tr_n2["rel_energy_errors"][1:])
    dE_n2_128 = _pos_or_nan(tr_n2["rel_energy_errors_f128"][1:])
    t_lag = np.array(tr_lag["times"][1:], dtype=float) / (2.0 * tr_lag["t_coll_analytical"])
    dE_lag_64 = _pos_or_nan(tr_lag["rel_energy_errors_f64"][1:])
    dE_lag_128 = _pos_or_nan(tr_lag["rel_energy_errors_f128"][1:])

    ax.semilogy(t_n2, dE_n2_64, "-", color="#1D4ED8", lw=1.8, label=rf"$N=2$ Sundman Bounce (Float64, $\max={tr_n2['max_dE_f64']:.1e}$)")
    ax.semilogy(t_lag, dE_lag_64, "--", color="#D97706", lw=1.8, label=rf"$N=3$ Lagrange McGehee (Float64, $\max={tr_lag['max_dE_f64']:.1e}$)")
    ax.semilogy(t_n2, dE_n2_128, "-", color="#059669", lw=1.8, label=rf"$N=2$ Sundman Bounce (Float128, $\max={tr_n2['max_dE_f128']:.1e}$)")
    ax.semilogy(t_lag, dE_lag_128, ":", color="#0891B2", lw=2.0, label=rf"$N=3$ Lagrange McGehee (Float128, $\max={tr_lag['max_dE_f128']:.1e}$)")
    ax.axvline(0.50, color="#DC2626", lw=1.5, ls=":", alpha=0.85, label=rf"Exact Collision $t = t_{{\mathrm{{coll}}}}$ ($N=2: {tr_n2['t_coll_analytical']:.3f}$, $L_3: {tr_lag['t_coll_analytical']:.3f}$)")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(1e-37, 1e-10)
    ax.set_xlabel(r"Normalized Collision Cycle Time $t / (2\,t_{\mathrm{coll}})$ (Bounce at $0.50$)")
    ax.set_ylabel(r"Relative Energy Error $|\Delta E(t)| / |E_0|$")
    ax.set_title(r"(a) L2.1: Exact Energy Preservation Across $N=2$ & $N=3$ Collisions")
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.5)
    _clean_log_axes(ax, False, True)

    # Panel (b): L2.2 Asymptotic Regularity of PBM4 (Z_4^{BM4}) Corrector Along Free-Fall w(r) ~ r^{-1/2}
    ax = axes[0, 1]
    sc = l2["c4_kick_scaling"]
    r_arr = np.array(sc["separations"], dtype=float)
    k4_arr = np.array(sc["k4_generator_abs"], dtype=float)
    dr_pbm4 = np.array(sc["pbm4_position_kick_norms"], dtype=float)
    dw_pbm4 = np.array(sc["pbm4_momentum_kick_norms"], dtype=float)
    dw_rel_pbm4 = np.array(sc["pbm4_relative_momentum_kick"], dtype=float)

    ax.loglog(r_arr, dr_pbm4, "o-", color="#1D4ED8", lw=1.8, ms=4.5, label=rf"Corrector Position Shift $\|\Delta\mathbf{{r}}_{{\mathrm{{corr}}}}\| \propto r^{{+{sc['alpha_dr_pbm4']:.2f}}}$ ($\mathcal{{O}}(r^4) \to 0$)")
    ax.loglog(r_arr, k4_arr, "s-", color="#059669", lw=1.7, ms=4.5, label=rf"Corrector Generating Potential $|K_4(\mathbf{{r}},\mathbf{{w}})| \propto r^{{+{sc['alpha_k4']:.2f}}}$ ($\mathcal{{O}}(r^{{7/2}}) \to 0$)")
    ax.loglog(r_arr, dw_rel_pbm4, "^--", color="#0891B2", lw=1.6, ms=4.5, label=rf"Relative Momentum Shift $\|\Delta\mathbf{{w}}_{{\mathrm{{corr}}}}\|/\|\mathbf{{w}}\| \propto r^{{+{sc['alpha_dw_rel_pbm4']:.2f}}}$ ($\mathcal{{O}}(r^3) \to 0$)")
    ax.loglog(r_arr, dw_pbm4, "h-", color="#D97706", lw=1.7, ms=4.5, label=rf"Absolute Momentum Shift $\|\Delta\mathbf{{w}}_{{\mathrm{{corr}}}}\| \propto r^{{+{sc['alpha_dw_pbm4']:.2f}}}$ ($\mathcal{{O}}(r^{{5/2}}) \to 0$)")

    ax.set_xlim(5e-13, 2.0)
    ax.set_ylim(1e-65, 1e-10)
    ax.set_xlabel(r"Inner Binary Separation $r = \|\mathbf{q}_1 - \mathbf{q}_0\|$ Along Free-Fall $w(r) \propto r^{-1/2}$ ($h_{\mathrm{sub}}(r) \propto r^{3/2}$)")
    ax.set_ylabel(r"Corrector Potential & Phase-Space Shifts ($h_0 = 0.01$)")
    ax.set_title(r"(b) L2: Asymptotic Free-Fall Regularity of $\mathrm{PBM}_4$ ($r \in [10^{-12}, 1]$ AU)")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=7.2)
    _clean_log_axes(ax, True, True)

    # Panel (c): L3a & L3b Lie-Bargmann Invariants & D_4 Quadruple Bounce Conservation
    ax = axes[1, 0]
    lb64 = l3["burrau"]["lie_bargmann_drift_f64"]
    lb128 = l3["burrau"]["lie_bargmann_drift_f128"]
    qd = l3["quadruple_d4"]

    categories = [
        r"Burrau $\|\Delta\mathbf{P}\|_2$",
        r"Burrau $\|\Delta\mathbf{L}\|_2$",
        r"Burrau $\|\Delta\mathbf{G}\|_2$",
        r"$D_4$ Bounce $|\Delta E|/|E_0|$",
        r"$D_4$ Bounce FB $\|\Delta\mathbf{q}\|$",
    ]
    vals64 = [
        float(lb64["abs_P"]),
        float(lb64["abs_L"]),
        float(lb64["abs_G"]),
        float(qd["max_dE_f64"]),
        float(qd["float64_fb"][-1]["delta_q_norm"]),
    ]
    vals128 = [
        float(lb128["abs_P"]),
        float(lb128["abs_L"]),
        float(lb128["abs_G"]),
        float(qd["max_dE_f128"]),
        float(qd["float128_fb"][-1]["delta_q_norm"]),
    ]

    x_pos = np.arange(len(categories))
    width = 0.35
    ax.bar(x_pos - width / 2, vals64, width, color="#1D4ED8", alpha=0.88, label=r"Float64 (IEEE 754 Double Precision)")
    ax.bar(x_pos + width / 2, vals128, width, color="#059669", alpha=0.88, label=r"Float128 (35-Digit Extended Precision)")
    ax.set_yscale("log")
    ax.set_ylim(1e-37, 1e-8)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(categories, rotation=10, ha="right")
    ax.set_ylabel(r"Maximum Residual Across Full Horizon")
    ax.set_title(r"(c) L3a & L3b: Lie–Bargmann Generators & $D_4$ Collision Invariants")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=8.0)
    _clean_log_axes(ax, False, True)

    # Panel (d): L4 Kozai-Lidov: Exact Noether L_{z,tot} vs Semisecular Inner Delaunay Action L_{z,in}
    ax = axes[1, 1]
    t_kl64 = np.array(l4["times"], dtype=float)
    dLz_tot64 = _pos_or_nan(l4["total_Lz_drifts_f64"])
    t_kl_del = np.array(l4["times"][1:], dtype=float)
    delaunay_dev = _pos_or_nan(l4["kozai_integral_deviations"][1:])
    t_kl128 = np.array(l4["times_f128"], dtype=float)
    dLz_tot128 = _pos_or_nan(l4["total_Lz_drifts_f128"])

    ax.semilogy(
        t_kl_del,
        delaunay_dev,
        "-",
        color="#D97706",
        lw=1.8,
        label=rf"Semisecular Inner Delaunay $|\Delta(\sqrt{{1-e_{{\mathrm{{in}}}}^2}}\cos i_{{\mathrm{{in}}}})|$ ($\max = {l4['max_kozai_int_dev']:.2e}$)",
    )
    ax.semilogy(
        t_kl64,
        dLz_tot64,
        "-",
        color="#1D4ED8",
        lw=1.8,
        label=rf"Exact Noether Total $|L_{{z,\mathrm{{tot}}}}(t) - L_{{z,\mathrm{{tot}}}}(0)|$ (Float64, $\max = {l4['max_total_Lz_drift_f64']:.2e}$)",
    )
    ax.semilogy(
        t_kl128,
        dLz_tot128,
        "--",
        color="#059669",
        lw=1.8,
        label=rf"Exact Noether Total $|L_{{z,\mathrm{{tot}}}}(t) - L_{{z,\mathrm{{tot}}}}(0)|$ (Float128, $\max = {l4['max_total_Lz_drift_f128']:.2e}$)",
    )

    ax.set_xlim(0.0, 305.0)
    ax.set_ylim(1e-37, 1e1)
    ax.set_xlabel(r"Physical Time $t$ ($a_{\mathrm{out}}/a_{\mathrm{in}} = 4.0, i_0 = 84.0^\circ$)")
    ax.set_ylabel(r"Vertical Angular Momentum / Action Variation")
    ax.set_title(r"(d) L4: Exact Noether Conservation of $L_{z,\mathrm{tot}}$ vs. Semisecular Action")
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.5)
    _clean_log_axes(ax, False, True)

    for out_dir in (CH07_FIG_DIR, ARTIFACT_DIR):
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "gen03_plate_B_algebraic_invariants.pdf")
        fig.savefig(out_dir / "gen03_plate_B_algebraic_invariants.png", dpi=300)
    plt.close(fig)


def plot_volet_C(data: Dict[str, Any]) -> None:
    l1 = data["L1_eccentricity_homotopy"]
    l2 = data["L2_mcgehee_and_corrector"]
    l3 = data["L3_burrau_and_quadruple_d4"]
    l4 = data["L4_kozai_lidov"]

    fig, axes = plt.subplots(2, 2, figsize=(14.6, 10.4), dpi=300)
    fig.subplots_adjust(wspace=0.28, hspace=0.35)
    fig.suptitle(
        r"Physical Trajectories, McGehee Exceptional Divisor Portraits & Kozai–Lidov Cycles",
        fontweight="bold",
        y=0.98,
    )

    # Panel (a): L1 Near-Parabolic & Rectilinear Keplerian Homotopy Distance r(t) down to 1e-12 AU (Symlog around t = 0.5 T_orb)
    ax = axes[0, 0]
    palette = ["#64748B", "#0891B2", "#059669", "#D97706", "#1D4ED8"]
    for idx, (k_str, orb) in enumerate(l1["sample_orbits"].items()):
        if "delta_u" in orb:
            du_arr = np.array(orb["delta_u"], dtype=float)
        else:
            du_arr = np.array(orb["times"], dtype=float) / l1["T_orb"] - 0.5
        r_arr = _pos_or_nan(orb["r_norm"])
        col = palette[idx % len(palette)]
        ax.plot(du_arr, r_arr, "-", color=col, lw=1.8, label=rf"$1 - e = {k_str}$ ($r_{{\min}} = {np.nanmin(r_arr):.1e}$ AU)")

    du_step64 = 0.5 / 129.0
    ax.axvspan(-du_step64, +du_step64, color="#FEF3C7", alpha=0.35, label=r"Single Macro-Step 64 ($\Delta\tilde{t} \in [-\frac{1}{258}, +\frac{1}{258}]$)")
    ax.axvline(0.0, color="#DC2626", lw=1.3, ls=":", label=r"Exact Pericenter Passage ($t = 0.5\,T_{\mathrm{orb}}$)")
    ax.set_xscale("symlog", linthresh=1e-20, linscale=0.35)
    ax.set_yscale("log")
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(1e-13, 5.0)
    xt = [-1e-1, -1e-6, -1e-12, -1e-18, 0.0, 1e-18, 1e-12, 1e-6, 1e-1]
    ax.set_xticks(xt)
    ax.set_xticklabels([r"$-10^{-1}$", r"$-10^{-6}$", r"$-10^{-12}$", r"$-10^{-18}$", r"$0$", r"$+10^{-18}$", r"$+10^{-12}$", r"$+10^{-6}$", r"$+10^{-1}$"], fontsize=7.5)
    ax.set_xlabel(r"Pericenter-Centered Time $(t - 0.5\,T_{\mathrm{orb}}) / T_{\mathrm{orb}}$ [symlog]")
    ax.set_ylabel(r"Binary Separation $r(t) = \|\mathbf{q}_1(t) - \mathbf{q}_0(t)\|$ [AU]")
    ax.set_title(r"(a) L1: Symlog Pericenter Traversal Down to $r_{\min} = 10^{-12}$ AU")
    ax.legend(loc="lower left", framealpha=0.92, fontsize=7.2)
    _clean_log_axes(ax, False, True)

    # Panel (b): L2.1 McGehee Blown-Up Phase Portrait (rho, v_rad = drho/dt * sqrt(rho)) Across rho = 0
    ax = axes[0, 1]
    tr_n2 = l2["traversal"]
    tr_lag = l2["lagrange_n3_mcgehee"]
    r_n2 = np.array(tr_n2["separations"], dtype=float)
    v_n2 = np.array(tr_n2["v_rad_mcgehee"], dtype=float)
    r_lag = np.array(tr_lag["hyper_radius_R"], dtype=float)
    v_lag = np.array(tr_lag["v_rad_mcgehee"], dtype=float)

    # Plot smooth analytical McGehee manifold curves alongside discrete integration points
    r_grid_n2 = np.linspace(0.0, 2.0, 300)
    v_th_n2 = np.sqrt(np.maximum(4.0 - 2.0 * r_grid_n2, 0.0))
    ax.plot(r_grid_n2, -v_th_n2, "--", color="#93C5FD", lw=1.5)
    ax.plot(r_grid_n2, +v_th_n2, "--", color="#93C5FD", lw=1.5)
    ax.plot(r_n2, v_n2, "o-", color="#1D4ED8", lw=1.6, ms=3.5, label=r"$N=2$ Sundman–Levi-Civita ($\mathcal{B}_C: -2.0 \mapsto +2.0$ at $r=0$)")

    gm_eff_lag = 1.0 / np.sqrt(3.0)
    r_grid_lag = np.linspace(0.0, tr_lag["R0"], 300)
    v_th_lag = np.sqrt(np.maximum(2.0 * gm_eff_lag * (1.0 - r_grid_lag / tr_lag["R0"]), 0.0))
    ax.plot(r_grid_lag, -v_th_lag, "--", color="#FCD34D", lw=1.5)
    ax.plot(r_grid_lag, +v_th_lag, "--", color="#FCD34D", lw=1.5)
    ax.plot(r_lag, v_lag, "s-", color="#D97706", lw=1.6, ms=3.5, label=r"$N=3$ Equilateral Lagrange McGehee ($\mathcal{B}_C: -1.075 \mapsto +1.075$ at $R=0$)")

    ax.annotate(
        "",
        xy=(0.02, +1.95),
        xytext=(0.02, -1.95),
        arrowprops=dict(arrowstyle="->", color="#DC2626", lw=2.2, ls="-"),
    )
    ax.text(
        0.14,
        1.48,
        r"Conormal Bounce Involution $\mathcal{B}_C$ at $\rho = 0$",
        color="#DC2626",
        fontsize=8.3,
        fontweight="bold",
        va="center",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#DC2626", alpha=0.92),
    )

    ax.axhline(0.0, color="#64748B", lw=0.9, ls=":")
    ax.set_xlim(-0.04, 2.08)
    ax.set_ylim(-2.35, 2.35)
    ax.set_xlabel(r"Radial Separation $r$ ($N=2$) / Circumradius $R$ ($N=3$ Lagrange) [AU]")
    ax.set_ylabel(r"McGehee Regularized Radial Velocity $v_{\mathrm{McG}} = \dot{\rho}\sqrt{\rho}$")
    ax.set_title(r"(b) L2.1: McGehee Blown-Up Phase Portrait $(\rho, \dot{\rho}\sqrt{\rho})$ Crossing $\rho = 0$")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=7.6)

    # Panel (c): L3a Burrau 3-4-5 Scattering & L3b Full D_4 Quadruple Collision Bounce
    ax = axes[1, 0]
    t_b = np.array(l3["burrau"]["times"], dtype=float)
    r_min_b = _pos_or_nan(l3["burrau"]["min_separations"])
    rho_b = _pos_or_nan(l3["burrau"]["hyper_radii"])
    t_d4 = np.array(l3["quadruple_d4"]["times"], dtype=float)
    R_d4 = _pos_or_nan(l3["quadruple_d4"]["radius_R"])
    t_coll_d4 = l3["quadruple_d4"]["t_coll_analytical"]

    ax.semilogy(t_b, rho_b, "--", color="#64748B", lw=1.4, label=r"L3a Burrau Cluster Hyper-Radius $\bar{\rho}_S(t)$")
    ax.semilogy(t_b, r_min_b, "-", color="#1D4ED8", lw=1.6, label=rf"L3a Burrau $\min_{{i<j}} r_{{ij}}(t)$ ($r_{{\min}} = {l3['burrau']['global_min_separation']:.2e}$ AU)")
    ax.semilogy(t_d4 * 4.8, R_d4, "s-", color="#DC2626", lw=1.8, ms=3.2, label=rf"L3b $D_4$ Square Radius $R(t)$ ($t\times 4.8$, bounce at $t_{{\mathrm{{coll}}}}^{{(D_4)}} = {t_coll_d4:.4f}$)")
    ax.axvline(t_coll_d4 * 4.8, color="#DC2626", lw=1.3, ls=":", label=rf"$D_4$ Quadruple Collision Instant ($4.8\,t_{{\mathrm{{coll}}}}^{{(D_4)}} = {4.8*t_coll_d4:.2f}$)")

    ax.set_xlim(0.0, 12.0)
    ax.set_ylim(5e-4, 20.0)
    ax.set_xlabel(r"Physical Time $t$ (Burrau) / Rescaled Time $4.8\,t$ ($D_4$ Quadruple Collision)")
    ax.set_ylabel(r"Characteristic Separation / Radius [AU]")
    ax.set_title(r"(c) L3a Burrau 3-4-5 Scattering & L3b $D_4$ Quadruple Collision Bounce")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=7.4)
    _clean_log_axes(ax, False, True)

    # Panel (d): L4 3D Kozai-Lidov Secular Anti-Correlated Oscillations (e_in(t) <-> i_in(t))
    ax = axes[1, 1]
    t_kl = np.array(l4["times"], dtype=float)
    ecc_kl = np.array(l4["eccentricities"], dtype=float)
    inc_kl = np.array(l4["inclinations_deg"], dtype=float)

    ax.plot(t_kl, ecc_kl, "-", color="#1D4ED8", lw=1.8, label=rf"Inner Eccentricity $e_{{\mathrm{{in}}}}(t)$ ($e_{{\max}} = {l4['e_max']:.6f}$)")
    ax.set_xlim(0.0, 300.0)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel(r"Physical Time $t$ ($a_{\mathrm{in}}=1.0, a_{\mathrm{out}}=4.0, m_2=0.4, i_0=84.0^\circ$)")
    ax.set_ylabel(r"Inner Binary Eccentricity $e_{\mathrm{in}}(t)$", color="#1D4ED8")
    ax.tick_params(axis="y", labelcolor="#1D4ED8")

    ax2 = ax.twinx()
    ax2.plot(t_kl, inc_kl, "--", color="#DC2626", lw=1.7, label=rf"Mutual Inclination $i_{{\mathrm{{in}}}}(t)$ ($i \in [{l4['i_min']:.1f}^\circ, {l4['i_max']:.1f}^\circ]$)")
    ax2.set_ylim(30.0, 90.0)
    ax2.set_ylabel(r"Mutual Inclination $i_{\mathrm{in}}(t)$ [deg]", color="#DC2626")
    ax2.tick_params(axis="y", labelcolor="#DC2626")
    ax2.grid(False)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="lower right", framealpha=0.92, fontsize=7.8)
    ax.set_title(r"(d) L4: 3D Kozai–Lidov Eccentricity–Inclination Cycles ($1-e_{\min}=1.13\times 10^{-3}$)")

    for out_dir in (CH07_FIG_DIR, ARTIFACT_DIR):
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "gen03_plate_C_physical_dynamics.pdf")
        fig.savefig(out_dir / "gen03_plate_C_physical_dynamics.png", dpi=300)
    plt.close(fig)


def plot_plate_L1_pericenter_conditioning(data: Dict[str, Any]) -> None:
    """Render Publication Plate L1 (1x2) for §7.1.3: Analytic Pericenter Conditioning on N=2 across 1-e in [1e-12, 1e-1] & e=1."""
    l1 = data["L1_eccentricity_homotopy"]

    fig, axes = plt.subplots(1, 2, figsize=(14.6, 5.4), dpi=300)
    fig.subplots_adjust(top=0.84, bottom=0.14, left=0.07, right=0.97, wspace=0.25)
    fig.suptitle(
        r"Limit Case L1: Analytic Kepler–Sundman–Stumpff Pericenter Conditioning ($N=2$, $1 - e \in [10^{-12}, 10^{-1}]$ & $e = 1.0$)",
        fontweight="bold",
        y=0.96,
    )

    # Panel (a): Full-Orbit Pericenter Plunge r(t) Down to r_min = 1e-12 AU on Symmetric Log Time Scale Centered at t = 0.5 T_orb
    ax = axes[0]
    palette = ["#64748B", "#0891B2", "#059669", "#D97706", "#1D4ED8"]
    for idx, (k_str, orb) in enumerate(l1["sample_orbits"].items()):
        if "delta_u" in orb:
            du_arr = np.array(orb["delta_u"], dtype=float)
        else:
            du_arr = np.array(orb["times"], dtype=float) / l1["T_orb"] - 0.5
        r_arr = _pos_or_nan(orb["r_norm"])
        col = palette[idx % len(palette)]
        ax.plot(
            du_arr,
            r_arr,
            "-",
            color=col,
            lw=1.95,
            label=rf"$1 - e = {k_str}$ ($r_{{\min}} = {np.nanmin(r_arr):.1e}$ AU, $k_{{\mathrm{{sym}}}} \equiv 0$)",
        )

    du_step64 = 0.5 / 129.0
    ax.axvspan(
        -du_step64,
        +du_step64,
        color="#FEF3C7",
        alpha=0.38,
        label=r"Single Macro-Step 64 ($\Delta\tilde{t} \in [-\frac{1}{258}, +\frac{1}{258}]$, $k_{\mathrm{sym}} \equiv 0$)",
    )
    ax.axvline(0.0, color="#DC2626", lw=1.5, ls=":", label=r"Exact Pericenter Passage ($t = 0.5\,T_{\mathrm{orb}}$, $\Delta\tilde{t} = 0$)")
    ax.set_xscale("symlog", linthresh=1e-20, linscale=0.35)
    ax.set_yscale("log")
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(2e-16, 5.0)
    xt = [-1e-1, -1e-6, -1e-12, -1e-18, 0.0, 1e-18, 1e-12, 1e-6, 1e-1]
    ax.set_xticks(xt)
    ax.set_xticklabels(
        [r"$-10^{-1}$", r"$-10^{-6}$", r"$-10^{-12}$", r"$-10^{-18}$", r"$0$", r"$+10^{-18}$", r"$+10^{-12}$", r"$+10^{-6}$", r"$+10^{-1}$"],
        fontsize=7.8,
    )
    ax.set_xlabel(r"Pericenter-Centered Orbital Time $\Delta\tilde{t} = (t - 0.5\,T_{\mathrm{orb}}) / T_{\mathrm{orb}}$ [symlog, $N_{\mathrm{steps}} = 129$]")
    ax.set_ylabel(r"Binary Separation $r(t) = \|\mathbf{q}_1(t) - \mathbf{q}_0(t)\|$ [AU]")
    ax.set_title(r"(a) Pericenter-Centered Symlog Radial Plunge Down to $r_{\min} = 10^{-12}$ AU")
    ax.legend(loc="lower center", ncol=2, framealpha=0.95, fontsize=6.9, columnspacing=0.9)
    _clean_log_axes(ax, False, True)

    # Panel (b): Full-Orbit Energy Error Across 1 - e in [1e-12, 1e-1] (Small Values 1e-12 on Left, Pure N=2)
    ax = axes[1]
    ome = np.array(l1["one_minus_e"], dtype=float)
    verlet = _pos_or_nan(l1["cartesian_verlet_rel_dE"])
    legacy_uncond = _pos_or_nan(l1["legacy_unconditioned_f64_rel_dE"])
    leaf64 = _pos_or_nan(l1["float64_pbm4_rel_dE"])
    leaf128 = _pos_or_nan(l1["float128_pbm4_rel_dE"])

    ax.loglog(ome, verlet, "o--", color="#DC2626", lw=1.8, ms=5.5, label=r"Unregularized Cartesian Verlet ($\propto (1-e)^{-5/2}$ Divergence)")
    ax.loglog(ome, legacy_uncond, "d:", color="#7C3AED", lw=1.7, ms=5.0, label=r"Unconditioned Stumpff Subtraction ($\propto (1-e)^{-2}\epsilon_{\mathrm{mach}}$)")
    ax.loglog(
        ome,
        leaf64,
        "^-",
        color="#1D4ED8",
        lw=2.1,
        ms=6.0,
        label=r"Conditioned Sundman–Stumpff Leaf (Float64, $e=1$: " + f"{l1['rectilinear_e1_f64_rel_dE']:.2e})",
    )
    ax.loglog(
        ome,
        leaf128,
        "h-",
        color="#059669",
        lw=2.1,
        ms=6.0,
        label=r"Conditioned Sundman–Stumpff Leaf (Float128, $e=1$: " + f"{l1['rectilinear_e1_f128_rel_dE']:.2e})",
    )

    ax.set_xlim(5e-13, 2e-1)
    ax.set_ylim(1e-38, 1e6)
    ax.set_xlabel(r"Pericenter Proximity Parameter $1 - e = r_{\min} / a$ ($10^{-12}$ on Left $\to 10^{-1}$ on Right)")
    ax.set_ylabel(r"Full-Orbit Relative Energy Error $\max_{t\in[0,T_{\mathrm{orb}}]} |\Delta H(t)|/|H_0|$")
    ax.set_title(r"(b) Cancellation-Free Pericenter Conditioning Across 12 Decades")
    ax.legend(loc="center right", bbox_to_anchor=(0.98, 0.36), framealpha=0.94, fontsize=7.8)
    _clean_log_axes(ax, True, True)

    for out_dir in (CH07_FIG_DIR, ARTIFACT_DIR):
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "fig_ch07_l1_pericenter_conditioning.pdf")
        fig.savefig(out_dir / "fig_ch07_l1_pericenter_conditioning.png", dpi=300)
    plt.close(fig)


def plot_plate_L2_pbm4_free_fall_regularity(data: Dict[str, Any]) -> None:
    """Render Publication Plate L2 (1x2) for §7.1.3: Asymptotic Free-Fall Regularity of PBM4 (Z_4^{BM4}) over 12 Decades r in [1e-12, 1] AU."""
    l2 = data["L2_mcgehee_and_corrector"]
    sc = l2["c4_kick_scaling"]
    r_arr = np.array(sc["separations"], dtype=float)
    k4_arr = np.array(sc["k4_generator_abs"], dtype=float)
    dr_pbm4 = np.array(sc["pbm4_position_kick_norms"], dtype=float)
    dw_pbm4 = np.array(sc["pbm4_momentum_kick_norms"], dtype=float)
    dw_rel_pbm4 = np.array(sc["pbm4_relative_momentum_kick"], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(14.6, 5.4), dpi=300)
    fig.subplots_adjust(top=0.84, bottom=0.14, left=0.07, right=0.97, wspace=0.25)
    fig.suptitle(
        r"Limit Case L2: Asymptotic Free-Fall Regularity of the $\mathrm{PBM}_4$ ($\mathcal{Z}_4^{\mathrm{BM4}}$) Cotangent Corrector ($r \in [10^{-12}, 10^0]$ AU)",
        fontweight="bold",
        y=0.96,
    )

    # Panel (a): Extinction of Corrector Generating Potential |K_4(r,w)| ~ r^{+3.5} and Position Shift ||Delta r_corr|| ~ r^{+4.0}
    ax = axes[0]
    ax.loglog(
        r_arr,
        k4_arr,
        "s-",
        color="#059669",
        lw=2.0,
        ms=4.5,
        label=rf"Corrector Generating Potential $|K_4(\mathbf{{r}},\mathbf{{w}})| \propto r^{{+{sc['alpha_k4']:.4f}}}$ ($\mathcal{{O}}(r^{{7/2}}) \to 0$)",
    )
    ax.loglog(
        r_arr,
        dr_pbm4,
        "o-",
        color="#1D4ED8",
        lw=2.0,
        ms=4.5,
        label=rf"Corrector Position Shift $\|\Delta\mathbf{{r}}_{{\mathrm{{corr}}}}\| = \|\nabla_{{\mathbf{{w}}}}K_4\| \propto r^{{+{sc['alpha_dr_pbm4']:.4f}}}$ ($\mathcal{{O}}(r^4) \to 0$)",
    )
    ax.set_xlim(5e-13, 2.0)
    ax.set_ylim(1e-65, 1e-10)
    ax.set_xlabel(r"Inner Binary Separation $r = \|\mathbf{q}_1 - \mathbf{q}_0\|$ [AU] ($10^{-12}$ on Left $\to 10^0$ on Right)")
    ax.set_ylabel(r"Generating Potential $|K_4|$ & Position Shift $\|\Delta\mathbf{r}_{\mathrm{corr}}\|$")
    ax.set_title(r"(a) Extinction of Generating Potential $|K_4| = \mathcal{O}(r^{7/2})$ & Position Shift $\|\Delta\mathbf{r}_{\mathrm{corr}}\| = \mathcal{O}(r^4)$")
    ax.legend(loc="lower right", framealpha=0.94, fontsize=7.8)
    _clean_log_axes(ax, True, True)

    # Panel (b): Extinction of Absolute Momentum Shift ||Delta w_corr|| ~ r^{+2.5} & Relative Momentum Shift ||Delta w_corr||/||w|| ~ r^{+3.0}
    ax = axes[1]
    ax.loglog(
        r_arr,
        dw_pbm4,
        "h-",
        color="#D97706",
        lw=2.0,
        ms=4.5,
        label=rf"Absolute Momentum Shift $\|\Delta\mathbf{{w}}_{{\mathrm{{corr}}}}\| = \|\nabla_{{\mathbf{{r}}}}K_4\| \propto r^{{+{sc['alpha_dw_pbm4']:.4f}}}$ ($\mathcal{{O}}(r^{{5/2}}) \to 0$)",
    )
    ax.loglog(
        r_arr,
        dw_rel_pbm4,
        "^--",
        color="#0891B2",
        lw=2.0,
        ms=4.5,
        label=rf"Relative Momentum Shift $\|\Delta\mathbf{{w}}_{{\mathrm{{corr}}}}\|/\|\mathbf{{w}}\| \propto r^{{+{sc['alpha_dw_rel_pbm4']:.4f}}}$ ($\mathcal{{O}}(r^3) \to 0$)",
    )
    ax.set_xlim(5e-13, 2.0)
    ax.set_ylim(1e-53, 1e-10)
    ax.set_xlabel(r"Inner Binary Separation $r = \|\mathbf{q}_1 - \mathbf{q}_0\|$ Along Free-Fall Branch $w(r) \propto r^{-1/2}$ [AU]")
    ax.set_ylabel(r"Momentum Shifts $\|\Delta\mathbf{w}_{\mathrm{corr}}\|$ & $\|\Delta\mathbf{w}_{\mathrm{corr}}\|/\|\mathbf{w}\|$")
    ax.set_title(r"(b) Extinction of Absolute ($\mathcal{O}(r^{5/2}) \to 0$) & Relative ($\mathcal{O}(r^3) \to 0$) Momentum Shifts")
    ax.legend(loc="lower right", framealpha=0.94, fontsize=7.8)
    _clean_log_axes(ax, True, True)

    for out_dir in (CH07_FIG_DIR, ARTIFACT_DIR):
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "fig_ch07_l2_pbm4_free_fall_regularity.pdf")
        fig.savefig(out_dir / "fig_ch07_l2_pbm4_free_fall_regularity.png", dpi=300)
    plt.close(fig)


def plot_plate_L3_mcgehee_central_collisions(data: Dict[str, Any]) -> None:
    """Render Publication Plate L3 (2x3) for §7.1.3: McGehee Exceptional Divisor Traversals Split into 3 Columns (k=2, k=3 Lagrange, k=4 D4)."""
    l2 = data["L2_mcgehee_and_corrector"]
    l3 = data["L3_burrau_and_quadruple_d4"]
    tr_n2 = l2["traversal"]
    tr_lag = l2["lagrange_n3_mcgehee"]
    qd = l3["quadruple_d4"]

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.6), dpi=300)
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.06, right=0.98, wspace=0.26, hspace=0.34)
    fig.suptitle(
        r"Limit Case L3: Exact McGehee Exceptional Divisor Traversals ($\rho = 0$) & Conormal Bounce Involution $\mathcal{B}_C^2 = \mathrm{Id}$ ($k = 2, 3, 4$)",
        fontweight="bold",
        y=0.97,
    )

    # ------------------------------------------------------------------
    # COLUMN 1: k = 2 Binary Sundman Collision (r0 = 2.0 AU, t_coll = 2.2214)
    # ------------------------------------------------------------------
    ax = axes[0, 0]
    r_n2 = np.array(tr_n2["separations"], dtype=float)
    v_n2 = np.array(tr_n2["v_rad_mcgehee"], dtype=float)
    r_grid_n2 = np.linspace(0.0, 2.0, 300)
    v_th_n2 = np.sqrt(np.maximum(4.0 - 2.0 * r_grid_n2, 0.0))

    ax.plot(r_grid_n2, -v_th_n2, "--", color="#93C5FD", lw=1.6, label=r"Analytical Infall / Rebound Branches")
    ax.plot(r_grid_n2, +v_th_n2, "--", color="#93C5FD", lw=1.6)
    ax.plot(r_n2, v_n2, "o-", color="#1D4ED8", lw=1.7, ms=3.5, label=r"Discrete Trajectory ($r_0 = 2.0$ AU)")
    ax.annotate("", xy=(0.75, -1.58), xytext=(1.25, -1.22), arrowprops=dict(arrowstyle="->", color="#1D4ED8", lw=2.0))
    ax.annotate("", xy=(1.25, +1.22), xytext=(0.75, +1.58), arrowprops=dict(arrowstyle="->", color="#1D4ED8", lw=2.0))
    ax.annotate("", xy=(0.02, +1.96), xytext=(0.02, -1.96), arrowprops=dict(arrowstyle="->", color="#DC2626", lw=2.4))
    ax.text(
        0.14,
        0.0,
        r"$\mathcal{B}_C: -2.000 \mapsto +2.000$" + "\n" + r"at $r = 0$ ($t_{\mathrm{coll}} = 2.2214$)",
        color="#DC2626",
        fontsize=8.0,
        fontweight="bold",
        va="center",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#DC2626", alpha=0.94),
    )
    ax.axhline(0.0, color="#64748B", lw=0.9, ls=":")
    ax.set_xlim(-0.05, 2.10)
    ax.set_ylim(-2.35, 2.35)
    ax.set_xlabel(r"Binary Separation $r = \|\mathbf{q}_1 - \mathbf{q}_0\|$ [AU]")
    ax.set_ylabel(r"McGehee Velocity $v_{\mathrm{McG}} = \dot{r}\sqrt{r}$")
    ax.set_title(r"(a) $k = 2$ Binary Sundman Phase Portrait")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=7.6)

    ax = axes[1, 0]
    t_n2 = np.array(tr_n2["times"][1:], dtype=float)
    dE_n2_64 = _pos_or_nan(tr_n2["rel_energy_errors"][1:])
    dE_n2_128 = _pos_or_nan(tr_n2["rel_energy_errors_f128"][1:])
    t_coll_n2 = tr_n2["t_coll_analytical"]
    fb_n2_dict = l2["forward_backward"]
    t_fb_n2 = np.array([x["t_max"] for x in fb_n2_dict["float64_fb"]], dtype=float)
    fb_n2_64 = _pos_or_nan([x["delta_E_rel"] for x in fb_n2_dict["float64_fb"]])
    fb_n2_128 = _pos_or_nan([x["delta_E_rel"] for x in fb_n2_dict["float128_fb"]])

    ax.semilogy(t_n2, dE_n2_64, "-", color="#1D4ED8", lw=1.8, label=rf"Float64 $|\Delta H(t)|/|H_0|$ ($\max = {tr_n2['max_dE_f64']:.2e}$)")
    ax.semilogy(t_fb_n2, fb_n2_64, "o", color="#1D4ED8", ms=5.5, label=rf"Float64 Round-Trip $|\Delta H_{{\mathrm{{FB}}}}|/|H_0|$ (${fb_n2_64[-1]:.2e}$)")
    ax.semilogy(t_n2, dE_n2_128, "-", color="#059669", lw=1.8, label=rf"Float128 $|\Delta H(t)|/|H_0|$ ($\max = {tr_n2['max_dE_f128']:.2e}$)")
    ax.semilogy(t_fb_n2, fb_n2_128, "^", color="#059669", ms=5.5, label=rf"Float128 Round-Trip $|\Delta H_{{\mathrm{{FB}}}}|/|H_0|$ (${fb_n2_128[-1]:.2e}$)")
    ax.axvline(t_coll_n2, color="#DC2626", lw=1.6, ls=":", label=rf"Exact Collision $t_{{\mathrm{{coll}}}} = {t_coll_n2:.4f}$ ($r = 0$)")
    ax.axvspan(t_coll_n2, float(t_n2[-1]), color="#FEE2E2", alpha=0.28)
    ax.set_xlim(0.0, float(t_n2[-1]))
    ax.set_ylim(1e-38, 1e-10)
    ax.set_xlabel(r"Physical Time $t$ ($\Delta t = 0.025, N_{\mathrm{steps}} = 177$)")
    ax.set_ylabel(r"Relative Hamiltonian Error $|\Delta H(t)| / |H_0|$")
    ax.set_title(r"(d) $k = 2$ Energy & Round-Trip Across Binary Collision")
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.0)
    _clean_log_axes(ax, False, True)

    # ------------------------------------------------------------------
    # COLUMN 2: k = 3 Equilateral Lagrange Triple Collision (R0 = 1.20 AU, t_coll = 1.9221)
    # ------------------------------------------------------------------
    ax = axes[0, 1]
    r_lag = np.array(tr_lag["hyper_radius_R"], dtype=float)
    v_lag = np.array(tr_lag["v_rad_mcgehee"], dtype=float)
    gm_eff_lag = 1.0 / np.sqrt(3.0)
    r_grid_lag = np.linspace(0.0, tr_lag["R0"], 300)
    v_th_lag = np.sqrt(np.maximum(2.0 * gm_eff_lag * (1.0 - r_grid_lag / tr_lag["R0"]), 0.0))

    ax.plot(r_grid_lag, -v_th_lag, "--", color="#FCD34D", lw=1.6, label=r"Analytical McGehee Homothetic Manifold")
    ax.plot(r_grid_lag, +v_th_lag, "--", color="#FCD34D", lw=1.6)
    ax.plot(r_lag, v_lag, "s-", color="#D97706", lw=1.7, ms=3.5, label=r"Discrete Trajectory ($R_0 = 1.20$ AU)")
    ax.annotate("", xy=(0.45, -0.85), xytext=(0.78, -0.64), arrowprops=dict(arrowstyle="->", color="#D97706", lw=2.0))
    ax.annotate("", xy=(0.78, +0.64), xytext=(0.45, +0.85), arrowprops=dict(arrowstyle="->", color="#D97706", lw=2.0))
    ax.annotate("", xy=(0.015, +1.05), xytext=(0.015, -1.05), arrowprops=dict(arrowstyle="->", color="#DC2626", lw=2.4))
    ax.text(
        0.10,
        0.0,
        r"$\mathcal{B}_C: -1.075 \mapsto +1.075$" + "\n" + r"at $R = 0$ ($t_{\mathrm{coll}}^{(L_3)} = 1.9216$)",
        color="#DC2626",
        fontsize=8.0,
        fontweight="bold",
        va="center",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#DC2626", alpha=0.94),
    )
    ax.axhline(0.0, color="#64748B", lw=0.9, ls=":")
    ax.set_xlim(-0.03, 1.28)
    ax.set_ylim(-1.35, 1.35)
    ax.set_xlabel(r"Equilateral Circumradius $R(t) = \|\mathbf{q}_i(t)\|$ [AU]")
    ax.set_ylabel(r"McGehee Velocity $v_{\mathrm{McG}} = \dot{R}\sqrt{R}$")
    ax.set_title(r"(b) $k = 3$ Equilateral Lagrange Phase Portrait")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=7.6)

    ax = axes[1, 1]
    t_lag = np.array(tr_lag["times"][1:], dtype=float)
    dE_lag_64 = _pos_or_nan(tr_lag["rel_energy_errors_f64"][1:])
    dE_lag_128 = _pos_or_nan(tr_lag["rel_energy_errors_f128"][1:])
    t_coll_lag = tr_lag["t_coll_analytical"]
    t_fb_lag = np.array([x["t_max"] for x in fb_n2_dict["lagrange_n3_float64_fb"]], dtype=float)
    fb_lag_64 = _pos_or_nan([x["delta_E_rel"] for x in fb_n2_dict["lagrange_n3_float64_fb"]])
    fb_lag_128 = _pos_or_nan([x["delta_E_rel"] for x in fb_n2_dict["lagrange_n3_float128_fb"]])

    ax.semilogy(t_lag, dE_lag_64, "--", color="#D97706", lw=1.8, label=rf"Float64 $|\Delta H(t)|/|H_0|$ ($\max = {tr_lag['max_dE_f64']:.2e}$)")
    ax.semilogy(t_fb_lag, fb_lag_64, "o", color="#1D4ED8", ms=5.5, label=rf"Float64 Round-Trip $|\Delta H_{{\mathrm{{FB}}}}|/|H_0|$ (${fb_lag_64[-1]:.2e}$)")
    ax.semilogy(t_lag, dE_lag_128, "-", color="#0891B2", lw=1.8, label=rf"Float128 $|\Delta H(t)|/|H_0|$ ($\max = {tr_lag['max_dE_f128']:.2e}$)")
    ax.semilogy(t_fb_lag, fb_lag_128, "^", color="#059669", ms=5.5, label=rf"Float128 Round-Trip $|\Delta H_{{\mathrm{{FB}}}}|/|H_0|$ (${fb_lag_128[-1]:.2e}$)")
    ax.axvline(t_coll_lag, color="#DC2626", lw=1.6, ls=":", label=rf"Exact Triple Collision $t_{{\mathrm{{coll}}}}^{{(L_3)}} = {t_coll_lag:.4f}$")
    ax.axvspan(t_coll_lag, float(t_lag[-1]), color="#FEE2E2", alpha=0.28)
    ax.set_xlim(0.0, float(t_lag[-1]))
    ax.set_ylim(1e-38, 1e-10)
    ax.set_xlabel(r"Physical Time $t$ ($\Delta t = 0.025, N_{\mathrm{steps}} = 153$)")
    ax.set_ylabel(r"Relative Hamiltonian Error $|\Delta H(t)| / |H_0|$")
    ax.set_title(r"(e) $k = 3$ Energy & Round-Trip Across Triple Collision")
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.0)
    _clean_log_axes(ax, False, True)

    # ------------------------------------------------------------------
    # COLUMN 3: k = 4 Symmetric Square D_4 Quadruple Collision (R0 = 1.00 AU, t_coll = 1.1353)
    # ------------------------------------------------------------------
    ax = axes[0, 2]
    r_d4 = np.array(qd["radius_R"], dtype=float)
    v_d4 = np.array(qd["v_rad_mcgehee"], dtype=float)
    gm_eff_d4 = 0.25 + 1.0 / np.sqrt(2.0)
    r_grid_d4 = np.linspace(0.0, 1.0, 300)
    v_th_d4 = np.sqrt(np.maximum(2.0 * gm_eff_d4 * (1.0 - r_grid_d4), 0.0))

    ax.plot(r_grid_d4, -v_th_d4, "--", color="#FCA5A5", lw=1.6, label=r"Analytical $D_4$ Homothetic Manifold")
    ax.plot(r_grid_d4, +v_th_d4, "--", color="#FCA5A5", lw=1.6)
    ax.plot(r_d4, v_d4, "^-", color="#DC2626", lw=1.7, ms=3.5, label=r"Discrete Trajectory ($R_0 = 1.00$ AU)")
    ax.annotate("", xy=(0.38, -1.09), xytext=(0.66, -0.81), arrowprops=dict(arrowstyle="->", color="#DC2626", lw=2.0))
    ax.annotate("", xy=(0.66, +0.81), xytext=(0.38, +1.09), arrowprops=dict(arrowstyle="->", color="#DC2626", lw=2.0))
    ax.annotate("", xy=(0.015, +1.35), xytext=(0.015, -1.35), arrowprops=dict(arrowstyle="->", color="#7C3AED", lw=2.4))
    ax.text(
        0.09,
        0.0,
        r"$\mathcal{B}_C: -1.384 \mapsto +1.384$" + "\n" + r"at $R = 0$ ($t_{\mathrm{coll}}^{(D_4)} = 1.1353$)",
        color="#7C3AED",
        fontsize=8.0,
        fontweight="bold",
        va="center",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#7C3AED", alpha=0.94),
    )
    ax.axhline(0.0, color="#64748B", lw=0.9, ls=":")
    ax.set_xlim(-0.03, 1.06)
    ax.set_ylim(-1.68, 1.68)
    ax.set_xlabel(r"Square Circumradius $R(t) = \|\mathbf{q}_i(t)\|$ [AU]")
    ax.set_ylabel(r"McGehee Velocity $v_{\mathrm{McG}} = \dot{R}\sqrt{R}$")
    ax.set_title(r"(c) $k = 4$ Symmetric Square $D_4$ Phase Portrait")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=7.6)

    ax = axes[1, 2]
    t_d4 = np.array(qd["times"][1:], dtype=float)
    dE_d4_64 = _pos_or_nan(qd["rel_energy_errors_f64"][1:])
    dE_d4_128 = _pos_or_nan(qd["rel_energy_errors_f128"][1:])
    t_coll_d4 = qd["t_coll_analytical"]
    t_fb = np.array([x["t_max"] for x in qd["float64_fb"]], dtype=float)
    fb_dE64 = _pos_or_nan([x["delta_E_rel"] for x in qd["float64_fb"]])
    fb_dE128 = _pos_or_nan([x["delta_E_rel"] for x in qd["float128_fb"]])

    ax.semilogy(t_d4, dE_d4_64, "-.", color="#DC2626", lw=1.8, label=rf"Float64 $|\Delta H(t)|/|H_0|$ ($\max = {qd['max_dE_f64']:.2e}$)")
    ax.semilogy(t_fb, fb_dE64, "o", color="#1D4ED8", ms=5.5, label=rf"Float64 Round-Trip $|\Delta H_{{\mathrm{{FB}}}}|/|H_0|$ (${fb_dE64[-1]:.2e}$)")
    ax.semilogy(t_d4, dE_d4_128, "-", color="#7C3AED", lw=1.8, label=rf"Float128 $|\Delta H(t)|/|H_0|$ ($\max = {qd['max_dE_f128']:.2e}$)")
    ax.semilogy(t_fb, fb_dE128, "^", color="#059669", ms=5.5, label=rf"Float128 Round-Trip $|\Delta H_{{\mathrm{{FB}}}}|/|H_0|$ (${fb_dE128[-1]:.2e}$)")
    ax.axvline(t_coll_d4, color="#DC2626", lw=1.6, ls=":", label=rf"Exact $D_4$ Collision $t_{{\mathrm{{coll}}}}^{{(D_4)}} = {t_coll_d4:.4f}$")
    ax.axvspan(t_coll_d4, float(t_d4[-1]), color="#FEE2E2", alpha=0.28)
    ax.set_xlim(0.0, float(t_d4[-1]))
    ax.set_ylim(1e-38, 1e-10)
    ax.set_xlabel(r"Physical Time $t$ ($\Delta t = 0.025, N_{\mathrm{steps}} = 100, \Delta_{D_4} = 0.00$)")
    ax.set_ylabel(r"Relative Hamiltonian Error $|\Delta H(t)| / |H_0|$")
    ax.set_title(r"(f) $k = 4$ Energy & Round-Trip Across $D_4$ Collision")
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.0)
    _clean_log_axes(ax, False, True)

    for out_dir in (CH07_FIG_DIR, ARTIFACT_DIR):
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "fig_ch07_l3_mcgehee_central_collisions.pdf")
        fig.savefig(out_dir / "fig_ch07_l3_mcgehee_central_collisions.png", dpi=300)
    plt.close(fig)


def plot_plate_L4_adaptive_multiscale_burrau_kozai(data: Dict[str, Any]) -> None:
    """Render Publication Plate L4 (2x2) for §7.1.3: Adaptive Multiscale Resolution L4a Burrau & L4b Kozai-Lidov with True Time Series."""
    l3 = data["L3_burrau_and_quadruple_d4"]
    l4 = data["L4_kozai_lidov"]
    bur = l3["burrau"]

    fig, axes = plt.subplots(2, 2, figsize=(15.0, 10.8), dpi=300)
    fig.subplots_adjust(top=0.91, bottom=0.08, left=0.08, right=0.92, wspace=0.34, hspace=0.35)
    fig.suptitle(
        r"Limit Case L4: Adaptive Multiscale Resolution Without Exact Collision — Burrau 3-4-5 Scattering (L4a) & 3D Kozai–Lidov Cycles (L4b)",
        fontweight="bold",
        y=0.97,
    )

    # Panel (a): L4a Burrau 3-4-5 Chaotic Scattering & Dyadic Level k_sym(t) over t in [0, 12.0]
    ax = axes[0, 0]
    t_b = np.array(bur["times"], dtype=float)
    r_min_b = _pos_or_nan(bur["min_separations"])
    rho_b = _pos_or_nan(bur["hyper_radii"])
    k_b = np.array(bur["dyadic_levels"], dtype=float)[-len(t_b):]

    ax.semilogy(t_b, rho_b, "--", color="#64748B", lw=1.4, label=r"Cluster Hyper-Radius $\bar{\rho}_S(t)$ [AU]")
    ax.semilogy(t_b, r_min_b, "-", color="#1D4ED8", lw=1.7, label=rf"Min Pairwise Distance $\min_{{i<j}} r_{{ij}}(t)$ ($r_{{\min}} = {bur['global_min_separation']:.2e}$ AU)")
    ax.set_xlim(0.0, 12.0)
    ax.set_ylim(5e-3, 120.0)
    ax.set_xlabel(r"Physical Time $t$ (Pythagorean 3-4-5 Problem, $m = \{3,4,5\}$)")
    ax.set_ylabel(r"Pairwise Separation / Hyper-Radius [AU]", color="#1D4ED8")
    ax.tick_params(axis="y", labelcolor="#1D4ED8")
    _clean_log_axes(ax, False, True)

    ax_k = ax.twinx()
    ax_k.step(t_b, k_b, where="post", color="#DC2626", lw=1.4, alpha=0.82, label=rf"Symmetric Dyadic Rung $k_{{\mathrm{{sym}}}}(t) \in \{{0,\dots,{int(np.max(k_b))}\}}$")
    ax_k.set_ylim(-0.5, 13.5)
    ax_k.set_ylabel(r"Symmetric Dyadic Level $k_{\mathrm{sym}}(t)$", color="#DC2626")
    ax_k.tick_params(axis="y", labelcolor="#DC2626")
    ax_k.grid(False)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax_k.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper right", framealpha=0.94, fontsize=7.3)
    ax.set_title(rf"(a) L4a: Burrau 3-4-5 Chaotic Scattering & Dyadic Adaptivity ($k_{{\max}}={int(np.max(k_b))}$)")

    # Panel (b): L4a Burrau 3-4-5 True Time Series of Energy & Noether Invariants (Float64 & Float128)
    ax = axes[0, 1]
    dE_b64 = _pos_or_nan(bur["rel_energy_errors"])
    dE_b128 = _pos_or_nan(bur["rel_energy_errors_f128"])
    abs_L_b64 = _pos_or_nan(bur.get("abs_L_history_f64", [bur["lie_bargmann_drift_f64"]["abs_L"]] * len(t_b)))
    abs_P_b64 = _pos_or_nan(bur.get("abs_P_history_f64", [bur["lie_bargmann_drift_f64"]["abs_P"]] * len(t_b)))
    abs_L_b128 = _pos_or_nan(bur.get("abs_L_history_f128", [bur["lie_bargmann_drift_f128"]["abs_L"]] * len(t_b)))
    abs_P_b128 = _pos_or_nan(bur.get("abs_P_history_f128", [bur["lie_bargmann_drift_f128"]["abs_P"]] * len(t_b)))

    ax.semilogy(t_b, dE_b64, "-", color="#1D4ED8", lw=1.7, label=rf"Float64/128 Shadow Energy $|\Delta H(t)|/|H_0|$ ($\mathcal{{O}}(h_{{\mathrm{{sub}}}}^6)$, $\max = {np.nanmax(dE_b64):.2e}$)")
    ax.semilogy(t_b, dE_b128, ":", color="#0EA5E9", lw=1.3, alpha=0.88)
    ax.semilogy(t_b, abs_L_b64, "-", color="#D97706", lw=1.7, label=rf"Float64 Noether $\|\Delta\mathbf{{L}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\max = {np.nanmax(abs_L_b64):.2e}$)")
    ax.semilogy(t_b, abs_P_b64, "--", color="#0891B2", lw=1.5, label=rf"Float64 Linear Momentum $\|\Delta\mathbf{{P}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\max = {np.nanmax(abs_P_b64):.2e}$)")
    ax.semilogy(t_b, abs_L_b128, "-", color="#059669", lw=1.8, label=rf"Float128 Noether $\|\Delta\mathbf{{L}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\max = {np.nanmax(abs_L_b128):.2e}$)")
    ax.semilogy(t_b, abs_P_b128, ":", color="#7C3AED", lw=1.7, label=rf"Float128 Linear Momentum $\|\Delta\mathbf{{P}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\max = {np.nanmax(abs_P_b128):.2e}$)")

    ax.set_xlim(0.0, 12.0)
    ax.set_ylim(5e-38, 1e-5)
    ax.set_xlabel(r"Physical Time $t$ (Burrau 3-4-5 Chaotic Scattering)")
    ax.set_ylabel(r"Time-Evolving Hamiltonian & Noether Residuals")
    ax.set_title(r"(b) L4a: Burrau 3-4-5 Shadow Energy ($\mathcal{O}(h_{\mathrm{sub}}^6)$) & Lie–Bargmann Noether Invariants")
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.0)
    _clean_log_axes(ax, False, True)

    # Panel (c): L4b 3D Kozai-Lidov Secular Cycles (e_in(t) <-> i_in(t)) & Dyadic Level k_sym(t)
    ax = axes[1, 0]
    t_kl = np.array(l4["times"], dtype=float)
    ecc_kl = np.array(l4["eccentricities"], dtype=float)
    inc_kl = np.array(l4["inclinations_deg"], dtype=float)
    k_kl = np.array(l4["dyadic_levels"], dtype=float)

    ax.plot(t_kl, ecc_kl, "-", color="#1D4ED8", lw=1.9, label=rf"Inner Eccentricity $e_{{\mathrm{{in}}}}(t)$ ($e_{{\max}} = {l4['e_max']:.6f}$, $1-e_{{\min}} = {l4['one_minus_e_min']:.2e}$)")
    ax.fill_between(
        t_kl,
        0.0,
        1.05,
        where=(k_kl >= 6),
        color="#A7F3D0",
        alpha=0.35,
        step="mid",
        label=rf"Peak Dyadic Refinement Windows ($k_{{\mathrm{{sym}}}}(t) \in [6, {int(np.max(k_kl))}]$)",
    )
    ax.set_xlim(0.0, 300.0)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel(r"Physical Time $t$ ($a_{\mathrm{in}}=1.0, a_{\mathrm{out}}=4.0, m_2=0.4, i_0=84.0^\circ$)")
    ax.set_ylabel(r"Inner Binary Eccentricity $e_{\mathrm{in}}(t)$", color="#1D4ED8")
    ax.tick_params(axis="y", labelcolor="#1D4ED8")

    ax2 = ax.twinx()
    ax2.plot(t_kl, inc_kl, "--", color="#DC2626", lw=1.7, label=rf"Mutual Inclination $i_{{\mathrm{{in}}}}(t) \in [{l4['i_min']:.1f}^\circ, {l4['i_max']:.1f}^\circ]$")
    ax2.set_ylim(30.0, 90.0)
    ax2.set_ylabel(r"Mutual Inclination $i_{\mathrm{in}}(t)$ [deg]", color="#DC2626")
    ax2.tick_params(axis="y", labelcolor="#DC2626")
    ax2.grid(False)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="center left", bbox_to_anchor=(0.02, 0.56), framealpha=0.94, fontsize=7.2)
    ax.set_title(r"(c) L4b: 3D Kozai–Lidov Eccentricity–Inclination Cycles & Dyadic Rung $k_{\mathrm{sym}}(t)$")

    # Panel (d): L4b Kozai-Lidov True Time Series of Semisecular Action, Energy & Exact Noether L_{z,tot}(t)
    ax = axes[1, 1]
    t_kl64 = np.array(l4["times"], dtype=float)
    dLz_tot64 = _pos_or_nan(l4["total_Lz_drifts_f64"])
    t_kl_del = np.array(l4["times"][1:], dtype=float)
    delaunay_dev = _pos_or_nan(l4["kozai_integral_deviations"][1:])
    t_kl128 = np.array(l4["times_f128"], dtype=float)
    dLz_tot128 = _pos_or_nan(l4["total_Lz_drifts_f128"])
    dE_kl64 = _pos_or_nan(l4["rel_energy_errors"])
    dE_kl128 = _pos_or_nan(l4["rel_energy_errors_f128"])

    ax.semilogy(
        t_kl_del,
        delaunay_dev,
        "-",
        color="#D97706",
        lw=1.7,
        label=rf"Semisecular Quadrupole Action $|\Delta(\sqrt{{1-e_{{\mathrm{{in}}}}^2}}\cos i_{{\mathrm{{in}}}})(t)|$ ($\max = {l4['max_kozai_int_dev']:.2e}$)",
    )
    ax.semilogy(
        t_kl64,
        dE_kl64,
        "-",
        color="#7C3AED",
        lw=1.4,
        alpha=0.80,
        label=rf"Float64 Shadow Energy $|\Delta H(t)|/|H_0|$ ($\mathcal{{O}}(h_{{\mathrm{{sub}}}}^6)$, $\max = {l4['max_rel_dE']:.2e}$)",
    )
    ax.semilogy(
        t_kl128,
        dE_kl128,
        "-.",
        color="#0891B2",
        lw=1.6,
        label=rf"Float128 Shadow Energy $|\Delta H(t)|/|H_0|$ ($\mathcal{{O}}(h_{{\mathrm{{sub}}}}^6)$, $\min = {np.nanmin(dE_kl128):.2e}$)",
    )
    ax.semilogy(
        t_kl64,
        dLz_tot64,
        "-",
        color="#1D4ED8",
        lw=1.7,
        label=rf"Exact Noether $|L_{{z,\mathrm{{tot}}}}(t) - L_{{z,\mathrm{{tot}}}}(0)|$ (Float64, $\max = {l4['max_total_Lz_drift_f64']:.2e}$)",
    )
    ax.semilogy(
        t_kl128,
        dLz_tot128,
        "--",
        color="#059669",
        lw=1.8,
        label=rf"Exact Noether $|L_{{z,\mathrm{{tot}}}}(t) - L_{{z,\mathrm{{tot}}}}(0)|$ (Float128, $\max = {l4['max_total_Lz_drift_f128']:.2e}$)",
    )

    ax.set_xlim(0.0, 305.0)
    ax.set_ylim(1e-37, 1e1)
    ax.set_xlabel(r"Physical Time $t$ ($a_{\mathrm{out}}/a_{\mathrm{in}} = 4.0, i_0 = 84.0^\circ$)")
    ax.set_ylabel(r"Time-Evolving Action / Energy / Noether Residual")
    ax.set_title(r"(d) L4b: Shadow Energy ($\mathcal{O}(h_{\mathrm{sub}}^6)$) & Exact Noether $L_{z,\mathrm{tot}}(t)$ vs. Semisecular Action")
    ax.legend(loc="center right", bbox_to_anchor=(0.98, 0.40), framealpha=0.94, fontsize=6.8)
    _clean_log_axes(ax, False, True)

    for out_dir in (CH07_FIG_DIR, ARTIFACT_DIR):
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "fig_ch07_l4_adaptive_multiscale_burrau_kozai.pdf")
        fig.savefig(out_dir / "fig_ch07_l4_adaptive_multiscale_burrau_kozai.png", dpi=300)
    plt.close(fig)


def main() -> None:
    assert JSON_PATH.exists(), f"Missing certified JSON: {JSON_PATH}"
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)
    _configure_matplotlib()
    plot_plate_L1_pericenter_conditioning(data)
    plot_plate_L2_pbm4_free_fall_regularity(data)
    plot_plate_L3_mcgehee_central_collisions(data)
    plot_plate_L4_adaptive_multiscale_burrau_kozai(data)
    print("=== [Plotter Generator 03 COMPLETE] Rendered 4 Dedicated Paper Plates L1, L2, L3, L4 (PDF & PNG) ===", flush=True)


if __name__ == "__main__":
    main()
