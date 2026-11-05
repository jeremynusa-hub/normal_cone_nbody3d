#!/usr/bin/env python3
r"""
Pure 3-Volet Laskar Plotter for Benchmark Generator 02:
Universal Obstruction of P_0, Shadow Manifold Z_4^{BM4}, Quadratic Floor &
Dynamic Atlas Gauge 1-Cocycle (§7.1.2, Stress-Tests ST-1 to ST-4).

Reads ONLY the certified JSON output:
  results/campaigns/generator_02_universal_obstruction_and_shadow_gauge.json

Generates 3 standalone spacious PDF/PNG plates (zero np.clip, zero overlapping insets,
zero synthetic curves):
  - Volet A: gen02_plate_A_computational.pdf / .png (2x2 grid)
  - Volet B: gen02_plate_B_algebraic_invariants.pdf / .png (2x2 grid)
  - Volet C: gen02_plate_C_physical_dynamics.pdf / .png (1x2 grid)
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
import matplotlib.ticker as mticker

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = WORKSPACE_ROOT / "results" / "campaigns" / "generator_02_universal_obstruction_and_shadow_gauge.json"
CH07_FIG_DIR = WORKSPACE_ROOT / "results" / "figures"
ARTIFACT_DIR = WORKSPACE_ROOT / "results" / "figures"

# Canonical color & marker palette (Section 2.3 of 09_pillar2_campaign_audit_and_laskar_plate_doctrine.md)
COL_BM4 = "#D97706"      # Amber / Orange
COL_VIOLET = "#7C3AED"   # Violet
COL_CYAN = "#0891B2"     # Cyan / Teal
COL_Z4 = "#1D4ED8"       # Royal Blue (Canonical Homological Z_4^{BM4})
COL_EPS_1E2 = "#DC2626"  # Crimson Red
COL_EPS_1E3 = "#0891B2"  # Cyan / Teal
COL_EPS_1E4 = "#1D4ED8"  # Royal Blue
COL_EMERALD = "#059669"  # Emerald Green
COL_SLATE = "#475569"    # Slate Grey


def _setup_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10.5,
        "axes.titlesize": 11.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 10.5,
        "legend.fontsize": 8.8,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "axes.grid": True,
        "grid.alpha": 0.28,
        "grid.linestyle": "--",
        "figure.dpi": 300,
    })


def _pos_or_nan(arr: Any) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    return np.where(a > 0.0, a, np.nan)


def _save_and_export(fig: plt.Figure, stem: str, legacy_names: list[str] | None = None) -> None:
    CH07_FIG_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = CH07_FIG_DIR / f"{stem}.pdf"
    png_path = CH07_FIG_DIR / f"{stem}.png"
    fig.savefig(pdf_path, bbox_inches="tight", dpi=300)
    fig.savefig(png_path, bbox_inches="tight", dpi=300)
    if ARTIFACT_DIR.exists() and ARTIFACT_DIR.resolve() != CH07_FIG_DIR.resolve():
        shutil.copy2(pdf_path, ARTIFACT_DIR / f"{stem}.pdf")
        shutil.copy2(png_path, ARTIFACT_DIR / f"{stem}.png")
    if legacy_names:
        for leg in legacy_names:
            shutil.copy2(png_path, CH07_FIG_DIR / leg)
            if leg.endswith(".png"):
                shutil.copy2(pdf_path, CH07_FIG_DIR / (leg[:-4] + ".pdf"))
            if ARTIFACT_DIR.exists() and ARTIFACT_DIR.resolve() != CH07_FIG_DIR.resolve():
                shutil.copy2(png_path, ARTIFACT_DIR / leg)
                if leg.endswith(".png"):
                    shutil.copy2(pdf_path, ARTIFACT_DIR / (leg[:-4] + ".pdf"))
    plt.close(fig)
    print(f"[RENDERED] {pdf_path.name} & {png_path.name}", flush=True)


def plot_volet_A_computational(data: Dict[str, Any]) -> None:
    """
    Publication Plate for Subsection 7.1.2 (2x2):
      (a) 35-Digit Float128 Energy Error Crossover (p=6 -> p=4) on 3D Eccentric Triple
      (b) Unclipped Empirical Convergence Order p_emp(h) vs Two-Channel Model
      (c) Two-Regime Tidal Scaling dE(eps) for BM4 vs PBM4 (O(eps^2*h^4) vs O(eps*h^6))
      (d) 35-Digit Float128 Trajectory Position Error Crossover ||q(T) - q_ref||_2 vs h_eff
    """
    fig, axes = plt.subplots(2, 2, figsize=(14.4, 10.6))
    fig.subplots_adjust(wspace=0.30, hspace=0.34)

    c_f128 = data["crossover_3d_f128"]["sweeps"]
    eps_meta = [
        ("eps_1e-02", r"\varepsilon = 10^{-2}", r"$R_{\mathrm{out}}=2.92$", COL_EPS_1E2, "s"),
        ("eps_1e-03", r"\varepsilon = 10^{-3}", r"$R_{\mathrm{out}}=6.30$", COL_EPS_1E3, "^"),
        ("eps_1e-04", r"\varepsilon = 10^{-4}", r"$R_{\mathrm{out}}=13.57$", COL_EPS_1E4, "o"),
    ]

    # --- Panel (a): Float128 Energy Crossover on 3D Eccentric Triple ---
    ax = axes[0, 0]
    for key, eps_tex, rout_tex, col, mk in eps_meta:
        entry = c_f128[key]
        h_eff = np.array(entry["tau_over_s"], dtype=float)
        dE_b = np.array(entry["dE_bm4"], dtype=float)
        dE_p = np.array(entry["dE_pbm4"], dtype=float)
        ax.loglog(h_eff, dE_b, linestyle="--", marker=mk, color=col, alpha=0.55, linewidth=1.4,
                  markersize=5.0, label=rf"$\mathrm{{BM}}_4^{{(7)}}$ (${eps_tex}$)")
        ax.loglog(h_eff, dE_p, linestyle="-", marker=mk, color=col, linewidth=2.2,
                  markersize=6.2, label=rf"$\mathrm{{PBM}}_4$ (${eps_tex}$, {rout_tex})")

    # Reference slopes p=6 and p=4
    h_ref = np.array(c_f128["eps_1e-04"]["tau_over_s"], dtype=float)
    y_ref6 = c_f128["eps_1e-04"]["dE_pbm4"][0] * (h_ref / h_ref[0]) ** 6
    y_ref4 = c_f128["eps_1e-02"]["dE_pbm4"][-1] * (h_ref / h_ref[-1]) ** 4
    ax.loglog(h_ref, y_ref6 * 0.45, ":", color="#1E293B", linewidth=1.3, label=r"Order-6 channel $\propto \varepsilon h_{\mathrm{eff}}^6$")
    ax.loglog(h_ref, y_ref4 * 1.8, "-.", color="#64748B", linewidth=1.3, label=r"Quadratic floor $\propto \varepsilon^2 h_{\mathrm{eff}}^4$")

    max_gain_1e4 = max(c_f128["eps_1e-04"]["gain_pbm4_vs_bm4"])
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_title(r"(a) Float128 Energy Error Crossover ($e_{\mathrm{in}}=0.30, i=28^\circ$)")
    ax.set_xlabel(r"Normalized Inverse Cost $h_{\mathrm{eff}} = h / s_{\mathrm{eval}}$ ($s_{\mathrm{eval}} = 6$)")
    ax.set_ylabel(r"Relative Hamiltonian Error $|\Delta H(T)| / |H_0|$")
    ax.text(0.04, 0.91, f"Max $\\mathrm{{PBM}}_4$ energy gain ($\\varepsilon=10^{{-4}}$): ${max_gain_1e4:,.0f}\\times$",
            transform=ax.transAxes, fontsize=9.0, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#EFF6FF", edgecolor=COL_Z4, alpha=0.92))
    ax.legend(loc="lower right", ncol=1, framealpha=0.92, fontsize=7.6)

    # --- Panel (b): Unclipped Empirical Convergence Order p_emp(h) ---
    ax = axes[0, 1]
    # Universal two-channel parameters C_6 and C_4 from the asymptotic endpoints
    c6_univ = c_f128["eps_1e-04"]["dq_pbm4"][1] / (1e-4 * (c_f128["eps_1e-04"]["h_steps"][1] ** 6))
    c4_univ = c_f128["eps_1e-02"]["dq_pbm4"][-1] / ((1e-2 ** 2) * (c_f128["eps_1e-02"]["h_steps"][-1] ** 4))
    h_dense = np.geomspace(0.012, 0.25, 200)
    for key, eps_tex, _rout_tex, col, mk in eps_meta:
        entry = c_f128[key]
        eps_val = float(entry["epsilon"])
        hm = np.array(entry["h_midpoints"], dtype=float)
        sl_q = np.array(entry["local_slopes_q_pbm4"], dtype=float)
        sl_b = np.array(entry["local_slopes_q_bm4"], dtype=float)
        term6_d = c6_univ * eps_val * (h_dense ** 6)
        term4_d = c4_univ * (eps_val ** 2) * (h_dense ** 4)
        sl_th_dense = (6.0 * term6_d + 4.0 * term4_d) / (term6_d + term4_d)
        ax.semilogx(hm, sl_q, linestyle="-", marker=mk, color=col, linewidth=2.2, markersize=6.2,
                    label=rf"$\mathrm{{PBM}}_4$ empirical (${eps_tex}$)")
        ax.semilogx(h_dense, sl_th_dense, linestyle=":", color=col, linewidth=1.6, alpha=0.88,
                    label=rf"Two-channel law $p_{{\mathrm{{th}}}}(h;\varepsilon)$ (${eps_tex}$)")
        if key == "eps_1e-04":
            ax.semilogx(hm, sl_b, linestyle="--", marker="d", color=COL_BM4, linewidth=1.7, markersize=5.2,
                        label=r"$\mathrm{BM}_4^{(7)}$ empirical ($p \approx 4.00$)")

    ax.axhline(6.0, color="#0F172A", linestyle="--", linewidth=1.1, alpha=0.7)
    ax.axhline(4.0, color="#475569", linestyle="-.", linewidth=1.1, alpha=0.7)
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_ylim(3.6, 6.85)
    ax.set_title(r"(b) Unclipped Local Empirical Order $p_{\mathrm{emp}}(h) : 6 \to 4$")
    ax.set_xlabel(r"Macro-Step Size $h$ (Geometric Midpoint $\sqrt{h_k h_{k+1}}$)")
    ax.set_ylabel(r"Local Empirical Order $p_{\mathrm{emp}}(h)$ (Unclipped)")
    ax.legend(loc="upper left", fontsize=7.4, framealpha=0.92)

    # --- Panel (c): Two-Regime Tidal Perturbation Scaling Across eps in [1e-5, 3e-2] ---
    ax = axes[1, 0]
    st1 = data["st1_obstruction_sweep"]
    eps_arr = np.array(st1["eps_values"], dtype=float)
    dE_bm4 = np.array(st1["dE_bm4"], dtype=float)
    dE_z4 = np.array(st1["dE_pbm4_homological_z4"], dtype=float)
    max_g_st1 = max(st1["gain_homological_z4_vs_bm4"])

    ax.loglog(eps_arr, dE_bm4, "--s", color=COL_BM4, linewidth=2.1, markersize=6.0,
              label=r"$\mathrm{BM}_4^{(7)}$ ($\mathcal{O}(\varepsilon h^4) \propto \varepsilon^1$, $10\times$/decade)")
    ax.loglog(eps_arr, dE_z4, "-D", color=COL_Z4, linewidth=2.4, markersize=6.4,
              label=rf"$\mathrm{{PBM}}_4$ (up to ${max_g_st1:.0f}\times$ suppression vs $\mathrm{{BM}}_4^{{(7)}}$)")

    # Separate reference lines for the two asymptotic regimes of PBM4:
    # 1. Quadratic floor regime eps in [1e-3, 3e-2] >= h^2 = 6.4e-3
    eps_high = eps_arr[:4]
    ref_quad = dE_z4[1] * ((eps_high / eps_high[1]) ** 2) * 0.55
    ax.loglog(eps_high, ref_quad, "-.", color="#0F172A", linewidth=1.4,
              label=rf"Quadratic floor $\mathcal{{O}}(\varepsilon^2 h^4) \propto \varepsilon^2$ (slope ${st1['slope_wrt_eps_quadratic_regime']:.2f}$)")
    # 2. Pre-asymptotic order-6 regime eps in [1e-5, 3e-4] << h^2
    eps_low = eps_arr[3:]
    ref_lin = dE_z4[-1] * ((eps_low / eps_low[-1]) ** 1) * 0.55
    ax.loglog(eps_low, ref_lin, ":", color=COL_EMERALD, linewidth=1.6,
              label=rf"Order-6 channel $\mathcal{{O}}(\varepsilon h^6) \propto \varepsilon^1$ (slope ${st1['slope_wrt_eps_order6_regime']:.2f}$)")

    ax.axvline(st1["h_squared"], color="#64748B", linestyle="--", linewidth=1.1, alpha=0.6)
    ax.text(st1["h_squared"] * 1.12, min(dE_z4) * 2.5, r"$\varepsilon_{\mathrm{cross}} \approx h^2 = 6.4\times 10^{-3}$",
            fontsize=8.6, color="#334155")
    ax.set_title(r"(c) Two-Regime Tidal Scaling $|\Delta H|/|H_0|$ vs $\varepsilon = (m_2/M_{01})(a_{\mathrm{in}}/R_{\mathrm{out}})^3$")
    ax.set_xlabel(r"Dimensionless Tidal Parameter $\varepsilon = \eta_S^3$ (at $h = 0.08, h_{\mathrm{eff}} = 0.0133$)")
    ax.set_ylabel(r"Max Relative Hamiltonian Error $\max_t |\Delta H(t)| / |H_0|$")
    ax.legend(loc="upper left", fontsize=7.8, framealpha=0.92)

    # --- Panel (d): Float128 Trajectory Position Error Crossover ||q(T) - q_ref||_2 ---
    ax = axes[1, 1]
    for key, eps_tex, rout_tex, col, mk in eps_meta:
        entry = c_f128[key]
        h_eff = np.array(entry["tau_over_s"], dtype=float)
        dq_b = np.array(entry["dq_bm4"], dtype=float)
        dq_p = np.array(entry["dq_pbm4"], dtype=float)
        ax.loglog(h_eff, dq_b, linestyle="--", marker=mk, color=col, alpha=0.55, linewidth=1.4,
                  markersize=5.0, label=rf"$\mathrm{{BM}}_4^{{(7)}}$ (${eps_tex}$)")
        ax.loglog(h_eff, dq_p, linestyle="-", marker=mk, color=col, linewidth=2.2,
                  markersize=6.2, label=rf"$\mathrm{{PBM}}_4$ (${eps_tex}$, {rout_tex})")

    yq_ref6 = c_f128["eps_1e-04"]["dq_pbm4"][0] * (h_ref / h_ref[0]) ** 6
    yq_ref4 = c_f128["eps_1e-02"]["dq_pbm4"][-1] * (h_ref / h_ref[-1]) ** 4
    ax.loglog(h_ref, yq_ref6 * 0.45, ":", color="#1E293B", linewidth=1.3, label=r"Order-6 channel $\propto \varepsilon h_{\mathrm{eff}}^6$")
    ax.loglog(h_ref, yq_ref4 * 1.8, "-.", color="#64748B", linewidth=1.3, label=r"Quadratic floor $\propto \varepsilon^2 h_{\mathrm{eff}}^4$")

    max_q_gain_1e4 = max(c_f128["eps_1e-04"]["gain_q_pbm4_vs_bm4"])
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_title(r"(d) Float128 Position Error Crossover $\|\mathbf{q}(T) - \mathbf{q}_{\mathrm{ref}}^{\mathrm{f128}}(T)\|_2$")
    ax.set_xlabel(r"Normalized Inverse Cost $h_{\mathrm{eff}} = h / s_{\mathrm{eval}}$ ($s_{\mathrm{eval}} = 6$)")
    ax.set_ylabel(r"RMS Trajectory Position Error $\|\mathbf{q} - \mathbf{q}_{\mathrm{ref}}^{\mathrm{f128}}\|_2$")
    ax.text(0.04, 0.91, f"Max $\\mathrm{{PBM}}_4$ position gain ($\\varepsilon=10^{{-4}}$): ${max_q_gain_1e4:,.0f}\\times$",
            transform=ax.transAxes, fontsize=9.0, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#EFF6FF", edgecolor=COL_Z4, alpha=0.92))
    ax.legend(loc="lower right", ncol=1, framealpha=0.92, fontsize=7.6)

    _save_and_export(fig, "fig_ch07_pbm4_quadratic_floor_crossover")


def plot_volet_B_algebraic_invariants(data: Dict[str, Any]) -> None:
    """
    Volet B (2x2):
      (a) Suppression Gain Factor of PBM4 over BM4 Across eps
      (b) 100-Orbit 3D Eccentric Secular Endurance & Continuous NAFF
      (c) Shadow Manifold Telescoping Amortization (s_eff(M) = 6 + 16/M -> 6.0)
      (d) Unclamped 3-Chart Cech 1-Cocycle & Lie-Bargmann Conservation of PBM4
    """
    fig, axes = plt.subplots(2, 2, figsize=(14.4, 10.6))
    fig.subplots_adjust(wspace=0.32, hspace=0.35)

    # --- Panel (a): Suppression Gain of PBM4 over BM4 ---
    ax = axes[0, 0]
    st1 = data["st1_obstruction_sweep"]
    eps_arr = np.array(st1["eps_values"], dtype=float)
    g_z4 = np.array(st1["gain_homological_z4_vs_bm4"], dtype=float)

    ax.loglog(eps_arr, g_z4, "-D", color=COL_Z4, linewidth=2.4, markersize=6.8,
              label=f"$\\mathrm{{PBM}}_4$ vs $\\mathrm{{BM}}_4^{{(7)}}$ (max ${max(g_z4):.1f}\\times$)")
    ax.axhline(1.0, color=COL_BM4, linestyle="-.", linewidth=1.3, label=r"$\mathrm{BM}_4^{(7)}$ baseline ($1.0\times$)")

    ax.set_title(r"(a) $\mathrm{PBM}_4$ Energy Suppression Gain Across $\varepsilon$")
    ax.set_xlabel(r"Dimensionless Tidal Parameter $\varepsilon = \eta_S^3$ ($e_{\mathrm{in}}=0.30, i=28^\circ$)")
    ax.set_ylabel(r"Energy Error Suppression Gain $\Delta E_{\mathrm{BM4}} / \Delta E_{\mathrm{PBM4}}$ ($\times$)")
    ax.legend(loc="lower left", fontsize=8.0, framealpha=0.92)

    # --- Panel (b): 100-Orbit Secular Endurance & NAFF ---
    ax = axes[0, 1]
    st3 = data["st3_secular_endurance"]
    orbs = np.array(st3["orbit_numbers"], dtype=float)
    dE_b = np.array(st3["orbit_max_dE_bm4"], dtype=float)
    dE_p = np.array(st3["orbit_max_dE_pbm4"], dtype=float)
    dL_p = np.array(st3["orbit_dL_pbm4"], dtype=float)

    ax.semilogy(orbs, _pos_or_nan(dE_b), "-", color=COL_BM4, linewidth=1.6, alpha=0.85,
                label=f"$\\mathrm{{BM}}_4^{{(7)}}$ ($\\max|\\Delta H|/|H_0| = {st3['max_dE_bm4']:.2e}$)")
    ax.semilogy(orbs, _pos_or_nan(dE_p), "-", color=COL_Z4, linewidth=2.1,
                label=f"$\\mathrm{{PBM}}_4$ ($\\max|\\Delta H|/|H_0| = {st3['max_dE_pbm4']:.2e}$, ${st3['energy_suppression_factor']:.1f}\\times$)")
    ax.semilogy(orbs, _pos_or_nan(dL_p), ":", color=COL_EMERALD, linewidth=1.5,
                label=r"Angular momentum drift $\|\Delta\mathbf{L}(t)\|_2 / \|\mathbf{L}_0\|_2$ ($\sim 10^{-15}$)")

    ax.text(
        0.04, 0.50,
        f"Continuous Laskar NAFF ($N={st3['n_samples']}$ samples, $h_{{\\mathrm{{ref}}}}=h/4$):\n"
        f"• $\\delta\\nu/\\nu_0(\\mathrm{{BM}}_4) = {st3['naff_drift_bm4']:.2e}$\n"
        f"• $\\delta\\nu/\\nu_0(\\mathrm{{PBM}}_4) = {st3['naff_drift_pbm4']:.2e}$ "
        f"($\\mathbf{{{st3['naff_gain_factor']:.1f}\\times}}$ gain)",
        transform=ax.transAxes, fontsize=8.6,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#F8FAFC", edgecolor="#94A3B8", alpha=0.94),
    )
    ax.set_title(r"(b) 100-Orbit Secular Endurance & Laskar NAFF")
    ax.set_xlabel(r"Elapsed Inner Binary Orbits ($t / T_{\mathrm{bin}}$, $e_{\mathrm{in}}=0.30, i=28^\circ$)")
    ax.set_ylabel(r"Per-Orbit Max Relative Error")
    ax.legend(loc="lower right", fontsize=7.9, framealpha=0.92)

    # --- Panel (c): Shadow Manifold Telescoping Amortization ---
    ax = axes[1, 0]
    tel = data["shadow_telescoping"]
    M_arr = np.array(tel["M_steps"], dtype=float)
    s_tel = np.array(tel["s_eff_telescoped"], dtype=float)
    spd = np.array(tel["wall_clock_speedup"], dtype=float)

    ax.semilogx(M_arr, s_tel, "-o", color=COL_Z4, linewidth=2.3, markersize=6.2,
                label=r"Shadow telescoped $\mathrm{PBM}_4$ ($s_{\mathrm{eff}}(M) = 6 + 16/M$)")
    ax.axhline(6.0, color=COL_BM4, linestyle=":", linewidth=1.5,
               label=r"Asymptotic interior cost $s_{\mathrm{eval}} = 6.0$ ($\mathrm{BM}_4^{(7)}$)")

    ax2 = ax.twinx()
    ax2.semilogx(M_arr, spd, "-.^", color=COL_EMERALD, linewidth=1.8, markersize=5.5,
                 label=r"Measured wall-clock speedup ($\times$)")
    ax2.set_ylabel(r"Wall-Clock Telescoping Speedup ($\times$)", color=COL_EMERALD)
    ax2.tick_params(axis="y", labelcolor=COL_EMERALD)
    ax2.grid(False)

    ax.set_title(r"(c) Shadow Manifold Telescoping Amortization ($\|\Delta\mathbf{q}_{\mathrm{tele}}\|_\infty = 0.0$)")
    ax.set_xlabel(r"Number of Interior Shadow Macro-Steps $M$ Between Readouts")
    ax.set_ylabel(r"Effective Apex Evaluations per Step $s_{\mathrm{eff}}(M)$")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="center right", fontsize=7.8, framealpha=0.92)

    # --- Panel (d): Unclamped 3-Chart Cech 1-Cocycle on PBM4 ---
    ax = axes[1, 1]
    st4_coc = data["st4_dynamic_gauge_and_cocycle"]["cech_cocycle_homological"]
    h_coc = np.array(st4_coc["h_steps"], dtype=float)
    coc_64 = _pos_or_nan(st4_coc["cocycle_defect_f64"])
    coc_128 = _pos_or_nan(st4_coc["cocycle_defect_f128"])
    lb_64 = _pos_or_nan(st4_coc["lie_bargmann_L_drift_f64"])
    lb_128 = _pos_or_nan(st4_coc["lie_bargmann_L_drift_f128"])

    ax.loglog(h_coc, coc_64, "-o", color=COL_Z4, linewidth=2.0, markersize=5.8,
              label=r"Float64 3-chart cocycle $\|(T_{C\to A}T_{B\to C}T_{A\to B} - \mathrm{Id})\tilde{\mathbf{z}}\|_\infty$")
    ax.loglog(h_coc, lb_64, "--^", color=COL_CYAN, linewidth=1.8, markersize=5.5,
              label=r"Float64 angular momentum drift $\|\Delta\mathbf{L}_{A\to B}\|_2$")
    ax.loglog(h_coc, coc_128, "-s", color=COL_EMERALD, linewidth=2.2, markersize=6.0,
              label=r"Float128 3-chart cocycle on $\mathrm{PBM}_4$ (unclamped)")
    ax.loglog(h_coc, lb_128, ":d", color=COL_BM4, linewidth=1.9, markersize=5.8,
              label=r"Float128 angular momentum drift $\|\Delta\mathbf{L}_{A\to B}\|_2$")

    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_ylim(1e-36, 1e-11)
    ax.set_title(r"(d) 3-Chart Čech 1-Cocycle & Lie–Bargmann Conservation")
    ax.set_xlabel(r"Macro-Step Size $h$ ($\mathcal{T}_A \to \mathcal{T}_B \to \mathcal{T}_C \to \mathcal{T}_A$)")
    ax.set_ylabel(r"Exact Unclamped Residual (Float64 & Float128)")
    ax.legend(loc="center right", fontsize=7.8, framealpha=0.92)

    _save_and_export(fig, "gen02_plate_B_algebraic_invariants", legacy_names=["fig_ch07_pbm4_vs_bm4_hierarchical_stress_test.png"])


def plot_volet_C_physical_dynamics(data: Dict[str, Any]) -> None:
    """
    Volet C (1x2):
      (a) 3D Inclined Eccentric Orbit Projection & Radial Momentum w_rad(t) = (r . w)/r
      (b) Canonical Offset ||q_shadow(t) - q_phys(t)||_2 = O(h^4) & Physical Energy Recovery
    """
    fig, axes = plt.subplots(1, 2, figsize=(14.4, 5.8))
    fig.subplots_adjust(wspace=0.32)

    phys = data["physical_3d_and_shadow_offset"]
    t_orb = np.array(phys["t_over_Tbin"], dtype=float)
    rz = np.array(phys["r_rel_z"], dtype=float)
    r_norm = np.array(phys["r_norm"], dtype=float)
    w_rad = np.array(phys["w_rad"], dtype=float)
    off_q = _pos_or_nan(phys["shadow_offset_q"])
    dE_sh = _pos_or_nan(phys["shadow_dE_raw"])
    dE_ph = _pos_or_nan(phys["phys_dE_recovered"])

    # --- Panel (a): Radial Distance r(t), Vertical Excursion z(t) & Radial Momentum w_rad(t) ---
    ax = axes[0]
    ax.plot(t_orb, r_norm, "-", color=COL_Z4, linewidth=2.2,
            label=r"Inner binary separation $r(t) = \|\mathbf{q}_1 - \mathbf{q}_0\|_2$ ($e_{\mathrm{in}}=0.30$)")
    ax.plot(t_orb, rz, "--", color=COL_CYAN, linewidth=1.8,
            label=r"Out-of-plane 3D excursion $r_z(t)$ ($i = 28^\circ$)")
    ax.plot(t_orb, w_rad, "-.", color=COL_EPS_1E2, linewidth=2.0,
            label=r"Radial momentum $w_{\mathrm{rad}}(t) = (\mathbf{r}\cdot\mathbf{w})/r \ne 0$")
    ax.axhline(0.0, color="#64748B", linestyle=":", linewidth=1.0)
    ax.set_ylim(-0.56, 1.62)
    ax.set_title(r"(a) 3D Eccentric Binary & Radial Momentum $w_{\mathrm{rad}}(t)$")
    ax.set_xlabel(r"Normalized Time $t / T_{\mathrm{bin}}$ (2 Complete Inner Binary Periods)")
    ax.set_ylabel(r"Orbital Coordinates & Canonical Radial Momentum")
    ax.legend(loc="upper right", fontsize=8.0, framealpha=0.92)

    # --- Panel (b): Shadow-to-Physical Offset ||q_tilde - q||_2 & Energy Recovery ---
    ax = axes[1]
    ax.semilogy(t_orb, dE_sh, "--", color=COL_BM4, linewidth=1.8,
                label=r"Shadow state energy oscillation $|H(\tilde{\mathbf{z}}(t)) - H_0|/|H_0| = \mathcal{O}(\varepsilon h^4)$")
    ax.semilogy(t_orb, dE_ph, "-", color=COL_Z4, linewidth=2.3,
                label=r"Recovered physical energy $|H(\mathcal{Z}_4^{\mathrm{BM4}}(\tilde{\mathbf{z}}(t))) - H_0|/|H_0| = \mathcal{O}(\varepsilon h^6)$")
    ax.semilogy(t_orb, off_q, ":", color=COL_VIOLET, linewidth=1.8,
                label=r"Canonical chart offset $\|\tilde{\mathbf{q}}(t) - \mathbf{q}(t)\|_2 = \mathcal{O}(\varepsilon h^4)$")

    ax.set_title(r"(b) Shadow Offset $\|\tilde{\mathbf{q}}(t)-\mathbf{q}(t)\|_2$ & Physical Recovery")
    ax.set_xlabel(r"Normalized Time $t / T_{\mathrm{bin}}$ ($\varepsilon = 10^{-2}, h = 0.0555$)")
    ax.set_ylabel(r"Relative Energy Error & Canonical Position Offset")
    ax.legend(loc="lower right", fontsize=8.0, framealpha=0.92)

    _save_and_export(fig, "gen02_plate_C_physical_dynamics")


def main() -> None:
    if not JSON_PATH.exists():
        raise FileNotFoundError(f"Certified JSON not found: {JSON_PATH}")
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)
    _setup_style()
    plot_volet_A_computational(data)
    print("[DONE] Canonical paper plate fig_ch07_pbm4_quadratic_floor_crossover rendered successfully.", flush=True)


if __name__ == "__main__":
    main()
