#!/usr/bin/env python3
r"""
Pure Plotter for Benchmark Generator 05: Long-Term Secular Stability,
9-Generator Lie-Bargmann Metrology, Laplace 4:2:1 Resonance & Laskar NAFF
Frequency Confinement (S1 & S2, §7.2.2).

Reads ONLY from:
    results/campaigns/generator_05_secular_stability_laplace_and_laskar_naff.json

Generates 3 standalone Laskar-segmented plates (PDF + PNG) in:
    - results/figures/
    - Artifact directory (for direct visual inspection)
    - Legacy composite figure fig_ch07_secular_brouwer_and_laplace_resonance.pdf/.png

Strict Plotting Rules:
    - Zero synthetic proxy z_act derived from energy error.
    - Zero np.clip on slopes or physical invariants.
    - Exact Float64 vs Float128 overlays from identical numerical parameters.
"""

from __future__ import annotations
import json
import shutil
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = (
    WORKSPACE_ROOT
    / "results"
    / "campaigns"
    / "generator_05_secular_stability_laplace_and_laskar_naff.json"
)
FIG_DIR = (
    WORKSPACE_ROOT / "results" / "figures"
)
ARTIFACT_DIR = WORKSPACE_ROOT / "results" / "figures"

INTEGRATOR_ORDER: List[str] = ["S2_Strang", "Y4_Yoshida", "BM4_Bare", "PBM4_Z4", "BM6_Bare"]
INTEGRATOR_COLORS: Dict[str, str] = {
    "S2_Strang": "#7f7f7f",
    "Y4_Yoshida": "#1f77b4",
    "BM4_Bare": "#ff7f0e",
    "PBM4_Z4": "#d62728",
    "BM6_Bare": "#2ca02c",
}
INTEGRATOR_LABELS: Dict[str, str] = {
    "S2_Strang": r"$\mathcal{S}_2$ (Strang, $p=2, s=1$)",
    "Y4_Yoshida": r"$\mathcal{Y}_4$ (Yoshida, $p=4, s=3$)",
    "BM4_Bare": r"$\mathrm{BM}_4^{(7)}$ Bare ($p=4, s=6$)",
    "PBM4_Z4": r"$\mathrm{PBM}_4(\mathcal{Z}_4^{\mathrm{BM4}})$ ($p_{\mathrm{eff}}=6, s=6$)",
    "BM6_Bare": r"$\mathrm{BM}_6^{(11)}$ Bare ($p=6, s=10$)",
}


def _configure_matplotlib() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.0,
        "figure.titlesize": 12.0,
        "figure.titleweight": "bold",
        "axes.grid": True,
        "grid.alpha": 0.28,
        "grid.linestyle": "--",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.12,
    })


def _pos_or_nan(arr: Any) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    return np.where(a > 0.0, a, np.nan)


def _pos_scalar_or_nan(val: float) -> float:
    v = float(val)
    return v if v > 0.0 else float("nan")


def _running_rms_or_nan(arr: Any) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    idx = np.arange(1, len(a) + 1, dtype=np.float64)
    rms = np.sqrt(np.cumsum(a * a) / idx)
    return np.where(rms > 0.0, rms, np.nan)


def _save_and_copy(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = FIG_DIR / f"{stem}.pdf"
    png_path = FIG_DIR / f"{stem}.png"
    fig.savefig(pdf_path, dpi=300)
    fig.savefig(png_path, dpi=300)
    if ARTIFACT_DIR.exists() and ARTIFACT_DIR.resolve() != FIG_DIR.resolve():
        shutil.copy2(pdf_path, ARTIFACT_DIR / f"{stem}.pdf")
        shutil.copy2(png_path, ARTIFACT_DIR / f"{stem}.png")
    print(f"[Plotter 05] Saved {pdf_path.name} and {png_path.name}")


# ==============================================================================
# VOLET A: Computational, Secular Energy & Dual-Precision Parity Plate (2x2)
# ==============================================================================
def plot_plate_A_computational(data: Dict[str, Any]) -> None:
    s1 = data["S1_laplace"]
    s2 = data["S2_solar_system"]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.6))
    fig.suptitle(
        "Generator 05 — Plate A (Computational & Secular Energy Scaling): "
        "Shadow Telescoping, Bounded Shadow Plateau & Equal-Parameter Dual Precision (§7.2.2)",
        y=0.98,
    )

    # Panel (a): S1 Jovian Laplace Resonance (N=4) Secular Energy Error (Float64 & Float128)
    ax = axes[0, 0]
    for key in INTEGRATOR_ORDER:
        entry = s1["integrators"][key]
        t_orb = np.asarray(entry["series"]["times_orbits"], dtype=np.float64)
        rel_de = _pos_or_nan(entry["series"]["rel_dE"])
        lw = 1.9 if key in ("PBM4_Z4", "BM6_Bare") else 1.2
        ax.semilogy(
            t_orb,
            rel_de,
            color=INTEGRATOR_COLORS[key],
            lw=lw,
            ls="-",
            label=f"{INTEGRATOR_LABELS[key].split(' (')[0]} F64 (max={entry['max_rel_dE']:.2e})",
        )
        if "rel_dE_f128" in entry["series"]:
            rel_de_128 = _pos_or_nan(entry["series"]["rel_dE_f128"])
            ax.semilogy(
                t_orb,
                rel_de_128,
                color=INTEGRATOR_COLORS[key],
                lw=lw,
                ls="--",
                alpha=0.88,
                label=f"{INTEGRATOR_LABELS[key].split(' (')[0]} F128 (max={entry.get('max_rel_dE_f128', np.nanmax(rel_de_128)):.2e})",
            )
    gain_s1 = s1["shadow_telescoping_gain_max"]
    gain_s1_128 = s1.get("shadow_telescoping_gain_max_f128", gain_s1)
    ax.set_title(
        f"(a) S1 Jovian Laplace Resonance ($N=4$): Gain F64={gain_s1:.1f}$\\times$, F128={gain_s1_128:.1f}$\\times$"
    )
    ax.set_xlabel("Time $t / T_{\\mathrm{Io}}$ [Io Orbits]")
    ax.set_ylabel("Relative Energy Error $|\\Delta H(t)| / |H_0|$")
    ax.legend(loc="lower right", ncol=2, framealpha=0.92, fontsize=6.4)

    # Panel (b): S2 Inner Solar System (N=5) Running RMS Energy (Float64 & Float128)
    ax = axes[0, 1]
    for key in INTEGRATOR_ORDER:
        entry = s2["integrators"][key]
        t_yr = _pos_or_nan(entry["series"]["times_yr"])
        rms_de = _pos_or_nan(entry["series"]["running_rms_dE"])
        lw = 1.9 if key in ("PBM4_Z4", "BM6_Bare") else 1.2
        ax.loglog(
            t_yr,
            rms_de,
            color=INTEGRATOR_COLORS[key],
            lw=lw,
            ls="-",
            label=f"{INTEGRATOR_LABELS[key].split(' (')[0]} F64 (RMS={entry['final_rms_dE']:.2e})",
        )
        if "running_rms_dE_f128" in entry["series"]:
            rms_de_128 = _pos_or_nan(entry["series"]["running_rms_dE_f128"])
            ax.loglog(
                t_yr,
                rms_de_128,
                color=INTEGRATOR_COLORS[key],
                lw=lw,
                ls="--",
                alpha=0.88,
                label=f"{INTEGRATOR_LABELS[key].split(' (')[0]} F128 (RMS={entry.get('final_rms_dE_f128', np.nanmax(rms_de_128)):.2e})",
            )
    gain_s2_rms = s2["shadow_telescoping_gain_rms"]
    gain_s2_rms_128 = s2.get("shadow_telescoping_gain_rms_f128", gain_s2_rms)
    ax.set_title(
        f"(b) S2 Inner Solar System ($N=5$): Gain F64={gain_s2_rms:.1f}$\\times$, F128={gain_s2_rms_128:.1f}$\\times$"
    )
    ax.set_xlabel("Integration Time $t$ [yr]")
    ax.set_ylabel("Running RMS Energy Error $\\mathrm{RMS}_{[0,t]}(|\\Delta H / H_0|)$")
    ax.legend(loc="lower right", ncol=2, framealpha=0.92, fontsize=6.4)

    # Panel (c): Equal-Parameter Dual-Precision Float64 vs Float128 Comparison (S1 & S2)
    ax = axes[1, 0]
    dp_s1 = s1["dual_precision_equal_params"]["integrators"]
    dp_s2 = s2["dual_precision_equal_params"]["integrators"]
    dp_keys = ["BM4_Bare", "PBM4_Z4", "BM6_Bare"]
    x_dp = np.arange(len(dp_keys))
    width = 0.11

    s1_f64_de = [_pos_scalar_or_nan(dp_s1[k]["float64"]["max_dE_rel"]) for k in dp_keys]
    s1_f128_de = [_pos_scalar_or_nan(dp_s1[k]["float128"]["max_dE_rel"]) for k in dp_keys]
    s2_f64_de = [_pos_scalar_or_nan(dp_s2[k]["float64"]["max_dE_rel"]) for k in dp_keys]
    s2_f128_de = [_pos_scalar_or_nan(dp_s2[k]["float128"]["max_dE_rel"]) for k in dp_keys]
    s1_f64_dl = [_pos_scalar_or_nan(dp_s1[k]["float64"]["max_dL_rel"]) for k in dp_keys]
    s1_f128_dl = [_pos_scalar_or_nan(dp_s1[k]["float128"]["max_dL_rel"]) for k in dp_keys]
    s2_f64_dl = [_pos_scalar_or_nan(dp_s2[k]["float64"]["max_dL_rel"]) for k in dp_keys]
    s2_f128_dl = [_pos_scalar_or_nan(dp_s2[k]["float128"]["max_dL_rel"]) for k in dp_keys]

    ax.bar(x_dp - 3.5 * width, s1_f64_de, width, color="#1f77b4", alpha=0.90, label="S1 Float64 $|\\Delta H|/|H_0|$")
    ax.bar(x_dp - 2.5 * width, s1_f128_de, width, color="#aec7e8", edgecolor="#1f77b4", hatch="//", label="S1 Float128 $|\\Delta H|/|H_0|$")
    ax.bar(x_dp - 1.5 * width, s2_f64_de, width, color="#d62728", alpha=0.90, label="S2 Float64 $|\\Delta H|/|H_0|$")
    ax.bar(x_dp - 0.5 * width, s2_f128_de, width, color="#ff9896", edgecolor="#d62728", hatch="//", label="S2 Float128 $|\\Delta H|/|H_0|$")
    ax.bar(x_dp + 0.5 * width, s1_f64_dl, width, color="#2ca02c", alpha=0.90, label="S1 Float64 $\\|\\Delta\\mathbf{L}\\|_2/\\|\\mathbf{L}_0\\|_2$")
    ax.bar(x_dp + 1.5 * width, s1_f128_dl, width, color="#98df8a", edgecolor="#2ca02c", hatch="//", label="S1 Float128 $\\|\\Delta\\mathbf{L}\\|_2/\\|\\mathbf{L}_0\\|_2$")
    ax.bar(x_dp + 2.5 * width, s2_f64_dl, width, color="#9467bd", alpha=0.90, label="S2 Float64 $\\|\\Delta\\mathbf{L}\\|_2/\\|\\mathbf{L}_0\\|_2$")
    ax.bar(x_dp + 3.5 * width, s2_f128_dl, width, color="#c5b0d5", edgecolor="#9467bd", hatch="//", label="S2 Float128 $\\|\\Delta\\mathbf{L}\\|_2/\\|\\mathbf{L}_0\\|_2$")
    ax.set_yscale("log")
    ax.set_xticks(x_dp)
    ax.set_xticklabels(["$\\mathrm{BM}_4^{(7)}$ Bare", "$\\mathrm{PBM}_4(\\mathcal{Z}_4^{\\mathrm{BM4}})$", "$\\mathrm{BM}_6^{(11)}$ Bare"])
    ax.axhline(2.22e-16, color="#333333", ls="--", lw=1.0, label="Float64 $\\epsilon_{\\mathrm{mach}}$")
    ax.axhline(1.0e-34, color="#2ca02c", ls=":", lw=1.1, label="Float128 $\\epsilon_{\\mathrm{mach}}$")
    ax.set_title("(c) Dual-Precision Parity (Float64 vs Float128, Equal Parameters)")
    ax.set_ylabel("Max Relative Error ($|\\Delta H|/|H_0|$ & $\\|\\Delta\\mathbf{L}\\|/\\|\\mathbf{L}_0\\|$)")
    ax.legend(loc="upper right", ncol=2, framealpha=0.92, fontsize=6.6)

    # Panel (d): Physical Secular Energy Slopes |beta_E| [yr^-1] & Cost-Accuracy Pareto
    ax = axes[1, 1]
    x_idx = np.arange(len(INTEGRATOR_ORDER))
    slopes_yr = [_pos_scalar_or_nan(abs(s2["integrators"][k]["secular_slope_per_yr"])) for k in INTEGRATOR_ORDER]
    max_de_s2 = [_pos_scalar_or_nan(s2["integrators"][k]["max_rel_dE"]) for k in INTEGRATOR_ORDER]
    bar_w = 0.34
    ax.bar(
        x_idx - 0.5 * bar_w,
        max_de_s2,
        bar_w,
        color=[INTEGRATOR_COLORS[k] for k in INTEGRATOR_ORDER],
        alpha=0.88,
        label="S2 Max Energy Error $\\max_t |\\Delta H(t)|/|H_0|$",
    )
    ax.bar(
        x_idx + 0.5 * bar_w,
        slopes_yr,
        bar_w,
        color=[INTEGRATOR_COLORS[k] for k in INTEGRATOR_ORDER],
        alpha=0.42,
        edgecolor="#111111",
        hatch="\\\\",
        label="S2 Secular Slope $|\\beta_E|$ [$\\mathrm{yr}^{-1}$] ($2\\pi\\,\\beta_{E,\\mathrm{rad}^{-1}}$)",
    )
    ax.set_yscale("log")
    ax.set_xticks(x_idx)
    ax.set_xticklabels([
        f"$\\mathcal{{S}}_2$\n($s=1$)",
        f"$\\mathcal{{Y}}_4$\n($s=3$)",
        f"$\\mathrm{{BM}}_4^{{(7)}}$\n($s=6$)",
        f"$\\mathrm{{PBM}}_4$\n($s=6$)",
        f"$\\mathrm{{BM}}_6^{{(11)}}$\n($s=10$)",
    ])
    ax.set_title("(d) S2 Energy Bound & Physical Drift Rate $|\\beta_E|$ [$\\mathrm{yr}^{-1}$]")
    ax.set_ylabel("Magnitude (Dimensionless / [$\\mathrm{yr}^{-1}$])")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.8)

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.955])
    _save_and_copy(fig, "gen05_plate_A_computational")
    plt.close(fig)

# ==============================================================================
# VOLET B: Algebraic Invariants, Lie-Bargmann & Authentic Laskar NAFF Plate (2x2)
# ==============================================================================
def plot_plate_B_algebraic_invariants(data: Dict[str, Any]) -> None:
    s1 = data["S1_laplace"]
    s2 = data["S2_solar_system"]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.6))
    fig.suptitle(
        "Generator 05 — Plate B (Algebraic Structure, 9-Generator Lie–Bargmann & Laskar NAFF): "
        "Multi-Integrator Conservation, Brouwer Diffusion & Order-3 Hanning Frequency Metrology (§7.2.2)",
        y=0.98,
    )

    x_idx = np.arange(len(INTEGRATOR_ORDER))
    xtick_names = ["$\\mathcal{S}_2$", "$\\mathcal{Y}_4$", "$\\mathrm{BM}_4^{(7)}$", "$\\mathrm{PBM}_4$", "$\\mathrm{BM}_6^{(11)}$"]

    # Panel (a): 9-Generator Lie-Bargmann Metrology Across All 5 Integrators (S1 & S2)
    ax = axes[0, 0]
    w = 0.14
    s1_dp = [_pos_scalar_or_nan(s1["integrators"][k]["lie_bargmann"]["max_dP_rel"]) for k in INTEGRATOR_ORDER]
    s1_dl = [_pos_scalar_or_nan(s1["integrators"][k]["lie_bargmann"]["max_dL_rel"]) for k in INTEGRATOR_ORDER]
    s1_dg = [_pos_scalar_or_nan(s1["integrators"][k]["lie_bargmann"]["max_dG_rel"]) for k in INTEGRATOR_ORDER]
    s2_dp = [_pos_scalar_or_nan(s2["integrators"][k]["lie_bargmann"]["max_dP_rel"]) for k in INTEGRATOR_ORDER]
    s2_dl = [_pos_scalar_or_nan(s2["integrators"][k]["lie_bargmann"]["max_dL_rel"]) for k in INTEGRATOR_ORDER]
    s2_dg = [_pos_scalar_or_nan(s2["integrators"][k]["lie_bargmann"]["max_dG_rel"]) for k in INTEGRATOR_ORDER]

    ax.bar(x_idx - 2.5 * w, s1_dp, w, color="#1f77b4", label="S1 $\\|\\Delta\\mathbf{P}\\|_2/M_{\\mathrm{tot}}$")
    ax.bar(x_idx - 1.5 * w, s1_dl, w, color="#d62728", label="S1 $\\|\\Delta\\mathbf{L}\\|_2/\\|\\mathbf{L}_0\\|_2$")
    ax.bar(x_idx - 0.5 * w, s1_dg, w, color="#2ca02c", label="S1 $\\|\\Delta\\mathbf{G}\\|_2/M_{\\mathrm{tot}}$")
    ax.bar(x_idx + 0.5 * w, s2_dp, w, color="#aec7e8", edgecolor="#1f77b4", hatch="//", label="S2 $\\|\\Delta\\mathbf{P}\\|_2/M_{\\mathrm{tot}}$")
    ax.bar(x_idx + 1.5 * w, s2_dl, w, color="#ff9896", edgecolor="#d62728", hatch="//", label="S2 $\\|\\Delta\\mathbf{L}\\|_2/\\|\\mathbf{L}_0\\|_2$")
    ax.bar(x_idx + 2.5 * w, s2_dg, w, color="#98df8a", edgecolor="#2ca02c", hatch="//", label="S2 $\\|\\Delta\\mathbf{G}\\|_2/M_{\\mathrm{tot}}$")
    ax.set_yscale("log")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(xtick_names)
    ax.axhline(2.22e-16, color="#333333", ls="--", lw=1.0, label="Float64 $\\epsilon_{\\mathrm{mach}}$")
    ax.set_title("(a) 9-Generator Lie–Bargmann Metrology Across All 5 Integrators (S1 & S2)")
    ax.set_ylabel("Max Relative Drift Over Full Horizon")
    ax.legend(loc="upper right", ncol=2, framealpha=0.92, fontsize=7.0)

    # Panel (b): Unclamped Lie-Bargmann Norms & Brouwer Stochastic t^{1/2} Diffusion on S2
    ax = axes[0, 1]
    pbm4_s2_ser = s2["integrators"]["PBM4_Z4"]["series"]
    bm6_s2_ser = s2["integrators"]["BM6_Bare"]["series"]
    t_yr = _pos_or_nan(pbm4_s2_ser["times_yr"])

    # Column 1 of legend (4 items): PBM4 (s=6) solid lines for L, G, P + Brouwer Law
    ax.loglog(
        t_yr,
        _pos_or_nan(pbm4_s2_ser["running_rms_dL_rel"]),
        color="#d62728",
        lw=1.9,
        label=r"$\mathrm{PBM}_4$: $\mathrm{RMS}(\|\Delta\mathbf{L}\|_2 / \|\mathbf{L}_0\|_2)$",
    )
    ax.loglog(
        t_yr,
        _pos_or_nan(pbm4_s2_ser["running_rms_dG_rel"]),
        color="#2ca02c",
        lw=1.9,
        label=r"$\mathrm{PBM}_4$: $\mathrm{RMS}(\|\Delta\mathbf{G}\|_2 / M_{\mathrm{tot}})$",
    )
    ax.loglog(
        t_yr,
        _running_rms_or_nan(pbm4_s2_ser["dP_norm_rel"]),
        color="#1f77b4",
        lw=1.9,
        label=r"$\mathrm{PBM}_4$: $\mathrm{RMS}(\|\Delta\mathbf{P}\|_2 / M_{\mathrm{tot}})$",
    )
    t_valid = np.asarray(pbm4_s2_ser["times_yr"][1:], dtype=np.float64)
    dl_valid = np.asarray(pbm4_s2_ser["running_rms_dL_rel"][1:], dtype=np.float64)
    c_br_l = dl_valid[-1] / np.sqrt(t_valid[-1])
    ax.loglog(
        t_valid,
        c_br_l * np.sqrt(t_valid),
        color="#111111",
        ls=":",
        lw=1.7,
        label=r"Brouwer Law $\propto t^{1/2}$ ($\sigma \propto \epsilon_{\mathrm{mach}}\sqrt{s\,n_{\mathrm{steps}}}$)",
    )

    # Column 2 of legend (3 items aligned with rows 1-3): BM6 (s=10) dashed lines of identical colors for L, G, P
    ax.loglog(
        t_yr,
        _pos_or_nan(bm6_s2_ser["running_rms_dL_rel"]),
        color="#d62728",
        lw=1.5,
        ls="--",
        alpha=0.85,
        label=r"$\mathrm{BM}_6^{(11)}$: $\mathrm{RMS}(\|\Delta\mathbf{L}\|_2 / \|\mathbf{L}_0\|_2)$",
    )
    ax.loglog(
        t_yr,
        _pos_or_nan(bm6_s2_ser["running_rms_dG_rel"]),
        color="#2ca02c",
        lw=1.5,
        ls="--",
        alpha=0.85,
        label=r"$\mathrm{BM}_6^{(11)}$: $\mathrm{RMS}(\|\Delta\mathbf{G}\|_2 / M_{\mathrm{tot}})$",
    )
    ax.loglog(
        t_yr,
        _running_rms_or_nan(bm6_s2_ser["dP_norm_rel"]),
        color="#1f77b4",
        lw=1.5,
        ls="--",
        alpha=0.85,
        label=r"$\mathrm{BM}_6^{(11)}$: $\mathrm{RMS}(\|\Delta\mathbf{P}\|_2 / M_{\mathrm{tot}})$",
    )
    ax.set_ylim(2e-22, 5e-13)
    ax.set_title("(b) S2 Unclamped Lie–Bargmann Norms & Brouwer $\\propto t^{1/2}$ Roundoff Law")
    ax.set_xlabel("Integration Time $t$ [yr]")
    ax.set_ylabel("Running RMS Relative Noether Drift")
    ax.legend(loc="upper left", ncol=2, framealpha=0.92, fontsize=6.5, columnspacing=0.9)

    # Panel (c): Authentic Laskar NAFF Frequency Error & Two-Window Shift on Mercury (z_orb, z_kepl, z_ecc)
    ax = axes[1, 0]
    w_naff = 0.24
    err_orb = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["rel_freq_err_vs_ref"]) for k in INTEGRATOR_ORDER]
    err_num_12 = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["numerical_two_window_drift_vs_ref"]) for k in INTEGRATOR_ORDER]
    raw_12 = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["raw_two_window_drift_orb"]) for k in INTEGRATOR_ORDER]

    ax.bar(
        x_idx - w_naff,
        err_orb,
        w_naff,
        color="#d62728",
        alpha=0.88,
        label=r"Fundamental Frequency Error $|\nu_{\mathrm{NAFF}} - \nu_{\mathrm{ref}}|/\nu_{\mathrm{ref}}$",
    )
    ax.bar(
        x_idx,
        err_num_12,
        w_naff,
        color="#1f77b4",
        alpha=0.88,
        label=r"Numerical 2-Window Drift $|\Delta\nu_{1\to 2} - \Delta\nu_{1\to 2}^{\mathrm{ref}}|/\nu_{\mathrm{ref}}$",
    )
    ax.bar(
        x_idx + w_naff,
        raw_12,
        w_naff,
        color="#2ca02c",
        alpha=0.55,
        edgecolor="#111111",
        hatch="//",
        label=r"Raw Physical + Numerical Shift $|\nu_{[T/2,T]} - \nu_{[0,T/2]}|/\nu_1$",
    )
    phys_shift = s2["reference_bm6_dt001"]["physical_two_window_shift_rel"]
    ax.axhline(
        phys_shift,
        color="#111111",
        ls="--",
        lw=1.2,
        label=f"Physical Secular Precession Shift $\\Delta\\nu_{{1\\to 2}}^{{\\mathrm{{phys}}}} = {phys_shift:.2e}$",
    )
    ax.set_yscale("log")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(xtick_names)
    ax.set_title("(c) S2 Authentic Laskar NAFF on Mercury's Signal $z_{\\mathrm{orb}}(t)$")
    ax.set_ylabel("Relative Frequency Deviation (Hanning $p=3$)")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.0)

    # Panel (d): Multi-Signal L2 Fidelity vs Reference
    ax = axes[1, 1]
    w_sig = 0.24
    l2_orb = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["rel_l2_z_orb_vs_ref"]) for k in INTEGRATOR_ORDER]
    l2_kepl = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["rel_l2_z_kepl_vs_ref"]) for k in INTEGRATOR_ORDER]
    l2_ecc = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["rel_l2_z_ecc_vs_ref"]) for k in INTEGRATOR_ORDER]
    ax.bar(
        x_idx - w_sig,
        l2_orb,
        w_sig,
        color="#1f77b4",
        alpha=0.88,
        label=r"Cartesian Signal $\|z_{\mathrm{orb}} - z_{\mathrm{orb}}^{\mathrm{ref}}\|_2 / \|z_{\mathrm{orb}}^{\mathrm{ref}}\|_2$",
    )
    ax.bar(
        x_idx,
        l2_kepl,
        w_sig,
        color="#ff7f0e",
        alpha=0.88,
        label=r"Keplerian Signal $\|a_1 e^{i\lambda_1} - (a_1 e^{i\lambda_1})^{\mathrm{ref}}\|_2 / \|z_{\mathrm{kepl}}^{\mathrm{ref}}\|_2$",
    )
    ax.bar(
        x_idx + w_sig,
        l2_ecc,
        w_sig,
        color="#2ca02c",
        alpha=0.88,
        label=r"Secular Eccentricity Signal $\|e_1 e^{i\varpi_1} - (e_1 e^{i\varpi_1})^{\mathrm{ref}}\|_2 / \|z_{\mathrm{ecc}}^{\mathrm{ref}}\|_2$",
    )
    ax.set_yscale("log")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(xtick_names)
    ax.set_title("(d) S2 Complex Orbital & Secular Signal Fidelity ($z_{\\mathrm{orb}}, z_{\\mathrm{kepl}}, z_{\\mathrm{ecc}}$)")
    ax.set_ylabel("Relative $L^2([0,T])$ Waveform Error vs Reference")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.0)

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.955])
    _save_and_copy(fig, "gen05_plate_B_algebraic_invariants")
    plt.close(fig)

# ==============================================================================
# VOLET C: Physical & Astrodynamical Reality Plate (2x2)
# ==============================================================================
def plot_plate_C_physical_dynamics(data: Dict[str, Any]) -> None:
    s1 = data["S1_laplace"]
    s2 = data["S2_solar_system"]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.6))
    fig.suptitle(
        "Generator 05 — Plate C (Physical & Astrodynamical Reality): "
        "Forced Laplace $4:2:1$ Libration, Mean-Longitude Ablation & Secular Hodographs (§7.2.2)",
        y=0.98,
    )

    # Panel (a): S1 Forced Laplace Libration in Mean Longitudes Across All 5 Integrators
    ax = axes[0, 0]
    for key in ["Y4_Yoshida", "BM4_Bare", "PBM4_Z4", "BM6_Bare"]:
        entry = s1["integrators"][key]
        t_orb = np.asarray(entry["series"]["times_orbits"], dtype=np.float64)
        th_m = np.asarray(entry["series"]["theta_lap_mean_deg"], dtype=np.float64)
        lw = 1.9 if key == "PBM4_Z4" else 1.15
        ls = "--" if key == "BM6_Bare" else "-"
        ax.plot(
            t_orb,
            th_m,
            color=INTEGRATOR_COLORS[key],
            lw=lw,
            ls=ls,
            label=(
                f"{INTEGRATOR_LABELS[key]} "
                f"($\\theta={entry['theta_lap_mean_median_deg']:.2f}^\\circ \\pm {entry['theta_lap_mean_half_amp_deg']:.2f}^\\circ$)"
            ),
        )
    # Inset or secondary reference showing S2_Strang resonance breaking
    s2_ent = s1["integrators"]["S2_Strang"]
    ax.axhline(
        180.0,
        color="#111111",
        ls=":",
        lw=1.2,
        label=r"Exact Laplace Equilibrium $\theta_{\mathrm{Lap}}^\star = 180^\circ$ ($\mathcal{S}_2$ breaks: $\pm 179.8^\circ$)",
    )
    ax.set_ylim(176.5, 183.5)
    ax.set_title("(a) S1 Forced Laplace $4:2:1$ Libration in Mean Longitudes $\\theta_{\\mathrm{Lap}}(t) = \\lambda_1 - 3\\lambda_2 + 2\\lambda_3$")
    ax.set_xlabel("Time $t / T_{\\mathrm{Io}}$ [Io Orbits]")
    ax.set_ylabel("Resonant Angle $\\theta_{\\mathrm{Lap}}(t)$ [deg]")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.2)

    # Panel (b): S1 Methodological Ablation: Mean vs True Polar Longitudes & Forced vs Circular ICs
    ax = axes[0, 1]
    pbm4_ser = s1["integrators"]["PBM4_Z4"]["series"]
    abl = s1["methodological_ablation"]
    t_orb = np.asarray(pbm4_ser["times_orbits"], dtype=np.float64)
    th_forced_mean = np.asarray(pbm4_ser["theta_lap_mean_deg"], dtype=np.float64)
    th_forced_true = np.asarray(pbm4_ser["theta_lap_true_deg"], dtype=np.float64)
    th_circ_mean = np.asarray(abl["unforced_theta_mean_deg"], dtype=np.float64)

    ax.plot(
        t_orb,
        th_circ_mean,
        color="#7f7f7f",
        lw=1.15,
        alpha=0.80,
        label=f"Unforced Circular ICs ($e_i(0)=0$, Mean $\\lambda_i$): $\\pm {abl['unforced_circular_mean_half_amp_deg']:.2f}^\\circ$",
    )
    ax.plot(
        t_orb,
        th_forced_true,
        color="#ff7f0e",
        lw=1.15,
        alpha=0.85,
        label=f"Forced Equilibrium + True Polar Angles $\\operatorname{{atan2}}(y_i,x_i)$: $\\pm {abl['forced_true_longitude_half_amp_deg']:.2f}^\\circ$",
    )
    ax.plot(
        t_orb,
        th_forced_mean,
        color="#d62728",
        lw=2.0,
        label=f"Forced Equilibrium + Mean Longitudes $\\lambda_i=\\varpi_i+M_i$: $\\pm {abl['forced_mean_longitude_half_amp_deg']:.2f}^\\circ$",
    )
    ax.axhline(180.0, color="#111111", ls=":", lw=1.1)
    ax.set_title("(b) S1 Methodological Ablation: Forced vs Unforced ICs & Mean vs True Longitudes")
    ax.set_xlabel("Time $t / T_{\\mathrm{Io}}$ [Io Orbits]")
    ax.set_ylabel("Resonant Angle $\\theta_{\\mathrm{Lap}}(t)$ [deg]")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.4)

    # Panel (c): S1 Forced Eccentricities e_1(t), e_2(t), e_3(t) & Resonant Phase Portrait
    ax = axes[1, 0]
    e_io = np.asarray(pbm4_ser["e_io"], dtype=np.float64)
    e_eu = np.asarray(pbm4_ser["e_eu"], dtype=np.float64)
    e_ga = np.asarray(pbm4_ser["e_ga"], dtype=np.float64)
    ax.plot(t_orb, 1e3 * e_eu, color="#1f77b4", lw=1.6, label=f"Europa $e_2(t)$ (forced $e_{{2,\\mathrm{{forced}}}} = 9.40\\times 10^{{-3}}$)")
    ax.plot(t_orb, 1e3 * e_io, color="#d62728", lw=1.6, label=f"Io $e_1(t)$ (forced $e_{{1,\\mathrm{{forced}}}} = 4.10\\times 10^{{-3}}$)")
    ax.plot(t_orb, 1e3 * e_ga, color="#2ca02c", lw=1.4, label=f"Ganymede $e_3(t)$ ($e_3(0) = 0.60\\times 10^{{-3}}$)")
    ax.axhline(9.4, color="#1f77b4", ls=":", lw=1.0)
    ax.axhline(4.1, color="#d62728", ls=":", lw=1.0)
    ax.set_title("(c) S1 Forced Eccentricity Confinement $e_{\\mathrm{Io}}(t), e_{\\mathrm{Europa}}(t), e_{\\mathrm{Ganymede}}(t)$")
    ax.set_xlabel("Time $t / T_{\\mathrm{Io}}$ [Io Orbits]")
    ax.set_ylabel("Osculating Eccentricity $10^3 \\times e_i(t)$")
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.6)

    # Panel (d): S2 Inner Solar System Secular Eccentricity Excursion (\Delta h_i, \Delta k_i) & 0 Synodic Flips
    ax = axes[1, 1]
    hodo = s2["integrators"]["PBM4_Z4"]["series"]["planetary_hodographs"]
    planet_specs = [
        ("1", "Mercury ($a_1=0.387$ AU, $e_1=0.2056$)", "#d62728"),
        ("2", "Venus ($a_2=0.723$ AU, $e_2=0.0068$)", "#ff7f0e"),
        ("3", "Earth ($a_3=1.000$ AU, $e_3=0.0167$)", "#1f77b4"),
        ("4", "Mars ($a_4=1.524$ AU, $e_4=0.0934$)", "#2ca02c"),
    ]
    for p_key, p_label, p_col in planet_specs:
        h_arr = np.asarray(hodo[p_key]["h_ecosw"], dtype=np.float64)
        k_arr = np.asarray(hodo[p_key]["k_esinw"], dtype=np.float64)
        dh = 1e5 * (h_arr - h_arr[0])
        dk = 1e5 * (k_arr - k_arr[0])
        ax.plot(dh, dk, color=p_col, lw=1.6, label=p_label)
        ax.scatter([dh[-1]], [dk[-1]], color=p_col, s=26, zorder=5)

    ax.axhline(0.0, color="#666666", lw=0.7, ls="--")
    ax.axvline(0.0, color="#666666", lw=0.7, ls="--")
    ax.set_title("(d) S2 Secular Eccentricity Drift $(\\Delta h_i, \\Delta k_i) \\times 10^5$ (0 Synodic Flips)")
    ax.set_xlabel("$10^5 \\times \\Delta(e_i \\cos\\varpi_i)$")
    ax.set_ylabel("$10^5 \\times \\Delta(e_i \\sin\\varpi_i)$")
    ax.legend(loc="best", framealpha=0.92, fontsize=7.4)

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.955])
    _save_and_copy(fig, "gen05_plate_C_physical_dynamics")
    plt.close(fig)


# ==============================================================================
# DEDICATED PUBLICATION PLATE S1: fig_ch07_s1_jovian_laplace_resonance (2x2)
# ==============================================================================
def plot_dedicated_plate_S1_laplace(data: Dict[str, Any]) -> None:
    s1 = data["S1_laplace"]

    fig = plt.figure(figsize=(14.8, 10.4))
    gs_main = fig.add_gridspec(
        2, 2,
        left=0.065, right=0.975, bottom=0.065, top=0.915,
        hspace=0.28, wspace=0.30,
    )
    fig.suptitle(
        "Protocol S1 — Resonant Jovian $4:2:1$ Laplace Satellite System ($N=4, k_{\\mathrm{sym}}=[3]$): "
        "Shadow Telescoping, Float64 Roundoff Ceiling vs Float128, Libration Hierarchy & Methodological Ablation",
        y=0.975,
        fontsize=11.5,
    )

    # ==========================================================================
    # Panel (a): S1 Secular Relative Energy Error — True Float64 (Solid) & Float128 (Dashed) Time Series
    # ==========================================================================
    ax_a = fig.add_subplot(gs_main[0, 0])
    for key in INTEGRATOR_ORDER:
        entry = s1["integrators"][key]
        t_orb = np.asarray(entry["series"]["times_orbits"], dtype=np.float64)
        rel_de_64 = _pos_or_nan(entry["series"]["rel_dE"])
        lw = 1.9 if key in ("PBM4_Z4", "BM6_Bare") else 1.2
        short_lbl = INTEGRATOR_LABELS[key].split(" (")[0]
        ax_a.semilogy(
            t_orb,
            rel_de_64,
            color=INTEGRATOR_COLORS[key],
            lw=lw,
            ls="-",
            alpha=0.90,
            label=f"{short_lbl} F64 (max={entry['max_rel_dE']:.2e})",
        )
        if "rel_dE_f128" in entry["series"]:
            rel_de_128 = _pos_or_nan(entry["series"]["rel_dE_f128"])
            max_128 = float(entry.get("max_rel_dE_f128", np.nanmax(rel_de_128)))
            ax_a.semilogy(
                t_orb,
                rel_de_128,
                color=INTEGRATOR_COLORS[key],
                lw=lw,
                ls="--",
                alpha=0.92,
                label=f"{short_lbl} F128 (max={max_128:.2e})",
            )

    gain_s1_64 = float(s1["shadow_telescoping_gain_max"])
    gain_s1_128 = float(s1.get("shadow_telescoping_gain_max_f128", gain_s1_64))
    k_s1 = s1["integrators"]["PBM4_Z4"].get("unique_dyadic_levels", [3])
    ax_a.set_title(
        f"(a) Secular Shadow Energy ($k_{{\\mathrm{{sym}}}}={k_s1}$): F64 ({gain_s1_64:.1f}$\\times$) vs F128 ({gain_s1_128:.1f}$\\times$)",
        fontsize=10.0,
    )
    ax_a.set_xlabel("Time $t / T_{\\mathrm{Io}}$ [Io Orbits]")
    ax_a.set_ylabel("Relative Energy Error $|\\Delta H(t)| / |H_0|$")
    ax_a.legend(loc="lower right", ncol=2, framealpha=0.93, fontsize=6.1, columnspacing=0.8)

    # ==========================================================================
    # Panel (b): 2-Row Split — (b1) Physical Laplace Libration (All 5 Integrators)
    #                        & (b2) Log Phase Error |theta(t) - theta_BM6(t)|
    # ==========================================================================
    gs_b = gs_main[0, 1].subgridspec(2, 1, height_ratios=[1.05, 1.05], hspace=0.14)
    ax_b1 = fig.add_subplot(gs_b[0, 0])
    ax_b2 = fig.add_subplot(gs_b[1, 0], sharex=ax_b1)

    s2_ent = s1["integrators"]["S2_Strang"]
    t_orb_s2 = np.asarray(s2_ent["series"]["times_orbits"], dtype=np.float64)
    th_s2 = np.asarray(s2_ent["series"]["theta_lap_mean_deg"], dtype=np.float64)
    ax_b1.plot(
        t_orb_s2,
        th_s2,
        color=INTEGRATOR_COLORS["S2_Strang"],
        lw=1.25,
        alpha=0.85,
        label=(
            f"$\\mathcal{{S}}_2$ "
            f"(${s2_ent['theta_lap_mean_median_deg']:.2f}^\\circ \\pm {s2_ent['theta_lap_mean_half_amp_deg']:.2f}^\\circ$)"
        ),
    )
    for key in ["Y4_Yoshida", "BM4_Bare", "PBM4_Z4", "BM6_Bare"]:
        entry = s1["integrators"][key]
        t_orb = np.asarray(entry["series"]["times_orbits"], dtype=np.float64)
        th_m = np.asarray(entry["series"]["theta_lap_mean_deg"], dtype=np.float64)
        lw = 2.1 if key == "PBM4_Z4" else 1.2
        ls = "--" if key == "BM6_Bare" else "-"
        ax_b1.plot(
            t_orb,
            th_m,
            color=INTEGRATOR_COLORS[key],
            lw=lw,
            ls=ls,
            label=(
                f"{INTEGRATOR_LABELS[key].split(' (')[0]} "
                f"(${entry['theta_lap_mean_median_deg']:.2f}^\\circ \\pm {entry['theta_lap_mean_half_amp_deg']:.2f}^\\circ$)"
            ),
        )
    ax_b1.axhline(180.0, color="#111111", ls=":", lw=1.1)
    ax_b1.set_ylabel("$\\theta_{\\mathrm{Lap}}(t)$ [deg]")
    ax_b1.set_title(
        "(b) Laplace $4:2:1$ Libration $\\theta_{\\mathrm{Lap}}(t)$ & Phase Error $|\\theta_{\\mathrm{Lap}} - \\theta_{\\mathrm{Lap}}^{\\mathrm{BM6}}|$",
        fontsize=10.0,
    )
    ax_b1.legend(loc="upper left", ncol=2, framealpha=0.92, fontsize=6.3, columnspacing=0.8)
    plt.setp(ax_b1.get_xticklabels(), visible=False)

    # (b2) Bottom: Log-scale resonant angle error |theta_Lap(t) - theta_Lap^(BM6)(t)| [deg]
    th_bm6 = np.asarray(s1["integrators"]["BM6_Bare"]["series"]["theta_lap_mean_deg"], dtype=np.float64)
    max_phase_errs: Dict[str, float] = {}
    for key in ["S2_Strang", "Y4_Yoshida", "BM4_Bare", "PBM4_Z4"]:
        entry = s1["integrators"][key]
        t_orb = np.asarray(entry["series"]["times_orbits"], dtype=np.float64)
        th_m = np.asarray(entry["series"]["theta_lap_mean_deg"], dtype=np.float64)
        err_deg = _pos_or_nan(np.abs(th_m - th_bm6))
        max_err = float(np.nanmax(err_deg))
        max_phase_errs[key] = max_err
        ls = "--" if key == "PBM4_Z4" else "-"
        lw = 1.8 if key in ("BM4_Bare", "PBM4_Z4") else 1.35
        ax_b2.semilogy(
            t_orb,
            err_deg,
            color=INTEGRATOR_COLORS[key],
            lw=lw,
            ls=ls,
            label=f"{INTEGRATOR_LABELS[key].split(' (')[0]}: max $|\\Delta\\theta| = {max_err:.2e}^\\circ$",
        )
    ax_b2.set_xlabel("Time $t / T_{\\mathrm{Io}}$ [Io Orbits]")
    ax_b2.set_ylabel("$|\\theta_{\\mathrm{Lap}} - \\theta_{\\mathrm{Lap}}^{\\mathrm{BM6}}|$ [deg]")
    ax_b2.legend(loc="lower right", ncol=2, framealpha=0.92, fontsize=6.3, columnspacing=0.8)

    # ==========================================================================
    # Panel (c): 2-Row Split — (c1) Zoomed Readout Ablation (Mean Longitudes vs True Polar Angles)
    #                        & (c2) Initialization Ablation (Unforced Circular ICs vs Forced Equilibrium + Eccentricities)
    # ==========================================================================
    gs_c = gs_main[1, 0].subgridspec(2, 1, height_ratios=[1.05, 1.05], hspace=0.14)
    ax_c1 = fig.add_subplot(gs_c[0, 0])
    ax_c2 = fig.add_subplot(gs_c[1, 0], sharex=ax_c1)

    pbm4_ser = s1["integrators"]["PBM4_Z4"]["series"]
    abl = s1["methodological_ablation"]
    t_orb = np.asarray(pbm4_ser["times_orbits"], dtype=np.float64)
    th_forced_mean = np.asarray(pbm4_ser["theta_lap_mean_deg"], dtype=np.float64)
    th_forced_true = np.asarray(pbm4_ser["theta_lap_true_deg"], dtype=np.float64)
    th_circ_mean = np.asarray(abl["unforced_theta_mean_deg"], dtype=np.float64)
    ratio_true_vs_mean = float(
        abl["forced_true_longitude_half_amp_deg"] / max(abl["forced_mean_longitude_half_amp_deg"], 1e-12)
    )
    ratio_circ_vs_forced = float(
        abl["unforced_circular_mean_half_amp_deg"] / max(abl["forced_mean_longitude_half_amp_deg"], 1e-12)
    )

    ax_c1.plot(
        t_orb,
        th_forced_true,
        color="#ff7f0e",
        lw=1.1,
        alpha=0.88,
        label=f"True Polar Angles $\\operatorname{{atan2}}(y_i,x_i)$: $\\pm {abl['forced_true_longitude_half_amp_deg']:.2f}^\\circ$ ($2e_i\\sin M_i$ ripple)",
    )
    ax_c1.plot(
        t_orb,
        th_forced_mean,
        color="#d62728",
        lw=2.2,
        label=f"Osculating Mean Longitudes $\\lambda_i=\\varpi_i+M_i$: $\\pm {abl['forced_mean_longitude_half_amp_deg']:.2f}^\\circ$ ({ratio_true_vs_mean:.1f}$\\times$ cleaner)",
    )
    ax_c1.axhline(180.0, color="#111111", ls=":", lw=1.1)
    ax_c1.set_ylim(174.8, 185.2)
    ax_c1.set_ylabel("Forced $\\theta_{\\mathrm{Lap}}$ [deg]")
    ax_c1.set_title(
        "(c) Ablation: Readout Coordinates (Top) & Initial Conditions (Bottom)",
        fontsize=10.0,
    )
    ax_c1.legend(loc="upper left", framealpha=0.92, fontsize=6.6)
    plt.setp(ax_c1.get_xticklabels(), visible=False)

    ax_c2.plot(
        t_orb,
        th_circ_mean,
        color="#555555",
        lw=1.5,
        label=f"Unforced Circular ICs ($e_i(0)=0$): $\\pm {abl['unforced_circular_mean_half_amp_deg']:.2f}^\\circ$ ({ratio_circ_vs_forced:.1f}$\\times$ larger)",
    )
    ax_c2.plot(
        t_orb,
        th_forced_mean,
        color="#d62728",
        lw=2.0,
        label=f"Forced Laplace Equilibrium: $\\pm {abl['forced_mean_longitude_half_amp_deg']:.2f}^\\circ$",
    )
    ax_c2.axhline(180.0, color="#111111", ls=":", lw=1.0)
    ax_c2.set_ylim(174.0, 249.0)
    ax_c2.set_xlabel("Time $t / T_{\\mathrm{Io}}$ [Io Orbits]")
    ax_c2.set_ylabel("Macro $\\theta_{\\mathrm{Lap}}$ [deg]")
    ax_c2.legend(loc="upper left", framealpha=0.92, fontsize=6.5)

    # Twin right axis on (c2) for Forced Eccentricities 10^3 e_i(t)
    ax_c2_r = ax_c2.twinx()
    e_io = np.asarray(pbm4_ser["e_io"], dtype=np.float64)
    e_eu = np.asarray(pbm4_ser["e_eu"], dtype=np.float64)
    e_ga = np.asarray(pbm4_ser["e_ga"], dtype=np.float64)
    ax_c2_r.plot(t_orb, 1e3 * e_eu, color="#1f77b4", lw=1.15, ls="--", alpha=0.85, label="Forced $10^3 e_{\\mathrm{Eu}}$ ($9.4$)")
    ax_c2_r.plot(t_orb, 1e3 * e_io, color="#2ca02c", lw=1.15, ls="--", alpha=0.85, label="Forced $10^3 e_{\\mathrm{Io}}$ ($4.1$)")
    ax_c2_r.plot(t_orb, 1e3 * e_ga, color="#9467bd", lw=1.15, ls=":", alpha=0.85, label="Forced $10^3 e_{\\mathrm{Ga}}$ ($0.6$)")
    ax_c2_r.set_ylim(0.0, 14.5)
    ax_c2_r.set_ylabel("Forced $10^3 e_i(t)$", fontsize=7.8, color="#1f77b4", labelpad=2)
    ax_c2_r.tick_params(axis="y", labelsize=7.2, colors="#1f77b4", pad=1)
    ax_c2_r.grid(False)
    ax_c2_r.legend(loc="center right", framealpha=0.90, fontsize=6.1)

    # ==========================================================================
    # Panel (d): S1 Energy & 9-Generator Lie-Bargmann Metrology — Paired Columns (Float64 Solid vs Float128 Hatched)
    # ==========================================================================
    ax = fig.add_subplot(gs_main[1, 1])
    x_idx = np.arange(len(INTEGRATOR_ORDER))
    xtick_names = ["$\\mathcal{S}_2$", "$\\mathcal{Y}_4$", "$\\mathrm{BM}_4^{(7)}$", "$\\mathrm{PBM}_4$", "$\\mathrm{BM}_6^{(11)}$"]
    dp_s1 = s1["dual_precision_equal_params"]["integrators"]

    # 4 observables x 2 precisions (Float64 solid vs Float128 hatched ///) under strictly equal parameters
    f64_de = [_pos_scalar_or_nan(dp_s1[k]["float64"]["max_dE_rel"]) for k in INTEGRATOR_ORDER]
    f64_dp = [_pos_scalar_or_nan(dp_s1[k]["float64"]["max_dP_rel"]) for k in INTEGRATOR_ORDER]
    f64_dl = [_pos_scalar_or_nan(dp_s1[k]["float64"]["max_dL_rel"]) for k in INTEGRATOR_ORDER]
    f64_dg = [_pos_scalar_or_nan(dp_s1[k]["float64"]["max_dG_rel"]) for k in INTEGRATOR_ORDER]

    f128_de = [_pos_scalar_or_nan(dp_s1[k]["float128"]["max_dE_rel"]) for k in INTEGRATOR_ORDER]
    f128_dp = [_pos_scalar_or_nan(dp_s1[k]["float128"]["max_dP_rel"]) for k in INTEGRATOR_ORDER]
    f128_dl = [_pos_scalar_or_nan(dp_s1[k]["float128"]["max_dL_rel"]) for k in INTEGRATOR_ORDER]
    f128_dg = [_pos_scalar_or_nan(dp_s1[k]["float128"]["max_dG_rel"]) for k in INTEGRATOR_ORDER]

    w = 0.092
    # Pair 1: Energy |dH|/|H0| (Purple)
    ax.bar(x_idx - 3.65 * w, f64_de, w, color="#9467bd", alpha=0.92, edgecolor="#111111", lw=0.5)
    ax.bar(x_idx - 2.65 * w, f128_de, w, color="#9467bd", alpha=0.35, edgecolor="#9467bd", hatch="///", lw=1.0)
    # Pair 2: Linear Momentum ||dP||_2 / M_tot (Blue)
    ax.bar(x_idx - 1.45 * w, f64_dp, w, color="#1f77b4", alpha=0.92, edgecolor="#111111", lw=0.5)
    ax.bar(x_idx - 0.45 * w, f128_dp, w, color="#1f77b4", alpha=0.35, edgecolor="#1f77b4", hatch="///", lw=1.0)
    # Pair 3: Angular Momentum ||dL||_2 / ||L0||_2 (Red)
    ax.bar(x_idx + 0.75 * w, f64_dl, w, color="#d62728", alpha=0.92, edgecolor="#111111", lw=0.5)
    ax.bar(x_idx + 1.75 * w, f128_dl, w, color="#d62728", alpha=0.35, edgecolor="#d62728", hatch="///", lw=1.0)
    # Pair 4: Galilean Boost ||dG||_2 / M_tot (Green)
    ax.bar(x_idx + 2.95 * w, f64_dg, w, color="#2ca02c", alpha=0.92, edgecolor="#111111", lw=0.5)
    ax.bar(x_idx + 3.95 * w, f128_dg, w, color="#2ca02c", alpha=0.35, edgecolor="#2ca02c", hatch="///", lw=1.0)

    ax.set_yscale("log")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(xtick_names)
    ax.axhline(2.22e-16, color="#333333", ls="--", lw=1.0)
    ax.axhline(1.0e-34, color="#2ca02c", ls=":", lw=1.1)
    ax.set_ylim(1e-41, 1e-1)
    ax.set_title(
        "(d) Energy Truncation vs Lie–Bargmann Conservation (F64 vs F128)",
        fontsize=10.0,
    )
    ax.set_ylabel("Max Relative Error / Drift (Float64 & Float128)")

    legend_handles = [
        Patch(facecolor="#9467bd", alpha=0.92, edgecolor="#111111", label=r"Float64 $|\Delta H|/|H_0|$"),
        Patch(facecolor="#1f77b4", alpha=0.92, edgecolor="#111111", label=r"Float64 $\|\Delta\mathbf{P}\|_2/M_{\mathrm{tot}}$"),
        Patch(facecolor="#d62728", alpha=0.92, edgecolor="#111111", label=r"Float64 $\|\Delta\mathbf{L}\|_2/\|\mathbf{L}_0\|_2$"),
        Patch(facecolor="#2ca02c", alpha=0.92, edgecolor="#111111", label=r"Float64 $\|\Delta\mathbf{G}\|_2/M_{\mathrm{tot}}$"),
        Line2D([0], [0], color="#333333", ls="--", lw=1.0, label=r"Float64 $\epsilon_{\mathrm{mach}}$ ($2.2\times 10^{-16}$)"),
        Patch(facecolor="#9467bd", alpha=0.35, edgecolor="#9467bd", hatch="///", label=r"Float128 $|\Delta H|/|H_0|$"),
        Patch(facecolor="#1f77b4", alpha=0.35, edgecolor="#1f77b4", hatch="///", label=r"Float128 $\|\Delta\mathbf{P}\|_2/M_{\mathrm{tot}}$"),
        Patch(facecolor="#d62728", alpha=0.35, edgecolor="#d62728", hatch="///", label=r"Float128 $\|\Delta\mathbf{L}\|_2/\|\mathbf{L}_0\|_2$"),
        Patch(facecolor="#2ca02c", alpha=0.35, edgecolor="#2ca02c", hatch="///", label=r"Float128 $\|\Delta\mathbf{G}\|_2/M_{\mathrm{tot}}$"),
        Line2D([0], [0], color="#2ca02c", ls=":", lw=1.1, label=r"Float128 $\epsilon_{\mathrm{mach}}$ ($10^{-34}$)"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", ncol=2, framealpha=0.92, fontsize=6.4, columnspacing=1.0)

    _save_and_copy(fig, "fig_ch07_s1_jovian_laplace_resonance")
    plt.close(fig)


# ==============================================================================
# DEDICATED PUBLICATION PLATE S2: fig_ch07_s2_solar_system_brouwer_naff (2x2)
# ==============================================================================
def plot_dedicated_plate_S2_solar_system(data: Dict[str, Any]) -> None:
    s2 = data["S2_solar_system"]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.6))
    fig.suptitle(
        "Protocol S2 — Inner Solar System ($N=5$: Sun, Mercury, Venus, Earth, Mars): "
        "Bounded Shadow Plateau ($\\propto t^0$), Unbiased Brouwer Noether Diffusion ($\\propto t^{1/2}$) & Authentic Laskar NAFF",
        y=0.98,
    )

    # Panel (a): S2 Running RMS Energy Error — True Float64 (Solid) & Float128 (Dashed) Time Series
    ax = axes[0, 0]
    for key in INTEGRATOR_ORDER:
        entry = s2["integrators"][key]
        t_yr = _pos_or_nan(entry["series"]["times_yr"])
        rms_de_64 = _running_rms_or_nan(entry["series"]["abs_dE"])
        lw = 2.0 if key in ("PBM4_Z4", "BM6_Bare") else 1.25
        short_lbl = INTEGRATOR_LABELS[key].split(" (")[0]
        ax.loglog(
            t_yr,
            rms_de_64,
            color=INTEGRATOR_COLORS[key],
            lw=lw,
            ls="-",
            alpha=0.90,
            label=f"{short_lbl} F64 (RMS={entry['final_rms_dE']:.2e})",
        )
        if "abs_dE_f128" in entry["series"]:
            rms_de_128 = _running_rms_or_nan(entry["series"]["abs_dE_f128"])
            final_rms_128 = float(entry.get("final_rms_dE_f128", rms_de_128[-1]))
            ax.loglog(
                t_yr,
                rms_de_128,
                color=INTEGRATOR_COLORS[key],
                lw=lw,
                ls="--",
                alpha=0.92,
                label=f"{short_lbl} F128 (RMS={final_rms_128:.2e})",
            )
    gain_s2_rms_64 = float(s2["shadow_telescoping_gain_rms"])
    gain_s2_rms_128 = float(s2.get("shadow_telescoping_gain_rms_f128", gain_s2_rms_64))
    k_s2 = s2["integrators"]["PBM4_Z4"].get("unique_dyadic_levels", [2, 3])
    n_junc_s2 = s2["integrators"]["PBM4_Z4"].get("n_dyadic_junctions", 0)
    ax.set_title(
        f"(a) Running RMS Shadow Energy ($k_{{\\mathrm{{sym}}}} \\in {k_s2}$, {n_junc_s2} Junc.): F64={gain_s2_rms_64:.1f}$\\times$, F128={gain_s2_rms_128:.1f}$\\times$",
        fontsize=9.8,
    )
    ax.set_xlabel("Integration Time $t$ [yr]")
    ax.set_ylabel("Running RMS Energy Error $\\mathrm{RMS}_{[0,t]}(|\\Delta H / H_0|)$")
    ax.legend(loc="lower right", ncol=2, framealpha=0.92, fontsize=6.1, columnspacing=0.8)

    # Panel (b): Unbiased Brouwer Stochastic Roundoff Diffusion (~ t^{1/2}) on Exact Lie-Bargmann Noether Invariants
    ax = axes[0, 1]
    pbm4_s2_ser = s2["integrators"]["PBM4_Z4"]["series"]
    bm6_s2_ser = s2["integrators"]["BM6_Bare"]["series"]
    t_yr = _pos_or_nan(pbm4_s2_ser["times_yr"])

    # Column 1 of legend (4 items): PBM4 (s=6) solid lines for L (red), G (green), P (blue) + Brouwer Law
    ax.loglog(
        t_yr,
        _pos_or_nan(pbm4_s2_ser["running_rms_dL_rel"]),
        color="#d62728",
        lw=1.9,
        label=r"$\mathrm{PBM}_4$: $\mathrm{RMS}(\|\Delta\mathbf{L}\|_2 / \|\mathbf{L}_0\|_2)$",
    )
    ax.loglog(
        t_yr,
        _pos_or_nan(pbm4_s2_ser["running_rms_dG_rel"]),
        color="#2ca02c",
        lw=1.9,
        label=r"$\mathrm{PBM}_4$: $\mathrm{RMS}(\|\Delta\mathbf{G}\|_2 / M_{\mathrm{tot}})$",
    )
    ax.loglog(
        t_yr,
        _running_rms_or_nan(pbm4_s2_ser["dP_norm_rel"]),
        color="#1f77b4",
        lw=1.9,
        label=r"$\mathrm{PBM}_4$: $\mathrm{RMS}(\|\Delta\mathbf{P}\|_2 / M_{\mathrm{tot}})$",
    )
    t_valid = np.asarray(pbm4_s2_ser["times_yr"][1:], dtype=np.float64)
    dl_valid = np.asarray(pbm4_s2_ser["running_rms_dL_rel"][1:], dtype=np.float64)
    c_br_l = dl_valid[-1] / np.sqrt(t_valid[-1])
    ax.loglog(
        t_valid,
        c_br_l * np.sqrt(t_valid),
        color="#111111",
        ls=":",
        lw=1.7,
        label=r"Brouwer Law $\propto t^{1/2}$ ($\sigma \propto \epsilon_{\mathrm{mach}}\sqrt{s\,n_{\mathrm{steps}}}$)",
    )

    # Column 2 of legend (3 items aligned with rows 1-3): BM6 (s=10) dashed lines of the exact same colors for L, G, P
    ax.loglog(
        t_yr,
        _pos_or_nan(bm6_s2_ser["running_rms_dL_rel"]),
        color="#d62728",
        lw=1.5,
        ls="--",
        alpha=0.85,
        label=r"$\mathrm{BM}_6^{(11)}$: $\mathrm{RMS}(\|\Delta\mathbf{L}\|_2 / \|\mathbf{L}_0\|_2)$",
    )
    ax.loglog(
        t_yr,
        _pos_or_nan(bm6_s2_ser["running_rms_dG_rel"]),
        color="#2ca02c",
        lw=1.5,
        ls="--",
        alpha=0.85,
        label=r"$\mathrm{BM}_6^{(11)}$: $\mathrm{RMS}(\|\Delta\mathbf{G}\|_2 / M_{\mathrm{tot}})$",
    )
    ax.loglog(
        t_yr,
        _running_rms_or_nan(bm6_s2_ser["dP_norm_rel"]),
        color="#1f77b4",
        lw=1.5,
        ls="--",
        alpha=0.85,
        label=r"$\mathrm{BM}_6^{(11)}$: $\mathrm{RMS}(\|\Delta\mathbf{P}\|_2 / M_{\mathrm{tot}})$",
    )
    ax.set_ylim(2e-22, 5e-13)
    ax.set_title("(b) Exact Lie–Bargmann Invariants & Brouwer $\\propto t^{1/2}$ Diffusion", fontsize=9.8)
    ax.set_xlabel("Integration Time $t$ [yr]")
    ax.set_ylabel("Running RMS Relative Noether Drift")
    ax.legend(loc="upper left", ncol=2, framealpha=0.92, fontsize=6.5, columnspacing=0.9)

    # Panel (c): Authentic Order-3 Hanning Laskar NAFF Frequency Metrology on Mercury's Orbital Signal
    ax = axes[1, 0]
    x_idx = np.arange(len(INTEGRATOR_ORDER))
    xtick_names = ["$\\mathcal{S}_2$", "$\\mathcal{Y}_4$", "$\\mathrm{BM}_4^{(7)}$", "$\\mathrm{PBM}_4$", "$\\mathrm{BM}_6^{(11)}$"]
    w_naff = 0.25
    err_orb = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["rel_freq_err_vs_ref"]) for k in INTEGRATOR_ORDER]
    err_num_12 = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["numerical_two_window_drift_vs_ref"]) for k in INTEGRATOR_ORDER]
    l2_orb = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["rel_l2_z_orb_vs_ref"]) for k in INTEGRATOR_ORDER]

    ax.bar(
        x_idx - w_naff,
        err_orb,
        w_naff,
        color="#d62728",
        alpha=0.90,
        label=r"Fundamental Frequency Error $|\nu_{\mathrm{NAFF}} - \nu_{\mathrm{ref}}|/\nu_{\mathrm{ref}}$ ($z_{\mathrm{orb}}$)",
    )
    ax.bar(
        x_idx,
        err_num_12,
        w_naff,
        color="#1f77b4",
        alpha=0.90,
        label=r"Numerical 2-Window Drift $|\Delta\nu_{1\to 2} - \Delta\nu_{1\to 2}^{\mathrm{ref}}|/\nu_{\mathrm{ref}}$",
    )
    ax.bar(
        x_idx + w_naff,
        l2_orb,
        w_naff,
        color="#2ca02c",
        alpha=0.72,
        edgecolor="#111111",
        hatch="//",
        label=r"Complex Waveform $L^2$ Error $\|z_{\mathrm{orb}} - z_{\mathrm{orb}}^{\mathrm{ref}}\|_2 / \|z_{\mathrm{orb}}^{\mathrm{ref}}\|_2$",
    )
    phys_shift = s2["reference_bm6_dt001"]["physical_two_window_shift_rel"]
    ax.axhline(
        phys_shift,
        color="#111111",
        ls="--",
        lw=1.2,
        label=f"Physical Secular Precession $\\Delta\\nu_{{1\\to 2}}^{{\\mathrm{{phys}}}} = {phys_shift:.2e}$",
    )
    ax.set_yscale("log")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(xtick_names)
    ax.set_title("(c) Authentic Laskar NAFF Frequency & Waveform Metrology on Mercury ($4{,}001$ Samples)")
    ax.set_ylabel("Relative Frequency / Waveform Error (Hanning $p=3$)")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=6.9)

    # Panel (d): Secular Planetary Eccentricity Hodographs (\Delta h_i, \Delta k_i) x 10^5 (0 Synodic Flips)
    ax = axes[1, 1]
    hodo = s2["integrators"]["PBM4_Z4"]["series"]["planetary_hodographs"]
    planet_specs = [
        ("1", "Mercury ($a_1=0.387$ AU, $e_1=0.2056$)", "#d62728"),
        ("2", "Venus ($a_2=0.723$ AU, $e_2=0.0068$)", "#ff7f0e"),
        ("3", "Earth ($a_3=1.000$ AU, $e_3=0.0167$)", "#1f77b4"),
        ("4", "Mars ($a_4=1.524$ AU, $e_4=0.0934$)", "#2ca02c"),
    ]
    for p_key, p_label, p_col in planet_specs:
        h_arr = np.asarray(hodo[p_key]["h_ecosw"], dtype=np.float64)
        k_arr = np.asarray(hodo[p_key]["k_esinw"], dtype=np.float64)
        dh = 1e5 * (h_arr - h_arr[0])
        dk = 1e5 * (k_arr - k_arr[0])
        ax.plot(dh, dk, color=p_col, lw=1.6, label=p_label)
        ax.scatter([dh[-1]], [dk[-1]], color=p_col, s=26, zorder=5)

    ax.axhline(0.0, color="#666666", lw=0.7, ls="--")
    ax.axvline(0.0, color="#666666", lw=0.7, ls="--")
    ax.set_title("(d) Secular Eccentricity Excursion $(\\Delta h_i, \\Delta k_i) \\times 10^5$ ($0$ Synodic Chart Flips)")
    ax.set_xlabel("$10^5 \\times \\Delta(e_i \\cos\\varpi_i)$")
    ax.set_ylabel("$10^5 \\times \\Delta(e_i \\sin\\varpi_i)$")
    ax.legend(loc="best", framealpha=0.92, fontsize=7.4)

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.955])
    _save_and_copy(fig, "fig_ch07_s2_solar_system_brouwer_naff")
    plt.close(fig)


# ==============================================================================
# COMPOSITE LEGACY FIGURE: fig_ch07_secular_brouwer_and_laplace_resonance
# ==============================================================================
def plot_composite_ch07_figure(data: Dict[str, Any]) -> None:
    s1 = data["S1_laplace"]
    s2 = data["S2_solar_system"]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.6))
    fig.suptitle(
        "Section 7.2.2 — Secular Brouwer Diffusion, Authentic Laskar NAFF & Forced Laplace $4:2:1$ Resonance",
        y=0.98,
    )

    # (a) S2 Running RMS Energy & Bounded Shadow Plateau
    ax = axes[0, 0]
    for key in INTEGRATOR_ORDER:
        entry = s2["integrators"][key]
        t_yr = _pos_or_nan(entry["series"]["times_yr"])
        rms_de = _pos_or_nan(entry["series"]["running_rms_dE"])
        ax.loglog(
            t_yr,
            rms_de,
            color=INTEGRATOR_COLORS[key],
            lw=2.0 if key in ("PBM4_Z4", "BM6_Bare") else 1.25,
            label=f"{INTEGRATOR_LABELS[key]} (RMS={entry['final_rms_dE']:.2e})",
        )
    ax.set_title(f"(a) S2 Secular Energy RMS & Shadow Telescoping ({s2['shadow_telescoping_gain_rms']:.1f}$\\times$ Gain)")
    ax.set_xlabel("Time $t$ [yr]")
    ax.set_ylabel("Running RMS Energy Error $\\mathrm{RMS}_{[0,t]}(|\\Delta H / H_0|)$")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=7.4)

    # (b) S2 Authentic Laskar NAFF Frequency Error
    ax = axes[0, 1]
    x_idx = np.arange(len(INTEGRATOR_ORDER))
    xtick_names = ["$\\mathcal{S}_2$", "$\\mathcal{Y}_4$", "$\\mathrm{BM}_4^{(7)}$", "$\\mathrm{PBM}_4$", "$\\mathrm{BM}_6^{(11)}$"]
    err_orb = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["rel_freq_err_vs_ref"]) for k in INTEGRATOR_ORDER]
    err_num_12 = [_pos_scalar_or_nan(s2["integrators"][k]["naff"]["numerical_two_window_drift_vs_ref"]) for k in INTEGRATOR_ORDER]
    ax.bar(x_idx - 0.17, err_orb, 0.34, color="#d62728", alpha=0.88, label=r"$|\nu_{\mathrm{NAFF}} - \nu_{\mathrm{ref}}|/\nu_{\mathrm{ref}}$")
    ax.bar(x_idx + 0.17, err_num_12, 0.34, color="#1f77b4", alpha=0.88, label=r"$|\Delta\nu_{1\to 2} - \Delta\nu_{1\to 2}^{\mathrm{ref}}|/\nu_{\mathrm{ref}}$")
    ax.set_yscale("log")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(xtick_names)
    ax.set_title("(b) S2 Authentic Laskar NAFF on Mercury's Orbital Signal $z_{\\mathrm{orb}}(t)$")
    ax.set_ylabel("Relative Frequency Error (Hanning $p=3$)")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.6)

    # (c) S1 Forced Laplace Libration
    ax = axes[1, 0]
    for key in ["Y4_Yoshida", "BM4_Bare", "PBM4_Z4", "BM6_Bare"]:
        entry = s1["integrators"][key]
        t_orb = np.asarray(entry["series"]["times_orbits"], dtype=np.float64)
        th_m = np.asarray(entry["series"]["theta_lap_mean_deg"], dtype=np.float64)
        ax.plot(
            t_orb,
            th_m,
            color=INTEGRATOR_COLORS[key],
            lw=1.9 if key == "PBM4_Z4" else 1.15,
            ls="--" if key == "BM6_Bare" else "-",
            label=f"{INTEGRATOR_LABELS[key]} ($\\pm {entry['theta_lap_mean_half_amp_deg']:.2f}^\\circ$)",
        )
    ax.axhline(180.0, color="#111111", ls=":", lw=1.2)
    ax.set_ylim(176.5, 183.5)
    ax.set_title("(c) S1 Forced Laplace $4:2:1$ Libration in Mean Longitudes $\\theta_{\\mathrm{Lap}}(t)$")
    ax.set_xlabel("Time $t / T_{\\mathrm{Io}}$ [Io Orbits]")
    ax.set_ylabel("Resonant Angle $\\theta_{\\mathrm{Lap}}(t)$ [deg]")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.4)

    # (d) 9-Generator Lie-Bargmann Metrology Across All 5 Integrators
    ax = axes[1, 1]
    w = 0.25
    s1_dl = [_pos_scalar_or_nan(s1["integrators"][k]["lie_bargmann"]["max_dL_rel"]) for k in INTEGRATOR_ORDER]
    s2_dl = [_pos_scalar_or_nan(s2["integrators"][k]["lie_bargmann"]["max_dL_rel"]) for k in INTEGRATOR_ORDER]
    s2_dg = [_pos_scalar_or_nan(s2["integrators"][k]["lie_bargmann"]["max_dG_rel"]) for k in INTEGRATOR_ORDER]
    ax.bar(x_idx - w, s1_dl, w, color="#1f77b4", label="S1 $\\|\\Delta\\mathbf{L}\\|_2/\\|\\mathbf{L}_0\\|_2$")
    ax.bar(x_idx, s2_dl, w, color="#d62728", label="S2 $\\|\\Delta\\mathbf{L}\\|_2/\\|\\mathbf{L}_0\\|_2$")
    ax.bar(x_idx + w, s2_dg, w, color="#2ca02c", label="S2 $\\|\\Delta\\mathbf{G}\\|_2/M_{\\mathrm{tot}}$")
    ax.set_yscale("log")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(xtick_names)
    ax.axhline(2.22e-16, color="#333333", ls="--", lw=1.0, label="Float64 $\\epsilon_{\\mathrm{mach}}$")
    ax.set_title("(d) Multi-Integrator Lie–Bargmann Conservation (`S1` & `S2`)")
    ax.set_ylabel("Max Relative Drift")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.4)

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.955])
    _save_and_copy(fig, "fig_ch07_secular_brouwer_and_laplace_resonance")
    plt.close(fig)


def main() -> None:
    if not JSON_PATH.exists():
        raise FileNotFoundError(f"Certified JSON not found: {JSON_PATH}")
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)

    _configure_matplotlib()
    plot_dedicated_plate_S1_laplace(data)
    plot_dedicated_plate_S2_solar_system(data)
    print("[Plotter 05] Canonical paper plates S1 and S2 generated successfully.")


if __name__ == "__main__":
    main()
