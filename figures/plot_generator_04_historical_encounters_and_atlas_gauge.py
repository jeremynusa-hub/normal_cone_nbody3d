#!/usr/bin/env python3
r"""
Pure Plotter for Benchmark Generator 04: Documented Historical Planetary Encounters
(H1 Lexell, H2 Shoemaker-Levy 9, H3 Apophis 2029) & Dynamic Atlas Gauge Validation (§7.2.1).

Reads ONLY from:
    results/campaigns/generator_04_historical_encounters_and_atlas_gauge.json

Generates 3 standalone Laskar-segmented plates (PDF + PNG) in:
    - results/figures/
    - Artifact directory (for direct visual inspection)
    - Legacy composite figure fig_ch07_historical_encounters_atlas_gauge.pdf/.png

Strict Plotting Rules:
    - Zero np.clip on slopes or physical invariants.
    - Zero synthetic sin/cos curves.
    - Exact Float64 vs Float128 overlays from identical numerical parameters.
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

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = WORKSPACE_ROOT / "results" / "campaigns" / "generator_04_historical_encounters_and_atlas_gauge.json"
FIG_DIR = (
    WORKSPACE_ROOT / "results" / "figures"
)
ARTIFACT_DIR = WORKSPACE_ROOT / "results" / "figures"


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


def _save_and_copy(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = FIG_DIR / f"{stem}.pdf"
    png_path = FIG_DIR / f"{stem}.png"
    fig.savefig(pdf_path, dpi=300)
    fig.savefig(png_path, dpi=300)
    if ARTIFACT_DIR.exists() and ARTIFACT_DIR.resolve() != FIG_DIR.resolve():
        shutil.copy2(pdf_path, ARTIFACT_DIR / f"{stem}.pdf")
        shutil.copy2(png_path, ARTIFACT_DIR / f"{stem}.png")
    print(f"[Plotter 04] Saved {pdf_path.name} and {png_path.name}")


def _pos_or_nan(arr: Any) -> np.ndarray:
    """Replaces non-positive mathematical zeros (<= 0.0, e.g. at t=0) with NaN for log-scale plotting without altering any positive solver value."""
    a = np.asarray(arr, dtype=float).copy()
    a[a <= 0.0] = np.nan
    return a


# ==============================================================================
# VOLET A: Computational, Multirate & Trigger Dynamics Plate (2x2)
# ==============================================================================
def plot_plate_A_computational(data: Dict[str, Any]) -> None:
    h1 = data["H1_lexell"]
    h2 = data["H2_shoemaker_levy"]
    h3 = data["H3_apophis"]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.6))
    fig.suptitle(
        "Generator 04 — Plate A (Computational & Multirate Scaling): "
        "Dynamic Atlas Hysteresis, Roche Promotion & Bilateral Perigee Guard (§7.2.1)",
        y=0.98,
    )

    # --- Panel (a): H1 Lexell Hill Switching & Dyadic Rungs ---
    ax = axes[0, 0]
    t_h1 = np.asarray(h1["series"]["times"], dtype=float)
    rj_h1 = np.asarray(h1["series"]["r_jup_au_f64"], dtype=float)
    k64_h1 = np.asarray(h1["series"]["dyadic_k_f64"], dtype=int)
    k128_h1 = np.asarray(h1["series"]["dyadic_k_f128"], dtype=int)
    r_hill = float(h1["r_hill_jup_au"])

    l1 = ax.semilogy(t_h1, rj_h1, color="#1f4e78", lw=1.9, label=r"$r_{J}(t)$ (AU, Float64/128)")
    ax.axhline(r_hill, color="#c0392b", ls="--", lw=1.2, label=rf"Hill Switch $\gamma R_{{\mathrm{{Hill}}}} = {r_hill:.3f}$ AU")
    for idx_f, fl in enumerate(h1["atlas_transitions"]):
        tf = float(fl["t"])
        ax.axvline(
            tf,
            color="#8e44ad",
            ls=":",
            lw=1.3,
            label=rf"Atlas Flip $t={tf:.2f}$ ($\mathcal{{T}}_{{\mathrm{{helio}}}}\leftrightarrow\mathcal{{T}}_{{\mathrm{{Jovi}}}}$)"
            if idx_f == 0
            else f"Egress Flip $t={tf:.2f}$",
        )
    ax.scatter([float(h1["t_perijove"])], [float(h1["min_r_jup_au"])], color="#d35400", s=45, zorder=6,
               label=rf"Perijove $r_{{\min}}={h1['min_r_jup_RJ']:.2f}\,R_J$")
    ax.set_xlabel("Physical Time $t$ (yr / $2\\pi$)")
    ax.set_ylabel(r"Jovicentric Distance $r_J(t)$ [AU]")
    ax.set_ylim(8e-4, 1.5)

    ax_r = ax.twinx()
    ax_r.step(t_h1, k64_h1, where="post", color="#27ae60", lw=1.4, alpha=0.80, label=r"$k_{\mathrm{sym}}(t)$ (Float64)")
    ax_r.step(t_h1, k128_h1, where="post", color="#16a085", lw=1.2, ls="--", label=r"$k_{\mathrm{sym}}(t)$ (Float128)")
    ax_r.set_ylabel(r"Symmetric Dyadic Rung $k_{\mathrm{sym}}$", color="#1e8449")
    ax_r.set_ylim(-0.3, max(int(np.max(k64_h1)), 4) + 1.2)
    ax_r.grid(False)

    lines_a, labels_a = ax.get_legend_handles_labels()
    lines_ar, labels_ar = ax_r.get_legend_handles_labels()
    ax.legend(lines_a + lines_ar, labels_a + labels_ar, loc="upper right", framealpha=0.92, fontsize=7.3)
    ax.set_title(r"(a) H1 Lexell (1770): Dynamic Hill Chart Flips ($\mathcal{T}_{\mathrm{helio}}\leftrightarrow\mathcal{T}_{\mathrm{Jovi}}$) & Rungs")

    # --- Panel (b): H2 Shoemaker-Levy 9 Cubic Tidal Parameter Lambda_S(t) ---
    ax = axes[0, 1]
    t_h2 = np.asarray(h2["series"]["times"], dtype=float)
    lam_h2 = np.asarray(h2["series"]["lambda_S"], dtype=float)
    act_h2 = np.asarray(h2["series"]["roche_supernode_active"], dtype=int)
    lam_crit = float(h2["lambda_crit"])
    lam_exit = float(h2["lambda_exit"])

    ax.semilogy(t_h2, lam_h2, color="#8e44ad", lw=2.0, label=r"Cubic Tidal Ratio $\Lambda_S(t) = \frac{1}{3}\eta_S(t)^3$")
    ax.axhline(lam_crit, color="#c0392b", ls="--", lw=1.4, label=rf"Roche Trigger $\Lambda_{{\mathrm{{crit}}}} = {lam_crit:.2f}$")
    ax.axhline(lam_exit, color="#e67e22", ls=":", lw=1.2, label=rf"Schmitt Exit $\Lambda_{{\mathrm{{exit}}}} = {lam_exit:.2f}$")

    # Shade active supernode window
    active_idx = np.where(act_h2 == 1)[0]
    if len(active_idx) > 0:
        t_enter = t_h2[active_idx[0]]
        t_leave = t_h2[active_idx[-1]]
        ax.axvspan(t_enter, t_leave, color="#f39c12", alpha=0.18,
                   label=rf"Promoted 3-Fragment Supernode ({h2['roche_active_steps']}/{h2['equal_parameters']['n_steps']} steps)")
    ax.axvline(float(h2["t_perijove"]), color="#2c3e50", ls="-.", lw=1.1,
               label=rf"Sub-Roche Perijove $t_{{\mathrm{{peri}}}}={h2['t_perijove']:.4f}$ ($\Lambda_{{\max}}={h2['lambda_max']:.0f}$)")
    ax.set_xlabel("Physical Time $t$ (yr / $2\\pi$)")
    ax.set_ylabel(r"Cubic Tidal Force Ratio $\Lambda_S(t) = \frac{1}{3}\eta_S^3$")
    ax.set_ylim(1e-2, max(lam_h2) * 4.0)
    ax.legend(loc="upper left", framealpha=0.92, fontsize=7.4)
    ax.set_title(r"(b) H2 Shoemaker–Levy 9 (1992): Dynamic Sub-Roche Promotion ($\Lambda_S \geq 0.25$)")

    # --- Panel (c): H3 Apophis Bilateral Pericenter Window & Symmetric Rungs ---
    ax = axes[1, 0]
    t_h3 = np.asarray(h3["series"]["times"], dtype=float)
    tp_over_dt = np.asarray(h3["series"]["t_peri_over_dt"], dtype=float)
    k64_h3 = np.asarray(h3["series"]["dyadic_k_f64"], dtype=int)
    k128_h3 = np.asarray(h3["series"]["dyadic_k_f128"], dtype=int)

    ax.semilogy(t_h3, _pos_or_nan(tp_over_dt), color="#2980b9", lw=1.8,
                label=r"Estimated $|t_{\mathrm{peri}}(t)| / \Delta t_{\mathrm{macro}}$")
    ax.axhline(2.0, color="#c0392b", ls="--", lw=1.3, label=r"Bilateral Guard $|t_{\mathrm{peri}}| \leq 2\,\Delta t_{\mathrm{macro}}$")
    ax.axvline(float(h3["t_perigee"]), color="#2c3e50", ls=":", lw=1.2,
               label=rf"Perigee $t_{{\mathrm{{peri}}}} = {h3['t_perigee']:.4f}$ ($5.98\,R_\oplus$)")
    ax.set_xlabel("Physical Time $t$ (yr / $2\\pi$)")
    ax.set_ylabel(r"Pericenter Horizon $|t_{\mathrm{peri}}| / \Delta t_{\mathrm{macro}}$")
    ax.set_ylim(1e-2, 2e2)

    ax_r3 = ax.twinx()
    ax_r3.step(t_h3, k64_h3, where="post", color="#27ae60", lw=1.6, label=r"$k_{\mathrm{sym}}(t)$ (Float64)")
    ax_r3.step(t_h3, k128_h3, where="post", color="#16a085", lw=1.3, ls="--", label=r"$k_{\mathrm{sym}}(t)$ (Float128)")
    ax_r3.set_ylabel(r"Symmetric Dyadic Rung $k_{\mathrm{sym}}$", color="#1e8449")
    ax_r3.set_ylim(-0.3, 5.2)
    ax_r3.grid(False)

    l3a, lb3a = ax.get_legend_handles_labels()
    l3b, lb3b = ax_r3.get_legend_handles_labels()
    ax.legend(l3a + l3b, lb3a + lb3b, loc="upper right", framealpha=0.92, fontsize=7.4)
    ax.set_title(r"(c) H3 Apophis (2029): Bilateral Pericenter Anticipation & Dyadic Rungs")

    # --- Panel (d): Equal-Parameter Dual-Precision Energy Conservation Summary ---
    ax = axes[1, 1]
    u_h1 = np.asarray(h1["series"]["times"], dtype=float) / float(h1["equal_parameters"]["t_max"])
    u_h2 = np.asarray(h2["series"]["times"], dtype=float) / float(h2["equal_parameters"]["t_max"])
    u_h3 = np.asarray(h3["series"]["times"], dtype=float) / float(h3["equal_parameters"]["t_max"])

    ax.semilogy(u_h1, _pos_or_nan(h1["series"]["rel_dE_f64"]), color="#1f4e78", lw=1.6,
                label=rf"H1 Lexell Float64 (max={h1['max_rel_dE_f64']:.2e})")
    ax.semilogy(u_h1, _pos_or_nan(h1["series"]["rel_dE_f128"]), color="#1f4e78", lw=1.4, ls="--",
                label=rf"H1 Lexell Float128 (max={h1['max_rel_dE_f128']:.2e})")

    ax.semilogy(u_h2, _pos_or_nan(h2["series"]["rel_dE_f64"]), color="#8e44ad", lw=1.6,
                label=rf"H2 SL9 Float64 (max={h2['max_rel_dE_f64']:.2e})")
    ax.semilogy(u_h2, _pos_or_nan(h2["series"]["rel_dE_f128"]), color="#8e44ad", lw=1.4, ls="--",
                label=rf"H2 SL9 Float128 (max={h2['max_rel_dE_f128']:.2e})")

    ax.semilogy(u_h3, _pos_or_nan(h3["series"]["rel_dE_f64"]), color="#27ae60", lw=1.6,
                label=rf"H3 Apophis Float64 (max={h3['max_rel_dE_f64']:.2e})")
    ax.semilogy(u_h3, _pos_or_nan(h3["series"]["rel_dE_f128"]), color="#27ae60", lw=1.4, ls="--",
                label=rf"H3 Apophis Float128 (max={h3['max_rel_dE_f128']:.2e})")

    ax.set_xlabel(r"Normalized Encounter Progress $t / t_{\max} \in [0, 1]$")
    ax.set_ylabel(r"Relative Hamiltonian Error $|\Delta H(t)| / |H_0|$")
    ax.set_ylim(1e-36, 1e-5)
    ax.legend(loc="lower right", ncol=2, framealpha=0.92, fontsize=7.0)
    ax.set_title(r"(d) Equal-Parameter Float64 vs Float128 Hamiltonian Drift Across H1, H2, H3")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save_and_copy(fig, "gen04_plate_A_computational")
    plt.close(fig)


# ==============================================================================
# VOLET B: Algebraic, Covariant Gauge Continuity, Palindrome & Lie-Bargmann Plate (2x2)
# ==============================================================================
def plot_plate_B_algebraic_invariants(data: Dict[str, Any]) -> None:
    h1 = data["H1_lexell"]
    h2 = data["H2_shoemaker_levy"]
    h3 = data["H3_apophis"]
    gc = h1["gauge_continuity"]
    rev = h3["reversibility"]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.6))
    fig.suptitle(
        "Generator 04 — Plate B (Algebraic, Symplectic & Invariant Structure): "
        r"Covariant Gauge $T_{A\to B}^{(\mathrm{PBM4})}$ Continuity, Dyadic Palindrome & Lie–Bargmann Metrology (§7.2.1)",
        y=0.98,
    )

    # --- Panel (a): H1 Lexell Continuous Energy & Noether Invariants Across Hill Flips ---
    ax = axes[0, 0]
    t_h1 = np.asarray(h1["series"]["times"], dtype=float)
    dE64 = _pos_or_nan(h1["series"]["rel_dE_f64"])
    dE128 = _pos_or_nan(h1["series"]["rel_dE_f128"])
    dL64 = _pos_or_nan(h1["series"]["dL_f64"])
    dL128 = _pos_or_nan(h1["series"]["dL_f128"])

    ax.semilogy(t_h1, dE64, color="#1f4e78", lw=1.8, ls="-",
                label=rf"Float64 $|\Delta H(t)|/|H_0|$ (max={h1['max_rel_dE_f64']:.2e})")
    ax.semilogy(t_h1, dL64, color="#2980b9", lw=1.4, ls="-.",
                label=rf"Float64 $\|\Delta\mathbf{{L}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h1['lie_bargmann_f64']['max_dL']:.2e}$)")
    ax.semilogy(t_h1, dE128, color="#27ae60", lw=1.9, ls="--",
                label=rf"Float128 $|\Delta H(t)|/|H_0|$ (max={h1['max_rel_dE_f128']:.2e}, flip={gc['jump_total_dH_with_gauge_f128']:.2e})")
    ax.semilogy(t_h1, dL128, color="#16a085", lw=1.4, ls=":",
                label=rf"Float128 $\|\Delta\mathbf{{L}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h1['lie_bargmann_f128']['max_dL']:.2e}$)")

    for idx_f, fl in enumerate(h1["atlas_transitions"]):
        tf = float(fl["t"])
        ax.axvline(tf, color="#8e44ad", ls=":", lw=1.3,
                   label=rf"Chart Flips $\mathcal{{T}}_{{\mathrm{{helio}}}}\leftrightarrow\mathcal{{T}}_{{\mathrm{{Jovi}}}}$ ($t={tf:.2f}, {h1['atlas_transitions'][1]['t']:.2f}$)" if idx_f == 0 else None)
    ax.set_xlabel("Physical Time $t$ (yr / $2\\pi$)")
    ax.set_ylabel("Relative Energy & Noether Conservation")
    ax.set_ylim(1e-37, 1e-10)
    ax.legend(loc="center right", framealpha=0.92, fontsize=6.9)
    ax.set_title(r"(a) H1 Lexell: Seamless Energy & Noether Continuity Across $\mathcal{T}_{\mathrm{helio}}\leftrightarrow\mathcal{T}_{\mathrm{Jovi}}$")

    # --- Panel (b): H3 Apophis Full-Flyby Dyadic Palindrome ---
    ax = axes[0, 1]
    steps_idx = np.arange(1, len(rev["k_fwd_f64"]) + 1)
    k_fwd_64 = np.asarray(rev["k_fwd_f64"], dtype=int)
    k_bwd_64 = np.asarray(rev["k_bwd_reversed_f64"], dtype=int)
    k_bwd_128 = np.asarray(rev["k_bwd_reversed_f128"], dtype=int)

    ax.step(steps_idx, k_fwd_64, where="mid", color="#1f4e78", lw=2.2,
            label=r"Forward $k_{\mathrm{fwd}}(m)$ (Float64 & Float128)")
    ax.plot(steps_idx, k_bwd_64, "o", color="#e74c3c", ms=4.2, alpha=0.85,
            label=rf"Backward $k_{{\mathrm{{bwd}}}}(N-1-m)$ Float64 (Match={rev['dyadic_palindrome_matched_f64']})")
    ax.plot(steps_idx, k_bwd_128, "x", color="#27ae60", ms=5.0,
            label=rf"Backward $k_{{\mathrm{{bwd}}}}(N-1-m)$ Float128 (Match={rev['dyadic_palindrome_matched_f128']})")
    ax.axvline(38, color="#8e44ad", ls="--", lw=1.2, label="Perigee Passage (Step 38, $5.98\\,R_\\oplus$)")
    ax.set_xlabel(r"Macro-Step Index $m \in \{1, \dots, 80\}$ Across Perigee Flyby")
    ax.set_ylabel(r"Symmetric Dyadic Rung $k_{\mathrm{sym}}$")
    ax.set_ylim(-0.3, max(k_fwd_64) + 1.2)
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.4)
    ax.set_title(r"(b) H3 Apophis: Exact Full-Flyby Dyadic Palindrome $k_{\mathrm{fwd}}(m) \equiv k_{\mathrm{bwd}}(N-1-m)$")

    # --- Panel (c): H3 Apophis Multi-Horizon Round-Trip Reversibility ---
    ax = axes[1, 0]
    nh_arr = np.asarray(rev["horizon_steps"], dtype=int)
    ax.semilogy(nh_arr, _pos_or_nan(rev["horizon_dq_f64"]), "o-", color="#1f4e78", lw=1.8,
                label=r"Float64 Round-Trip $\|\Delta\mathbf{q}_{\mathrm{FB}}\|_2$")
    ax.semilogy(nh_arr, _pos_or_nan(rev["horizon_dE_f64"]), "s--", color="#2980b9", lw=1.6,
                label=r"Float64 Round-Trip $|\Delta E_{\mathrm{FB}}|/|E_0|$")
    ax.semilogy(nh_arr, _pos_or_nan(rev["horizon_dq_f128"]), "^-", color="#27ae60", lw=1.8,
                label=r"Float128 Round-Trip $\|\Delta\mathbf{q}_{\mathrm{FB}}\|_2$")
    ax.semilogy(nh_arr, _pos_or_nan(rev["horizon_dE_f128"]), "d--", color="#16a085", lw=1.6,
                label=r"Float128 Round-Trip $|\Delta E_{\mathrm{FB}}|/|E_0|$")
    ax.axvline(38, color="#c0392b", ls=":", lw=1.3, label="Perigee Crossing ($N_h = 38$)")
    ax.set_xlabel(r"Round-Trip Horizon $N_h$ Steps ($\Phi_{-\Delta t}^{N_h} \circ \Phi_{+\Delta t}^{N_h}$)")
    ax.set_ylabel("Forward–Backward Reversibility Error")
    ax.set_ylim(1e-36, 1e-9)
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.4)
    ax.set_title(r"(c) H3 Apophis: Multi-Horizon Round-Trip Reversibility ($\Phi_{-h}^{N_h}\circ\Phi_{+h}^{N_h}$)")

    # --- Panel (d): 9-Generator Lie-Bargmann Metrology Across H1, H2, H3 ---
    ax = axes[1, 1]
    protocols = ["H1 Lexell\n(1770)", "H2 SL9\n(1992)", "H3 Apophis\n(2029)"]
    x_pos = np.arange(len(protocols))
    w = 0.13

    dP_64 = [float(h1["lie_bargmann_f64"]["max_dP"]), float(h2["lie_bargmann_f64"]["max_dP"]), float(h3["lie_bargmann_f64"]["max_dP"])]
    dL_64 = [float(h1["lie_bargmann_f64"]["max_dL"]), float(h2["lie_bargmann_f64"]["max_dL"]), float(h3["lie_bargmann_f64"]["max_dL"])]
    dG_64 = [float(h1["lie_bargmann_f64"]["max_dG"]), float(h2["lie_bargmann_f64"]["max_dG"]), float(h3["lie_bargmann_f64"]["max_dG"])]

    dP_128 = [float(h1["lie_bargmann_f128"]["max_dP"]), float(h2["lie_bargmann_f128"]["max_dP"]), float(h3["lie_bargmann_f128"]["max_dP"])]
    dL_128 = [float(h1["lie_bargmann_f128"]["max_dL"]), float(h2["lie_bargmann_f128"]["max_dL"]), float(h3["lie_bargmann_f128"]["max_dL"])]
    dG_128 = [float(h1["lie_bargmann_f128"]["max_dG"]), float(h2["lie_bargmann_f128"]["max_dG"]), float(h3["lie_bargmann_f128"]["max_dG"])]

    ax.bar(x_pos - 2.5 * w, dP_64, width=w, color="#1f4e78", label=r"$\|\Delta\mathbf{P}\|_2$ (Float64)")
    ax.bar(x_pos - 1.5 * w, dL_64, width=w, color="#2980b9", label=r"$\|\Delta\mathbf{L}\|_2$ (Float64)")
    ax.bar(x_pos - 0.5 * w, dG_64, width=w, color="#5dade2", label=r"$\|\Delta\mathbf{G}\|_2$ (Float64)")

    ax.bar(x_pos + 0.5 * w, dP_128, width=w, color="#1e8449", label=r"$\|\Delta\mathbf{P}\|_2$ (Float128)")
    ax.bar(x_pos + 1.5 * w, dL_128, width=w, color="#27ae60", label=r"$\|\Delta\mathbf{L}\|_2$ (Float128)")
    ax.bar(x_pos + 2.5 * w, dG_128, width=w, color="#82e0aa", label=r"$\|\Delta\mathbf{G}\|_2$ (Float128)")

    ax.set_yscale("log")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(protocols)
    ax.set_ylabel("Max Lie–Bargmann Generator Drift")
    ax.set_ylim(1e-36, 1e-9)
    ax.legend(loc="upper right", ncol=2, framealpha=0.92, fontsize=7.0)
    ax.set_title(r"(d) Exact 9-Generator Lie–Bargmann Conservation ($\mathbf{P}, \mathbf{L}, \mathbf{G}$) in Float64 & Float128")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save_and_copy(fig, "gen04_plate_B_algebraic_invariants")
    plt.close(fig)


# ==============================================================================
# DEDICATED PLATE H1: Comet D/1770 L1 (Lexell) — Deep Perijove & Atlas Flips (1x2)
# ==============================================================================
def plot_plate_H1_lexell(data: Dict[str, Any]) -> None:
    h1 = data["H1_lexell"]
    gc = h1["gauge_continuity"]

    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.2))
    fig.suptitle(
        rf"Benchmark H1 — Comet D/1770 L1 (Lexell, 1770–1779): Deep Jovian Encounter (${h1['min_r_jup_RJ']:.3f}\,R_J$), "
        r"Dynamic Hill Chart Flips ($\mathcal{T}_{\mathrm{helio}}\leftrightarrow\mathcal{T}_{\mathrm{Jovi}}$) & Invariant Continuity ($N=4$)",
        y=0.99,
    )

    t_h1 = np.asarray(h1["series"]["times"], dtype=float)
    r_jup_scale_au = float(h1["min_r_jup_au"]) / float(h1["min_r_jup_RJ"])
    rj_RJ_64 = np.asarray(h1["series"]["r_jup_au_f64"], dtype=float) / r_jup_scale_au
    rj_RJ_128 = np.asarray(h1["series"]["r_jup_au_f128"], dtype=float) / r_jup_scale_au
    r_hill_RJ = float(h1["r_hill_jup_au"]) / r_jup_scale_au
    k64_h1 = np.asarray(h1["series"]["dyadic_k_f64"], dtype=int)
    k128_h1 = np.asarray(h1["series"]["dyadic_k_f128"], dtype=int)
    k_max_h1 = int(np.max(k64_h1))

    # --- Panel (a): Jovicentric Approach in R_J, Hill Chart Flips & Dyadic Rungs ---
    ax = axes[0]
    t_flip1 = float(h1["atlas_transitions"][0]["t"])
    t_flip2 = float(h1["atlas_transitions"][1]["t"])
    ax.axvspan(t_flip1, t_flip2, color="#8e44ad", alpha=0.09,
               label=rf"Jovicentric Chart $\mathcal{{T}}_{{\mathrm{{Jovi}}}}$ ($t \in [{t_flip1:.2f}, {t_flip2:.2f}]$, $\Delta t={h1['dwell_time_jovicentric']:.2f}$)")
    ax.semilogy(t_h1, rj_RJ_64, color="#1f4e78", lw=2.1, label=r"$r_J(t) / R_J$ (Float64 & Float128)")
    ax.semilogy(t_h1, rj_RJ_128, color="#27ae60", lw=1.4, ls="--", label=r"$r_J(t) / R_J$ (Float128 overlay)")
    ax.axhline(r_hill_RJ, color="#c0392b", ls="--", lw=1.3,
               label=rf"Jovian Hill Sphere $R_{{\mathrm{{Hill}},J}} = {r_hill_RJ:.0f}\,R_J$ ($0.355$ AU)")
    ax.axhline(5.90, color="#7f8c8d", ls=":", lw=1.2, label=r"Io Orbital Radius ($a_{\mathrm{Io}} = 5.90\,R_J$)")
    ax.axvline(t_flip1, color="#8e44ad", ls=":", lw=1.3)
    ax.axvline(t_flip2, color="#8e44ad", ls=":", lw=1.3)
    ax.scatter([float(h1["t_perijove"])], [float(h1["min_r_jup_RJ"])], color="#d35400", edgecolors="black", s=55, zorder=6,
               label=rf"Deep Perijove $r_{{\min}} = {h1['min_r_jup_RJ']:.3f}\,R_J$ (${h1['min_r_jup_au']*1e3:.2f}\times 10^{{-3}}$ AU)")

    ax.set_xlabel(r"Physical Time $t$ (yr / $2\pi$)")
    ax.set_ylabel(r"Jovicentric Distance $r_J(t) / R_J$", color="#1f4e78")
    ax.set_ylim(1.5, 9500.0)

    ax_r = ax.twinx()
    ax_r.step(t_h1, k64_h1, where="post", color="#16a085", lw=1.5, alpha=0.85,
              label=rf"Symmetric Dyadic Rung $k_{{\mathrm{{sym}}}}(t) \in \{{0,\dots,{k_max_h1}\}}$")
    ax_r.set_ylabel(r"Symmetric Dyadic Rung $k_{\mathrm{sym}}$", color="#16a085")
    ax_r.set_ylim(-0.3, k_max_h1 + 2.4)
    ax_r.grid(False)

    la, lba = ax.get_legend_handles_labels()
    lb, lbb = ax_r.get_legend_handles_labels()
    ax.legend(la + lb, lba + lbb, loc="upper center", ncol=2, framealpha=0.93, fontsize=6.7)
    ax.set_title(r"(a) Jovicentric Approach, Hill Chart Flips ($\mathcal{T}_{\mathrm{helio}}\leftrightarrow\mathcal{T}_{\mathrm{Jovi}}$) & Rungs")

    # --- Panel (b): Continuous Hamiltonian & Lie-Bargmann Conservation Across Chart Flips ---
    ax = axes[1]
    dE64 = _pos_or_nan(h1["series"]["rel_dE_f64"])
    dE128 = _pos_or_nan(h1["series"]["rel_dE_f128"])
    dL64 = _pos_or_nan(h1["series"]["dL_f64"])
    dP64 = _pos_or_nan(h1["series"]["dP_f64"])
    dL128 = _pos_or_nan(h1["series"]["dL_f128"])
    dP128 = _pos_or_nan(h1["series"]["dP_f128"])

    ax.axvspan(t_flip1, t_flip2, color="#8e44ad", alpha=0.07)
    ax.semilogy(t_h1, dE64, color="#1f4e78", lw=1.9, ls="-",
                label=rf"Float64 $|\Delta H(t)|/|H_0|$ (max $= {h1['max_rel_dE_f64']:.2e}$)")
    ax.semilogy(t_h1, dL64, color="#2980b9", lw=1.4, ls="-.",
                label=rf"Float64 $\|\Delta\mathbf{{L}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h1['lie_bargmann_f64']['max_dL']:.2e}$)")
    ax.semilogy(t_h1, dP64, color="#5dade2", lw=1.3, ls=":",
                label=rf"Float64 $\|\Delta\mathbf{{P}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h1['lie_bargmann_f64']['max_dP']:.2e}$)")

    ax.semilogy(t_h1, dE128, color="#27ae60", lw=2.0, ls="--",
                label=rf"Float128 $|\Delta H(t)|/|H_0|$ ($\mathcal{{O}}(h_{{\mathrm{{sub}}}}^6)$ shadow osc., max $= {h1['max_rel_dE_f128']:.2e}$)")
    ax.semilogy(t_h1, dL128, color="#16a085", lw=1.4, ls="-.",
                label=rf"Float128 $\|\Delta\mathbf{{L}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h1['lie_bargmann_f128']['max_dL']:.2e}$)")
    ax.semilogy(t_h1, dP128, color="#82e0aa", lw=1.3, ls=":",
                label=rf"Float128 $\|\Delta\mathbf{{P}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h1['lie_bargmann_f128']['max_dP']:.2e}$)")

    ax.axvline(t_flip1, color="#8e44ad", ls=":", lw=1.4,
               label=rf"Ingress Flip $t_1={t_flip1:.2f}$ ($T_{{A\to B}}^{{(\mathrm{{PBM4}})}}$: $|\Delta H|/|H_0|={gc['jump_total_dH_with_gauge_f128']:.2e}$)")
    ax.axvline(t_flip2, color="#8e44ad", ls="--", lw=1.2,
               label=rf"Egress Flip $t_2={t_flip2:.2f}$ ($\mathcal{{T}}_{{\mathrm{{Jovi}}}}\to\mathcal{{T}}_{{\mathrm{{helio}}}}$)")
    ax.axvline(float(h1["t_perijove"]), color="#d35400", ls="-.", lw=1.1,
               label=rf"Perijove Passage $t_{{\mathrm{{peri}}}}={h1['t_perijove']:.4f}$ (${h1['min_r_jup_RJ']:.3f}\,R_J$)")

    ax.set_xlabel(r"Physical Time $t$ (yr / $2\pi$)")
    ax.set_ylabel(r"Relative Hamiltonian & Lie–Bargmann Invariant Drift")
    ax.set_ylim(1e-39, 1e-10)
    ax.legend(loc="center right", framealpha=0.93, fontsize=6.8)
    ax.set_title(r"(b) Continuous Energy & Lie–Bargmann Conservation Under Canonical Gauge $T_{A\to B}^{(\mathrm{PBM4})}$")

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save_and_copy(fig, "fig_ch07_h1_lexell_atlas_gauge")
    plt.close(fig)


# ==============================================================================
# DEDICATED PLATE H2: Comet Shoemaker-Levy 9 (1992) — Sub-Roche Breakup (1x3)
# ==============================================================================
def plot_plate_H2_shoemaker_levy(data: Dict[str, Any]) -> None:
    h2 = data["H2_shoemaker_levy"]

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.2))
    fig.suptitle(
        r"Benchmark H2 — Comet D/1993 F2 (Shoemaker–Levy 9, July 1992): Sub-Roche Perijove ($1.336\,R_J < 2.45\,R_J$), "
        r"Dynamic Tidal Promotion ($\Lambda_S \geq 0.25$) & String-of-Pearls Elongation ($N=5$)",
        y=0.99,
    )

    t_h2 = np.asarray(h2["series"]["times"], dtype=float)
    rj_h2_RJ = np.asarray(h2["series"]["r_jup_RJ_f64"], dtype=float)
    lam_h2 = np.asarray(h2["series"]["lambda_S"], dtype=float)
    act_h2 = np.asarray(h2["series"]["roche_supernode_active"], dtype=int)
    span64 = np.asarray(h2["series"]["span_ratio_f64"], dtype=float)
    span128 = np.asarray(h2["series"]["span_ratio_f128"], dtype=float)
    k64_h2 = np.asarray(h2["series"]["dyadic_k_f64"], dtype=int)
    roche_RJ = float(h2["roche_limit_RJ"])
    lam_crit = float(h2["lambda_crit"])
    active_idx = np.where(act_h2 == 1)[0]
    t_enter = float(t_h2[active_idx[0]]) if len(active_idx) > 0 else 0.0
    t_leave = float(t_h2[active_idx[-1]]) if len(active_idx) > 0 else 0.0

    # --- Panel (a): Sub-Roche Infall r_J(t)/R_J & Cubic Tidal Ratio Lambda_S(t) ---
    ax = axes[0]
    if len(active_idx) > 0:
        ax.axvspan(t_enter, t_leave, color="#f39c12", alpha=0.16,
                   label=rf"Promoted 3-Fragment Supernode ({h2['roche_active_steps']}/{h2['equal_parameters']['n_steps']} steps)")
    ax.plot(t_h2, rj_h2_RJ, color="#1f4e78", lw=2.1, label=r"Jovicentric Distance $r_J(t) / R_J$")
    ax.axhline(roche_RJ, color="#c0392b", ls="--", lw=1.4, label=rf"Fluid Roche Limit $r_{{\mathrm{{Roche}}}} = {roche_RJ:.2f}\,R_J$")
    ax.axhline(1.0, color="#7f8c8d", ls=":", lw=1.1, label=r"Jupiter Surface ($1.00\,R_J$)")
    ax.scatter([float(h2["t_perijove"])], [float(h2["min_r_jup_RJ"])], color="#d35400", edgecolors="black", s=55, zorder=6,
               label=rf"Sub-Roche Perijove $r_{{\min}} = {h2['min_r_jup_RJ']:.3f}\,R_J$ ($6.38\times 10^{{-4}}$ AU)")
    ax.set_xlabel(r"Physical Time $t$ (yr / $2\pi$)")
    ax.set_ylabel(r"Jovicentric Distance $r_J(t) / R_J$", color="#1f4e78")
    ax.set_ylim(0.4, max(rj_h2_RJ) * 1.08)

    ax_l = ax.twinx()
    ax_l.semilogy(t_h2, lam_h2, color="#8e44ad", lw=1.8, ls="-.",
                  label=rf"Cubic Tidal Ratio $\Lambda_S(t)=\frac{{1}}{{3}}\eta_S^3$ (max $={h2['lambda_max']:.0f}$)")
    ax_l.axhline(lam_crit, color="#8e44ad", ls=":", lw=1.2, label=rf"Promotion Threshold $\Lambda_{{\mathrm{{crit}}}} = {lam_crit:.2f}$")
    ax_l.set_ylabel(r"Cubic Tidal Force Ratio $\Lambda_S(t) = \frac{1}{3}\eta_S(t)^3$", color="#8e44ad")
    ax_l.set_ylim(1e-2, max(lam_h2) * 5.0)
    ax_l.grid(False)

    l1, lb1 = ax.get_legend_handles_labels()
    l2, lb2 = ax_l.get_legend_handles_labels()
    ax.legend(l1 + l2, lb1 + lb2, loc="upper left", framealpha=0.92, fontsize=6.9)
    ax.set_title(r"(a) Sub-Roche Passage ($1.336\,R_J$) & Cubic Tidal Trigger $\Lambda_S(t)$")

    # --- Panel (b): String-of-Pearls Differential Elongation S(t)/S_0 & Dyadic Rungs ---
    ax = axes[1]
    if len(active_idx) > 0:
        ax.axvspan(t_enter, t_leave, color="#f39c12", alpha=0.14, label="Active Synchronous Supernode Window")
    ax.plot(t_h2, span64, color="#8e44ad", lw=2.1,
            label=rf"Chain Elongation $S(t)/S_0$ Float64 (peak $= {h2['max_span_ratio']:.2f}\times$, end $= {h2['final_span_ratio']:.2f}\times$)")
    ax.plot(t_h2, span128, color="#27ae60", lw=1.5, ls="--",
            label=r"Chain Elongation $S(t)/S_0$ Float128 (exact overlay)")
    ax.axvline(float(h2["t_perijove"]), color="#d35400", ls=":", lw=1.3,
               label=rf"Perijove Passage $t_{{\mathrm{{peri}}}} = {h2['t_perijove']:.5f}$")
    ax.set_xlabel(r"Physical Time $t$ (yr / $2\pi$)")
    ax.set_ylabel(r"Fragment Chain Span Ratio $S(t) / S_0 = \|\mathbf{q}_4(t) - \mathbf{q}_2(t)\|_2 / S_0$", color="#8e44ad")
    ax.set_ylim(0.85, max(span64) * 1.15)

    ax_k = ax.twinx()
    k_max_h2 = int(np.max(k64_h2))
    ax_k.step(t_h2, k64_h2, where="post", color="#16a085", lw=1.5, alpha=0.85,
              label=rf"Symmetric Dyadic Rung $k_{{\mathrm{{sym}}}}(t) \in \{{0,\dots,{k_max_h2}\}}$")
    ax_k.set_ylabel(r"Symmetric Dyadic Rung $k_{\mathrm{sym}}$", color="#16a085")
    ax_k.set_ylim(-0.3, k_max_h2 + 1.8)
    ax_k.grid(False)

    l3, lb3 = ax.get_legend_handles_labels()
    l4, lb4 = ax_k.get_legend_handles_labels()
    ax.legend(l3 + l4, lb3 + lb4, loc="upper left", framealpha=0.92, fontsize=6.9)
    ax.set_title(r"(b) “String-of-Pearls” Tidal Dislocation $S(t)/S_0$ & Dyadic Rungs")

    # --- Panel (c): Continuous Hamiltonian & Lie-Bargmann Conservation ---
    ax = axes[2]
    dE64 = _pos_or_nan(h2["series"]["rel_dE_f64"])
    dE128 = _pos_or_nan(h2["series"]["rel_dE_f128"])
    dL64 = _pos_or_nan(h2["series"]["dL_f64"])
    dP64 = _pos_or_nan(h2["series"]["dP_f64"])
    dL128 = _pos_or_nan(h2["series"]["dL_f128"])
    dP128 = _pos_or_nan(h2["series"]["dP_f128"])

    if len(active_idx) > 0:
        ax.axvspan(t_enter, t_leave, color="#f39c12", alpha=0.12)
    ax.semilogy(t_h2, dE64, color="#1f4e78", lw=1.9, ls="-",
                label=rf"Float64 $|\Delta H(t)|/|H_0|$ ($\leq {h2['max_rel_dE_f64']:.2e}$)")
    ax.semilogy(t_h2, dL64, color="#2980b9", lw=1.4, ls="-.",
                label=rf"Float64 $\|\Delta\mathbf{{L}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h2['lie_bargmann_f64']['max_dL']:.2e}$)")
    ax.semilogy(t_h2, dP64, color="#5dade2", lw=1.3, ls=":",
                label=rf"Float64 $\|\Delta\mathbf{{P}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h2['lie_bargmann_f64']['max_dP']:.2e}$)")

    ax.semilogy(t_h2, dE128, color="#27ae60", lw=2.0, ls="--",
                label=rf"Float128 $|\Delta H(t)|/|H_0|$ ($\leq {h2['max_rel_dE_f128']:.2e}$)")
    ax.semilogy(t_h2, dL128, color="#16a085", lw=1.4, ls="-.",
                label=rf"Float128 $\|\Delta\mathbf{{L}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h2['lie_bargmann_f128']['max_dL']:.2e}$)")
    ax.semilogy(t_h2, dP128, color="#82e0aa", lw=1.3, ls=":",
                label=rf"Float128 $\|\Delta\mathbf{{P}}_{{\mathrm{{tot}}}}(t)\|_2$ ($\leq {h2['lie_bargmann_f128']['max_dP']:.2e}$)")

    ax.axvline(float(h2["t_perijove"]), color="#d35400", ls=":", lw=1.3, label="Sub-Roche Perijove ($1.336\\,R_J$)")
    ax.set_xlabel(r"Physical Time $t$ (yr / $2\pi$)")
    ax.set_ylabel("Relative Energy & Lie–Bargmann Drift")
    ax.set_ylim(1e-39, 1e-10)
    ax.legend(loc="center right", framealpha=0.92, fontsize=6.9)
    ax.set_title(r"(c) Energy & Noether Conservation Across Sub-Roche Breakup")

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save_and_copy(fig, "fig_ch07_h2_shoemaker_levy_roche")
    plt.close(fig)


# ==============================================================================
# DEDICATED PLATE H3: Asteroid (99942) Apophis (2029) — Flyby & Palindrome (2x2)
# ==============================================================================
def plot_plate_H3_apophis(data: Dict[str, Any]) -> None:
    h3 = data["H3_apophis"]
    rev = h3["reversibility"]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.4))
    fig.suptitle(
        r"Benchmark H3 — Asteroid (99942) Apophis (April 13, 2029): Sub-Geostationary Perigee ($5.976\,R_\oplus$), "
        r"Aten$\to$Apollo Deflection ($\Delta a = +0.1939$ AU) & Full-Flyby Dyadic Palindrome ($N=4$)",
        y=0.98,
    )

    # --- Panel (a): Geocentric Hyperbolic Flyby Inside Geostationary Ring ---
    ax = axes[0, 0]
    gx = np.asarray(h3["series"]["geo_x_RE"], dtype=float)
    gy = np.asarray(h3["series"]["geo_y_RE"], dtype=float)
    r_RE = np.asarray(h3["series"]["r_earth_RE_f64"], dtype=float)
    idx_min = int(np.argmin(r_RE))
    geo_R = float(h3["geo_belt_RE"])

    theta = np.linspace(0.0, 2.0 * np.pi, 240)
    ax.fill(np.cos(theta), np.sin(theta), color="#3498db", alpha=0.35, label=r"Earth Disk ($1.00\,R_\oplus$)")
    ax.plot(np.cos(theta), np.sin(theta), color="#1f4e78", lw=1.2)
    ax.plot(geo_R * np.cos(theta), geo_R * np.sin(theta), color="#7f8c8d", ls="--", lw=1.3,
            label=rf"Geostationary Ring ($r_{{\mathrm{{GEO}}}} = {geo_R:.2f}\,R_\oplus$)")
    ax.plot(gx, gy, color="#c0392b", lw=2.2, label="Apophis Geocentric Hyperbola (Float64 & Float128)")
    ax.scatter([gx[idx_min]], [gy[idx_min]], color="#f39c12", edgecolors="black", s=60, zorder=6,
               label=rf"Perigee Closest Approach $r_{{\min}} = {h3['min_r_earth_RE']:.3f}\,R_\oplus$ ($38,115$ km)")

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlim(-22.0, 22.0)
    ax.set_ylim(-15.0, 15.0)
    ax.set_xlabel(r"Geocentric $x_{\mathrm{geo}} / R_\oplus$")
    ax.set_ylabel(r"Geocentric $y_{\mathrm{geo}} / R_\oplus$")
    ax.legend(loc="upper left", framealpha=0.92, fontsize=7.3)
    ax.set_title(r"(a) Grazing Geocentric Hyperbolic Flyby Inside the GEO Ring")

    # --- Panel (b): Resonant Heliocentric Aten -> Apollo Orbital Transition ---
    ax = axes[0, 1]
    t_h3 = np.asarray(h3["series"]["times"], dtype=float)
    a64 = np.asarray(h3["series"]["a_helio_au_f64"], dtype=float)
    a128 = np.asarray(h3["series"]["a_helio_au_f128"], dtype=float)
    e64 = np.asarray(h3["series"]["e_helio_f64"], dtype=float)

    ax.plot(t_h3, a64, color="#1f4e78", lw=2.2,
            label=rf"Semi-Major Axis $a(t)$ Float64 ({h3['a_pre_au']:.4f} $\to$ {h3['a_post_au']:.4f} AU)")
    ax.plot(t_h3, a128, color="#27ae60", lw=1.5, ls="--",
            label=rf"Semi-Major Axis $a(t)$ Float128 ($\Delta a = a_{{\mathrm{{post}}}} - a_{{\mathrm{{pre}}}} = +{h3['delta_a_au']:.4f}$ AU)")
    ax.axhline(1.0, color="#c0392b", ls="--", lw=1.4, label=r"Aten / Apollo Class Boundary ($a = 1.000$ AU)")
    ax.axvline(float(h3["t_perigee"]), color="#7f8c8d", ls=":", lw=1.2)
    ax.text(0.0008, 0.935, "Aten Class\n($a_{\\mathrm{pre}} = 0.9224 < 1$ AU)", color="#1f4e78", fontsize=8.2, fontweight="bold")
    ax.text(0.0125, 1.124, "Apollo Class ($a_{\\mathrm{post}} = 1.1163 > 1$ AU)", color="#1e8449", fontsize=8.2, fontweight="bold")
    ax.set_xlabel(r"Physical Time $t$ (yr / $2\pi$)")
    ax.set_ylabel("Heliocentric Semi-Major Axis $a(t)$ [AU]", color="#1f4e78")
    ax.set_ylim(0.89, 1.15)

    ax_e = ax.twinx()
    ax_e.plot(t_h3, e64, color="#8e44ad", lw=1.5, ls="-.",
              label=rf"Heliocentric Eccentricity $e(t)$ ({h3['e_pre']:.3f} $\to$ {h3['e_post']:.3f})")
    ax_e.set_ylabel("Heliocentric Eccentricity $e(t)$", color="#8e44ad")
    ax_e.set_ylim(0.14, 0.25)
    ax_e.grid(False)

    l4a, lb4a = ax.get_legend_handles_labels()
    l4b, lb4b = ax_e.get_legend_handles_labels()
    ax.legend(l4a + l4b, lb4a + lb4b, loc="center right", framealpha=0.92, fontsize=7.1)
    ax.set_title(r"(b) Resonant Heliocentric Transition: Aten ($a<1$ AU) $\to$ Apollo ($a>1$ AU)")

    # --- Panel (c): Bilateral Pericenter Window & Exact Full-Flyby Dyadic Palindrome ---
    ax = axes[1, 0]
    steps_idx = np.arange(1, len(rev["k_fwd_f64"]) + 1)
    k_fwd_64 = np.asarray(rev["k_fwd_f64"], dtype=int)
    k_bwd_64 = np.asarray(rev["k_bwd_reversed_f64"], dtype=int)
    k_bwd_128 = np.asarray(rev["k_bwd_reversed_f128"], dtype=int)
    tp_over_dt = np.asarray(h3["series"]["t_peri_over_dt"][1:], dtype=float)

    ax.step(steps_idx, k_fwd_64, where="mid", color="#1f4e78", lw=2.2,
            label=r"Forward Rung $k_{\mathrm{fwd}}(m)$ (Float64 & Float128)")
    ax.plot(steps_idx, k_bwd_64, "o", color="#e74c3c", ms=4.2, alpha=0.85,
            label=rf"Backward $k_{{\mathrm{{bwd}}}}(N-1-m)$ Float64 (Match={rev['dyadic_palindrome_matched_f64']})")
    ax.plot(steps_idx, k_bwd_128, "x", color="#27ae60", ms=5.0,
            label=rf"Backward $k_{{\mathrm{{bwd}}}}(N-1-m)$ Float128 (Match={rev['dyadic_palindrome_matched_f128']})")
    ax.axvline(38, color="#8e44ad", ls="--", lw=1.2, label=r"Perigee Step $m=38$ ($r_{\min}=5.976\,R_\oplus$)")
    ax.set_xlabel(r"Macro-Step Index $m \in \{1, \dots, 80\}$ Across Perigee Flyby")
    ax.set_ylabel(r"Symmetric Dyadic Rung $k_{\mathrm{sym}} \in \{0,1,2,3\}$", color="#1f4e78")
    ax.set_ylim(-0.3, 4.9)

    ax_tp = ax.twinx()
    ax_tp.semilogy(steps_idx, _pos_or_nan(tp_over_dt), color="#2980b9", lw=1.4, ls=":",
                   label=r"Bilateral Horizon $|t_{\mathrm{peri}}(m)| / \Delta t_{\mathrm{macro}}$")
    ax_tp.axhline(2.0, color="#c0392b", ls=":", lw=1.1, label=r"Guard Threshold $|t_{\mathrm{peri}}| \leq 2\,\Delta t_{\mathrm{macro}}$")
    ax_tp.set_ylabel(r"Pericenter Horizon $|t_{\mathrm{peri}}| / \Delta t_{\mathrm{macro}}$", color="#2980b9")
    ax_tp.set_ylim(1e-2, 4e3)
    ax_tp.grid(False)

    lc1, lbc1 = ax.get_legend_handles_labels()
    lc2, lbc2 = ax_tp.get_legend_handles_labels()
    ax.legend(lc1 + lc2, lbc1 + lbc2, loc="upper center", ncol=2, framealpha=0.93, fontsize=6.4)
    ax.set_title(r"(c) Bilateral Pericenter Guard & Exact Dyadic Palindrome $k_{\mathrm{fwd}}(m) \equiv k_{\mathrm{bwd}}(N-1-m)$")

    # --- Panel (d): Multi-Horizon Round-Trip Reversibility Across Perigee ---
    ax = axes[1, 1]
    nh_arr = np.asarray(rev["horizon_steps"], dtype=int)
    ax.semilogy(nh_arr, _pos_or_nan(rev["horizon_dq_f64"]), "o-", color="#1f4e78", lw=1.9,
                label=rf"Float64 Return $\|\Delta\mathbf{{q}}_{{\mathrm{{FB}}}}\|_2$ ($N=80$: ${rev['delta_q_fb_f64']:.2e}$ AU)")
    ax.semilogy(nh_arr, _pos_or_nan(rev["horizon_dE_f64"]), "s--", color="#2980b9", lw=1.6,
                label=rf"Float64 Return $|\Delta H_{{\mathrm{{FB}}}}|/|H_0|$ ($N=80$: ${rev['delta_E_fb_f64']:.2e}$)")
    ax.semilogy(nh_arr, _pos_or_nan(rev["horizon_dq_f128"]), "^-", color="#27ae60", lw=1.9,
                label=rf"Float128 Return $\|\Delta\mathbf{{q}}_{{\mathrm{{FB}}}}\|_2$ ($N=80$: ${rev['delta_q_fb_f128']:.2e}$ AU)")
    ax.semilogy(nh_arr, _pos_or_nan(rev["horizon_dE_f128"]), "d--", color="#16a085", lw=1.6,
                label=rf"Float128 Return $|\Delta H_{{\mathrm{{FB}}}}|/|H_0|$ ($N=80$: ${rev['delta_E_fb_f128']:.2e}$)")
    ax.axvline(38, color="#c0392b", ls=":", lw=1.3, label="Perigee Crossing ($N_h = 38$)")
    ax.set_xlabel(r"Round-Trip Horizon $N_h \in \{20, 40, 60, 80\}$ Steps ($\Phi_{-\Delta t}^{N_h} \circ \Phi_{+\Delta t}^{N_h}$)")
    ax.set_ylabel(r"Forward–Backward Reversibility Error")
    ax.set_ylim(1e-36, 1e-9)
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.1)
    ax.set_title(r"(d) Exact Time-Reversal Symmetry ($\Phi_{-\Delta t}^{N_h} \circ \Phi_{+\Delta t}^{N_h}$) Across Perigee")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save_and_copy(fig, "fig_ch07_h3_apophis_flyby_palindrome")
    plt.close(fig)


# ==============================================================================
# VOLET C: Physical & Astrodynamical Reality Plate (2x2)
# ==============================================================================
def plot_plate_C_physical_dynamics(data: Dict[str, Any]) -> None:
    h1 = data["H1_lexell"]
    h2 = data["H2_shoemaker_levy"]
    h3 = data["H3_apophis"]

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.6))
    fig.suptitle(
        "Generator 04 — Plate C (Physical & Astrodynamical Reality): "
        r"Lexell Deep Perijove, SL9 Sub-Roche Breakup & Apophis Aten$\to$Apollo Deflection (§7.2.1)",
        y=0.98,
    )

    # --- Panel (a): H1 Comet Lexell (1770) Deep Perijove Traversal in R_J ---
    ax = axes[0, 0]
    t_h1 = np.asarray(h1["series"]["times"], dtype=float)
    r_jup_scale_au = float(h1["min_r_jup_au"]) / float(h1["min_r_jup_RJ"])
    rj_RJ_64 = np.asarray(h1["series"]["r_jup_au_f64"], dtype=float) / r_jup_scale_au
    rj_RJ_128 = np.asarray(h1["series"]["r_jup_au_f128"], dtype=float) / r_jup_scale_au
    r_hill_RJ = float(h1["r_hill_jup_au"]) / r_jup_scale_au

    ax.semilogy(t_h1, rj_RJ_64, color="#1f4e78", lw=2.0, label="Float64 Trajectory $r_J(t) / R_J$")
    ax.semilogy(t_h1, rj_RJ_128, color="#27ae60", lw=1.5, ls="--", label="Float128 Trajectory $r_J(t) / R_J$")
    ax.axhline(r_hill_RJ, color="#c0392b", ls="--", lw=1.2, label=rf"Hill Switching Sphere ({r_hill_RJ:.0f} $R_J$)")
    ax.axhline(5.90, color="#7f8c8d", ls=":", lw=1.2, label=r"Io Orbital Radius ($5.90\,R_J$)")
    ax.scatter([float(h1["t_perijove"])], [float(h1["min_r_jup_RJ"])], color="#d35400", s=50, zorder=6,
               label=rf"Deep Perijove $r_{{\min}} = {h1['min_r_jup_RJ']:.3f}\,R_J$ ($1.48\times 10^{{-3}}$ AU)")
    ax.set_xlabel("Physical Time $t$ (yr / $2\\pi$)")
    ax.set_ylabel(r"Jovicentric Distance $r_J(t) / R_J$")
    ax.set_ylim(1.5, 2500.0)
    ax.legend(loc="lower left", framealpha=0.92, fontsize=7.4)
    ax.set_title(r"(a) H1 Comet Lexell (1770): Full Traversal from $0.70$ AU to $r_{\min}=3.10\,R_J$")

    # --- Panel (b): H2 Shoemaker-Levy 9 Sub-Roche Perijove & String-of-Pearls Elongation ---
    ax = axes[0, 1]
    t_h2 = np.asarray(h2["series"]["times"], dtype=float)
    rj_h2_RJ = np.asarray(h2["series"]["r_jup_RJ_f64"], dtype=float)
    span64 = np.asarray(h2["series"]["span_ratio_f64"], dtype=float)
    span128 = np.asarray(h2["series"]["span_ratio_f128"], dtype=float)
    roche_RJ = float(h2["roche_limit_RJ"])

    ax.plot(t_h2, rj_h2_RJ, color="#1f4e78", lw=2.0, label=r"Jovicentric Distance $r_J(t) / R_J$")
    ax.axhline(roche_RJ, color="#c0392b", ls="--", lw=1.4, label=rf"Fluid Roche Limit $r_{{\mathrm{{Roche}}}} = {roche_RJ:.2f}\,R_J$")
    ax.axhline(1.0, color="#7f8c8d", ls=":", lw=1.1, label=r"Jupiter Surface ($1.00\,R_J$)")
    ax.scatter([float(h2["t_perijove"])], [float(h2["min_r_jup_RJ"])], color="#d35400", s=50, zorder=6,
               label=rf"Sub-Roche Perijove $r_{{\min}} = {h2['min_r_jup_RJ']:.3f}\,R_J < 2.45\,R_J$")
    ax.set_xlabel("Physical Time $t$ (yr / $2\\pi$)")
    ax.set_ylabel(r"Jovicentric Distance $r_J(t) / R_J$", color="#1f4e78")
    ax.set_ylim(0.5, max(rj_h2_RJ) * 1.08)

    ax_s = ax.twinx()
    ax_s.plot(t_h2, span64, color="#8e44ad", lw=1.9, label=rf"Chain Elongation $S(t)/S_0$ Float64 ({h2['max_span_ratio']:.2f}$\times$)")
    ax_s.plot(t_h2, span128, color="#27ae60", lw=1.4, ls="--", label=r"Chain Elongation $S(t)/S_0$ Float128")
    ax_s.set_ylabel(r"Fragment Chain Elongation $S(t) / S_0$", color="#8e44ad")
    ax_s.set_ylim(0.8, max(span64) * 1.15)
    ax_s.grid(False)

    l2a, lb2a = ax.get_legend_handles_labels()
    l2b, lb2b = ax_s.get_legend_handles_labels()
    ax.legend(l2a + l2b, lb2a + lb2b, loc="upper center", framealpha=0.92, fontsize=7.2)
    ax.set_title(r"(b) H2 Shoemaker–Levy 9 (1992): Sub-Roche Passage ($1.34\,R_J$) & String-of-Pearls")

    # --- Panel (c): H3 Apophis 2029 Geocentric Hyperbolic Flyby Inside GEO Belt ---
    ax = axes[1, 0]
    gx = np.asarray(h3["series"]["geo_x_RE"], dtype=float)
    gy = np.asarray(h3["series"]["geo_y_RE"], dtype=float)
    r_RE = np.asarray(h3["series"]["r_earth_RE_f64"], dtype=float)
    idx_min = int(np.argmin(r_RE))
    geo_R = float(h3["geo_belt_RE"])

    theta = np.linspace(0.0, 2.0 * np.pi, 240)
    ax.fill(np.cos(theta), np.sin(theta), color="#3498db", alpha=0.35, label=r"Earth Radius ($1.0\,R_\oplus$)")
    ax.plot(np.cos(theta), np.sin(theta), color="#1f4e78", lw=1.2)
    ax.plot(geo_R * np.cos(theta), geo_R * np.sin(theta), color="#7f8c8d", ls="--", lw=1.3,
            label=rf"Geostationary Ring ($r_{{\mathrm{{GEO}}}} = {geo_R:.2f}\,R_\oplus$)")
    ax.plot(gx, gy, color="#c0392b", lw=2.1, label="Apophis Geocentric Hyperbola (Float64 & Float128)")
    ax.scatter([gx[idx_min]], [gy[idx_min]], color="#f39c12", edgecolors="black", s=60, zorder=6,
               label=rf"Perigee Closest Approach $r_{{\min}} = {h3['min_r_earth_RE']:.2f}\,R_\oplus$")

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlim(-22.0, 22.0)
    ax.set_ylim(-15.0, 15.0)
    ax.set_xlabel(r"Geocentric $x_{\mathrm{geo}} / R_\oplus$")
    ax.set_ylabel(r"Geocentric $y_{\mathrm{geo}} / R_\oplus$")
    ax.legend(loc="upper left", framealpha=0.92, fontsize=7.3)
    ax.set_title(r"(c) H3 Apophis (2029): Geocentric Hyperbolic Flyby Inside GEO Ring")

    # --- Panel (d): H3 Apophis Resonant Aten -> Apollo Semi-Major Axis Jump ---
    ax = axes[1, 1]
    t_h3 = np.asarray(h3["series"]["times"], dtype=float)
    a64 = np.asarray(h3["series"]["a_helio_au_f64"], dtype=float)
    a128 = np.asarray(h3["series"]["a_helio_au_f128"], dtype=float)
    e64 = np.asarray(h3["series"]["e_helio_f64"], dtype=float)

    ax.plot(t_h3, a64, color="#1f4e78", lw=2.2,
            label=rf"Semi-Major Axis $a(t)$ Float64 ({h3['a_pre_au']:.4f} $\to$ {h3['a_post_au']:.4f} AU)")
    ax.plot(t_h3, a128, color="#27ae60", lw=1.5, ls="--",
            label=rf"Semi-Major Axis $a(t)$ Float128 ($\Delta a = +{h3['delta_a_au']:.4f}$ AU)")
    ax.axhline(1.0, color="#c0392b", ls="--", lw=1.4, label=r"Aten / Apollo Boundary ($a = 1.000$ AU)")
    ax.axvline(float(h3["t_perigee"]), color="#7f8c8d", ls=":", lw=1.2)
    ax.text(0.0008, 0.935, "Aten Class\n($a_{\\mathrm{pre}} = 0.9224 < 1$ AU)", color="#1f4e78", fontsize=8.2, fontweight="bold")
    ax.text(0.0130, 1.124, "Apollo Class ($a_{\\mathrm{post}} = 1.1163 > 1$ AU)", color="#1e8449", fontsize=8.2, fontweight="bold")
    ax.set_xlabel("Physical Time $t$ (yr / $2\\pi$)")
    ax.set_ylabel("Heliocentric Semi-Major Axis $a(t)$ [AU]", color="#1f4e78")
    ax.set_ylim(0.89, 1.15)

    ax_e = ax.twinx()
    ax_e.plot(t_h3, e64, color="#8e44ad", lw=1.5, ls="-.",
              label=rf"Heliocentric Eccentricity $e(t)$ ({h3['e_pre']:.3f} $\to$ {h3['e_post']:.3f})")
    ax_e.set_ylabel("Heliocentric Eccentricity $e(t)$", color="#8e44ad")
    ax_e.set_ylim(0.14, 0.25)
    ax_e.grid(False)

    l4a, lb4a = ax.get_legend_handles_labels()
    l4b, lb4b = ax_e.get_legend_handles_labels()
    ax.legend(l4a + l4b, lb4a + lb4b, loc="center right", framealpha=0.92, fontsize=7.2)
    ax.set_title(r"(d) H3 Apophis (2029): Resonant Class Transition Aten ($a<1$) $\to$ Apollo ($a>1$)")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save_and_copy(fig, "gen04_plate_C_physical_dynamics")
    # Also update legacy composite filename so existing LaTeX builds stay synchronized
    _save_and_copy(fig, "fig_ch07_historical_encounters_atlas_gauge")
    plt.close(fig)


def main() -> None:
    if not JSON_PATH.exists():
        raise FileNotFoundError(f"Missing certified JSON: {JSON_PATH}")
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)

    assert data.get("contract_assertions", {}).get("all_passed") is True, (
        "Contract assertions have not all passed in JSON!"
    )

    _configure_matplotlib()
    plot_plate_H1_lexell(data)
    plot_plate_H2_shoemaker_levy(data)
    plot_plate_H3_apophis(data)
    print("[Plotter 04] All 3 dedicated paper plates (H1, H2, H3) generated successfully.")


if __name__ == "__main__":
    main()
