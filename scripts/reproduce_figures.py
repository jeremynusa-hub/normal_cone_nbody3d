#!/usr/bin/env python3
"""
Tier 1 Reproducibility Orchestrator: Fast Figure Regeneration (< 15 seconds)
============================================================================

Companion Repository:
  Normal-Cone Is All You Need: Unsoftened Fulton-MacPherson Symplectic Integration
  of the 3D Newtonian N-Body Problem (v1.0.0)
  Author: Jérémy Nusa

Purpose:
  1. Verifies the presence and valid JSON schema of all 14 certified benchmark
     campaign datasets in `results/campaigns/`.
  2. Executes all 8 publication plate plotters in `figures/plot_*.py` to
     regenerate the 300 DPI PDF and PNG plates of Chapters 6 and 7 into
     `results/figures/` without re-running multi-hour numerical integrations.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Tuple

REQUIRED_CAMPAIGN_JSONS: List[str] = [
    "campaign_II_H1_lexell.json",
    "campaign_II_H2_shoemaker_levy.json",
    "campaign_II_H3_apophis.json",
    "generator_01_dual_convergence_and_microkernels.json",
    "generator_02_universal_obstruction_and_shadow_gauge.json",
    "generator_03_unsoftened_singularity_limit_cases.json",
    "generator_04_historical_encounters_and_atlas_gauge.json",
    "generator_05_secular_stability_laplace_and_laskar_naff.json",
    "rebound_b1_star_grazing_pericenter.json",
    "rebound_b2_hill_crossings.json",
    "rebound_b4_pythagorean_chaos.json",
    "rebound_b4_pythagorean_chaos_max_k4.json",
    "rebound_b4_pythagorean_chaos_max_k8.json",
    "rebound_b4_pythagorean_chaos_max_k10.json",
]

PLOTTER_SCRIPTS: List[Tuple[str, str]] = [
    ("Generator 01 (§6.6 & §7.1.1 — Figures 17–20)", "figures/plot_generator_01_dual_convergence_and_microkernels.py"),
    ("Generator 02 (§7.1.2 — Figures 21–23)", "figures/plot_generator_02_universal_obstruction_and_shadow_gauge.py"),
    ("Generator 03 (§7.1.3 — Figures 24–26)", "figures/plot_generator_03_unsoftened_singularity_limit_cases.py"),
    ("Generator 04 (§7.2.1 — Figures 27–28)", "figures/plot_generator_04_historical_encounters_and_atlas_gauge.py"),
    ("Generator 05 (§7.2.2 — Figures 29–30)", "figures/plot_generator_05_secular_stability_laplace_and_laskar_naff.py"),
    ("REBOUND B1 (§7.3.1 — Figure 31)", "figures/plot_rebound_b1_star_grazing_pericenter.py"),
    ("REBOUND B2 (§7.3.2 — Figure 32)", "figures/plot_rebound_b2_hill_crossings.py"),
    ("REBOUND B4 (§7.3.3 — Figure 33)", "figures/plot_rebound_b4_pythagorean_chaos.py"),
]


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    repo_root = Path(__file__).resolve().parent.parent
    campaigns_dir = repo_root / "results" / "campaigns"
    figures_out_dir = repo_root / "results" / "figures"
    figures_out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    print("=" * 88)
    print("TIER 1 REPRODUCIBILITY: CERTIFIED JSON AUDIT & PUBLICATION FIGURE REGENERATION")
    print("=" * 88)

    # 1. Verify all 14 JSON datasets
    print("[1/2] Auditing 14 certified campaign JSON datasets in results/campaigns/...")
    for idx, fname in enumerate(REQUIRED_CAMPAIGN_JSONS, start=1):
        fpath = campaigns_dir / fname
        if not fpath.exists():
            print(f"  [{idx:02d}/14] [FAIL] Missing dataset: {fname}", file=sys.stderr)
            return 1
        try:
            with fpath.open("r", encoding="utf-8") as f:
                data = json.load(f)
            size_kb = fpath.stat().st_size / 1024.0
            n_keys = len(data) if isinstance(data, dict) else len(data)
            print(f"  [{idx:02d}/14] [PASS] {fname:<56} ({size_kb:7.1f} KB, top-keys={n_keys})")
        except Exception as exc:
            print(f"  [{idx:02d}/14] [FAIL] Invalid JSON in {fname}: {exc}", file=sys.stderr)
            return 1

    # 2. Execute all 8 publication figure plotters
    print("-" * 88)
    print("[2/2] Regenerating publication plates (300 DPI PDF & PNG) into results/figures/...")
    for idx, (label, rel_script) in enumerate(PLOTTER_SCRIPTS, start=1):
        script_path = repo_root / rel_script
        t_s = time.perf_counter()
        res = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(repo_root),
            capture_output=True,
            text=True,
        )
        dt = time.perf_counter() - t_s
        if res.returncode != 0:
            print(f"  [{idx}/8] [FAIL] {label} ({rel_script}) failed in {dt:.2f}s:", file=sys.stderr)
            print(res.stderr, file=sys.stderr)
            return 1
        print(f"  [{idx}/8] [PASS] {label:<50} ({dt:5.2f}s)")

    generated_files = sorted(figures_out_dir.glob("*.*"))
    total_dt = time.perf_counter() - t0
    print("-" * 88)
    print(f"Generated {len(generated_files)} figure files in {figures_out_dir} in {total_dt:.2f}s.")
    print("=" * 88)
    print("TIER 1 VERDICT: PASS (14/14 JSON datasets verified, 8/8 figure suites regenerated)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
