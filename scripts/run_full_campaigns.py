#!/usr/bin/env python3
"""
Tier 3 Reproducibility Orchestrator: Full Numerical Benchmark Campaigns
=======================================================================

Companion Repository:
  Normal-Cone Is All You Need: Unsoftened Fulton-MacPherson Symplectic Integration
  of the 3D Newtonian N-Body Problem (v1.0.0)
  Author: Jérémy Nusa

Purpose:
  Provides a unified command-line entry point to re-execute any or all of the 8
  numerical benchmark generators (Chapters 6 & 7) that produce the certified JSON
  datasets in `results/campaigns/`, followed by automatic regeneration of all
  publication figures.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Tuple

CAMPAIGNS: Dict[str, Tuple[str, str]] = {
    "gen01": (
        "§6.6 & §7.1.1 — Dual-Precision Convergence & Micro-Kernel Scaling (Tables 3, 4, 6)",
        "benchmarks/run_generator_01_dual_convergence_and_microkernels.py",
    ),
    "gen02": (
        "§5.5 & §7.1.2 — Universal Obstruction & PBM4 Shadow Manifold (ST-1..ST-4)",
        "benchmarks/run_generator_02_universal_obstruction_and_shadow_gauge.py",
    ),
    "gen03": (
        "§7.1.3 — Unsoftened Singularity Traversals L1, L2, L3a, L3b, L4 (Table 7)",
        "benchmarks/run_generator_03_unsoftened_singularity_limit_cases.py",
    ),
    "gen04": (
        "§7.2.1 — Historical Planetary Encounters H1 Lexell, H2 SL9, H3 Apophis (Table 8)",
        "benchmarks/run_generator_04_historical_encounters_and_atlas_gauge.py",
    ),
    "gen05": (
        "§7.2.2 — Secular Stability S1 Laplace 4:2:1 & S2 264 Dyadic Junctions (Table 9)",
        "benchmarks/run_generator_05_secular_stability_laplace_and_laskar_naff.py",
    ),
    "rebound_b1": (
        "§7.3.1 — Head-to-Head vs REBOUND 5.2.2 B1 Star-Grazing Pericenter (Table 10)",
        "benchmarks/run_rebound_b1_star_grazing_pericenter.py",
    ),
    "rebound_b2": (
        "§7.3.2 — Head-to-Head vs REBOUND 5.2.2 B2 Deep Hill-Sphere Crossings (Table 10)",
        "benchmarks/run_rebound_b2_hill_crossings.py",
    ),
    "rebound_b4": (
        "§7.3.3 — Head-to-Head vs REBOUND 5.2.2 B4 Burrau 3-4-5 Reversibility (Table 11)",
        "benchmarks/run_rebound_b4_pythagorean_chaos.py",
    ),
}


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    parser = argparse.ArgumentParser(
        description="Tier 3 Reproducibility: Re-run full numerical benchmark campaigns (Chapters 6 & 7)."
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List all available benchmark campaigns and exit.",
    )
    parser.add_argument(
        "--campaign",
        type=str,
        default="all",
        choices=["all", *CAMPAIGNS.keys()],
        help="Select a specific campaign to execute or 'all' (default: all).",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent

    if args.list:
        print("Available Tier 3 Benchmark Campaigns:")
        for key, (desc, script) in CAMPAIGNS.items():
            print(f"  {key:<12} : {desc}\n                 ({script})")
        return 0

    selected = list(CAMPAIGNS.keys()) if args.campaign == "all" else [args.campaign]
    print("=" * 88)
    print(f"TIER 3 REPRODUCIBILITY: RUNNING {len(selected)} BENCHMARK CAMPAIGN(S)")
    print("=" * 88)

    for idx, key in enumerate(selected, start=1):
        desc, rel_script = CAMPAIGNS[key]
        print(f"[{idx}/{len(selected)}] Running {key}: {desc}")
        t0 = time.perf_counter()
        res = subprocess.run([sys.executable, rel_script], cwd=str(repo_root))
        dt = time.perf_counter() - t0
        if res.returncode != 0:
            print(f"[FAIL] Campaign {key} exited with code {res.returncode} after {dt:.1f}s.", file=sys.stderr)
            return res.returncode
        print(f"[PASS] Campaign {key} completed in {dt:.1f}s.")

    print("-" * 88)
    print("Regenerating all publication figures via scripts/reproduce_figures.py...")
    res_fig = subprocess.run([sys.executable, "scripts/reproduce_figures.py"], cwd=str(repo_root))
    return res_fig.returncode


if __name__ == "__main__":
    sys.exit(main())
