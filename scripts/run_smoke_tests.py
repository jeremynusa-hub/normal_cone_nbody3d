#!/usr/bin/env python3
"""
Tier 2 Reproducibility Orchestrator: Full Smoke Test & Formal Audit (< 2 minutes)
=================================================================================

Companion Repository:
  Normal-Cone Is All You Need: Unsoftened Fulton-MacPherson Symplectic Integration
  of the 3D Newtonian N-Body Problem (v1.0.0)
  Author: Jérémy Nusa

Purpose:
  Executes in sequence:
    1. The Lean 4 Formal Verification & Axiom Hygiene Audit (`lean4/scripts/verify_lean4_suite.py`
       across all 35 Appendix A modules and 5 certified chapter compilation suites).
    2. The Appendix A Standalone Showcase Interface Test Suite
       (`tests/unit/test_paper_showcase_interfaces.py`, including Arenstorf PCR3BP).
    3. The Full IEEE-754 Float64 Unit Test Suite (`tests/unit/`).
    4. The Full 116-bit / 35-digit Float128 Quadruple-Precision Mirror Unit Test Suite
       (`tests/unit_float128/`).
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path


def run_step(step_num: int, total_steps: int, title: str, cmd: list[str], cwd: Path) -> bool:
    print("-" * 88)
    print(f"[{step_num}/{total_steps}] {title}")
    print(f"  Command: {' '.join(cmd)}")
    t0 = time.perf_counter()
    res = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    dt = time.perf_counter() - t0
    if res.returncode != 0:
        print(res.stdout)
        print(res.stderr, file=sys.stderr)
        print(f"  -> [FAIL] ({dt:.2f}s)", file=sys.stderr)
        return False
    # Print concise tail summary
    lines = [line.strip() for line in res.stdout.splitlines() if line.strip()]
    tail = lines[-1] if lines else "OK"
    print(f"  -> [PASS] ({dt:.2f}s) | {tail}")
    return True


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    repo_root = Path(__file__).resolve().parent.parent
    t_start = time.perf_counter()

    print("=" * 88)
    print("TIER 2 REPRODUCIBILITY: LEAN 4 FORMAL AUDIT & DUAL-PRECISION PYTHON TEST SUITE")
    print("=" * 88)

    steps = [
        (
            "Lean 4 Formal Proof Suite & 0-sorryAx Audit (35 Appendix A Modules, §A.1–§A.20)",
            [sys.executable, "lean4/scripts/verify_lean4_suite.py"],
        ),
        (
            "Appendix A Standalone Reference Interfaces & PCR3BP Arenstorf Orbit (§3.3 & §6.1–§6.6)",
            [sys.executable, "-m", "pytest", "tests/unit/test_paper_showcase_interfaces.py", "-q"],
        ),
        (
            "Production IEEE-754 Float64 Symplectic Library Unit Test Suite (38 suites)",
            [sys.executable, "-m", "pytest", "tests/unit", "-q"],
        ),
        (
            "Mirror 116-bit / 35-digit Float128 Quadruple-Precision Library Unit Test Suite (8 suites)",
            [sys.executable, "-m", "pytest", "tests/unit_float128", "-q"],
        ),
    ]

    for idx, (title, cmd) in enumerate(steps, start=1):
        if not run_step(idx, len(steps), title, cmd, repo_root):
            print("=" * 88)
            print("TIER 2 VERDICT: FAIL", file=sys.stderr)
            return 1

    total_dt = time.perf_counter() - t_start
    print("=" * 88)
    print(f"TIER 2 VERDICT: PASS (All 4 verification stages succeeded in {total_dt:.2f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
