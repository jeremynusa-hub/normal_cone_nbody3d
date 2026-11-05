#!/usr/bin/env python3
"""
Formal Verification & Axiom Hygiene Auditor for the Lean 4 Proof Suite
======================================================================

Companion Repository:
  Unsoftened Normal-Cone Symplectic Integration of the 3D Newtonian N-Body Problem (v1.0.0)
  Author: Jérémy Nusa

Purpose:
  1. Verifies the presence and SHA-256 integrity of all 35 formal Lean 4 modules
     corresponding 1-to-1 with Appendix A (§A.1 to §A.20) of the monograph.
  2. Verifies the 5 self-contained chapter compilation suites under
     `NormalConeNBody/CertifiedSuites/` (`Sec01_Introduction.lean` to
     `Sec05_SymplecticSchemes.lean`).
  3. Audits every `.lean` module for complete absence of `sorry`, `sorryAx`,
     `admit`, or ad-hoc `axiom` declarations outside comments.
  4. Counts all formal declarations (`def`, `abbrev`, `structure`, `inductive`,
     `lemma`, `theorem`) and `#print axioms` audit statements.
  5. Optionally runs `lake build` or `lake env lean` when `--compile` is passed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Canonical mapping of the 35 Appendix A (§A.1 - §A.20) Lean 4 modules:
# (appendix_section, repo_rel_path, manuscript_rel_path, description)
LEAN4_MODULES_MANIFEST: List[Tuple[str, str, str, str]] = [
    # Chapter 1: §A.1 - §A.7 (7 modules)
    (
        "§A.1",
        "NormalConeNBody/Ch01_Introduction/01_ConfigurationSpace_Defs.lean",
        "chapters/ch01_introduction/sec01_background/code/01_ConfigurationSpace_Defs.lean",
        "Euclidean Configuration Space and Center-of-Mass Reduction",
    ),
    (
        "§A.2",
        "NormalConeNBody/Ch01_Introduction/02_CollisionDivisor_Defs.lean",
        "chapters/ch01_introduction/sec01_background/code/02_CollisionDivisor_Defs.lean",
        "Pairwise Collision Diagonals and Singular Locus Definitions",
    ),
    (
        "§A.3",
        "NormalConeNBody/Ch01_Introduction/03_CollisionDivisor_Theorems.lean",
        "chapters/ch01_introduction/sec01_background/code/03_CollisionDivisor_Theorems.lean",
        "Irreducibility and Complement Duality of the Collision Locus",
    ),
    (
        "§A.4",
        "NormalConeNBody/Ch01_Introduction/04_SymplecticForm_Defs.lean",
        "chapters/ch01_introduction/sec01_background/code/04_SymplecticForm_Defs.lean",
        "Canonical Cotangent Bundle and Darboux 2-Form Definitions",
    ),
    (
        "§A.5",
        "NormalConeNBody/Ch01_Introduction/05_SymplecticForm_Theorems.lean",
        "chapters/ch01_introduction/sec01_background/code/05_SymplecticForm_Theorems.lean",
        "Skew-Symmetry and Global Non-Degeneracy of the Darboux 2-Form",
    ),
    (
        "§A.6",
        "NormalConeNBody/Ch01_Introduction/06_BargmannAlgebra_Defs.lean",
        "chapters/ch01_introduction/sec01_background/code/06_BargmannAlgebra_Defs.lean",
        "Ten Classical Integrals and Galilean-Bargmann Generators",
    ),
    (
        "§A.7",
        "NormalConeNBody/Ch01_Introduction/07_BargmannAlgebra_Theorems.lean",
        "chapters/ch01_introduction/sec01_background/code/07_BargmannAlgebra_Theorems.lean",
        "Exact Lie-Poisson Commutation Relations of the Bargmann Algebra",
    ),
    # Chapter 2: §A.8 - §A.10 (6 modules)
    (
        "§A.8",
        "NormalConeNBody/Ch02_NormalCone/01_ReesDeformation_Defs.lean",
        "chapters/ch02_conormal_foundations/sec01_rees_deformation/code/01_ReesDeformation_Defs.lean",
        "Deformation to the Normal Cone and Shifted Cotangent Lift Definitions",
    ),
    (
        "§A.8",
        "NormalConeNBody/Ch02_NormalCone/02_ReesDeformation_Theorems.lean",
        "chapters/ch02_conormal_foundations/sec01_rees_deformation/code/02_ReesDeformation_Theorems.lean",
        "Exact Symplecticity and Non-Degeneracy of the Normal-Cone Lift",
    ),
    (
        "§A.9",
        "NormalConeNBody/Ch02_NormalCone/03_CanonicalSplitting_Defs.lean",
        "chapters/ch02_conormal_foundations/sec02_algorithmic_recipe/code/03_CanonicalSplitting_Defs.lean",
        "Four-Step Normal-Cone Symplectic Splitting Definitions",
    ),
    (
        "§A.9",
        "NormalConeNBody/Ch02_NormalCone/04_CanonicalSplitting_Theorems.lean",
        "chapters/ch02_conormal_foundations/sec02_algorithmic_recipe/code/04_CanonicalSplitting_Theorems.lean",
        "Exact Casimir Invariance and Time-Reversibility of the Splitting",
    ),
    (
        "§A.10",
        "NormalConeNBody/Ch02_NormalCone/05_HeavyTop_Defs.lean",
        "chapters/ch02_conormal_foundations/sec03_retrospective_cases/code/05_HeavyTop_Defs.lean",
        "Heavy Top Lie-Poisson Normal-Cone Fiber Bundle Definitions",
    ),
    (
        "§A.10",
        "NormalConeNBody/Ch02_NormalCone/06_HeavyTop_Theorems.lean",
        "chapters/ch02_conormal_foundations/sec03_retrospective_cases/code/06_HeavyTop_Theorems.lean",
        "Exact Preservation of Both Heavy Top Casimirs under Rodrigues Flow",
    ),
    # Chapter 3: §A.11 - §A.12 (2 modules)
    (
        "§A.11",
        "NormalConeNBody/Ch03_ThreeBody/01_PCR3BPHamiltonian_Theorems.lean",
        "chapters/ch03_three_body_paradigm/sec01_pcr3bp_canonical/code/01_PCR3BPHamiltonian_Theorems.lean",
        "Synodic PCR3BP Jacobi Conservation and Symplectic Coriolis Propagator",
    ),
    (
        "§A.12",
        "NormalConeNBody/Ch03_ThreeBody/02_PCR3BPRegularization_Theorems.lean",
        "chapters/ch03_three_body_paradigm/sec02_classical_regularizations/code/02_PCR3BPRegularization_Theorems.lean",
        "Levi-Civita Quadratic Normal-Cone Covering and Conformal Scaling",
    ),
    # Chapter 4: §A.13 - §A.15 (10 modules)
    (
        "§A.13",
        "NormalConeNBody/Ch04_FultonMacPherson/01_PartitionLattice_Defs.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec01_whitney_lattice/code/01_PartitionLattice_Defs.lean",
        "Partition Lattice Pi(N), Collision Codimension, and Non-Transversality",
    ),
    (
        "§A.13",
        "NormalConeNBody/Ch04_FultonMacPherson/02_ClusterHyperRadius_Theorems.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec01_whitney_lattice/code/02_ClusterHyperRadius_Theorems.lean",
        "Cluster Hyper-Radius and McGehee Blow-Up Potential Factorization",
    ),
    (
        "§A.13",
        "NormalConeNBody/Ch04_FultonMacPherson/03_MaximalFlagJacobiTree_Theorems.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec01_whitney_lattice/code/03_MaximalFlagJacobiTree_Theorems.lean",
        "Binary Jacobi Tree Hierarchy and Nested Fulton-MacPherson Flags",
    ),
    (
        "§A.13",
        "NormalConeNBody/Ch04_FultonMacPherson/04_DisjointCotangentCommutation_Theorems.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec01_whitney_lattice/code/04_DisjointCotangentCommutation_Theorems.lean",
        "Exact Cotangent Commutation of Disjoint Normal-Cone Leaves",
    ),
    (
        "§A.14",
        "NormalConeNBody/Ch04_FultonMacPherson/01_ApexUniversalLipschitzBound_Theorems.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec02_screened_apex_potential/code/01_ApexUniversalLipschitzBound_Theorems.lean",
        "Universal C^{1,1} Lipschitz and Hessian Bounds of the Screened Apex Potential",
    ),
    (
        "§A.14",
        "NormalConeNBody/Ch04_FultonMacPherson/02_LieBargmannExactInvariance_Theorems.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec02_screened_apex_potential/code/02_LieBargmannExactInvariance_Theorems.lean",
        "Exact Galilean-Bargmann Invariance of the Screened Apex Splitting",
    ),
    (
        "§A.14",
        "NormalConeNBody/Ch04_FultonMacPherson/03_ExtendedPalindromicSymplecticity_Theorems.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec02_screened_apex_potential/code/03_ExtendedPalindromicSymplecticity_Theorems.lean",
        "Extended Poincare-Sundman Symplecticity and Palindromic Time-Reversibility",
    ),
    (
        "§A.15",
        "NormalConeNBody/Ch04_FultonMacPherson/01_UnifiedSuperNodeDef_Theorems.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec03_synchronous_supernodes/code/01_UnifiedSuperNodeDef_Theorems.lean",
        "Synchronous Super-Node Collective Sundman Time-Step Positivity",
    ),
    (
        "§A.15",
        "NormalConeNBody/Ch04_FultonMacPherson/02_SuperNodeAntiZenoRegularization_Theorems.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec03_synchronous_supernodes/code/02_SuperNodeAntiZenoRegularization_Theorems.lean",
        "Anti-Zeno Physical Time Advance and Keplerian Normal-Cone Scaling",
    ),
    (
        "§A.15",
        "NormalConeNBody/Ch04_FultonMacPherson/03_AtlasTransitionSymplecticAntiChattering_Theorems.lean",
        "chapters/ch04_fulton_macpherson_splitting/sec03_synchronous_supernodes/code/03_AtlasTransitionSymplecticAntiChattering_Theorems.lean",
        "Symplectic SL(2,R) Chart Transitions and Hysteresis Anti-Chattering",
    ),
    # Chapter 5: §A.16 - §A.20 (10 modules)
    (
        "§A.16",
        "NormalConeNBody/Ch05_SymplecticSchemes/01_StrangOrderTwoSymplectic_Theorems.lean",
        "chapters/ch05_higher_order_symplectic/sec01_strang_propagator/code/01_StrangOrderTwoSymplectic_Theorems.lean",
        "Second-Order Palindromic Symplecticity of the Strang S2 Propagator",
    ),
    (
        "§A.16",
        "NormalConeNBody/Ch05_SymplecticSchemes/02_LieBargmannConservation_Theorems.lean",
        "chapters/ch05_higher_order_symplectic/sec01_strang_propagator/code/02_LieBargmannConservation_Theorems.lean",
        "Exact Conservation of the Ten Lie-Bargmann Generators under S2",
    ),
    (
        "§A.17",
        "NormalConeNBody/Ch05_SymplecticSchemes/01_YoshidaOrderFourSuzuki_Theorems.lean",
        "chapters/ch05_higher_order_symplectic/sec02_suzuki_barrier_mcgehee/code/01_YoshidaOrderFourSuzuki_Theorems.lean",
        "Sheng-Suzuki Barrier and Order-Four Yoshida Triple Jump Algebraic Equations",
    ),
    (
        "§A.18",
        "NormalConeNBody/Ch05_SymplecticSchemes/01_BlanesMoanShadowConfinement_Theorems.lean",
        "chapters/ch05_higher_order_symplectic/sec03_bm4_seven_stages/code/01_BlanesMoanShadowConfinement_Theorems.lean",
        "Seven-Stage Blanes-Moan BM4 FSAL Consistency and Shadow Confinement",
    ),
    (
        "§A.19",
        "NormalConeNBody/Ch05_SymplecticSchemes/01_BlanesMoanOrderSix_Theorems.lean",
        "chapters/ch05_higher_order_symplectic/sec04_bm6_eleven_stages/code/01_BlanesMoanOrderSix_Theorems.lean",
        "Eleven-Stage Blanes-Moan BM6 Order-Six Consistency and RKN Commutators",
    ),
    (
        "§A.20",
        "NormalConeNBody/Ch05_SymplecticSchemes/01_CorrectorBM4_Theorems.lean",
        "chapters/ch05_higher_order_symplectic/sec05_symplectic_correctors/code/01_CorrectorBM4_Theorems.lean",
        "Symplectic Corrector Conjugation C_h^{-1} o BM4 o C_h and Centered Equivalence",
    ),
    (
        "§A.20",
        "NormalConeNBody/Ch05_SymplecticSchemes/03_ConormalCorrectorPBM4_Order6.lean",
        "chapters/ch05_higher_order_symplectic/03_ConormalCorrectorPBM4_Order6.lean",
        "Single-Chart Sixth-Order Cancellation of [K,[K,V]] and [V,[V,K]] via PBM4",
    ),
    (
        "§A.20",
        "NormalConeNBody/Ch05_SymplecticSchemes/06_ConormalGaugeTransitionPBM4_Order6.lean",
        "chapters/ch05_higher_order_symplectic/06_ConormalGaugeTransitionPBM4_Order6.lean",
        "Inter-Chart Symplectic Gauge Connection T_{A->B}^{(PBM4)} and Cocycle Identity",
    ),
    (
        "§A.20",
        "NormalConeNBody/Ch05_SymplecticSchemes/07_GlobalConormalCorrectorPBM4_Order6.lean",
        "chapters/ch05_higher_order_symplectic/07_GlobalConormalCorrectorPBM4_Order6.lean",
        "Global Atlas Gauge-Invariance Theorem on the P_0-Dominated Regime",
    ),
    (
        "§A.20",
        "NormalConeNBody/Ch05_SymplecticSchemes/08_UniversalOrderSixObstruction.lean",
        "chapters/ch05_higher_order_symplectic/08_UniversalOrderSixObstruction.lean",
        "Bidegree-(3,2) Hall Lie Algebra L(5,2) Cohomological Obstruction Beyond P_0",
    ),
]

CERTIFIED_CHAPTER_SUITES: List[Tuple[str, str, str]] = [
    ("Ch01 (§A.1-§A.7)", "NormalConeNBody/CertifiedSuites/Sec01_Introduction.lean", "Sec01_Round2_Audit.lean"),
    ("Ch02 (§A.8-§A.10)", "NormalConeNBody/CertifiedSuites/Sec02_NormalCone.lean", "Sec02_Round2_Audit.lean"),
    ("Ch03 (§A.11-§A.12)", "NormalConeNBody/CertifiedSuites/Sec03_ThreeBody.lean", "Sec03_Round2_Audit.lean"),
    ("Ch04 (§A.13-§A.15)", "NormalConeNBody/CertifiedSuites/Sec04_FultonMacPherson.lean", "Sec04_Round2_Audit.lean"),
    ("Ch05 (§A.16-§A.20)", "NormalConeNBody/CertifiedSuites/Sec05_SymplecticSchemes.lean", "Sec05_Round2_Audit.lean"),
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def strip_lean_comments(text: str) -> str:
    """Remove block comments /- ... -/ and line comments -- ... for syntactic auditing."""
    # Remove nested or non-nested /- ... -/ comments iteratively
    prev = None
    cleaned = text
    while prev != cleaned:
        prev = cleaned
        cleaned = re.sub(r"/-((?!/-)[\s\S])*?-/", "", cleaned)
    cleaned = re.sub(r"--[^\n]*", "", cleaned)
    return cleaned


def count_declarations(text: str) -> Dict[str, int]:
    cleaned = strip_lean_comments(text)
    counts = {
        "def": len(re.findall(r"^\s*(?:noncomputable\s+|@\[[^\]]+\]\s*)*def\s+\S+", cleaned, flags=re.M)),
        "abbrev": len(re.findall(r"^\s*(?:noncomputable\s+)?abbrev\s+\S+", cleaned, flags=re.M)),
        "structure": len(re.findall(r"^\s*(?:@\[[^\]]+\]\s*)*structure\s+\S+", cleaned, flags=re.M)),
        "inductive": len(re.findall(r"^\s*inductive\s+\S+", cleaned, flags=re.M)),
        "lemma": len(re.findall(r"^\s*(?:@\[[^\]]+\]\s*)*lemma\s+\S+", cleaned, flags=re.M)),
        "theorem": len(re.findall(r"^\s*(?:@\[[^\]]+\]\s*)*theorem\s+\S+", cleaned, flags=re.M)),
    }
    counts["total"] = sum(counts.values())
    return counts


def audit_forbidden_tokens(text: str) -> List[str]:
    cleaned = strip_lean_comments(text)
    violations: List[str] = []
    for token in ("sorry", "sorryAx", "admit"):
        if re.search(rf"\b{token}\b", cleaned):
            violations.append(token)
    if re.search(r"^\s*axiom\s+\S+", cleaned, flags=re.M):
        violations.append("axiom")
    return violations


def sync_from_workspace(lean4_root: Path, manuscript_root: Path, audit_suite_dir: Path | None) -> None:
    for _, repo_rel, ms_rel, _ in LEAN4_MODULES_MANIFEST:
        src = manuscript_root / ms_rel
        dst = lean4_root / repo_rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)

    if audit_suite_dir and audit_suite_dir.exists():
        for _, suite_rel, audit_name in CERTIFIED_CHAPTER_SUITES:
            src = audit_suite_dir / audit_name
            dst = lean4_root / suite_rel
            if src.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="Verify Lean 4 formal proof suite (35 modules, 0 sorryAx).")
    parser.add_argument("--sync", action="store_true", help="Sync 35 modules from manuscript directory before verifying.")
    parser.add_argument("--manuscript-dir", type=str, default=None, help="Path to redaction/unsoftened_conormal_nbody3d.")
    parser.add_argument("--audit-suite-dir", type=str, default=None, help="Path to directory containing Sec01..05_Round2_Audit.lean.")
    parser.add_argument("--compile", action="store_true", help="Invoke `lake build` to compile the CertifiedSuites.")
    args = parser.parse_args()

    lean4_root = Path(__file__).resolve().parent.parent
    default_ms = lean4_root.parent.parent / "redaction" / "unsoftened_conormal_nbody3d"
    manuscript_root = Path(args.manuscript_dir) if args.manuscript_dir else default_ms
    audit_suite_dir = Path(args.audit_suite_dir) if args.audit_suite_dir else None

    if args.sync:
        if not manuscript_root.exists():
            print(f"[ERROR] Manuscript directory not found for --sync: {manuscript_root}", file=sys.stderr)
            return 1
        sync_from_workspace(lean4_root, manuscript_root, audit_suite_dir)
        print(f"[SYNC] Synchronized 35 Appendix A modules from {manuscript_root}")

    print("=" * 88)
    print("LEAN 4 FORMAL VERIFICATION & AXIOM HYGIENE AUDIT (APPENDIX A: §A.1 - §A.20)")
    print("=" * 88)

    total_decls = {"def": 0, "abbrev": 0, "structure": 0, "inductive": 0, "lemma": 0, "theorem": 0, "total": 0}
    all_passed = True
    ms_available = manuscript_root.exists()
    sha256_matches = 0

    for idx, (sec, repo_rel, ms_rel, desc) in enumerate(LEAN4_MODULES_MANIFEST, start=1):
        mod_path = lean4_root / repo_rel
        if not mod_path.exists():
            print(f"  [{idx:02d}/35] [FAIL] Missing module: {repo_rel}")
            all_passed = False
            continue

        text = mod_path.read_text(encoding="utf-8")
        violations = audit_forbidden_tokens(text)
        counts = count_declarations(text)
        for k in total_decls:
            total_decls[k] += counts[k]

        sha_ok = True
        if ms_available:
            ms_path = manuscript_root / ms_rel
            if ms_path.exists():
                if sha256_file(mod_path) == sha256_file(ms_path):
                    sha256_matches += 1
                else:
                    sha_ok = False
                    all_passed = False

        status = "PASS" if (not violations and sha_ok) else "FAIL"
        if violations:
            all_passed = False
        print(
            f"  [{idx:02d}/35] [{status}] {sec:<6} | {Path(repo_rel).name:<52} | "
            f"decls={counts['total']:2d} (thm/lem={counts['theorem'] + counts['lemma']:2d}) | "
            f"sorry/axiom={len(violations)}"
        )

    print("-" * 88)
    print("CERTIFIED CHAPTER COMPILATION SUITES (WITH #print axioms):")
    total_print_axioms = 0
    for label, suite_rel, _ in CERTIFIED_CHAPTER_SUITES:
        suite_path = lean4_root / suite_rel
        if not suite_path.exists():
            print(f"  [FAIL] Missing chapter suite: {suite_rel}")
            all_passed = False
            continue
        stext = suite_path.read_text(encoding="utf-8")
        sviol = audit_forbidden_tokens(stext)
        p_axioms = len(re.findall(r"^\s*#print\s+axioms\s+\S+", stext, flags=re.M))
        total_print_axioms += p_axioms
        if sviol:
            all_passed = False
        print(f"  [{'PASS' if not sviol else 'FAIL'}] {label:<20} | {suite_rel:<54} | #print axioms={p_axioms:2d}")

    print("-" * 88)
    print(f"SUMMARY OF FORMAL DECLARATIONS ACROSS THE 35 MODULES:")
    print(
        f"  Total Declarations : {total_decls['total']} "
        f"(def={total_decls['def']}, abbrev={total_decls['abbrev']}, "
        f"structure={total_decls['structure']}, inductive={total_decls['inductive']}, "
        f"lemma={total_decls['lemma']}, theorem={total_decls['theorem']})"
    )
    print(f"  Total #print axioms checks in CertifiedSuites : {total_print_axioms}")
    if ms_available:
        print(f"  Bit-for-bit SHA-256 identity with Appendix A  : {sha256_matches}/35 modules")

    if args.compile:
        print("-" * 88)
        print("RUNNING `lake build` ON NormalConeNBody...")
        res = subprocess.run(["lake", "build"], cwd=str(lean4_root), capture_output=True, text=True)
        print(res.stdout)
        if res.returncode != 0:
            print(res.stderr, file=sys.stderr)
            all_passed = False
        elif "sorryAx" in res.stdout:
            print("[FAIL] `sorryAx` detected in `#print axioms` output!", file=sys.stderr)
            all_passed = False
        else:
            print("[PASS] `lake build` succeeded with 0 errors and 0 `sorryAx`.")

    print("=" * 88)
    if all_passed:
        print("FINAL VERDICT: PASS (35/35 Lean 4 modules certified, 0 sorry, 0 sorryAx, 0 custom axioms)")
        return 0
    else:
        print("FINAL VERDICT: FAIL", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
