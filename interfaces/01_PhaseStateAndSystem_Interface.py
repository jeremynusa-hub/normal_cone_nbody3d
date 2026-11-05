"""
Unified Normal-Cone N-Body Integration: Core State and System Facade Interface.
Mathematical Monograph Alignment: Chapter 6, Section 6.1 (Appendix A.22).
Reference Codebase: normal_cone_nbody_float64/core/state.py
                   normal_cone_nbody_float64/system/nbody_system.py

Reference Implementation Interface:
Re-exports and delegates directly to the canonical reference classes:
1. ExtendedPhaseState: Immutable state in P_ext = T*Q[N] x R_t x R_E with
   exact barycentric reduction R_cm = 0, P_tot = 0.
2. NBodySystem: Unified 6-phase normal-cone simulation orchestrator with:
   - Shadow-manifold evolution active by default (in_shadow_mode=True);
   - Unobstructed palindromic Lie conjugator Z_4^{BM4} in P_1 (mode="homological")
     when processed integration is active;
   - 4-stage causal macro-boundary pipeline (Proposition 6.2):
     (1) Delaunay-Hill atlas adaptation (gamma=0.60, gamma_build=0.70),
     (2) Synchronous supernode detection (Hill/Roche hysteresis),
     (3) Symmetric dyadic scale evaluation k_sym = max(k_start, k_end),
     (4) Non-commutative covariant gauge connection T_{(A,h_A)->(B,h_B)}^{(PBM4)}.
"""

from __future__ import annotations
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.system.nbody_system import NBodySystem

__all__ = [
    "ExtendedPhaseState",
    "NBodySystem",
]


