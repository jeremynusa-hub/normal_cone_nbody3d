"""
03_DyadicScheduler_Interface.py
Standalone Python Interface for the Discrete Dyadic Multi-Rate Scheduler (Phase 2).

Colocalized with Chapter 6, Section 6.3 (Appendix A.24):
    "Macro-Time Loop, Dynamical Timescales, and Skeel--Gear Dyadic Scheduler"

Reference Implementation Interface:
Re-exports and delegates directly to the canonical reference scheduler in
normal_cone_nbody_float64/system/dyadic_scheduler.py:
1. Discrete dyadic time lattice subdivision: dt_k = Delta t_macro / 2^k.
2. Symmetric endpoint level evaluation:
       k_sym(z_n, dt) = max(k(z_n), k(Phi_dt^{(0)}(z_n))),
   guaranteeing exact time-reversal symmetry k_fwd(n) == k_bwd(N-1-n) across
   pericenter flybys and encounters.
3. Operational realization of the Leaf/Apex Decoupling Theorem (Theorems 5.11 & 6.2)
   with bilateral pericenter guard |t_peri| <= 2|dt| active on all unbound/hyperbolic
   encounters.
4. Bit-exact timestamp generation eliminating IEEE 754 summation drift.
"""

from __future__ import annotations
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from normal_cone_nbody_float64.system.dyadic_scheduler import (
    DEFAULT_MAX_DYADIC_LEVEL,
    ColoredInterval,
    DyadicScheduler,
)

DEFAULT_SAFETY_FACTOR: float = 0.05

__all__ = [
    "DEFAULT_MAX_DYADIC_LEVEL",
    "DEFAULT_SAFETY_FACTOR",
    "ColoredInterval",
    "DyadicScheduler",
]

