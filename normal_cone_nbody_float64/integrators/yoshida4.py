"""
Compatibility Re-Export Module for the 4th-Order 3-Stage Symmetric Yoshida Composition S_4.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.2 (The Sheng--Suzuki Barrier and High-Order Symmetric Compositions),
  Equations (5.43)-(5.60).
- Re-exports Yoshida4Integrator and certified stage weights from .yoshida.
"""

from .yoshida import (
    Yoshida4Integrator,
    YOSHIDA_W1,
    YOSHIDA_W0,
    YOSHIDA_WEIGHTS,
    YOSHIDA_DRIFT_WEIGHTS,
    YOSHIDA_IMPULSE_WEIGHTS,
)

__all__ = [
    "Yoshida4Integrator",
    "YOSHIDA_W1",
    "YOSHIDA_W0",
    "YOSHIDA_WEIGHTS",
    "YOSHIDA_DRIFT_WEIGHTS",
    "YOSHIDA_IMPULSE_WEIGHTS",
]
