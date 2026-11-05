"""
Modular Symplectic Integrator Suite on the Fulton--MacPherson Configuration Space.

Mathematical Monograph Alignment:
- Chapter 5: Higher-Order Symplectic Integration and Baker--Campbell--Hausdorff
  Palindromic Compositions (Sections 5.1–5.5).
- Formal Lean 4 Certification:
  * PaperUnifiedFoundations.lean
  * BM6Order6SymmetricComposition.lean
  * 02_LieBargmannExactInvariance_Theorems.lean
  * 01_BlanesMoanShadowConfinement_Theorems.lean
  * 01_CorrectorBM4_Theorems.lean
  * 02_CorrectorBM6_Theorems.lean.
"""

from .protocol import SymplecticIntegrator, IntegrationStepResult
from .strang import StrangIntegrator
from .yoshida import Yoshida4Integrator
from .blanes_moan_4 import BlanesMoan4Integrator
from .blanes_moan_6 import BlanesMoan6Integrator
from .correctors import (
    SymplecticCorrector,
    CorrectorC4,
    GlobalCorrectorC4,
    compute_analytical_directional_force_derivative,
    compute_conormal_corrector_kick,
    compute_global_conormal_c4_kick,
)
from .processed import (
    ProcessedSymplecticIntegrator,
    PBM4GaugeTransitionOperator,
    apply_gauge_connection,
    create_pbm4_integrator,
)

__all__ = [
    "SymplecticIntegrator",
    "IntegrationStepResult",
    "StrangIntegrator",
    "Yoshida4Integrator",
    "BlanesMoan4Integrator",
    "BlanesMoan6Integrator",
    "SymplecticCorrector",
    "CorrectorC4",
    "GlobalCorrectorC4",
    "compute_analytical_directional_force_derivative",
    "compute_conormal_corrector_kick",
    "compute_global_conormal_c4_kick",
    "ProcessedSymplecticIntegrator",
    "PBM4GaugeTransitionOperator",
    "apply_gauge_connection",
    "create_pbm4_integrator",
]


