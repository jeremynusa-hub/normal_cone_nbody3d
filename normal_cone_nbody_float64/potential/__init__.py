"""
Layer 3: Potential Package (`normal_cone_nbody_float64.potential`).

Provides regularized conormal potential fields, smooth scale-adapted screening,
intrinsic Hill tidal supernode detection, and unsoftened Newtonian baselines:
1. `screening.py`: Scale-adapted conormal screening function Theta(u_v) (Hermite C^2 & Cauchy C^infty).
2. `screened_apex.py`: Regularized C^{1,1} hierarchical apex potential and symplectic kick operator.
3. `supernode_detector.py`: Autonomous Hill tidal supernode detector and Marchal-Sundman thresholds.
4. `direct.py`: Unsoftened Newtonian potential and direct baseline force evaluator.

Theoretical Citations:
- Monograph: Chapter 1 and Sections 4.1, 4.2, 4.3.
- Formal Lean 4: 01_ApexUniversalLipschitzBound_Theorems.lean, 01_PartitionLattice_Defs.lean,
  and PaperUnifiedFoundations.lean.
"""

from .screening import (
    eval_screening_theta,
    eval_screening_theta_vectorized,
    get_screening_mode_for_order,
)
from .screened_apex import (
    ScreenedApexPotential,
    compute_quadrupole_tensor,
    compute_quadrupole_potential,
)
from .optimized_screened_apex import (
    OptimizedScreenedApexPotential,
)
from .supernode_detector import (
    compute_tidal_ratio,
    compute_natural_threshold,
    detect_supernodes,
)
from .direct import (
    compute_direct_newtonian_potential,
    compute_direct_newtonian_forces,
    direct_newtonian_step,
)

__all__ = [
    "eval_screening_theta",
    "eval_screening_theta_vectorized",
    "get_screening_mode_for_order",
    "ScreenedApexPotential",
    "OptimizedScreenedApexPotential",
    "compute_quadrupole_tensor",
    "compute_quadrupole_potential",
    "compute_tidal_ratio",
    "compute_natural_threshold",
    "detect_supernodes",
    "compute_direct_newtonian_potential",
    "compute_direct_newtonian_forces",
    "direct_newtonian_step",
]
