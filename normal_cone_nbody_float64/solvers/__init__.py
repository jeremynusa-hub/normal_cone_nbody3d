"""
Layer 2: Solvers Package (`normal_cone_nbody_float64.solvers`).

Provides universal, unsoftened analytical solvers and conormal regularizations:
1. `stumpff.py`: Universal Stumpff functions c_0(z) through c_5(z) and recurrence derivatives.
2. `kepler_sundman.py`: Universal two-body Sundman-Stumpff Kepler solver with cubic Halley iteration.
3. `supernode.py`: Anti-Zeno synchronous multi-body solver under unfolded Sundman clock dt_S = chi_S^2 dtau_S.

Theoretical Citations:
- Monograph: Sections 3.3, 3.4, 4.4 and Appendix A.
- Formal Lean 4: StumpffDynamics.lean, SuperNodeAntiZenoRegularization_Theorems.lean.
"""

from .stumpff import stumpff_c, stumpff_c_derivatives
from .kepler_sundman import solve_sundman_time_equation, step_kepler_sundman, compute_wronskian
from .supernode import step_supernode_synchronous, SupernodeCluster
from .optimized_kepler_sundman import (
    OptimizedKeplerSolver,
    step_kepler_sundman_fast,
    stumpff_c0123_scalar,
)
from .optimized_supernode import (
    OptimizedSupernodeSolver,
    step_supernode_synchronous_fast,
)
from .jit_supernode import (
    JitSupernodeSolver,
    step_supernode_synchronous_jit,
    NUMBA_AVAILABLE,
)

__all__ = [
    "stumpff_c",
    "stumpff_c_derivatives",
    "solve_sundman_time_equation",
    "step_kepler_sundman",
    "compute_wronskian",
    "step_supernode_synchronous",
    "SupernodeCluster",
    "OptimizedKeplerSolver",
    "step_kepler_sundman_fast",
    "stumpff_c0123_scalar",
    "OptimizedSupernodeSolver",
    "step_supernode_synchronous_fast",
    "JitSupernodeSolver",
    "step_supernode_synchronous_jit",
    "NUMBA_AVAILABLE",
]

