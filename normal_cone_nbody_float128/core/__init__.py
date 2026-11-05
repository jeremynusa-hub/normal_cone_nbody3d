"""
Core Layer: Canonical Algebra, Extended Phase Space, Invariants & Symplectic Foliations.

Mathematical Monograph Alignment:
- Chapter 2: The Deformation to the Normal Cone and Graded Rees Algebra (Sections 2.1 - 2.3).
- Chapter 3: The Three-Body Paradigm and Singular Desingularization (Sections 3.1 - 3.4).
- Chapter 4: The Partition Lattice Pi(N) and Fulton--MacPherson Manifold Q[N] (Section 4.1).
- Certified Formal Foundations: PaperUnifiedFoundations.lean.

Exports:
- Constants: G, ETA_CRIT, GAMMA_HYSTERESIS, DELTA_LAMBDA_HYSTERESIS, BM6 constants.
- ExtendedPhaseState: Immutable point on P_ext = T*Q[N] x R_t x R_E.
- LieBargmannInvariants: 10 generators of the extended Galilean Lie algebra.
- JacobiTree, JacobiTopologyNode: Immutable combinatorial hierarchy and mass metrics.
- JacobiCoordinates: Canonical relative coordinates and momenta.
- SL2TransitionMatrix: Symplectic lifts of Jacobi tree rotations in Sp(6N-6, R).
"""

from .constants import (
    G_DEFAULT,
    ETA_CRIT,
    ETA_ENTER_HILL,
    ETA_EXIT_HILL,
    GAMMA_HYSTERESIS,
    DELTA_LAMBDA_HYSTERESIS,
    EPS_MACH,
    DEFAULT_MAX_DYADIC_LEVEL,
    BM6_DRIFT_WEIGHTS,
    BM6_IMPULSE_WEIGHTS,
)
from .state import ExtendedPhaseState
from .invariants import (
    LieBargmannInvariants,
    compute_linear_momentum,
    compute_angular_momentum,
    compute_center_of_mass,
    compute_galilean_boost,
    compute_kinetic_energy,
    compute_newtonian_potential_energy,
    compute_hamiltonian,
    compute_cartesian_hamiltonian,
    compute_invariants,
    compute_truncation_residual,
)
from .jacobi import (
    JacobiTopologyNode,
    JacobiTree,
    build_canonical_tree,
    enumerate_all_jacobi_trees,
)
from .coordinates import (
    JacobiCoordinates,
    cartesian_to_jacobi,
    jacobi_to_cartesian,
    compute_jacobi_kinetic_energy,
    darboux_pairing_error,
    rees_deformation_retraction,
)
from .transitions import (
    SL2TransitionMatrix,
    build_symplectic_matrix_12,
    build_tree_transition_symplectic_12,
    rotate_jacobi_tree_coordinates,
    verify_symplectic_matrix_12,
)

__all__ = [
    "G_DEFAULT",
    "ETA_CRIT",
    "ETA_ENTER_HILL",
    "ETA_EXIT_HILL",
    "GAMMA_HYSTERESIS",
    "DELTA_LAMBDA_HYSTERESIS",
    "EPS_MACH",
    "DEFAULT_MAX_DYADIC_LEVEL",
    "BM6_DRIFT_WEIGHTS",
    "BM6_IMPULSE_WEIGHTS",
    "ExtendedPhaseState",
    "LieBargmannInvariants",
    "compute_linear_momentum",
    "compute_angular_momentum",
    "compute_center_of_mass",
    "compute_galilean_boost",
    "compute_kinetic_energy",
    "compute_newtonian_potential_energy",
    "compute_hamiltonian",
    "compute_cartesian_hamiltonian",
    "compute_invariants",
    "compute_truncation_residual",
    "JacobiTopologyNode",
    "JacobiTree",
    "build_canonical_tree",
    "enumerate_all_jacobi_trees",
    "JacobiCoordinates",
    "cartesian_to_jacobi",
    "jacobi_to_cartesian",
    "compute_jacobi_kinetic_energy",
    "darboux_pairing_error",
    "rees_deformation_retraction",
    "SL2TransitionMatrix",
    "build_symplectic_matrix_12",
    "build_tree_transition_symplectic_12",
    "rotate_jacobi_tree_coordinates",
    "verify_symplectic_matrix_12",
]
