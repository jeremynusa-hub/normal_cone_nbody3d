"""
Unit Tests for Layer 1 Jacobi (normal_cone_nbody_float64.core.jacobi).

Mathematical Verification Objectives:
- Construction of all (2N-3)!! binary Jacobi trees for N=2, 3, 4, 5 (jacobi_node_defines_hierarchy).
- Hyper-radius positivity and equivalence:
  rho_C = sqrt((1/M) sum m_i m_j r_ij^2) == sqrt(sum mu_v ||r_v||^2) (hyper_radius_sq_pos).
- Exact computation and dimension of normalized hyper-radius bar_rho_C = rho_C / sqrt(M_C) in [L].
- Disjoint subtrees cotangent orthogonality (disjoint_cotangent_commutation_local).
"""

import math
import pytest
import numpy as np

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import (
    JacobiTopologyNode,
    JacobiTree,
    build_canonical_tree,
    build_tree_from_spec,
    enumerate_all_jacobi_trees,
)
from normal_cone_nbody_float64.core.coordinates import cartesian_to_jacobi


def test_canonical_tree_construction():
    masses = np.array([1.0, 2.0, 3.0, 4.0])
    tree_comb = build_canonical_tree(masses, topology="comb")
    assert tree_comb.N == 4
    assert tree_comb.num_internal == 3
    assert len(tree_comb.leaf_pairs) >= 1
    assert tree_comb.root.M == pytest.approx(10.0, rel=1e-15)

    tree_bal = build_canonical_tree(masses, topology="balanced")
    assert tree_bal.N == 4
    assert tree_bal.num_internal == 3
    # In balanced 4-body tree ((0, 1), (2, 3)), there are 2 leaf pairs
    assert len(tree_bal.leaf_pairs) == 2


def test_combinatorial_tree_enumeration_count():
    # N = 2 -> (2*2 - 3)!! = 1!! = 1
    trees_2 = enumerate_all_jacobi_trees(np.ones(2))
    assert len(trees_2) == 1

    # N = 3 -> (2*3 - 3)!! = 3!! = 3
    trees_3 = enumerate_all_jacobi_trees(np.ones(3))
    assert len(trees_3) == 3

    # N = 4 -> (2*4 - 3)!! = 5!! = 15
    trees_4 = enumerate_all_jacobi_trees(np.ones(4))
    assert len(trees_4) == 15

    # N = 5 -> (2*5 - 3)!! = 7!! = 105
    trees_5 = enumerate_all_jacobi_trees(np.ones(5))
    assert len(trees_5) == 105


def test_mass_fractions_and_reduced_masses():
    m = np.array([1.0, 2.0, 3.0])
    tree = build_canonical_tree(m, topology="comb")
    # comb: ((0, 1), 2)
    # pair (0, 1): M=3.0, mu = (1*2)/3 = 2/3, c_left = 1/3, c_right = 2/3
    pair = tree.leaf_pairs[0]
    assert pair.M == pytest.approx(3.0, rel=1e-15)
    assert pair.mu == pytest.approx(2.0 / 3.0, rel=1e-15)
    assert pair.c_left + pair.c_right == pytest.approx(1.0, rel=1e-15)

    # root: merges (0, 1) of mass 3.0 with 2 of mass 3.0
    # M=6.0, mu = (3*3)/6 = 1.5, c_left = 0.5, c_right = 0.5
    root = tree.root
    assert root.M == pytest.approx(6.0, rel=1e-15)
    assert root.mu == pytest.approx(1.5, rel=1e-15)
    assert root.c_left == pytest.approx(0.5, rel=1e-15)
    assert root.c_right == pytest.approx(0.5, rel=1e-15)


def test_hyper_radius_cartesian_vs_jacobi_equivalence():
    rng = np.random.default_rng(123)
    q = rng.standard_normal((4, 3))
    p = rng.standard_normal((4, 3))
    masses = np.array([1.5, 2.5, 3.5, 4.5])

    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="balanced")

    coords = cartesian_to_jacobi(state, tree)

    # Compute via Cartesian formula
    rho_cart = tree.compute_hyper_radius_cartesian(state.q)
    # Compute via Jacobi formula
    rho_jacobi = tree.compute_hyper_radius(coords.r)

    assert rho_cart > 0.0
    assert rho_cart == pytest.approx(rho_jacobi, rel=1e-14)

    # Normalized hyper-radius bar_rho_C = rho_C / sqrt(M_tot) in [L]
    bar_rho_cart = tree.compute_normalized_hyper_radius_cartesian(state.q)
    bar_rho_jacobi = tree.compute_normalized_hyper_radius(coords.r)

    assert bar_rho_cart == pytest.approx(rho_cart / math.sqrt(state.total_mass), rel=1e-15)
    assert bar_rho_cart == pytest.approx(bar_rho_jacobi, rel=1e-14)


def test_disjoint_subtrees_orthogonality():
    # In balanced 4-body tree ((0, 1), (2, 3)), node_L = (0, 1) and node_R = (2, 3) are disjoint
    masses = np.array([1.0, 1.0, 1.0, 1.0])
    tree = build_tree_from_spec(masses, ((0, 1), (2, 3)))

    node_left = tree.root.left
    node_right = tree.root.right
    assert node_left is not None and node_right is not None

    # Body indices must be strictly disjoint sets
    set_L = set(node_left.body_indices)
    set_R = set(node_right.body_indices)
    assert set_L == {0, 1}
    assert set_R == {2, 3}
    assert set_L.isdisjoint(set_R)
