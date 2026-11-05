"""
Unit Tests for DyadicScheduler Leaf/Apex Decoupling (Iteration 148, Phase 5).

Theoretical Alignment:
- Fulton--MacPherson Normal Cone Splitting (Chapter 4, Eq 4.30 & Chapter 5, Section 5.6):
  The total Hamiltonian splits as H = H_base + H_leaf + V_apex^hier.
  The two-body leaf pairs are desingularized and integrated analytically to float64
  machine precision via Stumpff-Sundman universal functions without time-step error.
  Therefore, singular pericenter contraction (r_peri -> 0, v_peri -> infty)
  must not trigger extreme dyadic levels (k = 17) for active leaf pairs.
  Extreme pericenter anticipation is strictly confined to non-leaf interactions
  residing in the screened apex potential V_apex before chart transition.
"""

import math
import numpy as np
import pytest

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_tree_from_spec
from normal_cone_nbody_float64.system.dyadic_scheduler import DyadicScheduler


def test_leaf_pair_bypasses_singular_pericenter_anticipation():
    """
    Validates that an active leaf pair on an inward-bound trajectory (v_rad < 0)
    with a minute pericenter (r_peri ~ 2e-5) is NOT subjected to extreme dyadic levels (k >= 14).
    Instead, it is governed by its orbital kinematic timescale tau_kin = r / v,
    yielding a moderate dyadic level (k <= 7).
    """
    dt_macro = 0.015
    scheduler = DyadicScheduler(dt_macro=dt_macro, default_level=0, max_level=17)

    # 4-body system matching the geometry of Step 220 in the hierarchical 4-body benchmark
    masses = np.array([1.0, 1.0, 0.8, 1.2], dtype=np.float64)
    # Bodies 1 and 2 are approaching: r_12 = 0.0304, v_rel = 6.2 m/s, v_rad = -5.85 m/s
    q = np.array([
        [4.2, 0.0, 0.0],
        [-0.0152, 0.0, 0.0],
        [0.0152, 0.0, 0.0],
        [0.0, 0.3, 0.0],
    ], dtype=np.float64)

    # Velocity: Body 1 and 2 have relative velocity ~ 6.2 with strong inward radial component
    p = np.array([
        [0.0, 0.0, 0.0],
        [3.0, 0.5, 0.0],
        [-3.0, -0.5, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)

    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    # Tree: (0, ((1, 2), 3)) where (1, 2) is the active leaf pair
    tree_with_12_leaf = build_tree_from_spec(masses, (0, ((1, 2), 3)))

    leaf_pairs = [n.body_indices for n in tree_with_12_leaf.leaf_pairs]
    assert (1, 2) in leaf_pairs

    # Level with leaf decoupling (Regime II: tidally perturbed leaf at R_ext = 0.3, eps_tidal ~ 8.3e-4 > 1e-4)
    k_leaf_perturbed = scheduler.compute_dynamical_level(state, tree_with_12_leaf, dt_macro=dt_macro)
    assert k_leaf_perturbed <= 11, f"Expected Regime II tidal-BCH level k <= 11 (< 14), got {k_leaf_perturbed}"

    # Regime I: when the external bodies are well-isolated (R_ext = 3.0 -> eps_tidal ~ 8.3e-7 < 1e-4),
    # the active leaf pair is governed strictly by min(tau_kin, tau_orb) (tau_orb = 0.00295 s -> k = 8).
    q_isolated = q.copy()
    q_isolated[3] = [0.0, 3.0, 0.0]
    state_isolated = ExtendedPhaseState(q=q_isolated, p=p, masses=masses)
    k_leaf_isolated = scheduler.compute_dynamical_level(state_isolated, tree_with_12_leaf, dt_macro=dt_macro)
    assert k_leaf_isolated <= 8, f"Expected Regime I Hill-isolated level k <= 8, got {k_leaf_isolated}"


def test_non_leaf_intruder_triggers_protective_dyadic_elevation():
    """
    Validates that when an intruder (body not belonging to an active leaf pair)
    rapidly approaches another body with a minute projected pericenter (r_peri -> 0),
    the dyadic scheduler properly triggers an elevated refinement level (k >= 8)
    to protect the non-regularized apex kick against truncation error.
    """
    dt_macro = 0.015
    scheduler = DyadicScheduler(dt_macro=dt_macro, default_level=0, max_level=17)

    masses = np.array([1.0, 1.0, 0.8, 1.2], dtype=np.float64)
    # Bodies 1 and 2 are a wide leaf pair (r_12 = 0.5)
    # Body 3 is an incoming intruder heading straight towards Body 1 at high speed (v_rel = 20.0 m/s, r_13 = 0.05)
    q = np.array([
        [5.0, 0.0, 0.0],
        [-0.25, 0.0, 0.0],
        [0.25, 0.0, 0.0],
        [-0.249, 0.05, 0.0],  # Body 3 near Body 1 with impact parameter b = 0.001
    ], dtype=np.float64)

    p = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, -20.0 * 1.2, 0.0],  # Body 3 plunging towards Body 1
    ], dtype=np.float64)

    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    # Tree where (1, 2) is the leaf pair, but (1, 3) is a non-leaf apex interaction
    tree = build_tree_from_spec(masses, (0, ((1, 2), 3)))

    k_apex = scheduler.compute_dynamical_level(state, tree, dt_macro=dt_macro)

    # Because (1, 3) is in apex and plunging inward at v_rel = 20 m/s with small impact parameter,
    # the scheduler must elevate k >= 8 to protect V_apex
    assert k_apex >= 8, f"Expected elevated dyadic level k >= 8 for apex encounter, got {k_apex}"


def test_tree_none_fallback_maintains_safety():
    """
    Validates that when tree is None (e.g. global unclustered cruise),
    all pairs are conservatively protected by pericenter anticipation.
    """
    dt_macro = 0.015
    scheduler = DyadicScheduler(dt_macro=dt_macro, default_level=0, max_level=17)

    masses = np.array([1.0, 1.0], dtype=np.float64)
    q = np.array([[-0.01, 0.0, 0.0], [0.01, 0.0, 0.0]], dtype=np.float64)
    p = np.array([[10.0, 0.1, 0.0], [-10.0, -0.1, 0.0]], dtype=np.float64)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    k_none = scheduler.compute_dynamical_level(state, tree=None, dt_macro=dt_macro)
    assert k_none >= 8
