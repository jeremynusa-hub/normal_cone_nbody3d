"""
Unit Tests for Layer 3: Autonomous Intrinsic Hill Tidal Supernode Detector.

Theoretical Alignments:
- Monograph: Section 4.3 (Intrinsic Compactness and Marchal-Sundman Energy Threshold),
  Proposition 4.4 and Equations (4.181)-(4.227).
- Formal Lean 4: `01_PartitionLattice_Defs.lean` (`collision_diagonal_non_transversality`)
  and `01_UnifiedSuperNodeDef_Theorems.lean`.
- Invariants: Exact spatial scale invariance eta_S(lambda * q) = eta_S(q),
  transitive collision non-transversality coalescence, and critical threshold eta_crit = 3^{-1/3}.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.constants import ETA_CRITICAL_HILL
from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.jacobi import build_canonical_tree
from normal_cone_nbody_float64.potential.supernode_detector import (
    compute_tidal_ratio,
    compute_natural_threshold,
    detect_supernodes,
)


def test_tidal_ratio_scale_invariance():
    """
    Validates that the dimensionless Hill tidal ratio eta_S is strictly invariant
    under spatial homotheties q |-> lambda * q:
        eta_S(lambda * q) == eta_S(q).
    """
    bar_rho_S = 0.05
    r_P = 2.0
    M_S = 3.0
    M_ext = 1.0

    eta_base = compute_tidal_ratio(bar_rho_S, r_P, M_S, M_ext)

    for scale in [1e-4, 0.01, 0.5, 2.0, 100.0, 1e5]:
        eta_scaled = compute_tidal_ratio(scale * bar_rho_S, scale * r_P, M_S, M_ext)
        assert np.isclose(eta_scaled, eta_base, atol=1e-15), (
            f"Scale invariance failed at scale={scale}: base={eta_base}, scaled={eta_scaled}"
        )


def test_natural_threshold_scaling():
    """Validates the Marchal-Sundman natural clustering length formula."""
    G = 1.0
    M_tot = 10.0
    M_S = 2.0
    E0 = -5.0

    delta_nat = compute_natural_threshold(G, M_tot, M_S, E0, eta_crit=ETA_CRITICAL_HILL)
    expected = ETA_CRITICAL_HILL * (G * (10.0 ** (5.0 / 3.0)) * (2.0 ** (5.0 / 6.0))) / (2.0 * 5.0)

    assert np.isclose(delta_nat, expected, atol=1e-15)
    assert delta_nat > 0.0


def test_supernode_detection_compact_hierarchy():
    """
    Validates detection of a tight 3-body cluster in a 4-body system
    where eta_S <= 3^{-1/3}.
    """
    masses = np.array([1.0, 1.0, 1.0, 10.0], dtype=np.float64)
    # Bodies 0, 1, 2 form a tight equilateral cluster of size ~ 0.01
    # Body 3 is far away at distance ~ 10.0
    q = np.array([
        [0.01, 0.0, 0.0],
        [-0.005, 0.00866, 0.0],
        [-0.005, -0.00866, 0.0],
        [10.0, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    supernodes = detect_supernodes(state, tree)
    assert len(supernodes) >= 1
    # Cluster {0, 1, 2} must be detected as a supernode
    assert (0, 1, 2) in supernodes


def test_transitive_non_transversality_coalescence():
    """
    Validates Lean 4 theorem 'collision_diagonal_non_transversality':
    Transitive overlap of tight binary encounters merges into a unified supernode.
    """
    masses = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    # Triple {0, 1, 2} tightly coupled, body 3 far away
    q = np.array([
        [0.02, 0.0, 0.0],
        [0.0, 0.02, 0.0],
        [-0.02, 0.0, 0.0],
        [50.0, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    clusters = detect_supernodes(state, tree)
    # Must yield a merged supernode containing {0, 1, 2}
    found_triple = any(set([0, 1, 2]).issubset(set(c)) for c in clusters)
    assert found_triple, f"Transitive coalescence failed, found clusters: {clusters}"


def test_supernode_detection_simultaneous_triple_collision_n3():
    """
    Validates detection of a simultaneous compact triple collision for N=3
    without hierarchical binary dominance (all mutual distances ~ 0.1 <= 0.4).
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    # Equilateral triangle of size 0.1
    q = np.array([
        [0.05, 0.0, 0.0],
        [-0.025, 0.0433, 0.0],
        [-0.025, -0.0433, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    clusters = detect_supernodes(state, tree)
    assert len(clusters) == 1
    assert clusters[0] == (0, 1, 2)


def test_supernode_detection_tetrahedral_symmetry_n4():
    """
    Validates that a compact regular tetrahedral 4-body configuration is
    symmetrically qualified as a unified 4-body supernode (0, 1, 2, 3).
    """
    masses = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    scale = 0.05
    q = scale * np.array([
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ], dtype=np.float64) / np.sqrt(3)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    clusters = detect_supernodes(state, tree)
    assert len(clusters) == 1
    assert clusters[0] == (0, 1, 2, 3), f"Expected unified 4-body supernode, got {clusters}"


def test_schmitt_trigger_hysteresis_anti_chattering():
    """
    Validates that the canonical spectral Schmitt trigger band:
        [eta_enter, eta_exit] = [2^(-1/3), 2^(-1/6)] ~= [0.7937, 0.8909]
    eliminates topological chattering (flickering) strictly above eta_crit = 3^(-1/3) ~= 0.6934:
    1. Inactive cluster with eta in (2^(-1/3), 2^(-1/6)) stays INACTIVE (requires eta <= 2^(-1/3) to enter).
    2. Enters supernode status when eta <= 2^(-1/3) (radial equigravity anticipation).
    3. Remains ACTIVE when eta rises to 0.84 in (2^(-1/3), 2^(-1/6)) (requires eta > 2^(-1/6) to exit).
    4. Exits supernode status only when eta > 2^(-1/6) (spherical RMS disruption).
    """
    masses = np.array([1.0, 1.0, 1.0, 8.0], dtype=np.float64)
    # Bodies 0, 1, 2 form cluster of size ~ 0.5
    # M_S = 3.0, bar_rho_S ~ 0.5 / sqrt(3) ~ 0.288675
    # M_ext = 8.0 -> (M_ext / M_S)^(1/3) = (8/3)^(1/3) ~ 1.3867
    # eta = (bar_rho_S / r_P) * 1.3867
    scale = 10.0
    q_base = scale * np.array([
        [0.05, 0.0, 0.0],
        [-0.025, 0.04330127, 0.0],
        [-0.025, -0.04330127, 0.0],
        [1.0, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q_base)
    tree = build_canonical_tree(masses, topology="comb")

    # Measure exact bar_rho_S and factor
    m_sub = masses[:3]
    cdm_S = np.sum(m_sub[:, None] * q_base[:3], axis=0) / 3.0
    rho_S = np.sqrt(np.sum(m_sub * np.sum((q_base[:3] - cdm_S)**2, axis=1)))
    bar_rho_S = rho_S / np.sqrt(3.0)
    mass_factor = (8.0 / 3.0) ** (1.0 / 3.0)

    # 1. Choose r_P such that eta = 0.84 in (2^(-1/3) ~= 0.7937, 2^(-1/6) ~= 0.8909)
    r_P_mid = (bar_rho_S * mass_factor) / 0.84
    q1 = np.copy(q_base)
    q1[3] = [r_P_mid, 0.0, 0.0]
    state1 = ExtendedPhaseState(q=q1, p=p, masses=masses)

    # When initially inactive: must NOT trigger supernode (since 0.84 > eta_enter ~= 0.7937)
    sn_inactive = detect_supernodes(state1, tree, active_supernodes=[])
    assert len(sn_inactive) == 0, f"Expected no supernode, got {sn_inactive}"

    # 2. Choose r_P such that eta = 0.75 <= 2^(-1/3) ~= 0.7937
    r_P_close = (bar_rho_S * mass_factor) / 0.75
    q2 = np.copy(q_base)
    q2[3] = [r_P_close, 0.0, 0.0]
    state2 = ExtendedPhaseState(q=q2, p=p, masses=masses)

    sn_enter = detect_supernodes(state2, tree, active_supernodes=[])
    assert len(sn_enter) == 1 and (0, 1, 2) in sn_enter, f"Expected supernode to activate, got {sn_enter}"

    # 3. Now return to r_P_mid (eta = 0.84), but with active_supernodes=[(0, 1, 2)]
    # Under Schmitt trigger hysteresis, it must REMAIN ACTIVE because eta = 0.84 <= eta_exit ~= 0.8909!
    sn_stay = detect_supernodes(state1, tree, active_supernodes=[(0, 1, 2)])
    assert len(sn_stay) == 1 and (0, 1, 2) in sn_stay, (
        f"Hysteresis failed: supernode should have stayed active at eta=0.84 <= 0.8909, got {sn_stay}"
    )

    # 4. Choose r_P such that eta = 0.95 > 2^(-1/6) ~= 0.8909
    r_P_far = (bar_rho_S * mass_factor) / 0.95
    q3 = np.copy(q_base)
    q3[3] = [r_P_far, 0.0, 0.0]
    state3 = ExtendedPhaseState(q=q3, p=p, masses=masses)

    # Even with active_supernodes=[(0, 1, 2)], it must DEACTIVATE because eta = 0.95 > eta_exit ~= 0.8909!
    sn_exit = detect_supernodes(state3, tree, active_supernodes=[(0, 1, 2)])
    assert len(sn_exit) == 0, f"Expected supernode to exit at eta=0.95 > 0.8909, got {sn_exit}"


def test_internal_hierarchy_filter_rejects_hierarchical_subsystems():
    """
    Validates that a subsystem containing a tight binary sub-pair and a distant spectator
    (chi_S = r_max^{int} / r_min^{int} >= 3.0) is rejected by detect_supernodes,
    preserving analytical hierarchical Jacobi tree resolution (Definition 6.1 & Theorem 6.1).
    """
    from normal_cone_nbody_float64.potential.supernode_detector import (
        compute_internal_hierarchy_ratio,
    )

    masses = np.array([1.0, 1.0, 1.0, 10.0], dtype=np.float64)
    # Bodies 0 and 1 are a very tight binary (d = 0.001)
    # Body 2 is in hierarchical orbit around {0, 1} at d = 0.05
    # Body 3 is external perturber at d = 5.0
    q = np.array([
        [0.0, 0.0, 0.0],
        [0.001, 0.0, 0.0],
        [0.05, 0.0, 0.0],
        [5.0, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    chi = compute_internal_hierarchy_ratio(q, (0, 1, 2))
    assert chi >= 3.0, f"Expected chi_S >= 3.0 for hierarchical triple, got {chi}"

    # Must be filtered out because it is internally hierarchical
    supernodes = detect_supernodes(state, tree, max_hierarchy_ratio=3.0, enforce_non_hierarchical=True)
    assert (0, 1, 2) not in supernodes, (
        f"Hierarchical triple (chi_S={chi:.2f} >= 3.0) should not be flagged as supernode, got {supernodes}"
    )

    # When non-hierarchy enforcement is disabled, it would qualify under pure tidal Hill ratio
    supernodes_unfiltered = detect_supernodes(state, tree, enforce_non_hierarchical=False)
    assert (0, 1, 2) in supernodes_unfiltered


def test_dominant_tidal_perturber_selection_over_closer_dust_particle():
    """
    Validates D-144-07 (Part 1):
    When a candidate cluster S = {0, 1, 2} is surrounded by multiple external bodies,
    eta_S must be governed by the external perturber that maximizes the tidal ratio
    eta_S(k) = (bar_rho_S / r_k) * (m_k / M_S)**(1/3), NOT the one that merely minimizes
    Euclidean distance r_k (which fails when a nearly massless dust particle is closer than
    a massive star).
    """
    # Bodies 0, 1, 2: mass 1.0 each (M_S = 3.0), equilateral triangle with bar_rho_S = 0.5
    # Body 3: tiny dust particle (m_3 = 1e-12) at distance r_3 = 3.0 along +z
    # Body 4: massive star (m_4 = 3000.0) at distance r_4 = 4.0 along -z
    # Tidal ratio from body 3: (0.5 / 3.0) * (1e-12 / 3.0)**(1/3) ~= 1.15e-5 << eta_crit
    # Tidal ratio from body 4: (0.5 / 4.0) * (3000.0 / 3.0)**(1/3) = 0.125 * 10 = 1.25 > eta_exit
    masses = np.array([1.0, 1.0, 1.0, 1e-12, 3000.0], dtype=np.float64)
    q_disrupted = np.array([
        [0.5, 0.0, 0.0],
        [-0.25, 0.4330127018922193, 0.0],
        [-0.25, -0.4330127018922193, 0.0],
        [0.0, 0.0, 3.0],
        [0.0, 0.0, -4.0],
    ], dtype=np.float64)
    p = np.zeros_like(q_disrupted)
    state_disrupted = ExtendedPhaseState(q=q_disrupted, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    # Because the massive star at r_4 = 4.0 exerts eta_S(4) = 1.25 > eta_crit,
    # {0, 1, 2} is tidally disrupted and must NOT be classified as a self-bound Hill supernode
    sn_disrupted = detect_supernodes(state_disrupted, tree)
    assert (0, 1, 2) not in sn_disrupted, (
        f"Cluster {{0,1,2}} is tidally disrupted by massive body 4 (eta=1.25 > eta_crit), "
        f"but closer dust particle 3 masked the tidal shear: got {sn_disrupted}"
    )

    # When massive star 4 is moved to r_4 = 10.0, eta_S(4) = (0.5 / 10.0) * 10 = 0.50 <= eta_crit,
    # so {0, 1, 2} becomes tidally isolated and MUST be detected as a supernode
    q_isolated = np.copy(q_disrupted)
    q_isolated[4] = [0.0, 0.0, -10.0]
    state_isolated = ExtendedPhaseState(q=q_isolated, p=p, masses=masses)
    sn_isolated = detect_supernodes(state_isolated, tree)
    assert (0, 1, 2) in sn_isolated, f"Expected {{0,1,2}} supernode when massive perturber is at r=10, got {sn_isolated}"


def test_supernode_detection_tetrahedral_unequal_masses_n4():
    """
    Validates D-144-07 (Parts 2 & 3):
    1. Transitive non-transversality coalescence on a compact regular tetrahedron (N=4)
       is purely geometric on Q[N] and must unify all 4 bodies into (0, 1, 2, 3)
       even when all 4 masses are distinct (e.g. [1.0, 2.0, 3.0, 4.0]).
    2. An explicitly unbound hyperbolic system (E0 > 0) or two well-separated binaries
       must not be promoted to a global 4-body supernode by the fallback.
    """
    masses = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    scale = 0.05
    q = scale * np.array([
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ], dtype=np.float64) / np.sqrt(3.0)
    p = np.zeros_like(q)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    tree = build_canonical_tree(masses, topology="comb")

    clusters = detect_supernodes(state, tree)
    assert len(clusters) == 1
    assert clusters[0] == (0, 1, 2, 3), (
        f"Expected unified 4-body supernode (0, 1, 2, 3) for unequal masses, got {clusters}"
    )

    # Two well-separated tight binaries: (0, 1) at d=0.01 and (2, 3) at d=0.01 separated by D=10.0
    q_two_binaries = np.array([
        [0.0, 0.0, 0.0],
        [0.01, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [10.01, 0.0, 0.0],
    ], dtype=np.float64)
    state_bins = ExtendedPhaseState(q=q_two_binaries, p=p, masses=masses)
    tree_bins = build_canonical_tree(masses, topology="balanced")
    assert detect_supernodes(state_bins, tree_bins) == []


def test_schmitt_trigger_not_short_circuited_by_negative_E0():
    """
    Validates [ENG-FALSE-07] (Part 1):
    When state.E0 < 0 is set on a bound system with external perturbers, the global
    Marchal-Saari bound delta_nat must NOT short-circuit the instantaneous Hill tidal
    ratio eta_S and the spectral Schmitt trigger band [2^(-1/3), 2^(-1/6)].
    """
    masses = np.array([1.0, 1.0, 1.0, 8.0], dtype=np.float64)
    scale = 10.0
    q_base = scale * np.array([
        [0.05, 0.0, 0.0],
        [-0.025, 0.04330127, 0.0],
        [-0.025, -0.04330127, 0.0],
        [1.0, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q_base)
    tree = build_canonical_tree(masses, topology="comb")

    m_sub = masses[:3]
    cdm_S = np.sum(m_sub[:, None] * q_base[:3], axis=0) / 3.0
    rho_S = np.sqrt(np.sum(m_sub * np.sum((q_base[:3] - cdm_S) ** 2, axis=1)))
    bar_rho_S = rho_S / np.sqrt(3.0)
    mass_factor = (8.0 / 3.0) ** (1.0 / 3.0)

    # Set r_P so eta = 0.84 in (eta_enter = 0.7937, eta_exit = 0.8909) and E0 = -5.0
    r_P_mid = (bar_rho_S * mass_factor) / 0.84
    q1 = np.copy(q_base)
    q1[3] = [r_P_mid, 0.0, 0.0]
    state1 = ExtendedPhaseState(q=q1, p=p, masses=masses, E0=-5.0)

    # Must remain INACTIVE when initially inactive, even with E0 = -5.0!
    assert detect_supernodes(state1, tree, active_supernodes=[]) == []
    # Must remain ACTIVE when initially active!
    assert detect_supernodes(state1, tree, active_supernodes=[(0, 1, 2)]) == [(0, 1, 2)]


def test_roche_tidal_disruption_mode_shoemaker_levy_9():
    """
    Validates [ENG-FALSE-07] (Part 2, §7.2.1 Eq. 7.34, Benchmark H2 Shoemaker-Levy 9):
    1. compute_cubic_tidal_force_ratio evaluates Lambda_S = (1/3)*eta_S^3 = (M_ext/(3*M_S))*(bar_rho_S/r_P)^3.
    2. For a tenuous 3-fragment comet cluster (m_i = 1e-12, spread ~ 1.65e-5 AU) approaching Jupiter
       (M_J = 9.547919e-4) with negative total energy E0 < 0:
       - At initial distance r_0 = 0.012 AU (outside Roche limit, Lambda_S ~= 0.035 < 0.25):
         detect_supernodes(..., mode="roche_disruption") returns [] (no premature activation at t=0).
       - At perijove r_min = 6.33e-4 AU = 1.32 R_J (deep inside Roche limit 2.45 R_J, Lambda_S >= 0.25):
         detect_supernodes(..., mode="roche_disruption") returns [(1, 2, 3)].
    """
    from normal_cone_nbody_float64.potential.supernode_detector import (
        compute_cubic_tidal_force_ratio,
        LAMBDA_CRITICAL_ROCHE,
    )

    assert np.isclose(LAMBDA_CRITICAL_ROCHE, 0.25, atol=1e-15)

    m_jup = 9.547919e-4
    m_frag = 1e-12
    masses = np.array([m_jup, m_frag, m_frag, m_frag], dtype=np.float64)

    # Compact equilateral fragment cluster of radius d_frag = 8.165e-6 AU (bar_rho_S of H2)
    d_frag = 8.165e-6
    offsets = np.array([
        [d_frag, 0.0, 0.0],
        [-0.5 * d_frag, 0.5 * np.sqrt(3.0) * d_frag, 0.0],
        [-0.5 * d_frag, -0.5 * np.sqrt(3.0) * d_frag, 0.0],
    ], dtype=np.float64)

    # 1. At initial distance r_0 = 0.012 AU (outside Roche limit)
    q_init = np.zeros((4, 3), dtype=np.float64)
    q_init[1:] = np.array([0.012, 0.0, 0.0]) + offsets
    state_init = ExtendedPhaseState(q=q_init, p=np.zeros_like(q_init), masses=masses, E0=-1e-10)
    tree_init = build_canonical_tree(masses, topology="comb")
    # Build tree grouping comet fragments (1, 2, 3) vs Jupiter (0)
    from normal_cone_nbody_float64.system.atlas_manager import AtlasManager
    tree_sl9 = AtlasManager().build_optimal_tree(state_init)

    lam_init = compute_cubic_tidal_force_ratio(d_frag, 0.012, 3.0 * m_frag, m_jup)
    assert lam_init < LAMBDA_CRITICAL_ROCHE, f"Expected Lambda_S(0) < 0.25, got {lam_init}"
    assert detect_supernodes(state_init, tree_sl9) == []
    assert detect_supernodes(state_init, tree_sl9, mode="roche_disruption") == []
    assert detect_supernodes(state_init, tree_sl9, mode="hill_bound") == [(1, 2, 3)]

    # 2. At perijove r_min = 6.33e-4 AU (inside Roche limit, Lambda_S >> 0.25)
    q_peri = np.zeros((4, 3), dtype=np.float64)
    q_peri[1:] = np.array([6.33e-4, 0.0, 0.0]) + offsets
    state_peri = ExtendedPhaseState(q=q_peri, p=np.zeros_like(q_peri), masses=masses, E0=-1e-10)
    tree_peri = AtlasManager().build_optimal_tree(state_peri)

    lam_peri = compute_cubic_tidal_force_ratio(d_frag, 6.33e-4, 3.0 * m_frag, m_jup)
    assert lam_peri > 0.25, f"Expected Lambda_S(peri) > 0.25, got {lam_peri}"
    assert detect_supernodes(state_peri, tree_peri) == [(1, 2, 3)]
    assert detect_supernodes(state_peri, tree_peri, mode="roche_disruption") == [(1, 2, 3)]


def test_canonical_rule1_rejects_figure8_barycenter_crossing_and_burrau():
    """
    Validates Canonical Rule 1 (§4.3.2 & §6.4.2):
    Even when a body sits at the barycenter q_i = 0 with L_tot = 0 (as in Chenciner-Montgomery's
    Figure-8 choreography at t=0 where q_2 = (0,0,0), v_2 != 0, and x_0 x v_0 != 0),
    or in Burrau's 3-4-5 non-central right triangle at rest, detect_supernodes(..., mode="canonical")
    returns [] (no false global 3-body supernode promotion).
    """
    from normal_cone_nbody_float64.benchmarks.registry import (
        create_figure8_choreography_benchmark,
        create_burrau_pythagorean_benchmark,
    )
    from normal_cone_nbody_float64.system.atlas_manager import AtlasManager

    st_fig8 = create_figure8_choreography_benchmark()
    tree_fig8 = AtlasManager().build_optimal_tree(st_fig8)
    assert detect_supernodes(st_fig8, tree_fig8) == []

    st_burrau = create_burrau_pythagorean_benchmark()
    tree_burrau = AtlasManager().build_optimal_tree(st_burrau)
    assert detect_supernodes(st_burrau, tree_burrau) == []


