r"""
Comprehensive Unit Test Suite for High-Precision Float128 / Mpmath Core Subsystem.

Validates:
1. Multi-precision constants (mp.dps = 35, EPS_MACH_F128 < 2e-34, Blanes--Moan weights, homological constants).
2. ExtendedPhaseState immutability, bidirectionality (f64 <-> f128), and barycentering (< 1e-32).
3. Lie--Bargmann invariants: exact Burrau Hamiltonian (-769/60) to 35 digits, total momentum/angular momentum/boost.
4. Jacobi binary tree topology and coordinates: exact involutivity J^-1 o J = Id (< 1e-32) and Darboux pairing.
5. SL(2, R) symplectic chart transitions: |det(M) - 1| < 1e-32, Sp(12, R) certification, and Darboux 2-form conservation.
"""

import pytest
import numpy as np
import mpmath as mp

from normal_cone_nbody_float128.core.constants import (
    DEFAULT_DPS,
    EPS_MACH_F128,
    G_DEFAULT_F128,
    ETA_CRIT_F128,
    ETA_ENTER_HILL_F128,
    ETA_EXIT_HILL_F128,
    GAMMA_HYSTERESIS_F128,
    DELTA_LAMBDA_HYSTERESIS_F128,
    BM4_DRIFT_WEIGHTS_F128,
    BM4_IMPULSE_WEIGHTS_F128,
    BM6_DRIFT_WEIGHTS_F128,
    BM6_IMPULSE_WEIGHTS_F128,
    KAPPA4_BM4_F128,
    to_f128,
    to_f64,
)
from normal_cone_nbody_float128.core.state import ExtendedPhaseState
from normal_cone_nbody_float128.core.invariants import (
    compute_invariants,
    compute_kinetic_energy,
    compute_newtonian_potential_energy,
    compute_hamiltonian,
    compute_truncation_residual,
)
from normal_cone_nbody_float128.core.jacobi import (
    build_canonical_tree,
    build_tree_from_spec,
    JacobiTree,
)
from normal_cone_nbody_float128.core.coordinates import (
    cartesian_to_jacobi,
    jacobi_to_cartesian,
    compute_jacobi_kinetic_energy,
    darboux_pairing_error,
)
from normal_cone_nbody_float128.core.transitions import (
    SL2TransitionMatrix,
    build_symplectic_matrix_12,
    verify_symplectic_matrix_12,
)


def test_constants_precision_and_tolerances():
    """Validates that constants maintain 35 decimal places and exact mathematical relations."""
    assert mp.mp.dps >= 35
    assert EPS_MACH_F128 < mp.mpf("2e-34")

    # eta_crit = 3^(-1/3)
    expected_eta_crit = mp.power(mp.mpf("3.0"), -mp.mpf("1.0") / mp.mpf("3.0"))
    assert abs(ETA_CRIT_F128 - expected_eta_crit) < mp.mpf("1e-34")

    # eta_enter = 2^(-1/3), eta_exit = 2^(-1/6)
    expected_enter = mp.power(mp.mpf("2.0"), -mp.mpf("1.0") / mp.mpf("3.0"))
    expected_exit = mp.power(mp.mpf("2.0"), -mp.mpf("1.0") / mp.mpf("6.0"))
    assert abs(ETA_ENTER_HILL_F128 - expected_enter) < mp.mpf("1e-34")
    assert abs(ETA_EXIT_HILL_F128 - expected_exit) < mp.mpf("1e-34")

    # eta_enter == (eta_exit)^2
    assert abs(ETA_ENTER_HILL_F128 - ETA_EXIT_HILL_F128 ** 2) < mp.mpf("1e-34")

    # Sum of BM4 drift weights = 1.0 and cubic order condition 2*(a1^3+a2^3+a3^3)+a4^3 = 0
    sum_bm4_drift = sum(BM4_DRIFT_WEIGHTS_F128)
    assert abs(sum_bm4_drift - mp.mpf("1.0")) < mp.mpf("1e-34")
    cubic_bm4 = sum(w ** 3 for w in BM4_DRIFT_WEIGHTS_F128)
    assert abs(cubic_bm4) < mp.mpf("1e-34")

    # Sum of BM4 impulse weights = 1.0
    sum_bm4_impulse = sum(BM4_IMPULSE_WEIGHTS_F128)
    assert abs(sum_bm4_impulse - mp.mpf("1.0")) < mp.mpf("1e-34")

    # BM6 palindromic sum: 2 * (a1 + a2 + a3 + a4 + a5) + a6 = 1.0 to 35 digits
    sum_bm6_drift = 2 * sum(BM6_DRIFT_WEIGHTS_F128[:5]) + BM6_DRIFT_WEIGHTS_F128[5]
    assert abs(sum_bm6_drift - mp.mpf("1.0")) < mp.mpf("1e-34")
    sum_bm6_impulse = 2 * sum(BM6_IMPULSE_WEIGHTS_F128[:5]) + BM6_IMPULSE_WEIGHTS_F128[5]
    assert abs(sum_bm6_impulse - mp.mpf("1.0")) < mp.mpf("1e-34")

    # BM6 RKN order-3 quadrature condition: sum_k B_k c_k^2 == 1/3 to 35 digits
    A_seq = list(BM6_DRIFT_WEIGHTS_F128[:5]) + [BM6_DRIFT_WEIGHTS_F128[5]] + list(reversed(BM6_DRIFT_WEIGHTS_F128[:5]))
    B_seq = (
        list(BM6_IMPULSE_WEIGHTS_F128[:5])
        + [mp.mpf("0.5") * BM6_IMPULSE_WEIGHTS_F128[5], mp.mpf("0.5") * BM6_IMPULSE_WEIGHTS_F128[5]]
        + list(reversed(BM6_IMPULSE_WEIGHTS_F128[:5]))
    )
    c_nodes = [mp.mpf("0.0")]
    for ak in A_seq:
        c_nodes.append(c_nodes[-1] + ak)
    rkn_order3_res = sum(B_seq[k] * (c_nodes[k] ** 2) for k in range(12)) - mp.mpf("1.0") / mp.mpf("3.0")
    assert abs(rkn_order3_res) < mp.mpf("1e-34")

    # Yoshida 4th-order weights in 35-digit precision: 2*w1 + w0 = 1 and 2*w1^3 + w0^3 = 0
    from normal_cone_nbody_float128.integrators.yoshida import (
        YOSHIDA_DRIFT_WEIGHTS_F128,
    )
    assert abs(sum(YOSHIDA_DRIFT_WEIGHTS_F128) - mp.mpf("1.0")) < mp.mpf("1e-34")
    assert abs(sum(w ** 3 for w in YOSHIDA_DRIFT_WEIGHTS_F128)) < mp.mpf("1e-34")

    # Homological coefficients (Theorem 5.14 compliance: no KAPPA6_* in core.constants)
    import normal_cone_nbody_float128.core.constants as const_f128_mod

    assert abs(KAPPA4_BM4_F128 - mp.mpf("1.0") / mp.mpf("2880.0")) < mp.mpf("1e-34")
    assert not hasattr(const_f128_mod, "KAPPA6_BM6")
    assert not hasattr(const_f128_mod, "KAPPA6_BM6_F128")
    assert not hasattr(const_f128_mod, "KAPPA6_CONORMAL_LEAF")
    assert not hasattr(const_f128_mod, "KAPPA6_CONORMAL_LEAF_F128")


def test_extended_phase_state_f128():
    """Validates ExtendedPhaseState conversion, immutability, and barycentering in f128."""
    q_f64 = np.array([[1.0, 2.0, 0.0], [-2.0, -1.0, 0.0], [3.0, -4.0, 0.0]], dtype=np.float64)
    p_f64 = np.array([[0.1, 0.2, 0.0], [-0.3, 0.1, 0.0], [0.2, -0.3, 0.0]], dtype=np.float64)
    m_f64 = np.array([3.0, 4.0, 5.0], dtype=np.float64)

    s_f64 = ExtendedPhaseState(q=q_f64, p=p_f64, masses=m_f64, t=0.0)
    assert not s_f64.is_f128

    # Convert to f128
    s_f128 = s_f64.to_f128()
    assert s_f128.is_f128
    assert isinstance(s_f128.q[0, 0], mp.mpf)

    # Immutability
    with pytest.raises(ValueError):
        s_f128.q[0, 0] = mp.mpf("99.0")

    # Barycentering
    s_bary = s_f128.to_barycentric()
    r_cm = s_bary.center_of_mass
    p_tot = s_bary.total_momentum
    assert max(abs(x) for x in r_cm) < mp.mpf("1e-32")
    assert max(abs(x) for x in p_tot) < mp.mpf("1e-32")

    # Min pairwise distance and collision predicate
    min_d = s_bary.min_pairwise_distance()
    assert isinstance(min_d, mp.mpf)
    assert min_d > mp.mpf("1.0")
    assert s_bary.is_collision_free()

    # Edge case N = 1 with skip_validation=True in f128
    s_single_f128 = ExtendedPhaseState(
        q=np.array([[mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")]], dtype=object),
        p=np.array([[mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")]], dtype=object),
        masses=np.array([mp.mpf("1.0")], dtype=object),
        skip_validation=True,
    )
    min_d_single = s_single_f128.min_pairwise_distance()
    assert isinstance(min_d_single, mp.mpf)
    assert min_d_single == mp.mpf("inf")

    # Convert back to f64
    s_back = s_bary.to_f64()
    assert not s_back.is_f128
    assert isinstance(s_back.q[0, 0], np.float64)


def test_lie_bargmann_invariants_burrau_f128():
    """Validates exact analytical Hamiltonian H = -769/60 on Burrau Pythagorean triple to 35 decimals."""
    # Burrau initial configuration: right triangle with sides 3, 4, 5
    q = np.array([
        [mp.mpf("1.0"), mp.mpf("3.0"), mp.mpf("0.0")],
        [mp.mpf("-2.0"), mp.mpf("-1.0"), mp.mpf("0.0")],
        [mp.mpf("1.0"), mp.mpf("-1.0"), mp.mpf("0.0")],
    ], dtype=object)
    p = np.array([
        [mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")],
        [mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")],
        [mp.mpf("0.0"), mp.mpf("0.0"), mp.mpf("0.0")],
    ], dtype=object)
    masses = np.array([mp.mpf("3.0"), mp.mpf("4.0"), mp.mpf("5.0")], dtype=object)

    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    inv = compute_invariants(state)

    assert inv.is_f128
    # Exact analytical Hamiltonian: H0 = -769/60 = -12.816666666666666666666666666666667
    expected_H = -mp.mpf("769.0") / mp.mpf("60.0")
    assert abs(inv.H - expected_H) < mp.mpf("1e-34")

    # Invariants: P_tot == 0, L_tot == 0, G_boost == 0 at t=0 for rest state
    assert max(abs(x) for x in inv.P) == mp.mpf("0.0")
    assert max(abs(x) for x in inv.L) == mp.mpf("0.0")

    # Residual between identical states is bitwise zero
    res = compute_truncation_residual(inv, inv)
    assert all(x == mp.mpf("0.0") for x in res)


def test_jacobi_and_coordinates_f128():
    """Validates that Jacobi transformations achieve J^-1 o J = Id with error < 1e-32."""
    masses = np.array([mp.mpf("3.0"), mp.mpf("4.0"), mp.mpf("5.0")], dtype=object)
    tree = build_canonical_tree(masses)
    assert tree.is_f128

    # Verify cluster mass metrics
    assert abs(tree.root.M - mp.mpf("12.0")) < mp.mpf("1e-34")

    q = np.array([
        [mp.mpf("1.5"), mp.mpf("2.5"), mp.mpf("0.5")],
        [mp.mpf("-1.2"), mp.mpf("-0.8"), mp.mpf("0.3")],
        [mp.mpf("0.7"), mp.mpf("-1.4"), mp.mpf("-0.6")],
    ], dtype=object)
    p = np.array([
        [mp.mpf("0.4"), mp.mpf("0.1"), mp.mpf("-0.2")],
        [mp.mpf("-0.3"), mp.mpf("-0.2"), mp.mpf("0.1")],
        [mp.mpf("-0.1"), mp.mpf("0.1"), mp.mpf("0.1")],
    ], dtype=object)

    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    coords = cartesian_to_jacobi(state, tree)

    # Reconstruct Cartesian state
    state_rec = jacobi_to_cartesian(
        coords,
        tree,
        masses,
        center_of_mass=state.center_of_mass,
        total_momentum=state.total_momentum,
    )

    err_q = max(abs(state.q[i, c] - state_rec.q[i, c]) for i in range(3) for c in range(3))
    err_p = max(abs(state.p[i, c] - state_rec.p[i, c]) for i in range(3) for c in range(3))
    assert err_q < mp.mpf("1e-32")
    assert err_p < mp.mpf("1e-32")

    # Darboux pairing preservation
    err_darboux = darboux_pairing_error(state, coords, tree)
    assert err_darboux < mp.mpf("1e-32")

    # Hyper-radius consistency between Jacobi and Cartesian definitions
    rho_jacobi = tree.compute_hyper_radius(coords.r)
    rho_cart = tree.compute_hyper_radius_cartesian(state.q)
    assert abs(rho_jacobi - rho_cart) < mp.mpf("1e-32")

    rho_norm_jacobi = tree.compute_normalized_hyper_radius(coords.r)
    rho_norm_cart = tree.compute_normalized_hyper_radius_cartesian(state.q)
    assert isinstance(rho_norm_jacobi, mp.mpf)
    assert isinstance(rho_norm_cart, mp.mpf)
    assert abs(rho_norm_jacobi - rho_norm_cart) < mp.mpf("1e-32")


def test_sl2_transitions_f128():
    """Validates SL(2, R) unimodularity |det(M) - 1| < 1e-32 and Sp(12, R) symplecticity."""
    mA = mp.mpf("3.0")
    mB = mp.mpf("4.0")
    mC = mp.mpf("5.0")

    trans = SL2TransitionMatrix.from_masses(mA, mB, mC)
    assert trans.is_f128

    # Unimodular condition: det(M) == 1.0 exactly
    assert abs(trans.determinant - mp.mpf("1.0")) < mp.mpf("1e-34")

    # Sp(12, R) matrix certification
    Psi = build_symplectic_matrix_12(trans)
    assert verify_symplectic_matrix_12(Psi, tol=mp.mpf("1e-32"))

    # Darboux 2-form invariance under canonical position/momentum action
    r1 = np.array([mp.mpf("1.2"), mp.mpf("0.8"), mp.mpf("-0.4")], dtype=object)
    r2 = np.array([mp.mpf("-0.5"), mp.mpf("1.7"), mp.mpf("0.3")], dtype=object)
    w1 = np.array([mp.mpf("0.3"), mp.mpf("-0.2"), mp.mpf("0.5")], dtype=object)
    w2 = np.array([mp.mpf("-0.1"), mp.mpf("0.4"), mp.mpf("-0.2")], dtype=object)

    r1_p, r2_p = trans.transform_positions(r1, r2)
    w1_p, w2_p = trans.transform_momenta(w1, w2)

    pairing_before = np.sum(r1 * w1) + np.sum(r2 * w2)
    pairing_after = np.sum(r1_p * w1_p) + np.sum(r2_p * w2_p)
    assert abs(pairing_before - pairing_after) < mp.mpf("1e-34")

    # Tree-to-tree commutation in 35-digit precision:
    # rotate_jacobi_tree_coordinates(cartesian_to_jacobi(state, tree_A), tree_A, tree_B) == cartesian_to_jacobi(state, tree_B)
    from normal_cone_nbody_float128.core.jacobi import build_tree_from_spec
    from normal_cone_nbody_float128.core.transitions import rotate_jacobi_tree_coordinates

    masses = np.array([mA, mB, mC], dtype=object)
    tree_A = build_tree_from_spec(masses, ((0, 1), 2))
    tree_B = build_tree_from_spec(masses, (0, (1, 2)))
    q = np.array([
        [mp.mpf("1.5"), mp.mpf("2.5"), mp.mpf("0.5")],
        [mp.mpf("-1.2"), mp.mpf("-0.8"), mp.mpf("0.3")],
        [mp.mpf("0.7"), mp.mpf("-1.4"), mp.mpf("-0.6")],
    ], dtype=object)
    p = np.array([
        [mp.mpf("0.4"), mp.mpf("0.1"), mp.mpf("-0.2")],
        [mp.mpf("-0.3"), mp.mpf("-0.2"), mp.mpf("0.1")],
        [mp.mpf("-0.1"), mp.mpf("0.1"), mp.mpf("0.1")],
    ], dtype=object)
    state = ExtendedPhaseState(q=q, p=p, masses=masses)
    coords_A = cartesian_to_jacobi(state, tree_A)
    coords_B_direct = cartesian_to_jacobi(state, tree_B)

    pair_A = [n for n in tree_A.internal_nodes if n.body_indices == (0, 1)][0]
    pair_B = [n for n in tree_B.internal_nodes if n.body_indices == (1, 2)][0]
    coords_B_rot = rotate_jacobi_tree_coordinates(
        coords=coords_A,
        tree_A=tree_A,
        tree_B=tree_B,
        trans=trans,
        node_ids_A=(pair_A.id, tree_A.root.id),
        node_ids_B=(pair_B.id, tree_B.root.id),
    )
    for node in tree_B.internal_nodes:
        for d in range(3):
            assert abs(coords_B_rot.r[node.id][d] - coords_B_direct.r[node.id][d]) < mp.mpf("1e-33")
            assert abs(coords_B_rot.w[node.id][d] - coords_B_direct.w[node.id][d]) < mp.mpf("1e-33")
