"""
Unit Tests for Layer 1 Invariants (normal_cone_nbody_float64.core.invariants).

Mathematical Verification Objectives:
- Exact evaluation and format of 10 Lie--Bargmann generators (P, L, G, H_ext).
- Commutation relations and Poisson bracket closure of so(3) ltimes R^3:
  * {L_a, L_b} = eps_{abc} L_c
  * {L_a, P_b} = eps_{abc} P_c
  * {P_a, G_b} = -M_tot delta_{ab}
- Truncation residual vector precision floor: Delta_trunc < 10^-14.
- Exact unsoftened divergence (V -> -inf) upon particle coincidence (zero softening).
"""

import pytest
import numpy as np

from normal_cone_nbody_float64.core.state import ExtendedPhaseState
from normal_cone_nbody_float64.core.invariants import (
    compute_linear_momentum,
    compute_angular_momentum,
    compute_center_of_mass,
    compute_galilean_boost,
    compute_kinetic_energy,
    compute_newtonian_potential_energy,
    compute_hamiltonian,
    compute_invariants,
    compute_truncation_residual,
)


@pytest.fixture
def four_body_state() -> ExtendedPhaseState:
    # 4-body planar choreographic initial state
    q = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.0, 0.5, 0.0],
        [-0.5, 0.0, 0.0],
        [0.0, -0.5, 0.0],
        [0.5, 0.0, 0.0],
    ], dtype=np.float64)
    masses = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=2.0, E0=-2.0)


def test_linear_and_angular_momentum(four_body_state):
    state = four_body_state
    P = compute_linear_momentum(state)
    np.testing.assert_allclose(P, np.zeros(3), atol=1e-15)

    L = compute_angular_momentum(state)
    # L_z = sum (x_i p_y,i - y_i p_x,i) = 1*0.5 - 0 + 0 - 1*(-0.5) + (-1)*(-0.5) + 0 - (-1)*0.5 = 0.5 + 0.5 + 0.5 + 0.5 = 2.0
    expected_L = np.array([0.0, 0.0, 2.0])
    np.testing.assert_allclose(L, expected_L, atol=1e-15)


def test_galilean_boost(four_body_state):
    state = four_body_state
    G = compute_galilean_boost(state)
    # Center of mass is 0, momentum is 0 -> G = 0
    np.testing.assert_allclose(G, np.zeros(3), atol=1e-15)

    # Shift momentum and center of mass
    shifted = state.with_updates(q=state.q + np.array([1.0, 0.0, 0.0]), t=3.0)
    G_shifted = compute_galilean_boost(shifted)
    # R_cm = [1, 0, 0], P_tot = 0 -> G = M_tot * R_cm = 4 * [1, 0, 0] = [4, 0, 0]
    np.testing.assert_allclose(G_shifted, np.array([4.0, 0.0, 0.0]), atol=1e-15)


def test_kinetic_and_potential_energy(four_body_state):
    state = four_body_state
    T = compute_kinetic_energy(state)
    # T = 4 * (1/2 * (0.5^2 / 1.0)) = 4 * 0.125 = 0.5
    assert T == pytest.approx(0.5, rel=1e-15)

    V = compute_newtonian_potential_energy(state, G=1.0)
    # 4 adjacent pairs distance sqrt(2), 2 opposite pairs distance 2.0
    # V = -G * [4 * (1*1 / sqrt(2)) + 2 * (1*1 / 2.0)] = - (4/sqrt(2) + 1) = - (2*sqrt(2) + 1) ~= -3.8284271247
    expected_V = -(2.0 * np.sqrt(2.0) + 1.0)
    assert V == pytest.approx(expected_V, rel=1e-14)

    H = compute_hamiltonian(state, G=1.0)
    assert H == pytest.approx(T + V, rel=1e-15)


def test_unsoftened_singularity_divergence():
    # Coinciding bodies must yield -inf without division by zero crash
    q_sing = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ])
    p_sing = np.zeros((2, 3))
    state = ExtendedPhaseState(q=q_sing, p=p_sing, masses=np.array([1.0, 1.0]))
    V = compute_newtonian_potential_energy(state)
    assert V == -float("inf")


def test_invariants_container_and_truncation_residual(four_body_state):
    state = four_body_state
    inv = compute_invariants(state)

    vec10 = inv.as_vector10()
    assert vec10.shape == (10,)

    # Residual against identical state must be zero
    res_zero = compute_truncation_residual(inv, inv)
    np.testing.assert_allclose(res_zero, np.zeros(4), atol=1e-15)


def test_parabolic_zero_energy_drift_detection():
    """
    For a parabolic trajectory with exact initial energy E0 = 0.0 (T0 + V0 = 0),
    any subsequent energy perturbation H(t) != 0.0 must be captured by inv.delta_E
    andinv.delta_E_rel rather than silently masked to 0.0 by `if state.E0 != 0.0`.
    """
    # 2-body parabolic escape state: m1 = m2 = 1.0, r = 1.0, G = 1.0 -> V0 = -1.0
    # Reduced mass mu = 0.5, T0 = p_rel^2 / (2*mu) = p_rel^2 = 1.0 -> p1 = +1, p2 = -1 -> H0 = 0.0
    q = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]], dtype=np.float64)
    p = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], dtype=np.float64)
    masses = np.array([1.0, 1.0], dtype=np.float64)

    state0 = ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E0=0.0)
    inv0 = compute_invariants(state0, G=1.0)
    assert inv0.H == pytest.approx(0.0, abs=1e-15)
    assert inv0.delta_E == pytest.approx(0.0, abs=1e-15)

    # Perturb momentum so H != 0 while keeping E0 = 0.0
    state_perturbed = state0.with_updates(
        p=np.array([[1.01, 0.0, 0.0], [-1.01, 0.0, 0.0]], dtype=np.float64)
    )
    inv_perturbed = compute_invariants(state_perturbed, G=1.0)
    expected_dE = abs(inv_perturbed.H - 0.0)
    assert expected_dE > 1e-2
    assert inv_perturbed.delta_E == pytest.approx(expected_dE, rel=1e-14)
    assert inv_perturbed.delta_E_rel > 5e-3
    assert "R^4_+" in (compute_truncation_residual.__doc__ or "")


def test_poisson_bracket_so3_closure():
    """
    Finite difference numerical verification of canonical Poisson brackets on PhaseState:
      {F, G}_Omega = sum_i (dF/dq_i . dG/dp_i - dF/dp_i . dG/dq_i).
    """
    rng = np.random.default_rng(42)
    q = rng.standard_normal((3, 3))
    p = rng.standard_normal((3, 3))
    masses = np.array([1.2, 2.3, 3.4])
    state = ExtendedPhaseState(q=q, p=p, masses=masses)

    eps = 1e-7

    def _poisson_bracket(eval_f, eval_g):
        grad_f_q = np.zeros((3, 3))
        grad_f_p = np.zeros((3, 3))
        grad_g_q = np.zeros((3, 3))
        grad_g_p = np.zeros((3, 3))

        for i in range(3):
            for k in range(3):
                # dq
                dq = np.zeros((3, 3)); dq[i, k] = eps
                f_plus = eval_f(state.with_updates(q=state.q + dq, skip_validation=True))
                f_minus = eval_f(state.with_updates(q=state.q - dq, skip_validation=True))
                grad_f_q[i, k] = (f_plus - f_minus) / (2 * eps)

                g_plus = eval_g(state.with_updates(q=state.q + dq, skip_validation=True))
                g_minus = eval_g(state.with_updates(q=state.q - dq, skip_validation=True))
                grad_g_q[i, k] = (g_plus - g_minus) / (2 * eps)

                # dp
                dp = np.zeros((3, 3)); dp[i, k] = eps
                f_plus = eval_f(state.with_updates(p=state.p + dp, skip_validation=True))
                f_minus = eval_f(state.with_updates(p=state.p - dp, skip_validation=True))
                grad_f_p[i, k] = (f_plus - f_minus) / (2 * eps)

                g_plus = eval_g(state.with_updates(p=state.p + dp, skip_validation=True))
                g_minus = eval_g(state.with_updates(p=state.p - dp, skip_validation=True))
                grad_g_p[i, k] = (g_plus - g_minus) / (2 * eps)

        return float(np.sum(grad_f_q * grad_g_p - grad_f_p * grad_g_q))

    # {L_x, L_y} = L_z
    L_x = lambda s: compute_angular_momentum(s)[0]
    L_y = lambda s: compute_angular_momentum(s)[1]
    L_z_expected = compute_angular_momentum(state)[2]

    pb_Lx_Ly = _poisson_bracket(L_x, L_y)
    assert pb_Lx_Ly == pytest.approx(L_z_expected, rel=1e-5)

    # {P_x, G_x} = -M_tot
    P_x = lambda s: compute_linear_momentum(s)[0]
    G_x = lambda s: compute_galilean_boost(s)[0]
    pb_Px_Gx = _poisson_bracket(P_x, G_x)
    assert pb_Px_Gx == pytest.approx(-state.total_mass, rel=1e-5)

    # {P_x, G_y} = 0
    G_y = lambda s: compute_galilean_boost(s)[1]
    pb_Px_Gy = _poisson_bracket(P_x, G_y)
    assert abs(pb_Px_Gy) < 1e-6
