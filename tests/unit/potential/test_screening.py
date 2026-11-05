"""
Unit Tests for Layer 3: Smooth Scale-Adapted Conormal Screening Cutoff Functions Theta(u_v).

Theoretical Alignments:
- Monograph: Section 4.2 (Conormal Screening and Apex Potential Regularization),
  Equations (4.58)-(4.71).
- Formal Lean 4: `01_ApexUniversalLipschitzBound_Theorems.lean` (`apex_universal_lipschitz_bound`).
- Invariants: Exact zero core Theta(u <= 1) = 0, exact asymptotic unity Theta(u >= 2) = 1,
  monotonicity Theta' >= 0, and C^1/C^2 smoothness.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.potential.screening import (
    eval_screening_theta,
    eval_screening_theta_vectorized,
    get_screening_mode_for_order,
)


def test_screening_modes_order_adaptation():
    """Validates that get_screening_mode_for_order maps integrator orders properly."""
    assert get_screening_mode_for_order(2) == "c2"
    assert get_screening_mode_for_order(4) == "c4"
    assert get_screening_mode_for_order(6) == "c6"
    assert get_screening_mode_for_order(8) == "c8"
    assert get_screening_mode_for_order(10) == "cauchy"


def test_screening_boundaries_and_values():
    """Validates the core and asymptotic boundary conditions for all supported modes."""
    for mode in ["cauchy", "hermite", "c2", "c4", "c6", "c8"]:
        # Core region: u <= 1.0 => Theta == 0.0, Theta' == 0.0
        for u in [0.0, 0.5, 0.999, 1.0]:
            th, dth = eval_screening_theta(u, mode=mode)
            assert th == 0.0, f"Failed at u={u}, mode={mode}"
            assert dth == 0.0, f"Failed deriv at u={u}, mode={mode}"

        # Asymptotic region: u >= 2.0 => Theta == 1.0, Theta' == 0.0
        for u in [2.0, 2.001, 5.0, 100.0]:
            th, dth = eval_screening_theta(u, mode=mode)
            assert th == 1.0, f"Failed at u={u}, mode={mode}"
            assert dth == 0.0, f"Failed deriv at u={u}, mode={mode}"

        # Midpoint: u = 1.5 => Theta(1.5) == 0.5 by symmetry
        th_mid, _ = eval_screening_theta(1.5, mode=mode)
        assert np.isclose(th_mid, 0.5, atol=1e-14)


def test_screening_derivatives_vs_finite_differences():
    """Validates analytical derivatives Theta'(u) against central finite differences."""
    eps = 1e-6
    test_u = np.linspace(1.05, 1.95, 20)

    for mode in ["cauchy", "hermite", "c2", "c4", "c6", "c8"]:
        for u in test_u:
            th_plus, _ = eval_screening_theta(u + eps, mode=mode)
            th_minus, _ = eval_screening_theta(u - eps, mode=mode)
            _, dth_analytical = eval_screening_theta(u, mode=mode)

            fd_deriv = (th_plus - th_minus) / (2.0 * eps)
            assert np.isclose(fd_deriv, dth_analytical, rtol=1e-5, atol=1e-6), (
                f"Derivative mismatch at u={u}, mode={mode}: FD={fd_deriv}, Analytical={dth_analytical}"
            )


def test_screening_monotonicity():
    """Validates that Theta'(u) >= 0 everywhere on [1.0, 2.0] (strictly non-decreasing)."""
    grid = np.linspace(1.0, 2.0, 1000)
    for mode in ["cauchy", "hermite", "c2", "c4", "c6", "c8"]:
        th_prev = -1.0
        for u in grid:
            th, dth = eval_screening_theta(float(u), mode=mode)
            assert 0.0 <= th <= 1.0
            assert dth >= -1e-15, f"Negative derivative at u={u}, mode={mode}: {dth}"
            assert th >= th_prev - 1e-15, f"Non-monotonic at u={u}, mode={mode}"
            th_prev = th


def test_screening_vectorized():
    """Validates vectorized evaluation against scalar evaluations."""
    u_arr = np.linspace(0.0, 3.0, 151)
    for mode in ["cauchy", "hermite", "c2", "c4", "c6", "c8"]:
        th_vec, dth_vec = eval_screening_theta_vectorized(u_arr, mode=mode)
        for i, u in enumerate(u_arr):
            th_s, dth_s = eval_screening_theta(float(u), mode=mode)
            assert np.isclose(th_vec[i], th_s, atol=1e-15)
            assert np.isclose(dth_vec[i], dth_s, atol=1e-15)


def test_c8_centered_spline_precision_and_symmetry():
    """
    Validates the centered formulation of the C^8 Hermite spline:
    1. Condition number kappa ~ 1 around s = 0.5 (u = 1.5) with exact value 0.5.
    2. Strict anti-symmetry around midpoint: Theta(1.5 + xi) + Theta(1.5 - xi) == 1.0 (< 1e-15).
    3. Seamless transition between centered form (|xi| <= 0.25) and boundary factored form (|xi| > 0.25).
    4. Exact value of peak derivative at u = 1.5: Theta'(1.5) == 109395 / 32768 == 3.338470458984375.
    5. Monotonicity without subtractive cancellation on a dense grid of 10,000 points.
    """
    # 1. Exact value at midpoint u = 1.5 (s = 0.5, xi = 0.0)
    th_mid, dth_mid = eval_screening_theta(1.5, mode="c8")
    assert th_mid == 0.5
    assert dth_mid == 3.338470458984375  # 109395 / 32768 exact dyadic fraction

    # 2. Strict anti-symmetry around midpoint across [0, 0.5]
    for xi in np.linspace(0.001, 0.499, 100):
        th_plus, _ = eval_screening_theta(1.5 + xi, mode="c8")
        th_minus, _ = eval_screening_theta(1.5 - xi, mode="c8")
        assert np.isclose(th_plus + th_minus, 1.0, atol=1e-15), (
            f"Anti-symmetry violated at xi={xi}: sum={th_plus + th_minus}"
        )

    # 3. Transition continuity at |xi| = 0.25 (u = 1.25 and u = 1.75)
    eps = 1e-14
    th_left, _ = eval_screening_theta(1.25 - eps, mode="c8")
    th_right, _ = eval_screening_theta(1.25 + eps, mode="c8")
    assert np.isclose(th_left, th_right, atol=1e-13), f"Transition jump at u=1.25: {abs(th_left - th_right)}"

    th_left_hi, _ = eval_screening_theta(1.75 - eps, mode="c8")
    th_right_hi, _ = eval_screening_theta(1.75 + eps, mode="c8")
    assert np.isclose(th_left_hi, th_right_hi, atol=1e-13), f"Transition jump at u=1.75: {abs(th_left_hi - th_right_hi)}"

    # 4. Dense grid monotonicity
    grid = np.linspace(1.0, 2.0, 10001)
    th_vals = [eval_screening_theta(float(u), mode="c8")[0] for u in grid]
    diffs = np.diff(th_vals)
    assert np.all(diffs >= 0.0), f"Negative slope detected in C8: min diff = {np.min(diffs)}"


def test_screening_constants_and_no_pbm6_reference():
    """
    Validates D-144-02:
    1. Canonical transition window constants SCREENING_U0 == 1.0, SCREENING_U1 == 2.0, SCREENING_SPAN == 1.0
       are defined in core.constants and potential.screening.
    2. Obsolete 'PBM6' reference is purged from potential/screening.py docstrings.
    """
    import inspect
    import normal_cone_nbody_float64.core.constants as const_f64
    import normal_cone_nbody_float128.core.constants as const_f128
    import normal_cone_nbody_float64.potential.screening as scr_f64
    import normal_cone_nbody_float128.potential.screening as scr_f128

    assert getattr(const_f64, "SCREENING_U0", None) == 1.0
    assert getattr(const_f64, "SCREENING_U1", None) == 2.0
    assert getattr(const_f64, "SCREENING_SPAN", None) == 1.0

    assert getattr(const_f128, "SCREENING_U0", None) == 1.0
    assert getattr(const_f128, "SCREENING_U1", None) == 2.0
    assert getattr(const_f128, "SCREENING_SPAN", None) == 1.0

    assert getattr(scr_f64, "SCREENING_U0", None) == 1.0
    assert getattr(scr_f64, "SCREENING_U1", None) == 2.0
    assert getattr(scr_f128, "SCREENING_U0", None) == 1.0
    assert getattr(scr_f128, "SCREENING_U1", None) == 2.0

    assert "PBM6" not in inspect.getsource(scr_f64)
    assert "PBM6" not in inspect.getsource(scr_f128)

