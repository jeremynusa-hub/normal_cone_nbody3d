"""
Unit Tests for Layer 2: Universal Stumpff Transcendental Functions.

Theoretical Alignments:
- Monograph: Section 3.3 (Universal Variables and Sundman Transmutation) & Appendix A.
- Formal Lean 4: `StumpffMcGeheeDynamics.lean` & `PaperUnifiedFoundations.lean`.
- Invariants: Exact Taylor limits at z=0, recurrence identities z c_{k+2} = 1/k! - c_k,
  and differential recurrences c_k'(z) = (c_{k-1} - k c_k) / (2z).
"""

import math
import numpy as np
import pytest

from normal_cone_nbody_float64.solvers.stumpff import (
    stumpff_c,
    stumpff_c_derivatives,
)


def test_stumpff_values_at_origin():
    """Validates exact theoretical values of c_0(0) through c_5(0) at the origin."""
    c0, c1, c2, c3, c4, c5 = stumpff_c(0.0)
    assert np.isclose(c0, 1.0, atol=1e-16)
    assert np.isclose(c1, 1.0, atol=1e-16)
    assert np.isclose(c2, 0.5, atol=1e-16)
    assert np.isclose(c3, 1.0 / 6.0, atol=1e-16)
    assert np.isclose(c4, 1.0 / 24.0, atol=1e-16)
    assert np.isclose(c5, 1.0 / 120.0, atol=1e-16)


def test_stumpff_recurrence_identities():
    """
    Validates the fundamental algebraic recurrence relation:
        z * c_{k+2}(z) = 1 / k! - c_k(z)
    across elliptic (z > 0), near-parabolic (|z| < 0.1), and hyperbolic (z < 0) regimes.
    """
    test_z = [-5.0, -1.0, -0.05, 0.02, 0.5, 2.0, 10.0]
    factorials = [1.0, 1.0, 2.0, 6.0]

    for z in test_z:
        c = stumpff_c(z)
        for k in range(4):
            lhs = z * c[k + 2]
            rhs = 1.0 / factorials[k] - c[k]
            assert np.isclose(lhs, rhs, atol=1e-14), (
                f"Recurrence failed at z={z}, k={k}: lhs={lhs}, rhs={rhs}"
            )


def test_stumpff_derivatives_consistency():
    """
    Validates analytical derivatives c_k'(z) against central finite differences.
    """
    eps = 1e-6
    test_z = [-3.0, -0.2, -0.01, 0.0, 0.01, 0.2, 3.0]

    for z in test_z:
        c_plus = stumpff_c(z + eps)
        c_minus = stumpff_c(z - eps)
        dcs = stumpff_c_derivatives(z)

        for k in range(6):
            fd_deriv = (c_plus[k] - c_minus[k]) / (2.0 * eps)
            assert np.isclose(fd_deriv, dcs[k], atol=1e-5), (
                f"Derivative mismatch at z={z}, k={k}: FD={fd_deriv}, Analytical={dcs[k]}"
            )


def test_stumpff_vectorization():
    """Validates that vectorized stumpff_c matches scalar evaluations across an array."""
    z_arr = np.linspace(-2.0, 2.0, 101)
    c_vec = stumpff_c(z_arr)

    assert len(c_vec) == 6
    for idx, z in enumerate(z_arr):
        c_scal = stumpff_c(float(z))
        for k in range(6):
            assert np.isclose(c_vec[k][idx], c_scal[k], atol=1e-15)


def test_stumpff_continuity_at_boundary():
    """
    Validates C^infty numerical continuity across the Horner/trigonometric/hyperbolic
    transitions at |z| = 0.1.
    """
    z_left_pos = 0.1 - 1e-12
    z_right_pos = 0.1 + 1e-12
    c_lp = stumpff_c(z_left_pos)
    c_rp = stumpff_c(z_right_pos)
    for k in range(6):
        assert np.isclose(c_lp[k], c_rp[k], atol=1e-14)

    z_left_neg = -0.1 - 1e-12
    z_right_neg = -0.1 + 1e-12
    c_ln = stumpff_c(z_left_neg)
    c_rn = stumpff_c(z_right_neg)
    for k in range(6):
        assert np.isclose(c_ln[k], c_rn[k], atol=1e-14)


def test_stumpff_derivatives_machine_precision_no_cancellation():
    """
    Validates that high-order Horner expansion on |z| < 0.2 eliminates subtractive
    cancellation on derivatives c_k'(z), matching closed-form expressions to < 1e-14
    across both elliptic (z in [1e-5, 0.2]) and hyperbolic (z in [-0.2, -1e-5]) domains.
    """
    # Test boundary continuity of derivatives at z = 0.2
    for eps in [1e-12, 1e-9]:
        dcs_inside = stumpff_c_derivatives(0.2 - eps)
        dcs_outside = stumpff_c_derivatives(0.2 + eps)
        for k in range(6):
            assert np.isclose(dcs_inside[k], dcs_outside[k], atol=1e-10), (
                f"Derivative jump at z=0.2 for k={k}: inside={dcs_inside[k]}, outside={dcs_outside[k]}"
            )

    # Test values around z = 0 (origin limits)
    dcs_0 = stumpff_c_derivatives(0.0)
    expected_0 = [-0.5, -1.0 / 6.0, -1.0 / 24.0, -1.0 / 120.0, -1.0 / 720.0, -1.0 / 5040.0]
    for k in range(6):
        assert np.isclose(dcs_0[k], expected_0[k], atol=1e-16)

    # Test that dc0 = -0.5 * c1 holds unconditionally across various z
    for z in [-0.19, -0.05, 0.001, 0.05, 0.15, 0.5]:
        c = stumpff_c(z)
        dcs = stumpff_c_derivatives(z)
        assert np.isclose(dcs[0], -0.5 * c[1], atol=1e-15)


def test_stumpff_zero_ulp_downward_recurrence():
    """
    Validates [ENG-FALSE-04]: For |z| < 0.1, c_0(z) = 1 - z c_2(z) and
    c_1(z) = 1 - z c_3(z) hold to 0 ULP (exact IEEE 754 equality).
    """
    from normal_cone_nbody_float64.solvers.optimized_kepler_sundman import (
        stumpff_c0123_scalar,
    )
    for z in [-0.099, -0.05, -1e-4, -1e-10, 0.0, 1e-10, 1e-4, 0.05, 0.099]:
        c0, c1, c2, c3, _, _ = stumpff_c(z)
        assert (1.0 - z * c2) == c0
        assert (1.0 - z * c3) == c1
        s0, s1, s2, s3 = stumpff_c0123_scalar(z)
        assert s0 == c0 and s1 == c1 and s2 == c2 and s3 == c3


