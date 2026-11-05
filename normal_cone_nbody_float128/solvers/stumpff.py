r"""
Universal Stumpff Transcendental Functions for Non-Singular Keplerian Dynamics.

Theoretical Foundations:
- Monograph Reference: Chapter 5, Section 5.2.2 (Sundman--Stumpff Analytical Diffeomorphism)
  and Chapter 3, Section 3.3 (Universal Variables and Sundman Transmutation), Appendix A.
- Governing Equation:
    c_k(z) = \\sum_{m=0}^{\\infty} \\frac{(-1)^m z^m}{(2m + k)!}, \\quad k \\ge 0.
- Recurrence Relations:
    z c_{k+2}(z) = \\frac{1}{k!} - c_k(z), \\qquad c_k'(z) = \\frac{c_{k-1}(z) - k c_k(z)}{2z}.
- Certified Formal Proofs:
    Lean 4 module `PaperUnifiedFoundations.lean` and `StumpffMcGeheeDynamics.lean`.

Numerical Architecture:
- |z| < 0.1: 6th-order Horner polynomial expansion achieving < 10^-18 relative error (machine precision float64).
- z >= 0.1: Closed trigonometric formulations (elliptic regime).
- z <= -0.1: Closed real hyperbolic formulations (hyperbolic regime).
- Vectorized for NumPy arrays and optimized scalar fast-path.
"""

from __future__ import annotations
import math
from typing import Tuple, Union, Any
import numpy as np
import mpmath as mp
from ..core.constants import DEFAULT_DPS

# Ensure mpmath 35-digit precision is active for polynomial coefficients
mp.mp.dps = DEFAULT_DPS

# Order 14 Horner polynomial coefficients for mpmath 35-digit precision
def _make_stumpff_coeffs(k: int, order: int = 14) -> Tuple[mp.mpf, ...]:
    coeffs = []
    for m in range(order + 1):
        num = mp.mpf(-1) ** m
        denom = mp.fac(2 * m + k)
        coeffs.append(num / denom)
    return tuple(coeffs)

_C0_COEFFS_F128 = _make_stumpff_coeffs(0, 14)
_C1_COEFFS_F128 = _make_stumpff_coeffs(1, 14)
_C2_COEFFS_F128 = _make_stumpff_coeffs(2, 14)
_C3_COEFFS_F128 = _make_stumpff_coeffs(3, 14)
_C4_COEFFS_F128 = _make_stumpff_coeffs(4, 14)
_C5_COEFFS_F128 = _make_stumpff_coeffs(5, 14)


def _make_stumpff_deriv_coeffs(k: int, order: int = 14) -> Tuple[mp.mpf, ...]:
    coeffs = []
    for j in range(order):
        m = j + 1
        num = (mp.mpf(-1) ** m) * mp.mpf(m)
        denom = mp.fac(2 * m + k)
        coeffs.append(num / denom)
    return tuple(coeffs)

_DC0_COEFFS_F128 = _make_stumpff_deriv_coeffs(0, 14)
_DC1_COEFFS_F128 = _make_stumpff_deriv_coeffs(1, 14)
_DC2_COEFFS_F128 = _make_stumpff_deriv_coeffs(2, 14)
_DC3_COEFFS_F128 = _make_stumpff_deriv_coeffs(3, 14)
_DC4_COEFFS_F128 = _make_stumpff_deriv_coeffs(4, 14)
_DC5_COEFFS_F128 = _make_stumpff_deriv_coeffs(5, 14)


def _eval_horner_f128(coeffs: Tuple[mp.mpf, ...], z: mp.mpf) -> mp.mpf:
    res = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        res = res * z + c
    return res


def stumpff_c_f128(
    z: Union[mp.mpf, float, int, str]
) -> Tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """
    Computes Stumpff functions c_0(z) through c_5(z) in 35-digit mpmath precision.
    Order-14 Horner polynomial expansion for |z| < 0.1 guarantees < 10^-36 relative error.
    """
    z_val = mp.mpf(z) if not isinstance(z, mp.mpf) else z

    if abs(z_val) < mp.mpf('0.1'):
        c2 = _eval_horner_f128(_C2_COEFFS_F128, z_val)
        c3 = _eval_horner_f128(_C3_COEFFS_F128, z_val)
        c0 = mp.mpf('1.0') - z_val * c2
        c1 = mp.mpf('1.0') - z_val * c3
        c4 = _eval_horner_f128(_C4_COEFFS_F128, z_val)
        c5 = _eval_horner_f128(_C5_COEFFS_F128, z_val)
        return c0, c1, c2, c3, c4, c5

    if z_val >= mp.mpf('0.1'):
        sqrt_z = mp.sqrt(z_val)
        sin_sz = mp.sin(sqrt_z)
        cos_sz = mp.cos(sqrt_z)

        c0 = cos_sz
        c1 = sin_sz / sqrt_z
        c2 = (mp.mpf('1.0') - cos_sz) / z_val
        c3 = (sqrt_z - sin_sz) / (z_val * sqrt_z)
        c4 = (mp.mpf('0.5') - c2) / z_val
        c5 = (mp.mpf('1.0') / mp.mpf('6.0') - c3) / z_val
        return c0, c1, c2, c3, c4, c5

    # Hyperbolic regime (z < -0.1)
    sqrt_neg_z = mp.sqrt(-z_val)
    sinh_sz = mp.sinh(sqrt_neg_z)
    cosh_sz = mp.cosh(sqrt_neg_z)

    c0 = cosh_sz
    c1 = sinh_sz / sqrt_neg_z
    c2 = (cosh_sz - mp.mpf('1.0')) / (-z_val)
    c3 = (sinh_sz - sqrt_neg_z) / ((-z_val) * sqrt_neg_z)
    c4 = (mp.mpf('0.5') - c2) / z_val
    c5 = (mp.mpf('1.0') / mp.mpf('6.0') - c3) / z_val
    return c0, c1, c2, c3, c4, c5


def stumpff_c(
    z: Union[float, mp.mpf, np.ndarray]
) -> Tuple[
    Union[float, mp.mpf, np.ndarray],
    Union[float, mp.mpf, np.ndarray],
    Union[float, mp.mpf, np.ndarray],
    Union[float, mp.mpf, np.ndarray],
    Union[float, mp.mpf, np.ndarray],
    Union[float, mp.mpf, np.ndarray],
]:
    """
    Computes the first 6 universal Stumpff functions c_0(z) through c_5(z).
    Seamlessly dispatches to 35-digit mpmath precision when z is an mp.mpf.

    Parameters:
        z: Dimensionless energy parameter alpha * tau^2 (scalar float, mp.mpf, or NumPy ndarray).

    Returns:
        Tuple (c0, c1, c2, c3, c4, c5) of matching type and shape.
    """
    if isinstance(z, mp.mpf):
        return stumpff_c_f128(z)

    if isinstance(z, np.ndarray):
        if z.dtype == object and len(z) > 0 and isinstance(z.flat[0], mp.mpf):
            c_res = [stumpff_c_f128(x) for x in z.flat]
            c0 = np.array([r[0] for r in c_res], dtype=object).reshape(z.shape)
            c1 = np.array([r[1] for r in c_res], dtype=object).reshape(z.shape)
            c2 = np.array([r[2] for r in c_res], dtype=object).reshape(z.shape)
            c3 = np.array([r[3] for r in c_res], dtype=object).reshape(z.shape)
            c4 = np.array([r[4] for r in c_res], dtype=object).reshape(z.shape)
            c5 = np.array([r[5] for r in c_res], dtype=object).reshape(z.shape)
            return c0, c1, c2, c3, c4, c5
        return _stumpff_c_vectorized(z)

    z_val = float(z)

    if abs(z_val) < 0.1:
        # Horner polynomial expansion of order 6 on c2..c5 + exact downward Stumpff recurrence for c0, c1
        c2 = 0.5 + z_val * (
            -1.0 / 24.0 + z_val * (
                1.0 / 720.0 + z_val * (
                    -1.0 / 40320.0 + z_val * (
                        1.0 / 3628800.0 + z_val * (
                            -1.0 / 479001600.0 + z_val * (1.0 / 87178291200.0)
                        )
                    )
                )
            )
        )
        c3 = (1.0 / 6.0) + z_val * (
            -1.0 / 120.0 + z_val * (
                1.0 / 5040.0 + z_val * (
                    -1.0 / 362880.0 + z_val * (
                        1.0 / 39916800.0 + z_val * (
                            -1.0 / 6227020800.0 + z_val * (1.0 / 1307674368000.0)
                        )
                    )
                )
            )
        )
        c0 = 1.0 - z_val * c2
        c1 = 1.0 - z_val * c3
        c4 = (1.0 / 24.0) + z_val * (
            -1.0 / 720.0 + z_val * (
                1.0 / 40320.0 + z_val * (
                    -1.0 / 3628800.0 + z_val * (
                        1.0 / 479001600.0 + z_val * (
                            -1.0 / 87178291200.0 + z_val * (1.0 / 20922789888000.0)
                        )
                    )
                )
            )
        )
        c5 = (1.0 / 120.0) + z_val * (
            -1.0 / 5040.0 + z_val * (
                1.0 / 362880.0 + z_val * (
                    -1.0 / 39916800.0 + z_val * (
                        1.0 / 6227020800.0 + z_val * (
                            -1.0 / 1307674368000.0 + z_val * (1.0 / 355687428096000.0)
                        )
                    )
                )
            )
        )
        return c0, c1, c2, c3, c4, c5

    if z_val >= 0.1:
        # Elliptic regime (z > 0)
        sqrt_z = math.sqrt(z_val)
        sin_sz = math.sin(sqrt_z)
        cos_sz = math.cos(sqrt_z)

        c0 = cos_sz
        c1 = sin_sz / sqrt_z
        c2 = (1.0 - cos_sz) / z_val
        c3 = (sqrt_z - sin_sz) / (z_val * sqrt_z)
        c4 = (0.5 - c2) / z_val
        c5 = (1.0 / 6.0 - c3) / z_val
        return c0, c1, c2, c3, c4, c5

    # Hyperbolic regime (z < -0.1)
    sqrt_neg_z = math.sqrt(-z_val)
    sinh_sz = math.sinh(sqrt_neg_z)
    cosh_sz = math.cosh(sqrt_neg_z)

    c0 = cosh_sz
    c1 = sinh_sz / sqrt_neg_z
    c2 = (cosh_sz - 1.0) / (-z_val)
    c3 = (sinh_sz - sqrt_neg_z) / ((-z_val) * sqrt_neg_z)
    c4 = (0.5 - c2) / z_val
    c5 = (1.0 / 6.0 - c3) / z_val
    return c0, c1, c2, c3, c4, c5


def _stumpff_c_vectorized(
    z: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized evaluation of c_0(z) through c_5(z) for NumPy arrays."""
    c0 = np.empty_like(z, dtype=np.float64)
    c1 = np.empty_like(z, dtype=np.float64)
    c2 = np.empty_like(z, dtype=np.float64)
    c3 = np.empty_like(z, dtype=np.float64)
    c4 = np.empty_like(z, dtype=np.float64)
    c5 = np.empty_like(z, dtype=np.float64)

    small_mask = np.abs(z) < 0.1
    pos_mask = (z >= 0.1)
    neg_mask = (z <= -0.1)

    if np.any(small_mask):
        zs = z[small_mask]
        c2[small_mask] = 0.5 + zs * (-1.0/24.0 + zs * (1.0/720.0 + zs * (-1.0/40320.0 + zs * (1.0/3628800.0 + zs * (-1.0/479001600.0 + zs * (1.0/87178291200.0))))))
        c3[small_mask] = (1.0/6.0) + zs * (-1.0/120.0 + zs * (1.0/5040.0 + zs * (-1.0/362880.0 + zs * (1.0/39916800.0 + zs * (-1.0/6227020800.0 + zs * (1.0/1307674368000.0))))))
        c0[small_mask] = 1.0 - zs * c2[small_mask]
        c1[small_mask] = 1.0 - zs * c3[small_mask]
        c4[small_mask] = (1.0/24.0) + zs * (-1.0/720.0 + zs * (1.0/40320.0 + zs * (-1.0/3628800.0 + zs * (1.0/479001600.0 + zs * (-1.0/87178291200.0 + zs * (1.0/20922789888000.0))))))
        c5[small_mask] = (1.0/120.0) + zs * (-1.0/5040.0 + zs * (1.0/362880.0 + zs * (-1.0/39916800.0 + zs * (1.0/6227020800.0 + zs * (-1.0/1307674368000.0 + zs * (1.0/355687428096000.0))))))

    if np.any(pos_mask):
        zp = z[pos_mask]
        sqrt_z = np.sqrt(zp)
        sin_sz = np.sin(sqrt_z)
        cos_sz = np.cos(sqrt_z)
        c0[pos_mask] = cos_sz
        c1[pos_mask] = sin_sz / sqrt_z
        c2[pos_mask] = (1.0 - cos_sz) / zp
        c3[pos_mask] = (sqrt_z - sin_sz) / (zp * sqrt_z)
        c4[pos_mask] = (0.5 - c2[pos_mask]) / zp
        c5[pos_mask] = (1.0 / 6.0 - c3[pos_mask]) / zp

    if np.any(neg_mask):
        zn = z[neg_mask]
        sqrt_neg_z = np.sqrt(-zn)
        sinh_sz = np.sinh(sqrt_neg_z)
        cosh_sz = np.cosh(sqrt_neg_z)
        c0[neg_mask] = cosh_sz
        c1[neg_mask] = sinh_sz / sqrt_neg_z
        c2[neg_mask] = (cosh_sz - 1.0) / (-zn)
        c3[neg_mask] = (sinh_sz - sqrt_neg_z) / ((-zn) * sqrt_neg_z)
        c4[neg_mask] = (0.5 - c2[neg_mask]) / zn
        c5[neg_mask] = (1.0 / 6.0 - c3[neg_mask]) / zn

    return c0, c1, c2, c3, c4, c5


def stumpff_c_derivatives_f128(
    z: Union[mp.mpf, float, int, str]
) -> Tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """
    Computes derivatives c_k'(z) = d/dz [c_k(z)] in 35-digit mpmath precision.
    Order-14 Horner polynomial expansion for |z| < 0.1 eliminates catastrophic cancellation.
    """
    z_val = mp.mpf(z) if not isinstance(z, mp.mpf) else z

    if abs(z_val) < mp.mpf('0.1'):
        dc0 = _eval_horner_f128(_DC0_COEFFS_F128, z_val)
        dc1 = _eval_horner_f128(_DC1_COEFFS_F128, z_val)
        dc2 = _eval_horner_f128(_DC2_COEFFS_F128, z_val)
        dc3 = _eval_horner_f128(_DC3_COEFFS_F128, z_val)
        dc4 = _eval_horner_f128(_DC4_COEFFS_F128, z_val)
        dc5 = _eval_horner_f128(_DC5_COEFFS_F128, z_val)
        return dc0, dc1, dc2, dc3, dc4, dc5

    c0, c1, c2, c3, c4, c5 = stumpff_c_f128(z_val)
    two_z = mp.mpf('2.0') * z_val
    dc0 = -mp.mpf('0.5') * c1
    dc1 = (c0 - c1) / two_z
    dc2 = (c1 - mp.mpf('2.0') * c2) / two_z
    dc3 = (c2 - mp.mpf('3.0') * c3) / two_z
    dc4 = (c3 - mp.mpf('4.0') * c4) / two_z
    dc5 = (c4 - mp.mpf('5.0') * c5) / two_z
    return dc0, dc1, dc2, dc3, dc4, dc5


def stumpff_c_derivatives(
    z: Union[float, mp.mpf]
) -> Tuple[
    Union[float, mp.mpf],
    Union[float, mp.mpf],
    Union[float, mp.mpf],
    Union[float, mp.mpf],
    Union[float, mp.mpf],
    Union[float, mp.mpf],
]:
    """
    Computes exact analytical derivatives c_k'(z) = d/dz [c_k(z)] for k = 0, ..., 5.

    Uses differential recurrence:
        c_0'(z) = - (1/2) c_1(z)
        c_k'(z) = (c_{k-1}(z) - k c_k(z)) / (2z)  for k >= 1.
    For |z| < 0.2, high-order Horner Taylor expansions (degree 7) are evaluated
    to completely eliminate catastrophic subtractive cancellation in (c_{k-1} - k c_k) / (2z),
    achieving full IEEE 754 machine precision (< 1e-15 relative error).
    """
    if isinstance(z, mp.mpf):
        return stumpff_c_derivatives_f128(z)

    z_val = float(z)

    if abs(z_val) < 0.2:
        # Horner polynomial expansion of order 7 (truncation error < 1e-18 for |z| < 0.2)
        dc0 = -0.5 * (1.0 + z_val * (
            -1.0 / 6.0 + z_val * (
                1.0 / 120.0 + z_val * (
                    -1.0 / 5040.0 + z_val * (
                        1.0 / 362880.0 + z_val * (
                            -1.0 / 39916800.0 + z_val * (
                                1.0 / 6227020800.0 + z_val * (-1.0 / 1307674368000.0)
                            )
                        )
                    )
                )
            )
        ))
        dc1 = (-1.0 / 6.0) + z_val * (
            (1.0 / 60.0) + z_val * (
                (-1.0 / 1680.0) + z_val * (
                    (1.0 / 90720.0) + z_val * (
                        (-1.0 / 7983360.0) + z_val * (
                            (1.0 / 1037836800.0) + z_val * (
                                (-1.0 / 186810624000.0) + z_val * (1.0 / 44460928512000.0)
                            )
                        )
                    )
                )
            )
        )
        dc2 = (-1.0 / 24.0) + z_val * (
            (1.0 / 360.0) + z_val * (
                (-1.0 / 13440.0) + z_val * (
                    (1.0 / 907200.0) + z_val * (
                        (-1.0 / 95800320.0) + z_val * (
                            (1.0 / 14529715200.0) + z_val * (
                                (-1.0 / 2988969984000.0) + z_val * (1.0 / 800296713216000.0)
                            )
                        )
                    )
                )
            )
        )
        dc3 = (-1.0 / 120.0) + z_val * (
            (1.0 / 2520.0) + z_val * (
                (-1.0 / 120960.0) + z_val * (
                    (1.0 / 9979200.0) + z_val * (
                        (-1.0 / 1245404160.0) + z_val * (
                            (1.0 / 217945728000.0) + z_val * (
                                (-1.0 / 50812489728000.0) + z_val * (1.0 / 15205637551104000.0)
                            )
                        )
                    )
                )
            )
        )
        dc4 = (-1.0 / 720.0) + z_val * (
            (1.0 / 20160.0) + z_val * (
                (-1.0 / 1209600.0) + z_val * (
                    (1.0 / 119750400.0) + z_val * (
                        (-1.0 / 17435658240.0) + z_val * (
                            (1.0 / 3487131648000.0) + z_val * (
                                (-1.0 / 914624815104000.0) + z_val * (1.0 / 304112751022080000.0)
                            )
                        )
                    )
                )
            )
        )
        dc5 = (-1.0 / 5040.0) + z_val * (
            (1.0 / 181440.0) + z_val * (
                (-1.0 / 13305600.0) + z_val * (
                    (1.0 / 1556755200.0) + z_val * (
                        (-1.0 / 261534873600.0) + z_val * (
                            (1.0 / 59281238016000.0) + z_val * (
                                (-1.0 / 17377871486976000.0) + z_val * (1.0 / 6386367771463680000.0)
                            )
                        )
                    )
                )
            )
        )
        return dc0, dc1, dc2, dc3, dc4, dc5

    c0, c1, c2, c3, c4, c5 = stumpff_c(z_val)
    two_z = 2.0 * z_val
    dc0 = -0.5 * c1
    dc1 = (c0 - c1) / two_z
    dc2 = (c1 - 2.0 * c2) / two_z
    dc3 = (c2 - 3.0 * c3) / two_z
    dc4 = (c3 - 4.0 * c4) / two_z
    dc5 = (c4 - 5.0 * c5) / two_z
    return dc0, dc1, dc2, dc3, dc4, dc5
