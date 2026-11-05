"""
Smooth Scale-Adapted Conormal Screening Cutoff Functions Theta(u_v).

Theoretical Foundations:
- Monograph Reference: Section 4.2 (Conormal Screening and Apex Potential Regularization),
  Equations (4.58)-(4.71).
- Governing Variable:
    u_v = ||r_v|| / (bar_rho_{L(v)} + bar_rho_{R(v)}),
    where bar_rho_C = rho_C / sqrt{M_C} in [L] is the mass-normalized hyper-radius,
    rendering u_v strictly dimensionless and scale-invariant under homotheties.
- Canonical Regimes:
    - u <= 1.0: Theta(u) = 0.0, Theta'(u) = 0.0 (Compact conormal core: ghost pole eradicated).
    - 1.0 < u < 2.0: Smooth monotonic transition from 0 to 1.
    - u >= 2.0: Theta(u) = 1.0, Theta'(u) = 0.0 (Pure asymptotic multipolar regime).
- Supported Modes:
    1. "cauchy" (C^infty smooth bump mollifier):
       g1(x) = exp(-1/x), g2(x) = exp(-1/(1-x)), Theta = g1 / (g1 + g2).
    2. "hermite" (C^2 quintic polynomial spline):
       sigma(u) = 10 s^3 - 15 s^4 + 6 s^5, where s = u - 1 in [0, 1].
- Certified Formal Proofs:
    Lean 4 modules `PaperUnifiedFoundations.lean` and `01_ApexUniversalLipschitzBound_Theorems.lean`.
"""

from __future__ import annotations
import math
from typing import Tuple, Union, Any
import numpy as np
import mpmath as mp
from ..core.constants import (
    DEFAULT_DPS,
    SCREENING_U0,
    SCREENING_U1,
    SCREENING_SPAN,
    SCREENING_U0_F128,
    SCREENING_U1_F128,
    SCREENING_SPAN_F128,
)

mp.mp.dps = DEFAULT_DPS


def get_screening_mode_for_order(order: int) -> str:
    """
    Returns the mathematically adapted screening mode for a given symplectic integrator order.
    Ensures backward error analysis (BCH) is not truncated by derivative jump discontinuities.
    """
    if order <= 2:
        return "c2"
    elif order <= 4:
        return "c4"
    elif order <= 6:
        return "c6"
    elif order <= 8:
        return "c8"
    else:
        return "cauchy"


# Centered C8 polynomial coefficients for Q(z) where z = xi^2 = (s - 0.5)^2:
# sigma8(s) = 0.5 + xi * Q(z), with xi = s - 0.5 in [-0.25, 0.25].
# Horner form: Q(z) = ((((((((c8*z + c7)*z + c6)*z + c5)*z + c4)*z + c3)*z + c2)*z + c1)*z + c0)
# All coefficients are exact binary rationals (power of 2 denominators) with condition number kappa ~ 1.
_C8_CENTERED_COEFFS: Tuple[float, ...] = (
    12870.0,
    -29172.0,
    29452.5,
    -17403.75,
    6647.265625,
    -1709.296875,
    299.126953125,
    -35.6103515625,
    3.338470458984375,
)

_C8_CENTERED_COEFFS_F128: Tuple[mp.mpf, ...] = tuple(mp.mpf(str(c)) for c in _C8_CENTERED_COEFFS)


def eval_screening_theta_f128(
    u: Union[mp.mpf, float, str, int],
    mode: str = "cauchy",
) -> Tuple[mp.mpf, mp.mpf]:
    """
    Evaluates the scale-adapted conormal screening function Theta(u) and its derivative Theta'(u)
    in 35-digit mpmath precision.
    """
    u_val = mp.mpf(u) if not isinstance(u, mp.mpf) else u

    if u_val <= SCREENING_U0_F128:
        return mp.mpf('0.0'), mp.mpf('0.0')
    if u_val >= SCREENING_U1_F128:
        return mp.mpf('1.0'), mp.mpf('0.0')

    s = (u_val - SCREENING_U0_F128) / SCREENING_SPAN_F128

    if mode in ("c2", "hermite", "order_2"):
        p = s if s <= mp.mpf('0.5') else (mp.mpf('1.0') - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        val = mp.mpf('10.0') * p3 - mp.mpf('15.0') * p4 + mp.mpf('6.0') * p5
        theta = val if s <= mp.mpf('0.5') else (mp.mpf('1.0') - val)
        d_theta = mp.mpf('30.0') * ((s * (mp.mpf('1.0') - s)) ** 2)
        return theta, d_theta

    if mode in ("c4", "order_4"):
        p = s if s <= mp.mpf('0.5') else (mp.mpf('1.0') - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        poly = mp.mpf('126.0') - mp.mpf('420.0') * p + mp.mpf('540.0') * p2 - mp.mpf('315.0') * p3 + mp.mpf('70.0') * p4
        val = p5 * poly
        theta = val if s <= mp.mpf('0.5') else (mp.mpf('1.0') - val)
        d_theta = mp.mpf('630.0') * ((s * (mp.mpf('1.0') - s)) ** 4)
        return theta, d_theta

    if mode in ("c6", "order_6"):
        p = s if s <= mp.mpf('0.5') else (mp.mpf('1.0') - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        p6 = p5 * p
        p7 = p6 * p
        poly = (
            mp.mpf('1716.0')
            - mp.mpf('9009.0') * p
            + mp.mpf('20020.0') * p2
            - mp.mpf('24024.0') * p3
            + mp.mpf('16380.0') * p4
            - mp.mpf('6006.0') * p5
            + mp.mpf('924.0') * p6
        )
        val = p7 * poly
        theta = val if s <= mp.mpf('0.5') else (mp.mpf('1.0') - val)
        d_theta = mp.mpf('12012.0') * ((s * (mp.mpf('1.0') - s)) ** 6)
        return theta, d_theta

    if mode in ("c8", "order_8"):
        xi = s - mp.mpf('0.5')
        if abs(xi) <= mp.mpf('0.25'):
            z = xi * xi
            q = _C8_CENTERED_COEFFS_F128[0]
            for c in _C8_CENTERED_COEFFS_F128[1:]:
                q = q * z + c
            theta = mp.mpf('0.5') + xi * q
        else:
            p = s if s < mp.mpf('0.5') else (mp.mpf('1.0') - s)
            p2 = p * p
            p3 = p2 * p
            p4 = p3 * p
            p5 = p4 * p
            p6 = p5 * p
            p7 = p6 * p
            p8 = p7 * p
            p9 = p8 * p
            poly = (
                mp.mpf('24310.0')
                - mp.mpf('175032.0') * p
                + mp.mpf('556920.0') * p2
                - mp.mpf('1021020.0') * p3
                + mp.mpf('1178100.0') * p4
                - mp.mpf('875160.0') * p5
                + mp.mpf('408408.0') * p6
                - mp.mpf('109395.0') * p7
                + mp.mpf('12870.0') * p8
            )
            base = p9 * poly
            theta = base if s < mp.mpf('0.5') else (mp.mpf('1.0') - base)
        d_theta = mp.mpf('218790.0') * ((s * (mp.mpf('1.0') - s)) ** 8)
        return theta, d_theta

    if mode in ("cauchy", "smooth"):
        x = s
        inv_x = mp.mpf('1.0') / x
        inv_1mx = mp.mpf('1.0') / (mp.mpf('1.0') - x)
        e1 = mp.exp(-inv_x)
        e2 = mp.exp(-inv_1mx)
        sum_e = e1 + e2
        if sum_e <= mp.mpf('0.0'):
            return mp.mpf('0.5'), mp.mpf('0.0')
        theta = e1 / sum_e
        d_theta = (e1 * e2 * (inv_x * inv_x + inv_1mx * inv_1mx)) / (sum_e * sum_e)
        return theta, d_theta

    raise ValueError(
        f"Unknown screening mode '{mode}'. Choose 'c2', 'c4', 'c6', 'c8', or 'cauchy'."
    )


def eval_screening_theta(
    u: Union[float, mp.mpf],
    mode: str = "cauchy",
) -> Tuple[Union[float, mp.mpf], Union[float, mp.mpf]]:
    """
    Evaluates the scale-adapted conormal screening function Theta(u) and its derivative Theta'(u).

    Parameters:
        u: Dimensionless scale ratio ||r_v|| / (bar_rho_L + bar_rho_R).
        mode: Screening profile:
              - "c2" / "hermite": C^2 quintic Hermite spline (minimal for Strang order 2).
              - "c4": C^4 nonic Hermite spline (derivatives 1..4 vanish at edges, adapted for BM4/Yoshida).
              - "c6": C^6 tridecic Hermite spline (derivatives 1..6 vanish at edges, adapted for BM6/PBM4).
              - "c8": C^8 heptadecic Hermite spline (derivatives 1..8 vanish at edges, adapted for order-8 compositions).
              - "cauchy" / "smooth": C^infty Cauchy mollifier (canonical theoretical baseline).

    Returns:
        (theta, d_theta): Value of cutoff in [0, 1] and its first derivative.
    """
    if isinstance(u, mp.mpf):
        return eval_screening_theta_f128(u, mode=mode)

    if u <= SCREENING_U0:
        return 0.0, 0.0
    if u >= SCREENING_U1:
        return 1.0, 0.0

    s = (u - SCREENING_U0) / SCREENING_SPAN

    if mode in ("c2", "hermite", "order_2"):
        p = s if s <= 0.5 else (1.0 - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        val = 10.0 * p3 - 15.0 * p4 + 6.0 * p5
        theta = val if s <= 0.5 else (1.0 - val)
        d_theta = 30.0 * (s * (1.0 - s)) ** 2
        return theta, d_theta

    if mode in ("c4", "order_4"):
        p = s if s <= 0.5 else (1.0 - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        poly = 126.0 - 420.0 * p + 540.0 * p2 - 315.0 * p3 + 70.0 * p4
        val = p5 * poly
        theta = val if s <= 0.5 else (1.0 - val)
        d_theta = 630.0 * ((s * (1.0 - s)) ** 4)
        return theta, d_theta

    if mode in ("c6", "order_6"):
        p = s if s <= 0.5 else (1.0 - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        p6 = p5 * p
        p7 = p6 * p
        poly = (
            1716.0
            - 9009.0 * p
            + 20020.0 * p2
            - 24024.0 * p3
            + 16380.0 * p4
            - 6006.0 * p5
            + 924.0 * p6
        )
        val = p7 * poly
        theta = val if s <= 0.5 else (1.0 - val)
        d_theta = 12012.0 * ((s * (1.0 - s)) ** 6)
        return theta, d_theta

    if mode in ("c8", "order_8"):
        xi = s - 0.5
        if abs(xi) <= 0.25:
            # Centered form around s = 0.5 (kappa ~ 1): xi in [-0.25, 0.25], z in [0, 0.0625]
            z = xi * xi
            q = _C8_CENTERED_COEFFS[0]
            for c in _C8_CENTERED_COEFFS[1:]:
                q = q * z + c
            theta = 0.5 + xi * q
        else:
            # Factored boundary form near s = 0 and s = 1 (p < 0.25, zero cancellation)
            p = s if s < 0.5 else (1.0 - s)
            p2 = p * p
            p3 = p2 * p
            p4 = p3 * p
            p5 = p4 * p
            p6 = p5 * p
            p7 = p6 * p
            p8 = p7 * p
            p9 = p8 * p
            poly = (
                24310.0
                - 175032.0 * p
                + 556920.0 * p2
                - 1021020.0 * p3
                + 1178100.0 * p4
                - 875160.0 * p5
                + 408408.0 * p6
                - 109395.0 * p7
                + 12870.0 * p8
            )
            base = p9 * poly
            theta = base if s < 0.5 else (1.0 - base)
        d_theta = 218790.0 * ((s * (1.0 - s)) ** 8)
        return theta, d_theta

    if mode in ("cauchy", "smooth"):
        x = s
        inv_x = 1.0 / x
        inv_1mx = 1.0 / (1.0 - x)

        e1 = math.exp(-inv_x) if inv_x < 700.0 else 0.0
        e2 = math.exp(-inv_1mx) if inv_1mx < 700.0 else 0.0
        sum_e = e1 + e2

        if sum_e <= 0.0:
            return 0.5, 0.0

        theta = e1 / sum_e
        d_theta = (e1 * e2 * (inv_x * inv_x + inv_1mx * inv_1mx)) / (sum_e * sum_e)
        return theta, d_theta

    raise ValueError(
        f"Unknown screening mode '{mode}'. Choose 'c2', 'c4', 'c6', 'c8', or 'cauchy'."
    )


def eval_screening_theta_vectorized(
    u_arr: np.ndarray,
    mode: str = "cauchy",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized evaluation of Theta(u) and Theta'(u) for NumPy arrays.
    """
    if u_arr.dtype == object and len(u_arr) > 0 and isinstance(u_arr.flat[0], mp.mpf):
        results = [eval_screening_theta_f128(x, mode=mode) for x in u_arr.flat]
        theta = np.array([r[0] for r in results], dtype=object).reshape(u_arr.shape)
        d_theta = np.array([r[1] for r in results], dtype=object).reshape(u_arr.shape)
        return theta, d_theta

    theta = np.zeros_like(u_arr, dtype=np.float64)
    d_theta = np.zeros_like(u_arr, dtype=np.float64)

    theta[u_arr >= 2.0] = 1.0

    trans_mask = (u_arr > 1.0) & (u_arr < 2.0)
    if not np.any(trans_mask):
        return theta, d_theta

    u_trans = u_arr[trans_mask]
    s = u_trans - 1.0

    if mode in ("c2", "hermite", "order_2"):
        p = np.where(s <= 0.5, s, 1.0 - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        val = 10.0 * p3 - 15.0 * p4 + 6.0 * p5
        theta[trans_mask] = np.where(s <= 0.5, val, 1.0 - val)
        d_theta[trans_mask] = 30.0 * ((s * (1.0 - s)) ** 2)
    elif mode in ("c4", "order_4"):
        p = np.where(s <= 0.5, s, 1.0 - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        poly = 126.0 - 420.0 * p + 540.0 * p2 - 315.0 * p3 + 70.0 * p4
        val = p5 * poly
        theta[trans_mask] = np.where(s <= 0.5, val, 1.0 - val)
        d_theta[trans_mask] = 630.0 * ((s * (1.0 - s)) ** 4)
    elif mode in ("c6", "order_6"):
        p = np.where(s <= 0.5, s, 1.0 - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        p6 = p5 * p
        p7 = p6 * p
        poly = (
            1716.0
            - 9009.0 * p
            + 20020.0 * p2
            - 24024.0 * p3
            + 16380.0 * p4
            - 6006.0 * p5
            + 924.0 * p6
        )
        val = p7 * poly
        theta[trans_mask] = np.where(s <= 0.5, val, 1.0 - val)
        d_theta[trans_mask] = 12012.0 * ((s * (1.0 - s)) ** 6)
    elif mode in ("c8", "order_8"):
        xi = s - 0.5
        abs_xi = np.abs(xi)
        center_mask = abs_xi <= 0.25

        # 1. Centered evaluation on |xi| <= 0.25 (kappa ~ 1)
        z = xi * xi
        q = np.full_like(z, _C8_CENTERED_COEFFS[0])
        for c in _C8_CENTERED_COEFFS[1:]:
            q = q * z + c
        val_centered = 0.5 + xi * q

        # 2. Factored boundary evaluation on |xi| > 0.25 (zero cancellation)
        p = np.where(s <= 0.5, s, 1.0 - s)
        p2 = p * p
        p3 = p2 * p
        p4 = p3 * p
        p5 = p4 * p
        p6 = p5 * p
        p7 = p6 * p
        p8 = p7 * p
        p9 = p8 * p
        poly = (
            24310.0
            - 175032.0 * p
            + 556920.0 * p2
            - 1021020.0 * p3
            + 1178100.0 * p4
            - 875160.0 * p5
            + 408408.0 * p6
            - 109395.0 * p7
            + 12870.0 * p8
        )
        base = p9 * poly
        val_boundary = np.where(s <= 0.5, base, 1.0 - base)

        theta[trans_mask] = np.where(center_mask, val_centered, val_boundary)
        d_theta[trans_mask] = 218790.0 * ((s * (1.0 - s)) ** 8)
    elif mode in ("cauchy", "smooth"):
        for idx in np.flatnonzero(trans_mask):
            th, dth = eval_screening_theta(float(u_arr.ravel()[idx]), mode="cauchy")
            theta.ravel()[idx] = th
            d_theta.ravel()[idx] = dth
    else:
        raise ValueError(
            f"Unknown screening mode '{mode}'. Choose 'c2', 'c4', 'c6', 'c8', or 'cauchy'."
        )

    return theta, d_theta
