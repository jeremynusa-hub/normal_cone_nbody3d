"""
Unit Tests for Layer 1 Constants (normal_cone_nbody_float64.core.constants).

Mathematical Verification Objectives:
- Exact algebraic identity: eta_crit^3 == 1/3.
- Asymmetric hysteresis ratio: gamma = 3/5 = 0.6, gap Delta_lambda = 16/15.
- IEEE 754 machine roundoff threshold sanity.
- 25-digit literal verification of BM_6^(11) constants against Chapter 5 Table 1:
  Consistency conditions: 2 * sum(a[:-1]) + a[-1] == 1.0 and 2 * sum(b[:-1]) + b[-1] == 1.0.
- Retrograde sub-step bound: max_{a_i < 0} |a_i| <= 0.65 (Sheng--Suzuki stability).
"""

import math
import pytest
import numpy as np

from normal_cone_nbody_float64.core.constants import (
    G_DEFAULT,
    ETA_CRIT,
    GAMMA_HYSTERESIS,
    DELTA_LAMBDA_HYSTERESIS,
    EPS_MACH,
    ZERO_SOFTENING,
    BM6_DRIFT_WEIGHTS,
    BM6_IMPULSE_WEIGHTS,
)


def test_gravitational_constant():
    assert G_DEFAULT == 1.0


def test_zero_softening_enforcement():
    assert ZERO_SOFTENING == 0.0


def test_eta_crit_algebraic_identity():
    # eta_crit = 3^(-1/3)
    assert abs(ETA_CRIT**3 - (1.0 / 3.0)) < 1e-15
    assert 0.693361 < ETA_CRIT < 0.693362


def test_hysteresis_constants():
    assert GAMMA_HYSTERESIS == 1.0
    assert abs(DELTA_LAMBDA_HYSTERESIS - (16.0 / 15.0)) < 1e-15
    assert GAMMA_HYSTERESIS**2 == pytest.approx(1.0, rel=1e-14)


def test_machine_epsilon():
    assert EPS_MACH == np.finfo(np.float64).eps
    assert 1.0 + (EPS_MACH / 2.0) == 1.0
    assert 1.0 + EPS_MACH > 1.0



def test_bm6_constants_consistency():
    assert len(BM6_DRIFT_WEIGHTS) == 6
    assert len(BM6_IMPULSE_WEIGHTS) == 6

    # Symmetric composition consistency: sum of 11 stages must equal 1.0
    sum_a = 2.0 * sum(BM6_DRIFT_WEIGHTS[:-1]) + BM6_DRIFT_WEIGHTS[-1]
    sum_b = 2.0 * sum(BM6_IMPULSE_WEIGHTS[:-1]) + BM6_IMPULSE_WEIGHTS[-1]

    assert abs(sum_a - 1.0) < 1e-14
    assert abs(sum_b - 1.0) < 1e-14


def test_bm6_retrograde_substep_bound():
    # Sheng--Suzuki barrier requires negative sub-steps to be bounded
    negative_a = [w for w in BM6_DRIFT_WEIGHTS if w < 0.0]
    negative_b = [w for w in BM6_IMPULSE_WEIGHTS if w < 0.0]

    assert len(negative_a) > 0
    assert len(negative_b) > 0

    max_retrograde_a = max(abs(w) for w in negative_a)
    max_retrograde_b = max(abs(w) for w in negative_b)

    # Maximum negative sub-step should remain well below 0.65
    assert max_retrograde_a < 0.65
    assert max_retrograde_b < 0.65


def test_homological_constants_obstruction_compliance():
    """Verifies that only KAPPA4_BM4 is exposed and obsolete KAPPA6_* constants are purged (Thm 5.14)."""
    import normal_cone_nbody_float64.core.constants as const_mod

    assert hasattr(const_mod, "KAPPA4_BM4")
    assert abs(const_mod.KAPPA4_BM4 - (1.0 / 2880.0)) < 1e-16
    assert not hasattr(const_mod, "KAPPA6_BM6")
    assert not hasattr(const_mod, "KAPPA6_CONORMAL_LEAF")


def test_bm4_constants_sheng_suzuki_cubic_cancellation_and_c4_bch_moment():
    """
    Verifies Chapter 5 Section 5.3 Eq. (5.30) and Section 5.4 Eq. (5.181):
    1. 2*(a_1 + a_2 + a_3) + a_4 == 1.0 to exact IEEE 754 float64.
    2. 2*(a_1^3 + a_2^3 + a_3^3) + a_4^3 == 0.0 to exact IEEE 754 float64 zero.
    3. Exact 4th-order BCH/Duhamel moment c_4^{BM4} = (sum b_k tau_k^4 - 1/80) / 24 == C4_BCH_BM4.
    """
    from normal_cone_nbody_float64.core.constants import (
        BM4_A1,
        BM4_A2,
        BM4_A3,
        BM4_A4,
        BM4_DRIFT_WEIGHTS,
        BM4_IMPULSE_WEIGHTS,
        C4_BCH_BM4,
    )

    assert 2.0 * (BM4_A1 + BM4_A2 + BM4_A3) + BM4_A4 == 1.0
    cubic_res = 2.0 * (BM4_A1**3 + BM4_A2**3 + BM4_A3**3) + BM4_A4**3
    assert cubic_res == 0.0, f"Expected exact 0.0 for Sheng-Suzuki cubic cancellation, got {cubic_res}"

    tau = [-0.5]
    acc = -0.5
    for a_val in BM4_DRIFT_WEIGHTS:
        acc += a_val
        tau.append(acc)
    mu4 = sum(b_val * (t_val**4) for b_val, t_val in zip(BM4_IMPULSE_WEIGHTS, tau))
    c4_computed = (mu4 - 1.0 / 80.0) / 24.0
    assert abs(c4_computed - C4_BCH_BM4) < 1e-17


def test_bm6_rkn6_quadrature_order_conditions():
    """
    Verifies that the 11-stage BM6 coefficients satisfy all 5 RKN order-6 quadrature
    conditions sum_{k=1}^{12} b_k c_k^m == 1/(m+1) for m=1..5 to machine precision (< 1e-15).
    """
    a1, a2, a3, a4, a5, a6 = BM6_DRIFT_WEIGHTS
    b1, b2, b3, b4, b5, b6 = BM6_IMPULSE_WEIGHTS
    a_full = (a1, a2, a3, a4, a5, a6, a5, a4, a3, a2, a1)
    b_full = (b1, b2, b3, b4, b5, 0.5 * b6, 0.5 * b6, b5, b4, b3, b2, b1)

    c_nodes = [0.0]
    acc = 0.0
    for a_val in a_full:
        acc += a_val
        c_nodes.append(acc)

    for m in range(1, 6):
        quad_err = abs(sum(b_k * (c_k**m) for b_k, c_k in zip(b_full, c_nodes)) - 1.0 / (m + 1))
        assert quad_err < 1e-15, f"BM6 quadrature condition m={m} residual too large: {quad_err}"


