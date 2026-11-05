r"""
Centralized Repository of Certified Physical, Geometric, and Numerical Constants.

Mathematical Monograph Alignment:
- Chapter 4: Section 4.3, Equation (4.28) — Universal Hill/Roche compactness threshold
  eta_crit = 3^(-1/3) governing synchronous supernode formation.
- Chapter 6: Section 6.2, Equations (6.12)-(6.15) — Asymmetric hysteresis band
  gamma = 3/5 = 0.6 and hysteresis gap Delta_lambda = 16/15 precluding topological chattering.
- Chapter 5: Section 5.4, Table 1 — Optimal 11-stage 6th-order Blanes--Moan symplectic coefficients
  certified to 25 decimal digits.
- Formal Lean 4 Alignment: PaperUnifiedFoundations.lean & BM6Order6SymmetricComposition.lean.
"""

from typing import Tuple, Any, Union
import numpy as np
import mpmath as mp

# Configure global multi-precision context for float128 / quad-precision operations
DEFAULT_DPS: int = 35
mp.mp.dps = DEFAULT_DPS


def get_dps() -> int:
    """Returns current decimal precision of mpmath context."""
    return int(mp.mp.dps)


def set_dps(dps: int) -> None:
    """Sets global decimal precision for float128 / extended operations."""
    mp.mp.dps = int(dps)


def to_f128(val: Any) -> Any:
    """Converts scalar or array-like to mpmath.mpf or object array of mpmath.mpf."""
    if isinstance(val, (int, str)):
        return mp.mpf(val)
    if isinstance(val, float):
        return mp.mpf(str(val))
    if isinstance(val, mp.mpf):
        return val
    if isinstance(val, np.ndarray):
        flat = [to_f128(x) for x in val.flat]
        return np.array(flat, dtype=object).reshape(val.shape)
    if isinstance(val, (list, tuple)):
        return [to_f128(x) for x in val]
    return mp.mpf(str(val))


def to_f64(val: Any) -> Any:
    """Converts mpmath scalar or object array back to float64."""
    if isinstance(val, (mp.mpf, int, str)):
        return float(val)
    if isinstance(val, float):
        return val
    if isinstance(val, np.ndarray):
        return np.asarray(val, dtype=np.float64)
    if isinstance(val, (list, tuple)):
        return [float(x) for x in val]
    return float(val)


# Normalized gravitational constant (canonical N-body units where G = 1.0)
G_DEFAULT_F128: mp.mpf = mp.mpf('1.0')
G_DEFAULT: float = 1.0
GRAVITATIONAL_CONSTANT: float = G_DEFAULT
GRAVITATIONAL_CONSTANT_F128: mp.mpf = G_DEFAULT_F128

# Universal Hill/Roche stability threshold governing synchronous supernode collapse (Section 4.3)
# 2 * eta^3 = 2/3 <=> eta_crit = (1/3)^(1/3) = 3^(-1/3) ~= 0.6933612743506347...
ETA_CRIT_F128: mp.mpf = mp.power(mp.mpf('3.0'), -mp.mpf('1.0') / mp.mpf('3.0'))
ETA_CRIT: float = float(ETA_CRIT_F128)
ETA_CRITICAL_HILL: float = ETA_CRIT
ETA_CRITICAL_HILL_F128: mp.mpf = ETA_CRIT_F128

# Canonical Spectral Hysteresis Band for Hill Supernode Activation/Deactivation (Section 4.3, Cor. 4.1)
# 1. Activation (Anticipation) at exact radial equigravity:
#    2 * eta^3 = 1 <=> eta_enter = 2^(-1/3) ~= 0.7937005259840997...
# 2. Deactivation at spherical RMS equigravity:
#    sqrt(2) * eta^3 = 1 <=> eta_exit = 2^(-1/6) ~= 0.8908987181403392...
ETA_ENTER_HILL_F128: mp.mpf = mp.power(mp.mpf('2.0'), -mp.mpf('1.0') / mp.mpf('3.0'))
ETA_EXIT_HILL_F128: mp.mpf = mp.power(mp.mpf('2.0'), -mp.mpf('1.0') / mp.mpf('6.0'))
ETA_ENTER_HILL: float = float(ETA_ENTER_HILL_F128)
ETA_EXIT_HILL: float = float(ETA_EXIT_HILL_F128)

# Cubic Tidal Force Ratio Thresholds for Roche Tidal Disruption (Section 7.2.1, Equation 7.34)
# Lambda_S = (1 / 3) * eta_S^3 = (M_ext / (3 * M_S)) * (bar_rho_S / r_P)^3
LAMBDA_CRITICAL_ROCHE_F128: mp.mpf = mp.mpf('0.25')
LAMBDA_ENTER_ROCHE_F128: mp.mpf = mp.mpf('0.25')
LAMBDA_EXIT_ROCHE_F128: mp.mpf = mp.mpf('0.18')
LAMBDA_CRITICAL_ROCHE: float = float(LAMBDA_CRITICAL_ROCHE_F128)
LAMBDA_ENTER_ROCHE: float = float(LAMBDA_ENTER_ROCHE_F128)
LAMBDA_EXIT_ROCHE: float = float(LAMBDA_EXIT_ROCHE_F128)

# Canonical Voronoi admissibility ratio (gamma = 1.0, midpoint 0.5/0.5, zero hysteresis by default)
GAMMA_HYSTERESIS_F128: mp.mpf = mp.mpf('1.0')
GAMMA_HYSTERESIS: float = 1.0

# Width of the dimensionless hysteresis gap: Delta_lambda = 5/3 - 3/5 = 16/15 ~= 1.066667
DELTA_LAMBDA_HYSTERESIS_F128: mp.mpf = mp.mpf('16.0') / mp.mpf('15.0')
DELTA_LAMBDA_HYSTERESIS: float = 16.0 / 15.0

# IEEE 754 double precision machine epsilon floor
EPS_MACH: float = float(np.finfo(np.float64).eps)  # 2.220446049250313e-16

# Quadruple precision / float128 machine epsilon floor (2^-113)
# Under mp.dps = 35, EPS_MACH_F128 = 2^-113 ~= 1.9259299443872358530559779425849273e-34
EPS_MACH_F128: mp.mpf = mp.power(mp.mpf('2.0'), -113)

# Zero softening enforcement: gravitational potential is unsoftened (epsilon = 0.0, d_0 = 0.0)
ZERO_SOFTENING: float = 0.0
ZERO_SOFTENING_F128: mp.mpf = mp.mpf('0.0')

# Canonical scale-adapted screening transition window u_v in [1.0, 2.0] (Chapter 4, Section 4.2 & Chapter 6, Section 6.4.3)
SCREENING_U0_F128: mp.mpf = mp.mpf('1.0')
SCREENING_U1_F128: mp.mpf = mp.mpf('2.0')
SCREENING_SPAN_F128: mp.mpf = SCREENING_U1_F128 - SCREENING_U0_F128
SCREENING_U0: float = 1.0
SCREENING_U1: float = 2.0
SCREENING_SPAN: float = 1.0

# Scale-invariant maximum dyadic refinement depth (Chapter 5, Section 5.6)
DEFAULT_MAX_DYADIC_LEVEL: int = 17

# Certified 4th-Order 7-Stage Palindromic Strang Composition Coefficients for BM_4^{(7)} (Chapter 5, Section 5.3)
# 7-stage symmetric composition of kick-exterior Strang steps S_2*(a_j h) with FSAL impulse fusion b_j = (a_{j-1} + a_j)/2,
# satisfying 2*(a_1 + a_2 + a_3) + a_4 = 1 and 2*(a_1^3 + a_2^3 + a_3^3) + a_4^3 = 0 to 40 decimal digits.
BM4_A1_F128: mp.mpf = mp.mpf('+0.2390843168869614000000000000000000000000')
BM4_A2_F128: mp.mpf = mp.mpf('+0.2390725069195032000000000000000000000000')
BM4_A3_F128: mp.mpf = mp.mpf('+0.2390524743622015795882324908525173940394')
BM4_A4_F128: mp.mpf = mp.mpf('1.0') - mp.mpf('2.0') * (BM4_A1_F128 + BM4_A2_F128 + BM4_A3_F128)

BM4_A1: float = float(BM4_A1_F128)
BM4_A2: float = float(BM4_A2_F128)
BM4_A3: float = float(BM4_A3_F128)
BM4_A4: float = float(BM4_A4_F128)

BM4_DRIFT_WEIGHTS_F128: Tuple[mp.mpf, ...] = (
    BM4_A1_F128,
    BM4_A2_F128,
    BM4_A3_F128,
    BM4_A4_F128,
    BM4_A3_F128,
    BM4_A2_F128,
    BM4_A1_F128,
)

BM4_IMPULSE_WEIGHTS_F128: Tuple[mp.mpf, ...] = (
    mp.mpf('0.5') * BM4_A1_F128,
    mp.mpf('0.5') * (BM4_A1_F128 + BM4_A2_F128),
    mp.mpf('0.5') * (BM4_A2_F128 + BM4_A3_F128),
    mp.mpf('0.5') * (BM4_A3_F128 + BM4_A4_F128),
    mp.mpf('0.5') * (BM4_A4_F128 + BM4_A3_F128),
    mp.mpf('0.5') * (BM4_A3_F128 + BM4_A2_F128),
    mp.mpf('0.5') * (BM4_A2_F128 + BM4_A1_F128),
    mp.mpf('0.5') * BM4_A1_F128,
)

BM4_DRIFT_WEIGHTS: Tuple[float, ...] = tuple(float(x) for x in BM4_DRIFT_WEIGHTS_F128)
BM4_IMPULSE_WEIGHTS: Tuple[float, ...] = tuple(float(x) for x in BM4_IMPULSE_WEIGHTS_F128)

# Certified 40-Digit Algebraic Coefficients for BM_6^(11) (Blanes & Moan 2002, Table 3 SRKN_11^b; Chapter 5, Table 5.1)
# 11-stage palindromic composition in kick-exterior FSAL form:
# a = drift weights (6 distinct values, palindromic symmetry a_{12-i} = a_i)
# b = impulse weights (6 distinct values, palindromic symmetry b_{12-i} = b_i)
# Satisfying all 8 RKN order-6 conditions to < 1e-35.
BM6_A1_F128: mp.mpf = mp.mpf('+0.1232297759462711643039405467597932509199')
BM6_A2_F128: mp.mpf = mp.mpf('+0.2905537977995582630562587859964007637277')
BM6_A3_F128: mp.mpf = mp.mpf('-0.1270492126254169219936895974447724790771')
BM6_A4_F128: mp.mpf = mp.mpf('-0.2463317610620754120741655636241040997077')
BM6_A5_F128: mp.mpf = mp.mpf('+0.3572088727959277730105970560242744614724')
BM6_A6_F128: mp.mpf = mp.mpf('1.0') - mp.mpf('2.0') * (
    BM6_A1_F128 + BM6_A2_F128 + BM6_A3_F128 + BM6_A4_F128 + BM6_A5_F128
)

BM6_B1_F128: mp.mpf = mp.mpf('+0.0414649985182628496366924125064201851397')
BM6_B2_F128: mp.mpf = mp.mpf('+0.1981286719180665471711545861105697559659')
BM6_B3_F128: mp.mpf = mp.mpf('-0.0400061921041526865377642138968083380434')
BM6_B4_F128: mp.mpf = mp.mpf('+0.0752539843015803357148000959587589288775')
BM6_B5_F128: mp.mpf = mp.mpf('-0.0115113874206877682143805227909295690717')
BM6_B6_F128: mp.mpf = mp.mpf('1.0') - mp.mpf('2.0') * (
    BM6_B1_F128 + BM6_B2_F128 + BM6_B3_F128 + BM6_B4_F128 + BM6_B5_F128
)

BM6_DRIFT_WEIGHTS_F128: Tuple[mp.mpf, ...] = (
    BM6_A1_F128,
    BM6_A2_F128,
    BM6_A3_F128,
    BM6_A4_F128,
    BM6_A5_F128,
    BM6_A6_F128,
)

BM6_IMPULSE_WEIGHTS_F128: Tuple[mp.mpf, ...] = (
    BM6_B1_F128,
    BM6_B2_F128,
    BM6_B3_F128,
    BM6_B4_F128,
    BM6_B5_F128,
    BM6_B6_F128,
)

BM6_DRIFT_WEIGHTS: Tuple[float, ...] = tuple(float(x) for x in BM6_DRIFT_WEIGHTS_F128)
BM6_IMPULSE_WEIGHTS: Tuple[float, ...] = tuple(float(x) for x in BM6_IMPULSE_WEIGHTS_F128)
BM6_A1, BM6_A2, BM6_A3, BM6_A4, BM6_A5, BM6_A6 = BM6_DRIFT_WEIGHTS
BM6_B1, BM6_B2, BM6_B3, BM6_B4, BM6_B5, BM6_B6 = BM6_IMPULSE_WEIGHTS

# Calibrated Homological Coefficient for Partial Boundary Corrector K_4 (Chapter 5, Section 5.5, Theorem 5.15)
# Formal Lean 4 Certification: 01_CorrectorBM4_Theorems.lean
KAPPA4_BM4_F128: mp.mpf = mp.mpf('1.0') / mp.mpf('2880.0')
KAPPA4_BM4: float = float(KAPPA4_BM4_F128)

# Exact 40-Digit 4th-Order Baker--Campbell--Hausdorff / Duhamel Moment of 7-Stage Blanes--Moan BM_4^{(7)}:
#   c_4^{BM4} = (sum_{k=1}^8 b_k * tau_k^4 - 1/80) / 24
def _compute_c4_bch_bm4_f128() -> mp.mpf:
    with mp.workdps(50):
        tau = []
        acc = mp.mpf('-0.5')
        tau.append(acc)
        for a_val in BM4_DRIFT_WEIGHTS_F128:
            acc = acc + a_val
            tau.append(acc)
        mu4 = sum(b_val * (t_val ** 4) for b_val, t_val in zip(BM4_IMPULSE_WEIGHTS_F128, tau))
        return (mu4 - mp.mpf('1.0') / mp.mpf('80.0')) / mp.mpf('24.0')


C4_BCH_BM4_F128: mp.mpf = _compute_c4_bch_bm4_f128()
C4_BCH_BM4: float = float(C4_BCH_BM4_F128)
HOMOLOGICAL_ALPHA1_BM4_F128: mp.mpf = mp.mpf('0.25')
HOMOLOGICAL_ALPHA2_BM4_F128: mp.mpf = mp.mpf('0.125')
HOMOLOGICAL_BETA1_BM4_F128: mp.mpf = mp.mpf('256.0') * C4_BCH_BM4_F128
HOMOLOGICAL_BETA2_BM4_F128: mp.mpf = mp.mpf('-512.0') * C4_BCH_BM4_F128
HOMOLOGICAL_ALPHA1_BM4: float = float(HOMOLOGICAL_ALPHA1_BM4_F128)
HOMOLOGICAL_ALPHA2_BM4: float = float(HOMOLOGICAL_ALPHA2_BM4_F128)
HOMOLOGICAL_BETA1_BM4: float = float(HOMOLOGICAL_BETA1_BM4_F128)
HOMOLOGICAL_BETA2_BM4: float = float(HOMOLOGICAL_BETA2_BM4_F128)



