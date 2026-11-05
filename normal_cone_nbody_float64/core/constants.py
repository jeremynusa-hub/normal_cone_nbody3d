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

from typing import Tuple
import numpy as np

# Normalized gravitational constant (canonical N-body units where G = 1.0)
G_DEFAULT: float = 1.0
GRAVITATIONAL_CONSTANT: float = G_DEFAULT

# Universal Hill/Roche stability threshold governing synchronous supernode collapse (Section 4.3)
# Derived from the collinear L_1 Lagrange point where radial tidal acceleration equals 2/3 of internal gravity:
# 2 * eta^3 = 2/3 <=> eta_crit = (1/3)^(1/3) = 3^(-1/3) ~= 0.6933612743506347
ETA_CRIT: float = float(3.0 ** (-1.0 / 3.0))
ETA_CRITICAL_HILL: float = ETA_CRIT

# Canonical Spectral Hysteresis Band for Hill Supernode Activation/Deactivation (Section 4.3, Cor. 4.1)
# Derived from the spectral decomposition of the symmetric trace-free Newtonian tidal tensor:
#   T(r_P) = I_3 - 3 * r_hat_P * r_hat_P^T  with eigenvalues (-2, +1, +1)
# 1. Activation (Anticipation) at exact radial equigravity (maximum tidal pull equals internal binding):
#    2 * eta^3 = 1 <=> eta_enter = 2^(-1/3) ~= 0.7937005259840998
# 2. Deactivation at spherical RMS equigravity (uniform mean tidal pull over S^2 equals internal binding):
#    sqrt(2) * eta^3 = 1 <=> eta_exit = 2^(-1/6) ~= 0.8908987181403393
# Note: Both thresholds are strictly above eta_crit (3^(-1/3) < 2^(-1/3) < 2^(-1/6)), guaranteeing
# proactive activation before the Hill boundary is reached and satisfying eta_enter = (eta_exit)^2.
ETA_ENTER_HILL: float = float(2.0 ** (-1.0 / 3.0))
ETA_EXIT_HILL: float = float(2.0 ** (-1.0 / 6.0))

# Cubic Tidal Force Ratio Thresholds for Roche Tidal Disruption (Section 7.2.1, Equation 7.34)
# Lambda_S = (1 / 3) * eta_S^3 = (M_ext / (3 * M_S)) * (bar_rho_S / r_P)^3
LAMBDA_CRITICAL_ROCHE: float = 0.25
LAMBDA_ENTER_ROCHE: float = 0.25
LAMBDA_EXIT_ROCHE: float = 0.18

# Canonical Voronoi admissibility ratio (gamma = 1.0, midpoint 0.5/0.5, zero hysteresis by default;
# optional gamma < 1.0 e.g. 0.60 supported for hysteresis experiments)
GAMMA_HYSTERESIS: float = 1.0

# Width of the dimensionless hysteresis gap: Delta_lambda = 5/3 - 3/5 = 16/15 ~= 1.066667
DELTA_LAMBDA_HYSTERESIS: float = 16.0 / 15.0

# IEEE 754 double precision machine epsilon floor
EPS_MACH: float = float(np.finfo(np.float64).eps)  # 2.220446049250313e-16

# Zero softening enforcement: gravitational potential is unsoftened (epsilon = 0.0, d_0 = 0.0)
ZERO_SOFTENING: float = 0.0

# Canonical scale-adapted screening transition window u_v in [1.0, 2.0] (Chapter 4, Section 4.2 & Chapter 6, Section 6.4.3)
SCREENING_U0: float = 1.0
SCREENING_U1: float = 2.0
SCREENING_SPAN: float = SCREENING_U1 - SCREENING_U0

# Scale-invariant maximum dyadic refinement depth (Chapter 5, Section 5.6)
# k_max = 17 permits 2^17 = 131,072 substeps, resolving extreme pericenters (r ~ 10^-4, v ~ 100 m/s)
# ensuring h_sub / tau_kin <= 0.05 without artificial clipping
DEFAULT_MAX_DYADIC_LEVEL: int = 17

# Certified 4th-Order 7-Stage Blanes--Moan Coefficients for BM_4^{(7)} (Chapter 5, Section 5.3)
# Satisfying 2*(a_1 + a_2 + a_3) + a_4 = 1.0 and 2*(a_1^3 + a_2^3 + a_3^3) + a_4^3 = 0.0 to exact IEEE 754 precision.
BM4_A1: float = +0.2390843168869614
BM4_A2: float = +0.2390725069195032
BM4_A3: float = +0.23905247436220158
BM4_A4: float = 1.0 - 2.0 * (BM4_A1 + BM4_A2 + BM4_A3)  # -0.4344185963373323

BM4_DRIFT_WEIGHTS: Tuple[float, ...] = (
    BM4_A1,
    BM4_A2,
    BM4_A3,
    BM4_A4,
    BM4_A3,
    BM4_A2,
    BM4_A1,
)

BM4_IMPULSE_WEIGHTS: Tuple[float, ...] = (
    0.5 * BM4_A1,
    0.5 * (BM4_A1 + BM4_A2),
    0.5 * (BM4_A2 + BM4_A3),
    0.5 * (BM4_A3 + BM4_A4),
    0.5 * (BM4_A4 + BM4_A3),
    0.5 * (BM4_A3 + BM4_A2),
    0.5 * (BM4_A2 + BM4_A1),
    0.5 * BM4_A1,
)

# Certified 40-Digit Algebraic Coefficients for BM_6^(11) (Blanes & Moan 2002, Table 3 SRKN_11^b; Chapter 5, Table 1)
# 11-stage palindromic composition in kick-exterior FSAL form:
# a = drift weights (6 distinct values, palindromic symmetry a_{12-i} = a_i)
# b = impulse weights (6 distinct values, palindromic symmetry b_{12-i} = b_i)
BM6_A1: float = +0.1232297759462711643039405467597932509199
BM6_A2: float = +0.2905537977995582630562587859964007637277
BM6_A3: float = -0.1270492126254169219936895974447724790771
BM6_A4: float = -0.2463317610620754120741655636241040997077
BM6_A5: float = +0.3572088727959277730105970560242744614724
BM6_A6: float = +0.2047770542914702673941175445768161973296

BM6_B1: float = +0.0414649985182628496366924125064201851397
BM6_B2: float = +0.1981286719180665471711545861105697559659
BM6_B3: float = -0.0400061921041526865377642138968083380434
BM6_B4: float = +0.0752539843015803357148000959587589288775
BM6_B5: float = -0.0115113874206877682143805227909295690717
BM6_B6: float = +0.4733398495738614444590002842239780742640

BM6_DRIFT_WEIGHTS: Tuple[float, ...] = (
    BM6_A1,
    BM6_A2,
    BM6_A3,
    BM6_A4,
    BM6_A5,
    BM6_A6,
)

BM6_IMPULSE_WEIGHTS: Tuple[float, ...] = (
    BM6_B1,
    BM6_B2,
    BM6_B3,
    BM6_B4,
    BM6_B5,
    BM6_B6,
)

# Calibrated Homological Coefficient for Partial Boundary Corrector K_4 (Chapter 5, Section 5.5, Theorem 5.15)
# Formal Lean 4 Certification: 01_CorrectorBM4_Theorems.lean
KAPPA4_BM4: float = 1.0 / 2880.0

# Exact 4th-Order Baker--Campbell--Hausdorff / Duhamel Moment of 7-Stage Blanes--Moan BM_4^{(7)}:
#   c_4^{BM4} = (sum_{k=1}^8 b_k * tau_k^4 - 1/80) / 24 = +1.7600050100410129e-4 (~= 1 / 5681.802008026074)
# where tau_k = sum_{j < k} a_j - 1/2 are the centered impulse abscissae.
C4_BCH_BM4: float = +1.76000501004101292305e-4
HOMOLOGICAL_ALPHA1_BM4: float = 0.25
HOMOLOGICAL_ALPHA2_BM4: float = 0.125
HOMOLOGICAL_BETA1_BM4: float = 256.0 * C4_BCH_BM4
HOMOLOGICAL_BETA2_BM4: float = -512.0 * C4_BCH_BM4



