r"""
Universal Kepler-Sundman Two-Body Solver with Cubic Halley Convergence.

Theoretical Foundations:
- Monograph Reference: Chapter 5, Section 5.2.2 (Sundman--Stumpff Analytical Diffeomorphism
  under Retrograde Steps, Proposition 5.2.2, Equations 5.17-5.23).
- Sundman Regularized Time Equation:
    dt = r(tau) dtau = (rho_0 c_0(z) + sigma_0 tau c_1(z) + GM tau^2 c_2(z)) dtau.
    F(tau) = rho_0 tau c_1(z) + sigma_0 tau^2 c_2(z) + GM tau^3 c_3(z) - dt = 0,
    where z = alpha tau^2 and alpha = 2 GM / rho_0 - ||v_0||^2.
- Exact Halley Iteration Derivatives:
    F'(tau) = rho(tau) = rho_0 c_0(z) + sigma_0 tau c_1(z) + GM tau^2 c_2(z),
    F''(tau) = sigma_0 c_0(z) + (GM - alpha rho_0) tau c_1(z).
- Wronskian Unimodularity (Darboux Symplectic Invariance):
    f dot{g} - dot{f} g == +1.0 (exact to machine precision IEEE 754 < 1e-15).
- Certified Formal Proofs:
    Lean 4 modules `01_YoshidaOrderFourSuzuki_Theorems.lean` and `PaperUnifiedFoundations.lean`.
"""

from __future__ import annotations
import math
from typing import Tuple, Optional, Union, Any
import numpy as np
import mpmath as mp

from .stumpff import stumpff_c, stumpff_c_f128
from ..core.profiler import get_current_profiler

HALLEY_TOL_F128: mp.mpf = mp.mpf('1e-32')
HALLEY_TOL_F64: float = 1e-15
BATTIN_PROJ_TOL_F128: mp.mpf = mp.mpf('1e-32')
BATTIN_PROJ_TOL_F64: float = 1e-12

# ==============================================================================
# IEEE 754 Hardware Hyperbolic Guard Constants (DISC-NUM-03)
# ==============================================================================
# In float64 (IEEE 754 binary64), sinh(x) overflows for x > ~709.7827 = ln(2^1024),
# so -z <= (709.78)^2 ~= 503,791. Z_MIN_HYPERBOLIC = -500_000.0 is the exact safe guard.
# In float128 / quad precision (IEEE 754-2008 binary128), exponent reaches 16383,
# sinh(x) overflows for x > 16384 * ln(2) ~= 11356.52.
# So -z <= (11356.52)^2 ~= 1.2897e+8. Z_MIN_HYPERBOLIC_F128 = mp.mpf('-1.28e8')
# provides the authentic machine hardware bound for quadruple precision.
Z_MIN_HYPERBOLIC = -500_000.0
Z_MIN_HYPERBOLIC_F128: mp.mpf = mp.mpf('-1.28e8')


def solve_sundman_time_equation_f128(
    r0: Union[mp.mpf, float],
    sigma0: Union[mp.mpf, float],
    GM: Union[mp.mpf, float],
    alpha: Union[mp.mpf, float],
    dt_target: Union[mp.mpf, float],
    tol: Union[mp.mpf, float] = HALLEY_TOL_F128,
    max_iter: int = 40,
) -> mp.mpf:
    """
    Numerically inverts the universal Sundman time equation F(tau) = dt_target
    via cubic Halley iteration in 35-digit mpmath precision (tol = 1e-32).
    """
    r0 = mp.mpf(r0)
    sigma0 = mp.mpf(sigma0)
    GM = mp.mpf(GM)
    alpha = mp.mpf(alpha)
    dt_target = mp.mpf(dt_target)
    tol = mp.mpf(tol)

    if dt_target == mp.mpf('0.0'):
        return mp.mpf('0.0')

    sgn_dt = mp.mpf('1.0') if dt_target > mp.mpf('0.0') else mp.mpf('-1.0')
    abs_dt = abs(dt_target)
    eff_sigma0 = sgn_dt * sigma0
    abs_tol = tol * abs_dt

    # 1. Physics-based ballistic and collision initial guess
    v0_sq = max((mp.mpf('2.0') * GM / r0) - alpha, mp.mpf('0.0'))
    r_collision_scale = (GM * (abs_dt ** 2)) ** (mp.mpf('1.0') / mp.mpf('3.0'))

    if r0 < mp.mpf('0.1') * r_collision_scale and v0_sq * (abs_dt ** 2) < (r_collision_scale ** 2):
        tau = (mp.mpf('6.0') * abs_dt / GM) ** (mp.mpf('1.0') / mp.mpf('3.0'))
    elif alpha < mp.mpf('0.0'):
        sqrt_neg_alpha = mp.sqrt(-alpha)
        arg = max(sqrt_neg_alpha * abs_dt / max(r0, mp.mpf('1e-60')), mp.mpf('0.0'))
        tau = mp.asinh(arg) / sqrt_neg_alpha
        tau_max_hyp = mp.sqrt(mp.fabs(Z_MIN_HYPERBOLIC_F128) / mp.fabs(alpha))
        tau = min(tau, tau_max_hyp)
    else:
        r_floor_sq = max(mp.mpf('1e-48'), mp.mpf('0.01') * (r_collision_scale ** 2))
        r_mid_sq = max(r0 * r0 + eff_sigma0 * abs_dt + mp.mpf('0.25') * v0_sq * (abs_dt ** 2), r_floor_sq)
        r_end_sq = max(r0 * r0 + mp.mpf('2.0') * eff_sigma0 * abs_dt + v0_sq * (abs_dt ** 2), r_floor_sq)
        r_mid = mp.sqrt(r_mid_sq)
        r_end = mp.sqrt(r_end_sq)
        tau = (abs_dt / mp.mpf('6.0')) * (mp.mpf('1.0') / max(r0, mp.sqrt(r_floor_sq)) + mp.mpf('4.0') / r_mid + mp.mpf('1.0') / r_end)
        if alpha > mp.mpf('0.0'):
            tau_ell_max = (alpha * abs_dt / max(GM, mp.mpf('1e-60')) + mp.pi) / mp.sqrt(alpha)
            tau = min(tau, tau_ell_max)

    tau_floor = min(mp.mpf('1e-32'), mp.mpf('1e-32') * abs_dt / max(r0, mp.mpf('1.0')))
    tau = max(tau, tau_floor)

    def _eval_sundman(t: mp.mpf) -> Tuple[mp.mpf, mp.mpf, mp.mpf]:
        zt = max(alpha * t * t, Z_MIN_HYPERBOLIC_F128)
        c0, c1, c2, c3, _, _ = stumpff_c_f128(zt)
        t_sq = t * t
        t_cb = t_sq * t
        f = r0 * t * c1 + eff_sigma0 * t_sq * c2 + GM * t_cb * c3 - abs_dt
        f1 = r0 * c0 + eff_sigma0 * t * c1 + GM * t_sq * c2
        f2 = eff_sigma0 * c0 + (GM - alpha * r0) * t * c1
        return f, f1, f2

    # 2. Cubic Halley iteration with trust-region step bounding
    converged = False
    halley_iters = 0
    for _ in range(max_iter):
        halley_iters += 1
        f_val, f_prime, f_prime2 = _eval_sundman(tau)
        if abs(f_val) <= abs_tol:
            converged = True
            break

        f_prime = max(f_prime, mp.mpf('1e-60'))
        denom = mp.mpf('2.0') * (f_prime ** 2) - f_val * f_prime2

        if denom <= mp.mpf('0.1') * (f_prime ** 2):
            step = f_val / f_prime
        else:
            step = (mp.mpf('2.0') * f_prime * f_val) / denom

        if abs(step) > mp.mpf('0.5') * tau:
            step = mp.sign(step) * min(abs(step), mp.mpf('0.5') * tau)

        tau -= step
        if alpha < mp.mpf('0.0'):
            tau = min(tau, tau_max_hyp)
        tau = max(tau, tau_floor)

    # 3. Guaranteed monotone bracket fallback (Deterministic bisection, 115 iterations, 10^-35 resolution)
    used_bisection = False
    if not converged:
        used_bisection = True
        a = mp.mpf('0.0')
        b = max(tau, abs_dt / r0, mp.mpf('1.0'))
        if alpha < mp.mpf('0.0'):
            tau_max_hyp = mp.sqrt(mp.fabs(Z_MIN_HYPERBOLIC_F128) / mp.fabs(alpha))
            b = min(b, tau_max_hyp)

        fb = _eval_sundman(b)[0]
        while fb < mp.mpf('0.0') and b < mp.mpf('1e16'):
            b *= mp.mpf('2.0')
            if alpha < mp.mpf('0.0'):
                tau_max_hyp = mp.sqrt(mp.fabs(Z_MIN_HYPERBOLIC_F128) / mp.fabs(alpha))
                if b > tau_max_hyp:
                    b = tau_max_hyp
                    break
            fb = _eval_sundman(b)[0]

        for _ in range(115):
            mid = mp.mpf('0.5') * (a + b)
            f_mid, _, _ = _eval_sundman(mid)
            if abs(f_mid) <= abs_tol or (b - a) < mp.mpf('1e-34') * mid:
                tau = mid
                converged = True
                break
            if f_mid < mp.mpf('0.0'):
                a = mid
            else:
                b = mid
        tau = mid

    prof = get_current_profiler()
    if prof is not None:
        prof.record_kepler(
            halley_iters=halley_iters,
            stumpff_evals=halley_iters + (115 if used_bisection else 0),
            is_bisection=used_bisection,
        )

    return sgn_dt * tau


def solve_sundman_time_equation(
    r0: float,
    sigma0: float,
    GM: float,
    alpha: float,
    dt_target: float,
    tol: float = 1e-15,
    max_iter: int = 25,
) -> float:
    """
    Numerically inverts the universal Sundman time equation F(tau) = dt_target via cubic Halley iteration.

    Parameters:
        r0: Initial separation distance ||r_0|| > 0.
        sigma0: Initial radial velocity inner product r_0 . v_0.
        GM: Standard gravitational parameter G (m_1 + m_2).
        alpha: Energy parameter 2 GM / r_0 - ||v_0||^2.
        dt_target: Physical time step to cross (can be negative for retrograde stages).
        tol: Relative convergence tolerance.
        max_iter: Maximum number of Halley iterations (cubic convergence typically reaches tolerance in <= 4 steps).

    Returns:
        tau: Fictitious Sundman time parameter satisfying F(tau) = dt_target.
    """
    if any(isinstance(x, mp.mpf) for x in (r0, sigma0, GM, alpha, dt_target, tol)) or tol < 1e-16:
        return float(solve_sundman_time_equation_f128(r0, sigma0, GM, alpha, dt_target, tol=tol, max_iter=max_iter))

    if dt_target == 0.0:
        return 0.0

    # ==========================================================================
    # Time-Reversal Symmetry and Retrograde Step Factorization (DISC-NUM-02)
    # ==========================================================================
    # High-order symplectic partitioned Runge-Kutta-Nyström schemes (such as Yoshida S4
    # and Blanes-Moan BM6) inherently require negative intermediate step weights
    # (a_i < 0, b_i < 0) to bypass the Suzuki (1991) order barrier.
    # To guarantee exact microscopic time-reversibility S(-dt) o S(dt) = Id to
    # IEEE 754 machine precision (< 1e-15) under retrograde steps (dt_target < 0),
    # the Sundman time equation F(tau) = dt is solved via momentum involution:
    #   T : (r, w) |-> (r, -w),  yielding sigma0_eff = sgn(dt) * sigma0.
    # This maps the backward trajectory onto a forward motion under inverted radial
    # momentum, ensuring identical Halley convergence basins and preventing any
    # secular phase drift or spurious asymmetry across negative substeps.
    sgn_dt = 1.0 if dt_target > 0.0 else -1.0
    abs_dt = abs(dt_target)
    eff_sigma0 = sgn_dt * sigma0
    abs_tol = tol * abs_dt

    # 1. Physics-based ballistic and collision initial guess
    v0_sq = max((2.0 * GM / r0) - alpha, 0.0)
    r_collision_scale = (GM * (abs_dt**2)) ** (1.0 / 3.0)

    if r0 < 0.1 * r_collision_scale and v0_sq * (abs_dt**2) < (r_collision_scale**2):
        # Gravitationally dominated free-fall from near rest
        tau = (6.0 * abs_dt / GM) ** (1.0 / 3.0)
    elif alpha < 0.0:
        # High-energy hyperbolic regime: dt ~ (r0 / sqrt(-alpha)) sinh(sqrt(-alpha) tau)
        sqrt_neg_alpha = math.sqrt(-alpha)
        arg = max(sqrt_neg_alpha * abs_dt / max(r0, 1e-30), 0.0)
        tau = math.asinh(arg) / sqrt_neg_alpha
    else:
        # High-velocity / rectilinear ballistic Simpson quadrature for elliptic/parabolic
        r_floor_sq = max(1e-24, 0.01 * (r_collision_scale**2))
        r_mid_sq = max(r0 * r0 + eff_sigma0 * abs_dt + 0.25 * v0_sq * (abs_dt**2), r_floor_sq)
        r_end_sq = max(r0 * r0 + 2.0 * eff_sigma0 * abs_dt + v0_sq * (abs_dt**2), r_floor_sq)
        r_mid = math.sqrt(r_mid_sq)
        r_end = math.sqrt(r_end_sq)
        tau = (abs_dt / 6.0) * (1.0 / max(r0, math.sqrt(r_floor_sq)) + 4.0 / r_mid + 1.0 / r_end)
        if alpha > 0.0:
            tau_ell_max = (alpha * abs_dt / max(GM, 1e-30) + math.pi) / math.sqrt(alpha)
            tau = min(tau, tau_ell_max)

    if alpha < 0.0:
        tau_max_hyp = math.sqrt(abs(Z_MIN_HYPERBOLIC) / abs(alpha))
        tau = min(tau, tau_max_hyp)

    tau_floor = min(1e-15, 1e-15 * abs_dt / max(r0, 1.0))
    tau = max(tau, tau_floor)

    def _eval_sundman(t: float) -> Tuple[float, float, float]:
        zt = max(alpha * t * t, Z_MIN_HYPERBOLIC)
        c0, c1, c2, c3, _, _ = stumpff_c(zt)
        t_sq = t * t
        t_cb = t_sq * t
        f = r0 * t * c1 + eff_sigma0 * t_sq * c2 + GM * t_cb * c3 - abs_dt
        f1 = r0 * c0 + eff_sigma0 * t * c1 + GM * t_sq * c2
        f2 = eff_sigma0 * c0 + (GM - alpha * r0) * t * c1
        return f, f1, f2

    # 2. Cubic Halley iteration with trust-region step bounding
    converged = False
    halley_iters = 0
    for _ in range(max_iter):
        halley_iters += 1
        f_val, f_prime, f_prime2 = _eval_sundman(tau)
        if abs(f_val) <= abs_tol:
            converged = True
            break

        f_prime = max(f_prime, 1e-30)
        denom = 2.0 * (f_prime**2) - f_val * f_prime2

        # Trust region: if denominator degenerates or approaches zero, fall back to damped Newton step
        if denom <= 0.1 * (f_prime**2):
            step = f_val / f_prime
        else:
            step = (2.0 * f_prime * f_val) / denom

        # Guard against wild leaps overshooting the convergence basin
        if not np.isfinite(step) or abs(step) > 0.5 * tau:
            step = math.copysign(min(abs(step), 0.5 * tau), step)

        tau -= step
        if alpha < 0.0:
            tau = min(tau, tau_max_hyp)
        tau = max(tau, tau_floor)

    # 3. Guaranteed monotone bracket fallback (Deterministic bisection, 60 iterations, 10^-18 resolution)
    # Strictly valid because F'(tau) = r(tau) > 0 universally (F is strictly increasing)
    used_bisection = False
    if not converged:
        used_bisection = True
        a = 0.0
        b = max(tau, abs_dt / r0, 1.0)
        if alpha < 0.0:
            tau_max_hyp = math.sqrt(500_000.0 / abs(alpha))
            b = min(b, tau_max_hyp)

        fb = _eval_sundman(b)[0]
        while fb < 0.0 and b < 1e8:
            b *= 2.0
            if alpha < 0.0:
                tau_max_hyp = math.sqrt(500_000.0 / abs(alpha))
                if b >= tau_max_hyp:
                    b = tau_max_hyp
                    fb = _eval_sundman(b)[0]
                    break
            fb = _eval_sundman(b)[0]

        for _ in range(60):
            mid = 0.5 * (a + b)
            f_mid, _, _ = _eval_sundman(mid)
            if abs(f_mid) <= abs_tol or (b - a) < 1e-15 * mid:
                tau = mid
                converged = True
                break
            if f_mid < 0.0:
                a = mid
            else:
                b = mid
        tau = mid

    prof = get_current_profiler()
    if prof is not None:
        prof.record_kepler(
            halley_iters=halley_iters,
            stumpff_evals=halley_iters + (60 if used_bisection else 0),
            is_bisection=used_bisection,
        )

    return sgn_dt * tau


def step_kepler_sundman_f128(
    r0: np.ndarray,
    w0: np.ndarray,
    mu: Union[mp.mpf, float],
    GM: Union[mp.mpf, float],
    dt_target: Union[mp.mpf, float],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Analytically integrates an unsoftened two-body Keplerian channel across physical step dt_target
    in 35-digit mpmath precision.
    """
    r0_arr = np.array([mp.mpf(x) for x in r0], dtype=object)
    w0_arr = np.array([mp.mpf(x) for x in w0], dtype=object)
    mu_mp = mp.mpf(mu)
    GM_mp = mp.mpf(GM)
    dt_mp = mp.mpf(dt_target)

    if dt_mp == mp.mpf('0.0'):
        return r0_arr.copy(), w0_arr.copy()

    if dt_mp < mp.mpf('0.0'):
        r_rev, w_rev = step_kepler_sundman_f128(r0_arr, -w0_arr, mu_mp, GM_mp, -dt_mp)
        return r_rev, -w_rev

    v0 = w0_arr / mu_mp
    rx, ry, rz = r0_arr[0], r0_arr[1], r0_arr[2]
    vx, vy, vz = v0[0], v0[1], v0[2]

    rho0_sq = rx * rx + ry * ry + rz * rz
    rho0 = mp.sqrt(rho0_sq)
    sigma0 = rx * vx + ry * vy + rz * vz
    v0_sq = vx * vx + vy * vy + vz * vz

    # Specific orbital angular momentum L = r0 x v0
    lx = ry * vz - rz * vy
    ly = rz * vx - rx * vz
    lz = rx * vy - ry * vx
    L_sq = lx * lx + ly * ly + lz * lz

    if rho0 < mp.mpf('1e-32'):
        v_norm = mp.sqrt(v0_sq)
        if v_norm > mp.mpf('1e-32'):
            inv_v = mp.mpf('1.0') / v_norm
            s_out = np.array([vx * inv_v, vy * inv_v, vz * inv_v], dtype=object)
        else:
            s_out = np.array([mp.mpf('1.0'), mp.mpf('0.0'), mp.mpf('0.0')], dtype=object)

        tau = (mp.mpf('6.0') * dt_mp / GM_mp) ** (mp.mpf('1.0') / mp.mpf('3.0'))
        c0, c1, c2, c3, _, _ = stumpff_c_f128(mp.mpf('0.0'))
        rho_final = GM_mp * (tau ** 2) * c2
        v_phys = c1 / (c2 * tau)
        r_final = rho_final * s_out
        w_final = (mu_mp * v_phys) * s_out
        return r_final, w_final

    # Energy parameter alpha = 2 GM / rho0 - v0^2
    alpha = (mp.mpf('2.0') * GM_mp / rho0) - v0_sq

    # 1. Invert Sundman clock via cubic Halley iteration
    tau = solve_sundman_time_equation_f128(rho0, sigma0, GM_mp, alpha, dt_mp)

    # 2. Universal Stumpff functions
    z = alpha * tau * tau
    c0, c1, c2, c3, _, _ = stumpff_c_f128(z)

    tau_sq = tau * tau
    rho_raw = rho0 * c0 + sigma0 * tau * c1 + GM_mp * tau_sq * c2
    sigma_end = sigma0 * c0 + (GM_mp - alpha * rho0) * tau * c1

    # 3. Head-on radial collision handling (L_sq ~ 0, e = 1):
    # Pure analytical Levi-Civita / Sundman-Stumpff derivative v_radial = F''(tau) / F'(tau)
    scale_L2 = GM_mp * rho0 + v0_sq * rho0_sq
    is_head_on = bool(L_sq < mp.mpf('1e-48') * max(scale_L2, mp.mpf('1e-60')))
    if is_head_on:
        inv_rho0 = mp.mpf('1.0') / rho0
        s0_x = rx * inv_rho0
        s0_y = ry * inv_rho0
        s0_z = rz * inv_rho0

        denom_peri = mp.mpf('2.0') * GM_mp - alpha * rho_raw
        if abs(rho_raw) < mp.mpf('0.25') * rho0 and denom_peri > mp.mpf('0.5') * GM_mp:
            rho_final_abs = (sigma_end * sigma_end + L_sq) / denom_peri
            rho_final_abs = (sigma_end * sigma_end + L_sq) / (mp.mpf('2.0') * GM_mp - alpha * rho_final_abs)
        else:
            rho_final_abs = abs(rho_raw)
        if rho_final_abs > mp.mpf('1e-32'):
            v_radial = sigma_end / rho_final_abs
        else:
            v_radial = (c1 / (c2 * tau)) if abs(c2 * tau) > mp.mpf('1e-60') else mp.mpf('0.0')

        r_final = np.array([rho_final_abs * s0_x, rho_final_abs * s0_y, rho_final_abs * s0_z], dtype=object)
        w_final = np.array([mu_mp * v_radial * s0_x, mu_mp * v_radial * s0_y, mu_mp * v_radial * s0_z], dtype=object)
        return r_final, w_final

    # 4. Regular 3D orbit (L > 0)
    denom_peri = mp.mpf('2.0') * GM_mp - alpha * rho_raw
    if rho_raw < mp.mpf('0.25') * rho0 and denom_peri > mp.mpf('0.5') * GM_mp and L_sq > mp.mpf('0.0'):
        rho_final = (sigma_end * sigma_end + L_sq) / denom_peri
        rho_final = (sigma_end * sigma_end + L_sq) / (mp.mpf('2.0') * GM_mp - alpha * rho_final)
        f_prime = - (GM_mp / rho0) * tau * c1
        g_prime = rho0 * c0 + sigma0 * tau * c1
        rpx = f_prime * rx + g_prime * vx
        rpy = f_prime * ry + g_prime * vy
        rpz = f_prime * rz + g_prime * vz
        # L x r'
        lcr_x = ly * rpz - lz * rpy
        lcr_y = lz * rpx - lx * rpz
        lcr_z = lx * rpy - ly * rpx
        sx = sigma_end * rpx - lcr_x
        sy = sigma_end * rpy - lcr_y
        sz_v = sigma_end * rpz - lcr_z
        s_norm = mp.sqrt(sx * sx + sy * sy + sz_v * sz_v)
        if s_norm > mp.mpf('0.0'):
            shx, shy, shz = sx / s_norm, sy / s_norm, sz_v / s_norm
            r_final = np.array([rho_final * shx, rho_final * shy, rho_final * shz], dtype=object)
            ux = ly * shz - lz * shy
            uy = lz * shx - lx * shz
            uz = lx * shy - ly * shx
            v_final = np.array([
                (sigma_end * shx + ux) / rho_final,
                (sigma_end * shy + uy) / rho_final,
                (sigma_end * shz + uz) / rho_final,
            ], dtype=object)
            return r_final, mu_mp * v_final
    else:
        rho_final = rho_raw

    # Incremental Gauss-Stumpff formulation in 35-digit mpmath precision
    df = - (GM_mp / rho0) * tau_sq * c2
    f_coeff = mp.mpf('1.0') + df
    g_coeff = rho0 * tau * c1 + sigma0 * tau_sq * c2
    f_dot = - (GM_mp / (rho0 * rho_final)) * tau * c1
    dg_dot = - (GM_mp / rho_final) * tau_sq * c2

    if abs(f_coeff) >= mp.mpf('0.5'):
        dg_dot = (f_dot * g_coeff - df) / f_coeff
        g_dot = mp.mpf('1.0') + dg_dot
    else:
        g_dot = mp.mpf('1.0') + dg_dot
        max_diag = max(abs(f_coeff), abs(g_dot))
        max_off = max(abs(g_coeff), abs(f_dot))
        if max_diag >= max_off:
            if abs(f_coeff) >= abs(g_dot) and abs(f_coeff) > mp.mpf('1e-24'):
                dg_dot = (f_dot * g_coeff - df) / f_coeff
                g_dot = mp.mpf('1.0') + dg_dot
            elif abs(g_dot) > mp.mpf('1e-24'):
                df = (f_dot * g_coeff - dg_dot) / g_dot
                f_coeff = mp.mpf('1.0') + df
        else:
            if abs(g_coeff) >= abs(f_dot) and abs(g_coeff) > mp.mpf('1e-24'):
                f_dot = (df + dg_dot + df * dg_dot) / g_coeff
            elif abs(f_dot) > mp.mpf('1e-24'):
                g_coeff = (df + dg_dot + df * dg_dot) / f_dot

    r_final = np.array([
        rx + (df * rx + g_coeff * vx),
        ry + (df * ry + g_coeff * vy),
        rz + (df * rz + g_coeff * vz),
    ], dtype=object)

    mu_f_dot = mu_mp * f_dot
    w_final = np.array([
        w0_arr[0] + (mu_f_dot * rx + dg_dot * w0_arr[0]),
        w0_arr[1] + (mu_f_dot * ry + dg_dot * w0_arr[1]),
        w0_arr[2] + (mu_f_dot * rz + dg_dot * w0_arr[2]),
    ], dtype=object)
    return r_final, w_final


def step_kepler_sundman(
    r0: np.ndarray,
    w0: np.ndarray,
    mu: Union[float, mp.mpf],
    GM: Union[float, mp.mpf],
    dt_target: Union[float, mp.mpf],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Analytically integrates an unsoftened two-body Keplerian channel across physical step dt_target.
    Dispatches transparently to 35-digit mpmath precision when inputs are extended precision types.

    Parameters:
        r0: Relative separation vector r_0 in R^3.
        w0: Relative conjugate Jacobi momentum w_0 in R^3 (where v_0 = w_0 / mu).
        mu: Reduced mass m_1 m_2 / (m_1 + m_2).
        GM: Standard gravitational parameter G (m_1 + m_2).
        dt_target: Physical time step (positive or negative).

    Returns:
        (r_final, w_final): Propagated relative separation and conjugate momentum.
    """
    if (r0.dtype == object or w0.dtype == object or
        isinstance(mu, mp.mpf) or isinstance(GM, mp.mpf) or isinstance(dt_target, mp.mpf) or
        (len(r0) > 0 and isinstance(r0.flat[0], mp.mpf)) or
        (len(w0) > 0 and isinstance(w0.flat[0], mp.mpf))):
        return step_kepler_sundman_f128(r0, w0, mu, GM, dt_target)
    if dt_target == 0.0:
        return np.copy(r0), np.copy(w0)

    # Time-reversal symmetry handling for retrograde steps
    if dt_target < 0.0:
        r_rev, w_rev = step_kepler_sundman(r0, -w0, mu, GM, -dt_target)
        return r_rev, -w_rev

    v0 = w0 / mu
    rx, ry, rz = float(r0[0]), float(r0[1]), float(r0[2])
    vx, vy, vz = float(v0[0]), float(v0[1]), float(v0[2])

    rho0_sq = rx * rx + ry * ry + rz * rz
    rho0 = math.sqrt(rho0_sq)
    sigma0 = rx * vx + ry * vy + rz * vz
    v0_sq = vx * vx + vy * vy + vz * vz

    # Specific orbital angular momentum L = r0 x v0
    lx = ry * vz - rz * vy
    ly = rz * vx - rx * vz
    lz = rx * vy - ry * vx
    L_sq = lx * lx + ly * ly + lz * lz

    if rho0 < 1e-15:
        # Radial analytical continuation from exact collision singularity (Sundman 1912)
        v_norm = math.sqrt(v0_sq)
        if v_norm > 1e-15:
            inv_v = 1.0 / v_norm
            s_out = np.array([vx * inv_v, vy * inv_v, vz * inv_v], dtype=np.float64)
        else:
            s_out = np.array([1.0, 0.0, 0.0], dtype=np.float64)

        tau = (6.0 * dt_target / GM) ** (1.0 / 3.0)
        c0, c1, c2, c3, _, _ = stumpff_c(0.0)
        rho_final = GM * (tau**2) * c2
        v_phys = c1 / (c2 * tau)
        r_final = rho_final * s_out
        w_final = (mu * v_phys) * s_out
        return r_final, w_final

    # Energy parameter alpha = 2 GM / rho0 - v0^2
    alpha = (2.0 * GM / rho0) - v0_sq

    # 1. Invert Sundman clock via cubic Halley iteration
    tau = solve_sundman_time_equation(rho0, sigma0, GM, alpha, dt_target)

    # 2. Universal Stumpff functions
    z = alpha * tau * tau
    if z < Z_MIN_HYPERBOLIC:
        z = Z_MIN_HYPERBOLIC
    c0, c1, c2, c3, _, _ = stumpff_c(z)

    tau_sq = tau * tau
    rho_raw = rho0 * c0 + sigma0 * tau * c1 + GM * tau_sq * c2
    sigma_end = sigma0 * c0 + (GM - alpha * rho0) * tau * c1

    # 3. Head-on radial collision handling (L_sq ~ 0, e = 1)
    scale_L2 = GM * rho0 + v0_sq * rho0_sq
    is_head_on = bool(L_sq < 1e-24 * max(scale_L2, 1e-30))
    if is_head_on:
        inv_rho0 = 1.0 / rho0
        s0_x = rx * inv_rho0
        s0_y = ry * inv_rho0
        s0_z = rz * inv_rho0

        denom_peri = 2.0 * GM - alpha * rho_raw
        if abs(rho_raw) < 0.25 * rho0 and denom_peri > 0.5 * GM:
            rho_final_abs = (sigma_end * sigma_end + L_sq) / denom_peri
            rho_final_abs = (sigma_end * sigma_end + L_sq) / (2.0 * GM - alpha * rho_final_abs)
        else:
            rho_final_abs = abs(rho_raw)
        if rho_final_abs > 1e-15:
            v_radial = sigma_end / rho_final_abs
        else:
            v_radial = (c1 / (c2 * tau)) if abs(c2 * tau) > 1e-30 else 0.0

        r_final = np.array([rho_final_abs * s0_x, rho_final_abs * s0_y, rho_final_abs * s0_z], dtype=np.float64)
        w_final = np.array([mu * v_radial * s0_x, mu * v_radial * s0_y, mu * v_radial * s0_z], dtype=np.float64)
        return r_final, w_final

    # 4. Regular 3D orbit (L > 0)
    denom_peri = 2.0 * GM - alpha * rho_raw
    if rho_raw < 0.25 * rho0 and denom_peri > 0.5 * GM and L_sq > 0.0:
        rho_final = (sigma_end * sigma_end + L_sq) / denom_peri
        rho_final = (sigma_end * sigma_end + L_sq) / (2.0 * GM - alpha * rho_final)
        f_prime = - (GM / rho0) * tau * c1
        g_prime = rho0 * c0 + sigma0 * tau * c1
        r_prime = np.array([f_prime * rx + g_prime * vx, f_prime * ry + g_prime * vy, f_prime * rz + g_prime * vz], dtype=np.float64)
        L_vec = np.array([lx, ly, lz], dtype=np.float64)
        s_unnorm = sigma_end * r_prime - np.cross(L_vec, r_prime)
        s_norm = float(np.linalg.norm(s_unnorm))
        if s_norm > 0.0:
            s_hat = s_unnorm / s_norm
            r_final = rho_final * s_hat
            r_actual = float(np.linalg.norm(r_final))
            u_vec = np.cross(L_vec, s_hat)
            v_final = (sigma_end * s_hat + u_vec) / r_actual
            return r_final, mu * v_final
    else:
        rho_final = rho_raw

    df = - (GM / rho0) * tau_sq * c2
    f_coeff = 1.0 + df
    g_coeff = rho0 * tau * c1 + sigma0 * tau_sq * c2
    f_dot = - (GM / (rho0 * rho_final)) * tau * c1
    dg_dot = - (GM / rho_final) * tau_sq * c2

    if abs(f_coeff) >= 0.5:
        dg_dot = (f_dot * g_coeff - df) / f_coeff
        g_dot = 1.0 + dg_dot
    else:
        g_dot = 1.0 + dg_dot
        max_diag = max(abs(f_coeff), abs(g_dot))
        max_off = max(abs(g_coeff), abs(f_dot))
        if max_diag >= max_off:
            if abs(f_coeff) >= abs(g_dot) and abs(f_coeff) > 1e-12:
                dg_dot = (f_dot * g_coeff - df) / f_coeff
                g_dot = 1.0 + dg_dot
            elif abs(g_dot) > 1e-12:
                df = (f_dot * g_coeff - dg_dot) / g_dot
                f_coeff = 1.0 + df
        else:
            if abs(g_coeff) >= abs(f_dot) and abs(g_coeff) > 1e-12:
                f_dot = (df + dg_dot + df * dg_dot) / g_coeff
            elif abs(f_dot) > 1e-12:
                g_coeff = (df + dg_dot + df * dg_dot) / f_dot

    r_final = np.array([
        rx + (df * rx + g_coeff * vx),
        ry + (df * ry + g_coeff * vy),
        rz + (df * rz + g_coeff * vz),
    ], dtype=np.float64)

    mu_f_dot = mu * f_dot
    w_final = np.array([
        w0[0] + (mu_f_dot * rx + dg_dot * w0[0]),
        w0[1] + (mu_f_dot * ry + dg_dot * w0[1]),
        w0[2] + (mu_f_dot * rz + dg_dot * w0[2]),
    ], dtype=np.float64)
    return r_final, w_final


def compute_wronskian(
    r0: np.ndarray,
    w0: np.ndarray,
    mu: Union[float, mp.mpf],
    GM: Union[float, mp.mpf],
    dt_target: Union[float, mp.mpf],
) -> Union[float, mp.mpf]:
    """
    Evaluates the Wronskian determinant f * dot{g} - dot{f} * g of the Sundman-Stumpff flow.
    Must equal 1.0 to within machine precision (IEEE 754 < 1e-15 or float128 < 1e-32).
    """
    if (r0.dtype == object or w0.dtype == object or
        isinstance(mu, mp.mpf) or isinstance(GM, mp.mpf) or isinstance(dt_target, mp.mpf) or
        (len(r0) > 0 and isinstance(r0.flat[0], mp.mpf)) or
        (len(w0) > 0 and isinstance(w0.flat[0], mp.mpf))):
        r0_arr = np.array([mp.mpf(x) for x in r0], dtype=object)
        w0_arr = np.array([mp.mpf(x) for x in w0], dtype=object)
        mu_mp = mp.mpf(mu)
        GM_mp = mp.mpf(GM)
        dt_mp = mp.mpf(dt_target)
        if dt_mp == mp.mpf('0.0'):
            return mp.mpf('1.0')
        v0 = w0_arr / mu_mp
        rho0 = mp.sqrt(r0_arr[0]**2 + r0_arr[1]**2 + r0_arr[2]**2)
        sigma0 = r0_arr[0]*v0[0] + r0_arr[1]*v0[1] + r0_arr[2]*v0[2]
        v0_sq = v0[0]**2 + v0[1]**2 + v0[2]**2
        alpha = (mp.mpf('2.0') * GM_mp / rho0) - v0_sq
        tau = solve_sundman_time_equation_f128(rho0, sigma0, GM_mp, alpha, dt_mp)
        z = alpha * tau * tau
        c0, c1, c2, c3, _, _ = stumpff_c_f128(z)
        tau_sq = tau * tau
        rho_final = rho0 * c0 + sigma0 * tau * c1 + GM_mp * tau_sq * c2
        df = - (GM_mp / rho0) * tau_sq * c2
        g_coeff = rho0 * tau * c1 + sigma0 * tau_sq * c2
        f_dot = - (GM_mp / (rho0 * rho_final)) * tau * c1
        dg_dot = - (GM_mp / rho_final) * tau_sq * c2
        return mp.mpf('1.0') + (df + dg_dot + df * dg_dot - f_dot * g_coeff)

    if dt_target == 0.0:
        return 1.0

    v0 = w0 / mu
    rho0 = float(np.linalg.norm(r0))
    sigma0 = float(np.dot(r0, v0))
    v0_sq = float(np.dot(v0, v0))
    alpha = (2.0 * GM / rho0) - v0_sq

    tau = solve_sundman_time_equation(rho0, sigma0, GM, alpha, dt_target)
    z = alpha * tau * tau
    if z < Z_MIN_HYPERBOLIC:
        z = Z_MIN_HYPERBOLIC
    c0, c1, c2, c3, _, _ = stumpff_c(z)

    tau_sq = tau * tau
    rho_final = rho0 * c0 + sigma0 * tau * c1 + GM * tau_sq * c2

    df = - (GM / rho0) * tau_sq * c2
    g_coeff = rho0 * tau * c1 + sigma0 * tau_sq * c2
    f_dot = - (GM / (rho0 * rho_final)) * tau * c1
    dg_dot = - (GM / rho_final) * tau_sq * c2

    return float(1.0 + (df + dg_dot + df * dg_dot - f_dot * g_coeff))
