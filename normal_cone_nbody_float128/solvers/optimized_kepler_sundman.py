r"""
Optimized Universal Kepler-Sundman Diffeomorphic Leaf Solver & Fast Conormal Drift.

Mathematical & Theoretical Foundations:
- Monograph Reference: Chapter 3, Section 3.2 & Section 3.3 (Universal Kepler-Sundman Diffeomorphisms)
  and Chapter 5, Section 5.2.2 (Sundman-Stumpff Analytical Diffeomorphism).
- Strict Zero-Softening: \varepsilon = 0, d_0 = 0 preserved unconditionally.
- Lie-Bargmann & Noether Invariants: Total linear momentum, angular momentum, and
  unimodular transfer matrix determinant det(M) = f \dot{g} - \dot{f} g = 1 (Battin 1999)
  strictly preserved to machine precision float64 (< 10^-15).

Algorithmic Optimizations:
1. Micro-kernel Stumpff Evaluation (`stumpff_c0123_scalar`):
   - Horner order 6 polynomials for |z| < 0.1, closed-form trig/hyp for |z| >= 0.1.
   - Truncates computation to strictly needed functions c_0(z), c_1(z), c_2(z), c_3(z),
     completely bypassing unnecessary evaluation of c_4(z) and c_5(z).
2. Inlined Halley Iteration with Single-Pass Stumpff Evaluation:
   - Eliminates Python closure allocations (`_eval_sundman`) and stack frame overhead.
   - At convergence, directly reuses c_0, c_1, c_2, c_3 and \rho_{final} = f'(\tau),
     eliminating the redundant second call to Stumpff functions at the end of each step.
3. Fast Jacobi Coordinate Projection for N=3:
   - Eliminates recursive tree traversals and dictionary allocations during conormal drift.
   - Direct linear projection between Cartesian phase space and Jacobi tree chart.
"""

from __future__ import annotations
import math
from typing import Tuple, Optional, Dict, Set, List
import numpy as np

from ..core.constants import G_DEFAULT
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree, JacobiTopologyNode
from ..core.coordinates import JacobiCoordinates, cartesian_to_jacobi, jacobi_to_cartesian
from ..core.profiler import get_current_profiler
from .kepler_sundman import step_kepler_sundman, solve_sundman_time_equation, Z_MIN_HYPERBOLIC


def stumpff_c0123_scalar(z_val: float) -> Tuple[float, float, float, float]:
    """
    Evaluates only the first 4 universal Stumpff functions c_0(z) through c_3(z)
    for a scalar float z, avoiding unused computation of c_4 and c_5.

    Bit-exact identical to stumpff_c(z)[:4] across all z in R.
    """
    if abs(z_val) < 0.1:
        # Horner polynomial expansion of order 6 on c2, c3 + exact downward Stumpff recurrence for c0, c1
        c2 = 0.5 + z_val * (
            (-1.0 / 24.0) + z_val * (
                (1.0 / 720.0) + z_val * (
                    (-1.0 / 40320.0) + z_val * (
                        (1.0 / 3628800.0) + z_val * (
                            (-1.0 / 479001600.0) + z_val * (1.0 / 87178291200.0)
                        )
                    )
                )
            )
        )
        c3 = (1.0 / 6.0) + z_val * (
            (-1.0 / 120.0) + z_val * (
                (1.0 / 5040.0) + z_val * (
                    (-1.0 / 362880.0) + z_val * (
                        (1.0 / 39916800.0) + z_val * (
                            (-1.0 / 6227020800.0) + z_val * (1.0 / 1307674368000.0)
                        )
                    )
                )
            )
        )
        c0 = 1.0 - z_val * c2
        c1 = 1.0 - z_val * c3
        return c0, c1, c2, c3

    if z_val >= 0.1:
        # Elliptic regime (z > 0)
        sqrt_z = math.sqrt(z_val)
        sin_sz = math.sin(sqrt_z)
        cos_sz = math.cos(sqrt_z)
        c0 = cos_sz
        c1 = sin_sz / sqrt_z
        c2 = (1.0 - cos_sz) / z_val
        c3 = (sqrt_z - sin_sz) / (z_val * sqrt_z)
        return c0, c1, c2, c3

    # Hyperbolic regime (z < -0.1)
    sqrt_neg_z = math.sqrt(-z_val)
    sinh_sz = math.sinh(sqrt_neg_z)
    cosh_sz = math.cosh(sqrt_neg_z)
    c0 = cosh_sz
    c1 = sinh_sz / sqrt_neg_z
    c2 = (cosh_sz - 1.0) / (-z_val)
    c3 = (sinh_sz - sqrt_neg_z) / ((-z_val) * sqrt_neg_z)
    return c0, c1, c2, c3


def step_kepler_sundman_fast(
    r0: np.ndarray,
    w0: np.ndarray,
    mu: float,
    GM: float,
    dt_target: float,
    abs_tol: Optional[float] = None,
    max_iter: int = 25,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Optimized two-body Keplerian channel integration across physical step dt_target.
    Yields bit-exact numerical results matching step_kepler_sundman (< 1e-15 in float64,
    or 35-digit mpmath.mpf when inputs are object arrays).
    """
    if getattr(r0, "dtype", None) == object or getattr(w0, "dtype", None) == object:
        return step_kepler_sundman(r0, w0, mu, GM, dt_target)

    if dt_target == 0.0:
        return np.copy(r0), np.copy(w0)

    # Time-reversal symmetry handling for retrograde steps
    if dt_target < 0.0:
        r_rev, w_rev = step_kepler_sundman_fast(r0, -w0, mu, GM, -dt_target, abs_tol, max_iter)
        return r_rev, -w_rev

    inv_mu = 1.0 / mu
    vx = float(w0[0]) * inv_mu
    vy = float(w0[1]) * inv_mu
    vz = float(w0[2]) * inv_mu

    rx = float(r0[0])
    ry = float(r0[1])
    rz = float(r0[2])

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
        # Radial continuation fallback to baseline
        return step_kepler_sundman(r0, w0, mu, GM, dt_target)

    alpha = (2.0 * GM / rho0) - v0_sq
    abs_dt = dt_target
    eff_tol = (1e-15 * abs_dt) if abs_tol is None else (abs_tol * abs_dt if abs_dt < 1.0 else abs_tol)
    tau_floor = min(1e-15, 1e-15 * abs_dt / max(rho0, 1.0))

    # 1. Physics-based ballistic and collision initial guess tau_0
    r_collision_scale = (GM * (abs_dt**2)) ** (1.0 / 3.0)
    if rho0 < 0.1 * r_collision_scale and v0_sq * (abs_dt**2) < (r_collision_scale**2):
        tau = (6.0 * abs_dt / GM) ** (1.0 / 3.0)
    elif alpha < 0.0:
        sqrt_neg_alpha = math.sqrt(-alpha)
        arg = max((sqrt_neg_alpha * abs_dt) / max(rho0, 1e-30), 0.0)
        tau = math.asinh(arg) / sqrt_neg_alpha
    else:
        r_floor_sq = max(1e-24, 0.01 * (r_collision_scale**2))
        r_mid_sq = max(rho0_sq + sigma0 * abs_dt + 0.25 * v0_sq * (abs_dt**2), r_floor_sq)
        r_end_sq = max(rho0_sq + 2.0 * sigma0 * abs_dt + v0_sq * (abs_dt**2), r_floor_sq)
        tau = (abs_dt / 6.0) * (1.0 / max(rho0, math.sqrt(r_floor_sq)) + 4.0 / math.sqrt(r_mid_sq) + 1.0 / math.sqrt(r_end_sq))
        if alpha > 0.0:
            tau_ell_max = (alpha * abs_dt / max(GM, 1e-30) + math.pi) / math.sqrt(alpha)
            tau = min(tau, tau_ell_max)

    if alpha < 0.0:
        tau_max_hyp = math.sqrt(abs(Z_MIN_HYPERBOLIC) / abs(alpha))
        tau = min(tau, tau_max_hyp)

    tau = max(tau, tau_floor)

    # 2. Inlined cubic Halley iteration
    GM_minus_alpha_r0 = GM - alpha * rho0
    c0 = c1 = c2 = c3 = 0.0
    f1 = 0.0
    converged = False
    halley_iters = 0

    for _ in range(max_iter):
        halley_iters += 1
        zt = max(alpha * tau * tau, Z_MIN_HYPERBOLIC)
        c0, c1, c2, c3 = stumpff_c0123_scalar(zt)
        t_sq = tau * tau
        t_cb = t_sq * tau
        f_val = rho0 * tau * c1 + sigma0 * t_sq * c2 + GM * t_cb * c3 - abs_dt

        f_prime = rho0 * c0 + sigma0 * tau * c1 + GM * t_sq * c2

        if abs(f_val) <= eff_tol:
            f1 = f_prime
            converged = True
            break

        f_prime2 = sigma0 * c0 + GM_minus_alpha_r0 * tau * c1
        f_prime_safe = max(f_prime, 1e-30)
        denom = 2.0 * (f_prime_safe**2) - f_val * f_prime2

        if denom <= 0.1 * (f_prime_safe**2):
            step = f_val / f_prime_safe
        else:
            step = (2.0 * f_prime_safe * f_val) / denom

        if not np.isfinite(step) or abs(step) > 0.5 * tau:
            step = math.copysign(min(abs(step), 0.5 * tau), step)

        tau -= step
        if alpha < 0.0:
            tau = min(tau, tau_max_hyp)
        tau = max(tau, tau_floor)

    if not converged:
        # Fallback to robust bisection routine in baseline
        return step_kepler_sundman(r0, w0, mu, GM, dt_target)

    # 3. Final quantities: REUSE c0, c1, c2, c3 and f1 from the converged Halley step!
    tau_sq = tau * tau
    rho_raw = f1
    sigma_end = sigma0 * c0 + GM_minus_alpha_r0 * tau * c1

    # Record profiling metrics
    prof = get_current_profiler()
    if prof is not None:
        prof.record_kepler(
            halley_iters=halley_iters,
            stumpff_evals=halley_iters,
            is_bisection=False,
        )

    # 4. Head-on radial collision handling (L_sq ~ 0, e = 1)
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

    # 5. Regular 3D orbit (L > 0)
    denom_peri = 2.0 * GM - alpha * rho_raw
    if rho_raw < 0.25 * rho0 and denom_peri > 0.5 * GM and L_sq > 0.0:
        rho_final = (sigma_end * sigma_end + L_sq) / denom_peri
        rho_final = (sigma_end * sigma_end + L_sq) / (2.0 * GM - alpha * rho_final)
        f_prime = - (GM / rho0) * tau * c1
        g_prime = rho0 * c0 + sigma0 * tau * c1
        rpx = f_prime * rx + g_prime * vx
        rpy = f_prime * ry + g_prime * vy
        rpz = f_prime * rz + g_prime * vz
        lcr_x = ly * rpz - lz * rpy
        lcr_y = lz * rpx - lx * rpz
        lcr_z = lx * rpy - ly * rpx
        sx = sigma_end * rpx - lcr_x
        sy = sigma_end * rpy - lcr_y
        sz_v = sigma_end * rpz - lcr_z
        s_norm = math.sqrt(sx * sx + sy * sy + sz_v * sz_v)
        if s_norm > 0.0:
            inv_sn = 1.0 / s_norm
            shx, shy, shz = sx * inv_sn, sy * inv_sn, sz_v * inv_sn
            r_final = np.array([rho_final * shx, rho_final * shy, rho_final * shz], dtype=np.float64)
            r_actual = math.sqrt(r_final[0]**2 + r_final[1]**2 + r_final[2]**2)
            ux = ly * shz - lz * shy
            uy = lz * shx - lx * shz
            uz = lx * shy - ly * shx
            inv_r = mu / r_actual
            w_final = np.array([
                (sigma_end * shx + ux) * inv_r,
                (sigma_end * shy + uy) * inv_r,
                (sigma_end * shz + uz) * inv_r,
            ], dtype=np.float64)
            return r_final, w_final
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


class OptimizedKeplerSolver:
    """
    Decoupled Optimized Kepler-Sundman Solver & Fast Conormal Drift Engine.

    Features:
    - Inlined fast Kepler-Sundman solver (`step_kepler_sundman_fast`).
    - Caching of Jacobi linear projection structures per tree topology.
    - Specialized high-speed drift implementation for N=3 three-body problems.
    """

    def __init__(self, G: float = G_DEFAULT) -> None:
        self.G = G
        self._cached_tree_id: Optional[int] = None
        self._cached_tree_key: Optional[Tuple[int, Tuple[Tuple[int, int], ...]]] = None
        self._leaf_node_id: Optional[int] = None
        self._comp_node_id: Optional[int] = None
        self._i1: int = 0
        self._i2: int = 1
        self._i3: int = 2
        self._m1: float = 0.0
        self._m2: float = 0.0
        self._m3: float = 0.0
        self._m12: float = 0.0
        self._M_tot: float = 0.0
        self._mu1: float = 0.0
        self._mu2: float = 0.0
        self._c1: float = 0.0
        self._c2: float = 0.0

    def step_kepler(
        self,
        r0: np.ndarray,
        w0: np.ndarray,
        mu: float,
        GM: float,
        dt_target: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Solves a single Keplerian leaf channel."""
        return step_kepler_sundman_fast(r0, w0, mu, GM, dt_target)

    def _setup_n3_cache(self, tree: JacobiTree) -> None:
        """Pre-caches particle indices and mass fractions for N=3 tree topology."""
        if len(tree.leaf_pairs) != 1 or len(tree.composite_nodes) != 1:
            return

        leaf_node = tree.leaf_pairs[0]
        comp_node = tree.composite_nodes[0]

        assert leaf_node.left is not None and leaf_node.right is not None
        i1 = leaf_node.left.body_index
        i2 = leaf_node.right.body_index
        assert i1 is not None and i2 is not None

        # The third body is the one not in leaf_node
        all_bodies = {0, 1, 2}
        i3 = list(all_bodies - {i1, i2})[0]

        self._i1, self._i2, self._i3 = i1, i2, i3
        m = tree.masses
        self._m1, self._m2, self._m3 = float(m[i1]), float(m[i2]), float(m[i3])
        self._m12 = self._m1 + self._m2
        self._M_tot = self._m12 + self._m3
        self._mu1 = (self._m1 * self._m2) / self._m12
        self._mu2 = (self._m12 * self._m3) / self._M_tot
        self._c1 = self._m1 / self._m12
        self._c2 = self._m2 / self._m12
        self._leaf_node_id = leaf_node.id
        self._comp_node_id = comp_node.id
        self._cached_tree_id = id(tree)
        self._cached_tree_key = (id(tree), tuple(lp.body_indices for lp in tree.leaf_pairs))

    def step_drift_and_kepler(
        self,
        state: ExtendedPhaseState,
        tree: JacobiTree,
        dt: float,
        G: Optional[float] = None,
        supernode_bodies: Optional[Set[int]] = None,
        supernode_clusters: Optional[List[Tuple[int, ...]]] = None,
    ) -> ExtendedPhaseState:
        """
        Executes the composite drift and Kepler leaf step using the optimized solver.
        """
        if dt == 0.0:
            return state

        if getattr(state, "is_f128", False):
            from ..integrators.common import step_drift_and_kepler
            grav_const_f128 = self.G if G is None else G
            return step_drift_and_kepler(
                state,
                tree,
                dt,
                G=grav_const_f128,
                supernode_bodies=supernode_bodies,
                supernode_clusters=supernode_clusters,
            )

        grav_const = float(self.G) if G is None else float(G)
        from ..potential.screened_apex import _normalize_supernode_clusters
        clusters = _normalize_supernode_clusters(supernode_bodies, supernode_clusters)
        cluster_sets = [set(c) for c in clusters]
        sn_set = set().union(*cluster_sets) if cluster_sets else set()

        # Fast path for N=3 without supernodes
        if state.N == 3 and not sn_set and len(tree.leaf_pairs) == 1:
            tree_key = (id(tree), tuple(lp.body_indices for lp in tree.leaf_pairs))
            if self._cached_tree_key != tree_key:
                self._setup_n3_cache(tree)

            i1, i2, i3 = self._i1, self._i2, self._i3
            m1, m2, m3 = self._m1, self._m2, self._m3
            m12, M_tot = self._m12, self._M_tot
            mu1, mu2 = self._mu1, self._mu2
            c1, c2 = self._c1, self._c2

            q = state.q
            p = state.p

            # 1. Direct Cartesian to Jacobi for N=3
            r1 = q[i1] - q[i2]
            X12 = c1 * q[i1] + c2 * q[i2]
            r2 = X12 - q[i3]

            P12 = p[i1] + p[i2]
            w1 = mu1 * (p[i1] / m1 - p[i2] / m2)
            w2 = mu2 * (P12 / m12 - p[i3] / m3)

            # 2. Center of mass drift
            P_tot = state.total_momentum
            R_cm_next = state.center_of_mass + dt * (P_tot / M_tot)

            # 3. Composite node drift: r2 <- r2 + dt * (w2 / mu2)
            r2_next = r2 + dt * (w2 / mu2)
            w2_next = w2

            # 4. Leaf pair Kepler solve
            prof = get_current_profiler()
            mu_param = grav_const * m12
            if prof is not None:
                with prof.time_stage("kepler"):
                    r1_next, w1_next = step_kepler_sundman_fast(
                        r0=r1,
                        w0=w1,
                        mu=mu1,
                        GM=mu_param,
                        dt_target=dt,
                    )
            else:
                r1_next, w1_next = step_kepler_sundman_fast(
                    r0=r1,
                    w0=w1,
                    mu=mu1,
                    GM=mu_param,
                    dt_target=dt,
                )

            # 5. Direct Jacobi to Cartesian for N=3
            # X12_next = R_cm_next + (m3 / M_tot) * r2_next
            # q3_next  = R_cm_next - (m12 / M_tot) * r2_next
            c_comp_right = m3 / M_tot
            c_comp_left = m12 / M_tot
            X12_next = R_cm_next + c_comp_right * r2_next

            q_out = np.empty_like(q)
            q_out[i3] = R_cm_next - c_comp_left * r2_next
            q_out[i1] = X12_next + c2 * r1_next
            q_out[i2] = X12_next - c1 * r1_next

            P12_next = c_comp_left * P_tot + w2_next
            p_out = np.empty_like(p)
            p_out[i3] = c_comp_right * P_tot - w2_next
            p_out[i1] = c1 * P12_next + w1_next
            p_out[i2] = c2 * P12_next - w1_next

            return state.with_updates(
                q=q_out,
                p=p_out,
                t=state.t + dt,
            )

        # General path for arbitrary N or when supernodes are active
        coords = cartesian_to_jacobi(state, tree)
        r_dict = dict(coords.r)
        w_dict = dict(coords.w)

        R_cm_next = state.center_of_mass + dt * (state.total_momentum / tree.root.M)
        P_tot = state.total_momentum

        for node in tree.composite_nodes:
            if cluster_sets and any(set(node.body_indices).issubset(c_set) for c_set in cluster_sets):
                continue
            w_v = w_dict[node.id]
            r_dict[node.id] = r_dict[node.id] + dt * (w_v / node.mu)

        prof = get_current_profiler()
        for node in tree.leaf_pairs:
            if cluster_sets and any(set(node.body_indices).issubset(c_set) for c_set in cluster_sets):
                continue
            r0 = r_dict[node.id]
            w0 = w_dict[node.id]
            mu_param = grav_const * node.M

            if prof is not None:
                with prof.time_stage("kepler"):
                    r_next, w_next = step_kepler_sundman_fast(
                        r0=r0,
                        w0=w0,
                        mu=node.mu,
                        GM=mu_param,
                        dt_target=dt,
                    )
            else:
                r_next, w_next = step_kepler_sundman_fast(
                    r0=r0,
                    w0=w0,
                    mu=node.mu,
                    GM=mu_param,
                    dt_target=dt,
                )
            r_dict[node.id] = r_next
            w_dict[node.id] = w_next

        new_coords = JacobiCoordinates(r=r_dict, w=w_dict)
        updated_state = jacobi_to_cartesian(
            new_coords,
            tree,
            state.masses,
            t=state.t + dt,
            E=state.E,
            E0=state.E0,
            center_of_mass=R_cm_next,
            total_momentum=P_tot,
        )

        from ..integrators.common import step_supernode_flow
        updated_state = step_supernode_flow(
            updated_state,
            dt,
            G=grav_const,
            supernode_clusters=clusters,
            supernode_bodies=sn_set,
            drift_center_of_mass=False,
            kepler_solver=self,
        )
        return updated_state
