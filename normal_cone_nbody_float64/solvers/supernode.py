r"""
Multi-Body Synchronous Supernode Solver with Anti-Zeno Regularization.

Theoretical Foundations:
- Monograph Reference: Chapter 4, Section 4.3 (Synchronous Supernodes and Excluded Middle Dilemma),
  Chapter 4, Section 4.4 (Anti-Zeno Regularization via Sundman Time Unfolding), and Chapter 5, Section 5.1.
- Excluded Middle Dilemma:
    Asynchronous sub-stepping inside a compact multi-body cluster (k >= 3) induces a Poisson
    bracket divergence in O(rho_S^{-5/2}). Synchronous sub-stepping is geometrically mandated.
- Unfolded Sundman Time Clock:
    dt_S = chi_S^2 dtau_S = rho_S dtau_S, \\quad \\chi_S = \\sqrt{\\rho_S} > 0,
    where rho_S = (\\sum_{i \\in I_S} m_i ||q_i - q_S^{cdm}||^2)^{1/2}.
- Anti-Zeno Finite Traversal Property (Saari's Theorem):
    Almost everywhere on configuration space, rho_min > 0, bounding the required sub-steps:
    N_{sub} <= ceil(Delta t / (rho_min Delta tau_S)) < infty.
- Exact Conservation Guarantees:
    - Internal linear momentum P_S = sum_{i in I_S} p_i == cste (< 1e-15).
    - Internal angular momentum L_S = sum_{i in I_S} q_i x p_i == cste (< 1e-14).
    - Decoupled center of mass trajectory (zero drift from internal forces).
- Certified Formal Proofs:
    Lean 4 modules `PaperUnifiedFoundations.lean`, `01_UnifiedSuperNodeDef_Theorems.lean`,
    and `02_SuperNodeAntiZenoRegularization_Theorems.lean`.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, List, Optional, Any, Callable, Sequence
import math
import numpy as np

from .kepler_sundman import step_kepler_sundman
from ..core.constants import BM4_DRIFT_WEIGHTS, BM4_IMPULSE_WEIGHTS
from ..core.profiler import get_current_profiler


def _resolve_integrator_stages(
    drift_weights: Optional[Sequence[float]] = None,
    impulse_weights: Optional[Sequence[float]] = None,
) -> Tuple[Sequence[float], Sequence[float]]:
    """
    Resolves the active symplectic composition weights (drift and impulse stages).
    If not explicitly passed, queries the current execution context (set by the
    parent integrator, e.g. BM4, BM6, or Strang).
    If no active context exists, defaults to certified 7-stage 4th-order Blanes--Moan (BM4).
    """
    if drift_weights is not None and impulse_weights is not None:
        return drift_weights, impulse_weights
    try:
        from ..integrators.common import get_current_integrator_stages
        stages = get_current_integrator_stages()
        if stages is not None:
            return stages[0], stages[1]
    except ImportError:
        pass
    return BM4_DRIFT_WEIGHTS, BM4_IMPULSE_WEIGHTS


def _resolve_kepler_step(kepler_solver: Optional[Any] = None) -> Callable[..., Tuple[np.ndarray, np.ndarray]]:
    """
    Dynamically resolves the active Kepler-Sundman step function.
    If an explicit solver is passed or set in the ContextVar, delegates
    to its step method/callable or `step_kepler_sundman_fast`.
    Otherwise, defaults unconditionally to canonical baseline `step_kepler_sundman`.
    """
    if kepler_solver is None:
        try:
            from ..integrators.common import get_current_kepler_solver
            kepler_solver = get_current_kepler_solver()
        except ImportError:
            kepler_solver = None

    if kepler_solver is not None:
        from .optimized_kepler_sundman import step_kepler_sundman_fast
        cls_name = type(kepler_solver).__name__
        if cls_name in ("OptimizedKeplerSolver", "OptimizedKeplerSundmanSolver"):
            if getattr(kepler_solver, "tol", None) is None and getattr(kepler_solver, "max_iter", 50) == 50:
                return step_kepler_sundman_fast
        if hasattr(kepler_solver, "step_kepler"):
            return kepler_solver.step_kepler
        if callable(kepler_solver):
            return kepler_solver
        return step_kepler_sundman_fast
    return step_kepler_sundman


@dataclass(frozen=True)
class SupernodeCluster:
    """
    Representation of a detected synchronous supernode cluster of k >= 3 bodies.

    Attributes:
        body_indices: Tuple of integer indices belonging to the supernode.
        mass: Total mass M_S = sum_{i in I_S} m_i.
        hyper_radius: Instantaneous hyper-radius rho_S.
        center_of_mass_q: Center of mass position in R^3.
        center_of_mass_p: Center of mass total momentum in R^3.
    """
    body_indices: Tuple[int, ...]
    mass: float
    hyper_radius: float
    center_of_mass_q: np.ndarray
    center_of_mass_p: np.ndarray

    def __post_init__(self) -> None:
        if len(self.body_indices) < 2:
            raise ValueError(f"Supernode must contain at least 2 bodies, got {len(self.body_indices)}.")
        if self.mass <= 0.0:
            raise ValueError(f"Total cluster mass must be strictly positive, got {self.mass}.")
        if self.hyper_radius < 0.0:
            raise ValueError(f"Hyper-radius must be non-negative, got {self.hyper_radius}.")


def compute_cluster_hyper_radius(
    q: np.ndarray,
    masses: np.ndarray,
    body_indices: Tuple[int, ...],
) -> Tuple[float, np.ndarray]:
    """
    Computes the cluster hyper-radius rho_S and center of mass:
        rho_S = (sum_{i in I_S} m_i ||q_i - q_S^{cdm}||^2)^{1/2}.
    """
    idx = list(body_indices)
    m_sub = masses[idx]
    q_sub = q[idx]
    m_tot = float(np.sum(m_sub))
    if m_tot <= 0.0:
        raise ValueError(f"Cluster mass must be positive, got {m_tot}.")

    q_cdm = np.sum(m_sub[:, None] * q_sub, axis=0) / m_tot
    diff = q_sub - q_cdm
    sq_dists = np.sum(diff ** 2, axis=1)
    rho_sq = float(np.sum(m_sub * sq_dists))
    rho = float(np.sqrt(max(0.0, rho_sq)))
    return rho, q_cdm


def _compute_internal_cluster_acc(
    x_arr: np.ndarray,
    m_sub: np.ndarray,
    G: float,
) -> Tuple[np.ndarray, float]:
    """
    Computes exact unsoftened internal pairwise Newtonian accelerations and minimal separation.
    Zero softening (d_safe = 0.0), self-interactions strictly masked to zero.
    Strictly protected against division by zero and NaN even during exact collision (d_ij = 0.0).
    """
    dr = x_arr[None, :, :] - x_arr[:, None, :]  # shape (k, k, 3), dr[i, j] = x[j] - x[i]
    d2 = np.sum(dr ** 2, axis=-1)
    np.fill_diagonal(d2, 1.0)
    d = np.sqrt(d2)

    nonzero_mask = d2 > 0.0
    inv_d3 = np.zeros_like(d2)
    inv_d3[nonzero_mask] = G / (d2[nonzero_mask] * d[nonzero_mask])
    np.fill_diagonal(inv_d3, 0.0)

    acc = np.sum(m_sub[None, :, None] * inv_d3[:, :, None] * dr, axis=1)
    np.fill_diagonal(d, np.inf)
    r_min = float(np.min(d)) if len(m_sub) >= 2 else 1.0
    return acc, r_min


def _compute_internal_apex_forces(
    x_arr: np.ndarray,
    m_sub: np.ndarray,
    G: float,
    pair_idx: Optional[Tuple[int, int]] = None,
) -> np.ndarray:
    """
    Computes exact unsoftened internal gravitational forces on all cluster bodies.
    Under the recursive normal cone splitting (Equation 6.72):
        H_S^{int} = H_{leaf}^{int} + H_{base}^{int} + V_{apex}^{int}
    If pair_idx = (i_min, j_min) is provided, the singular binary interaction between
    i_min and j_min is omitted from V_{apex}^{int} because it is integrated analytically
    by the universal Kepler-Sundman solver.
    The resulting cross-forces strictly satisfy exact linear momentum conservation:
        sum_i F_i = 0,
    and angular momentum conservation:
        sum_i x_i x F_i = 0.
    """
    k = len(m_sub)
    F = np.zeros_like(x_arr)
    for i in range(k):
        for j in range(i + 1, k):
            if pair_idx is not None and (i, j) == pair_idx:
                continue
            disp = x_arr[j] - x_arr[i]
            d2_ij = float(np.sum(disp ** 2))
            if d2_ij > 0.0:
                dist = math.sqrt(d2_ij)
                f = (G * m_sub[i] * m_sub[j] / (d2_ij * dist)) * disp
                F[i] += f
                F[j] -= f
    return F


def _step_cluster_internal_multistage(
    x: np.ndarray,
    v: np.ndarray,
    m_sub: np.ndarray,
    h_step: float,
    G: float,
    pair_idx: Optional[Tuple[int, int]],
    kepler_solver: Optional[Any],
    drift_weights: Sequence[float],
    impulse_weights: Sequence[float],
) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Executes one elementary substep h_step using a palindromic partitioned Runge-Kutta-Nyström
    splitting (e.g. BM4, BM6, or Strang) with recursive normal cone decomposition:
    - Kick phase b_s * h_step: advance velocities under internal apex force F^{apex} = - grad V_{apex}^{int}
    - Drift / Kepler phase a_s * h_step:
        * Analytical Kepler-Sundman flow for pair_idx = (i*, j*)
        * Free linear kinetic drift for spectator bodies in S \ {i*, j*}
    """
    k = len(m_sub)
    x_out = np.copy(x)
    v_out = np.copy(v)
    kepler_step_fn = _resolve_kepler_step(kepler_solver)

    n_drifts = len(drift_weights)
    n_kicks = len(impulse_weights)

    for s in range(n_kicks):
        # 1. Apex kick stage
        b_s = impulse_weights[s]
        if b_s != 0.0:
            dt_kick = b_s * h_step
            F_apex = _compute_internal_apex_forces(x_out, m_sub, G, pair_idx=pair_idx)
            v_out += dt_kick * (F_apex / m_sub[:, None])

        # 2. Drift / Kepler stage
        if s < n_drifts:
            a_s = drift_weights[s]
            if a_s != 0.0:
                dt_drift = a_s * h_step
                if pair_idx is not None:
                    i_min, j_min = pair_idx
                    m1, m2 = m_sub[i_min], m_sub[j_min]
                    M_p = m1 + m2
                    mu_p = (m1 * m2) / M_p

                    r_p = x_out[j_min] - x_out[i_min]
                    w_p = mu_p * (v_out[j_min] - v_out[i_min])
                    X_p = (m1 * x_out[i_min] + m2 * x_out[j_min]) / M_p
                    V_p = (m1 * v_out[i_min] + m2 * v_out[j_min]) / M_p

                    prof = get_current_profiler()
                    if prof is not None:
                        with prof.time_stage("kepler"):
                            r_p_new, w_p_new = kepler_step_fn(
                                r_p, w_p, mu=mu_p, GM=G * M_p, dt_target=dt_drift
                            )
                    else:
                        r_p_new, w_p_new = kepler_step_fn(
                            r_p, w_p, mu=mu_p, GM=G * M_p, dt_target=dt_drift
                        )
                    X_p_new = X_p + dt_drift * V_p

                    x_out[i_min] = X_p_new - (m2 / M_p) * r_p_new
                    x_out[j_min] = X_p_new + (m1 / M_p) * r_p_new
                    v_out[i_min] = V_p - (w_p_new / m1)
                    v_out[j_min] = V_p + (w_p_new / m2)

                    # Spectators linear kinetic drift
                    for l in range(k):
                        if l != i_min and l != j_min:
                            x_out[l] += dt_drift * v_out[l]
                else:
                    # Full linear kinetic drift of all cluster bodies
                    x_out += dt_drift * v_out

    return x_out, v_out


def _identify_closest_pair_and_distances(
    x: np.ndarray,
    k: int,
) -> Tuple[Tuple[int, int], float, float]:
    """
    Identifies the closest pair (i_min, j_min) with i_min < j_min, minimal distance r_min,
    and second minimal distance r_second.
    """
    dr = x[None, :, :] - x[:, None, :]  # shape (k, k, 3)
    d2 = np.sum(dr ** 2, axis=-1)
    np.fill_diagonal(d2, np.inf)
    d = np.sqrt(d2)

    i_min, j_min = np.unravel_index(np.argmin(d), d.shape)
    i_min, j_min = int(i_min), int(j_min)
    if i_min > j_min:
        i_min, j_min = j_min, i_min
    r_min = float(d[i_min, j_min])

    d_copy = np.copy(d)
    d_copy[i_min, j_min] = np.inf
    d_copy[j_min, i_min] = np.inf
    r_second = float(np.min(d_copy)) if k > 2 else float(r_min)
    return (i_min, j_min), r_min, r_second


def _compute_local_substep_bound(
    x: np.ndarray,
    v: np.ndarray,
    m_sub: np.ndarray,
    m_tot: float,
    G: float,
) -> float:
    """
    Computes the local adaptive time-step bound h_bound(x, v), invariant under time reversal
    (x, v) -> (x, -v) and quantized to float32 to eliminate trajectory roundoff sensitivity
    at dyadic bisection boundaries.
    """
    k = len(m_sub)
    _, _, r_second = _identify_closest_pair_and_distances(x, k)
    rho_S = float(np.sqrt(max(0.0, np.sum(x ** 2))))
    v_max = float(np.max(np.linalg.norm(v, axis=1)))

    tau_smooth = 2.0 * math.pi * math.sqrt((max(r_second, 1e-12) ** 3) / max(G * m_tot, 1e-30))
    tau_kin = max(r_second, 1e-12) / (v_max + 1e-12)
    tau_sundman = max(rho_S, 1e-12) / (v_max + 1e-12)
    h_raw = min(0.05 * tau_smooth, 0.05 * tau_kin, 0.05 * tau_sundman)
    return float(np.float32(max(h_raw, 1e-16)))


def _step_cluster_internal_substep(
    x: np.ndarray,
    v: np.ndarray,
    m_sub: np.ndarray,
    h_step: float,
    G: float,
    kepler_solver: Optional[Any] = None,
    drift_weights: Optional[Sequence[float]] = None,
    impulse_weights: Optional[Sequence[float]] = None,
) -> Tuple[np.ndarray, np.ndarray, float, float, bool]:
    """
    Performs one multi-stage symplectic sub-step of internal cluster motion under the
    recursive normal cone splitting (Equation 6.72) with a Symmetric Half-Step Pair Bridge
    if the closest pair switches during the substep, guaranteeing exact time-reversal
    symmetry R o Phi_{-h} o R = Phi_{+h}^{-1}.

    Returns:
        (x_next, v_next, r_min, r_second, has_close_pair)
    """
    k = len(m_sub)
    pair_start, r_min, r_second = _identify_closest_pair_and_distances(x, k)
    has_close_pair = True
    dw, iw = _resolve_integrator_stages(drift_weights, impulse_weights)

    x_out, v_out = _step_cluster_internal_multistage(
        x, v, m_sub, h_step, G,
        pair_idx=pair_start,
        kepler_solver=kepler_solver,
        drift_weights=dw,
        impulse_weights=iw,
    )

    if k > 2:
        pair_end, _, _ = _identify_closest_pair_and_distances(x_out, k)
        if pair_end != pair_start:
            # Symmetric Half-Step Pair Bridge: Phi_{h/2}^{(pair_end)} o Phi_{h/2}^{(pair_start)}
            x_mid, v_mid = _step_cluster_internal_multistage(
                x, v, m_sub, 0.5 * h_step, G,
                pair_idx=pair_start,
                kepler_solver=kepler_solver,
                drift_weights=dw,
                impulse_weights=iw,
            )
            x_out, v_out = _step_cluster_internal_multistage(
                x_mid, v_mid, m_sub, 0.5 * h_step, G,
                pair_idx=pair_end,
                kepler_solver=kepler_solver,
                drift_weights=dw,
                impulse_weights=iw,
            )

    return x_out, v_out, r_min, r_second, has_close_pair


def step_supernode_synchronous(
    q: np.ndarray,
    p: np.ndarray,
    masses: np.ndarray,
    body_indices: Tuple[int, ...],
    dt: float,
    G: float = 1.0,
    n_substeps: Optional[int] = None,
    drift_center_of_mass: bool = True,
    max_steps: int = 100_000,
    kepler_solver: Optional[Any] = None,
    drift_weights: Optional[Sequence[float]] = None,
    impulse_weights: Optional[Sequence[float]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrates the internal degrees of freedom of a synchronous supernode cluster
    under an unfolded Sundman clock with finite dyadic sub-stepping (Anti-Zeno regularized).
    Close binary encounters inside the cluster (k >= 3) are regularized via analytical
    Stumpff-Sundman isolation, preventing b^{-2} energy blowup and temporal stagnation.

    Parameters:
        q: Cartesian positions of all N bodies (shape (N, 3)).
        p: Cartesian momenta of all N bodies (shape (N, 3)).
        masses: Masses of all N bodies (shape (N,)).
        body_indices: Indices of the bodies belonging to the supernode.
        dt: Macro physical time step.
        G: Gravitational constant.
        n_substeps: Optional explicit sub-step count. If None, dynamically determined from
                    the internal dynamical timescale tau_dyn ~ sqrt(r_min^3 / (G M_S)).
        drift_center_of_mass: Whether the supernode center of mass drifts under Galilean inertia.
        max_steps: Maximum allowable sub-steps before raising RuntimeError for incomplete integration.
        kepler_solver: Optional Kepler solver instance or token; if None, resolves via ContextVar or baseline.

    Returns:
        (q_new, p_new): Updated positions and momenta with internal evolution.
    """
    if dt == 0.0:
        return np.copy(q), np.copy(p)

    prof = get_current_profiler()
    if prof is not None:
        prof.record_supernode(activations=1)

    kepler_step_fn = _resolve_kepler_step(kepler_solver)

    idx = list(body_indices)
    k = len(idx)
    m_sub = masses[idx]
    q_sub = q[idx]
    p_sub = p[idx]

    m_tot = float(np.sum(m_sub))
    q_cdm = np.sum(m_sub[:, None] * q_sub, axis=0) / m_tot
    p_cdm = np.sum(p_sub, axis=0)

    # Relative internal coordinates and velocities
    x = q_sub - q_cdm
    v = p_sub / m_sub[:, None] - (p_cdm / m_tot)
    q_cdm_target = q_cdm + dt * (p_cdm / m_tot) if drift_center_of_mass else q_cdm

    # Pure binary supernode (k = 2)
    if k == 2:
        m1, m2 = m_sub[0], m_sub[1]
        M_p = m1 + m2
        mu_p = (m1 * m2) / M_p
        r_p = x[1] - x[0]
        w_p = mu_p * (v[1] - v[0])
        X_p = (m1 * x[0] + m2 * x[1]) / M_p
        V_p = (m1 * v[0] + m2 * v[1]) / M_p

        if prof is not None:
            with prof.time_stage("kepler"):
                r_p_new, w_p_new = kepler_step_fn(r_p, w_p, mu=mu_p, GM=G * M_p, dt_target=dt)
        else:
            r_p_new, w_p_new = kepler_step_fn(r_p, w_p, mu=mu_p, GM=G * M_p, dt_target=dt)
        X_p_new = X_p + dt * V_p

        x_cur = np.zeros_like(x)
        v_cur = np.zeros_like(v)
        x_cur[0] = X_p_new - (m2 / M_p) * r_p_new
        x_cur[1] = X_p_new + (m1 / M_p) * r_p_new
        v_cur[0] = V_p - (w_p_new / m1)
        v_cur[1] = V_p + (w_p_new / m2)

        q_new = np.copy(q)
        p_new = np.copy(p)
        q_new[idx] = x_cur + q_cdm_target
        p_new[idx] = m_sub[:, None] * (v_cur + p_cdm / m_tot)
        return q_new, p_new

    r_norms = np.linalg.norm(x, axis=1)
    r_mean = float(np.mean(r_norms))
    bar_rho_sq = float(np.sum(m_sub * (r_norms ** 2))) / m_tot

    # Internal angular momentum pre-filter (Eq. 6.94 / 6.107: L^2 < 1e-20 * (bar_rho_S^2 + 1e-24))
    L_int = np.sum([np.cross(x[i], m_sub[i] * v[i]) for i in range(k)], axis=0)
    L_sq = float(np.dot(L_int, L_int))

    is_radial_symmetric = False
    lam_mean = 0.0
    if L_sq < 1e-20 * (bar_rho_sq + 1e-24):
        # Check for homothety central configuration collapse (exact radial collapse)
        acc_int, _ = _compute_internal_cluster_acc(x, m_sub, G)

        lambdas = []
        for i in range(k):
            if r_norms[i] > 1e-12:
                lam = - float(np.dot(acc_int[i], x[i])) / (r_norms[i] ** 2)
                lambdas.append(lam)

        if len(lambdas) == k and all(l > 0.0 for l in lambdas):
            lam_mean = float(np.mean(lambdas))
            lam_std = float(np.std(lambdas))
            if lam_std / lam_mean < 0.02:
                is_radial_symmetric = True

    if is_radial_symmetric:
        r0 = max(r_mean, 1e-15)
        u_dirs = x / r_norms[:, None]
        v_rad = float(np.mean([np.dot(u_dirs[i], v[i]) for i in range(k)]))
        GM_eff = float(lam_mean * (r0 ** 3))

        r_1d = np.array([r0, 0.0, 0.0], dtype=np.float64)
        w_1d = np.array([v_rad, 0.0, 0.0], dtype=np.float64)  # mu_eff = 1.0

        if prof is not None:
            with prof.time_stage("kepler"):
                r_next_1d, w_next_1d = kepler_step_fn(
                    r_1d, w_1d, mu=1.0, GM=GM_eff, dt_target=dt
                )
        else:
            r_next_1d, w_next_1d = kepler_step_fn(
                r_1d, w_1d, mu=1.0, GM=GM_eff, dt_target=dt
            )
        r_next = float(r_next_1d[0])
        v_next = float(w_next_1d[0])

        x_new = r_next * u_dirs
        v_new = v_next * u_dirs

        q_new = np.copy(q)
        p_new = np.copy(p)
        q_new[idx] = x_new + q_cdm_target
        p_new[idx] = (v_new + p_cdm / m_tot) * m_sub[:, None]
        return q_new, p_new

    # General multi-body cluster: conormal synchronous sub-stepping with close-pair regularizer
    x_cur = np.copy(x)
    v_cur = np.copy(v)

    dw, iw = _resolve_integrator_stages(drift_weights, impulse_weights)

    if n_substeps is not None:
        # Explicit uniform sub-step count (e.g. from unit tests)
        n_steps = max(1, int(n_substeps))
        h = dt / n_steps
        for _ in range(n_steps):
            x_cur, v_cur, _, _, _ = _step_cluster_internal_substep(
                x_cur, v_cur, m_sub, h, G,
                kepler_solver=kepler_solver,
                drift_weights=dw,
                impulse_weights=iw,
            )
    else:
        # Symmetric dyadic bisection stack under unfolded Sundman clock with close-pair desingularization.
        # Because the acceptance criterion |h_try| <= min(h_start, h_end) is invariant under time reversal
        # (x_cur, v_cur) <-> (x_next, -v_next) and 0.5 * h_try is bit-exact in binary floating-point,
        # forward and backward passes traverse the exact same dyadic tree of sub-steps in reverse order.
        abs_dt = abs(dt)
        t_accum = 0.0
        steps = 0
        stack: List[float] = [float(dt)]
        max_stack_depth = 64

        while stack and steps < max_steps:
            h_try = stack.pop()
            abs_h_try = abs(h_try)

            h_start = _compute_local_substep_bound(x_cur, v_cur, m_sub, m_tot, G)
            if abs_h_try > h_start and (len(stack) + 2 <= max_stack_depth) and abs_h_try > 1e-15:
                half_h = 0.5 * h_try
                stack.append(half_h)
                stack.append(half_h)
                continue

            x_next, v_next, _, _, has_close_pair = _step_cluster_internal_substep(
                x_cur, v_cur, m_sub, h_try, G,
                kepler_solver=kepler_solver,
                drift_weights=dw,
                impulse_weights=iw,
            )

            h_end = _compute_local_substep_bound(x_next, v_next, m_sub, m_tot, G)
            if abs_h_try > h_end and (len(stack) + 2 <= max_stack_depth) and abs_h_try > 1e-15:
                half_h = 0.5 * h_try
                stack.append(half_h)
                stack.append(half_h)
                continue

            x_cur = x_next
            v_cur = v_next

            if prof is not None:
                prof.record_supernode(substeps=1, acc_evals=len(iw), is_close_pair=has_close_pair)

            t_accum += abs_h_try
            steps += 1

        if stack or (abs_dt - t_accum > 1e-12 * abs_dt):
            raise RuntimeError(
                f"Supernode internal integration reached max_steps ({max_steps}) with incomplete physical time: "
                f"t_accum={t_accum:.6e} < abs_dt={abs_dt:.6e}."
            )

    # Canonical Dirac projection of machine-precision roundoff residual onto constraint manifold sum m_i v_i = 0
    P_rel_cdm = np.sum(m_sub[:, None] * v_cur, axis=0)
    if np.linalg.norm(P_rel_cdm) > 0.0:
        v_cur -= (P_rel_cdm / m_tot)

    q_new = np.copy(q)
    p_new = np.copy(p)
    q_new[idx] = x_cur + q_cdm_target
    p_new[idx] = m_sub[:, None] * (v_cur + p_cdm / m_tot)
    return q_new, p_new
