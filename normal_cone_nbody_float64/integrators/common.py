"""
Common Symplectic Sub-Step Operators: Kinetic Drift, Stumpff-Sundman Kepler, and Apex Kick.

Mathematical Monograph Alignment:
- Chapter 5: Section 5.1, Equations (5.10)-(5.24) — Elementary Lie group flows
  (center-of-mass drift, composite-node kinematics, analytic Kepler-Sundman leaf pairs,
   synchronous supernodes, and screened apex cotangent kicks).
- Formal Lean 4 Alignment: PaperUnifiedFoundations.lean
  (sundman_stumpff_diffeomorphism, lie_bargmann_closure_dim10).
"""

from __future__ import annotations
from typing import Optional, Set, Dict, List, Tuple, Any, Sequence
import numpy as np

from contextvars import ContextVar, Token
from ..core.constants import G_DEFAULT
from ..core.state import ExtendedPhaseState
from ..core.jacobi import JacobiTree
from ..core.coordinates import cartesian_to_jacobi, jacobi_to_cartesian, JacobiCoordinates
from ..solvers.kepler_sundman import step_kepler_sundman
from ..solvers.supernode import step_supernode_synchronous
from ..potential.screened_apex import ScreenedApexPotential, _normalize_supernode_clusters
from ..core.profiler import get_current_profiler


_CURRENT_KEPLER_SOLVER: ContextVar[Optional[Any]] = ContextVar("current_kepler_solver", default=None)
_CURRENT_SUPERNODE_SOLVER: ContextVar[Optional[Any]] = ContextVar("current_supernode_solver", default=None)
_CURRENT_APEX_EVALUATOR: ContextVar[Optional[Any]] = ContextVar("current_apex_evaluator", default=None)
_CURRENT_BYPASS_KEPLER: ContextVar[bool] = ContextVar("current_bypass_kepler", default=False)
_CURRENT_INTEGRATOR_STAGES: ContextVar[Optional[Tuple[Sequence[float], Sequence[float]]]] = ContextVar(
    "current_integrator_stages", default=None
)
_CURRENT_SUPERNODE_N_SUBSTEPS: ContextVar[Optional[int]] = ContextVar(
    "current_supernode_n_substeps", default=None
)


def set_current_supernode_n_substeps(n_substeps: Optional[int]) -> Token:
    """Sets the forced n_substeps override for synchronous supernodes (Canonical Rule 3)."""
    return _CURRENT_SUPERNODE_N_SUBSTEPS.set(n_substeps)


def get_current_supernode_n_substeps() -> Optional[int]:
    """Retrieves the forced n_substeps override for synchronous supernodes, or None if adaptive."""
    return _CURRENT_SUPERNODE_N_SUBSTEPS.get()


def reset_current_supernode_n_substeps(token: Token) -> None:
    """Resets the forced n_substeps override ContextVar using the provided token."""
    _CURRENT_SUPERNODE_N_SUBSTEPS.reset(token)


def set_current_kepler_solver(solver: Optional[Any]) -> Token:
    """Sets the active Kepler solver for the current execution context."""
    return _CURRENT_KEPLER_SOLVER.set(solver)


def get_current_kepler_solver() -> Optional[Any]:
    """Retrieves the active Kepler solver, or None if using canonical baseline."""
    return _CURRENT_KEPLER_SOLVER.get()


def reset_current_kepler_solver(token: Token) -> None:
    """Resets the active Kepler solver ContextVar using the provided token."""
    _CURRENT_KEPLER_SOLVER.reset(token)


def set_current_supernode_solver(solver: Optional[Any]) -> Token:
    """Sets the active Supernode solver for the current execution context."""
    return _CURRENT_SUPERNODE_SOLVER.set(solver)


def get_current_supernode_solver() -> Optional[Any]:
    """Retrieves the active Supernode solver, or None if using canonical baseline."""
    return _CURRENT_SUPERNODE_SOLVER.get()


def reset_current_supernode_solver(token: Token) -> None:
    """Resets the active Supernode solver ContextVar using the provided token."""
    _CURRENT_SUPERNODE_SOLVER.reset(token)


def set_current_apex_evaluator(evaluator: Optional[Any]) -> Token:
    """Sets the active Screened Apex potential evaluator for the current execution context (§6.6.4 Eq. 6.132)."""
    return _CURRENT_APEX_EVALUATOR.set(evaluator)


def get_current_apex_evaluator() -> Optional[Any]:
    """Retrieves the active Screened Apex potential evaluator, or None if using canonical baseline."""
    return _CURRENT_APEX_EVALUATOR.get()


def reset_current_apex_evaluator(token: Token) -> None:
    """Resets the active Screened Apex potential evaluator ContextVar using the provided token."""
    _CURRENT_APEX_EVALUATOR.reset(token)


def set_current_bypass_kepler(bypass: bool) -> Token:
    """Sets whether to bypass the Kepler-Sundman solver and execute pure linear drift."""
    return _CURRENT_BYPASS_KEPLER.set(bypass)


def get_current_bypass_kepler() -> bool:
    """Retrieves whether Kepler-Sundman solver bypass is active."""
    return _CURRENT_BYPASS_KEPLER.get()


def reset_current_bypass_kepler(token: Token) -> None:
    """Resets the Kepler-Sundman solver bypass ContextVar using the provided token."""
    _CURRENT_BYPASS_KEPLER.reset(token)


def set_current_integrator_stages(stages: Optional[Tuple[Sequence[float], Sequence[float]]]) -> Token:
    """Sets the active parent integrator stages (drift_weights, impulse_weights) in execution context."""
    return _CURRENT_INTEGRATOR_STAGES.set(stages)


def get_current_integrator_stages() -> Optional[Tuple[Sequence[float], Sequence[float]]]:
    """Retrieves the active parent integrator stages, or None if none set."""
    return _CURRENT_INTEGRATOR_STAGES.get()


def reset_current_integrator_stages(token: Token) -> None:
    """Resets the active parent integrator stages ContextVar using the provided token."""
    _CURRENT_INTEGRATOR_STAGES.reset(token)



def step_drift_and_kepler(
    state: ExtendedPhaseState,
    tree: JacobiTree,
    dt: float,
    G: float = G_DEFAULT,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[List[Tuple[int, ...]]] = None,
    kepler_solver: Optional[Any] = None,
    drift_weights: Optional[Sequence[float]] = None,
    impulse_weights: Optional[Sequence[float]] = None,
) -> ExtendedPhaseState:
    """
    Executes the exact analytic decoupled Hamiltonian flow for duration dt:
        Phi_{dt} = Phi_{dt}^{T_0} o Phi_{dt}^{T_{comp}} o Phi_{dt}^{H_{leaf}} o Phi_{dt}^{H_{super}}.

    Properties:
    - Center-of-mass momentum P_tot is conserved identically.
    - Wronskian unimodularity |f g_dot - f_dot g - 1| < 1e-15 on all leaf pairs.
    - Exact time reversibility under retrograde sub-steps (dt < 0).
    - Preserves Darboux pairing and Lie-Bargmann linear and angular momentum.
    """
    if dt == 0.0 or state.N <= 1:
        return state.with_updates(t=state.t + dt)

    active_kepler = kepler_solver or get_current_kepler_solver()
    if active_kepler is not None and hasattr(active_kepler, "step_drift_and_kepler"):
        return active_kepler.step_drift_and_kepler(
            state,
            tree,
            dt,
            G=G,
            supernode_bodies=supernode_bodies,
            supernode_clusters=supernode_clusters,
            drift_weights=drift_weights,
            impulse_weights=impulse_weights,
        )
    step_kepler_fn = active_kepler if callable(active_kepler) else step_kepler_sundman

    clusters = _normalize_supernode_clusters(supernode_bodies, supernode_clusters)
    cluster_sets = [set(c) for c in clusters]
    sn_set: Set[int] = set().union(*cluster_sets) if cluster_sets else set()

    # Convert Cartesian state to canonical Jacobi tree coordinates
    coords = cartesian_to_jacobi(state, tree)
    r_dict = dict(coords.r)
    w_dict = dict(coords.w)
    dr_dict = {nid: np.zeros(3, dtype=np.float64) for nid in r_dict}
    dw_dict = {nid: np.zeros(3, dtype=np.float64) for nid in w_dict}

    # 1. Center-of-mass drift: d_cm = dt * (P_tot / M_tot)
    d_cm = dt * (state.total_momentum / tree.root.M)

    # 2. Composite nodes drift: dr_v <- dt * (w_v / mu_v) for v in V_comp^free(T)
    for node in tree.composite_nodes:
        # Exclude composite node iff all its bodies belong to a single active supernode (§4.2 Def. 4.3)
        if cluster_sets and any(set(node.body_indices).issubset(c_set) for c_set in cluster_sets):
            continue
        w_v = w_dict[node.id]
        dr_dict[node.id] = dt * (w_v / node.mu)

    # 3. Leaf pairs Keplerian step in Stumpff-Sundman universal series for v in V_leaf^free(T)
    bypass_kepler = get_current_bypass_kepler()
    prof = get_current_profiler()
    for node in tree.leaf_pairs:
        # Exclude leaf pair iff both bodies belong to a single active supernode
        if cluster_sets and any(set(node.body_indices).issubset(c_set) for c_set in cluster_sets):
            continue
        r0 = r_dict[node.id]
        w0 = w_dict[node.id]
        if bypass_kepler:
            dr_dict[node.id] = dt * (w0 / node.mu)
            continue
        mu_param = G * node.M

        if prof is not None:
            with prof.time_stage("kepler"):
                r_next, w_next = step_kepler_fn(
                    r0=r0,
                    w0=w0,
                    mu=node.mu,
                    GM=mu_param,
                    dt_target=dt,
                )
        else:
            r_next, w_next = step_kepler_fn(
                r0=r0,
                w0=w0,
                mu=node.mu,
                GM=mu_param,
                dt_target=dt,
            )
        dr_dict[node.id] = r_next - r0
        dw_dict[node.id] = w_next - w0

    # 4. Invert Jacobi increments back to Cartesian configuration space (loss-free incremental update)
    delta_state = jacobi_to_cartesian(
        JacobiCoordinates(r=dr_dict, w=dw_dict),
        tree,
        state.masses,
        t=0.0,
        E=0.0,
        E0=0.0,
        center_of_mass=d_cm,
        total_momentum=np.zeros(3, dtype=np.float64),
    )
    updated_state = state.with_updates(
        q=state.q + delta_state.q,
        p=state.p + delta_state.p,
        t=state.t + dt,
    )

    # Preserve pre-drift internal relative coordinates of each supernode cluster S_k
    # around its drifted cluster barycenter (R_{S_k}^{after}, P_{S_k}^{after}) so that
    # steps 1-4 translate only the cluster barycenter even if S_k is not a single node of tree.
    if clusters:
        q_adj = np.copy(updated_state.q)
        p_adj = np.copy(updated_state.p)
        for c_idx in clusters:
            if len(c_idx) >= 3:
                idx_list = list(c_idx)
                m_c = state.masses[idx_list]
                M_c = float(np.sum(m_c))
                R_c_before = np.sum(m_c[:, None] * state.q[idx_list], axis=0) / M_c
                P_c_before = np.sum(state.p[idx_list], axis=0)
                R_c_after = np.sum(m_c[:, None] * updated_state.q[idx_list], axis=0) / M_c
                P_c_after = np.sum(updated_state.p[idx_list], axis=0)
                for b_i in idx_list:
                    q_adj[b_i] = R_c_after + (state.q[b_i] - R_c_before)
                    p_adj[b_i] = (state.masses[b_i] / M_c) * P_c_after + (
                        state.p[b_i] - (state.masses[b_i] / M_c) * P_c_before
                    )
        updated_state = updated_state.with_updates(q=q_adj, p=p_adj)

    # 5. Advance synchronous supernode clusters individually if present
    # Note: Center of mass translation was already applied in jacobi_to_cartesian, so drift_center_of_mass=False
    active_kepler = kepler_solver if kepler_solver is not None else get_current_kepler_solver()
    updated_state = step_supernode_flow(
        updated_state,
        dt,
        G=G,
        supernode_clusters=clusters,
        supernode_bodies=sn_set,
        drift_center_of_mass=False,
        kepler_solver=active_kepler,
        drift_weights=drift_weights,
        impulse_weights=impulse_weights,
    )

    return updated_state


def step_supernode_flow(
    state: ExtendedPhaseState,
    dt: float,
    G: float = G_DEFAULT,
    supernode_clusters: Optional[List[Tuple[int, ...]]] = None,
    supernode_bodies: Optional[Set[int]] = None,
    drift_center_of_mass: bool = False,
    kepler_solver: Optional[Any] = None,
    drift_weights: Optional[Sequence[float]] = None,
    impulse_weights: Optional[Sequence[float]] = None,
) -> ExtendedPhaseState:
    """
    Executes the coupled symplectic Hamiltonian flow of synchronous supernode clusters (k >= 3)
    under local unfolded conormal Sundman time:
        dt_S = chi_S^2 dtau_S = rho_S dtau_S.

    Mathematical Monograph Alignment:
    - Chapter 4: Section 4.3 (Synchronous Supernodes and Anti-Zeno Regularization),
      Theorem 4.3 and Equations (4.55)-(4.65).
    - Preserves total cluster linear momentum (< 1e-15) and angular momentum.
    - Note on Center-of-Mass Kinematics:
      When called during the drift phase of step_drift_and_kepler, the global barycenter
      drift has already been executed during jacobi_to_cartesian; drift_center_of_mass is
      therefore set to False to avoid duplicate barycentric translation.
    """
    if dt == 0.0 or state.N <= 2:
        return state

    prof = get_current_profiler()

    def _run_supernode_update(st: ExtendedPhaseState) -> ExtendedPhaseState:
        updated_state = st
        clusters = _normalize_supernode_clusters(supernode_bodies, supernode_clusters)
        if not clusters:
            return updated_state

        active_kepler = kepler_solver if kepler_solver is not None else get_current_kepler_solver()
        active_supernode = get_current_supernode_solver()

        if active_supernode is not None and hasattr(active_supernode, "step_supernode_synchronous"):
            step_sn_fn = active_supernode.step_supernode_synchronous
        elif callable(active_supernode):
            step_sn_fn = active_supernode
        else:
            step_sn_fn = step_supernode_synchronous

        active_stages = get_current_integrator_stages()
        dw = drift_weights if drift_weights is not None else (active_stages[0] if active_stages is not None else None)
        iw = impulse_weights if impulse_weights is not None else (active_stages[1] if active_stages is not None else None)

        extra_kwargs: Dict[str, Any] = {}
        active_n_sub = get_current_supernode_n_substeps()
        if (dw is not None and iw is not None) or active_n_sub is not None:
            import inspect
            try:
                sig = inspect.signature(step_sn_fn)
                if dw is not None and iw is not None and "drift_weights" in sig.parameters:
                    extra_kwargs["drift_weights"] = dw
                    extra_kwargs["impulse_weights"] = iw
                if active_n_sub is not None and "n_substeps" in sig.parameters:
                    extra_kwargs["n_substeps"] = active_n_sub
            except (ValueError, TypeError):
                pass

        for sn_indices in clusters:
            if len(sn_indices) >= 3:
                sn_q, sn_p = step_sn_fn(
                    updated_state.q,
                    updated_state.p,
                    updated_state.masses,
                    sn_indices,
                    dt,
                    G=G,
                    drift_center_of_mass=drift_center_of_mass,
                    kepler_solver=active_kepler,
                    **extra_kwargs,
                )
                updated_state = updated_state.with_updates(q=sn_q, p=sn_p)
        return updated_state


    if prof is not None:
        with prof.time_stage("supernode"):
            return _run_supernode_update(state)
    return _run_supernode_update(state)


# Canonical alias for compatibility with operator-splitting kick/flow wrappers
step_supernode_kick = step_supernode_flow


def step_apex_kick(
    state: ExtendedPhaseState,
    tree: JacobiTree,
    apex: ScreenedApexPotential,
    dt: float,
    supernode_bodies: Optional[Set[int]] = None,
    supernode_clusters: Optional[List[Tuple[int, ...]]] = None,
) -> ExtendedPhaseState:
    """
    Executes the vertical cotangent shear kick:
        p_i <- p_i - dt * nabla_{q_i} W_{apex}^{screen}(q).

    Mathematical Monograph Alignment:
    - Chapter 4: Section 4.2 (Conormal Fulton-MacPherson Splitting), Theorem 4.1 & 4.2.
    - Chapter 5: Section 5.5, Equations (5.82)-(5.84).
    - Chapter 6: Section 6.6.4 (ContextVar Modular Dispatch), Equation (6.132).
    """
    if dt == 0.0 or state.N <= 1:
        return state
    active_apex = get_current_apex_evaluator() or apex
    prof = get_current_profiler()

    def _call_apex() -> ExtendedPhaseState:
        try:
            return active_apex.step_apex_kick(
                state,
                tree,
                dt=dt,
                supernode_bodies=supernode_bodies,
                supernode_clusters=supernode_clusters,
            )
        except TypeError:
            return active_apex.step_apex_kick(
                state,
                tree,
                dt=dt,
                supernode_bodies=supernode_bodies,
            )

    if prof is not None:
        with prof.time_stage("apex"):
            res = _call_apex()
            prof.record_apex(kicks=1, pairs=len(tree.leaf_pairs) if tree else 0)
        return res
    return _call_apex()


step_screened_apex_kick = step_apex_kick
