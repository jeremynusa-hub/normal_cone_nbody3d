"""
Unit tests for OptimizedScreenedApexPotential.

Verifies:
1. Exact bit-exact force concordance with baseline ScreenedApexPotential (< 1e-15 max diff).
2. Strict Lie-Bargmann linear and angular momentum conservation (< 1e-15).
3. Exact equivalence of symplectic kick operator step_apex_kick.
4. Compatibility across 3-body and N-body configurations.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.core.constants import GRAVITATIONAL_CONSTANT
from normal_cone_nbody_float64.benchmarks.registry import (
    create_burrau_pythagorean_benchmark,
    create_figure8_choreography_benchmark,
    create_radial_collision_benchmark,
)
from normal_cone_nbody_float64.system.atlas_manager import AtlasManager
from normal_cone_nbody_float64.potential.screened_apex import ScreenedApexPotential
from normal_cone_nbody_float64.potential.optimized_screened_apex import (
    OptimizedScreenedApexPotential,
)


@pytest.mark.parametrize("benchmark_factory", [
    create_burrau_pythagorean_benchmark,
    create_figure8_choreography_benchmark,
    create_radial_collision_benchmark,
])
def test_optimized_apex_forces_concordance(benchmark_factory):
    state = benchmark_factory()
    atlas = AtlasManager()
    tree = atlas.build_optimal_tree(state)

    base_apex = ScreenedApexPotential(G=GRAVITATIONAL_CONSTANT)
    opt_apex = OptimizedScreenedApexPotential(G=GRAVITATIONAL_CONSTANT)

    # 1. Pure kick force (include_screened_monopole=False) — 0.0 ULP exact match
    f_base, v_base = base_apex.compute_apex_forces(state, tree, include_screened_monopole=False)
    f_opt, v_opt = opt_apex.compute_apex_forces(state, tree, include_screened_monopole=False)

    assert np.array_equal(f_base, f_opt), f"Force mismatch: {f_base - f_opt}"
    assert v_base == v_opt, f"Potential mismatch: {v_base} vs {v_opt}"

    # 2. Lie-Bargmann linear momentum conservation: sum F_i = 0
    total_force = np.sum(f_opt, axis=0)
    assert np.all(np.abs(total_force) < 1e-15), f"Non-zero net force: {total_force}"

    # 3. Lie-Bargmann angular momentum conservation: sum q_i x F_i = 0
    total_torque = np.sum(np.cross(state.q, f_opt), axis=0)
    assert np.all(np.abs(total_torque) < 1e-15), f"Non-zero net torque: {total_torque}"


def test_optimized_step_apex_kick():
    state = create_burrau_pythagorean_benchmark()
    atlas = AtlasManager()
    tree = atlas.build_optimal_tree(state)

    base_apex = ScreenedApexPotential(G=GRAVITATIONAL_CONSTANT)
    opt_apex = OptimizedScreenedApexPotential(G=GRAVITATIONAL_CONSTANT)

    dt = 0.005
    st_base = base_apex.step_apex_kick(state, tree, dt=dt)
    st_opt = opt_apex.step_apex_kick(state, tree, dt=dt)

    assert np.array_equal(st_base.p, st_opt.p), f"Kick p mismatch: {st_base.p - st_opt.p}"


def test_optimized_composite_screening():
    state = create_burrau_pythagorean_benchmark()
    atlas = AtlasManager()
    tree = atlas.build_optimal_tree(state)

    base_apex = ScreenedApexPotential(G=GRAVITATIONAL_CONSTANT, screening_mode="c6")
    opt_apex = OptimizedScreenedApexPotential(G=GRAVITATIONAL_CONSTANT, screening_mode="c6")

    f_base, v_base = base_apex.compute_apex_forces(state, tree, include_screened_monopole=True)
    f_opt, v_opt = opt_apex.compute_apex_forces(state, tree, include_screened_monopole=True)

    assert np.array_equal(f_base, f_opt), f"Force mismatch with composite screening: {f_base - f_opt}"
    assert v_base == v_opt, f"Potential mismatch with composite screening: {v_base} vs {v_opt}"


def test_optimized_apex_n5_and_gc_cache_safety_zero_ulp():
    """
    Validates 0.0 ULP concordance for N=5 (both standard and composite screening)
    and verifies that _mass_cache and _topo_cache retain strong references so that
    ephemeral arrays/trees created in a loop never collide via id() reuse after GC.
    """
    from normal_cone_nbody_float64.core.state import ExtendedPhaseState
    from normal_cone_nbody_float64.core.jacobi import build_canonical_tree

    rng = np.random.default_rng(42)
    base_apex = ScreenedApexPotential(G=1.0, screening_mode="c6")
    opt_apex = OptimizedScreenedApexPotential(G=1.0, screening_mode="c6")

    for _ in range(10):
        masses_5 = rng.uniform(0.1, 2.0, size=5)
        pos_5 = rng.uniform(-5.0, 5.0, size=(5, 3))
        vel_5 = rng.uniform(-0.5, 0.5, size=(5, 3))
        state_5 = ExtendedPhaseState(q=pos_5, p=vel_5 * masses_5[:, None], masses=masses_5)
        tree_5 = build_canonical_tree(masses_5, topology="comb")

        f_b, v_b = base_apex.compute_apex_forces(state_5, tree_5, include_screened_monopole=True)
        f_o, v_o = opt_apex.compute_apex_forces(state_5, tree_5, include_screened_monopole=True)
        assert np.array_equal(f_b, f_o)
        assert v_b == v_o


def test_optimized_apex_coincident_coordinates_zero_division_guard():
    """
    Validates [ENG-FALSE-06]:
    OptimizedScreenedApexPotential.compute_apex_forces guards against r^2 == 0.0 in both
    the specialized N=3 branch and the general N > 3 vectorized branch without emitting
    RuntimeWarning or producing NaN / Inf.
    """
    import warnings
    from normal_cone_nbody_float64.core.state import ExtendedPhaseState
    from normal_cone_nbody_float64.core.jacobi import build_canonical_tree

    opt_apex = OptimizedScreenedApexPotential(G=1.0)

    # 1. N=3 coincident non-leaf pair (bodies 0 and 2 at identical position)
    m3 = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    q3 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float64)
    st3 = ExtendedPhaseState(q=q3, p=np.zeros_like(q3), masses=m3)
    tree3 = build_canonical_tree(m3, topology="comb")

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        f3, v3 = opt_apex.compute_apex_forces(st3, tree3)
    assert np.all(np.isfinite(f3)) and np.isfinite(v3)

    # 2. N=4 coincident non-leaf pair (bodies 0 and 3 at identical position)
    m4 = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    q4 = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)
    st4 = ExtendedPhaseState(q=q4, p=np.zeros_like(q4), masses=m4)
    tree4 = build_canonical_tree(m4, topology="comb")

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        f4, v4 = opt_apex.compute_apex_forces(st4, tree4)
    assert np.all(np.isfinite(f4)) and np.isfinite(v4)

