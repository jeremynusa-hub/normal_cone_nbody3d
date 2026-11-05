"""
Unit Tests for Layer 1 State (normal_cone_nbody_float64.core.state).

Mathematical Verification Objectives:
- Immutability / copy-on-write integrity of ExtendedPhaseState (frozen dataclass, read-only arrays).
- Validation of dimensional consistency ((N, 3), m_i > 0, N >= 2).
- Regular locus predicate: regular_iff_not_collision (NonCoincident q <-> not InCollisionLocus q).
- Exact center-of-mass pinning and zero momentum projection via to_barycentric().
"""

import pytest
import numpy as np

from normal_cone_nbody_float64.core.state import ExtendedPhaseState


@pytest.fixture
def three_body_state() -> ExtendedPhaseState:
    # 3-body Lagrange equilateral configuration
    q = np.array([
        [1.0, 0.0, 0.0],
        [-0.5, np.sqrt(3)/2, 0.0],
        [-0.5, -np.sqrt(3)/2, 0.0],
    ], dtype=np.float64)
    p = np.array([
        [0.0, 0.5, 0.0],
        [-np.sqrt(3)/4, -0.25, 0.0],
        [np.sqrt(3)/4, -0.25, 0.0],
    ], dtype=np.float64)
    masses = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E0=-1.5)


def test_state_basic_properties(three_body_state):
    state = three_body_state
    assert state.N == 3
    assert state.total_mass == pytest.approx(6.0, rel=1e-15)
    assert state.E == -1.5
    assert state.E0 == -1.5


def test_immutability(three_body_state):
    state = three_body_state
    # Dataclass is frozen
    with pytest.raises(Exception):
        state.t = 1.0  # type: ignore

    # Underlying numpy arrays have writeable=False
    with pytest.raises(ValueError):
        state.q[0, 0] = 999.0
    with pytest.raises(ValueError):
        state.p[0, 0] = 999.0


def test_validation_errors():
    # Negative mass error
    with pytest.raises(ValueError, match="strictly positive"):
        ExtendedPhaseState(
            q=np.zeros((2, 3)),
            p=np.zeros((2, 3)),
            masses=np.array([1.0, -0.5]),
        )

    # Shape mismatch error
    with pytest.raises(ValueError, match="shape"):
        ExtendedPhaseState(
            q=np.zeros((3, 2)),
            p=np.zeros((3, 3)),
            masses=np.array([1.0, 1.0, 1.0]),
        )

    # Body count < 2
    with pytest.raises(ValueError, match="at least N >= 2"):
        ExtendedPhaseState(
            q=np.zeros((1, 3)),
            p=np.zeros((1, 3)),
            masses=np.array([1.0]),
        )


def test_with_updates(three_body_state):
    state = three_body_state
    new_masses = np.array([2.0, 3.0, 5.0], dtype=np.float64)
    new_state = state.with_updates(t=2.5, E=-1.2, masses=new_masses)
    assert new_state.t == 2.5
    assert new_state.E == -1.2
    assert new_state.total_mass == pytest.approx(10.0, rel=1e-15)
    np.testing.assert_allclose(new_state.masses, new_masses)
    assert state.t == 0.0  # Original remains unchanged
    assert state.total_mass == pytest.approx(6.0, rel=1e-15)
    np.testing.assert_allclose(new_state.q, state.q)


def test_barycentric_projection(three_body_state):
    # Add arbitrary center-of-mass shift
    shifted_q = three_body_state.q + np.array([10.0, -5.0, 3.0])
    shifted_p = three_body_state.p + np.array([-2.0, 1.0, 0.5])
    shifted_state = three_body_state.with_updates(q=shifted_q, p=shifted_p)

    bary_state = shifted_state.to_barycentric()
    np.testing.assert_allclose(bary_state.center_of_mass, np.zeros(3), atol=1e-15)
    np.testing.assert_allclose(bary_state.total_momentum, np.zeros(3), atol=1e-15)


def test_regular_locus_predicate():
    # Coincident bodies
    q_collision = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ])
    state_col = ExtendedPhaseState(
        q=q_collision,
        p=np.zeros((3, 3)),
        masses=np.ones(3),
    )
    assert not state_col.is_collision_free()

    # Separated bodies
    q_reg = np.array([
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ])
    state_reg = ExtendedPhaseState(
        q=q_reg,
        p=np.zeros((3, 3)),
        masses=np.ones(3),
    )
    assert state_reg.is_collision_free()
    assert state_reg.is_collision_free(tol=1e-3)
    assert not state_reg.is_collision_free(tol=0.2)
