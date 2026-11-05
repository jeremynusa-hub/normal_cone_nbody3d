"""
Unit Tests for Layer 7: Astrophysical Benchmark Registry (registry.py).

Theoretical Alignments:
- Monograph: Chapter 6 and Chapter 7 (Canonical Experimental Benchmarks).
- Invariants: Exact dimensional compliance (N >= 2, 3D), positive mass vectors,
  and initial center-of-mass momentum conservation ||P_0|| < 1e-14.
"""

import numpy as np
import pytest

from normal_cone_nbody_float64.benchmarks.registry import (
    BENCHMARK_REGISTRY,
    get_benchmark,
)
from normal_cone_nbody_float64.core.state import ExtendedPhaseState


def test_registry_contains_all_eight_canonical_benchmarks():
    """
    Validates that the registry enumerates all 15 canonical benchmarks:
    1. radial_collision
    2. burrau_pythagorean
    3. arenstorf_lunar
    4. figure8_choreography
    5. sitnikov_oscillatory
    6. quadruple_collision
    7. jovian_solar_system
    8. pareto_test
    9. eccentric_binary
    10. jovian_laplace_resonance
    11. comet_lexell_1770
    12. shoemaker_levy_9_1992
    13. apophis_2029_flyby
    14. hierarchical_3d_eccentric_triple
    15. inner_solar_system_mercury
    """
    expected = {
        "radial_collision",
        "burrau_pythagorean",
        "arenstorf_lunar",
        "figure8_choreography",
        "sitnikov_oscillatory",
        "quadruple_collision",
        "jovian_solar_system",
        "pareto_test",
        "eccentric_binary",
        "jovian_laplace_resonance",
        "comet_lexell_1770",
        "shoemaker_levy_9_1992",
        "apophis_2029_flyby",
        "hierarchical_3d_eccentric_triple",
        "inner_solar_system_mercury",
    }
    assert set(BENCHMARK_REGISTRY.keys()) == expected


def test_historical_encounters_and_laskar_mercury_benchmarks():
    """
    Validates the 4 real documented astronomical encounter & Laskar planetary benchmarks:
    1. comet_lexell_1770: N=4 (Sun, Jupiter, Earth, Comet D/1770 L1 Lexell) starting outside
       Jupiter's Hill switching sphere at r0_j = 0.70 AU and dynamically flipping
       T_helio -> T_Jovi -> T_helio through deep perijove r_min ~ 3.10 R_J.
    2. shoemaker_levy_9_1992: N=5 (Sun, Jupiter + 3 self-gravitating comet progenitor
       fragments) approaching Roche-limit tidal disruption.
    3. apophis_2029_flyby: N=4 (Sun, Earth, Moon, Asteroid 99942 Apophis) grazing
       Earth encounter inside the geostationary belt deflecting Apophis from Aten
       (a_pre = 0.9224 AU < 1 AU) to Apollo (a_post ~ 1.116 AU > 1 AU, Delta_a > 0).
    4. inner_solar_system_mercury: N=5 (Sun, Mercury, Venus, Earth, Mars) with
       Mercury's high eccentricity e_Merc = 0.205630 (Laskar / Farrés 2013 bottleneck).
    """
    from normal_cone_nbody_float64.system.nbody_system import NBodySystem
    from normal_cone_nbody_float64.integrators.processed import create_pbm4_integrator
    from normal_cone_nbody_float64.diagnostics.metrics import compute_osculating_elements

    s_lex = get_benchmark("comet_lexell_1770")
    assert s_lex.N == 4
    # Comet (body 3) starts at r0_j = 0.70 AU outside Jupiter's Hill switching boundary (~0.581 AU)
    r_jup_lex = np.linalg.norm(s_lex.q[3] - s_lex.q[1])
    assert np.isclose(r_jup_lex, 0.70, rtol=1e-4)

    # Verify dynamic chart flip T_helio -> T_Jovi -> T_helio and deep perijove ~ 3.10 R_J
    sys_lex = NBodySystem(state=s_lex, integrator=create_pbm4_integrator(G=1.0), G=1.0)
    sys_lex.scheduler.max_level = 6
    assert sys_lex.tree.get_leaf_pair_indices() == [(0, 2)]  # Heliocentric start
    min_rj = r_jup_lex
    seen_jovi = False
    while sys_lex.t < 12.5:
        dt_step = 0.002 if 5.39 <= sys_lex.t < 5.46 else 0.10
        sys_lex.step(dt=dt_step)
        st_phys = sys_lex.get_physical_state()
        d_j = float(np.linalg.norm(st_phys.q[3] - st_phys.q[1]))
        min_rj = min(min_rj, d_j)
        if (1, 3) in sys_lex.tree.get_leaf_pair_indices():
            seen_jovi = True
    assert seen_jovi is True, "Comet Lexell failed to trigger dynamic flip to Jovicentric chart"
    assert sys_lex.tree.get_leaf_pair_indices() == [(0, 2)], "Comet Lexell failed to flip back to heliocentric chart after exit"
    assert len(sys_lex.atlas_manager.transition_history) >= 2
    assert min_rj < 0.0020, f"Expected deep perijove r_min ~ 0.00148 AU, got {min_rj}"

    s_sl9 = get_benchmark("shoemaker_levy_9_1992")
    assert s_sl9.N == 5
    # Progenitor fragments (2, 3, 4) form a compact cluster near Jupiter (body 1)
    r_frag_12 = np.linalg.norm(s_sl9.q[3] - s_sl9.q[2])
    r_jup_sl9 = np.linalg.norm(s_sl9.q[2] - s_sl9.q[1])
    assert r_frag_12 < 0.01 * r_jup_sl9

    s_apo = get_benchmark("apophis_2029_flyby")
    assert s_apo.N == 4
    r_earth_apo = np.linalg.norm(s_apo.q[3] - s_apo.q[1])
    assert np.isclose(r_earth_apo, 0.0020, rtol=1e-2)
    elem_pre = compute_osculating_elements(s_apo, primary_idx=0, body_idx=3, G=1.0)
    assert np.isclose(elem_pre["a"], 0.9224, atol=1e-3), f"Expected Aten a_pre ~ 0.9224 AU, got {elem_pre['a']}"
    assert elem_pre["a"] < 1.0

    # Integrate Apophis across the 2029 perigee flyby and verify transition to Apollo (a_post > 1.0 AU)
    sys_apo = NBodySystem(state=s_apo, integrator=create_pbm4_integrator(G=1.0), G=1.0)
    sys_apo.scheduler.max_level = 6
    sys_apo.evolve(t_span=0.02, dt=0.001)
    elem_post = compute_osculating_elements(sys_apo.get_physical_state(), primary_idx=0, body_idx=3, G=1.0)
    assert elem_post["a"] > 1.05, f"Expected Apollo a_post > 1.05 AU, got {elem_post['a']}"
    assert elem_post["a"] - elem_pre["a"] > 0.15

    s_inner = get_benchmark("inner_solar_system_mercury")
    assert s_inner.N == 5
    # Verify Mercury (body 1) starts at perihelion r_peri = a*(1-e) = 0.387098*(1-0.205630)
    r_merc = np.linalg.norm(s_inner.q[1] - s_inner.q[0])
    assert np.isclose(r_merc, 0.387098 * (1.0 - 0.205630), rtol=1e-6)


def test_jovian_laplace_resonance_phase_invariant():
    """
    Validates the 4-body Jovian Laplace resonance benchmark (Jupiter, Io, Europa, Ganymede):
    - N = 4 bodies.
    - Forced eccentricities e_1 = 0.0041, e_2 = 0.0094, e_3 = 0.0006 at pericenter.
    - Laplace resonant phase theta_Lap(t) = lambda_1 - 3*lambda_2 + 2*lambda_3 = pi (mod 2*pi)
      in osculating mean longitudes, librating tightly around 180 deg (< 1 deg amplitude).
    """
    from normal_cone_nbody_float64.diagnostics.metrics import (
        compute_osculating_elements,
        compute_laplace_resonant_angle,
    )
    from normal_cone_nbody_float64.system.nbody_system import NBodySystem
    from normal_cone_nbody_float64.integrators.blanes_moan_4 import BlanesMoan4Integrator

    state = get_benchmark("jovian_laplace_resonance")
    assert state.N == 4

    el1 = compute_osculating_elements(state, primary_idx=0, body_idx=1, G=1.0)
    el2 = compute_osculating_elements(state, primary_idx=0, body_idx=2, G=1.0)
    el3 = compute_osculating_elements(state, primary_idx=0, body_idx=3, G=1.0)
    assert np.isclose(el1["e"], 0.0041, atol=1e-4)
    assert np.isclose(el2["e"], 0.0094, atol=1e-4)
    assert np.isclose(el3["e"], 0.0006, atol=1e-4)

    lap0_rad = compute_laplace_resonant_angle(state, G=1.0)
    lap0_deg = float(np.degrees(lap0_rad))
    assert np.isclose(lap0_deg, 180.0, atol=1e-9)

    # Evolve for ~20 units and verify tight libration around 180 deg without circulation
    sys_lap = NBodySystem(state=state, integrator=BlanesMoan4Integrator(G=1.0), G=1.0)
    angles_deg = [lap0_deg]
    for _ in range(40):
        sys_lap.step(dt=0.5)
        ang_rad = compute_laplace_resonant_angle(sys_lap.get_physical_state(), G=1.0)
        angles_deg.append(float(np.degrees(ang_rad)))
    max_dev_deg = max(abs(a - 180.0) for a in angles_deg)
    assert max_dev_deg < 2.0, f"Laplace resonant angle drifted from 180 deg by {max_dev_deg:.3f} deg"



@pytest.mark.parametrize("bench_name", list(BENCHMARK_REGISTRY.keys()))
def test_benchmark_state_validity_and_barycenter_integrity(bench_name: str):
    """
    Validates that every registered benchmark creates a structurally flawless
    ExtendedPhaseState satisfying:
    - Positive masses m_i > 0.
    - Exact dimensions (N, 3).
    - Null initial total momentum ||P_0|| < 1e-14 (isolated system conservation).
    - Finite Hamiltonian energy.
    """
    state = get_benchmark(bench_name)

    assert isinstance(state, ExtendedPhaseState)
    state.validate()

    assert state.N >= 2
    assert state.q.shape == (state.N, 3)
    assert state.p.shape == (state.N, 3)
    assert np.all(state.masses > 0.0)

    # Initial momentum zero check
    p_norm = np.linalg.norm(state.total_momentum)
    assert p_norm < 1e-14, f"Benchmark '{bench_name}' has non-zero total momentum: {p_norm}"

    # Energy finite check
    assert np.isfinite(state.E)
    assert np.isfinite(state.E0)
