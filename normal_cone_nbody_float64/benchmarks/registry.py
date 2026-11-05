r"""
Canonical Astrophysical N-Body Benchmark Registry.

Theoretical Foundations:
- Monograph Reference: Chapter 7, Section 7.1 (\label{sec:unsoftened_benchmarks})
  and Section 7.2 (\label{sec:long_term_stability}); Chapter 3, Section 3.3 (\label{sec:arenstorf_validation}).
- The 8 Canonical Configurations:
    1. Radial Two-Body Collision (L = 0, e = 1.0, exact Stumpff-Sundman regularization).
    2. Burrau's Pythagorean 3-Body Problem (masses 3-4-5, extreme close encounters).
    3. Arenstorf Periodic Lunar Orbit (PCR3BP, close lunar flyby).
    4. Chenciner--Montgomery Figure-8 Choreography (N = 3, equal masses, zero momentum).
    5. Sitnikov Oscillatory 3-Body Problem (N = 3, vertical axis oscillation).
    6. Symmetric Quadruple Central Collision (N = 4, D_4 symmetry).
    7. Outer Solar System (Sun-Jupiter-Saturn, long-term quasi-secular stability).
    8. Pareto Work-Precision Convergence Benchmark (smooth non-colliding chaotic 3-body).
- Invariants: Zero initial total momentum ||P_0|| < 1e-15, center of mass R_cm = 0,
  and exact unsoftened geometry (eps = 0, d0 = 0).
"""

from __future__ import annotations
import math
from typing import Dict, Callable, Any, Optional
import numpy as np

from ..core.constants import GRAVITATIONAL_CONSTANT
from ..core.state import ExtendedPhaseState
from ..core.invariants import compute_cartesian_hamiltonian


def create_radial_collision_benchmark(
    r0: float = 2.0,
    m0: float = 1.0,
    m1: float = 1.0,
    G: float = GRAVITATIONAL_CONSTANT,
) -> ExtendedPhaseState:
    """
    Two-body head-on radial collision from rest (angular momentum L = 0, e = 1.0).
    Impact occurs at t_coll = pi * r0^{3/2} / (2 * sqrt(2 * G * (m0 + m1))).
    """
    masses = np.array([float(m0), float(m1)], dtype=np.float64)
    # Positions along the X-axis centered at barycenter
    # m0 * x0 + m1 * x1 = 0, x1 - x0 = r0
    x0 = -r0 * (m1 / (m0 + m1))
    x1 = r0 * (m0 / (m0 + m1))

    q = np.array([
        [x0, 0.0, 0.0],
        [x1, 0.0, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q)

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_burrau_pythagorean_benchmark(
    G: float = GRAVITATIONAL_CONSTANT,
) -> ExtendedPhaseState:
    """
    Burrau's Pythagorean 3-Body Problem:
    Masses 3, 4, 5 placed at the vertices of a 3-4-5 right triangle at rest.
    Undergoes multiple close triple encounters before asymptotic ejection.
    """
    masses = np.array([3.0, 4.0, 5.0], dtype=np.float64)
    # Standard Pythagorean vertices:
    # Body 0 (m=3): (1, 3, 0)
    # Body 1 (m=4): (-2, -1, 0)
    # Body 2 (m=5): (1, -1, 0)
    # Barycenter is already exactly (0, 0, 0)
    q = np.array([
        [1.0, 3.0, 0.0],
        [-2.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q)

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_arenstorf_lunar_benchmark(
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Arenstorf periodic lunar orbit (PCR3BP in barycentric inertial frame).
    Primary Earth (m0 = 1 - mu), Moon (m1 = mu), mu = 0.012277471.
    Third body: satellite of mass m2 = 1e-12 on periodic orbit.
    """
    mu = 0.012277471
    m0 = 1.0 - mu
    m1 = mu
    m2 = 1e-12
    masses = np.array([m0, m1, m2], dtype=np.float64)

    # Earth-Moon circular orbit on X-axis (separation 1.0, angular velocity omega = 1.0)
    x_earth = -mu
    x_moon = 1.0 - mu
    v_earth = -mu * 1.0
    v_moon = (1.0 - mu) * 1.0

    # Arenstorf satellite initial state (inertial frame at t=0, rotating x=0.994, y=0)
    # Rotating velocity (u=0, v=-2.0015851063790825224)
    # Inertial velocity: V_x = u - omega * y = 0, V_y = v + omega * x = -2.001585106379 + 0.994
    v_arenstorf_y = -2.0015851063790825 + 0.994

    q = np.array([
        [x_earth, 0.0, 0.0],
        [x_moon, 0.0, 0.0],
        [0.994, 0.0, 0.0],
    ], dtype=np.float64)

    v = np.array([
        [0.0, -mu, 0.0],
        [0.0, 1.0 - mu, 0.0],
        [0.0, v_arenstorf_y, 0.0],
    ], dtype=np.float64)

    p = v * masses[:, np.newaxis]
    # Neutralize small center of mass velocity
    p_cm = np.sum(p, axis=0) / np.sum(masses)
    p -= masses[:, np.newaxis] * p_cm

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_figure8_choreography_benchmark(
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Chenciner & Montgomery Figure-8 Choreography for three equal masses.
    High-precision initial conditions (Simó 2000, 16 digits).
    Period T approx 6.32591398.
    """
    masses = np.array([1.0, 1.0, 1.0], dtype=np.float64)

    # Initial positions (Simó 2000, 16 digits)
    x1 = 0.9700043566977735
    y1 = -0.2430875317772714
    q = np.array([
        [x1, y1, 0.0],
        [-x1, -y1, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)

    # Initial velocities (Simó 2000, 16 digits)
    vx3 = -0.9324073714441582
    vy3 = -0.8647314614271871
    vx1 = -0.5 * vx3
    vy1 = -0.5 * vy3
    vx2 = vx1
    vy2 = vy1

    v = np.array([
        [vx1, vy1, 0.0],
        [vx2, vy2, 0.0],
        [vx3, vy3, 0.0],
    ], dtype=np.float64)

    p = v * masses[:, np.newaxis]
    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_sitnikov_oscillatory_benchmark(
    z0: float = 0.5,
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Sitnikov 3-body problem:
    Two equal primary masses m0 = m1 = 0.5 orbiting circularly in the XY plane.
    Infinitesimal particle m2 = 1e-12 oscillating along the Z axis.
    """
    masses = np.array([0.5, 0.5, 1e-12], dtype=np.float64)

    # Primaries circular orbit with separation r = 1.0
    # Orbital speed v = sqrt(G * M_prim / (4 * r_prim)) = sqrt(1.0 * 1.0 / 4) = 0.5
    q = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 0.0, float(z0)],
    ], dtype=np.float64)

    v = np.array([
        [0.0, -0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float64)

    p = v * masses[:, np.newaxis]
    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_quadruple_collision_benchmark(
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Symmetric Quadruple Central Collision (N = 4):
    Four equal masses placed symmetrically on the axes (+-1, 0) and (0, +-1) at rest.
    Collapses simultaneously to origin with exact D_4 point group symmetry.
    """
    masses = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)

    q = np.array([
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
    ], dtype=np.float64)
    p = np.zeros_like(q)

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_jovian_solar_system_benchmark(
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Outer Solar System (Sun-Jupiter-Saturn, N = 3):
    Long-term quasi-secular resonance test.
    Sun (m0 = 1.0), Jupiter (m1 = 9.547919e-4), Saturn (m2 = 2.858856e-4).
    """
    m_sun = 1.0
    m_jup = 9.547919e-4
    m_sat = 2.858856e-4
    masses = np.array([m_sun, m_jup, m_sat], dtype=np.float64)

    # Simplified circular coplanar approximation (AU and yr / 2pi units)
    # Jupiter: a ~ 5.2 AU, v ~ 1 / sqrt(5.2)
    # Saturn: a ~ 9.5 AU, v ~ 1 / sqrt(9.5)
    r_jup = 5.204
    r_sat = 9.582
    v_jup = math.sqrt(G * m_sun / r_jup)
    v_sat = math.sqrt(G * m_sun / r_sat)

    q = np.array([
        [0.0, 0.0, 0.0],
        [r_jup, 0.0, 0.0],
        [-r_sat, 0.0, 0.0],
    ], dtype=np.float64)

    v = np.array([
        [0.0, 0.0, 0.0],
        [0.0, v_jup, 0.0],
        [0.0, -v_sat, 0.0],
    ], dtype=np.float64)

    # Adjust barycenter
    cdm = np.sum(masses[:, None] * q, axis=0) / np.sum(masses)
    q -= cdm
    p = v * masses[:, np.newaxis]
    p_cm = np.sum(p, axis=0) / np.sum(masses)
    p -= masses[:, np.newaxis] * p_cm

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_pareto_test_benchmark(
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Standard hierarchical non-colliding chaotic 3-body system
    optimized for Pareto work-precision convergence analysis across integrators.
    """
    masses = np.array([1.0, 1.0, 0.5], dtype=np.float64)
    q = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 5.0, 0.2],
    ], dtype=np.float64)
    p = np.array([
        [0.0, -0.65, 0.0],
        [0.0, 0.65, 0.0],
        [0.2, 0.0, -0.05],
    ], dtype=np.float64)

    # Exact barycenter correction
    cdm = np.sum(masses[:, None] * q, axis=0) / np.sum(masses)
    q -= cdm
    p_cm = np.sum(p, axis=0) / np.sum(masses)
    p -= masses[:, np.newaxis] * p_cm

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_eccentric_binary_benchmark(
    e: float = 0.9,
    a: float = 1.0,
    m0: float = 1.0,
    m1: float = 1.0,
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Unsoftened eccentric two-body Keplerian binary (e in [0, 1)).
    Initialized at apoapsis r_max = a * (1 + e).
    Relative velocity v_apo = sqrt(G * (m0 + m1) * (1 - e) / (a * (1 + e))).
    Exact zero momentum and barycenter at origin.
    """
    masses = np.array([float(m0), float(m1)], dtype=np.float64)
    m_tot = float(m0 + m1)
    r_max = float(a * (1.0 + e))
    v_apo = math.sqrt(float(G) * m_tot * (1.0 - e) / r_max)

    x0 = -r_max * (float(m1) / m_tot)
    x1 = r_max * (float(m0) / m_tot)

    q = np.array([
        [x0, 0.0, 0.0],
        [x1, 0.0, 0.0],
    ], dtype=np.float64)

    vy0 = -v_apo * (float(m1) / m_tot)
    vy1 = v_apo * (float(m0) / m_tot)

    v = np.array([
        [0.0, vy0, 0.0],
        [0.0, vy1, 0.0],
    ], dtype=np.float64)

    p = v * masses[:, np.newaxis]
    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_jovian_laplace_resonance_benchmark(
    G: float = 1.0,
    use_forced_equilibrium: bool = True,
) -> ExtendedPhaseState:
    """
    Resonant 4-body Jovian Galilean satellite system (Jupiter, Io, Europa, Ganymede)
    in 4:2:1 Laplace mean-motion resonance (Chapter 7, Section 7.2):
    - Jupiter (m0 = 1.0), Io (m1 = 4.70e-5), Europa (m2 = 2.53e-5), Ganymede (m3 = 7.80e-5).
    - When use_forced_equilibrium=True (default):
      * Initialized at the forced resonant Laplace equilibrium at pericenter
        (M_1 = M_2 = M_3 = 0, varpi_1 = 0, varpi_2 = pi, varpi_3 = 0 =>
        lambda_1 = 0, lambda_2 = pi, lambda_3 = 0 => theta_Lap = lambda_1 - 3*lambda_2 + 2*lambda_3 = pi)
        with physical forced eccentricities e_1 = 0.0041, e_2 = 0.0094, e_3 = 0.0006
        and resonant mean motions n_1 = 1.0, n_2 = 0.5*(1 - Delta_n),
        n_3 = 0.25*(1 - 3*Delta_n + delta_Lap) (Delta_n = 0.003634, delta_Lap = -0.000865)
        compensating mutual apsidal precession so theta_Lap(t) librates tightly around 180 deg.
    """
    m_jup = 1.0
    m_io = 4.70e-5
    m_eu = 2.53e-5
    m_ga = 7.80e-5
    masses = np.array([m_jup, m_io, m_eu, m_ga], dtype=np.float64)

    if use_forced_equilibrium:
        e1, e2, e3 = 0.0041, 0.0094, 0.0006
        delta_n = 0.003634
        delta_lap = -0.000865
        n1 = 1.0
        n2 = 0.5 * (1.0 - delta_n)
        n3 = 0.25 * (1.0 - 3.0 * delta_n + delta_lap)

        a1 = (float(G) * (m_jup + m_io) / (n1 * n1)) ** (1.0 / 3.0)
        a2 = (float(G) * (m_jup + m_eu) / (n2 * n2)) ** (1.0 / 3.0)
        a3 = (float(G) * (m_jup + m_ga) / (n3 * n3)) ** (1.0 / 3.0)

        r1 = a1 * (1.0 - e1)
        r2 = a2 * (1.0 - e2)
        r3 = a3 * (1.0 - e3)

        v1 = math.sqrt(float(G) * (m_jup + m_io) * (1.0 + e1) / r1)
        v2 = math.sqrt(float(G) * (m_jup + m_eu) * (1.0 + e2) / r2)
        v3 = math.sqrt(float(G) * (m_jup + m_ga) * (1.0 + e3) / r3)
    else:
        r1 = 1.0
        r2 = 2.0 ** (2.0 / 3.0)
        r3 = 4.0 ** (2.0 / 3.0)
        v1 = math.sqrt(float(G) * (m_jup + m_io) / r1)
        v2 = math.sqrt(float(G) * (m_jup + m_eu) / r2)
        v3 = math.sqrt(float(G) * (m_jup + m_ga) / r3)

    q = np.array([
        [0.0, 0.0, 0.0],
        [r1, 0.0, 0.0],
        [-r2, 0.0, 0.0],
        [r3, 0.0, 0.0],
    ], dtype=np.float64)

    v = np.array([
        [0.0, 0.0, 0.0],
        [0.0, v1, 0.0],
        [0.0, -v2, 0.0],
        [0.0, v3, 0.0],
    ], dtype=np.float64)

    cdm = np.sum(masses[:, None] * q, axis=0) / np.sum(masses)
    q -= cdm
    p = v * masses[:, np.newaxis]
    p_cm = np.sum(p, axis=0) / np.sum(masses)
    p -= masses[:, np.newaxis] * p_cm

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_comet_lexell_1770_benchmark(
    G: float = 1.0,
    r0_j: float = 0.70,
) -> ExtendedPhaseState:
    """
    Documented Historical Close Encounter 1: Comet D/1770 L1 (Lexell) 1779 Deep
    Perijove Ejection in the 4-body Sun-Jupiter-Earth-Comet system (Lexell 1776;
    Laplace 1805; Le Verrier 1844):
    - Body 0: Sun (m0 = 1.0)
    - Body 1: Jupiter (m1 = 9.5458e-4, a_J = 5.20336 AU)
    - Body 2: Earth (m2 = 3.0035e-6, a_E = 1.00000 AU)
    - Body 3: Comet D/1770 L1 (Lexell, m3 = 1.0e-14) starting outside Jupiter's
      Hill switching sphere at r0_j = 0.70 AU (> gamma_build * R_Hill,J ~ 0.581 AU)
      on a Coriolis-compensated hyperbolic approach aiming at deep perijove
      r_min ~ 0.00148 AU (~3.10 R_J, inside Io's orbit), dynamically triggering
      the chart transitions T_helio -> T_Jovi -> T_helio across the encounter.
    """
    m_sun = 1.0
    m_jup = 9.5458e-4
    m_earth = 3.0035e-6
    m_com = 1.0e-14
    masses = np.array([m_sun, m_jup, m_earth, m_com], dtype=np.float64)

    a_jup = 5.20336
    a_earth = 1.00000
    v_jup = math.sqrt(float(G) * (m_sun + m_jup) / a_jup)
    v_earth = math.sqrt(float(G) * (m_sun + m_earth) / a_earth)

    rp_j = 0.0015
    v_inf_j = 0.25 * v_jup
    mu_j = float(G) * m_jup
    v0_j = math.sqrt(v_inf_j * v_inf_j + 2.0 * mu_j / float(r0_j))
    vp_j = math.sqrt(v_inf_j * v_inf_j + 2.0 * mu_j / rp_j)
    if abs(float(r0_j) - 0.70) < 1e-9:
        # Coriolis/solar-tidal compensated aiming angle from r0_j = 0.70 AU to hit r_min = 0.0014815 AU (3.10 R_J)
        sin_gamma = -0.005767603157
    else:
        sin_gamma = min(0.99, (rp_j * vp_j) / (float(r0_j) * v0_j))
    cos_gamma = math.sqrt(1.0 - sin_gamma * sin_gamma)

    q = np.array([
        [0.0, 0.0, 0.0],
        [a_jup, 0.0, 0.0],
        [0.0, a_earth, 0.0],
        [a_jup + float(r0_j), 0.0, 0.001],
    ], dtype=np.float64)

    v = np.array([
        [0.0, 0.0, 0.0],
        [0.0, v_jup, 0.0],
        [-v_earth, 0.0, 0.0],
        [-v0_j * cos_gamma, v_jup + v0_j * sin_gamma, 0.0],
    ], dtype=np.float64)

    cdm = np.sum(masses[:, None] * q, axis=0) / np.sum(masses)
    q -= cdm
    p = v * masses[:, np.newaxis]
    p_cm = np.sum(p, axis=0) / np.sum(masses)
    p -= masses[:, np.newaxis] * p_cm

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_shoemaker_levy_9_1992_benchmark(
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Documented Historical Close Encounter 2: Comet D/1993 F2 (Shoemaker-Levy 9)
    July 7, 1992 Perijove Tidal Disruption inside Jupiter's Fluid Roche Limit
    (r_min ~ 1.33 R_J ~ 6.35e-4 AU, Sekanina, Chodas & Yeomans 1994; Asphaug & Benz 1996):
    - Body 0: Sun (m0 = 1.0)
    - Body 1: Jupiter (m1 = 9.5458e-4, a_J = 5.20336 AU)
    - Bodies 2, 3, 4: Self-gravitating 3-fragment rubble-pile comet progenitor
      (each m_f = 1.0e-11, initial inter-fragment separation delta_r = 2.0e-5 AU)
      approaching perijove r_p = 6.35e-4 AU from r0 = 0.012 AU on a highly eccentric
      bound Jovicentric orbit (e_J ~ 0.996).
    """
    m_sun = 1.0
    m_jup = 9.5458e-4
    m_frag = 1.0e-11
    masses = np.array([m_sun, m_jup, m_frag, m_frag, m_frag], dtype=np.float64)

    a_jup = 5.20336
    v_jup = math.sqrt(float(G) * (m_sun + m_jup) / a_jup)

    # Jovicentric eccentric bound orbit (a_sl9 = 0.16 AU, r_p = 6.35e-4 AU -> e = 1 - r_p/a_sl9)
    r0_j = 0.012
    rp_j = 6.35e-4
    a_sl9 = 0.16
    mu_j = float(G) * m_jup
    v0_j = math.sqrt(mu_j * (2.0 / r0_j - 1.0 / a_sl9))
    vp_j = math.sqrt(mu_j * (2.0 / rp_j - 1.0 / a_sl9))
    sin_gamma = min(0.99, (rp_j * vp_j) / (r0_j * v0_j))
    cos_gamma = math.sqrt(1.0 - sin_gamma * sin_gamma)

    delta_r = 2.0e-5
    v_bind = math.sqrt(float(G) * (3.0 * m_frag) / delta_r)

    xc = a_jup + r0_j
    yc = 0.0
    vxc = -v0_j * cos_gamma
    vyc = v_jup + v0_j * sin_gamma

    q = np.array([
        [0.0, 0.0, 0.0],
        [a_jup, 0.0, 0.0],
        [xc - delta_r, yc, 0.0],
        [xc, yc, 0.0],
        [xc + delta_r, yc, 0.0],
    ], dtype=np.float64)

    v = np.array([
        [0.0, 0.0, 0.0],
        [0.0, v_jup, 0.0],
        [vxc, vyc - 0.5 * v_bind, 0.0],
        [vxc, vyc, 0.0],
        [vxc, vyc + 0.5 * v_bind, 0.0],
    ], dtype=np.float64)

    cdm = np.sum(masses[:, None] * q, axis=0) / np.sum(masses)
    q -= cdm
    p = v * masses[:, np.newaxis]
    p_cm = np.sum(p, axis=0) / np.sum(masses)
    p -= masses[:, np.newaxis] * p_cm

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_apophis_2029_flyby_benchmark(
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Documented Real Celestial Encounter 3: Asteroid (99942) Apophis April 13, 2029
    Grazing Geocentric Flyby inside the Geostationary Belt in the 4-body
    Sun-Earth-Moon-Apophis system (Giorgini et al. 2008; Farnocchia et al. 2013):
    - Body 0: Sun (m0 = 1.0)
    - Body 1: Earth (m1 = 3.0035e-6, a_E = 1.00000 AU)
    - Body 2: Moon (m2 = 3.6943e-8, a_Moon = 0.0025695 AU around Earth)
    - Body 3: Asteroid (99942) Apophis (m3 = 4.0e-19) starting at r0 = 0.0020 AU
      from Earth on a trailing-side hyperbolic geocentric approach (phi_approach = 171.575217 deg)
      aiming at perigee r_min = 2.541e-4 AU (38,012 km ~ 5.96 R_Earth) with v_inf = 5.86 km/s
      (0.1967 v_Earth), deflecting Apophis outward from Aten class (a_pre = 0.9224 AU < 1 AU)
      to Apollo class (a_post ~ 1.116 AU > 1 AU, Delta_a > 0).
    """
    m_sun = 1.0
    m_earth = 3.0035e-6
    m_moon = 3.6943e-8
    m_apo = 4.0e-19
    masses = np.array([m_sun, m_earth, m_moon, m_apo], dtype=np.float64)

    a_earth = 1.00000
    a_moon = 0.0025695
    v_earth = math.sqrt(float(G) * (m_sun + m_earth + m_moon) / a_earth)
    v_moon_rel = math.sqrt(float(G) * (m_earth + m_moon) / a_moon)

    r0_e = 0.0020
    rp_e = 2.541e-4  # 38,012 km in AU
    v_inf_e = 0.1967 * v_earth
    mu_e = float(G) * m_earth
    v0_e = math.sqrt(v_inf_e * v_inf_e + 2.0 * mu_e / r0_e)
    vp_e = math.sqrt(v_inf_e * v_inf_e + 2.0 * mu_e / rp_e)
    sin_gamma = min(0.99, (rp_e * vp_e) / (r0_e * v0_e))
    cos_gamma = math.sqrt(1.0 - sin_gamma * sin_gamma)

    # Trailing-side approach rotation angle phi_approach = 171.575217 deg so a_pre = 0.9224 AU -> a_post > 1.10 AU
    phi_app = math.radians(171.575217)
    cos_phi = math.cos(phi_app)
    sin_phi = math.sin(phi_app)

    dr_x = r0_e * cos_phi
    dr_y = r0_e * sin_phi
    dv0_x = -v0_e * cos_gamma
    dv0_y = v0_e * sin_gamma
    dv_x = dv0_x * cos_phi - dv0_y * sin_phi
    dv_y = dv0_x * sin_phi + dv0_y * cos_phi

    q = np.array([
        [0.0, 0.0, 0.0],
        [a_earth, 0.0, 0.0],
        [a_earth, a_moon, 0.0],
        [a_earth + dr_x, dr_y, 1.0e-5],
    ], dtype=np.float64)

    v = np.array([
        [0.0, 0.0, 0.0],
        [0.0, v_earth, 0.0],
        [-v_moon_rel, v_earth, 0.0],
        [dv_x, v_earth + dv_y, 0.0],
    ], dtype=np.float64)

    cdm = np.sum(masses[:, None] * q, axis=0) / np.sum(masses)
    q -= cdm
    p = v * masses[:, np.newaxis]
    p_cm = np.sum(p, axis=0) / np.sum(masses)
    p -= masses[:, np.newaxis] * p_cm

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_hierarchical_3d_eccentric_triple_benchmark(
    e_in: float = 0.30,
    inc_deg: float = 28.0,
    a_in: float = 1.0,
    a_out: float = 10.0,
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Canonical 3D Eccentric Hierarchical Triple System (Chapter 5, Theorem 5.14 &
    Chapter 7, Section 7.1.2 Universal Obstruction Benchmark):
    - Inner binary: m0 = 0.5, m1 = 0.5 with semi-major axis a_in = 1.0,
      eccentricity e_in = 0.30, initialized away from apsides so <r_in, v_in> != 0.
    - Outer tertiary: m2 = 0.5 orbiting at a_out = 10.0 with mutual inclination
      inc_deg = 28.0 deg.
    """
    m0, m1, m2 = 0.5, 0.5, 0.5
    masses = np.array([m0, m1, m2], dtype=np.float64)
    m_in = m0 + m1
    m_tot = m_in + m2

    # Inner binary at eccentric anomaly E_0 = pi/4 (non-apsidal phase with radial velocity != 0)
    E_anom = 0.25 * math.pi
    cos_E = math.cos(E_anom)
    sin_E = math.sin(E_anom)
    r_in_norm = float(a_in) * (1.0 - float(e_in) * cos_E)
    x_in = float(a_in) * (cos_E - float(e_in))
    y_in = float(a_in) * math.sqrt(1.0 - float(e_in) * float(e_in)) * sin_E
    n_in = math.sqrt(float(G) * m_in / (float(a_in) ** 3))
    vx_in = -n_in * (float(a_in) ** 2 / r_in_norm) * sin_E
    vy_in = n_in * (float(a_in) ** 2 / r_in_norm) * math.sqrt(1.0 - float(e_in) * float(e_in)) * cos_E

    # Outer tertiary inclined by inc_deg at true anomaly f_out = pi/3
    inc_rad = math.radians(float(inc_deg))
    f_out = math.pi / 3.0
    v_out = math.sqrt(float(G) * m_tot / float(a_out))
    r_out_vec = np.array([
        float(a_out) * math.cos(f_out),
        float(a_out) * math.sin(f_out) * math.cos(inc_rad),
        float(a_out) * math.sin(f_out) * math.sin(inc_rad),
    ], dtype=np.float64)
    v_out_vec = np.array([
        -v_out * math.sin(f_out),
        v_out * math.cos(f_out) * math.cos(inc_rad),
        v_out * math.cos(f_out) * math.sin(inc_rad),
    ], dtype=np.float64)

    r_in_vec = np.array([x_in, y_in, 0.0], dtype=np.float64)
    v_in_vec = np.array([vx_in, vy_in, 0.0], dtype=np.float64)

    q_cm_in = -(m2 / m_tot) * r_out_vec
    v_cm_in = -(m2 / m_tot) * v_out_vec
    q2 = (m_in / m_tot) * r_out_vec
    v2 = (m_in / m_tot) * v_out_vec

    q0 = q_cm_in - (m1 / m_in) * r_in_vec
    q1 = q_cm_in + (m0 / m_in) * r_in_vec
    v0 = v_cm_in - (m1 / m_in) * v_in_vec
    v1 = v_cm_in + (m0 / m_in) * v_in_vec

    q = np.vstack([q0, q1, q2])
    v = np.vstack([v0, v1, v2])
    p = v * masses[:, np.newaxis]

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


def create_inner_solar_system_mercury_benchmark(
    G: float = 1.0,
) -> ExtendedPhaseState:
    """
    Jacques Laskar's Inner Planetary Hierarchy Benchmark (Sun + Mercury, Venus,
    Earth, Mars; N = 5, Farrés, Laskar et al. 2013, Table 1, eps_Jacobi = 4.54e-6):
    - Dominated by Mercury's short orbital period (a_1 = 0.387098 AU) and high
      eccentricity (e_1 = 0.205630) plus inclination (i_1 = 7.005 deg), which
      creates the classic O(eps * tau^{2p}) pericenter resolution bottleneck.
    - All 4 terrestrial planets initialized at perihelion r_p = a_i * (1 - e_i)
      with exact inclination and perihelion velocity v_p = sqrt(G*(M_sun+m_i)*(1+e_i)/(a_i*(1-e_i))).
    """
    m_sun = 1.0
    m_merc = 1.6601e-7
    m_ven = 2.4478e-6
    m_earth = 3.0404e-6
    m_mars = 3.2271e-7
    masses = np.array([m_sun, m_merc, m_ven, m_earth, m_mars], dtype=np.float64)

    a_list = [0.387098, 0.723332, 1.000000, 1.523679]
    e_list = [0.205630, 0.006772, 0.016708, 0.093400]
    inc_deg_list = [7.005, 3.395, 0.000, 1.850]

    q = np.zeros((5, 3), dtype=np.float64)
    v = np.zeros((5, 3), dtype=np.float64)

    for idx in range(4):
        m_p = masses[idx + 1]
        a_p = a_list[idx]
        e_p = e_list[idx]
        inc_rad = math.radians(inc_deg_list[idx])
        r_peri = a_p * (1.0 - e_p)
        v_peri = math.sqrt(float(G) * (m_sun + m_p) * (1.0 + e_p) / r_peri)
        # Alternate initial longitudes (0, pi/2, pi, 3pi/2) to avoid artificial alignment
        theta = 0.5 * math.pi * idx
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)
        cos_i = math.cos(inc_rad)
        sin_i = math.sin(inc_rad)

        q[idx + 1, 0] = r_peri * cos_t
        q[idx + 1, 1] = r_peri * sin_t * cos_i
        q[idx + 1, 2] = r_peri * sin_t * sin_i

        v[idx + 1, 0] = -v_peri * sin_t
        v[idx + 1, 1] = v_peri * cos_t * cos_i
        v[idx + 1, 2] = v_peri * cos_t * sin_i

    cdm = np.sum(masses[:, None] * q, axis=0) / np.sum(masses)
    q -= cdm
    p = v * masses[:, np.newaxis]
    p_cm = np.sum(p, axis=0) / np.sum(masses)
    p -= masses[:, np.newaxis] * p_cm

    E0 = compute_cartesian_hamiltonian(q, p, masses, G=G)
    return ExtendedPhaseState(q=q, p=p, masses=masses, t=0.0, E=E0, E0=E0)


BENCHMARK_REGISTRY: Dict[str, Callable[..., ExtendedPhaseState]] = {
    "radial_collision": create_radial_collision_benchmark,
    "burrau_pythagorean": create_burrau_pythagorean_benchmark,
    "arenstorf_lunar": create_arenstorf_lunar_benchmark,
    "figure8_choreography": create_figure8_choreography_benchmark,
    "sitnikov_oscillatory": create_sitnikov_oscillatory_benchmark,
    "quadruple_collision": create_quadruple_collision_benchmark,
    "jovian_solar_system": create_jovian_solar_system_benchmark,
    "pareto_test": create_pareto_test_benchmark,
    "eccentric_binary": create_eccentric_binary_benchmark,
    "jovian_laplace_resonance": create_jovian_laplace_resonance_benchmark,
    "comet_lexell_1770": create_comet_lexell_1770_benchmark,
    "shoemaker_levy_9_1992": create_shoemaker_levy_9_1992_benchmark,
    "apophis_2029_flyby": create_apophis_2029_flyby_benchmark,
    "hierarchical_3d_eccentric_triple": create_hierarchical_3d_eccentric_triple_benchmark,
    "inner_solar_system_mercury": create_inner_solar_system_mercury_benchmark,
}

BENCHMARK_PROFILES: Dict[str, str] = {
    "radial_collision": "chaotic_encounter",
    "burrau_pythagorean": "chaotic_encounter",
    "arenstorf_lunar": "smooth_periodic",
    "figure8_choreography": "smooth_periodic",
    "sitnikov_oscillatory": "chaotic_encounter",
    "quadruple_collision": "chaotic_encounter",
    "jovian_solar_system": "hierarchical_multiscale",
    "pareto_test": "smooth_periodic",
    "hierarchical_4body": "hierarchical_multiscale",
    "eccentric_binary": "smooth_periodic",
    "jovian_laplace_resonance": "hierarchical_multiscale",
    "comet_lexell_1770": "hierarchical_multiscale",
    "shoemaker_levy_9_1992": "hierarchical_multiscale",
    "apophis_2029_flyby": "hierarchical_multiscale",
    "hierarchical_3d_eccentric_triple": "hierarchical_multiscale",
    "inner_solar_system_mercury": "hierarchical_multiscale",
}


def get_benchmark(name: str, **kwargs: Any) -> ExtendedPhaseState:
    """Retrieves an instantiated benchmark ExtendedPhaseState by name."""
    if name not in BENCHMARK_REGISTRY:
        raise KeyError(f"Unknown benchmark '{name}'. Available: {list(BENCHMARK_REGISTRY.keys())}")
    return BENCHMARK_REGISTRY[name](**kwargs)


def get_benchmark_profile(name: str) -> str:
    """Retrieves the canonical KPI profile for a benchmark configuration."""
    return BENCHMARK_PROFILES.get(name, "smooth_periodic")
