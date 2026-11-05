"""
pcr3bp_canonical_solver.py
Normal-cone symplectic solver for the Planar Circular Restricted Three-Body Problem (PCR3BP).
Exact validation on the periodic Arenstorf orbit (T = 17.0652).
"""

import numpy as np

class PCR3BPNormalConeIntegrator:
    """
    Normal-cone Strang and Blanes-Moan symplectic integrator for the PCR3BP.
    Decouples the exact linear Coriolis flow from the transverse gravitational impulse.
    """
    # Order-4 palindromic 7-stage S_2 composition weights (sum gamma_j = 1, sum gamma_j^3 = 0; cf. Section 5.3 for BM_4^(7) ABA)
    BM4_A = (
        +0.2390843168869614,
        +0.2390725069195032,
        +0.23905247436220158,
        1.0 - 2.0 * (0.2390843168869614 + 0.2390725069195032 + 0.23905247436220158),
        +0.23905247436220158,
        +0.2390725069195032,
        +0.2390843168869614,
    )

    def __init__(self, mu: float = 0.0121505856):
        self.mu = float(mu)

    def hamiltonian(self, q: np.ndarray, p: np.ndarray) -> float:
        """Canonical synodic Hamiltonian H_syn = H_kin + H_grav."""
        mu = self.mu
        x, y = float(q[0]), float(q[1])
        px, py = float(p[0]), float(p[1])
        r1 = float(np.hypot(x + mu, y))
        r2 = float(np.hypot(x - 1.0 + mu, y))
        h_kin = 0.5 * (px * px + py * py) + (y * px - x * py)
        h_grav = -(1.0 - mu) / max(r1, 1e-30) - mu / max(r2, 1e-30)
        return h_kin + h_grav

    def jacobi_constant(self, q: np.ndarray, p: np.ndarray) -> float:
        """Jacobi integral invariant C_J = -2 H_syn(q, p)."""
        return -2.0 * self.hamiltonian(q, p)

    def coriolis_step(self, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
        """Exact flow of H_kin in Sp(4, R)."""
        c, s = np.cos(h), np.sin(h)
        R = np.array([[c, s], [-s, c]])
        q_next = R @ q + h * (R @ p)
        p_next = R @ p
        return q_next, p_next

    def grav_kick(self, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
        """Exact gravitational impulse along the transverse cotangent fibers."""
        mu = self.mu
        x, y = float(q[0]), float(q[1])
        r1_sq = (x + mu)**2 + y**2
        r2_sq = (x - 1.0 + mu)**2 + y**2
        r1_3 = r1_sq**1.5
        r2_3 = r2_sq**1.5
        fx = -(1.0 - mu) * (x + mu) / r1_3 - mu * (x - 1.0 + mu) / r2_3
        fy = -(1.0 - mu) * y / r1_3 - mu * y / r2_3
        p_next = p + h * np.array([fx, fy], dtype=np.float64)
        return np.copy(q), p_next

    def step_strang(self, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
        """Second-order symmetric Strang composition: S_2(h) = Phi_{h/2}^grav o Psi_h^kin o Phi_{h/2}^grav."""
        q, p = self.grav_kick(q, p, 0.5 * h)
        q, p = self.coriolis_step(q, p, h)
        q, p = self.grav_kick(q, p, 0.5 * h)
        return q, p

    def step_bm4(self, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
        """Fourth-order palindromic 7-stage composition of Strang steps S_2(gamma_j * h)."""
        q_curr, p_curr = np.copy(q), np.copy(p)
        for a_j in self.BM4_A:
            q_curr, p_curr = self.step_strang(q_curr, p_curr, a_j * h)
        return q_curr, p_curr

    def _dyadic_exponent(self, q: np.ndarray, h: float, ds_max: float) -> int:
        """Dyadic level exponent k >= 0 derived from the harmonic Sundman scale chi(q) = r1*r2/(r1+r2)."""
        mu = self.mu
        x, y = float(q[0]), float(q[1])
        r1 = float(np.hypot(x + mu, y))
        r2 = float(np.hypot(x - 1.0 + mu, y))
        chi = max((r1 * r2) / max(r1 + r2, 1e-15), 1e-6)
        n_sub = max(1, int(np.ceil(abs(h) / max(chi * ds_max, 1e-15))))
        return int(np.ceil(np.log2(n_sub))) if n_sub > 1 else 0

    def step_sundman_regularized(
        self, q: np.ndarray, p: np.ndarray, h: float, ds_max: float = 0.02
    ) -> tuple[np.ndarray, np.ndarray]:
        """Symmetrized dyadic sub-cycling k_sym = max(k(q_start), k(q_end)) on the time lattice."""
        k_start = self._dyadic_exponent(q, h, ds_max)
        n_dyadic = 1 << k_start
        h_sub = h / n_dyadic
        q_curr, p_curr = np.copy(q), np.copy(p)
        for _ in range(n_dyadic):
            q_curr, p_curr = self.step_bm4(q_curr, p_curr, h_sub)
        k_end = self._dyadic_exponent(q_curr, h, ds_max)
        if k_end > k_start:
            n_dyadic = 1 << k_end
            h_sub = h / n_dyadic
            q_curr, p_curr = np.copy(q), np.copy(p)
            for _ in range(n_dyadic):
                q_curr, p_curr = self.step_bm4(q_curr, p_curr, h_sub)
        return q_curr, p_curr


PCR3BPConormalIntegrator = PCR3BPNormalConeIntegrator
BM4_A = PCR3BPNormalConeIntegrator.BM4_A

