-- Dependencies: Mathlib.Tactic.Ring
import Mathlib.Tactic.Ring

namespace UnsoftenedNBody.PCR3BP

/-- Planar 2D vector in the rotating synodic frame over a commutative ring k. -/
structure Vec2 (k : Type) where
  x : k
  y : k

namespace PhaseState
variable {k : Type} [CommRing k]

/-- Fundamental Theorem: Exact Jacobi first integral conservation identity dC_J/dt ≡ 0 under the synodic flow. -/
theorem jacobi_conservation_identity
    (q p dq dp gradV : Vec2 k)
    (h_dq : dq = ⟨p.x + q.y, p.y - q.x⟩)
    (h_dp : dp = ⟨p.y - gradV.x, -p.x - gradV.y⟩)
    (dV : k)
    (h_dV : dV = gradV.x * dq.x + gradV.y * dq.y) :
    let d_ang := dq.x * p.y - dq.y * p.x + q.x * dp.y - q.y * dp.x
    let d_kin := 2 * (p.x * dp.x + p.y * dp.y)
    let d_CJ := 2 * d_ang - d_kin - 2 * dV
    d_CJ = 0 := by
  intro d_ang d_kin d_CJ
  dsimp [d_CJ, d_ang, d_kin]
  rw [h_dV, h_dq, h_dp]
  ring

/-- Canonical Darboux 2-form on planar phase space T*R² evaluated on tangent pairs (q1, p1) and (q2, p2). -/
def omegaPlanar (q1 p1 q2 p2 : Vec2 k) : k :=
  (q1.x * p2.x - p1.x * q2.x) + (q1.y * p2.y - p1.y * q2.y)

/-- Exact finite-step linear Coriolis-kinematic propagator Ψ_h^kin(q, p) = (R(-h)(q + h p), R(-h) p). -/
def coriolisStep (c s h : k) (q p : Vec2 k) : Vec2 k × Vec2 k :=
  (⟨c * (q.x + h * p.x) + s * (q.y + h * p.y),
    -s * (q.x + h * p.x) + c * (q.y + h * p.y)⟩,
   ⟨c * p.x + s * p.y, -s * p.x + c * p.y⟩)

/-- Fundamental Symplecticity Theorem: The discrete Coriolis-kinematic step Ψ_h^kin belongs to Sp(4, k),
    preserving the planar Darboux 2-form ω₀ identically for any rotation (c, s) on the unit circle c² + s² = 1. -/
theorem coriolis_step_symplectic (c s h : k) (h_rot : c^2 + s^2 = 1)
    (q1 p1 q2 p2 : Vec2 k) :
    let s1 := coriolisStep c s h q1 p1
    let s2 := coriolisStep c s h q2 p2
    omegaPlanar s1.1 s1.2 s2.1 s2.2 = omegaPlanar q1 p1 q2 p2 := by
  intro s1 s2
  dsimp [s1, s2, coriolisStep, omegaPlanar]
  calc _ = (c^2 + s^2) * ((q1.x * p2.x - p1.x * q2.x) + (q1.y * p2.y - p1.y * q2.y)) := by ring
    _ = 1 * ((q1.x * p2.x - p1.x * q2.x) + (q1.y * p2.y - p1.y * q2.y)) := by rw [h_rot]
    _ = (q1.x * p2.x - p1.x * q2.x) + (q1.y * p2.y - p1.y * q2.y) := by ring

end PhaseState
end UnsoftenedNBody.PCR3BP
