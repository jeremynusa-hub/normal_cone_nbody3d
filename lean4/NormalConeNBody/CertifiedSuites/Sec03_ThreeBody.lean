import Mathlib.Algebra.Ring.Defs
import Mathlib.Tactic.Ring

-- === 01_PCR3BPHamiltonian_Theorems.lean ===
namespace UnsoftenedNBody.PCR3BP

structure Vec2 (k : Type) where
  x : k
  y : k

namespace PhaseState
variable {k : Type} [CommRing k]

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

def omegaPlanar (q1 p1 q2 p2 : Vec2 k) : k :=
  (q1.x * p2.x - p1.x * q2.x) + (q1.y * p2.y - p1.y * q2.y)

def coriolisStep (c s h : k) (q p : Vec2 k) : Vec2 k × Vec2 k :=
  (⟨c * (q.x + h * p.x) + s * (q.y + h * p.y),
    -s * (q.x + h * p.x) + c * (q.y + h * p.y)⟩,
   ⟨c * p.x + s * p.y, -s * p.x + c * p.y⟩)

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

-- === 02_PCR3BPRegularization_Theorems.lean ===
namespace LeviCivita

universe u

variable {k : Type u} [CommRing k]

def levi_civita_map (u v : k) : k × k :=
  (u^2 - v^2, 2 * u * v)

theorem levi_civita_normal_cone_covering (u v : k) :
    let (x, y) := levi_civita_map u v
    x^2 + y^2 = (u^2 + v^2)^2 := by
  dsimp [levi_civita_map]
  ring

theorem levi_civita_antipodal_covering (u v : k) :
    levi_civita_map (-u) (-v) = levi_civita_map u v := by
  dsimp [levi_civita_map]
  ext <;> ring

theorem levi_civita_conformal_orthogonality (u v : k) :
    (2 * u) * (-2 * v) + (2 * v) * (2 * u) = 0 := by
  ring

theorem levi_civita_conformal_scaling (u v : k) :
    (2 * u) * (2 * u) - (-2 * v) * (2 * v) = 4 * (u^2 + v^2) := by
  ring

end LeviCivita

#print axioms UnsoftenedNBody.PCR3BP.PhaseState.jacobi_conservation_identity
#print axioms UnsoftenedNBody.PCR3BP.PhaseState.coriolis_step_symplectic
#print axioms LeviCivita.levi_civita_normal_cone_covering
#print axioms LeviCivita.levi_civita_antipodal_covering
#print axioms LeviCivita.levi_civita_conformal_orthogonality
#print axioms LeviCivita.levi_civita_conformal_scaling
