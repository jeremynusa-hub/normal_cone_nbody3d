-- Dependencies: 01_ConfigurationSpace_Defs, 05_HeavyTop_Defs
import Mathlib.Tactic.Ring

universe u

namespace HeavyTop

/-- Theorem: The transverse normal fiber torque kick strictly preserves Casimir C₁ = ||Γ||². -/
theorem torque_kick_preserves_casimir1 {k : Type u} [CommRing k]
    (c : k) (chi : Vec3 k) (h : k) (s : HeavyTopState k) :
    casimir1 (conormalTorqueKick c chi h s) = casimir1 s := by
  dsimp [conormalTorqueKick, casimir1]

/-- Theorem: Orthogonality of cross product: (Γ × χ) · Γ = 0. -/
theorem cross_dot_self {k : Type u} [CommRing k] (a b : Vec3 k) :
    Vec3.dot (Vec3.cross a b) a = 0 := by
  dsimp [Vec3.dot, Vec3.cross]
  ring

/-- Theorem: The transverse normal fiber torque kick strictly preserves Casimir C₂ = Π · Γ.
    This holds algebraically because the gravitational torque Γ × χ is strictly
    orthogonal to the direction of gravity Γ. -/
theorem torque_kick_preserves_casimir2 {k : Type u} [CommRing k]
    (c : k) (chi : Vec3 k) (h : k) (s : HeavyTopState k) :
    casimir2 (conormalTorqueKick c chi h s) = casimir2 s := by
  dsimp [conormalTorqueKick, casimir2, Vec3.dot, Vec3.add, Vec3.smul, Vec3.cross]
  ring

/-- Theorem: The groupoid rotation step strictly preserves Casimir C₁ = ||Γ||². -/
theorem groupoid_rot_preserves_casimir1 {k : Type u} [CommRing k]
    (R : SO3Isometry k) (s : HeavyTopState k) :
    casimir1 (groupoidRotationStep R s) = casimir1 s := by
  dsimp [groupoidRotationStep, casimir1]
  exact R.preserves_dot s.Gamma s.Gamma

/-- Theorem: The groupoid rotation step strictly preserves Casimir C₂ = Π · Γ. -/
theorem groupoid_rot_preserves_casimir2 {k : Type u} [CommRing k]
    (R : SO3Isometry k) (s : HeavyTopState k) :
    casimir2 (groupoidRotationStep R s) = casimir2 s := by
  dsimp [groupoidRotationStep, casimir2]
  exact R.preserves_dot s.Pi s.Gamma

/-- Fundamental Theorem: Exact Casimir Preservation of the McLachlan Normal Cone Top Splitting.
    Every step of the Strang splitting (kick and rotation) preserves both Casimirs C₁ and C₂
    identically and unconditionally. -/
theorem heavy_top_strang_casimirs_exact {k : Type u} [CommRing k]
    (c : k) (chi : Vec3 k) (R : SO3Isometry k) (h_half : k) (s : HeavyTopState k) :
    let s1 := conormalTorqueKick c chi h_half s
    let s2 := groupoidRotationStep R s1
    let s3 := conormalTorqueKick c chi h_half s2
    casimir1 s3 = casimir1 s ∧ casimir2 s3 = casimir2 s := by
  intro s1 s2 s3
  dsimp [s3, s2, s1]
  rw [torque_kick_preserves_casimir1]
  rw [groupoid_rot_preserves_casimir1]
  rw [torque_kick_preserves_casimir1]
  rw [torque_kick_preserves_casimir2]
  rw [groupoid_rot_preserves_casimir2]
  rw [torque_kick_preserves_casimir2]
  exact ⟨rfl, rfl⟩

/-- Concrete Instantiation: Exact Casimir Preservation under Explicit 3D Rodrigues Rotation.
    For any unit rotation axis `n` (`n · n = 1`) and trigonometric angle parameters `(c_ang, s_ang)`
    satisfying `c_ang^2 + s_ang^2 = 1`, the Strang splitting with explicit 3D Rodrigues rotation
    `rodriguesIsometry` preserves both Casimirs `C₁` and `C₂` identically. -/
theorem heavy_top_strang_rodrigues_casimirs_exact {k : Type u} [CommRing k]
    (c : k) (chi n : Vec3 k) (c_ang s_ang : k)
    (hn : Vec3.dot n n = 1) (h_circle : c_ang * c_ang + s_ang * s_ang = 1)
    (h_half : k) (s : HeavyTopState k) :
    let R := rodriguesIsometry n c_ang s_ang hn h_circle
    let s1 := conormalTorqueKick c chi h_half s
    let s2 := groupoidRotationStep R s1
    let s3 := conormalTorqueKick c chi h_half s2
    casimir1 s3 = casimir1 s ∧ casimir2 s3 = casimir2 s :=
  heavy_top_strang_casimirs_exact c chi (rodriguesIsometry n c_ang s_ang hn h_circle) h_half s

end HeavyTop
