-- Dependencies: 01_ConfigurationSpace_Defs
import Mathlib.Tactic.Ring

universe u

namespace HeavyTop

/-- State of the heavy top on the dual Lie algebra se(3)*.
    `Pi` is the angular momentum in the body frame.
    `Gamma` is the unit gravity vector in the body frame. -/
structure HeavyTopState (k : Type u) where
  Pi : Vec3 k
  Gamma : Vec3 k

/-- First Casimir invariant: C₁ = ||Γ||² = Γ · Γ. -/
def casimir1 {k : Type u} [CommRing k] (s : HeavyTopState k) : k :=
  Vec3.dot s.Gamma s.Gamma

/-- Second Casimir invariant: C₂ = Π · Γ. -/
def casimir2 {k : Type u} [CommRing k] (s : HeavyTopState k) : k :=
  Vec3.dot s.Pi s.Gamma

/-- Normal fiber torque kick in the transverse fiber:
    Under gravity, the torque is K(Γ) = c * (Γ × χ), where χ is the vector
    from the fixed point to the center of mass.
    The kick updates Π ↦ Π + h * K(Γ) while leaving Γ invariant. -/
def conormalTorqueKick {k : Type u} [CommRing k]
    (c : k) (chi : Vec3 k) (h : k) (s : HeavyTopState k) : HeavyTopState k :=
  let torque := Vec3.smul c (Vec3.cross s.Gamma chi)
  ⟨Vec3.add s.Pi (Vec3.smul h torque), s.Gamma⟩

/-- Orthogonal transformation on Vec3 (modeling an element of SO(3)). -/
structure SO3Isometry (k : Type u) [CommRing k] where
  rot : Vec3 k → Vec3 k
  preserves_dot : ∀ (u v : Vec3 k), Vec3.dot (rot u) (rot v) = Vec3.dot u v

/-- Explicit principal-axis rotation around the z-axis with cosine `c` and sine `s`. -/
def zAxisRot {k : Type u} [CommRing k] (c s : k) (v : Vec3 k) : Vec3 k :=
  ⟨c * v.x - s * v.y, s * v.x + c * v.y, v.z⟩

/-- Constructive SO(3) isometry witness for principal-axis rotation around the z-axis
    whenever `c^2 + s^2 = 1`. -/
def zAxisRotationIsometry {k : Type u} [CommRing k]
    (c s : k) (h_circle : c * c + s * s = 1) : SO3Isometry k where
  rot := zAxisRot c s
  preserves_dot := by
    intro u v
    calc Vec3.dot (zAxisRot c s u) (zAxisRot c s v)
      _ = (c * c + s * s) * (u.x * v.x + u.y * v.y) + u.z * v.z := by
          dsimp [zAxisRot, Vec3.dot]; ring
      _ = 1 * (u.x * v.x + u.y * v.y) + u.z * v.z := by rw [h_circle]
      _ = Vec3.dot u v := by dsimp [Vec3.dot]; ring

/-- Explicit 3D Rodrigues rotation around a unit axis `n` (`n · n = 1`) with cosine `c` and sine `s`:
    rot(v) = c v + s (n × v) + (1 - c) (n · v) n. -/
def rodriguesRot {k : Type u} [CommRing k] (n : Vec3 k) (c s : k) (v : Vec3 k) : Vec3 k :=
  Vec3.add (Vec3.add (Vec3.smul c v) (Vec3.smul s (Vec3.cross n v)))
           (Vec3.smul ((1 - c) * Vec3.dot n v) n)

/-- Constructive 3D Rodrigues SO(3) isometry witness: for any unit axis `n` (`n · n = 1`)
    and trigonometric pair `(c, s)` satisfying `c^2 + s^2 = 1`, Rodrigues' formula defines
    an exact Euclidean isometry on `Vec3 k`. -/
def rodriguesIsometry {k : Type u} [CommRing k]
    (n : Vec3 k) (c s : k)
    (hn : Vec3.dot n n = 1) (h_circle : c * c + s * s = 1) : SO3Isometry k where
  rot := rodriguesRot n c s
  preserves_dot := by
    intro u v
    calc Vec3.dot (rodriguesRot n c s u) (rodriguesRot n c s v)
      _ = (c * c + s * s * Vec3.dot n n) * Vec3.dot u v +
          ((1 - c) * (1 - c) * Vec3.dot n n + 2 * c * (1 - c) - s * s) *
            (Vec3.dot n u * Vec3.dot n v) := by
          dsimp [rodriguesRot, Vec3.dot, Vec3.add, Vec3.smul, Vec3.cross]; ring
      _ = (c * c + s * s * 1) * Vec3.dot u v +
          ((1 - c) * (1 - c) * 1 + 2 * c * (1 - c) - s * s) *
            (Vec3.dot n u * Vec3.dot n v) := by rw [hn]
      _ = (c * c + s * s) * Vec3.dot u v +
          (1 - (c * c + s * s)) * (Vec3.dot n u * Vec3.dot n v) := by ring
      _ = 1 * Vec3.dot u v + (1 - 1) * (Vec3.dot n u * Vec3.dot n v) := by rw [h_circle]
      _ = Vec3.dot u v := by ring

/-- Groupoid kinematic rotation step on se(3)*:
    Both Π and Γ rotate by the same body-frame rotation R ∈ SO(3). -/
def groupoidRotationStep {k : Type u} [CommRing k]
    (R : SO3Isometry k) (s : HeavyTopState k) : HeavyTopState k :=
  ⟨R.rot s.Pi, R.rot s.Gamma⟩

end HeavyTop
