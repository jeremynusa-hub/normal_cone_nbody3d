import Mathlib.Tactic.Ring
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

universe u

structure Vec3 (k : Type u) where
  x : k
  y : k
  z : k
  deriving Repr, DecidableEq

namespace Vec3
variable {k : Type u}

def zero [Zero k] : Vec3 k := ⟨0, 0, 0⟩
def add [Add k] (u v : Vec3 k) : Vec3 k := ⟨u.x + v.x, u.y + v.y, u.z + v.z⟩
def sub [Sub k] (u v : Vec3 k) : Vec3 k := ⟨u.x - v.x, u.y - v.y, u.z - v.z⟩
def neg [Neg k] (u : Vec3 k) : Vec3 k := ⟨-u.x, -u.y, -u.z⟩
def smul [Mul k] (s : k) (v : Vec3 k) : Vec3 k := ⟨s * v.x, s * v.y, s * v.z⟩
def dot [Add k] [Mul k] (u v : Vec3 k) : k := u.x * v.x + u.y * v.y + u.z * v.z
def cross [Sub k] [Mul k] (u v : Vec3 k) : Vec3 k :=
  ⟨u.y * v.z - u.z * v.y,
   u.z * v.x - u.x * v.z,
   u.x * v.y - u.y * v.x⟩

@[ext]
theorem ext {u v : Vec3 k} (hx : u.x = v.x) (hy : u.y = v.y) (hz : u.z = v.z) : u = v := by
  cases u; cases v; congr

end Vec3

-- === 01_ReesDeformation_Defs.lean ===
structure NormalConePoint (k : Type u) where
  q_parallel : Vec3 k
  n          : Vec3 k

def deformationPhysicalPosition {k : Type u} [CommRing k]
    (d : NormalConePoint k) (t : k) : Vec3 k :=
  Vec3.add d.q_parallel (Vec3.smul t d.n)

structure ConormalTangentVector (k : Type u) where
  dq : Vec3 k
  dp : Vec3 k
  dn : Vec3 k
  dw : Vec3 k

def omegaNormal {k : Type u} [CommRing k]
    (v1 v2 : ConormalTangentVector k) : k :=
  (Vec3.dot v1.dq v2.dp - Vec3.dot v2.dq v1.dp) +
  (Vec3.dot v1.dn v2.dw - Vec3.dot v2.dn v1.dw)

def transverseConormalTangentFlow {k : Type u} [CommRing k]
    (s : k) (dF : Vec3 k) (v : ConormalTangentVector k) : ConormalTangentVector k :=
  { dq := v.dq
    dp := v.dp
    dn := v.dn
    dw := Vec3.add v.dw (Vec3.smul s dF) }

-- === 02_ReesDeformation_Theorems.lean ===
section ReesTheorems
variable {k : Type u} [CommRing k]

theorem deformation_limit_zero (d : NormalConePoint k) :
    deformationPhysicalPosition d 0 = d.q_parallel := by
  dsimp [deformationPhysicalPosition, Vec3.smul, Vec3.add]
  cases d with
  | mk qp n =>
    dsimp
    ext <;> ring

theorem omegaNormal_skew (v1 v2 : ConormalTangentVector k) :
    omegaNormal v1 v2 = - omegaNormal v2 v1 := by
  dsimp [omegaNormal]
  ring

theorem omegaNormal_self (v : ConormalTangentVector k) :
    omegaNormal v v = 0 := by
  dsimp [omegaNormal]
  ring

theorem transverse_conormal_flow_symplectic
    (s : k) (dF1 dF2 : Vec3 k) (v1 v2 : ConormalTangentVector k)
    (h_hess_symm : Vec3.dot dF1 v2.dn = Vec3.dot dF2 v1.dn) :
    omegaNormal (transverseConormalTangentFlow s dF1 v1)
                (transverseConormalTangentFlow s dF2 v2) =
    omegaNormal v1 v2 := by
  dsimp [transverseConormalTangentFlow, omegaNormal, Vec3.dot, Vec3.add, Vec3.smul]
  calc (v1.dq.x * v2.dp.x + v1.dq.y * v2.dp.y + v1.dq.z * v2.dp.z -
        (v2.dq.x * v1.dp.x + v2.dq.y * v1.dp.y + v2.dq.z * v1.dp.z)) +
       (v1.dn.x * (v2.dw.x + s * dF2.x) + v1.dn.y * (v2.dw.y + s * dF2.y) + v1.dn.z * (v2.dw.z + s * dF2.z) -
        (v2.dn.x * (v1.dw.x + s * dF1.x) + v2.dn.y * (v1.dw.y + s * dF1.y) + v2.dn.z * (v1.dw.z + s * dF1.z)))
    = omegaNormal v1 v2 + s * (Vec3.dot v1.dn dF2 - Vec3.dot v2.dn dF1) := by
        dsimp [omegaNormal, Vec3.dot]; ring
    _ = omegaNormal v1 v2 + s * (Vec3.dot dF1 v2.dn - Vec3.dot v2.dn dF1) := by
        have h_comm1 : Vec3.dot v1.dn dF2 = Vec3.dot dF2 v1.dn := by
          dsimp [Vec3.dot]; ring
        have h_symm : Vec3.dot dF2 v1.dn = Vec3.dot dF1 v2.dn := h_hess_symm.symm
        rw [h_comm1, h_symm]
    _ = omegaNormal v1 v2 := by
        have h_comm2 : Vec3.dot dF1 v2.dn = Vec3.dot v2.dn dF1 := by
          dsimp [Vec3.dot]; ring
        rw [h_comm2]
        ring

theorem omegaNormal_nondegenerate
    (v1 : ConormalTangentVector k)
    (h : ∀ v2 : ConormalTangentVector k, omegaNormal v1 v2 = 0) :
    v1.dq = Vec3.zero ∧ v1.dp = Vec3.zero ∧ v1.dn = Vec3.zero ∧ v1.dw = Vec3.zero := by
  have hq_x := h ⟨Vec3.zero, ⟨1, 0, 0⟩, Vec3.zero, Vec3.zero⟩
  have hq_y := h ⟨Vec3.zero, ⟨0, 1, 0⟩, Vec3.zero, Vec3.zero⟩
  have hq_z := h ⟨Vec3.zero, ⟨0, 0, 1⟩, Vec3.zero, Vec3.zero⟩
  have hp_x := h ⟨⟨1, 0, 0⟩, Vec3.zero, Vec3.zero, Vec3.zero⟩
  have hp_y := h ⟨⟨0, 1, 0⟩, Vec3.zero, Vec3.zero, Vec3.zero⟩
  have hp_z := h ⟨⟨0, 0, 1⟩, Vec3.zero, Vec3.zero, Vec3.zero⟩
  have hn_x := h ⟨Vec3.zero, Vec3.zero, Vec3.zero, ⟨1, 0, 0⟩⟩
  have hn_y := h ⟨Vec3.zero, Vec3.zero, Vec3.zero, ⟨0, 1, 0⟩⟩
  have hn_z := h ⟨Vec3.zero, Vec3.zero, Vec3.zero, ⟨0, 0, 1⟩⟩
  have hw_x := h ⟨Vec3.zero, Vec3.zero, ⟨1, 0, 0⟩, Vec3.zero⟩
  have hw_y := h ⟨Vec3.zero, Vec3.zero, ⟨0, 1, 0⟩, Vec3.zero⟩
  have hw_z := h ⟨Vec3.zero, Vec3.zero, ⟨0, 0, 1⟩, Vec3.zero⟩
  dsimp [omegaNormal, Vec3.dot, Vec3.zero] at hq_x hq_y hq_z hp_x hp_y hp_z hn_x hn_y hn_z hw_x hw_y hw_z
  have eqx : v1.dq.x = 0 := by
    calc v1.dq.x = (v1.dq.x * 1 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hq_x
  have eqy : v1.dq.y = 0 := by
    calc v1.dq.y = (v1.dq.x * 0 + v1.dq.y * 1 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hq_y
  have eqz : v1.dq.z = 0 := by
    calc v1.dq.z = (v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 1 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hq_z
  have epx : v1.dp.x = 0 := by
    calc v1.dp.x = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (1 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hp_x]
    _ = 0 := by ring
  have epy : v1.dp.y = 0 := by
    calc v1.dp.y = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 1 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hp_y]
    _ = 0 := by ring
  have epz : v1.dp.z = 0 := by
    calc v1.dp.z = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 1 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hp_z]
    _ = 0 := by ring
  have enx : v1.dn.x = 0 := by
    calc v1.dn.x = (v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 1 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hn_x
  have eny : v1.dn.y = 0 := by
    calc v1.dn.y = (v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 1 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hn_y
  have enz : v1.dn.z = 0 := by
    calc v1.dn.z = (v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 1 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hn_z
  have ewx : v1.dw.x = 0 := by
    calc v1.dw.x = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (1 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hw_x]
    _ = 0 := by ring
  have ewy : v1.dw.y = 0 := by
    calc v1.dw.y = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 1 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hw_y]
    _ = 0 := by ring
  have ewz : v1.dw.z = 0 := by
    calc v1.dw.z = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 1 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hw_z]
    _ = 0 := by ring
  rcases v1 with ⟨dq, dp, dn, dw⟩
  rcases dq with ⟨qx, qy, qz⟩
  rcases dp with ⟨px, py, pz⟩
  rcases dn with ⟨nx, ny, nz⟩
  rcases dw with ⟨wx, wy, wz⟩
  dsimp at eqx eqy eqz epx epy epz enx eny enz ewx ewy ewz
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [eqx, eqy, eqz]; rfl
  · rw [epx, epy, epz]; rfl
  · rw [enx, eny, enz]; rfl
  · rw [ewx, ewy, ewz]; rfl

end ReesTheorems

-- === 03_CanonicalSplitting_Defs.lean & 04_CanonicalSplitting_Theorems.lean ===
namespace ConormalSplitting

structure SplittingState (Base : Type u) (Fiber : Type u) where
  base : Base
  fiber : Fiber

structure ConormalIntegrator (k : Type u) (Base : Type u) (Fiber : Type u) where
  stepNL : k → SplittingState Base Fiber → SplittingState Base Fiber
  stepG  : k → SplittingState Base Fiber → SplittingState Base Fiber
  stepNL_base_id : ∀ (h : k) (s : SplittingState Base Fiber), (stepNL h s).base = s.base

def strangComposition {k : Type u} {Base Fiber : Type u}
    (I : ConormalIntegrator k Base Fiber)
    (h_half h : k) (s : SplittingState Base Fiber) : SplittingState Base Fiber :=
  I.stepNL h_half (I.stepG h (I.stepNL h_half s))

def PreservesInvariant {k : Type u} {Base Fiber : Type u}
    (step : k → SplittingState Base Fiber → SplittingState Base Fiber)
    (C : SplittingState Base Fiber → k) : Prop :=
  ∀ (h : k) (s : SplittingState Base Fiber), C (step h s) = C s

structure IsOneParameterFlow {k : Type u} [AddGroup k] {α : Type u} (flow : k → α → α) : Prop where
  group_law : ∀ (h₁ h₂ : k) (x : α), flow h₁ (flow h₂ x) = flow (h₁ + h₂) x
  zero_id   : ∀ (x : α), flow 0 x = x

theorem stepNL_preserves_base_casimir {k : Type u} {Base Fiber : Type u}
    (I : ConormalIntegrator k Base Fiber)
    (c : Base → k) :
    PreservesInvariant I.stepNL (fun s => c s.base) := by
  intro h s
  dsimp
  rw [I.stepNL_base_id h s]

theorem strang_exact_casimir_preservation {k : Type u} {Base Fiber : Type u}
    (I : ConormalIntegrator k Base Fiber)
    (C : SplittingState Base Fiber → k)
    (hNL : PreservesInvariant I.stepNL C)
    (hG  : PreservesInvariant I.stepG C)
    (h_half h : k) (s : SplittingState Base Fiber) :
    C (strangComposition I h_half h s) = C s := by
  unfold strangComposition
  rw [hNL, hG, hNL]

theorem strang_preserves_base_casimir {k : Type u} {Base Fiber : Type u}
    (I : ConormalIntegrator k Base Fiber)
    (c : Base → k)
    (hG : PreservesInvariant I.stepG (fun s => c s.base))
    (h_half h : k) (s : SplittingState Base Fiber) :
    c (strangComposition I h_half h s).base = c s.base :=
  strang_exact_casimir_preservation I (fun s => c s.base)
    (stepNL_preserves_base_casimir I c) hG h_half h s

theorem flow_group_inverse {k : Type u} [AddGroup k] {α : Type u} {flow : k → α → α}
    (grp : IsOneParameterFlow flow) (h : k) (x : α) :
    flow (-h) (flow h x) = x := by
  rw [grp.group_law, neg_add_cancel, grp.zero_id]

theorem strang_time_reversibility {k : Type u} [AddGroup k] {Base Fiber : Type u}
    (I : ConormalIntegrator k Base Fiber)
    (grpNL : IsOneParameterFlow I.stepNL)
    (grpG  : IsOneParameterFlow I.stepG)
    (h_half h : k) (s : SplittingState Base Fiber) :
    strangComposition I (-h_half) (-h) (strangComposition I h_half h s) = s := by
  unfold strangComposition
  set s1 := I.stepNL h_half s
  set s2 := I.stepG h s1
  set s3 := I.stepNL h_half s2
  have h_kick1 : I.stepNL (-h_half) s3 = s2 := by
    dsimp [s3]
    exact flow_group_inverse grpNL h_half s2
  rw [h_kick1]
  have h_groupoid : I.stepG (-h) s2 = s1 := by
    dsimp [s2]
    exact flow_group_inverse grpG h s1
  rw [h_groupoid]
  have h_kick2 : I.stepNL (-h_half) s1 = s := by
    dsimp [s1]
    exact flow_group_inverse grpNL h_half s
  exact h_kick2

end ConormalSplitting

-- === 05_HeavyTop_Defs.lean & 06_HeavyTop_Theorems.lean ===
namespace HeavyTop

structure HeavyTopState (k : Type u) where
  Pi : Vec3 k
  Gamma : Vec3 k

def casimir1 {k : Type u} [CommRing k] (s : HeavyTopState k) : k :=
  Vec3.dot s.Gamma s.Gamma

def casimir2 {k : Type u} [CommRing k] (s : HeavyTopState k) : k :=
  Vec3.dot s.Pi s.Gamma

def conormalTorqueKick {k : Type u} [CommRing k]
    (c : k) (chi : Vec3 k) (h : k) (s : HeavyTopState k) : HeavyTopState k :=
  let torque := Vec3.smul c (Vec3.cross s.Gamma chi)
  ⟨Vec3.add s.Pi (Vec3.smul h torque), s.Gamma⟩

structure SO3Isometry (k : Type u) [CommRing k] where
  rot : Vec3 k → Vec3 k
  preserves_dot : ∀ (u v : Vec3 k), Vec3.dot (rot u) (rot v) = Vec3.dot u v

def zAxisRot {k : Type u} [CommRing k] (c s : k) (v : Vec3 k) : Vec3 k :=
  ⟨c * v.x - s * v.y, s * v.x + c * v.y, v.z⟩

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

def rodriguesRot {k : Type u} [CommRing k] (n : Vec3 k) (c s : k) (v : Vec3 k) : Vec3 k :=
  Vec3.add (Vec3.add (Vec3.smul c v) (Vec3.smul s (Vec3.cross n v)))
           (Vec3.smul ((1 - c) * Vec3.dot n v) n)

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

def groupoidRotationStep {k : Type u} [CommRing k]
    (R : SO3Isometry k) (s : HeavyTopState k) : HeavyTopState k :=
  ⟨R.rot s.Pi, R.rot s.Gamma⟩

theorem torque_kick_preserves_casimir1 {k : Type u} [CommRing k]
    (c : k) (chi : Vec3 k) (h : k) (s : HeavyTopState k) :
    casimir1 (conormalTorqueKick c chi h s) = casimir1 s := by
  dsimp [conormalTorqueKick, casimir1]

theorem cross_dot_self {k : Type u} [CommRing k] (a b : Vec3 k) :
    Vec3.dot (Vec3.cross a b) a = 0 := by
  dsimp [Vec3.dot, Vec3.cross]
  ring

theorem torque_kick_preserves_casimir2 {k : Type u} [CommRing k]
    (c : k) (chi : Vec3 k) (h : k) (s : HeavyTopState k) :
    casimir2 (conormalTorqueKick c chi h s) = casimir2 s := by
  dsimp [conormalTorqueKick, casimir2, Vec3.dot, Vec3.add, Vec3.smul, Vec3.cross]
  ring

theorem groupoid_rot_preserves_casimir1 {k : Type u} [CommRing k]
    (R : SO3Isometry k) (s : HeavyTopState k) :
    casimir1 (groupoidRotationStep R s) = casimir1 s := by
  dsimp [groupoidRotationStep, casimir1]
  exact R.preserves_dot s.Gamma s.Gamma

theorem groupoid_rot_preserves_casimir2 {k : Type u} [CommRing k]
    (R : SO3Isometry k) (s : HeavyTopState k) :
    casimir2 (groupoidRotationStep R s) = casimir2 s := by
  dsimp [groupoidRotationStep, casimir2]
  exact R.preserves_dot s.Pi s.Gamma

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

#print axioms deformation_limit_zero
#print axioms omegaNormal_skew
#print axioms omegaNormal_self
#print axioms transverse_conormal_flow_symplectic
#print axioms omegaNormal_nondegenerate
#print axioms ConormalSplitting.stepNL_preserves_base_casimir
#print axioms ConormalSplitting.strang_exact_casimir_preservation
#print axioms ConormalSplitting.strang_preserves_base_casimir
#print axioms ConormalSplitting.flow_group_inverse
#print axioms ConormalSplitting.strang_time_reversibility
#print axioms HeavyTop.zAxisRotationIsometry
#print axioms HeavyTop.rodriguesIsometry
#print axioms HeavyTop.torque_kick_preserves_casimir1
#print axioms HeavyTop.cross_dot_self
#print axioms HeavyTop.torque_kick_preserves_casimir2
#print axioms HeavyTop.groupoid_rot_preserves_casimir1
#print axioms HeavyTop.groupoid_rot_preserves_casimir2
#print axioms HeavyTop.heavy_top_strang_casimirs_exact
#print axioms HeavyTop.heavy_top_strang_rodrigues_casimirs_exact
