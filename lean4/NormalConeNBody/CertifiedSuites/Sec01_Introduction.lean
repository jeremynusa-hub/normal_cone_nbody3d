import Mathlib.Tactic.Ring
universe u

-- ============================================================================
-- 01_ConfigurationSpace_Defs.lean
-- ============================================================================
@[ext]
structure Vec3 (k : Type u) where
  x : k
  y : k
  z : k

namespace Vec3
variable {k : Type u} [CommRing k]
def zero : Vec3 k := ⟨0, 0, 0⟩
def add (v1 v2 : Vec3 k) : Vec3 k := ⟨v1.x + v2.x, v1.y + v2.y, v1.z + v2.z⟩
def smul (c : k) (v : Vec3 k) : Vec3 k := ⟨c * v.x, c * v.y, c * v.z⟩
def dot (v1 v2 : Vec3 k) : k := v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
def cross (v1 v2 : Vec3 k) : Vec3 k := ⟨v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x⟩
end Vec3

def AmbientConfiguration (N : Nat) (k : Type u) := Fin N → Vec3 k

structure PhaseState3D (N : Nat) (k : Type u) where
  q : Fin N → Vec3 k
  p : Fin N → Vec3 k
  m : Fin N → k

def SpatialTranslation {N : Nat} {k : Type u} [CommRing k] (q : AmbientConfiguration N k) (a : Vec3 k) : AmbientConfiguration N k :=
  fun i => Vec3.add (q i) a

def sumFin {k : Type u} [CommRing k] : (n : Nat) → (Fin n → k) → k
  | 0, _ => 0
  | n + 1, f => f 0 + sumFin n (fun i => f ⟨i.val + 1, Nat.succ_lt_succ i.isLt⟩)

def RelativeConfiguration (N : Nat) {k : Type u} [CommRing k] (m : Fin N → k) :=
  { q : AmbientConfiguration N k //
    (sumFin N (fun i => m i * (q i).x)) = 0 ∧
    (sumFin N (fun i => m i * (q i).y)) = 0 ∧
    (sumFin N (fun i => m i * (q i).z)) = 0 }

-- ============================================================================
-- 02_CollisionDivisor_Defs.lean
-- ============================================================================
def BinaryCollision {N : Nat} {k : Type u} (q : Fin N → Vec3 k) (i j : Fin N) : Prop :=
  q i = q j

def InCollisionLocus {N : Nat} {k : Type u} (q : Fin N → Vec3 k) : Prop :=
  ∃ (i j : Fin N), i ≠ j ∧ BinaryCollision q i j

def NonCoincident {N : Nat} {k : Type u} (q : Fin N → Vec3 k) : Prop :=
  ∀ (i j : Fin N), i ≠ j → q i ≠ q j

-- ============================================================================
-- 03_CollisionDivisor_Theorems.lean
-- ============================================================================
theorem non_coincident_pair_symm {N : Nat} {k : Type u} (q : Fin N → Vec3 k) (i j : Fin N) :
    q i ≠ q j ↔ q j ≠ q i := by
  constructor
  · intro h hij; exact h hij.symm
  · intro h hji; exact h hji.symm

theorem regular_iff_not_collision {N : Nat} {k : Type u} (q : Fin N → Vec3 k) :
    NonCoincident q ↔ ¬ InCollisionLocus q := by
  dsimp [NonCoincident, InCollisionLocus, BinaryCollision]
  constructor
  · intro h_reg ⟨i, j, h_ne, h_eq⟩; exact (h_reg i j h_ne) h_eq
  · intro h_ncol i j h_ne h_eq; apply h_ncol; exact ⟨i, j, h_ne, h_eq⟩

-- ============================================================================
-- 04_SymplecticForm_Defs.lean
-- ============================================================================
structure TangentVector3D (N : Nat) (k : Type u) where
  dq : Fin N → Vec3 k
  dp : Fin N → Vec3 k

variable {N : Nat} {k : Type u} [CommRing k]

def omega0At (v1 v2 : TangentVector3D N k) (i : Fin N) : k :=
  Vec3.dot (v1.dq i) (v2.dp i) - Vec3.dot (v2.dq i) (v1.dp i)

def omega0 (v1 v2 : TangentVector3D N k) : k := sumFin N (fun i => omega0At v1 v2 i)

-- ============================================================================
-- 05_SymplecticForm_Theorems.lean
-- ============================================================================
theorem sumFin_zero {k : Type u} [CommRing k] (n : Nat) (f : Fin n → k) (h : ∀ i, f i = 0) : sumFin n f = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => dsimp [sumFin]; rw [h 0, ih _ (fun i => h _)]; ring

theorem sumFin_neg {k : Type u} [CommRing k] (n : Nat) (f g : Fin n → k) (h : ∀ i, f i = - g i) : sumFin n f = - sumFin n g := by
  induction n with
  | zero => dsimp [sumFin]; ring
  | succ n ih => dsimp [sumFin]; rw [h 0, ih _ _ (fun i => h _)]; ring

theorem omega0_skew (v1 v2 : TangentVector3D N k) : omega0 v1 v2 = - omega0 v2 v1 := by
  dsimp [omega0]; apply sumFin_neg; intro i; dsimp [omega0At]; ring

theorem omega0_self (v : TangentVector3D N k) : omega0 v v = 0 := by
  dsimp [omega0]; apply sumFin_zero; intro i; dsimp [omega0At]; ring

theorem omega0At_nondeg (dq dp : Vec3 k)
    (h : ∀ (dq' dp' : Vec3 k), Vec3.dot dq dp' - Vec3.dot dq' dp = 0) :
    dq = Vec3.zero ∧ dp = Vec3.zero := by
  have hqx := h Vec3.zero ⟨1, 0, 0⟩
  have hqy := h Vec3.zero ⟨0, 1, 0⟩
  have hqz := h Vec3.zero ⟨0, 0, 1⟩
  have hpx := h ⟨1, 0, 0⟩ Vec3.zero
  have hpy := h ⟨0, 1, 0⟩ Vec3.zero
  have hpz := h ⟨0, 0, 1⟩ Vec3.zero
  have ex : dq.x = 0 := by
    calc dq.x = (Vec3.dot dq ⟨1, 0, 0⟩ - Vec3.dot Vec3.zero dp) := by dsimp [Vec3.dot, Vec3.zero]; ring
    _ = 0 := hqx
  have ey : dq.y = 0 := by
    calc dq.y = (Vec3.dot dq ⟨0, 1, 0⟩ - Vec3.dot Vec3.zero dp) := by dsimp [Vec3.dot, Vec3.zero]; ring
    _ = 0 := hqy
  have ez : dq.z = 0 := by
    calc dq.z = (Vec3.dot dq ⟨0, 0, 1⟩ - Vec3.dot Vec3.zero dp) := by dsimp [Vec3.dot, Vec3.zero]; ring
    _ = 0 := hqz
  have epx : dp.x = 0 := by
    calc dp.x = - (Vec3.dot dq Vec3.zero - Vec3.dot ⟨1, 0, 0⟩ dp) := by dsimp [Vec3.dot, Vec3.zero]; ring
    _ = - 0 := by rw [hpx]
    _ = 0 := by ring
  have epy : dp.y = 0 := by
    calc dp.y = - (Vec3.dot dq Vec3.zero - Vec3.dot ⟨0, 1, 0⟩ dp) := by dsimp [Vec3.dot, Vec3.zero]; ring
    _ = - 0 := by rw [hpy]
    _ = 0 := by ring
  have epz : dp.z = 0 := by
    calc dp.z = - (Vec3.dot dq Vec3.zero - Vec3.dot ⟨0, 0, 1⟩ dp) := by dsimp [Vec3.dot, Vec3.zero]; ring
    _ = - 0 := by rw [hpz]
    _ = 0 := by ring
  rcases dq with ⟨qx, qy, qz⟩
  rcases dp with ⟨px, py, pz⟩
  dsimp at ex ey ez epx epy epz
  constructor
  · rw [ex, ey, ez]; rfl
  · rw [epx, epy, epz]; rfl

theorem omega0_tangent_nondeg (v : TangentVector3D N k)
    (h : ∀ i (dq' dp' : Vec3 k), Vec3.dot (v.dq i) dp' - Vec3.dot dq' (v.dp i) = 0) :
    (∀ i, v.dq i = Vec3.zero) ∧ (∀ i, v.dp i = Vec3.zero) := by
  constructor
  · intro i; exact (omega0At_nondeg (v.dq i) (v.dp i) (h i)).1
  · intro i; exact (omega0At_nondeg (v.dq i) (v.dp i) (h i)).2

theorem sumFin_single {k : Type u} [CommRing k] (n : Nat) (i0 : Fin n) (f : Fin n → k)
    (h : ∀ j, j ≠ i0 → f j = 0) : sumFin n f = f i0 := by
  induction n with
  | zero => exact Fin.elim0 i0
  | succ n ih =>
    dsimp [sumFin]
    by_cases h0 : (0 : Fin (n + 1)) = i0
    · subst h0
      have h_rest : sumFin n (fun j => f ⟨j.val + 1, Nat.succ_lt_succ j.isLt⟩) = 0 := by
        apply sumFin_zero
        intro j
        apply h
        intro h_eq
        have h_val : j.val + 1 = 0 := congrArg Fin.val h_eq
        omega
      rw [h_rest, add_zero]
    · have hf0 : f 0 = 0 := h 0 h0
      have hi0_pos : 0 < i0.val := Nat.pos_of_ne_zero (fun h_val => h0 (Fin.ext h_val.symm))
      let i0_pred : Fin n := ⟨i0.val - 1, by omega⟩
      have hi0_eq : (⟨i0_pred.val + 1, Nat.succ_lt_succ i0_pred.isLt⟩ : Fin (n + 1)) = i0 := by
        ext; dsimp [i0_pred]; omega
      have h_ind := ih i0_pred (fun j => f ⟨j.val + 1, Nat.succ_lt_succ j.isLt⟩) (by
        intro j hj
        apply h
        intro h_eq
        apply hj
        ext
        have h_val : j.val + 1 = i0.val := congrArg Fin.val h_eq
        dsimp [i0_pred]; omega)
      rw [hf0, zero_add, h_ind, hi0_eq]

def localizedTangent (i : Fin N) (dq' dp' : Vec3 k) : TangentVector3D N k where
  dq j := if j = i then dq' else Vec3.zero
  dp j := if j = i then dp' else Vec3.zero

theorem omega0_localized (v : TangentVector3D N k) (i : Fin N) (dq' dp' : Vec3 k) :
    omega0 v (localizedTangent i dq' dp') = Vec3.dot (v.dq i) dp' - Vec3.dot dq' (v.dp i) := by
  dsimp [omega0]
  rw [sumFin_single N i (fun j => omega0At v (localizedTangent i dq' dp') j)]
  · simp [omega0At, localizedTangent]
  · intro j hj
    simp [omega0At, localizedTangent, hj]
    dsimp [Vec3.dot, Vec3.zero]
    ring

theorem omega0_global_nondeg (v : TangentVector3D N k)
    (h : ∀ (w : TangentVector3D N k), omega0 v w = 0) :
    (∀ i, v.dq i = Vec3.zero) ∧ (∀ i, v.dp i = Vec3.zero) := by
  apply omega0_tangent_nondeg
  intro i dq' dp'
  rw [← omega0_localized v i dq' dp']
  exact h (localizedTangent i dq' dp')

-- ============================================================================
-- 06_BargmannAlgebra_Defs.lean
-- ============================================================================
inductive Axis3 | X | Y | Z deriving DecidableEq

namespace Axis3
variable {k : Type u} [CommRing k]
def unitVec (a : Axis3) : Vec3 k := match a with | Axis3.X => ⟨1, 0, 0⟩ | Axis3.Y => ⟨0, 1, 0⟩ | Axis3.Z => ⟨0, 0, 1⟩
def proj (a : Axis3) (v : Vec3 k) : k := match a with | Axis3.X => v.x | Axis3.Y => v.y | Axis3.Z => v.z
end Axis3

structure FormallySmoothObservable (N : Nat) (k : Type u) where
  val : PhaseState3D N k → k
  grad_q : PhaseState3D N k → Fin N → Vec3 k
  grad_p : PhaseState3D N k → Fin N → Vec3 k

def poissonBracketAt {N : Nat} {k : Type u} [CommRing k] (F G : FormallySmoothObservable N k) (s : PhaseState3D N k) (i : Fin N) : k :=
  Vec3.dot (F.grad_q s i) (G.grad_p s i) - Vec3.dot (F.grad_p s i) (G.grad_q s i)

def poissonBracketVal {N : Nat} {k : Type u} [CommRing k] (F G : FormallySmoothObservable N k) (s : PhaseState3D N k) : k :=
  sumFin N (fun i => poissonBracketAt F G s i)

namespace BargmannAlgebra
variable {N : Nat} {k : Type u} [CommRing k]

def linearMomentumComponent (a : Axis3) : FormallySmoothObservable N k where
  val s := sumFin N (fun i => Axis3.proj a (s.p i))
  grad_q _ _ := Vec3.zero
  grad_p _ _ := Axis3.unitVec a

def centerOfMassComponent (a : Axis3) (t : k) : FormallySmoothObservable N k where
  val s := sumFin N (fun i => s.m i * Axis3.proj a (s.q i) - t * Axis3.proj a (s.p i))
  grad_q s i := Vec3.smul (s.m i) (Axis3.unitVec a)
  grad_p _ _ := Vec3.smul (-t) (Axis3.unitVec a)

def totalMass (s : PhaseState3D N k) : k := sumFin N (fun i => s.m i)

def angularMomentumComponent (a : Axis3) : FormallySmoothObservable N k where
  val s := sumFin N (fun i => Axis3.proj a (Vec3.cross (s.q i) (s.p i)))
  grad_q s i := Vec3.cross (s.p i) (Axis3.unitVec a)
  grad_p s i := Vec3.cross (Axis3.unitVec a) (s.q i)

-- ============================================================================
-- 07_BargmannAlgebra_Theorems.lean
-- ============================================================================
theorem poisson_P_G_same (a : Axis3) (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (linearMomentumComponent a) (centerOfMassComponent a t) s = - totalMass s := by
  dsimp [poissonBracketVal, totalMass]
  have h_rw : (fun i => poissonBracketAt (linearMomentumComponent a) (centerOfMassComponent a t) s i) = (fun i => - s.m i) := by
    funext i; dsimp [poissonBracketAt, linearMomentumComponent, centerOfMassComponent]
    cases a <;> { dsimp [Axis3.unitVec, Vec3.zero, Vec3.smul, Vec3.dot]; ring }
  rw [h_rw]; exact sumFin_neg N (fun i => - s.m i) (fun i => s.m i) (fun _ => rfl)

theorem poisson_P_G_diff {a b : Axis3} (h : a ≠ b) (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (linearMomentumComponent a) (centerOfMassComponent b t) s = 0 := by
  dsimp [poissonBracketVal]; apply sumFin_zero; intro i
  dsimp [poissonBracketAt, linearMomentumComponent, centerOfMassComponent]
  cases a <;> cases b <;> try contradiction
  all_goals { dsimp [Axis3.unitVec, Vec3.zero, Vec3.smul, Vec3.dot]; ring }

theorem poisson_L_self (a : Axis3) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent a) (angularMomentumComponent a) s = 0 := by
  dsimp [poissonBracketVal]; apply sumFin_zero; intro i
  dsimp [poissonBracketAt]
  have h_comm : Vec3.dot ((angularMomentumComponent a).grad_q s i) ((angularMomentumComponent a).grad_p s i) =
                Vec3.dot ((angularMomentumComponent a).grad_p s i) ((angularMomentumComponent a).grad_q s i) := by
    dsimp [Vec3.dot]; ring
  rw [h_comm, sub_self]

theorem poisson_L_X_L_Y (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.X) (angularMomentumComponent Axis3.Y) s =
      (angularMomentumComponent Axis3.Z).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Axis3.proj]
  ring

theorem poisson_L_Y_L_Z (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Y) (angularMomentumComponent Axis3.Z) s =
      (angularMomentumComponent Axis3.X).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Axis3.proj]
  ring

theorem poisson_L_Z_L_X (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Z) (angularMomentumComponent Axis3.X) s =
      (angularMomentumComponent Axis3.Y).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Axis3.proj]
  ring

theorem poisson_L_P_same (a : Axis3) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent a) (linearMomentumComponent a) s = 0 := by
  dsimp [poissonBracketVal]; apply sumFin_zero; intro i
  dsimp [poissonBracketAt, angularMomentumComponent, linearMomentumComponent]
  cases a <;> { dsimp [Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.zero]; ring }

theorem poisson_L_X_P_Y (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.X) (linearMomentumComponent Axis3.Y) s =
      (linearMomentumComponent Axis3.Z).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, linearMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.zero, Axis3.proj]
  ring

theorem poisson_L_Y_P_Z (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Y) (linearMomentumComponent Axis3.Z) s =
      (linearMomentumComponent Axis3.X).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, linearMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.zero, Axis3.proj]
  ring

theorem poisson_L_Z_P_X (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Z) (linearMomentumComponent Axis3.X) s =
      (linearMomentumComponent Axis3.Y).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, linearMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.zero, Axis3.proj]
  ring

theorem poisson_L_G_same (a : Axis3) (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent a) (centerOfMassComponent a t) s = 0 := by
  dsimp [poissonBracketVal]; apply sumFin_zero; intro i
  dsimp [poissonBracketAt, angularMomentumComponent, centerOfMassComponent]
  cases a <;> { dsimp [Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.smul]; ring }

theorem poisson_L_X_G_Y (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.X) (centerOfMassComponent Axis3.Y t) s =
      (centerOfMassComponent Axis3.Z t).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, centerOfMassComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.smul, Axis3.proj]
  ring

theorem poisson_L_Y_G_Z (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Y) (centerOfMassComponent Axis3.Z t) s =
      (centerOfMassComponent Axis3.X t).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, centerOfMassComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.smul, Axis3.proj]
  ring

theorem poisson_L_Z_G_X (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Z) (centerOfMassComponent Axis3.X t) s =
      (centerOfMassComponent Axis3.Y t).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, centerOfMassComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.smul, Axis3.proj]
  ring
end BargmannAlgebra

-- ============================================================================
-- AXIOM HYGIENE VERIFICATION (#print axioms on all 25 theorems)
-- ============================================================================
#print axioms non_coincident_pair_symm
#print axioms regular_iff_not_collision
#print axioms sumFin_zero
#print axioms sumFin_neg
#print axioms omega0_skew
#print axioms omega0_self
#print axioms omega0At_nondeg
#print axioms omega0_tangent_nondeg
#print axioms sumFin_single
#print axioms omega0_localized
#print axioms omega0_global_nondeg
#print axioms BargmannAlgebra.poisson_P_G_same
#print axioms BargmannAlgebra.poisson_P_G_diff
#print axioms BargmannAlgebra.poisson_L_self
#print axioms BargmannAlgebra.poisson_L_X_L_Y
#print axioms BargmannAlgebra.poisson_L_Y_L_Z
#print axioms BargmannAlgebra.poisson_L_Z_L_X
#print axioms BargmannAlgebra.poisson_L_P_same
#print axioms BargmannAlgebra.poisson_L_X_P_Y
#print axioms BargmannAlgebra.poisson_L_Y_P_Z
#print axioms BargmannAlgebra.poisson_L_Z_P_X
#print axioms BargmannAlgebra.poisson_L_G_same
#print axioms BargmannAlgebra.poisson_L_X_G_Y
#print axioms BargmannAlgebra.poisson_L_Y_G_Z
#print axioms BargmannAlgebra.poisson_L_Z_G_X
