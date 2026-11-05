-- Dependencies: 01_ConfigurationSpace_Defs, 04_SymplecticForm_Defs
import Mathlib.Tactic.Ring
universe u

theorem sumFin_zero {k : Type u} [CommRing k] (n : Nat) (f : Fin n → k) (h : ∀ i, f i = 0) : sumFin n f = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => dsimp [sumFin]; rw [h 0, ih _ (fun i => h _)]; ring

theorem sumFin_neg {k : Type u} [CommRing k] (n : Nat) (f g : Fin n → k) (h : ∀ i, f i = - g i) : sumFin n f = - sumFin n g := by
  induction n with
  | zero => dsimp [sumFin]; ring
  | succ n ih => dsimp [sumFin]; rw [h 0, ih _ _ (fun i => h _)]; ring

variable {N : Nat} {k : Type u} [CommRing k]

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

