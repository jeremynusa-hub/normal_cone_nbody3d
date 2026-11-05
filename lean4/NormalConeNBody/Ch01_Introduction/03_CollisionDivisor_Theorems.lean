-- Dependencies: 01_ConfigurationSpace_Defs, 02_CollisionDivisor_Defs
universe u

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
