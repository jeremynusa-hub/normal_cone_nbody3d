-- Dependencies: 01_ConfigurationSpace_Defs (Vec3, AmbientConfiguration)
universe u

def BinaryCollision {N : Nat} {k : Type u} (q : Fin N → Vec3 k) (i j : Fin N) : Prop :=
  q i = q j

def InCollisionLocus {N : Nat} {k : Type u} (q : Fin N → Vec3 k) : Prop :=
  ∃ (i j : Fin N), i ≠ j ∧ BinaryCollision q i j

def NonCoincident {N : Nat} {k : Type u} (q : Fin N → Vec3 k) : Prop :=
  ∀ (i j : Fin N), i ≠ j → q i ≠ q j
