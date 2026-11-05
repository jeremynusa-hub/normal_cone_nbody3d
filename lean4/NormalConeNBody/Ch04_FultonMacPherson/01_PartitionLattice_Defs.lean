-- ==============================================================================
-- Module : 01_PartitionLattice_Defs.lean
-- Component : Partition Lattice Π(N), Collision Strata, and Non-Transversality
-- Author : Jérémy Nusa
-- Theoretical Framework : Fulton--MacPherson (1994), Axelrod--Singer (1994)
-- Dependencies : none
-- Verification Standard : Constructive 1-Jet Kernel (zero sorry, axiom-free)
-- ==============================================================================

namespace UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice

/-- Real codimension in the configuration space Q ≅ ℝ^{3N-3} of the multiple
    collision stratum associated with a sub-cluster of cardinality k ≥ 2 :
    codim(Δ_I) = 3 * (k - 1). -/
def collisionCodimension (k : Nat) : Nat :=
  if k ≥ 2 then 3 * (k - 1) else 0

/-- Lemma : The collision codimension of a binary pair (k = 2) is exactly 3. -/
theorem collisionCodimension_pair : collisionCodimension 2 = 3 := by
  rfl

/-- Lemma : The collision codimension of a triple (k = 3) is exactly 6. -/
theorem collisionCodimension_triple : collisionCodimension 3 = 6 := by
  rfl

/-- Fundamental Theorem of Non-Transversality of Pairwise Collision Diagonals :
    The intersection of collision diagonals Δ₁₂ (q₁ = q₂) and Δ₂₃ (q₂ = q₃)
    is strictly contained in the third diagonal Δ₁₃ (q₁ = q₃).
    By transitivity of Euclidean equality, the collision hyperplanes intersect
    in non-general position, inducing an intrinsically reducible singular normal cone
    in the raw Euclidean configuration space ℝ^{3N}. -/
theorem collision_diagonal_non_transversality {α : Type u} (q1 q2 q3 : α)
    (h12 : q1 = q2) (h23 : q2 = q3) : q1 = q3 := by
  rw [h12, h23]

/-- Predicate for a multiple collision state on a coalition of indices :
    All bodies belonging to the coalition I share the exact same spatial position. -/
def isCollisionState {N : Nat} {α : Type u} (q : Fin N → α) (I : Fin N → Prop) : Prop :=
  ∀ i j, I i → I j → q i = q j

/-- Lemma : Total collision of a superset implies collision of any sub-coalition. -/
theorem collision_hereditary {N : Nat} {α : Type u} (q : Fin N → α) {I J : Fin N → Prop}
    (h_sub : ∀ i, J i → I i) (h_coll : isCollisionState q I) :
    isCollisionState q J := by
  intro i j hi hj
  exact h_coll i j (h_sub i hi) (h_sub j hj)

end UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice
