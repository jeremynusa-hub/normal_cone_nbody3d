-- Dependencies: 01_ConfigurationSpace_Defs, 01_ReesDeformation_Defs, 03_CanonicalSplitting_Defs
import Mathlib.Tactic.Ring

universe u

namespace ConormalSplitting

/-- Constructive Base-Casimir Invariance under Normal Fiber Kicks:
    Because the transverse normal fiber kick `stepNL` leaves the base coordinate invariant
    (`stepNL_base_id`), any observable depending purely on the base manifold `Z` is
    automatically and unconditionally preserved by `stepNL`. -/
theorem stepNL_preserves_base_casimir {k : Type u} {Base Fiber : Type u}
    (I : ConormalIntegrator k Base Fiber)
    (c : Base → k) :
    PreservesInvariant I.stepNL (fun s => c s.base) := by
  intro h s
  dsimp
  rw [I.stepNL_base_id h s]

/-- Fundamental Theorem: Unconditional Casimir Preservation under Normal Cone Strang Splitting.
    If the normal fiber momentum kick Φ^{NL} and the groupoid exponential Φ^G both preserve
    a scalar Casimir functional C, then their symmetric composition S₂(h) unconditionally
    preserves C without numerical drift or coordinate projection. -/
theorem strang_exact_casimir_preservation {k : Type u} {Base Fiber : Type u}
    (I : ConormalIntegrator k Base Fiber)
    (C : SplittingState Base Fiber → k)
    (hNL : PreservesInvariant I.stepNL C)
    (hG  : PreservesInvariant I.stepG C)
    (h_half h : k) (s : SplittingState Base Fiber) :
    C (strangComposition I h_half h s) = C s := by
  unfold strangComposition
  rw [hNL, hG, hNL]

/-- Corollary: Any base-manifold invariant `c(s.base)` preserved by the groupoid flow `stepG`
    is unconditionally preserved by the full Strang composition `S₂(h)` by virtue of
    `stepNL_base_id`, without requiring a separate invariance hypothesis on `stepNL`. -/
theorem strang_preserves_base_casimir {k : Type u} {Base Fiber : Type u}
    (I : ConormalIntegrator k Base Fiber)
    (c : Base → k)
    (hG : PreservesInvariant I.stepG (fun s => c s.base))
    (h_half h : k) (s : SplittingState Base Fiber) :
    c (strangComposition I h_half h s).base = c s.base :=
  strang_exact_casimir_preservation I (fun s => c s.base)
    (stepNL_preserves_base_casimir I c) hG h_half h s

/-- Inversion Lemma: Any 1-parameter flow satisfying the group law admits an exact time-inverse:
    Φ(-h) ∘ Φ(h) = id. -/
theorem flow_group_inverse {k : Type u} [AddGroup k] {α : Type u} {flow : k → α → α}
    (grp : IsOneParameterFlow flow) (h : k) (x : α) :
    flow (-h) (flow h x) = x := by
  rw [grp.group_law, neg_add_cancel, grp.zero_id]

/-- Fundamental Theorem: Exact Time-Reversibility of Symmetric Normal Cone Splitting.
    The Strang composition S₂(h) = Φ_{h/2}^{NL} ∘ Φ_h^G ∘ Φ_{h/2}^{NL} is an exact
    time-reversible involution: S₂(-h) ∘ S₂(h) = id. -/
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
