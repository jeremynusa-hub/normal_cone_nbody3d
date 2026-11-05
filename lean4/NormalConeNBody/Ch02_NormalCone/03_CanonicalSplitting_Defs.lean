-- Dependencies: 01_ConfigurationSpace_Defs, 01_ReesDeformation_Defs
import Mathlib.Tactic.Ring

universe u

namespace ConormalSplitting

/-- Phase state of a dynamical system decomposed along the normal cone fibration.
    `base` represents generalized coordinates on the symmetry/singular locus Z.
    `fiber` represents conjugate transverse momenta in N_{Z/M}^*. -/
structure SplittingState (Base : Type u) (Fiber : Type u) where
  base : Base
  fiber : Fiber

/-- The canonical normal cone splitting operators:
    - `stepNL`: exact abelian translation along the normal fiber (normal fiber kick).
    - `stepG`: exact symplectic groupoid exponential on the base manifold. -/
structure ConormalIntegrator (k : Type u) (Base : Type u) (Fiber : Type u) where
  stepNL : k → SplittingState Base Fiber → SplittingState Base Fiber
  stepG  : k → SplittingState Base Fiber → SplittingState Base Fiber
  -- Transverse kicks leave the base coordinate invariant:
  stepNL_base_id : ∀ (h : k) (s : SplittingState Base Fiber), (stepNL h s).base = s.base

/-- Second-order symmetric Strang composition operator:
    S₂(h) = Φ_{h/2}^{NL} ∘ Φ_h^G ∘ Φ_{h/2}^{NL}. -/
def strangComposition {k : Type u} {Base Fiber : Type u}
    (I : ConormalIntegrator k Base Fiber)
    (h_half h : k) (s : SplittingState Base Fiber) : SplittingState Base Fiber :=
  I.stepNL h_half (I.stepG h (I.stepNL h_half s))

/-- Predicate for exact scalar invariant preservation (e.g. Poisson Casimirs). -/
def PreservesInvariant {k : Type u} {Base Fiber : Type u}
    (step : k → SplittingState Base Fiber → SplittingState Base Fiber)
    (C : SplittingState Base Fiber → k) : Prop :=
  ∀ (h : k) (s : SplittingState Base Fiber), C (step h s) = C s

/-- Axiomatic 1-parameter flow group property for time-reversibility proofs. -/
structure IsOneParameterFlow {k : Type u} [AddGroup k] {α : Type u} (flow : k → α → α) : Prop where
  group_law : ∀ (h₁ h₂ : k) (x : α), flow h₁ (flow h₂ x) = flow (h₁ + h₂) x
  zero_id   : ∀ (x : α), flow 0 x = x

end ConormalSplitting
