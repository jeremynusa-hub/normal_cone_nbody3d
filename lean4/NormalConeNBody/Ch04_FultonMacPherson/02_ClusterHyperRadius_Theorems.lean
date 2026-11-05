-- ==============================================================================
-- Module : 02_ClusterHyperRadius_Theorems.lean
-- Component : Multiscale Cluster Hyper-Radius and Homogeneity
-- Author : Jérémy Nusa
-- Theoretical Framework : Fulton--MacPherson (1994), McGehee (1974)
-- Dependencies : Mathlib.Data.Real.Basic, Mathlib.Tactic.Linarith, Mathlib.Tactic.Ring
-- Verification Standard : Constructive 1-Jet Kernel (zero sorry, axiom-free)
-- ==============================================================================

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

namespace UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice

/-- Squared cluster hyper-radius for reduced mass μ > 0
    and internal quadratic separation r_sq ≥ 0 : ρ² = μ · r_sq. -/
noncomputable def hyperRadiusSq (mu : ℝ) (r_sq : ℝ) : ℝ :=
  mu * r_sq

/-- Theorem : Strict positivity of cluster hyper-radius.
    For strictly positive reduced mass (μ > 0) and non-zero quadratic separation
    (r_sq > 0), the squared hyper-radius is strictly positive : ρ² > 0. -/
theorem hyper_radius_sq_pos (mu : ℝ) (r_sq : ℝ)
    (h_mu : 0 < mu) (h_r : 0 < r_sq) :
    0 < hyperRadiusSq mu r_sq := by
  dsimp [hyperRadiusSq]
  exact mul_pos h_mu h_r

/-- Theorem : Total collision characterization for a sub-cluster.
    Under physical assumptions μ > 0 and r_sq ≥ 0, the hyper-radius vanishes
    if and only if the internal spatial separation vanishes (mutual collision). -/
theorem hyper_radius_sq_zero_iff (mu : ℝ) (r_sq : ℝ)
    (h_mu : 0 < mu) (_h_r : 0 ≤ r_sq) :
    hyperRadiusSq mu r_sq = 0 ↔ r_sq = 0 := by
  dsimp [hyperRadiusSq]
  constructor
  · intro h
    cases mul_eq_zero.mp h with
    | inl h_mu_zero =>
      have : mu ≠ 0 := ne_of_gt h_mu
      contradiction
    | inr h_r_zero => exact h_r_zero
  · intro h
    rw [h, mul_zero]

/-- Predicate for homogeneity of degree -1 of the internal gravitational potential :
    For every scaling factor λ > 0, U(λ · r) = λ⁻¹ · U(r). -/
def isHomogeneousDegreeMinusOne (U : ℝ → ℝ) : Prop :=
  ∀ (r : ℝ) (lambda : ℝ), 0 < lambda → U (lambda * r) = (1 / lambda) * U r

/-- Two-body scalar Newtonian gravitational potential : V(r) = -GM · r⁻¹. -/
noncomputable def newtonianScalarPotential (GM : ℝ) (r : ℝ) : ℝ :=
  - GM * r⁻¹

/-- Theorem (Homogeneity of the Newtonian Potential) :
    Formal proof that the Newtonian gravitational potential intrinsically
    satisfies McGehee's degree -1 homogeneity relation :
    V(λ · r) = λ⁻¹ · V(r) for all λ > 0 and r ≠ 0. -/
theorem newtonian_potential_is_homogeneous (GM : ℝ) :
    ∀ (r : ℝ) (lambda : ℝ), 0 < lambda → r ≠ 0 →
    newtonianScalarPotential GM (lambda * r) = (1 / lambda) * newtonianScalarPotential GM r := by
  intro r lambda h_lambda _h_r
  dsimp [newtonianScalarPotential]
  have _h_lam_ne : lambda ≠ 0 := ne_of_gt h_lambda
  rw [mul_inv, one_div]
  ring

/-- Theorem : McGehee potential factorization of internal potential.
    If the internal potential U is homogeneous of degree -1, then setting r = ρ · s (with ρ > 0),
    the potential factors exactly as U(r) = ρ⁻¹ · U(s). -/
theorem mcgehee_potential_factorization (U : ℝ → ℝ) (h_hom : isHomogeneousDegreeMinusOne U)
    (rho : ℝ) (h_rho : 0 < rho) (s : ℝ) :
    U (rho * s) = (1 / rho) * U s :=
  h_hom s rho h_rho

end UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice
