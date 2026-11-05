/-
  Project: Unsoftened Symplectic Integration of the 3D Newtonian N-Body Problem
  File: chapters/ch04_fulton_macpherson_splitting/sec03_synchronous_supernodes/code/03_AtlasTransitionSymplecticAntiChattering_Theorems.lean
  Author: Jérémy Nusa
  Theoretical Framework: Jacobi Symplectic Geometry Sp(6N-6, ℝ) / Hysteresis Anti-Chattering Criterion
  Verification Standard: Constructive 1-Jet Kernel (zero sorry, axiom-free)

  Component: Theorem 4.4 — Symplectic Atlas Transitions and Chattering Elimination
  Theoretical Derivation: Section 4.3 (Theorem 4.4 Symplectic Atlas Transitions)
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

namespace UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes

/-- 2x2 linear Jacobi coordinate transition matrix between adjacent clustering trees:
    M = [[alpha, beta], [gamma, delta]]. -/
structure JacobiTransitionMatrix2x2 where
  alpha : ℝ
  beta : ℝ
  gamma : ℝ
  delta : ℝ

/-- Determinant of the 2x2 transition matrix. -/
def det2x2 (M : JacobiTransitionMatrix2x2) : ℝ :=
  M.alpha * M.delta - M.beta * M.gamma

/-- Membership predicate in SL(2, ℝ): det(M) = 1. -/
def isInSL2 (M : JacobiTransitionMatrix2x2) : Prop :=
  det2x2 M = 1

/-- Explicit mass-dependent 2x2 Jacobi coordinate transition matrix M_coord (Theorem 4.4)
    associated with the elementary binary tree rotation ((S₁, S₂), S₃) ↔ (S₁, (S₂, S₃)):
      alpha = -m₁ / (m₁ + m₂),                        beta  = 1,
      gamma = -m₂ (m₁ + m₂ + m₃) / ((m₁+m₂)(m₂+m₃)),  delta = -m₃ / (m₂ + m₃). -/
noncomputable def jacobiTransitionMatrix (m1 m2 m3 : ℝ) : JacobiTransitionMatrix2x2 where
  alpha := -m1 / (m1 + m2)
  beta := 1
  gamma := -m2 * (m1 + m2 + m3) / ((m1 + m2) * (m2 + m3))
  delta := -m3 / (m2 + m3)

/-- Theorem (Unimodularity of the Jacobi Tree Rotation Matrix, Theorem 4.4).
    For any cluster masses m₁, m₂, m₃ with non-zero binary sub-masses m₁ + m₂ ≠ 0 and m₂ + m₃ ≠ 0,
    the Jacobi transition matrix M_coord belongs to SL(2, ℝ) (det M_coord = +1). -/
theorem jacobi_transition_matrix_in_sl2 (m1 m2 m3 : ℝ)
    (h12 : m1 + m2 ≠ 0) (h23 : m2 + m3 ≠ 0) :
    isInSL2 (jacobiTransitionMatrix m1 m2 m3) := by
  dsimp [isInSL2, det2x2, jacobiTransitionMatrix]
  field_simp [h12, h23]
  ring

/-- Lemma (Symplectic area preservation by SL(2, ℝ) transformations).
    For any vectors (x₁, x₂) and (y₁, y₂), transformation by M ∈ SL(2, ℝ)
    strictly preserves the canonical area form x₁ y₂ - x₂ y₁. -/
theorem sl2_preserves_symplectic_area (M : JacobiTransitionMatrix2x2) (hM : isInSL2 M)
    (x1 x2 y1 y2 : ℝ) :
    let x1' := M.alpha * x1 + M.beta * x2
    let x2' := M.gamma * x1 + M.delta * x2
    let y1' := M.alpha * y1 + M.beta * y2
    let y2' := M.gamma * y1 + M.delta * y2
    x1' * y2' - x2' * y1' = x1 * y2 - x2 * y1 := by
  intro x1' x2' y1' y2'
  dsimp [x1', x2', y1', y2']
  have h_det : M.alpha * M.delta - M.beta * M.gamma = 1 := hM
  calc
    (M.alpha * x1 + M.beta * x2) * (M.gamma * y1 + M.delta * y2) -
    (M.gamma * x1 + M.delta * x2) * (M.alpha * y1 + M.beta * y2)
      = (M.alpha * M.delta - M.beta * M.gamma) * (x1 * y2 - x2 * y1) := by ring
    _ = 1 * (x1 * y2 - x2 * y1) := by rw [h_det]
    _ = x1 * y2 - x2 * y1 := by ring

/-- Theorem (Cotangent Lift Symplectomorphism Ψ_rot = diag(M, (M⁻¹)ᵀ) ∈ Sp(4, ℝ), Theorem 4.4).
    When relative Jacobi positions transform by M ∈ SL(2, ℝ) and conjugate momenta transform
    by the inverse transpose (M⁻¹)ᵀ = [[delta, -gamma], [-beta, alpha]], the canonical
    Darboux 2-form ω₀((dr, dπ), (dr̄, dπ̄)) = (dr₁ dπ̄₁ - dr̄₁ dπ₁) + (dr₂ dπ̄₂ - dr̄₂ dπ₂)
    is identically preserved. -/
theorem cotangent_lift_preserves_darboux (M : JacobiTransitionMatrix2x2) (hM : isInSL2 M)
    (dr1 dr2 dpi1 dpi2 dr1_bar dr2_bar dpi1_bar dpi2_bar : ℝ) :
    let dr1' := M.alpha * dr1 + M.beta * dr2
    let dr2' := M.gamma * dr1 + M.delta * dr2
    let dpi1' := M.delta * dpi1 - M.gamma * dpi2
    let dpi2' := -M.beta * dpi1 + M.alpha * dpi2
    let dr1_bar' := M.alpha * dr1_bar + M.beta * dr2_bar
    let dr2_bar' := M.gamma * dr1_bar + M.delta * dr2_bar
    let dpi1_bar' := M.delta * dpi1_bar - M.gamma * dpi2_bar
    let dpi2_bar' := -M.beta * dpi1_bar + M.alpha * dpi2_bar
    (dr1' * dpi1_bar' - dr1_bar' * dpi1') + (dr2' * dpi2_bar' - dr2_bar' * dpi2') =
    (dr1 * dpi1_bar - dr1_bar * dpi1) + (dr2 * dpi2_bar - dr2_bar * dpi2) := by
  intro dr1' dr2' dpi1' dpi2' dr1_bar' dr2_bar' dpi1_bar' dpi2_bar'
  dsimp [dr1', dr2', dpi1', dpi2', dr1_bar', dr2_bar', dpi1_bar', dpi2_bar']
  have h_det : M.alpha * M.delta - M.beta * M.gamma = 1 := hM
  calc
    ((M.alpha * dr1 + M.beta * dr2) * (M.delta * dpi1_bar - M.gamma * dpi2_bar) -
     (M.alpha * dr1_bar + M.beta * dr2_bar) * (M.delta * dpi1 - M.gamma * dpi2)) +
    ((M.gamma * dr1 + M.delta * dr2) * (-M.beta * dpi1_bar + M.alpha * dpi2_bar) -
     (M.gamma * dr1_bar + M.delta * dr2_bar) * (-M.beta * dpi1 + M.alpha * dpi2))
      = (M.alpha * M.delta - M.beta * M.gamma) *
        ((dr1 * dpi1_bar - dr1_bar * dpi1) + (dr2 * dpi2_bar - dr2_bar * dpi2)) := by ring
    _ = 1 * ((dr1 * dpi1_bar - dr1_bar * dpi1) + (dr2 * dpi2_bar - dr2_bar * dpi2)) := by rw [h_det]
    _ = (dr1 * dpi1_bar - dr1_bar * dpi1) + (dr2 * dpi2_bar - dr2_bar * dpi2) := by ring

/-- Theorem: Absolute elimination of topological chattering under hysteresis γ = 0.6.
    For any mutual distance d > 0 and hysteresis threshold γ = 3/5 (γ² = 9/25 < 1),
    it is mathematically impossible to satisfy simultaneously the back-and-forth flip conditions:
    d₁₂ < γ * d_int  and  d_int < γ * d₁₂. -/
theorem hysteresis_anti_chattering_contradiction (d12 d_int gamma : ℝ)
    (h_d12 : 0 < d12) (h_gamma_pos : 0 ≤ gamma) (h_gamma_sq : gamma * gamma < 1)
    (h_flip_forward : d12 < gamma * d_int) (h_flip_backward : d_int < gamma * d12) :
    False := by
  have h_chain : d12 < (gamma * gamma) * d12 := by
    calc
      d12 < gamma * d_int := h_flip_forward
      _ ≤ gamma * (gamma * d12) := by
        have h_mul := mul_le_mul_of_nonneg_left (le_of_lt h_flip_backward) h_gamma_pos
        exact h_mul
      _ = (gamma * gamma) * d12 := by ring
  have h_contra : (1 - gamma * gamma) * d12 < 0 := by
    calc
      (1 - gamma * gamma) * d12 = d12 - (gamma * gamma) * d12 := by ring
      _ < 0 := sub_neg.mpr h_chain
  have h_pos_coeff : 0 < 1 - gamma * gamma := sub_pos.mpr h_gamma_sq
  have _h_prod_pos : 0 < (1 - gamma * gamma) * d12 := mul_pos h_pos_coeff h_d12
  linarith

/-- Theorem: Symplectic Atlas Transitions and Unconditional Chattering Elimination.
    1. Every unimodular Jacobi transition matrix in SL(2, ℝ) preserves the canonical symplectic structure.
    2. The asymmetric hysteresis criterion with γ = 0.6 (γ² = 0.36 < 1) strictly prevents any infinite
       topological chattering cycle (anti-chattering). -/
theorem atlas_transition_symplectic_anti_chattering (M : JacobiTransitionMatrix2x2) (hM : isInSL2 M)
    (x1 x2 y1 y2 : ℝ) (d12 d_int : ℝ) (h_d12 : 0 < d12) :
    -- 1. Symplectic transition preservation
    (let x1' := M.alpha * x1 + M.beta * x2
     let x2' := M.gamma * x1 + M.delta * x2
     let y1' := M.alpha * y1 + M.beta * y2
     let y2' := M.gamma * y1 + M.delta * y2
     x1' * y2' - x2' * y1' = x1 * y2 - x2 * y1) ∧
    -- 2. Non-existence of any chattering cycle for gamma = 3/5
    (¬ (d12 < (3 / 5 : ℝ) * d_int ∧ d_int < (3 / 5 : ℝ) * d12)) := by
  refine ⟨sl2_preserves_symplectic_area M hM x1 x2 y1 y2, ?_⟩
  intro ⟨hfwd, hback⟩
  have h_gamma_pos : (0 : ℝ) ≤ 3 / 5 := by linarith
  have h_gamma_sq : (3 / 5 : ℝ) * (3 / 5 : ℝ) < 1 := by
    norm_num
  exact hysteresis_anti_chattering_contradiction d12 d_int (3 / 5) h_d12 h_gamma_pos h_gamma_sq hfwd hback

end UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes
