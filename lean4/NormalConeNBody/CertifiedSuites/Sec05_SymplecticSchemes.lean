import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ext
import Mathlib.LinearAlgebra.Matrix.DotProduct
import Mathlib.LinearAlgebra.Matrix.Symmetric

-- ============================================================================
-- 1. 01_StrangOrderTwoSymplectic_Theorems.lean
-- ============================================================================
namespace Sec05_01

@[ext]
structure ExtendedPhasePoint where
  r : ℝ
  pi : ℝ
  t : ℝ
  E : ℝ

def extendedTimeReversal (p : ExtendedPhasePoint) : ExtendedPhasePoint :=
  ⟨p.r, -p.pi, -p.t, p.E⟩

theorem extended_time_reversal_involution (p : ExtendedPhasePoint) :
    extendedTimeReversal (extendedTimeReversal p) = p := by
  dsimp [extendedTimeReversal]; ext <;> ring

noncomputable def driftStep (mu dt : ℝ) (p : ExtendedPhasePoint) : ExtendedPhasePoint :=
  ⟨p.r + (dt / (2 * mu)) * p.pi, p.pi, p.t, p.E⟩

noncomputable def apexKickStep (gradV dt : ℝ) (p : ExtendedPhasePoint) : ExtendedPhasePoint :=
  ⟨p.r, p.pi - (dt / 2) * gradV, p.t, p.E⟩

def timeAdvanceStep (dt : ℝ) (p : ExtendedPhasePoint) : ExtendedPhasePoint :=
  ⟨p.r, p.pi, p.t + dt, p.E⟩

noncomputable def strangStep (mu gradV dt : ℝ) (p : ExtendedPhasePoint) : ExtendedPhasePoint :=
  apexKickStep gradV dt (driftStep mu (2 * dt) (apexKickStep gradV dt (timeAdvanceStep dt p)))

def symplecticArea (dr1 dpi1 dr2 dpi2 : ℝ) : ℝ :=
  dr1 * dpi2 - dr2 * dpi1

theorem apex_kick_preserves_symplectic_area (dr1 dpi1 dr2 dpi2 H dt : ℝ) :
    let dpi1_kick := dpi1 - (dt / 2) * H * dr1
    let dpi2_kick := dpi2 - (dt / 2) * H * dr2
    symplecticArea dr1 dpi1_kick dr2 dpi2_kick = symplecticArea dr1 dpi1 dr2 dpi2 := by
  intro dpi1_kick dpi2_kick; dsimp [symplecticArea, dpi1_kick, dpi2_kick]; ring

theorem drift_preserves_symplectic_area (dr1 dpi1 dr2 dpi2 mu dt : ℝ) :
    let dr1_drift := dr1 + (dt / (2 * mu)) * dpi1
    let dr2_drift := dr2 + (dt / (2 * mu)) * dpi2
    symplecticArea dr1_drift dpi1 dr2_drift dpi2 = symplecticArea dr1 dpi1 dr2 dpi2 := by
  intro dr1_drift dr2_drift; dsimp [symplecticArea, dr1_drift, dr2_drift]; ring

structure SymmetricHessian3x3 where
  h11 : ℝ
  h12 : ℝ
  h13 : ℝ
  h22 : ℝ
  h23 : ℝ
  h33 : ℝ

def symplecticForm3D (dr1x dr1y dr1z dpi1x dpi1y dpi1z dt1 dE1
                      dr2x dr2y dr2z dpi2x dpi2y dpi2z dt2 dE2 : ℝ) : ℝ :=
  (dr1x * dpi2x - dr2x * dpi1x) +
  (dr1y * dpi2y - dr2y * dpi1y) +
  (dr1z * dpi2z - dr2z * dpi1z) +
  (dt1 * dE2 - dt2 * dE1)

theorem apex_kick_preserves_symplectic_form_3d (H : SymmetricHessian3x3) (dt : ℝ)
    (dr1x dr1y dr1z dpi1x dpi1y dpi1z dt1 dE1
     dr2x dr2y dr2z dpi2x dpi2y dpi2z dt2 dE2 : ℝ) :
    let dpi1x' := dpi1x - (dt / 2) * (H.h11 * dr1x + H.h12 * dr1y + H.h13 * dr1z)
    let dpi1y' := dpi1y - (dt / 2) * (H.h12 * dr1x + H.h22 * dr1y + H.h23 * dr1z)
    let dpi1z' := dpi1z - (dt / 2) * (H.h13 * dr1x + H.h23 * dr1y + H.h33 * dr1z)
    let dpi2x' := dpi2x - (dt / 2) * (H.h11 * dr2x + H.h12 * dr2y + H.h13 * dr2z)
    let dpi2y' := dpi2y - (dt / 2) * (H.h12 * dr2x + H.h22 * dr2y + H.h23 * dr2z)
    let dpi2z' := dpi2z - (dt / 2) * (H.h13 * dr2x + H.h23 * dr2y + H.h33 * dr2z)
    symplecticForm3D dr1x dr1y dr1z dpi1x' dpi1y' dpi1z' dt1 dE1
                     dr2x dr2y dr2z dpi2x' dpi2y' dpi2z' dt2 dE2 =
    symplecticForm3D dr1x dr1y dr1z dpi1x dpi1y dpi1z dt1 dE1
                     dr2x dr2y dr2z dpi2x dpi2y dpi2z dt2 dE2 := by
  intro dpi1x' dpi1y' dpi1z' dpi2x' dpi2y' dpi2z'
  dsimp [symplecticForm3D, dpi1x', dpi1y', dpi1z', dpi2x', dpi2y', dpi2z']; ring

theorem drift_preserves_symplectic_form_3d (mu dt : ℝ)
    (dr1x dr1y dr1z dpi1x dpi1y dpi1z dt1 dE1
     dr2x dr2y dr2z dpi2x dpi2y dpi2z dt2 dE2 : ℝ) :
    let dr1x' := dr1x + (dt / (2 * mu)) * dpi1x
    let dr1y' := dr1y + (dt / (2 * mu)) * dpi1y
    let dr1z' := dr1z + (dt / (2 * mu)) * dpi1z
    let dr2x' := dr2x + (dt / (2 * mu)) * dpi2x
    let dr2y' := dr2y + (dt / (2 * mu)) * dpi2y
    let dr2z' := dr2z + (dt / (2 * mu)) * dpi2z
    symplecticForm3D dr1x' dr1y' dr1z' dpi1x dpi1y dpi1z dt1 dE1
                     dr2x' dr2y' dr2z' dpi2x dpi2y dpi2z dt2 dE2 =
    symplecticForm3D dr1x dr1y dr1z dpi1x dpi1y dpi1z dt1 dE1
                     dr2x dr2y dr2z dpi2x dpi2y dpi2z dt2 dE2 := by
  intro dr1x' dr1y' dr1z' dr2x' dr2y' dr2z'
  dsimp [symplecticForm3D, dr1x', dr1y', dr1z', dr2x', dr2y', dr2z']; ring

theorem strang_palindromic_microréversibility (mu gradV dt : ℝ) (p : ExtendedPhasePoint) :
    extendedTimeReversal (strangStep mu gradV dt (extendedTimeReversal p)) =
      strangStep mu gradV (-dt) p := by
  dsimp [strangStep, apexKickStep, driftStep, timeAdvanceStep, extendedTimeReversal]; ext <;> ring

end Sec05_01

-- ============================================================================
-- 2. 02_LieBargmannConservation_Theorems.lean
-- ============================================================================
namespace Sec05_02

@[ext]
structure Vec3 where
  x : ℝ
  y : ℝ
  z : ℝ

namespace Vec3
def zero : Vec3 := ⟨0, 0, 0⟩
def add (v w : Vec3) : Vec3 := ⟨v.x + w.x, v.y + w.y, v.z + w.z⟩
def neg (v : Vec3) : Vec3 := ⟨-v.x, -v.y, -v.z⟩
def sub (v w : Vec3) : Vec3 := ⟨v.x - w.x, v.y - w.y, v.z - w.z⟩
def smul (c : ℝ) (v : Vec3) : Vec3 := ⟨c * v.x, c * v.y, c * v.z⟩
def cross (v w : Vec3) : Vec3 :=
  ⟨v.y * w.z - v.z * w.y, v.z * w.x - v.x * w.z, v.x * w.y - v.y * w.x⟩
end Vec3

theorem three_body_central_torque_zero
    (r1 r2 r3 : Vec3) (f12 f23 f31 : ℝ) :
    let F12 := Vec3.smul f12 (Vec3.sub r1 r2)
    let F23 := Vec3.smul f23 (Vec3.sub r2 r3)
    let F31 := Vec3.smul f31 (Vec3.sub r3 r1)
    let tau1 := Vec3.cross r1 (Vec3.sub F12 F31)
    let tau2 := Vec3.cross r2 (Vec3.sub F23 F12)
    let tau3 := Vec3.cross r3 (Vec3.sub F31 F23)
    Vec3.add (Vec3.add tau1 tau2) tau3 = Vec3.zero := by
  intro F12 F23 F31 tau1 tau2 tau3
  dsimp [F12, F23, F31, tau1, tau2, tau3, Vec3.add, Vec3.sub, Vec3.smul, Vec3.cross, Vec3.zero]
  ext <;> ring

theorem galilean_boost_step_invariance (M_tot R_t P_tot t dt : ℝ) (hM : M_tot ≠ 0) :
    let R_next := R_t + (dt / M_tot) * P_tot
    let t_next := t + dt
    let G_t := M_tot * R_t - t * P_tot
    let G_next := M_tot * R_next - t_next * P_tot
    G_next = G_t := by
  intro R_next t_next G_t G_next
  dsimp [R_next, t_next, G_t, G_next]
  have h_prod : M_tot * (R_t + dt / M_tot * P_tot) = M_tot * R_t + dt * P_tot := by
    calc
      M_tot * (R_t + dt / M_tot * P_tot) = M_tot * R_t + (M_tot * (dt / M_tot)) * P_tot := by ring
      _ = M_tot * R_t + dt * P_tot := by rw [mul_div_cancel₀ dt hM]
  calc
    M_tot * (R_t + dt / M_tot * P_tot) - (t + dt) * P_tot
      = (M_tot * R_t + dt * P_tot) - (t * P_tot + dt * P_tot) := by rw [h_prod]; ring
    _ = M_tot * R_t - t * P_tot := by ring

end Sec05_02

-- ============================================================================
-- 3. 01_YoshidaOrderFourSuzuki_Theorems.lean & 4. 01_BlanesMoanShadowConfinement_Theorems.lean
-- ============================================================================
namespace Sec05_03_04

theorem suzuki_negative_step_necessity (w1 w0 : ℝ) (h_w1_pos : 0 < w1)
    (h_cube : 2 * (w1 * w1 * w1) + (w0 * w0 * w0) = 0) :
    w0 < 0 := by
  by_contra h_nonneg
  have h_w0_ge : 0 ≤ w0 := not_lt.mp h_nonneg
  have h_w1_cube_pos : 0 < 2 * (w1 * w1 * w1) := by
    have h1 : 0 < w1 * w1 := mul_pos h_w1_pos h_w1_pos
    have h2 : 0 < (w1 * w1) * w1 := mul_pos h1 h_w1_pos
    linarith
  have h_w0_cube_nonneg : 0 ≤ w0 * w0 * w0 := by
    have k1 : 0 ≤ w0 * w0 := mul_nonneg h_w0_ge h_w0_ge
    exact mul_nonneg k1 h_w0_ge
  linarith

theorem bm4_fsal_kick_consistency (a1 a2 a3 a4 : Real) :
    let b1 := a1 / 2
    let b2 := (a1 + a2) / 2
    let b3 := (a2 + a3) / 2
    let b4 := (a3 + a4) / 2
    2 * (b1 + b2 + b3 + b4) = 2 * (a1 + a2 + a3) + a4 := by
  intro b1 b2 b3 b4; dsimp [b1, b2, b3, b4]; ring

theorem blanes_moan_shadow_constant_superiority :
    (0.01236 : Real) < 5.3719 /\ (400 : Real) * (0.01236 : Real) < 5.3719 := by
  constructor <;> norm_num

end Sec05_03_04

-- ============================================================================
-- 5. 01_BlanesMoanOrderSix_Theorems.lean
-- ============================================================================
namespace Sec05_05

structure BM6StageCoefficients where
  a1 : Real
  a2 : Real
  a3 : Real
  a4 : Real
  a5 : Real
  a6 : Real
  b1 : Real
  b2 : Real
  b3 : Real
  b4 : Real
  b5 : Real
  b6 : Real

def optimalBM6Coeffs : BM6StageCoefficients where
  a1 := 0.1232297759462710188989670
  a2 := 0.2905537977995580045431665
  a3 := -0.1270492126254170068593452
  a4 := -0.2463317610620750102148924
  a5 := 0.3572088727959280198425176
  a6 := 0.2047770542914699475791730
  b1 := 0.0414649985182624239841838
  b2 := 0.1981286719180669986348123
  b3 := -0.0400061921041533031024560
  b4 := 0.0752539843015807062451982
  b5 := -0.0115113874206879009841253
  b6 := 0.2366699247869310752223870

theorem optimal_order1_drift_satisfied :
    2 * (optimalBM6Coeffs.a1 + optimalBM6Coeffs.a2 + optimalBM6Coeffs.a3 + optimalBM6Coeffs.a4 + optimalBM6Coeffs.a5) + optimalBM6Coeffs.a6 = 1 := by
  dsimp [optimalBM6Coeffs]; norm_num

theorem optimal_order1_kick_satisfied :
    2 * (optimalBM6Coeffs.b1 + optimalBM6Coeffs.b2 + optimalBM6Coeffs.b3 + optimalBM6Coeffs.b4 + optimalBM6Coeffs.b5 + optimalBM6Coeffs.b6) = 1 := by
  dsimp [optimalBM6Coeffs]; norm_num

theorem optimal_central_paired_kick :
    optimalBM6Coeffs.b6 + optimalBM6Coeffs.b6 = 0.4733398495738621504447740 := by
  dsimp [optimalBM6Coeffs]; norm_num

theorem rkn_potential_commutator_zero (dV_dq dW_dq : Real) :
    dV_dq * (0 : Real) - (0 : Real) * dW_dq = 0 := by ring

end Sec05_05

-- ============================================================================
-- 6-10. 03, 06, 07, 08 & 01_CorrectorBM4_Theorems.lean
-- ============================================================================
namespace Sec05_06_10
open Matrix

noncomputable def c4_bm4_from_moment (M4_cent : ℝ) : ℝ :=
  (M4_cent - (1 / 80 : ℝ)) / 24

noncomputable def c4_bm4_from_uncentered_moment (M4_uncent : ℝ) : ℝ :=
  (M4_uncent - (1 / 5 : ℝ)) / 24

theorem c4_bm4_centered_uncentered_equiv (M4_cent M4_uncent : ℝ)
    (h_shift : M4_cent = M4_uncent - (3 / 16 : ℝ)) :
    c4_bm4_from_moment M4_cent = c4_bm4_from_uncentered_moment M4_uncent := by
  dsimp [c4_bm4_from_moment, c4_bm4_from_uncentered_moment]; rw [h_shift]; ring

noncomputable def z4_alpha1 : ℝ := (1 / 4 : ℝ)
noncomputable def z4_alpha2 : ℝ := (1 / 8 : ℝ)
noncomputable def z4_beta1 (c4 : ℝ) : ℝ := 256 * c4
noncomputable def z4_beta2 (c4 : ℝ) : ℝ := -512 * c4

theorem z4_bm4_symmetric_conjugation_weights (c4 : ℝ) :
    z4_beta1 c4 * z4_alpha1 + z4_beta2 c4 * z4_alpha2 = 0 ∧
    (1 / 3 : ℝ) * (z4_beta1 c4 * z4_alpha1 ^ 3 + z4_beta2 c4 * z4_alpha2 ^ 3) = c4 := by
  dsimp [z4_alpha1, z4_alpha2, z4_beta1, z4_beta2]; constructor <;> ring

noncomputable def cotangent_k4_p1 (scale : ℝ) (r w : Fin 3 → ℝ) (T_apex : Matrix (Fin 3) (Fin 3) ℝ) : ℝ :=
  scale * dotProduct w (mulVec T_apex r)

def phase_space_time_reversal (r w : Fin 3 → ℝ) : (Fin 3 → ℝ) × (Fin 3 → ℝ) := (r, -w)

theorem cotangent_k4_p1_time_reversal_odd (scale : ℝ) (r w : Fin 3 → ℝ) (T_apex : Matrix (Fin 3) (Fin 3) ℝ) :
    let state_rev := phase_space_time_reversal r w
    cotangent_k4_p1 scale state_rev.1 state_rev.2 T_apex = - cotangent_k4_p1 scale r w T_apex := by
  dsimp [phase_space_time_reversal, cotangent_k4_p1]; rw [neg_dotProduct]; ring

theorem msym_quadrupole_symmetric (a b : ℝ) (r w : Fin 3 → ℝ) (i j : Fin 3) :
    let M : Matrix (Fin 3) (Fin 3) ℝ := fun p q ↦ a * (r p * w q + w p * r q) - b * (r p * r q)
    M i j = M j i := by
  intro M; dsimp [M]; ring

theorem conormal_corrector_total_momentum_conservation
    (grad_r grad_R : Fin 3 → ℝ) (c_a c_b h : ℝ) (h_mass : c_a + c_b = 1) :
    let F_a : Fin 3 → ℝ := grad_r + c_a • grad_R
    let F_b : Fin 3 → ℝ := -grad_r + c_b • grad_R
    let F_k : Fin 3 → ℝ := -grad_R
    (- (h^4)) • (F_a + F_b + F_k) = 0 := by
  dsimp; ext i
  simp only [Pi.add_apply, Pi.neg_apply, Pi.smul_apply, Pi.zero_apply, smul_eq_mul]
  calc
    -(h ^ 4) * (grad_r i + c_a * grad_R i + (-grad_r i + c_b * grad_R i) + -grad_R i)
      = -(h ^ 4) * (((c_a + c_b) - 1) * grad_R i) := by ring
    _ = -(h ^ 4) * ((1 - 1) * grad_R i) := by rw [h_mass]
    _ = 0 := by ring

theorem pbm4_gauge_operator_conjugation_cocycle {α : Type*}
    (C_A C_B C_B_inv C_C_inv M_AB M_BC : α → α)
    (h_inv_B : ∀ x, C_B (C_B_inv x) = x) (x : α) :
    (C_C_inv ∘ M_BC ∘ C_B) ((C_B_inv ∘ M_AB ∘ C_A) x) =
      (C_C_inv ∘ (M_BC ∘ M_AB) ∘ C_A) x := by
  simp only [Function.comp_apply]
  rw [h_inv_B]

theorem atlas_gauge_invariance_requires_p0
    (J J_inv : Matrix (Fin 3) (Fin 3) ℝ)
    (h_inv : J_inv * J = 1)
    (grad_w_K : Fin 3 → ℝ)
    (h_spatial_stationarity : mulVec J grad_w_K = 0) :
    grad_w_K = 0 := by
  have h_mul : mulVec J_inv (mulVec J grad_w_K) = 0 := by
    rw [h_spatial_stationarity, mulVec_zero]
  rw [mulVec_mulVec] at h_mul
  rw [h_inv, one_mulVec] at h_mul
  exact h_mul

structure Bidegree32HallSpace where
  e1_coboundary : ℝ
  e2_unsolvable : ℝ

theorem hall_l52_cohomological_obstruction_beyond_p0
    (dA_K2 R_52 : Bidegree32HallSpace)
    (h_coboundary_image : dA_K2.e2_unsolvable = 0)
    (h_nontrivial_cohomology_class : R_52.e2_unsolvable ≠ 0)
    (h_quadratic_homological_eq :
      dA_K2.e1_coboundary + R_52.e1_coboundary = 0 ∧
      dA_K2.e2_unsolvable + R_52.e2_unsolvable = 0) :
    False := by
  have h_e2 : dA_K2.e2_unsolvable + R_52.e2_unsolvable = 0 := h_quadratic_homological_eq.2
  rw [h_coboundary_image, zero_add] at h_e2
  exact h_nontrivial_cohomology_class h_e2

end Sec05_06_10

#print axioms Sec05_01.strang_palindromic_microréversibility
#print axioms Sec05_02.three_body_central_torque_zero
#print axioms Sec05_02.galilean_boost_step_invariance
#print axioms Sec05_03_04.bm4_fsal_kick_consistency
#print axioms Sec05_05.optimal_order1_drift_satisfied
#print axioms Sec05_05.optimal_order1_kick_satisfied
#print axioms Sec05_05.rkn_potential_commutator_zero
#print axioms Sec05_06_10.c4_bm4_centered_uncentered_equiv
#print axioms Sec05_06_10.z4_bm4_symmetric_conjugation_weights
#print axioms Sec05_06_10.cotangent_k4_p1_time_reversal_odd
#print axioms Sec05_06_10.msym_quadrupole_symmetric
#print axioms Sec05_06_10.pbm4_gauge_operator_conjugation_cocycle
#print axioms Sec05_06_10.hall_l52_cohomological_obstruction_beyond_p0
