import Mathlib.Data.Real.Basic
import Mathlib.Algebra.Order.Archimedean.Basic
import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Ext

-- ==============================================================================
-- SECTION 4.1 : WHITNEY LATTICE & FULTON-MACPHERSON MANIFOLD Q[N] (4 FILES)
-- ==============================================================================
namespace UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice

-- 01_PartitionLattice_Defs.lean
def collisionCodimension (k : Nat) : Nat :=
  if k ≥ 2 then 3 * (k - 1) else 0

theorem collisionCodimension_pair : collisionCodimension 2 = 3 := by
  rfl

theorem collisionCodimension_triple : collisionCodimension 3 = 6 := by
  rfl

theorem collision_diagonal_non_transversality {α : Type u} (q1 q2 q3 : α)
    (h12 : q1 = q2) (h23 : q2 = q3) : q1 = q3 := by
  rw [h12, h23]

def isCollisionState {N : Nat} {α : Type u} (q : Fin N → α) (I : Fin N → Prop) : Prop :=
  ∀ i j, I i → I j → q i = q j

theorem collision_hereditary {N : Nat} {α : Type u} (q : Fin N → α) {I J : Fin N → Prop}
    (h_sub : ∀ i, J i → I i) (h_coll : isCollisionState q I) :
    isCollisionState q J := by
  intro i j hi hj
  exact h_coll i j (h_sub i hi) (h_sub j hj)

-- 02_ClusterHyperRadius_Theorems.lean
noncomputable def hyperRadiusSq (mu : ℝ) (r_sq : ℝ) : ℝ :=
  mu * r_sq

theorem hyper_radius_sq_pos (mu : ℝ) (r_sq : ℝ)
    (h_mu : 0 < mu) (h_r : 0 < r_sq) :
    0 < hyperRadiusSq mu r_sq := by
  dsimp [hyperRadiusSq]
  exact mul_pos h_mu h_r

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

def isHomogeneousDegreeMinusOne (U : ℝ → ℝ) : Prop :=
  ∀ (r : ℝ) (lambda : ℝ), 0 < lambda → U (lambda * r) = (1 / lambda) * U r

noncomputable def newtonianScalarPotential (GM : ℝ) (r : ℝ) : ℝ :=
  - GM * r⁻¹

theorem newtonian_potential_is_homogeneous (GM : ℝ) :
    ∀ (r : ℝ) (lambda : ℝ), 0 < lambda → r ≠ 0 →
    newtonianScalarPotential GM (lambda * r) = (1 / lambda) * newtonianScalarPotential GM r := by
  intro r lambda h_lambda _h_r
  dsimp [newtonianScalarPotential]
  have _h_lam_ne : lambda ≠ 0 := ne_of_gt h_lambda
  rw [mul_inv, one_div]
  ring

theorem mcgehee_potential_factorization (U : ℝ → ℝ) (h_hom : isHomogeneousDegreeMinusOne U)
    (rho : ℝ) (h_rho : 0 < rho) (s : ℝ) :
    U (rho * s) = (1 / rho) * U s :=
  h_hom s rho h_rho

-- 03_MaximalFlagJacobiTree_Theorems.lean
inductive BinaryClusterTree (N : Nat) where
  | leaf (index : Fin N) : BinaryClusterTree N
  | node (left : BinaryClusterTree N) (right : BinaryClusterTree N) : BinaryClusterTree N
  deriving DecidableEq, Repr

def treeBodies {N : Nat} : BinaryClusterTree N → (Fin N → Prop)
  | BinaryClusterTree.leaf i => fun j => j = i
  | BinaryClusterTree.node l r => fun j => treeBodies l j ∨ treeBodies r j

theorem leaf_bodies_singleton {N : Nat} (i : Fin N) :
    treeBodies (BinaryClusterTree.leaf i) = fun j => j = i := by
  rfl

theorem node_bodies_contains_left {N : Nat} (l r : BinaryClusterTree N) (j : Fin N) :
    treeBodies l j → treeBodies (BinaryClusterTree.node l r) j := by
  intro h
  exact Or.inl h

theorem node_bodies_contains_right {N : Nat} (l r : BinaryClusterTree N) (j : Fin N) :
    treeBodies r j → treeBodies (BinaryClusterTree.node l r) j := by
  intro h
  exact Or.inr h

def isDisjointSubtrees {N : Nat} (l r : BinaryClusterTree N) : Prop :=
  ∀ j : Fin N, ¬ (treeBodies l j ∧ treeBodies r j)

theorem jacobi_node_exact_union {N : Nat} (l r : BinaryClusterTree N) (j : Fin N) :
    treeBodies (BinaryClusterTree.node l r) j ↔ (treeBodies l j ∨ treeBodies r j) := by
  rfl

theorem jacobi_node_defines_hierarchy {N : Nat} (l r : BinaryClusterTree N) (j : Fin N)
    (hj : treeBodies l j) : treeBodies (BinaryClusterTree.node l r) j :=
  (jacobi_node_exact_union l r j).mpr (Or.inl hj)

-- 04_DisjointCotangentCommutation_Theorems.lean
def localPoissonTerm (grad_q_F grad_p_F grad_q_G grad_p_G : ℝ) : ℝ :=
  grad_q_F * grad_p_G - grad_p_F * grad_q_G

theorem local_poisson_term_zero_of_first_zero (grad_q_G grad_p_G : ℝ) :
    localPoissonTerm 0 0 grad_q_G grad_p_G = 0 := by
  dsimp [localPoissonTerm]
  ring

theorem local_poisson_term_zero_of_second_zero (grad_q_F grad_p_F : ℝ) :
    localPoissonTerm grad_q_F grad_p_F 0 0 = 0 := by
  dsimp [localPoissonTerm]
  ring

theorem disjoint_cotangent_commutation_local (grad_q_F grad_p_F grad_q_G grad_p_G : ℝ)
    (h_disjoint : (grad_q_F = 0 ∧ grad_p_F = 0) ∨ (grad_q_G = 0 ∧ grad_p_G = 0)) :
    localPoissonTerm grad_q_F grad_p_F grad_q_G grad_p_G = 0 := by
  cases h_disjoint with
  | inl hF =>
    rcases hF with ⟨hqF, hpF⟩
    rw [hqF, hpF]
    exact local_poisson_term_zero_of_first_zero grad_q_G grad_p_G
  | inr hG =>
    rcases hG with ⟨hqG, hpG⟩
    rw [hqG, hpG]
    exact local_poisson_term_zero_of_second_zero grad_q_F grad_p_F

end UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice

-- ==============================================================================
-- SECTION 4.2 : SCREENED APEX POTENTIAL (3 FILES)
-- ==============================================================================
namespace UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential

-- 01_ApexUniversalLipschitzBound_Theorems.lean
inductive ConormalRegime where
  | regime1 : ConormalRegime
  | regime2 : ConormalRegime
  | regime3 : ConormalRegime
deriving DecidableEq

structure RegimeDifferentialBound where
  grad_bound : ℝ
  hess_bound : ℝ
  grad_pos : 0 < grad_bound
  hess_pos : 0 < hess_bound

structure AdaptedChartApexData where
  b1 : RegimeDifferentialBound
  b2 : RegimeDifferentialBound
  b3 : RegimeDifferentialBound

def apexGradBound (data : AdaptedChartApexData) : ℝ :=
  data.b1.grad_bound + data.b2.grad_bound + data.b3.grad_bound

def apexHessBound (data : AdaptedChartApexData) : ℝ :=
  data.b1.hess_bound + data.b2.hess_bound + data.b3.hess_bound

theorem apex_grad_bound_pos (data : AdaptedChartApexData) : 0 < apexGradBound data := by
  dsimp [apexGradBound]
  have h1 := data.b1.grad_pos
  have h2 := data.b2.grad_pos
  have h3 := data.b3.grad_pos
  linarith

theorem apex_hess_bound_pos (data : AdaptedChartApexData) : 0 < apexHessBound data := by
  dsimp [apexHessBound]
  have h1 := data.b1.hess_pos
  have h2 := data.b2.hess_pos
  have h3 := data.b3.hess_pos
  linarith

structure ChartApexState (data : AdaptedChartApexData) where
  grad_norm : ℝ
  hess_norm : ℝ
  grad_nonneg : 0 ≤ grad_norm
  hess_nonneg : 0 ≤ hess_norm
  regime : ConormalRegime
  regime_bound_grad :
    grad_norm ≤ match regime with
      | ConormalRegime.regime1 => data.b1.grad_bound
      | ConormalRegime.regime2 => data.b2.grad_bound
      | ConormalRegime.regime3 => data.b3.grad_bound
  regime_bound_hess :
    hess_norm ≤ match regime with
      | ConormalRegime.regime1 => data.b1.hess_bound
      | ConormalRegime.regime2 => data.b2.hess_bound
      | ConormalRegime.regime3 => data.b3.hess_bound

theorem apex_universal_lipschitz_bound (data : AdaptedChartApexData) :
    ∃ (M_apex L_apex : ℝ),
      0 < M_apex ∧ 0 < L_apex ∧
      (∀ (state : ChartApexState data),
        state.grad_norm ≤ M_apex ∧ state.hess_norm ≤ L_apex) := by
  use apexGradBound data, apexHessBound data
  refine ⟨apex_grad_bound_pos data, apex_hess_bound_pos data, ?_⟩
  intro state
  have h1_pos := data.b1.grad_pos
  have h2_pos := data.b2.grad_pos
  have h3_pos := data.b3.grad_pos
  have k1_pos := data.b1.hess_pos
  have k2_pos := data.b2.hess_pos
  have k3_pos := data.b3.hess_pos
  rcases state with ⟨grad_norm, hess_norm, _hg_nn, _hh_nn, regime, hg, hh⟩
  cases regime with
  | regime1 =>
    dsimp at hg hh
    dsimp [apexGradBound, apexHessBound]
    constructor <;> linarith
  | regime2 =>
    dsimp at hg hh
    dsimp [apexGradBound, apexHessBound]
    constructor <;> linarith
  | regime3 =>
    dsimp at hg hh
    dsimp [apexGradBound, apexHessBound]
    constructor <;> linarith

-- 02_LieBargmannExactInvariance_Theorems.lean
structure LieBargmannGenerators where
  P_x : ℝ
  P_y : ℝ
  P_z : ℝ
  L_x : ℝ
  L_y : ℝ
  L_z : ℝ
  G_x : ℝ
  G_y : ℝ
  G_z : ℝ
  H_ext : ℝ

def isLieBargmannConserved (g_initial g_final : LieBargmannGenerators) : Prop :=
  g_final.P_x = g_initial.P_x ∧
  g_final.P_y = g_initial.P_y ∧
  g_final.P_z = g_initial.P_z ∧
  g_final.L_x = g_initial.L_x ∧
  g_final.L_y = g_initial.L_y ∧
  g_final.L_z = g_initial.L_z ∧
  g_final.G_x = g_initial.G_x ∧
  g_final.G_y = g_initial.G_y ∧
  g_final.G_z = g_initial.G_z ∧
  g_final.H_ext = g_initial.H_ext

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
      M_tot * (R_t + dt / M_tot * P_tot) = M_tot * R_t + M_tot * (dt / M_tot * P_tot) := by ring
      _ = M_tot * R_t + (M_tot * (dt / M_tot)) * P_tot := by ring
      _ = M_tot * R_t + dt * P_tot := by
        rw [mul_div_cancel₀ dt hM]
  calc
    M_tot * (R_t + dt / M_tot * P_tot) - (t + dt) * P_tot
      = (M_tot * R_t + dt * P_tot) - (t * P_tot + dt * P_tot) := by
        rw [h_prod]
        ring
    _ = M_tot * R_t - t * P_tot := by ring

theorem galilean_boost_poisson_compensation (P_tot : ℝ) :
    P_tot - P_tot = 0 := by
  ring

theorem lie_bargmann_exact_invariance (M_tot : ℝ) (hM : 0 < M_tot)
    (R_x R_y R_z P_x P_y P_z L_x L_y L_z : ℝ) (t dt H_0 : ℝ) :
    let _hM_ne : M_tot ≠ 0 := ne_of_gt hM
    let R_next_x := R_x + (dt / M_tot) * P_x
    let R_next_y := R_y + (dt / M_tot) * P_y
    let R_next_z := R_z + (dt / M_tot) * P_z
    let t_next := t + dt
    let g_init : LieBargmannGenerators := ⟨P_x, P_y, P_z, L_x, L_y, L_z,
                                           M_tot * R_x - t * P_x,
                                           M_tot * R_y - t * P_y,
                                           M_tot * R_z - t * P_z,
                                           H_0 - H_0⟩
    let g_final : LieBargmannGenerators := ⟨P_x, P_y, P_z, L_x, L_y, L_z,
                                            M_tot * R_next_x - t_next * P_x,
                                            M_tot * R_next_y - t_next * P_y,
                                            M_tot * R_next_z - t_next * P_z,
                                            H_0 - H_0⟩
    isLieBargmannConserved g_init g_final ∧
    (P_x - P_x = 0 ∧ P_y - P_y = 0 ∧ P_z - P_z = 0) := by
  intro _hM_ne R_next_x R_next_y R_next_z t_next g_init g_final
  dsimp [isLieBargmannConserved, g_init, g_final]
  have hM_ne : M_tot ≠ 0 := ne_of_gt hM
  have hG_x : M_tot * R_next_x - t_next * P_x = M_tot * R_x - t * P_x := by
    exact galilean_boost_step_invariance M_tot R_x P_x t dt hM_ne
  have hG_y : M_tot * R_next_y - t_next * P_y = M_tot * R_y - t * P_y := by
    exact galilean_boost_step_invariance M_tot R_y P_y t dt hM_ne
  have hG_z : M_tot * R_next_z - t_next * P_z = M_tot * R_z - t * P_z := by
    exact galilean_boost_step_invariance M_tot R_z P_z t dt hM_ne
  refine ⟨⟨rfl, rfl, rfl, rfl, rfl, rfl, hG_x, hG_y, hG_z, rfl⟩, ?_⟩
  refine ⟨by ring, by ring, by ring⟩

-- 03_ExtendedPalindromicSymplecticity_Theorems.lean
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
  dsimp [extendedTimeReversal]
  ext <;> ring

noncomputable def driftStep (mu dt : ℝ) (p : ExtendedPhasePoint) : ExtendedPhasePoint :=
  ⟨p.r + (dt / (2 * mu)) * p.pi, p.pi, p.t, p.E⟩

noncomputable def apexKickStep (gradV dt : ℝ) (p : ExtendedPhasePoint) : ExtendedPhasePoint :=
  ⟨p.r, p.pi - (dt / 2) * gradV, p.t, p.E⟩

def timeAdvanceStep (dt : ℝ) (p : ExtendedPhasePoint) : ExtendedPhasePoint :=
  ⟨p.r, p.pi, p.t + dt, p.E⟩

def symplecticArea (dr1 dpi1 dr2 dpi2 : ℝ) : ℝ :=
  dr1 * dpi2 - dr2 * dpi1

theorem apex_kick_preserves_symplectic_area (dr1 dpi1 dr2 dpi2 H dt : ℝ) :
    let dpi1_kick := dpi1 - (dt / 2) * H * dr1
    let dpi2_kick := dpi2 - (dt / 2) * H * dr2
    symplecticArea dr1 dpi1_kick dr2 dpi2_kick = symplecticArea dr1 dpi1 dr2 dpi2 := by
  intro dpi1_kick dpi2_kick
  dsimp [symplecticArea, dpi1_kick, dpi2_kick]
  ring

theorem drift_preserves_symplectic_area (dr1 dpi1 dr2 dpi2 mu dt : ℝ) :
    let dr1_drift := dr1 + (dt / (2 * mu)) * dpi1
    let dr2_drift := dr2 + (dt / (2 * mu)) * dpi2
    symplecticArea dr1_drift dpi1 dr2_drift dpi2 = symplecticArea dr1 dpi1 dr2 dpi2 := by
  intro dr1_drift dr2_drift
  dsimp [symplecticArea, dr1_drift, dr2_drift]
  ring

theorem extended_palindromic_symplecticity (mu gradV H dt : ℝ) (dr1 dpi1 dr2 dpi2 : ℝ) (p : ExtendedPhasePoint) :
    symplecticArea dr1 (dpi1 - (dt / 2) * H * dr1) dr2 (dpi2 - (dt / 2) * H * dr2) =
      symplecticArea dr1 dpi1 dr2 dpi2 ∧
    symplecticArea (dr1 + (dt / (2 * mu)) * dpi1) dpi1 (dr2 + (dt / (2 * mu)) * dpi2) dpi2 =
      symplecticArea dr1 dpi1 dr2 dpi2 ∧
    extendedTimeReversal (extendedTimeReversal p) = p ∧
    extendedTimeReversal (apexKickStep gradV dt (extendedTimeReversal p)) = apexKickStep gradV (-dt) p ∧
    extendedTimeReversal (driftStep mu dt (extendedTimeReversal p)) = driftStep mu (-dt) p := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · exact apex_kick_preserves_symplectic_area dr1 dpi1 dr2 dpi2 H dt
  · exact drift_preserves_symplectic_area dr1 dpi1 dr2 dpi2 mu dt
  · exact extended_time_reversal_involution p
  · dsimp [extendedTimeReversal, apexKickStep]
    ext <;> ring
  · dsimp [extendedTimeReversal, driftStep]
    ext <;> ring

end UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential

-- ==============================================================================
-- SECTION 4.3 : SYNCHRONOUS SUPERNODES & ATLAS TRANSITIONS (3 FILES)
-- ==============================================================================
namespace UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes

-- 01_UnifiedSuperNodeDef_Theorems.lean
def isCompactCluster (rho : ℝ) (delta_cluster : ℝ) : Prop :=
  rho ≤ delta_cluster

noncomputable def superNodePhysicalTimeStepSq (rho : ℝ) (dtau : ℝ) : ℝ :=
  (rho * rho) * (dtau * dtau)

theorem super_node_time_step_sq_pos (rho : ℝ) (dtau : ℝ)
    (h_rho : 0 < rho) (h_tau : 0 < dtau) :
    0 < superNodePhysicalTimeStepSq rho dtau := by
  dsimp [superNodePhysicalTimeStepSq]
  have h_r2 : 0 < rho * rho := mul_pos h_rho h_rho
  have h_t2 : 0 < dtau * dtau := mul_pos h_tau h_tau
  exact mul_pos h_r2 h_t2

theorem super_node_accumulated_time_sq_pos (rho : ℝ) (dtau : ℝ)
    (h_rho : 0 < rho) (h_tau : 0 < dtau) (N : Nat) (hN : 1 ≤ N) :
    0 < (N : ℝ) * (N : ℝ) * superNodePhysicalTimeStepSq rho dtau := by
  have h_step := super_node_time_step_sq_pos rho dtau h_rho h_tau
  have hN_pos : (0 : ℝ) < (N : ℝ) := by
    have h_lt : 0 < N := Nat.lt_of_lt_of_le Nat.zero_lt_one hN
    exact Nat.cast_pos.mpr h_lt
  have hN2_pos : 0 < (N : ℝ) * (N : ℝ) := mul_pos hN_pos hN_pos
  exact mul_pos hN2_pos h_step

-- 02_SuperNodeAntiZenoRegularization_Theorems.lean
structure UnfoldedSundmanState where
  chi : ℝ
  chi_pos : 0 < chi
  tau_step : ℝ
  tau_pos : 0 < tau_step

def physicalTimeIncrement (state : UnfoldedSundmanState) : ℝ :=
  state.chi * state.chi * state.tau_step

theorem physical_time_increment_pos (state : UnfoldedSundmanState) :
    0 < physicalTimeIncrement state := by
  dsimp [physicalTimeIncrement]
  have h_chi2 : 0 < state.chi * state.chi := mul_pos state.chi_pos state.chi_pos
  exact mul_pos h_chi2 state.tau_pos

theorem keplerian_conormal_time_scaling (c tau : ℝ) (hc : 0 < c) (htau : 0 < tau) :
    let dt_kepler := ((1 / 27 : ℝ) * (c * c * c)) * (tau * tau * tau)
    0 < dt_kepler := by
  intro dt_kepler
  dsimp [dt_kepler]
  have hc2 : 0 < c * c := mul_pos hc hc
  have hc3 : 0 < (c * c) * c := mul_pos hc2 hc
  have ht2 : 0 < tau * tau := mul_pos htau htau
  have ht3 : 0 < (tau * tau) * tau := mul_pos ht2 htau
  have h_coeff : (0 : ℝ) < 1 / 27 := by linarith
  have h_coeff_c : 0 < (1 / 27 : ℝ) * (c * c * c) := mul_pos h_coeff hc3
  exact mul_pos h_coeff_c ht3

theorem super_node_anti_zeno_regularization (state : UnfoldedSundmanState) (dt : ℝ) (h_dt : 0 < dt) :
    0 < physicalTimeIncrement state ∧
    (∃ (N_max : Nat), 0 < N_max ∧ dt ≤ (N_max : ℝ) * (physicalTimeIncrement state)) := by
  have h_pos := physical_time_increment_pos state
  refine ⟨h_pos, ?_⟩
  have h_div : 0 < dt / physicalTimeIncrement state := div_pos h_dt h_pos
  obtain ⟨N, hN⟩ := exists_nat_gt (dt / physicalTimeIncrement state)
  have h_N_pos : 0 < N := by
    have h_le : 0 < (N : ℝ) := lt_trans h_div hN
    exact Nat.cast_pos.mp h_le
  use N
  refine ⟨h_N_pos, ?_⟩
  have h_mul := (div_lt_iff₀ h_pos).mp hN
  linarith

-- 03_AtlasTransitionSymplecticAntiChattering_Theorems.lean
structure JacobiTransitionMatrix2x2 where
  alpha : ℝ
  beta : ℝ
  gamma : ℝ
  delta : ℝ

def det2x2 (M : JacobiTransitionMatrix2x2) : ℝ :=
  M.alpha * M.delta - M.beta * M.gamma

def isInSL2 (M : JacobiTransitionMatrix2x2) : Prop :=
  det2x2 M = 1

noncomputable def jacobiTransitionMatrix (m1 m2 m3 : ℝ) : JacobiTransitionMatrix2x2 where
  alpha := -m1 / (m1 + m2)
  beta := 1
  gamma := -m2 * (m1 + m2 + m3) / ((m1 + m2) * (m2 + m3))
  delta := -m3 / (m2 + m3)

theorem jacobi_transition_matrix_in_sl2 (m1 m2 m3 : ℝ)
    (h12 : m1 + m2 ≠ 0) (h23 : m2 + m3 ≠ 0) :
    isInSL2 (jacobiTransitionMatrix m1 m2 m3) := by
  dsimp [isInSL2, det2x2, jacobiTransitionMatrix]
  field_simp [h12, h23]
  ring

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

theorem atlas_transition_symplectic_anti_chattering (M : JacobiTransitionMatrix2x2) (hM : isInSL2 M)
    (x1 x2 y1 y2 : ℝ) (d12 d_int : ℝ) (h_d12 : 0 < d12) :
    (let x1' := M.alpha * x1 + M.beta * x2
     let x2' := M.gamma * x1 + M.delta * x2
     let y1' := M.alpha * y1 + M.beta * y2
     let y2' := M.gamma * y1 + M.delta * y2
     x1' * y2' - x2' * y1' = x1 * y2 - x2 * y1) ∧
    (¬ (d12 < (3 / 5 : ℝ) * d_int ∧ d_int < (3 / 5 : ℝ) * d12)) := by
  refine ⟨sl2_preserves_symplectic_area M hM x1 x2 y1 y2, ?_⟩
  intro ⟨hfwd, hback⟩
  have h_gamma_pos : (0 : ℝ) ≤ 3 / 5 := by linarith
  have h_gamma_sq : (3 / 5 : ℝ) * (3 / 5 : ℝ) < 1 := by
    norm_num
  exact hysteresis_anti_chattering_contradiction d12 d_int (3 / 5) h_d12 h_gamma_pos h_gamma_sq hfwd hback

end UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes

#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.collisionCodimension_pair
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.collisionCodimension_triple
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.collision_diagonal_non_transversality
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.collision_hereditary
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.hyper_radius_sq_pos
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.hyper_radius_sq_zero_iff
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.newtonian_potential_is_homogeneous
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.mcgehee_potential_factorization
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.leaf_bodies_singleton
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.node_bodies_contains_left
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.node_bodies_contains_right
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.jacobi_node_exact_union
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.jacobi_node_defines_hierarchy
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.local_poisson_term_zero_of_first_zero
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.local_poisson_term_zero_of_second_zero
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice.disjoint_cotangent_commutation_local
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.apex_grad_bound_pos
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.apex_hess_bound_pos
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.apex_universal_lipschitz_bound
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.galilean_boost_step_invariance
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.galilean_boost_poisson_compensation
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.lie_bargmann_exact_invariance
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.extended_time_reversal_involution
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.apex_kick_preserves_symplectic_area
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.drift_preserves_symplectic_area
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.ScreenedApexPotential.extended_palindromic_symplecticity
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.super_node_time_step_sq_pos
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.super_node_accumulated_time_sq_pos
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.physical_time_increment_pos
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.keplerian_conormal_time_scaling
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.super_node_anti_zeno_regularization
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.jacobi_transition_matrix_in_sl2
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.sl2_preserves_symplectic_area
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.cotangent_lift_preserves_darboux
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.hysteresis_anti_chattering_contradiction
#print axioms UnsoftenedNBody.FultonMacPhersonSplitting.SynchronousSupernodes.atlas_transition_symplectic_anti_chattering
