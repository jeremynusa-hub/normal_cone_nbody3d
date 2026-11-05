-- Dependencies: 01_ConfigurationSpace_Defs, 04_SymplecticForm_Defs, 01_ReesDeformation_Defs
import Mathlib.Tactic.Ring

universe u

variable {k : Type u} [CommRing k]

/-- Fundamental Theorem of Calaque-Safronov Degeneration:
    At the central apex t = 0, the physical evaluation morphism retracts
    strictly and identically onto the singular base stratum Z (q_parallel). -/
theorem deformation_limit_zero (d : NormalConePoint k) :
    deformationPhysicalPosition d 0 = d.q_parallel := by
  dsimp [deformationPhysicalPosition, Vec3.smul, Vec3.add]
  cases d with
  | mk qp n =>
    dsimp
    ext <;> ring

/-- Skew-symmetry of the canonical normal cone symplectic 2-form ω_normal on T*(N_Z M). -/
theorem omegaNormal_skew (v1 v2 : ConormalTangentVector k) :
    omegaNormal v1 v2 = - omegaNormal v2 v1 := by
  dsimp [omegaNormal]
  ring

/-- Alternating property of ω_normal. -/
theorem omegaNormal_self (v : ConormalTangentVector k) :
    omegaNormal v v = 0 := by
  dsimp [omegaNormal]
  ring

/-- Exact Symplectomorphism Theorem of the Transverse Normal Fiber Flow:
    Under potential Hessian symmetry (dot dF1 v2.dn = dot dF2 v1.dn),
    the normal fiber kick Φ_s along the flat abelian fibers of N_Z* M
    strictly preserves the symplectic form ω_normal. -/
theorem transverse_conormal_flow_symplectic
    (s : k) (dF1 dF2 : Vec3 k) (v1 v2 : ConormalTangentVector k)
    (h_hess_symm : Vec3.dot dF1 v2.dn = Vec3.dot dF2 v1.dn) :
    omegaNormal (transverseConormalTangentFlow s dF1 v1)
                (transverseConormalTangentFlow s dF2 v2) =
    omegaNormal v1 v2 := by
  dsimp [transverseConormalTangentFlow, omegaNormal, Vec3.dot, Vec3.add, Vec3.smul]
  calc (v1.dq.x * v2.dp.x + v1.dq.y * v2.dp.y + v1.dq.z * v2.dp.z -
        (v2.dq.x * v1.dp.x + v2.dq.y * v1.dp.y + v2.dq.z * v1.dp.z)) +
       (v1.dn.x * (v2.dw.x + s * dF2.x) + v1.dn.y * (v2.dw.y + s * dF2.y) + v1.dn.z * (v2.dw.z + s * dF2.z) -
        (v2.dn.x * (v1.dw.x + s * dF1.x) + v2.dn.y * (v1.dw.y + s * dF1.y) + v2.dn.z * (v1.dw.z + s * dF1.z)))
    = omegaNormal v1 v2 + s * (Vec3.dot v1.dn dF2 - Vec3.dot v2.dn dF1) := by
        dsimp [omegaNormal, Vec3.dot]; ring
    _ = omegaNormal v1 v2 + s * (Vec3.dot dF1 v2.dn - Vec3.dot v2.dn dF1) := by
        have h_comm1 : Vec3.dot v1.dn dF2 = Vec3.dot dF2 v1.dn := by
          dsimp [Vec3.dot]; ring
        have h_symm : Vec3.dot dF2 v1.dn = Vec3.dot dF1 v2.dn := h_hess_symm.symm
        rw [h_comm1, h_symm]
    _ = omegaNormal v1 v2 := by
        have h_comm2 : Vec3.dot dF1 v2.dn = Vec3.dot v2.dn dF1 := by
          dsimp [Vec3.dot]; ring
        rw [h_comm2]
        ring

/-- Constructive Non-Degeneracy of the Normal Cone Symplectic 2-Form ω_normal:
    If ω_normal(v1, v2) = 0 for all variations v2, then every component of v1
    vanishes identically, certifying that ω_normal is a genuine
    non-degenerate symplectic form on T*(N_Z M). -/
theorem omegaNormal_nondegenerate
    (v1 : ConormalTangentVector k)
    (h : ∀ v2 : ConormalTangentVector k, omegaNormal v1 v2 = 0) :
    v1.dq = Vec3.zero ∧ v1.dp = Vec3.zero ∧ v1.dn = Vec3.zero ∧ v1.dw = Vec3.zero := by
  have hq_x := h ⟨Vec3.zero, ⟨1, 0, 0⟩, Vec3.zero, Vec3.zero⟩
  have hq_y := h ⟨Vec3.zero, ⟨0, 1, 0⟩, Vec3.zero, Vec3.zero⟩
  have hq_z := h ⟨Vec3.zero, ⟨0, 0, 1⟩, Vec3.zero, Vec3.zero⟩
  have hp_x := h ⟨⟨1, 0, 0⟩, Vec3.zero, Vec3.zero, Vec3.zero⟩
  have hp_y := h ⟨⟨0, 1, 0⟩, Vec3.zero, Vec3.zero, Vec3.zero⟩
  have hp_z := h ⟨⟨0, 0, 1⟩, Vec3.zero, Vec3.zero, Vec3.zero⟩
  have hn_x := h ⟨Vec3.zero, Vec3.zero, Vec3.zero, ⟨1, 0, 0⟩⟩
  have hn_y := h ⟨Vec3.zero, Vec3.zero, Vec3.zero, ⟨0, 1, 0⟩⟩
  have hn_z := h ⟨Vec3.zero, Vec3.zero, Vec3.zero, ⟨0, 0, 1⟩⟩
  have hw_x := h ⟨Vec3.zero, Vec3.zero, ⟨1, 0, 0⟩, Vec3.zero⟩
  have hw_y := h ⟨Vec3.zero, Vec3.zero, ⟨0, 1, 0⟩, Vec3.zero⟩
  have hw_z := h ⟨Vec3.zero, Vec3.zero, ⟨0, 0, 1⟩, Vec3.zero⟩
  dsimp [omegaNormal, Vec3.dot, Vec3.zero] at hq_x hq_y hq_z hp_x hp_y hp_z hn_x hn_y hn_z hw_x hw_y hw_z
  have eqx : v1.dq.x = 0 := by
    calc v1.dq.x = (v1.dq.x * 1 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hq_x
  have eqy : v1.dq.y = 0 := by
    calc v1.dq.y = (v1.dq.x * 0 + v1.dq.y * 1 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hq_y
  have eqz : v1.dq.z = 0 := by
    calc v1.dq.z = (v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 1 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hq_z
  have epx : v1.dp.x = 0 := by
    calc v1.dp.x = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (1 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hp_x]
    _ = 0 := by ring
  have epy : v1.dp.y = 0 := by
    calc v1.dp.y = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 1 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hp_y]
    _ = 0 := by ring
  have epz : v1.dp.z = 0 := by
    calc v1.dp.z = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 1 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hp_z]
    _ = 0 := by ring
  have enx : v1.dn.x = 0 := by
    calc v1.dn.x = (v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 1 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hn_x
  have eny : v1.dn.y = 0 := by
    calc v1.dn.y = (v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 1 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hn_y
  have enz : v1.dn.z = 0 := by
    calc v1.dn.z = (v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                   (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 1 - (0 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z)) := by ring
    _ = 0 := hn_z
  have ewx : v1.dw.x = 0 := by
    calc v1.dw.x = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (1 * v1.dw.x + 0 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hw_x]
    _ = 0 := by ring
  have ewy : v1.dw.y = 0 := by
    calc v1.dw.y = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 1 * v1.dw.y + 0 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hw_y]
    _ = 0 := by ring
  have ewz : v1.dw.z = 0 := by
    calc v1.dw.z = - ((v1.dq.x * 0 + v1.dq.y * 0 + v1.dq.z * 0 - (0 * v1.dp.x + 0 * v1.dp.y + 0 * v1.dp.z)) +
                      (v1.dn.x * 0 + v1.dn.y * 0 + v1.dn.z * 0 - (0 * v1.dw.x + 0 * v1.dw.y + 1 * v1.dw.z))) := by ring
    _ = - 0 := by rw [hw_z]
    _ = 0 := by ring
  rcases v1 with ⟨dq, dp, dn, dw⟩
  rcases dq with ⟨qx, qy, qz⟩
  rcases dp with ⟨px, py, pz⟩
  rcases dn with ⟨nx, ny, nz⟩
  rcases dw with ⟨wx, wy, wz⟩
  dsimp at eqx eqy eqz epx epy epz enx eny enz ewx ewy ewz
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [eqx, eqy, eqz]; rfl
  · rw [epx, epy, epz]; rfl
  · rw [enx, eny, enz]; rfl
  · rw [ewx, ewy, ewz]; rfl

