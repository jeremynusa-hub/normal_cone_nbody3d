-- Dependencies: 01_ConfigurationSpace_Defs, 04_SymplecticForm_Defs
import Mathlib.Tactic.Ring

universe u

/-- Point in the normal cone / deformation family bundle.
    `q_parallel` denotes the base position along the singular stratum Z.
    `n` denotes the transverse direction in the normal fiber. -/
structure NormalConePoint (k : Type u) where
  q_parallel : Vec3 k
  n          : Vec3 k

/-- Spatial evaluation morphism of the Calaque-Safronov Rees deformation family:
    ev_t(q_parallel, n) = q_parallel + t · n.
    - For t ≠ 0, reconstitutes the generic ambient configuration in M.
    - For t = 0, retracts identically onto the singular base stratum Z. -/
def deformationPhysicalPosition {k : Type u} [CommRing k]
    (d : NormalConePoint k) (t : k) : Vec3 k :=
  Vec3.add d.q_parallel (Vec3.smul t d.n)

/-- Tangent variation vector on the cotangent bundle of the normal cone T*(N_Z M).
    Comprises base variations (dq, dp) and transverse fiber variations (dn, dw). -/
structure ConormalTangentVector (k : Type u) where
  dq : Vec3 k
  dp : Vec3 k
  dn : Vec3 k
  dw : Vec3 k

/-- Canonical Darboux symplectic 2-form on T*(N_Z M):
    ω_normal(v1, v2) = (dq1 · dp2 - dq2 · dp1) + (dn1 · dw2 - dn2 · dw1).
    Inherits an abelian linear foliation along the normal fibers. -/
def omegaNormal {k : Type u} [CommRing k]
    (v1 v2 : ConormalTangentVector k) : k :=
  (Vec3.dot v1.dq v2.dp - Vec3.dot v2.dq v1.dp) +
  (Vec3.dot v1.dn v2.dw - Vec3.dot v2.dn v1.dw)

/-- Transverse normal fiber flow kick by potential force gradient dF along the fiber:
    Leaves base variations (dq, dp, dn) invariant and translates fiber momentum dw. -/
def transverseConormalTangentFlow {k : Type u} [CommRing k]
    (s : k) (dF : Vec3 k) (v : ConormalTangentVector k) : ConormalTangentVector k :=
  { dq := v.dq
    dp := v.dp
    dn := v.dn
    dw := Vec3.add v.dw (Vec3.smul s dF) }
