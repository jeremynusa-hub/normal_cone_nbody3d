-- Dependencies: 01_ConfigurationSpace_Defs (Vec3)
import Mathlib.Tactic.Ring
universe u

structure TangentVector3D (N : Nat) (k : Type u) where
  dq : Fin N → Vec3 k
  dp : Fin N → Vec3 k

variable {N : Nat} {k : Type u} [CommRing k]

def omega0At (v1 v2 : TangentVector3D N k) (i : Fin N) : k :=
  Vec3.dot (v1.dq i) (v2.dp i) - Vec3.dot (v2.dq i) (v1.dp i)

def omega0 (v1 v2 : TangentVector3D N k) : k := sumFin N (fun i => omega0At v1 v2 i)
