-- Dependencies: 01_ConfigurationSpace_Defs (Vec3, PhaseState3D), 04_SymplecticForm_Defs (sumFin)
import Mathlib.Tactic.Ring
universe u

inductive Axis3 | X | Y | Z deriving DecidableEq

namespace Axis3
variable {k : Type u} [CommRing k]
def unitVec (a : Axis3) : Vec3 k := match a with | Axis3.X => ⟨1, 0, 0⟩ | Axis3.Y => ⟨0, 1, 0⟩ | Axis3.Z => ⟨0, 0, 1⟩
def proj (a : Axis3) (v : Vec3 k) : k := match a with | Axis3.X => v.x | Axis3.Y => v.y | Axis3.Z => v.z
end Axis3

structure FormallySmoothObservable (N : Nat) (k : Type u) where
  val : PhaseState3D N k → k
  grad_q : PhaseState3D N k → Fin N → Vec3 k
  grad_p : PhaseState3D N k → Fin N → Vec3 k

def poissonBracketAt {N : Nat} {k : Type u} [CommRing k] (F G : FormallySmoothObservable N k) (s : PhaseState3D N k) (i : Fin N) : k :=
  Vec3.dot (F.grad_q s i) (G.grad_p s i) - Vec3.dot (F.grad_p s i) (G.grad_q s i)

def poissonBracketVal {N : Nat} {k : Type u} [CommRing k] (F G : FormallySmoothObservable N k) (s : PhaseState3D N k) : k :=
  sumFin N (fun i => poissonBracketAt F G s i)

namespace BargmannAlgebra
variable {N : Nat} {k : Type u} [CommRing k]

def linearMomentumComponent (a : Axis3) : FormallySmoothObservable N k where
  val s := sumFin N (fun i => Axis3.proj a (s.p i))
  grad_q _ _ := Vec3.zero
  grad_p _ _ := Axis3.unitVec a

def centerOfMassComponent (a : Axis3) (t : k) : FormallySmoothObservable N k where
  val s := sumFin N (fun i => s.m i * Axis3.proj a (s.q i) - t * Axis3.proj a (s.p i))
  grad_q s i := Vec3.smul (s.m i) (Axis3.unitVec a)
  grad_p _ _ := Vec3.smul (-t) (Axis3.unitVec a)

def totalMass (s : PhaseState3D N k) : k := sumFin N (fun i => s.m i)

def angularMomentumComponent (a : Axis3) : FormallySmoothObservable N k where
  val s := sumFin N (fun i => Axis3.proj a (Vec3.cross (s.q i) (s.p i)))
  grad_q s i := Vec3.cross (s.p i) (Axis3.unitVec a)
  grad_p s i := Vec3.cross (Axis3.unitVec a) (s.q i)
end BargmannAlgebra
