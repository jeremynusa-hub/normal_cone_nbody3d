-- Dependencies: 01_ConfigurationSpace_Defs, 04_SymplecticForm_Defs, 05_SymplecticForm_Theorems (sumFin_neg, sumFin_zero), 06_BargmannAlgebra_Defs
import Mathlib.Tactic.Ring
universe u

namespace BargmannAlgebra
variable {N : Nat} {k : Type u} [CommRing k]

theorem poisson_P_G_same (a : Axis3) (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (linearMomentumComponent a) (centerOfMassComponent a t) s = - totalMass s := by
  dsimp [poissonBracketVal, totalMass]
  have h_rw : (fun i => poissonBracketAt (linearMomentumComponent a) (centerOfMassComponent a t) s i) = (fun i => - s.m i) := by
    funext i; dsimp [poissonBracketAt, linearMomentumComponent, centerOfMassComponent]
    cases a <;> { dsimp [Axis3.unitVec, Vec3.zero, Vec3.smul, Vec3.dot]; ring }
  rw [h_rw]; exact sumFin_neg N (fun i => - s.m i) (fun i => s.m i) (fun _ => rfl)

theorem poisson_P_G_diff {a b : Axis3} (h : a ≠ b) (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (linearMomentumComponent a) (centerOfMassComponent b t) s = 0 := by
  dsimp [poissonBracketVal]; apply sumFin_zero; intro i
  dsimp [poissonBracketAt, linearMomentumComponent, centerOfMassComponent]
  cases a <;> cases b <;> try contradiction
  all_goals { dsimp [Axis3.unitVec, Vec3.zero, Vec3.smul, Vec3.dot]; ring }

theorem poisson_L_self (a : Axis3) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent a) (angularMomentumComponent a) s = 0 := by
  dsimp [poissonBracketVal]; apply sumFin_zero; intro i
  dsimp [poissonBracketAt]
  have h_comm : Vec3.dot ((angularMomentumComponent a).grad_q s i) ((angularMomentumComponent a).grad_p s i) =
                Vec3.dot ((angularMomentumComponent a).grad_p s i) ((angularMomentumComponent a).grad_q s i) := by
    dsimp [Vec3.dot]; ring
  rw [h_comm, sub_self]

theorem poisson_L_X_L_Y (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.X) (angularMomentumComponent Axis3.Y) s =
      (angularMomentumComponent Axis3.Z).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Axis3.proj]
  ring

theorem poisson_L_Y_L_Z (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Y) (angularMomentumComponent Axis3.Z) s =
      (angularMomentumComponent Axis3.X).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Axis3.proj]
  ring

theorem poisson_L_Z_L_X (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Z) (angularMomentumComponent Axis3.X) s =
      (angularMomentumComponent Axis3.Y).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Axis3.proj]
  ring

theorem poisson_L_P_same (a : Axis3) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent a) (linearMomentumComponent a) s = 0 := by
  dsimp [poissonBracketVal]; apply sumFin_zero; intro i
  dsimp [poissonBracketAt, angularMomentumComponent, linearMomentumComponent]
  cases a <;> { dsimp [Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.zero]; ring }

theorem poisson_L_X_P_Y (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.X) (linearMomentumComponent Axis3.Y) s =
      (linearMomentumComponent Axis3.Z).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, linearMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.zero, Axis3.proj]
  ring

theorem poisson_L_Y_P_Z (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Y) (linearMomentumComponent Axis3.Z) s =
      (linearMomentumComponent Axis3.X).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, linearMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.zero, Axis3.proj]
  ring

theorem poisson_L_Z_P_X (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Z) (linearMomentumComponent Axis3.X) s =
      (linearMomentumComponent Axis3.Y).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, linearMomentumComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.zero, Axis3.proj]
  ring

theorem poisson_L_G_same (a : Axis3) (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent a) (centerOfMassComponent a t) s = 0 := by
  dsimp [poissonBracketVal]; apply sumFin_zero; intro i
  dsimp [poissonBracketAt, angularMomentumComponent, centerOfMassComponent]
  cases a <;> { dsimp [Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.smul]; ring }

theorem poisson_L_X_G_Y (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.X) (centerOfMassComponent Axis3.Y t) s =
      (centerOfMassComponent Axis3.Z t).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, centerOfMassComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.smul, Axis3.proj]
  ring

theorem poisson_L_Y_G_Z (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Y) (centerOfMassComponent Axis3.Z t) s =
      (centerOfMassComponent Axis3.X t).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, centerOfMassComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.smul, Axis3.proj]
  ring

theorem poisson_L_Z_G_X (t : k) (s : PhaseState3D N k) :
    poissonBracketVal (angularMomentumComponent Axis3.Z) (centerOfMassComponent Axis3.X t) s =
      (centerOfMassComponent Axis3.Y t).val s := by
  dsimp [poissonBracketVal, angularMomentumComponent, centerOfMassComponent]
  congr 1; funext i
  dsimp [poissonBracketAt, Axis3.unitVec, Vec3.cross, Vec3.dot, Vec3.smul, Axis3.proj]
  ring
end BargmannAlgebra

