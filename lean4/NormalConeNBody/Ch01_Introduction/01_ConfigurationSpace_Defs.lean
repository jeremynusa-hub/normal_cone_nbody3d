-- Dependencies: none (foundational spatial primitives)
import Mathlib.Tactic.Ring
universe u

@[ext]
structure Vec3 (k : Type u) where
  x : k
  y : k
  z : k

namespace Vec3
variable {k : Type u} [CommRing k]
def zero : Vec3 k := ⟨0, 0, 0⟩
def add (v1 v2 : Vec3 k) : Vec3 k := ⟨v1.x + v2.x, v1.y + v2.y, v1.z + v2.z⟩
def smul (c : k) (v : Vec3 k) : Vec3 k := ⟨c * v.x, c * v.y, c * v.z⟩
def dot (v1 v2 : Vec3 k) : k := v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
def cross (v1 v2 : Vec3 k) : Vec3 k := ⟨v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x⟩
end Vec3

def AmbientConfiguration (N : Nat) (k : Type u) := Fin N → Vec3 k

structure PhaseState3D (N : Nat) (k : Type u) where
  q : Fin N → Vec3 k
  p : Fin N → Vec3 k
  m : Fin N → k

def SpatialTranslation {N : Nat} {k : Type u} [CommRing k] (q : AmbientConfiguration N k) (a : Vec3 k) : AmbientConfiguration N k :=
  fun i => Vec3.add (q i) a

def sumFin {k : Type u} [CommRing k] : (n : Nat) → (Fin n → k) → k
  | 0, _ => 0
  | n + 1, f => f 0 + sumFin n (fun i => f ⟨i.val + 1, Nat.succ_lt_succ i.isLt⟩)

def RelativeConfiguration (N : Nat) {k : Type u} [CommRing k] (m : Fin N → k) :=
  { q : AmbientConfiguration N k //
    (sumFin N (fun i => m i * (q i).x)) = 0 ∧
    (sumFin N (fun i => m i * (q i).y)) = 0 ∧
    (sumFin N (fun i => m i * (q i).z)) = 0 }
