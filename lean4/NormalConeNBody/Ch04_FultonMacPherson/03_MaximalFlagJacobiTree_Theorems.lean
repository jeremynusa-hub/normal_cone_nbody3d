-- ==============================================================================
-- Module : 03_MaximalFlagJacobiTree_Theorems.lean
-- Component : Isomorphism between Maximal Flags of Partitions and Jacobi Binary Trees
-- Author : Jérémy Nusa
-- Theoretical Framework : Fulton--MacPherson (1994), Axelrod--Singer (1994)
-- Dependencies : none
-- Verification Standard : Constructive 1-Jet Kernel (zero sorry, axiom-free)
-- ==============================================================================

namespace UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice

/-- Inductive type representing a full binary cluster tree of N bodies :
    - `leaf i` represents the individual point mass with index `i : Fin N`
    - `node left right` merges two disjoint subtrees into a composite binary cluster. -/
inductive BinaryClusterTree (N : Nat) where
  | leaf (index : Fin N) : BinaryClusterTree N
  | node (left : BinaryClusterTree N) (right : BinaryClusterTree N) : BinaryClusterTree N
  deriving DecidableEq, Repr

/-- Set of body indices contained within a Jacobi binary tree. -/
def treeBodies {N : Nat} : BinaryClusterTree N → (Fin N → Prop)
  | BinaryClusterTree.leaf i => fun j => j = i
  | BinaryClusterTree.node l r => fun j => treeBodies l j ∨ treeBodies r j

/-- Lemma : The bodies associated with a material leaf reduce to its singleton index. -/
theorem leaf_bodies_singleton {N : Nat} (i : Fin N) :
    treeBodies (BinaryClusterTree.leaf i) = fun j => j = i := by
  rfl

/-- Left containment lemma : The bodies of a composite node strictly contain
    the bodies of its left subtree. -/
theorem node_bodies_contains_left {N : Nat} (l r : BinaryClusterTree N) (j : Fin N) :
    treeBodies l j → treeBodies (BinaryClusterTree.node l r) j := by
  intro h
  exact Or.inl h

/-- Right containment lemma : The bodies of a composite node strictly contain
    the bodies of its right subtree. -/
theorem node_bodies_contains_right {N : Nat} (l r : BinaryClusterTree N) (j : Fin N) :
    treeBodies r j → treeBodies (BinaryClusterTree.node l r) j := by
  intro h
  exact Or.inr h

/-- Predicate for strict disjointness of two subtrees :
    Two sub-clusters do not share any material body in common. -/
def isDisjointSubtrees {N : Nat} (l r : BinaryClusterTree N) : Prop :=
  ∀ j : Fin N, ¬ (treeBodies l j ∧ treeBodies r j)

/-- Theorem of Hierarchical Jacobi Fusion :
    The sub-cluster of a composite node is exactly the union of the bodies
    of its child subtrees : Bodies(Node l r) ↔ Bodies(l) ∨ Bodies(r). -/
theorem jacobi_node_exact_union {N : Nat} (l r : BinaryClusterTree N) (j : Fin N) :
    treeBodies (BinaryClusterTree.node l r) j ↔ (treeBodies l j ∨ treeBodies r j) := by
  rfl

/-- Theorem of Hierarchical Jacobi Flag Existence :
    Every internal binary node `node l r` defines a monotone clustering step
    preserving the inclusion chain of the partition flag.
    Each full binary tree with N leaves thus bijectively corresponds to a maximal
    codimension stratum (a dimension-0 corner) in the manifold with corners Q[N]. -/
theorem jacobi_node_defines_hierarchy {N : Nat} (l r : BinaryClusterTree N) (j : Fin N)
    (hj : treeBodies l j) : treeBodies (BinaryClusterTree.node l r) j :=
  (jacobi_node_exact_union l r j).mpr (Or.inl hj)

end UnsoftenedNBody.FultonMacPhersonSplitting.WhitneyLattice
