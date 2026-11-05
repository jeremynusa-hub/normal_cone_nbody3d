/-
  Project: Unsoftened Normal-Cone Symplectic Integration of the 3D Newtonian N-Body Problem (v1.0.0)
  Author: Jérémy Nusa
  Companion Formal Verification Package: NormalConeNBody (Lean 4 + Mathlib 4)

  Structure:
  1. `NormalConeNBody/Ch01_Introduction/` .. `Ch05_SymplecticSchemes/`
     Contains the 35 formal Lean 4 modules corresponding 1-to-1 (bit-for-bit SHA-256 identity)
     with Appendix A (§A.1 to §A.20) of the monograph.
  2. `NormalConeNBody/CertifiedSuites/`
     Contains the 5 self-contained chapter compilation suites (`Sec01_Introduction.lean` to
     `Sec05_SymplecticSchemes.lean`) configured as the root targets of `lakefile.lean`,
     with explicit `#print axioms` verification at the bottom of every chapter suite.
-/

import NormalConeNBody.CertifiedSuites.Sec03_ThreeBody
import NormalConeNBody.CertifiedSuites.Sec04_FultonMacPherson
import NormalConeNBody.CertifiedSuites.Sec05_SymplecticSchemes
