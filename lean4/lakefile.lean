import Lake
open Lake DSL

package «normal_cone_nbody» where
  leanOptions := #[
    ⟨`pp.unicode.fun, true⟩,
    ⟨`autoImplicit, false⟩,
    ⟨`relaxedAutoImplicit, false⟩
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git"

@[default_target]
lean_lib «NormalConeNBody» where
  roots := #[
    `NormalConeNBody.CertifiedSuites.Sec01_Introduction,
    `NormalConeNBody.CertifiedSuites.Sec02_NormalCone,
    `NormalConeNBody.CertifiedSuites.Sec03_ThreeBody,
    `NormalConeNBody.CertifiedSuites.Sec04_FultonMacPherson,
    `NormalConeNBody.CertifiedSuites.Sec05_SymplecticSchemes
  ]
