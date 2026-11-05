# Lean 4 Formal Verification Suite (`NormalConeNBody` — Appendix A, `§A.1`, `Tables 10–11`)

This directory contains the self-contained **Lean 4 + Mathlib 4** formal verification package accompanying the monograph:

> **Normal-Cone Is All You Need: Unsoftened Fulton–MacPherson Symplectic Integration of the 3D Newtonian $N$-Body Problem (`v1.0.0`)**  
> **Author:** Jérémy Nusa  
> **Zenodo DOI:** [`10.5281/zenodo.23982968`](https://doi.org/10.5281/zenodo.23982968) | **GitHub:** [`jeremynusa-hub/normal_cone_nbody3d`](https://github.com/jeremynusa-hub/normal_cone_nbody3d)

---

## 1. Epistemic Scope & Zero-`sorryAx` Guarantee

- **35 Formal Modules (`§A.1`, `Tables 10–11`)**: Every `.lean` module under `NormalConeNBody/Ch01_Introduction/` through `NormalConeNBody/Ch05_SymplecticSchemes/` corresponds 1-to-1 with the module locator in **Appendix A (`§A.1`, `Tables 10–11`, `pp. 187–188`)** of the monograph.
- **416 Formal Declarations (`214` Theorems) & `0` `sorryAx`**: Across the 35 modules (`125` consolidated declarations in the 5 chapter compilation suites), every proof is closed constructively inside the Lean 4 kernel (`v4.35.0-rc2`) and Mathlib 4 without a single `sorry`, `sorryAx`, `admit`, or custom `axiom`.
- **Automated `#print axioms` Certification**: The 5 self-contained chapter compilation units in `NormalConeNBody/CertifiedSuites/` execute `99` explicit `#print axioms` checks, confirming that every theorem depends exclusively on Lean 4's standard foundational axioms (`propext`, `Quot.sound`, `Classical.choice`).

---

## 2. Directory Layout

```text
lean4/
├── lean-toolchain               # Pinned Lean 4 toolchain (leanprover/lean4:v4.35.0-rc2)
├── lakefile.lean                # Lake package configuration (NormalConeNBody)
├── lake-manifest.json           # Pinned Mathlib 4 dependency manifest
├── NormalConeNBody.lean         # Top-level library entry point
├── README.md                    # This documentation and traceability table
├── scripts/
│   └── verify_lean4_suite.py    # Automated SHA-256, declaration & 0-sorryAx auditor
└── NormalConeNBody/
    ├── Ch01_Introduction/       # 7 modules (§A.1, Table 10): Darboux 2-form & Bargmann algebra
    ├── Ch02_NormalCone/         # 6 modules (§A.1, Table 10): Normal-cone lift, splitting & Heavy Top
    ├── Ch03_ThreeBody/          # 2 modules (§A.1, Table 10): Synodic PCR3BP & Levi-Civita covering
    ├── Ch04_FultonMacPherson/   # 10 modules (§A.1, Tables 10–11): Q[N] tree, C^{1,1} apex & super-nodes
    ├── Ch05_SymplecticSchemes/  # 10 modules (§A.1, Table 11): S2, Suzuki, BM4, BM6, PBM4 & L(5,2)
    └── CertifiedSuites/         # 5 self-contained compilable chapter suites with #print axioms
        ├── Sec01_Introduction.lean
        ├── Sec02_NormalCone.lean
        ├── Sec03_ThreeBody.lean
        ├── Sec04_FultonMacPherson.lean
        └── Sec05_SymplecticSchemes.lean
```

---

## 3. Quick Verification Commands

### Fast Syntactic, Declaration & `SHA-256` Audit (`< 1 second`, Python 3 only)

```bash
python lean4/scripts/verify_lean4_suite.py
```

### Full Lean 4 + Mathlib 4 Kernel Compilation & `#print axioms` Audit

```bash
cd lean4
lake exe cache get
lake build
```

Or compile any individual chapter suite directly:

```bash
cd lean4
lake env lean NormalConeNBody/CertifiedSuites/Sec01_Introduction.lean
lake env lean NormalConeNBody/CertifiedSuites/Sec02_NormalCone.lean
lake env lean NormalConeNBody/CertifiedSuites/Sec03_ThreeBody.lean
lake env lean NormalConeNBody/CertifiedSuites/Sec04_FultonMacPherson.lean
lake env lean NormalConeNBody/CertifiedSuites/Sec05_SymplecticSchemes.lean
```

---

## 4. Complete Traceability Matrix (`Appendix A: §A.1, Tables 10–11` $\leftrightarrow$ `35 Lean 4 Modules`)

| # | Appendix Locator | Module Path (`NormalConeNBody/...`) | Decls (`Thm`) | Formalized Mathematical Result |
| :-: | :---: | :--- | :---: | :--- |
| 01 | `§A.1 (Table 10)` | `Ch01_Introduction/01_ConfigurationSpace_Defs.lean` | `11` (`0`) | Ambient configuration space $\mathbb{R}^{3N}$ and barycentric reduction |
| 02 | `§A.1 (Table 10)` | `Ch01_Introduction/02_CollisionDivisor_Defs.lean` | `3` (`0`) | Pairwise collision diagonals $\Delta_{ij}$ and singular locus $\Delta$ |
| 03 | `§A.1 (Table 10)` | `Ch01_Introduction/03_CollisionDivisor_Theorems.lean` | `2` (`2`) | Symmetry and complement duality of $Q_{\mathrm{reg}} = \mathbb{R}^{3N} \setminus \Delta$ |
| 04 | `§A.1 (Table 10)` | `Ch01_Introduction/04_SymplecticForm_Defs.lean` | `3` (`0`) | Canonical cotangent bundle $T^*Q$ and Darboux 2-form $\omega_0$ |
| 05 | `§A.1 (Table 10)` | `Ch01_Introduction/05_SymplecticForm_Theorems.lean` | `10` (`9`) | Skew-symmetry and global non-degeneracy of $\omega_0$ |
| 06 | `§A.1 (Table 10)` | `Ch01_Introduction/06_BargmannAlgebra_Defs.lean` | `10` (`0`) | Ten classical integrals $(\mathbf{P}, \mathbf{J}, \mathbf{G}, H)$ and Poisson bracket |
| 07 | `§A.1 (Table 10)` | `Ch01_Introduction/07_BargmannAlgebra_Theorems.lean` | `14` (`14`) | Exact Lie–Poisson commutation relations of the Bargmann algebra |
| 08 | `§A.1 (Table 10)` | `Ch02_NormalCone/01_ReesDeformation_Defs.lean` | `5` (`0`) | Deformation to the normal cone $\mathcal{D}(Q, \Delta)$ and shifted cotangent lift |
| 09 | `§A.1 (Table 10)` | `Ch02_NormalCone/02_ReesDeformation_Theorems.lean` | `5` (`5`) | Symplecticity and non-degeneracy of the normal-cone 2-form $\omega_{\mathcal{N}}$ |
| 10 | `§A.1 (Table 10)` | `Ch02_NormalCone/03_CanonicalSplitting_Defs.lean` | `5` (`0`) | Four-step normal-cone symplectic splitting recipe |
| 11 | `§A.1 (Table 10)` | `Ch02_NormalCone/04_CanonicalSplitting_Theorems.lean` | `5` (`5`) | Exact Casimir invariance and time-reversibility of the Strang splitting |
| 12 | `§A.1 (Table 10)` | `Ch02_NormalCone/05_HeavyTop_Defs.lean` | `10` (`0`) | Heavy top $\mathfrak{se}(3)^*$ normal-cone fiber bundle and Rodrigues rotation |
| 13 | `§A.1 (Table 10)` | `Ch02_NormalCone/06_HeavyTop_Theorems.lean` | `7` (`7`) | Exact preservation of both Casimirs $C_1 = \|\boldsymbol{\Gamma}\|^2$ and $C_2 = \mathbf{\Pi}\cdot\boldsymbol{\Gamma}$ |
| 14 | `§A.1 (Table 10)` | `Ch03_ThreeBody/01_PCR3BPHamiltonian_Theorems.lean` | `5` (`2`) | Exact Jacobi integral identity $\dot{C}_J \equiv 0$ and $\mathrm{Sp}(4)$ Coriolis step |
| 15 | `§A.1 (Table 10)` | `Ch03_ThreeBody/02_PCR3BPRegularization_Theorems.lean` | `5` (`4`) | Levi-Civita double covering $z = u^2$ and conformal scaling $\|J(u)\|^2 = 4|u|^2$ |
| 16 | `§A.1 (Table 10)` | `Ch04_FultonMacPherson/01_PartitionLattice_Defs.lean` | `6` (`4`) | Partition lattice $\Pi(N)$, codimension $3(k-1)$, and non-transversality |
| 17 | `§A.1 (Table 10)` | `Ch04_FultonMacPherson/02_ClusterHyperRadius_Theorems.lean` | `7` (`4`) | Cluster hyper-radius $\rho_S$ and McGehee potential factorization |
| 18 | `§A.1 (Table 10)` | `Ch04_FultonMacPherson/03_MaximalFlagJacobiTree_Theorems.lean` | `8` (`5`) | Binary Jacobi tree hierarchy and nested Fulton–MacPherson flags |
| 19 | `§A.1 (Table 10)` | `Ch04_FultonMacPherson/04_DisjointCotangentCommutation_Theorems.lean` | `4` (`3`) | Exact Poisson commutation $\{H_{\mathrm{leaf}}^{(S_1)}, H_{\mathrm{leaf}}^{(S_2)}\} = 0$ |
| 20 | `§A.1 (Table 11)` | `Ch04_FultonMacPherson/01_ApexUniversalLipschitzBound_Theorems.lean` | `9` (`3`) | Universal $C^{1,1}$ Lipschitz and Hessian bounds of $W_{\mathrm{apex}}^{\mathrm{screen}}$ |
| 21 | `§A.1 (Table 11)` | `Ch04_FultonMacPherson/02_LieBargmannExactInvariance_Theorems.lean` | `5` (`3`) | Exact Galilean–Bargmann invariance of the screened apex splitting |
| 22 | `§A.1 (Table 11)` | `Ch04_FultonMacPherson/03_ExtendedPalindromicSymplecticity_Theorems.lean` | `10` (`4`) | Extended Poincaré–Sundman symplecticity and palindromic reversibility |
| 23 | `§A.1 (Table 11)` | `Ch04_FultonMacPherson/01_UnifiedSuperNodeDef_Theorems.lean` | `4` (`2`) | Synchronous super-node collective Sundman step positivity |
| 24 | `§A.1 (Table 11)` | `Ch04_FultonMacPherson/02_SuperNodeAntiZenoRegularization_Theorems.lean` | `5` (`3`) | Anti-Zeno physical time increment and Keplerian $r \propto \tau^2$ scaling |
| 25 | `§A.1 (Table 11)` | `Ch04_FultonMacPherson/03_AtlasTransitionSymplecticAntiChattering_Theorems.lean` | `9` (`5`) | Symplectic $\mathrm{SL}(2,\mathbb{R})$ chart flips and hysteresis anti-chattering |
| 26 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/01_StrangOrderTwoSymplectic_Theorems.lean` | `17` (`7`) | Second-order palindromic symplecticity of the Strang $\mathrm{S}_2$ propagator |
| 27 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/02_LieBargmannConservation_Theorems.lean` | `22` (`12`) | Exact conservation of all 10 Lie–Bargmann generators under $\mathrm{S}_2$ |
| 28 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/01_YoshidaOrderFourSuzuki_Theorems.lean` | `7` (`4`) | Sheng–Suzuki barrier ($w_0 < 0$) and order-4 Yoshida triple jump |
| 29 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/01_BlanesMoanShadowConfinement_Theorems.lean` | `10` (`8`) | Seven-stage $\mathrm{BM}_4$ FSAL consistency and positive leaf steps |
| 30 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/01_BlanesMoanOrderSix_Theorems.lean` | `63` (`23`) | Eleven-stage $\mathrm{BM}_6$ order-6 consistency and RKN commutators |
| 31 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/01_CorrectorBM4_Theorems.lean` | `59` (`26`) | Symplectic corrector conjugation $\mathcal{Z}_4^{\mathrm{BM4}}(h)$ and centered equivalence |
| 32 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/03_ConormalCorrectorPBM4_Order6.lean` | `22` (`13`) | Single-chart sixth-order cancellation of $[K,[K,V]]$ and $[V,[V,K]]$ via $\mathrm{PBM}_4$ |
| 33 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/06_ConormalGaugeTransitionPBM4_Order6.lean` | `14` (`11`) | Inter-chart symplectic gauge connection $T_{A\to B}^{(\mathrm{PBM4})}$ and cocycle law |
| 34 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/07_GlobalConormalCorrectorPBM4_Order6.lean` | `17` (`13`) | Global atlas gauge-invariance theorem on the $\mathcal{P}_0$-dominated regime |
| 35 | `§A.1 (Table 11)` | `Ch05_SymplecticSchemes/08_UniversalOrderSixObstruction.lean` | `18` (`13`) | Bidegree-$(3,2)$ Hall Lie algebra $\mathcal{L}(5,2)$ cohomological obstruction beyond $\mathcal{P}_0$ |
