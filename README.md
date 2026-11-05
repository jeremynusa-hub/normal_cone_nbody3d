# Normal-Cone Is All You Need: Unsoftened Fulton–MacPherson Symplectic Integration of the 3D Newtonian $N$-Body Problem (`v1.0.0`)

[![Lean 4 Formal Proofs](https://img.shields.io/badge/Lean_4_(v4.35.0--rc2)-35_Modules_%7C_416_Decls_%7C_0_sorryAx-0072B2)](lean4/README.md)
[![Python Symplectic Suite](https://img.shields.io/badge/Python_3.10%2B-Float64_(Numba_%40njit)_%26_Float128_(116--bit)-009E73)](pyproject.toml)
[![Unit & Showcase Tests](https://img.shields.io/badge/Pytest_Suite-330%2F330_Passing-2CA02C)](scripts/run_smoke_tests.py)
[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0_(Non--Commercial)-D97706)](LICENSE)

**Author:** Jérémy Nusa  
**Monograph Version:** `v1.0.0` (`293` pages, `33` Figures, `11` Tables, `Appendix A: §A.1–§A.21`)

This repository is the official, self-contained **companion formal verification, dual-precision software, and reproducible benchmark archive** for the monograph:

> **Jérémy Nusa**, *Normal-Cone Is All You Need: Unsoftened Fulton–MacPherson Symplectic Integration of the 3D Newtonian $N$-Body Problem* (`v1.0.0`, 2026).

It unites **three certified pillars** under a single root directory:
1. **Pillar I — Formal Verification in `Lean 4` + `Mathlib 4` ([`lean4/`](lean4/README.md), Appendix `§A.1–§A.20`)**:
   - **35 formal `.lean` modules** (`bit-for-bit SHA-256` identical to the code listings injected in `Appendix A`, `pp. 187–286`) and **5 self-contained chapter compilation suites** (`NormalConeNBody/CertifiedSuites/Sec01..Sec05`), totaling **416 formal declarations (`214` theorems, `125` consolidated suite declarations)** and **99 explicit `#print axioms` checks** with **`0` `sorry`, `0` `sorryAx`, and `0` custom axioms**.
2. **Pillar II — Dual-Precision Symplectic $N$-Body Engines (`normal_cone_nbody_float64/` & `normal_cone_nbody_float128/`, Appendix `§A.21`)**:
   - Two strictly 1-to-1 mirrored 46-module Python packages implementing the unsoftened ($\epsilon = 0, d_0 = 0$) normal-cone symplectic integrator on the Fulton–MacPherson compactification $\mathcal{Q}[N]$:
     - **`normal_cone_nbody_float64/`**: Hardware IEEE-754 double-precision (`53`-bit) engine accelerated by LLVM `Numba @njit(fastmath=False)` micro-kernels (`62.54x` apex speedup at $N=64$, `59.83x` synchronous super-node speedup, `0.0 ULP` algebraic deviation);
     - **`normal_cone_nbody_float128/`**: Exact structural mirror in `116`-bit / `35`-significant-decimal-digit quadruple precision (`mpmath` / `gmpy2`) to expose asymptotic sixth-order convergence slopes ($\mathrm{BM}_6^{(11)}$ and $\mathrm{PBM}_4(\mathcal{Z}_4^{\mathrm{BM4}})$) down to $10^{-31}$ below the IEEE-754 `Float64` roundoff floor;
     - **`interfaces/`**: 7 standalone pedagogical reference facades (`pcr3bp_canonical_solver.py` and `01..06_*_Interface.py`);
     - **`tests/`**: **45 unit test suites (`330` tests: `265` Float64 + `65` Float128)** certifying all algebraic, symplectic, and Noether invariants.
3. **Pillar III — Reproducible Benchmark Campaigns, Certified JSON Datasets & 300 DPI Figure Plotters (`benchmarks/`, `results/campaigns/`, `figures/`, `scripts/`)**:
   - **8 benchmark campaign runners** (Chapters 6 & 7, including head-to-head confrontations against `REBOUND 5.2.2` `WHFast`, `MERCURIUS`, `TRACE`, and `IAS15` on protocols `B1`, `B2`, and `B4`), **14 certified JSON datasets**, and **8 publication plate plotters** (`30` PDF/PNG files at `300 DPI`).

---

## 1. Repository Architecture (`normal_cone_nbody3d_github/`)

```text
normal_cone_nbody3d_github/
├── README.md, CITATION.cff, .zenodo.json, LICENSE   # Academic metadata, Zenodo schema & CC BY-NC-SA 4.0 license
├── pyproject.toml, requirements.txt, environment.yml # PEP 621 packaging (pip & Conda/Mamba)
├── .github/workflows/ci.yml                         # Automated CI: Job 1 (Lean 4) & Job 2 (Python Float64/Float128)
├── scripts/                                         # 3-Tier Reproducibility Orchestrator for Referees
│   ├── reproduce_figures.py                         #   [Tier 1: <45s] Audits 14 JSONs & regenerates 30 plates (PDF/PNG)
│   ├── run_smoke_tests.py                           #   [Tier 2: <90s] Audits 35 Lean 4 modules + runs 330 pytest tests
│   └── run_full_campaigns.py                        #   [Tier 3: Full] Re-runs Generators 01–05 & REBOUND B1, B2, B4
├── lean4/                                           # [Pillar I] Lean 4 + Mathlib 4 Formal Proof Suite (§A.1–§A.20)
│   ├── lean-toolchain, lakefile.lean, lake-manifest.json, NormalConeNBody.lean, README.md
│   ├── scripts/verify_lean4_suite.py                #   Automated SHA-256 (35/35), declaration & 0-sorryAx auditor
│   └── NormalConeNBody/
│       ├── Ch01_Introduction/                       #   7 modules (§A.1–§A.7): Darboux 2-form & Bargmann algebra
│       ├── Ch02_NormalCone/                         #   6 modules (§A.8–§A.10): Normal-cone lift, splitting & Heavy Top
│       ├── Ch03_ThreeBody/                          #   2 modules (§A.11–§A.12): Synodic PCR3BP & Levi-Civita covering
│       ├── Ch04_FultonMacPherson/                   #   10 modules (§A.13–§A.15): Q[N] tree, C^{1,1} apex & super-nodes
│       ├── Ch05_SymplecticSchemes/                  #   10 modules (§A.16–§A.20): S2, Suzuki, BM4, BM6, PBM4 & L(5,2)
│       └── CertifiedSuites/                         #   5 self-contained chapter suites (Sec01..Sec05, 99 #print axioms)
├── normal_cone_nbody_float64/                       # [Pillar II-a] IEEE-754 Float64 + Numba @njit Engine (46 modules)
│   ├── core/                                        #   state.py, jacobi.py, coordinates.py, transitions.py, invariants.py
│   ├── potential/                                   #   screened_apex.py, optimized_screened_apex.py (@njit), screening.py
│   ├── solvers/                                     #   stumpff.py, kepler_sundman.py, optimized_kepler_sundman.py, jit_supernode.py
│   ├── integrators/                                 #   common.py (ContextVar), strang.py, yoshida.py, blanes_moan_{4,6}.py, processed.py
│   ├── system/                                      #   nbody_system.py, atlas_manager.py, dyadic_scheduler.py, history.py
│   ├── diagnostics/                                 #   metrics.py (Laskar NAFF, Brouwer law), collision_radar.py, visualizer.py
│   └── benchmarks/                                  #   registry.py (17 canonical N-body systems), runner.py
├── normal_cone_nbody_float128/                      # [Pillar II-b] 116-Bit / 35-Digit Float128 Mirror Engine (46 modules)
│   └── core/, potential/, solvers/, integrators/, system/, diagnostics/, benchmarks/
├── interfaces/                                      # Standalone Reference Facades (§3.3 & §6.1–§6.6, Appendix §A.21)
│   ├── pcr3bp_canonical_solver.py                   #   §3.3: Planar Synodic PCR3BP Normal-Cone Solver (Arenstorf Orbit)
│   └── 01..06_*_Interface.py                        #   §6.1–§6.6: Phase state, tree/Z4, dyadic scheduler, subcycle, atlas, @njit
├── benchmarks/                                      # Reproducible Benchmark Campaign Generators (Chapters 6 & 7, §A.21)
│   ├── contracts/00..05_*.md                        #   Formal benchmark verification contracts
│   ├── run_generator_01_dual_convergence_and_microkernels.py
│   ├── run_generator_02_universal_obstruction_and_shadow_gauge.py
│   ├── run_generator_03_unsoftened_singularity_limit_cases.py
│   ├── run_generator_04_historical_encounters_and_atlas_gauge.py
│   ├── run_generator_05_secular_stability_laplace_and_laskar_naff.py
│   ├── run_rebound_b1_star_grazing_pericenter.py
│   ├── run_rebound_b2_hill_crossings.py
│   └── run_rebound_b4_pythagorean_chaos.py
├── figures/                                         # Publication Plate Plotters (300 DPI PDF/PNG from certified JSONs)
│   └── plot_generator_01..05_*.py, plot_rebound_b1..b4_*.py
├── results/campaigns/                               # 14 Certified Reference JSON Datasets (Chapters 6 & 7)
└── tests/                                           # Automated Test Suites (330 tests across 45 test modules)
    ├── unit/                                        #   265 Float64 unit & Appendix A showcase interface tests
    └── unit_float128/                               #   65 Float128 (35-digit quadruple precision) unit tests
```

---

## 2. Three-Speed Reproducibility Quickstart (For Referees & Readers)

### Installation (`< 1 minute`)

```bash
# Option A: Standard Python virtual environment (pip)
pip install -e ".[dev]"

# Option B: Conda / Mamba environment
conda env create -f environment.yml
conda activate normal-cone-nbody3d
```

### Speed 1 — Instant Audit of Certified JSONs & Regeneration of All Publication Figures (`~ 45 seconds`)

Verifies the integrity of all **14 certified JSON campaign datasets** in `results/campaigns/` and executes all **8 figure plotters** (`figures/plot_*.py`) to regenerate the **30 publication plates (`300 DPI` PDF & PNG)** in `results/figures/`:

```bash
python scripts/reproduce_figures.py
```

### Speed 2 — Complete Formal `Lean 4` Audit & Dual-Precision Python Test Suite (`~ 85 seconds`)

Executes sequentially: (1) the syntactic, declaration, and `0-sorryAx` audit of all **35 `Lean 4` modules** (`§A.1–§A.20`), (2) the **Appendix A standalone interface tests** (including the Arenstorf synodic orbit), (3) all **265 `Float64` unit tests**, and (4) all **65 `Float128` (`35`-digit) unit tests**:

```bash
python scripts/run_smoke_tests.py
```

To also compile the 5 self-contained `Lean 4` chapter suites inside the `Lean 4` (`v4.35.0-rc2`) + `Mathlib 4` kernel and print all `99` `#print axioms` certificates:

```bash
cd lean4
lake exe cache get
lake build
```

### Speed 3 — Full Numerical Benchmark Campaign Re-Execution (`Chapters 6 & 7`)

```bash
# List all 8 available benchmark campaigns (gen01..gen05, rebound_b1, rebound_b2, rebound_b4)
python scripts/run_full_campaigns.py --list

# Re-run a specific campaign (e.g., Generator 01 or REBOUND B1) or all campaigns
python scripts/run_full_campaigns.py --campaign gen01
python scripts/run_full_campaigns.py --campaign all
```

---

## 3. Unified Bi-Modal Traceability Matrix (`Monograph Sections 1–8` $\leftrightarrow$ `Lean 4 §A.1–§A.20` $\leftrightarrow$ `Python & Benchmarks §A.21`)

| Monograph Section | Mathematical, Geometric, or Algorithmic Result | Formal `Lean 4` Proof Modules (`lean4/NormalConeNBody/...`, `§A.1–§A.20`) | Python Production Modules, Interfaces & Benchmark Suites (`§A.21`) | Figures & Tables in `v1.0.0` |
| :--- | :--- | :--- | :--- | :--- |
| **§1.1–§1.3** | Configuration space $\mathbb{R}^{3N}$, collision locus $\Delta$, Darboux 2-form $\omega_0$, and 10-generator Galilean–Bargmann Lie–Poisson algebra | **`§A.1–§A.7`** (`Ch01_Introduction/01..07_*.lean` & `CertifiedSuites/Sec01_Introduction.lean`, `53` decls, `25` thms) | `core/state.py`, `core/invariants.py` (`compute_lie_bargmann_invariants`) | `Fig. 1` |
| **§2.1–§2.3** | Deformation to the normal cone $\mathcal{D}(Q,\Delta)$, shifted cotangent lift $\omega_{\mathcal{N}}$, 4-step splitting recipe, and Heavy Top $\mathfrak{se}(3)^*$ Casimirs $(C_1, C_2)$ | **`§A.8–§A.10`** (`Ch02_NormalCone/01..06_*.lean` & `CertifiedSuites/Sec02_NormalCone.lean`, `37` decls, `17` thms) | `solvers/kepler_sundman.py` (Rodrigues / Sundman normal-cone fibers) | `Figs. 2–5` |
| **§3.1–§3.3** | Planar synodic PCR3BP ($\dot{C}_J \equiv 0$, $\mathrm{Sp}(4)$ Coriolis propagator), Levi-Civita covering $z=u^2$, and Arenstorf periodic orbit validation | **`§A.11–§A.12`** (`Ch03_ThreeBody/01..02_*.lean` & `CertifiedSuites/Sec03_ThreeBody.lean`, `10` decls, `6` thms) | `interfaces/pcr3bp_canonical_solver.py` (`PCR3BPNormalConeIntegrator`), `tests/unit/test_paper_showcase_interfaces.py` | `Figs. 6–8`, `Table 1` |
| **§4.1** | Partition lattice $\Pi(N)$, codimension $3(k-1)$, cluster hyper-radius $\rho_S$, binary Jacobi trees $\mathcal{T}$, and disjoint leaf commutation | **`§A.13`** (`Ch04_FultonMacPherson/01..04_*.lean`, `25` decls, `16` thms) | `core/jacobi.py`, `core/coordinates.py` (`JacobiTree`, `build_optimal_tree`) | `Figs. 9–10` |
| **§4.2** | Screened apex potential $W_{\mathrm{apex}}^{\mathrm{screen}}$, universal $C^{1,1}$ Lipschitz/Hessian bounds, Lie–Bargmann invariance, and extended Sundman symplecticity | **`§A.14`** (`Ch04_FultonMacPherson/01..03_Apex*.lean`, `24` decls, `10` thms) | `potential/screened_apex.py`, `potential/screening.py` (`ScreenedApexPotential`, `septic_screening_theta`) | `Fig. 11` |
| **§4.3** | Synchronous super-nodes ($k \ge 3$), collective Sundman clock $\mathrm{d}t_S = \bar{\rho}_S\mathrm{d}\tau_S$, anti-Zeno scaling, $\mathrm{SL}(2,\mathbb{R})$ chart flips, and hysteresis $\gamma=0.60$ | **`§A.15`** (`Ch04_FultonMacPherson/01..03_*SuperNode*.lean`, `18` decls, `10` thms) | `solvers/supernode.py`, `potential/supernode_detector.py`, `core/transitions.py`, `system/atlas_manager.py` | `Fig. 12` |
| **§5.1–§5.4** | Palindromic splitting hierarchy: Strang $\mathcal{S}_2$, Sheng–Suzuki barrier ($w_0<0$), Yoshida $\mathcal{S}_4$, 7-stage $\mathrm{BM}_4^{(7)}$ FSAL, and 11-stage $\mathrm{BM}_6^{(11)}$ | **`§A.16–§A.19`** (`Ch05_SymplecticSchemes/01..02_*.lean`, `119` decls, `54` thms) | `integrators/{strang, yoshida, blanes_moan_4, blanes_moan_6}.py`, `system/dyadic_scheduler.py` | `Figs. 13–14`, `Table 2` |
| **§5.5** | Universal Order-6 Obstruction Theorem in $\mathcal{L}(5,2)$, cotangent corrector $K_4 = c_4^{\mathrm{BM4}} h^4 D_A^3 B \in \mathcal{P}_1$, conjugator $\mathcal{Z}_4^{\mathrm{BM4}}(h)$, and gauge connection $T_{A\to B}^{(\mathrm{PBM4})}$ | **`§A.20`** (`Ch05_SymplecticSchemes/01,03,06,07,08_*.lean`, `130` decls, `76` thms) | `integrators/correctors.py` (`PalindromicCorrectorZ4`), `integrators/processed.py` (`ProcessedBM4Integrator`) | `Fig. 15` |
| **§6.1–§6.5** | Six-phase algorithmic pipeline, Delaunay–Hill tree initialization, symmetric dyadic scheduler ($k_{\mathrm{sym}}$), leaf/apex decoupling, and Voronoi atlas | *(Certified by `§A.13–§A.20`)* | `system/{nbody_system, atlas_manager, dyadic_scheduler}.py`, `interfaces/01..05_*_Interface.py` | `Fig. 16`, `Table 3` |
| **§6.6** | High-performance vectorization & LLVM `Numba @njit(fastmath=False)` micro-kernels ($N \in [3, 64]$, `62.54x` apex, `59.83x` super-node, `0.0 ULP`) | *(Bit-for-bit `0.0 ULP` algebraic parity with `§A.14–§A.15`)* | `potential/optimized_screened_apex.py`, `solvers/{optimized_kepler_sundman, optimized_supernode, jit_supernode}.py`, `integrators/common.py`, `interfaces/06_*_Interface.py` | `Fig. 17`, `Tables 4–5` |
| **§7.1.1** | Dual-precision symplectic convergence (`Float64` vs `35`-digit `Float128`) on Chenciner–Montgomery Figure-8 and 3D eccentric triple | *(Empirical realization of `§A.16–§A.20`)* | `benchmarks/run_generator_01_dual_convergence_and_microkernels.py`, `figures/plot_generator_01_*.py` | `Figs. 18–20`, `Table 6` |
| **§7.1.2** | Two-channel truncation law, order $6 \to 4$ crossover $h_{\times}(\varepsilon)$, and universal quadratic floor $\mathcal{O}(\varepsilon^2 h^4)$ (`ST-1..ST-4`) | *(Empirical verification of `§A.20`)* | `benchmarks/run_generator_02_universal_obstruction_and_shadow_gauge.py`, `figures/plot_generator_02_*.py` | `Figs. 21–23` |
| **§7.1.3** | Unsoftened singularity limit cases (`L1` $e\to 1$, `L2` $\mathcal{Z}_4^{\mathrm{BM4}}$ free-fall & McGehee $k\in\{2,3,4\}$, `L3` Burrau & $D_4$, `L4` Kozai–Lidov) | *(Empirical verification of `§A.13–§A.15`)* | `benchmarks/run_generator_03_unsoftened_singularity_limit_cases.py`, `figures/plot_generator_03_*.py` | `Figs. 24–27`, `Table 7` |
| **§7.2.1** | Historical planetary encounters & non-commutative atlas gauge connection (`H1` Lexell 1770, `H2` Shoemaker–Levy 9 1992, `H3` Apophis 2029) | *(Empirical verification of `§A.15` & `§A.20`)* | `benchmarks/run_generator_04_historical_encounters_and_atlas_gauge.py`, `figures/plot_generator_04_*.py` | `Figs. 28–30`, `Table 8` |
| **§7.2.2** | Secular stability, Brouwer's $\mathcal{O}(N_{\mathrm{steps}}^{1/2})$ law & Laskar NAFF (`S1` Jovian Laplace $4:2:1$, `S2` 264 dyadic junctions) | *(Empirical verification of `§A.16–§A.20`)* | `benchmarks/run_generator_05_secular_stability_laplace_and_laskar_naff.py`, `figures/plot_generator_05_*.py` | `Figs. 29–30`, `Table 9` |
| **§7.3** | Head-to-head confrontation against `REBOUND 5.2.2` (`WHFast`, `MERCURIUS`, `TRACE`, `IAS15`) on `B1` Star-Grazing, `B2` Hill Crossings & `B4` Burrau 3-4-5 | *(Cross-code validation of `§A.13–§A.20`)* | `benchmarks/run_rebound_{b1_star_grazing_pericenter, b2_hill_crossings, b4_pythagorean_chaos}.py`, `figures/plot_rebound_b1..b4_*.py` | `Figs. 31–33`, `Tables 10–11` |

---

## 4. License (`CC BY-NC-SA 4.0` — Academic & Non-Commercial Research Only)

Copyright © 2026 **Jérémy Nusa**. All Commercial Rights Reserved.

- This repository (including all Python source code, micro-kernels, 6-phase pipeline architecture, `Lean 4` formal proofs, benchmark datasets, and figures) is licensed under the **Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (`CC BY-NC-SA 4.0`)** supplemented by the explicit non-commercial software and derivative translation clauses in [`LICENSE`](LICENSE).
- **Academic & Peer-Review Use**: Researchers, universities, public research institutes, and referees are freely permitted to inspect, execute, verify, and build upon this work for non-commercial scientific research with proper citation.
- **Commercial Prohibition**: Any commercial use, corporate R&D deployment, commercial cloud/SaaS hosting, or commercial derivative translation/port into another language (`C`, `C++`, `Rust`, `Fortran`, `Julia`, `CUDA`, etc.) is strictly prohibited without a prior written commercial license from **Jérémy Nusa**.

---

## 5. Citation (`BibTeX`)

```bibtex
@misc{nusa2026normalcone,
  author       = {Nusa, J{\'e}r{\'e}my},
  title        = {Normal-Cone Is All You Need: Unsoftened {F}ulton--{M}ac{P}herson Symplectic Integration of the 3{D} {N}ewtonian {$N$}-Body Problem},
  year         = {2026},
  version      = {1.0.0},
  note         = {Monograph (293 pp.) and companion formal Lean 4 \& dual-precision Python repository}
}
```

