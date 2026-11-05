# Canonical Benchmark Generator Contract Template (`00_canonical_generator_contract_template.md`)

**Contract Governance Version**: Pillar II — Cycle V (Iteration 58, Arbitrages `A1`, `A2`, `A3`)  
**Reference Doctrine**: [`09_pillar2_campaign_audit_and_laskar_plate_doctrine.md`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/reviews/substantive_audit/09_pillar2_campaign_audit_and_laskar_plate_doctrine.md)  
**Applicable Generators (`Iterations 59–63`)**:
- `01_generator_dual_convergence_and_microkernels.md` (`Iteration 59`, `§6.6` & `§7.1.1`)
- `02_generator_universal_obstruction_and_shadow_gauge.md` (`Iteration 60`, `§7.1.2`)
- `03_generator_unsoftened_singularity_limit_cases.md` (`Iteration 61`, `§7.1.3`)
- `04_generator_historical_encounters_and_atlas_gauge.md` (`Iteration 62`, `§7.2.1`)
- `05_generator_secular_stability_laplace_and_laskar_naff.md` (`Iteration 63`, `§7.2.2`)

---

## Section 1 — Mathematical Objective, Target Theorems & Formal Lean 4 Anchors

1. **Generator Identifier & Title**: `Generator 0X — [Full Descriptive Title]`
2. **Target Manuscript Sections, Propositions & Tables**:
   - Primary LaTeX section(s): `§X.Y.Z` (`chapters/ch0X.../sec0Y.../content.tex`)
   - Target Table(s): `Table N`
   - Substantive Audit Findings Remediated: `[CRIT-X.Y-01]`, `[MAJ-X.Y-01]`, ... (cross-referenced to Reports `01–07`)
3. **Mathematical & Theoretical Propositions Tested**:
   - Explicit statement of the mathematical property under test (e.g., asymptotic order $p \in \{2, 4, 6\}$, Theorem 5.14 universal obstruction of $\mathcal{P}_0$, Theorem 5.15–5.17 cotangent lift $\mathcal{P}_1$ and palindromic Lie conjugator $\mathcal{Z}_4^{\mathrm{BM4}}$, Theorem 5.18 covariant gauge connection $T_{A\to B}^{(\mathrm{PBM4})}$, McGehee blow-up regularity at $d_0=0$, Brouwer unbiased random walk $\propto t^{1/2}$, Laskar NAFF frequency stability).
4. **Formal Lean 4 Proof Modules Linked (`Appendix A`)**:
   - Exact `.lean` file(s) and certified theorem identifier(s) corresponding to the numerical invariant verified by this generator.

---

## Section 2 — Paper Showcase Interfaces (`01..06_*.py`) & Dual-Precision Library Bindings

1. **Mandatory Paper Showcase Imports (`Float64`)**:
   - Every generator MUST import its `Float64` classes and functions through the colocated paper showcase interfaces certified in Iteration 57:
     - `pcr3bp_canonical_solver.py` (`Appendix A.13`)
     - `01_PhaseStateAndSystem_Interface.py` (`Appendix A.22`: `ExtendedPhaseState`, `NBodySystem`)
     - `02_TreeBuilderAndCorrector_Interface.py` (`Appendix A.23`: `build_optimal_tree`, `detect_supernodes`, `PalindromicCorrectorZ4`, `CotangentLiftCorrectorP1`, `SymplecticBoundaryCorrector`, `GlobalCorrectorC4`, `create_pbm4_integrator`)
     - `03_DyadicScheduler_Interface.py` (`Appendix A.24`: `DyadicScheduler`, `ColoredInterval`)
     - `04_SymplecticSubcycle_Interface.py` (`Appendix A.25`: `stumpff_c`, `step_kepler_sundman`, `step_supernode_synchronous`, `compute_intercluster_forces`, `StrangIntegrator`, `Yoshida4Integrator`, `BlanesMoan4Integrator`, `BlanesMoan6Integrator`)
     - `05_AtlasTransitions_Interface.py` (`Appendix A.26`: `SL2TransitionMatrix`, `AtlasManager`, `PBM4GaugeTransitionOperator`, `apply_gauge_connection`, `compute_lie_bargmann_invariants`)
     - `06_HighPerformanceOptimizations_Interface.py` (`Appendix A.27`: `OptimizedScreenedApexPotential`, `OptimizedKeplerSolver`, `OptimizedSupernodeSolver`, `JitSupernodeSolver`, `set_current_*`)
2. **Quadruple-Precision (`Float128`, 35 Digits `mpmath`) Mirror Bindings**:
   - Exact imports from `normal_cone_nbody_float128` with mandatory runtime assertion:
     ```python
     assert state_f128.is_f128 is True
     assert isinstance(state_f128.q[0, 0], mpmath.mpf)
     ```
3. **Explicit Configuration Flags**:
   - `in_shadow_mode`: `True` (default) vs `False` (when testing un-telescoped ablation)
   - `corrector.mode`: `"homological"` ($\mathcal{Z}_4^{\mathrm{BM4}} \in \mathcal{P}_1$, default), `"cotangent_lift"` ($\in \mathcal{P}_1$), or `"global"` ($\in \mathcal{P}_0$ obstructed baseline)
   - `use_optimized_kernels`: `True` / `False`

---

## Section 3 — Initial Conditions, Complete Physical Horizons & Numerical Protocol (`Zero Pre-Pericenter Cutoff`)

1. **Benchmark Registry Factory & Exact Physical Parameters**:
   - Canonical factory in `benchmarks/registry.py` (`create_*_benchmark`), body count $N$, masses $\{m_i\}$, initial orbital elements $(a_i, e_i, i_i, \Omega_i, \omega_i, M_i)$, and softening parameter $\varepsilon = 0, d_0 = 0$.
2. **Complete Physical Integration Horizon (`t_span`, `dt_macro`, Step Counts)**:
   - Explicit mathematical justification that $t_{\mathrm{span}}$ covers the **entire physical phenomenon of interest** (e.g., $t_{\mathrm{span}} \ge T_{\mathrm{orb}}$ across pericenter, $t_{\mathrm{span}} > t_{\mathrm{coll}}$ across total/quadruple collision, full ingress $\to$ perijove/perigee $\to$ egress for planetary encounters, or integer resonant/secular periods).
   - **Strict Prohibition**: Truncating `n_steps` or `t_span` in `Float128` or `Float64` before the pericenter, Roche limit, or collision bounce is a contract violation.
3. **Step-Size Grid & Farrés–Laskar Normalization**:
   - Step sizes $h \in \{h_1, \dots, h_K\}$, force evaluations per step $s_{\mathrm{eval}} \in \{1, 3, 6, 10\}$, and normalized inverse-cost abscissa $h_{\mathrm{eff}} = h / s_{\mathrm{eval}}$.

---

## Section 4 — Exact Measurement Protocol & Prohibited Proxies (`Zero Synthetic Proxy`)

1. **Exact Observable Definitions**:
   - Relative energy error: $\Delta E_{\mathrm{rel}}(t) = |H(\mathbf{q}(t), \mathbf{p}(t)) - H_0| / |H_0|$ evaluated on the **physical state** $\mathbf{z}(t) = \mathcal{C}_h(\tilde{\mathbf{z}}(t))$ via `sys.observe_physical_state()`.
   - Position / phase-space convergence error: $\|\mathbf{q}(T) - \mathbf{q}_{\mathrm{ref}}^{\mathrm{f128}}(T)\|_2$ against a genuine `Float128` `BM6` reference trajectory.
   - 10 Galilean Lie–Bargmann generators: $\mathbf{P}_{\mathrm{tot}}, \mathbf{L}_{\mathrm{tot}}, \mathbf{G}_{\mathrm{boost}}, H$ via `compute_lie_bargmann_invariants`.
   - Astrodynamical observables: osculating elements $(a, e, i, \Omega, \omega, M, \lambda)$ via `compute_osculating_elements`, Laplace angle $\theta_{\mathrm{Lap}}(t) = \lambda_1 - 3\lambda_2 + 2\lambda_3$ via `compute_laplace_resonant_angle`, and Laskar NAFF fundamental frequency $\nu$ via `extract_laskar_fundamental_frequency` on $z(t) = a(t)e^{i\lambda(t)}$ and $e(t)e^{i\varpi(t)}$.
2. **Explicit Anti-Proxy Declarations**:
   - Zero 0-step static algebraic fixture substituted for a dynamical trajectory metric;
   - Zero synthetic phase constructed from energy errors `abs(dE_arr)`;
   - Zero `max(err, 1e-35)` clamping in persisted JSON metrics.

---

## Section 5 — Automated Quantitative Assertion Gate (`PASS/FAIL` Contract Verification)

Every generator script MUST execute a dedicated `verify_contract_assertions(data: dict) -> None` function before writing its JSON output file. If any assertion fails, the generator MUST raise an `AssertionError` and refuse to mark the campaign as certified:

| Assertion ID | Metric / JSON Key | Mathematical Meaning | Mandatory Threshold (`PASS` Condition) |
| :---: | :--- | :--- | :--- |
| `ASSERT-01` | `is_f128_verified` | True 35-digit `mpmath` execution in `Float128` sweeps | `True` (`floor <= 1e-30`) |
| `ASSERT-02` | `event_traversed` | Pericenter / Roche limit / collision bounce actually crossed | $r_{\min} \le r_{\mathrm{target}}$ and $t_{\mathrm{end}} > t_{\mathrm{event}}$ |
| `ASSERT-03` | `order_or_gain` | Asymptotic slope $p_{\mathrm{emp}}$ or shadow corrector gain $\mathcal{E}(\mathrm{BM}_4)/\mathcal{E}(\mathrm{PBM}_4)$ | Matches theoretical order/gain bound |
| `ASSERT-04` | `lie_bargmann_conservation` | Conservation of $\mathbf{P}_{\mathrm{tot}}, \mathbf{L}_{\mathrm{tot}}, \mathbf{G}_{\mathrm{boost}}$ | $\le 10^{-14}$ (`Float64`), $\le 10^{-32}$ (`Float128`) |
| `ASSERT-05` | `reversibility_palindrome` | Round-trip $\Phi_{-h}^N \circ \Phi_h^N = \mathrm{Id}$ & dyadic palindrome | `dyadic_palindrome_matched is True` |

---

## Section 6 — Panel-by-Panel Specification of the Three Segmented Plates (`Volet A`, `Volet B`, `Volet C`)

### 6.1. Volet A — Computational Performance, Scaling & Farrés–Laskar Inverse-Cost Plate (`*_plate_A_computational.pdf/.png`)
- **Subplot Grid & Dimensions**: `$1\times 2$` or `$2\times 2$`, `figsize=(13.5, 9.8)`, `dpi=300`, zero overlapping insets.
- **Panel (a)**: [Exact $x$-axis variable, $y$-axis variable, JSON keys plotted, reference slopes shown]
- **Panel (b)**: [Exact $x$-axis variable, $y$-axis variable, JSON keys plotted, unclipped empirical slopes]
- **Panel (c)/(d)**: [Dyadic level profile $k_{\mathrm{sym}}(t)$, Pareto work-precision frontier, or micro-kernel speedup]

### 6.2. Volet B — Algebraic Structure, Shadow Manifold, Gauge Cocycle, Lie–Bargmann & Laskar NAFF Plate (`*_plate_B_algebraic_invariants.pdf/.png`)
- **Subplot Grid & Dimensions**: `$1\times 2$` or `$2\times 2$`, `figsize=(13.5, 9.8)`, `dpi=300`, zero overlapping insets.
- **Panel (a)**: [Shadow manifold $\mathcal{Z}_4^{\mathrm{BM4}} \in \mathcal{P}_1$ vs $\mathcal{P}_1$ cotangent lift vs obstructed $\mathcal{P}_0$ vs bare $\mathrm{BM}_4$]
- **Panel (b)**: [Dynamic gauge transition $T_{A\to B}^{(\mathrm{PBM4})}$, Čech 1-cocycle closure, or round-trip palindrome]
- **Panel (c)/(d)**: [9 Lie–Bargmann generators across all integrators, Brouwer $\propto t^{1/2}$ law, or Laskar NAFF frequency drift $\delta\nu/\nu$]

### 6.3. Volet C — Physical & Astrodynamical Reality, McGehee Exceptional Divisor & Phase Portraits Plate (`*_plate_C_physical_dynamics.pdf/.png`)
- **Subplot Grid & Dimensions**: `$1\times 2$` or `$2\times 2$`, `figsize=(13.5, 9.8)`, `dpi=300`, zero overlapping insets.
- **Panel (a)**: [Physical orbital trajectory / McGehee $(\rho, v_{\mathrm{rad}})$ exceptional divisor projection across $d_0=0$]
- **Panel (b)**: [Osculating elements $(a(t), e(t), i(t))$, Roche-limit tidal elongation $\Delta r_{\mathrm{tid}}(t)$, or Laplace resonant libration $\theta_{\mathrm{Lap}}(t)$]

---

## Section 7 — Bijective Traceability Table (`JSON Key <-> Plot Panel <-> Manuscript Table / Text`)

| JSON File & Exact Key Path | Volet Plate & Panel ID | Target LaTeX Table / Section (`§6.6`–`§7.3`) | Certified Numerical Value (`Float64` / `Float128`) |
| :--- | :--- | :--- | :--- |
| `results/campaigns/generator_0X_*.json :: [key_path]` | `Volet A — Panel (a)` | `Table N, Row M` / `§7.X.Y` | `[Populated & verified upon generator execution]` |

