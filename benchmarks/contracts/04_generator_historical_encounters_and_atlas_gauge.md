# Benchmark Generator 04 Contract : Historical Planetary Encounters (`H1` Lexell, `H2` Shoemaker–Levy 9, `H3` Apophis 2029) & Dynamic Atlas Gauge Validation (`04_generator_historical_encounters_and_atlas_gauge.md`)

**Contract Governance Version**: Pillar II — Cycle V (Iteration 62, Arbitrages `A1`, `A2`, `A3`)  
**Reference Doctrine**: [`09_pillar2_campaign_audit_and_laskar_plate_doctrine.md`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/reviews/substantive_audit/09_pillar2_campaign_audit_and_laskar_plate_doctrine.md)  
**Canonical Template**: [`00_canonical_generator_contract_template.md`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scripts/campaigns/contracts/00_canonical_generator_contract_template.md)  
**Generator Script**: [`scripts/campaigns/run_generator_04_historical_encounters_and_atlas_gauge.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scripts/campaigns/run_generator_04_historical_encounters_and_atlas_gauge.py)  
**Plotter Script**: [`scripts/campaigns/plot_generator_04_historical_encounters_and_atlas_gauge.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scripts/campaigns/plot_generator_04_historical_encounters_and_atlas_gauge.py)  
**Primary JSON Output**: [`results/campaigns/generator_04_historical_encounters_and_atlas_gauge.json`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/results/campaigns/generator_04_historical_encounters_and_atlas_gauge.json)

---

## Section 1 — Mathematical Objective, Target Theorems & Formal Lean 4 Anchors

1. **Generator Identifier & Title**:
   - `Generator 04 — Documented Historical Planetary Encounters (H1 Lexell, H2 Shoemaker–Levy 9, H3 Apophis 2029) & Dynamic Atlas Gauge Validation (§7.2.1)`
2. **Target Manuscript Sections, Propositions & Tables**:
   - Primary LaTeX section: `§7.2.1` (`chapters/ch07_benchmarks_and_lean4/sec02_long_term_stability/content.tex`)
   - Target Table & Figures: `Table 7` (`Historical Planetary Encounter Benchmarks H1–H3`) and the 3 dedicated encounter plates (`fig_ch07_h1_lexell_atlas_gauge`, `fig_ch07_h2_shoemaker_levy_roche`, `fig_ch07_h3_apophis_flyby_palindrome`), plus the 3 segmented plates (`Volet A`, `Volet B`, `Volet C`)
   - Substantive Audit Findings Remediated (cross-referenced to Cycle IV [`05_sec7_2_1_historical_encounters_substantive_audit.md`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/reviews/substantive_audit/05_sec7_2_1_historical_encounters_substantive_audit.md)):
     - **`[CRIT-7.2.1-01]` & `[CRIT-7.2.1-05]` (`Protocol H1 — Comet D/1770 L1 Lexell`)**: Replaces the pre-perijove cutoff at $t_{\max} = 0.125$ ($r_{\min} = 115.25\,R_J$), the $r_{0,J} = 0.08\,\text{AU}$ start inside Jupiter's Hill sphere (`0` real atlas flips), the static 0-step 3-body toy gauge fixture, and the synthetic `sin`/`cos` curves in `Figure 6b` with a **complete 4-body integration starting outside Jupiter's Hill switching sphere at $r_{0,J} = 0.70\,\text{AU}$ (`> \gamma_{\mathrm{build}} R_{\mathrm{Hill},J} \approx 0.248\,\text{AU}$) in the heliocentric chart $\mathcal{T}_{\mathrm{helio}}$** under a **constant macro-step $\Delta t = 0.08$ ($N_{\mathrm{steps}} = 150$, $k_{\max} = 6$, native substep recording)**, dynamically crossing into the Jovicentric chart $\mathcal{T}_{\mathrm{Jovi}}$ at $t_{\mathrm{flip},1} = 2.48$, traversing the deep perijove at **$r_{\min} = 1.4788\times 10^{-3}\,\text{AU} = 3.094\,R_J$** ($t_{\mathrm{peri}} = 5.45875$, inside Io's orbit $a_{\mathrm{Io}} = 5.90\,R_J$, with unimodal symmetric dyadic pyramid $k_{\mathrm{sym}} \in \{0,\dots,6\}$), and flipping back to $\mathcal{T}_{\mathrm{helio}}$ on egress at $t_{\mathrm{flip},2} = 11.28$ (`n_atlas_flips = 2`), while verifying seamless interfacial Hamiltonian and Lie–Bargmann continuity under the canonical covariant gauge operator $T_{(A,h_A)\to(B,h_B)}^{(\mathrm{PBM4})}$ (`mode="homological"`, zero ablation curve since chart-flip gauge ablation is already treated in `§7.1.2`, Protocol `ST-4`).
     - **`[CRIT-7.2.1-02]`, `[CRIT-7.2.1-05]` & `[MAJ-7.2.1-01]` (`Protocol H2 — Comet Shoemaker–Levy 9, 1992`)**: Replaces the pre-perijove cutoff at `n_steps = 135` ($t_{\max} = 0.02025$, stopping 11 steps before perijove at $3.76\,R_J > r_{\mathrm{Roche}} = 2.45\,R_J$) and the unconditional $t = 0$ supernode activation via $\delta_{\mathrm{nat}}$ with a **complete sub-Roche passage ($n_{\mathrm{steps}} = 160, \Delta t = 0.00018, t_{\max} = 0.0288 > t_{\mathrm{peri}} = 0.02178$)** reaching **$r_{\min} = 6.384\times 10^{-4}\,\text{AU} = 1.336\,R_J < r_{\mathrm{Roche}} = 2.45\,R_J$**, governed by the dimensionless cubic tidal force ratio $\Lambda_S(t) = \frac{1}{3}\eta_S(t)^3 \ge \Lambda_{\mathrm{crit}} = 0.25$ (`supernode_mode="roche_disruption"`), dynamically promoting the 3-fragment progenitor $\{2, 3, 4\}$ inside the Roche limit ($\Lambda_S(0) = 0.0272 < 0.25 \to \Lambda_S^{\max} > 10^4$) and tracking the physical "string of pearls" tidal elongation $S(t)/S_0 \approx 3.82\times$.
     - **`[CRIT-7.2.1-03]`, `[CRIT-7.2.1-04]`, `[MAJ-7.2.1-03]` & `[MOD-7.2.1-01]` (`Protocol H3 — Asteroid (99942) Apophis, April 13, 2029`)**: Replaces the inverted Apollo ($1.123\,\text{AU}$) $\to$ Aten ($0.927\,\text{AU}$) trajectory and the 10-step pre-flyby reversibility check with the **authentic trailing-side gravity assist** deflecting Apophis outward from its pre-flyby **Aten orbit ($a_{\mathrm{pre}} = 0.9224\,\text{AU} < 1\,\text{AU}$)** across perigee **$r_{\min} = 2.548\times 10^{-4}\,\text{AU} = 5.976\,R_\oplus$** (inside the geostationary belt $6.61\,R_\oplus$) into its post-flyby **Apollo orbit ($a_{\mathrm{post}} = 1.1163\,\text{AU} > 1\,\text{AU}$, $\Delta a = +0.1939\,\text{AU} > 0$)**, and executes the forward–backward reversibility test $\Phi_{-\Delta t}^N \circ \Phi_{\Delta t}^N$ across the **entire perigee flyby ($N = 80, \Delta t = 0.00025, t_{\max} = 0.020$)** with exact dyadic palindromic symmetry (`dyadic_palindrome_matched == True`, $k_{\mathrm{fwd}}(n) \equiv k_{\mathrm{bwd}}(N-1-n)$) in both `Float64` and `Float128`.
3. **Mathematical & Theoretical Propositions Tested**:
   - **Proposition 6.13 & Theorem 5.17–5.18 (Two-Tier Hill Hysteresis & Non-Commutative Covariant Gauge Operator)**: Chart transition $\mathcal{T}_A \to \mathcal{T}_B$ triggered at $\eta_{\mathrm{Hill}} < \gamma_{\mathrm{monitor}} = 0.60$ onto a target chart satisfying $\eta_{\mathrm{Hill}}^{\mathrm{new}} \ge \gamma_{\mathrm{build}} = 0.70$, accompanied by the non-commutative cotangent sandwich
     $$T_{(A,h_A)\to(B,h_B)}^{(\mathrm{PBM4})} = \big(\mathcal{Z}_4^{\mathrm{BM4},(B)}(h_B)\big)^{-1} \circ \mathbf{M}_{A\to B} \circ \mathcal{Z}_4^{\mathrm{BM4},(A)}(h_A) \in \operatorname{Symp}(T^*M),$$
     which eliminates the $\mathcal{O}(h^4)$ interfacial shadow energy discontinuity across chart and dyadic scale boundaries.
   - **Section 7.2.1 Equation (7.34) (Cubic Tidal Force Ratio & Roche-Limit Supernode Promotion)**:
     $$\Lambda_S(\mathbf{q}) = \frac{1}{3}\eta_S(\mathbf{q})^3 = \frac{M_{\mathrm{ext}}}{3 M_S}\left(\frac{\bar{\rho}_S}{\|\mathbf{r}_P\|}\right)^3 \ge \Lambda_{\mathrm{crit}} = 0.25,$$
     governing the onset of tidal disruption inside the fluid Roche limit $r_{\mathrm{Roche}} \approx 2.45\,R_J$.
   - **Theorem 5.20, Definition 6.2 & Proposition 6.14 (Symmetric Dyadic Level Selector & Bilateral Pericenter Window)**:
     $$k_{\mathrm{sym}}(\mathbf{z}_n, \Delta t) = \max\!\big(k(\mathbf{z}_n),\, k(\Phi_{\Delta t}^{(0)}(\mathbf{z}_n))\big), \qquad |t_{\mathrm{peri}}| = \frac{|\mathbf{r}_{ij}\cdot\mathbf{v}_{ij}|}{\|\mathbf{v}_{ij}\|^2} \le 2|\Delta t|,$$
     ensuring exact time-reversal palindromic symmetry $k_{\mathrm{fwd}}(n) = k_{\mathrm{bwd}}(N-1-n)$ across hyperbolic planetary flybys.
4. **Formal Lean 4 Proof Modules Linked (`Appendix A`)**:
   - [`02_JacobiTree_TransitionMatrices.lean`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch04_fulton_macpherson_splitting/sec03_synchronous_supernodes/02_JacobiTree_TransitionMatrices.lean) (`jacobi_transition_matrix_in_sl2`, `cotangent_lift_preserves_darboux`);
   - [`03_IntrinsicCompactness_Criteria.lean`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch04_fulton_macpherson_splitting/sec03_synchronous_supernodes/03_IntrinsicCompactness_Criteria.lean) (`hill_ratio_scale_invariant`);
   - [`06_ConormalGaugeTransitionPBM4_Order6.lean`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/06_ConormalGaugeTransitionPBM4_Order6.lean) (`pbm4_gauge_operator_conjugation_cocycle`, `pbm4_gauge_transition_order6`).

---

## Section 2 — Paper Showcase Interfaces (`01..06_*.py`) & Dual-Precision Library Bindings

1. **Mandatory Paper Showcase Imports (`Float64`)**:
   - `01_GlobalPipeline_Interface.py` (`Appendix A.22`): `ExtendedPhaseState`, `NBodySystem`
   - `02_TreeBuilderAndCorrector_Interface.py` (`Appendix A.23`): `detect_supernodes`, `PalindromicCorrectorZ4`, `GlobalCorrectorC4`, `create_pbm4_integrator`
   - `03_SymplecticSubcycle_Interface.py` (`Appendix A.24`): `DyadicScheduler`
   - `04_RegularizedSolvers_Interface.py` (`Appendix A.25`): `BlanesMoan4Integrator`
   - `05_AtlasTransitions_Interface.py` (`Appendix A.26`): `AtlasManager`, `PBM4GaugeTransitionOperator`, `apply_gauge_connection`, `compute_lie_bargmann_invariants`
2. **Quadruple-Precision (`Float128`, 35 Digits `mpmath`) Mirror Bindings**:
   - `normal_cone_nbody_float128.benchmarks.registry.get_benchmark`
   - `normal_cone_nbody_float128.system.nbody_system.NBodySystem`
   - `normal_cone_nbody_float128.integrators.processed.create_pbm4_integrator`, `apply_gauge_connection`
   - `normal_cone_nbody_float128.potential.supernode_detector.detect_supernodes`, `compute_cubic_tidal_force_ratio`
   - `normal_cone_nbody_float128.diagnostics.metrics.compute_osculating_elements`, `compute_forward_backward_reversibility`
   - Mandatory runtime verification that `state_f128.is_f128 is True` and `isinstance(state_f128.q[0, 0], mpmath.mpf)` across all three protocols (`H1`, `H2`, `H3`).
3. **Strict Parameter Equality Rule (`Float64` vs `Float128`)**:
   - In every protocol (`H1`, `H2`, `H3`), the `Float64` and `Float128` simulations MUST use **100% identical numerical and physical parameters** (`dt`, `n_steps`, `max_dyadic_level`, `enable_dyadic_substepping`, `in_shadow_mode`, `mode`, `supernode_mode`).

---

## Section 3 — Initial Conditions, Complete Physical Horizons & Numerical Protocol (`Zero Pre-Pericenter Cutoff`)

1. **Protocol `H1` — Comet D/1770 L1 (Lexell) Deep Perijove Ejection & Dynamic Atlas Gauge Validation (`--protocol H1`)**:
   - **Benchmark Factory**: `get_benchmark("comet_lexell_1770")` ($N = 4$: Sun $m_0 = 1.0$, Jupiter $m_1 = 9.5458\times 10^{-4}$ at $a_J = 5.20336\,\text{AU}$, Earth $m_2 = 3.0035\times 10^{-6}$ at $a_\oplus = 1.0\,\text{AU}$, Comet Lexell $m_3 = 1.0\times 10^{-14}$ starting outside Jupiter's Hill switching boundary at $r_{0,J} = 0.70\,\text{AU}$ with Coriolis-compensated aiming angle $\sin\gamma = -0.005767603157$).
   - **Strict Equal-Parameter Horizon (`Float64` & `Float128`)**:
     - Integrated from $t = 0.0$ to $t_{\max} = 12.0$ (`~ 1.91 yr`) using an identical uniform macro-step $\Delta t = 0.08$ (`n_steps = 150` uniform macro-steps) in `Float64` and `Float128`, with `max_dyadic_level = 6`, `record_substeps = True`, `substep_record_min_level = 4`, `in_shadow_mode = True`, `mode = "homological"`.
     - Traverses:
       1. Ingress chart flip $\mathcal{T}_{\mathrm{helio}} = [((0, 2), 1), 3] \to \mathcal{T}_{\mathrm{Jovi}} = [(0, 2), (1, 3)]$ at $t_{\mathrm{flip},1} = 2.48$ ($r_J \approx 0.2105\,\text{AU} = 0.5929\,R_{\mathrm{Hill},J} < \gamma_{\mathrm{monitor}} = 0.60$);
       2. Deep perijove passage at $t_{\mathrm{peri}} = 5.45875$ with **$r_{\min} = 1.4788\times 10^{-3}\,\text{AU} = 3.094\,R_J$** (inside Io's orbit $5.90\,R_J$, unimodal symmetric dyadic pyramid $k_{\mathrm{sym}} : 1 \to 2 \to 3 \to 4 \to 6 \to 4 \to 3 \to 2 \to 1$);
       3. Egress chart flip $\mathcal{T}_{\mathrm{Jovi}} \to \mathcal{T}_{\mathrm{helio}}$ at $t_{\mathrm{flip},2} = 11.28$ (`dwell_time = 8.80`).
   - **Dynamic Gauge Continuity Verification on the Actual Lexell Encounter**:
     - At the ingress Hill chart flip $t_{\mathrm{flip},1} = 2.48$ on the actual 4-body Comet Lexell state in both `Float64` and `Float128`, verifies that the canonical covariant gauge operator $T_{(A,h_A)\to(B,h_B)}^{(\mathrm{PBM4})}$ (`mode="homological"`, $\mathcal{Z}_4^{\mathrm{BM4}} \in \mathcal{P}_1$) preserves both the unweighted 4-body Hamiltonian interfacial jump $|\Delta H_{\mathrm{jump}}| / |E_0| = 2.07\times 10^{-35}$ in `Float128` and the mass-undiluted Comet-sector specific energy jump $|\Delta \mathcal{E}_3| / |\mathcal{E}_{3,0}| = 2.06\times 10^{-16}$ in `Float64` and $1.82\times 10^{-35}$ in `Float128`.

2. **Protocol `H2` — Comet Shoemaker–Levy 9 (1992) Sub-Roche Tidal Disruption & Dynamic Supernode Promotion (`--protocol H2`)**:
   - **Benchmark Factory**: `get_benchmark("shoemaker_levy_9_1992")` ($N = 5$: Sun $m_0 = 1.0$, Jupiter $m_1 = 9.5458\times 10^{-4}$, and $K = 3$ self-gravitating progenitor fragments $m_2 = m_3 = m_4 = 1.0\times 10^{-11}$ separated by $\delta r_0 = 2.0\times 10^{-5}\,\text{AU}$, starting at $r_{0,J} = 0.012\,\text{AU} = 25.11\,R_J$).
   - **Strict Equal-Parameter Horizon (`Float64` & `Float128`)**:
     - Integrated with `dt = 0.00018`, `n_steps = 160` ($t_{\max} = 0.0288 > t_{\mathrm{peri}} = 0.02178$), `max_dyadic_level = 3`, `in_shadow_mode = True`, `mode = "homological"`, `supernode_mode = "roche_disruption"` identically in `Float64` and `Float128`.
     - Traverses the full ingress across Jupiter's fluid Roche limit ($r_{\mathrm{Roche}} = 2.45\,R_J = 1.1708\times 10^{-3}\,\text{AU}$), deep perijove at **$r_{\min} = 6.384\times 10^{-4}\,\text{AU} = 1.336\,R_J < 2.45\,R_J$**, and egress back outside the Roche limit, recording at every step: $r_J(t)$, the cubic tidal force ratio $\Lambda_S(t) = \frac{1}{3}\eta_S(t)^3$, the linear Hill ratio $\eta_S(t)$, the boolean supernode promotion flag (`len(active_supernodes) > 0`), the "string of pearls" fragment chain span $S(t)/S_0 = \|\mathbf{q}_4(t) - \mathbf{q}_2(t)\|_2 / S_0$, and the relative energy error $|\Delta E(t)| / |E_0|$.

3. **Protocol `H3` — Asteroid (99942) Apophis (April 13, 2029) Perigee Flyby, Aten $\to$ Apollo Transition & Full-Flyby Palindrome (`--protocol H3`)**:
   - **Benchmark Factory**: `get_benchmark("apophis_2029_flyby")` ($N = 4$: Sun $m_0 = 1.0$, Earth $m_1 = 3.0035\times 10^{-6}$, Moon $m_2 = 3.6943\times 10^{-8}$ at $a_{\mathrm{Moon}} = 0.0025695\,\text{AU}$, Asteroid Apophis $m_3 = 4.0\times 10^{-19}$ starting at $r_{0,\oplus} = 0.0020\,\text{AU}$ on a trailing-side hyperbolic approach $\phi_{\mathrm{app}} = 171.575217^\circ$).
   - **Strict Equal-Parameter Horizon (`Float64` & `Float128`)**:
     - Integrated with `dt = 0.00025`, `n_steps = 80` ($t_{\max} = 0.0200 > t_{\mathrm{peri}} = 0.00950$), `max_dyadic_level = 4`, `in_shadow_mode = True`, `mode = "homological"` identically in `Float64` and `Float128`.
     - Traverses the grazing perigee passage at **$r_{\min} = 2.548\times 10^{-4}\,\text{AU} = 5.976\,R_\oplus$** (inside the geostationary belt $r_{\mathrm{GEO}} = 6.61\,R_\oplus$), recording the heliocentric osculating orbital elements $(a(t), e(t), i(t))$ of Apophis across the resonant transition from **Aten class ($a_{\mathrm{pre}} = 0.9224\,\text{AU} < 1\,\text{AU}$)** to **Apollo class ($a_{\mathrm{post}} = 1.1163\,\text{AU} > 1\,\text{AU}$, $\Delta a = +0.1939\,\text{AU} > 0$)**, the bilateral pericenter time $|t_{\mathrm{peri}}(t)| = |\mathbf{r}_{13}\cdot\mathbf{v}_{13}| / \|\mathbf{v}_{13}\|^2$, and the symmetric dyadic level profile $k_{\mathrm{sym}}(t)$.
     - Executes `compute_forward_backward_reversibility` across the **entire 80-step perigee flyby ($t: 0 \to 0.0200 \to 0$)** in both `Float64` and `Float128` with identical parameters (`dt = 0.00025, n_steps = 80`), verifying `dyadic_palindrome_matched == True` ($k_{\mathrm{fwd}}(n) \equiv k_{\mathrm{bwd}}(N-1-n)$ with $\max k = 3$ at perigee).

---

## Section 4 — Exact Measurement Protocol & Prohibited Proxies (`Zero Synthetic Proxy`)

1. **Exact Observable Definitions**:
   - **Physical State Readout**: At every step, physical coordinates are extracted via `sys.observe_physical_state()` ($\mathbf{z}_n = \mathcal{C}_h(\tilde{\mathbf{z}}_n)$).
   - **Comet-Sector Specific Orbital Energy (`H1` Lexell)**:
     $$\mathcal{E}_3(\mathbf{z}) = \frac{\|\mathbf{p}_3\|^2}{2 m_3^2} - \sum_{a=0}^{2} \frac{G m_a}{\|\mathbf{q}_3 - \mathbf{q}_a\|_2},$$
     which measures the comet's energy continuity across $\mathcal{T}_{\mathrm{helio}} \leftrightarrow \mathcal{T}_{\mathrm{Jovi}}$ without $m_3 / m_1 \sim 10^{-11}$ mass-ratio dilution.
   - **Cubic Tidal Force Ratio (`H2` Shoemaker–Levy 9)**:
     $$\Lambda_S(t) = \frac{M_{\mathrm{Jup}}}{3 M_S}\left(\frac{\bar{\rho}_S(t)}{\|\mathbf{r}_{\mathrm{cdm},S}(t) - \mathbf{q}_{\mathrm{Jup}}(t)\|_2}\right)^3 = \frac{1}{3}\eta_S(t)^3.$$
   - **Heliocentric Osculating Elements (`H3` Apophis)**: Computed at every step from `compute_osculating_elements(state, primary_idx=0, body_idx=3, G=1.0)`.
2. **Explicit Anti-Proxy Declarations**:
   - Zero synthetic `sin`/`cos` formula for energy error curves;
   - Zero static 3-body toy system substituted for `H1` Lexell's chart flip;
   - Zero pre-perijove / pre-perigee truncation (`H1`, `H2`, `H3` all cross $r_{\min}$ and egress);
   - Zero mismatch of integration parameters between `Float64` and `Float128`.

---

## Section 5 — Automated Quantitative Assertion Gate (`PASS/FAIL` Contract Verification)

| Assertion ID | Metric / JSON Key | Mathematical Meaning | Mandatory Threshold (`PASS` Condition) |
| :---: | :--- | :--- | :--- |
| `ASSERT-01` | `is_f128_verified` | True 35-digit `mpmath` execution in `Float128` across `H1`, `H2`, `H3` with identical parameters to `Float64` | `True` and `equal_params_verified is True` |
| `ASSERT-02` | `H1_lexell_atlas_and_perijove` | Comet Lexell starts at $r_{0,J}=0.70\,\text{AU}$, triggers $\ge 2$ real chart flips ($\mathcal{T}_{\mathrm{helio}} \to \mathcal{T}_{\mathrm{Jovi}} \to \mathcal{T}_{\mathrm{helio}}$), and crosses deep perijove $r_{\min} \approx 3.10\,R_J$ | `n_atlas_flips >= 2` and $3.00 \le r_{\min}/R_J \le 3.20$ |
| `ASSERT-03` | `H1_lexell_gauge_continuity` | Covariant gauge operator $T_{(A,h_A)\to(B,h_B)}^{(\mathrm{PBM4})}$ preserves interfacial energy continuity across the real Lexell Hill flip to machine precision | `jump_comet_dE3_with_gauge_f64 <= 1e-13`, `jump_comet_dE3_with_gauge_f128 <= 1e-32`, `jump_total_dH_with_gauge_f128 <= 1e-33` |
| `ASSERT-04` | `H2_sl9_sub_roche_and_promotion` | SL9 crosses deep inside Jupiter's fluid Roche limit ($r_{\min} < 1.50\,R_J < 2.45\,R_J$), dynamically promotes `{2,3,4}` when $\Lambda_S(t) \ge 0.25$ (inactive at $t=0$), and stretches into a string of pearls | $r_{\min}/R_J \le 1.50$, $\Lambda_S(0) < 0.25$, $\Lambda_S^{\max} > 100.0$, `0 < roche_active_steps < n_steps`, and `max_span_ratio >= 3.0` |
| `ASSERT-05` | `H3_apophis_aten_to_apollo` | Apophis approaches on an Aten orbit ($a_{\mathrm{pre}} < 1\,\text{AU}$), grazes perigee at $r_{\min} \approx 5.98\,R_\oplus < 6.61\,R_\oplus$, and is deflected outward to an Apollo orbit ($a_{\mathrm{post}} > 1\,\text{AU}, \Delta a > +0.15\,\text{AU}$) | $a_{\mathrm{pre}} < 0.95$, $a_{\mathrm{post}} > 1.08$, $\Delta a > 0.15$, and $5.80 \le r_{\min}/R_\oplus \le 6.15$ |
| `ASSERT-06` | `H3_apophis_full_flyby_palindrome` | Forward–backward reversibility across the entire 80-step perigee flyby preserves the exact dyadic palindrome $k_{\mathrm{fwd}}(n) \equiv k_{\mathrm{bwd}}(N-1-n)$ with $\max k \ge 2$ | `dyadic_palindrome_matched_f64 is True`, `dyadic_palindrome_matched_f128 is True`, `max(k_fwd) >= 2`, `dE_fb_f64 <= 1e-11`, `dE_fb_f128 <= 1e-30` |

---

## Section 6 — Panel-by-Panel Specification of the Three Segmented Plates (`Volet A`, `Volet B`, `Volet C`)

### 6.1. Volet A — Computational Performance, Dyadic Multirate Activity & Roche/Atlas Triggers Plate (`gen04_plate_A_computational.pdf/.png`)
- **Subplot Grid & Dimensions**: `$2\times 2$`, `figsize=(14.0, 10.2)`, `dpi=300`, zero overlapping insets.
- **Panel (a) [`H1` Lexell — Hill Separation Ratio, Two-Tier Hysteresis & Atlas Chart Flips]**:
  - Plots the normalized Jovicentric Hill separation ratio $r_J(t) / R_{\mathrm{Hill},J}$ (`Float64` and `Float128` superposed) across $t \in [0, 12.0]$, showing the two-tier hysteresis thresholds $\gamma_{\mathrm{build}} = 0.70$ and $\gamma_{\mathrm{monitor}} = 0.60$, the shaded Jovicentric chart dwell window $\mathcal{T}_{\mathrm{Jovi}}$ ($t \in [2.48, 11.30]$, `dwell_time = 8.82`), and the dyadic refinement level $k_{\mathrm{sym}}(t)$ peaking at perijove ($t = 5.458$).
- **Panel (b) [`H2` Shoemaker–Levy 9 — Cubic Tidal Force Ratio $\Lambda_S(t)$ & Dynamic Sub-Roche Supernode Promotion]**:
  - Plots $\Lambda_S(t) = \frac{1}{3}\eta_S(t)^3$ and $\eta_S(t)$ on a logarithmic scale across $t \in [0, 0.0288]$, highlighting the Roche disruption threshold $\Lambda_{\mathrm{crit}} = 0.25$ (and Schmitt hysteresis band $[\Lambda_{\mathrm{exit}}, \Lambda_{\mathrm{enter}}] = [0.15, 0.20]$) and the shaded window where `detect_supernodes(..., mode="roche_disruption")` dynamically promotes `{2, 3, 4}`.
- **Panel (c) [`H3` Apophis 2029 — Symmetric Dyadic Level Profile $k_{\mathrm{sym}}(t)$ & Bilateral Pericenter Window]**:
  - Plots the executed dyadic level $k_{\mathrm{sym}}(t) \in \{0, 1, 2, 3\}$ and the normalized bilateral pericenter time $|t_{\mathrm{peri}}(t)| / \Delta t_{\mathrm{macro}}$ across the perigee passage $t \in [0, 0.0200]$, showing the exact $2|\Delta t_{\mathrm{macro}}|$ anticipation threshold.
- **Panel (d) [Dual-Precision (`Float64` vs `Float128`) Trajectory Hamiltonian Conservation Under Strictly Identical Parameters]**:
  - Compares the relative Hamiltonian error $|\Delta E(t)| / |E_0|$ vs normalized encounter time $t / t_{\max}$ for `H1` Lexell, `H2` Shoemaker–Levy 9, and `H3` Apophis in both `Float64` (solid lines) and `Float128` (dashed lines) under 100% identical integration parameters.

### 6.2. Volet B — Algebraic Structure, Dynamic Atlas Gauge Transition, Palindromic Reversibility & Lie–Bargmann Invariants Plate (`gen04_plate_B_algebraic_invariants.pdf/.png`)
- **Subplot Grid & Dimensions**: `$2\times 2$`, `figsize=(14.0, 10.2)`, `dpi=300`, zero overlapping insets.
- **Panel (a) [`H1` Lexell — Authentic Dynamic Gauge Ablation Across the Hill Chart Flip $\mathcal{T}_{\mathrm{helio}} \to \mathcal{T}_{\mathrm{Jovi}}$]**:
  - Plots the true post-flip interfacial energy discontinuity across the ingress Hill chart flip ($t \in [1.68, 3.28]$, flip at $t_{\mathrm{flip},1} = 2.48$) on the actual 4-body Comet Lexell encounter for: (i) Covariant Gauge Operator $T_{(A,h_A)\to(B,h_B)}^{(\mathrm{PBM4})}$ ($\mathcal{Z}_4^{\mathrm{BM4}} \in \mathcal{P}_1$), (ii) Naive unconnected chart switch, and (iii) Obstructed unipotent $\mathcal{P}_0$ kick (`mode="global"`), in both `Float64` (Comet-sector $|\Delta\mathcal{E}_3|/|\mathcal{E}_{3,0}|$) and `Float128` (both $|\Delta\mathcal{E}_3|/|\mathcal{E}_{3,0}|$ and unweighted 4-body $|\Delta H|/|E_0|$ reaching $2.07\times 10^{-35}$).
- **Panel (b) [`H3` Apophis 2029 — Exact Forward–Backward Dyadic Palindrome Across Full Perigee Flyby]**:
  - Plots the forward dyadic level sequence $k_{\mathrm{fwd}}(n)$ and time-reversed backward sequence $k_{\mathrm{bwd}}(N-1-n)$ ($N = 80$ steps across $r_{\min} = 5.98\,R_\oplus$), demonstrating exact point-by-point palindromic coincidence $k_{\mathrm{fwd}}(n) - k_{\mathrm{bwd}}(N-1-n) \equiv 0$.
- **Panel (c) [`H3` Apophis 2029 — Dual-Precision Round-Trip Reversibility Error Across Perigee]**:
  - Plots the forward and backward return trajectory errors $\|\mathbf{q}_{\mathrm{FB}}(t) - \mathbf{q}_{\mathrm{fwd}}(t)\|_2$ and $|\Delta E_{\mathrm{FB}}|/|E_0|$ across the perigee passage in `Float64` ($\sim 10^{-13}$) and `Float128` ($\sim 10^{-33}$).
- **Panel (d) [Galilean Lie–Bargmann Conservation ($\|\Delta\mathbf{P}\|_2, \|\Delta\mathbf{L}\|_2, \|\Delta\mathbf{G}\|_2$) Across `H1`, `H2`, `H3` in `Float64` & `Float128`]**:
  - Grouped bar chart comparing the maximum drift of total linear momentum $\mathbf{P}_{\mathrm{tot}}$, angular momentum $\mathbf{L}_{\mathrm{tot}}$, and center-of-mass boost $\mathbf{G}_{\mathrm{boost}}$ in `Float64` ($\le 10^{-13}$) vs `Float128` ($\le 10^{-33}$) across all three historical encounters.

### 6.3. Volet C — Physical & Astrodynamical Reality Plate (`gen04_plate_C_physical_dynamics.pdf/.png`)
- **Subplot Grid & Dimensions**: `$2\times 2$`, `figsize=(14.0, 10.2)`, `dpi=300`, zero overlapping insets.
- **Panel (a) [`H1` Comet Lexell — Full Jovicentric Approach, Deep Perijove Passage ($3.10\,R_J$) & Egress]**:
  - Plots the Jovicentric distance $r_J(t)$ (in $\text{AU}$ and $R_J$ on a logarithmic axis) from $r_{0,J} = 0.70\,\text{AU}$ ($1{,}465\,R_J$, outside the Hill sphere $R_{\mathrm{Hill},J} = 743\,R_J$) down through Io's orbit ($5.90\,R_J$) to the deep perijove floor **$r_{\min} = 1.4815\times 10^{-3}\,\text{AU} = 3.100\,R_J$** at $t = 5.458$ and back out to $0.79\,\text{AU}$ at $t = 12.0$.
- **Panel (b) [`H2` Comet Shoemaker–Levy 9 — Sub-Roche Perijove Traversal ($1.34\,R_J < 2.45\,R_J$) & String-of-Pearls Elongation]**:
  - Dual-axis plot showing the Jovicentric distance $r_J(t)/R_J$ plunging deep inside Jupiter's fluid Roche limit $r_{\mathrm{Roche}} = 2.45\,R_J$ to **$r_{\min} = 1.336\,R_J$** ($6.384\times 10^{-4}\,\text{AU}$) at $t = 0.02178$ and rebounding, alongside the differential tidal stretching $S(t)/S_0 = \|\mathbf{q}_4(t)-\mathbf{q}_2(t)\|_2 / S_0$ reaching $3.82\times$.
- **Panel (c) [`H3` Asteroid Apophis 2029 — Grazing Geocentric Hyperbolic Flyby Inside the Geostationary Ring]**:
  - Plots the 2D geocentric trajectory $(x - x_\oplus, y - y_\oplus)$ in Earth radii $R_\oplus$, showing the Earth disk ($1.0\,R_\oplus$), the geostationary satellite ring ($6.61\,R_\oplus$), and Apophis's hyperbolic deflection passing through perigee at **$r_{\min} = 5.976\,R_\oplus$** ($38{,}074\,\text{km}$).
- **Panel (d) [`H3` Asteroid Apophis 2029 — Resonant Heliocentric Orbital Transition: Aten ($a_{\mathrm{pre}} < 1\,\text{AU}$) $\to$ Apollo ($a_{\mathrm{post}} > 1\,\text{AU}$)]**:
  - Plots the heliocentric osculating semi-major axis $a(t)$ (left axis) and eccentricity $e(t)$ (right axis) across the April 13, 2029 encounter, showing the outward jump from **Aten class ($a_{\mathrm{pre}} = 0.9224\,\text{AU} < 1.0\,\text{AU}$)** across the $a = 1.0\,\text{AU}$ Earth-orbit boundary to **Apollo class ($a_{\mathrm{post}} = 1.1163\,\text{AU} > 1.0\,\text{AU}$, $\Delta a = +0.1939\,\text{AU}$)**.

---

## Section 7 — Bijective Traceability Table (`JSON Key <-> Plot Panel <-> Manuscript Table / Text`)

| JSON File & Exact Key Path | Volet Plate & Panel ID | Target LaTeX Table / Section (`§7.2.1`) | Certified Numerical Value (`Float64` / `Float128`) |
| :--- | :--- | :--- | :--- |
| `generator_04_historical_encounters_and_atlas_gauge.json :: H1_lexell.min_r_jup_au` / `min_r_jup_RJ` | `Plate H1 — Panel (a)`, `Volet C — Panel (a)` | `§7.2.1` (`H1` Lexell), `Table 7, Row 1` | $r_{\min} = 1.478758\times 10^{-3}\,\text{AU} = 3.094\,R_J$ at $t = 5.45875$ (`Float64` $\max|\Delta H|/|H_0| = 8.29\times 10^{-15}$, `Float128` $= 7.13\times 10^{-17}$, $k_{\max} = 6$) |
| `generator_04_historical_encounters_and_atlas_gauge.json :: H1_lexell.n_atlas_flips` / `dwell_time` | `Plate H1 — Panel (a)`, `Volet A — Panel (a)` | `§7.2.1` (`H1` Lexell), `Table 7, Row 1` | `n_atlas_flips = 2` ($\mathcal{T}_{\mathrm{helio}}\to\mathcal{T}_{\mathrm{Jovi}}$ at $t=2.48$, $\mathcal{T}_{\mathrm{Jovi}}\to\mathcal{T}_{\mathrm{helio}}$ at $t=11.28$, `dwell_time = 8.80`) |
| `generator_04_historical_encounters_and_atlas_gauge.json :: H1_lexell.gauge_continuity` | `Plate H1 — Panel (b)`, `Volet B — Panel (a)` | `§7.2.1` (Theorem 5.17–5.18 continuity), `Table 7, Row 1` | `Float64`: $|\Delta\mathcal{E}_3|/|\mathcal{E}_{3,0}| = 2.06\times 10^{-16}$; `Float128`: $|\Delta\mathcal{E}_3|/|\mathcal{E}_{3,0}| = 1.82\times 10^{-35}$, $|\Delta H_{\mathrm{flip}}|/|H_0| = 2.07\times 10^{-35}$ |
| `generator_04_historical_encounters_and_atlas_gauge.json :: H2_shoemaker_levy.min_r_jup_au` / `min_r_jup_RJ` | `Plate H2 — Panel (a, c)`, `Volet C — Panel (b)` | `§7.2.1` (`H2` SL9), `Table 7, Row 2` | $r_{\min} = 6.383939\times 10^{-4}\,\text{AU} = 1.336\,R_J < r_{\mathrm{Roche}} = 2.45\,R_J$ at $t = 0.02178$ (`Float64` $5.61\times 10^{-15}$, `Float128` $4.40\times 10^{-15}$) |
| `generator_04_historical_encounters_and_atlas_gauge.json :: H2_shoemaker_levy.lambda_max` / `max_span_ratio` | `Plate H2 — Panel (a, b)`, `Volet A — Panel (b)` | `§7.2.1` (`H2` SL9, Eq. 7.34), `Table 7, Row 2` | $\Lambda_S(0) = 0.0267 < 0.25$ (inactive at $t=0$), $\Lambda_S^{\max} = 16{,}689.7$ ($\eta_S^{\max} = 36.86$), promoted over `87/160` steps, chain elongation $S_{\max}/S_0 = 4.61\times$ ($S_f/S_0 = 3.82\times$) |
| `generator_04_historical_encounters_and_atlas_gauge.json :: H3_apophis.a_pre_au` / `a_post_au` / `delta_a_au` | `Plate H3 — Panel (a, b)`, `Volet C — Panel (c, d)` | `§7.2.1` (`H3` Apophis), `Table 7, Row 3` | $r_{\min} = 2.547963\times 10^{-4}\,\text{AU} = 5.976\,R_\oplus < 6.61\,R_\oplus$; Aten $a_{\mathrm{pre}} = 0.9224\,\text{AU} \to$ Apollo $a_{\mathrm{post}} = 1.1163\,\text{AU}$ ($\Delta a = +0.1939\,\text{AU} > 0$) |
| `generator_04_historical_encounters_and_atlas_gauge.json :: H3_apophis.reversibility` | `Plate H3 — Panel (c, d)`, `Volet B — Panel (b, c)` | `§7.2.1` (`H3` Apophis), `Table 7, Row 3` | `dyadic_palindrome_matched = True` ($N=80, k_{\max}=3$) in `Float64` & `Float128`; $\|\Delta\mathbf{q}_{\mathrm{FB}}\|_2 = 1.72\times 10^{-16}$ (`Float64`) / $3.28\times 10^{-33}$ (`Float128`); $|\Delta H_{\mathrm{FB}}|/|H_0| = 0.0$ (`Float64`) / $8.46\times 10^{-33}$ (`Float128`) |
