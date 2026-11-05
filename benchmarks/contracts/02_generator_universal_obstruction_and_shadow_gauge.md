# Benchmark Generator Contract 02 — Universal Obstruction of $\mathcal{P}_0$, Shadow Manifold $\mathcal{Z}_4^{\mathrm{BM4}}$, Quadratic Floor & Dynamic Atlas Gauge 1-Cocycle (`02_generator_universal_obstruction_and_shadow_gauge.md`)

**Contract Governance Version**: Pillar II — Cycle V (Iteration 60, Arbitrages `A1`, `A2`, `A3`)  
**Generator Script**: [`scripts/campaigns/run_generator_02_universal_obstruction_and_shadow_gauge.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scripts/campaigns/run_generator_02_universal_obstruction_and_shadow_gauge.py)  
**Plotter Script**: [`scripts/campaigns/plot_generator_02_universal_obstruction_and_shadow_gauge.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scripts/campaigns/plot_generator_02_universal_obstruction_and_shadow_gauge.py)  
**Certified JSON Output**: [`results/campaigns/generator_02_universal_obstruction_and_shadow_gauge.json`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/results/campaigns/generator_02_universal_obstruction_and_shadow_gauge.json) (and synchronized [`data/pbm4_vs_bm4_stress_test.json`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/data/pbm4_vs_bm4_stress_test.json))

---

## Section 1 — Mathematical Objective, Target Theorems & Formal Lean 4 Anchors

1. **Generator Identifier & Title**:
   - `Generator 02 — Universal Obstruction of P_0, Shadow Manifold Z_4^{BM4}, Quadratic Floor & Dynamic Atlas Gauge 1-Cocycle`
2. **Target Manuscript Sections, Propositions & Figures**:
   - **Subsection 7.1.2 (`Experimental Proof of the Universal Obstruction, Cotangent-Lifted Processing, and the Quadratic Floor`)**: [`chapters/ch07_benchmarks_and_lean4/sec01_unsoftened_benchmarks/content.tex`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch07_benchmarks_and_lean4/sec01_unsoftened_benchmarks/content.tex#L69-L109), **`Figure 3`** (`\label{fig:pbm4_quadratic_floor}`), **`Figure 4`** (`\label{fig:pbm4_stress_test}`, Stress-Tests **ST-1** to **ST-4**).
   - **Substantive Audit Findings Remediated (`Report 03`)**:
     - `[CRIT-7.1.2-01]`: Replacement of the static 0-step tautology $\mathcal{C}_B \circ (\mathcal{C}_B^{-1} \circ \mathcal{C}_A) \circ \mathcal{C}_A^{-1}(\mathbf{z}_0) \equiv \mathbf{z}_0$ in Protocol `ST-4` by a **genuine dynamical multi-step integration across a Jacobi chart switch** $\Phi_{h_B}^{(B), N_B} \circ T_{(A,h_A)\to(B,h_B)}^{(\mathrm{PBM4})} \circ \Phi_{h_A}^{(A), N_A}$, and replacement of `c4_global` (`mode="global"` $\in \mathcal{P}_0$) clamped at `1e-35` by the true non-commutative 3-chart Čech 1-cocycle $T_{C\to A}^{(\mathrm{PBM4})} \circ T_{B\to C}^{(\mathrm{PBM4})} \circ T_{A\to B}^{(\mathrm{PBM4})} = \mathrm{Id}$ on `mode="homological"` ($\mathcal{Z}_4^{\mathrm{BM4}}$) without any artificial clamping.
     - `[CRIT-7.1.2-02]`: Execution of the `Float128` order $6 \to 4$ crossover sweep on the **authentic 3D eccentric hierarchical triple ($e_{\mathrm{in}}=0.30, i=28^\circ$)** over an integer number of binary periods (eliminating the 2D circular $e_{\mathrm{in}}=0, i=0^\circ$ substitution and phase-oscillation slope overshoot, with **zero `np.clip`** on empirical slopes).
     - `[CRIT-7.1.2-03]`: Verification that canonical factories (`create_pbm4_integrator`, `CorrectorC4`) and paper showcase interfaces (`02_*.py`, `05_*.py`) default to `mode="homological"` ($\mathcal{Z}_4^{\mathrm{BM4}} \in \mathcal{P}_1$) while keeping `mode="global"` ($\in \mathcal{P}_0$) and `mode="cotangent_lift"` ($\in \mathcal{P}_1$) as explicit ablation baselines.
     - `[CRIT-7.1.2-04]` & `[MAJ-7.1.2-03]`: Transparent accounting of the interior FSAL stage cost and boundary corrector evaluations ($2M \to 2$ telescoping amortization), and exact alignment of the reported `ST-2` efficiency gain.
     - `[MAJ-7.1.2-01]` & `[MAJ-7.1.2-02]`: Explicit separation in `ST-1` between the pre-asymptotic order-6 regime $\mathcal{O}(\varepsilon h^6) \propto \varepsilon^1$ ($\varepsilon \ll h^2$, constant gain $\propto h^{-2}$) and the quadratic floor regime $\mathcal{O}(\varepsilon^2 h^4) \propto \varepsilon^2$ ($\varepsilon \gtrsim h^2$), plus explanation of the static $l=2$ quadrupole scaling ($\propto \varepsilon^{4/3}$) of `CotangentLiftCorrectorC4` vs full non-perturbative `HomologicalLieCorrectorC4` ($\mathcal{Z}_4^{\mathrm{BM4}}$).
     - `[MOD-7.1.2-01]` & `[MOD-7.1.2-02]`: Upgrade of `ST-3` reference trajectory to $h_{\mathrm{ref}} = h/4$ with `sample_stride = 1` ($40$ points/orbit, $4{,}000$ samples over 100 orbits) and continuous Hann-windowed Laskar NAFF frequency extraction, plus update of `scripts/audit_pbm4_divergence_check.py`.
3. **Mathematical & Theoretical Propositions Tested**:
   - **Theorem 5.14 (`Universal Obstruction in P_0`)**: Any pure-configuration generator $K_4(\mathbf{r}) \in \mathcal{P}_0$ (`mode="global"`) cannot cancel the momentum-dependent terms $(\mathbf{r}\cdot\mathbf{w})^2$ of $D_A^3 B$ on eccentric/3D orbits ($w_{\mathrm{rad}} \ne 0$), yielding $\mathcal{O}(1)$ gain ($\approx 1.0\times\text{–}1.1\times$) over raw $\mathrm{BM}_4^{(7)}$.
   - **Theorems 5.15–5.18 (`Cotangent-Lifted Corrector P_1 & Palindromic Lie Conjugator Z_4^{BM4}`)**: Lifting $K_4(\mathbf{r},\mathbf{w}) \in \mathcal{P}_1$ or applying the palindromic Lie conjugator $\mathcal{Z}_4^{\mathrm{BM4}}(h)$ eliminates the $\mathcal{O}(\varepsilon h^4)$ error channel, leaving $\Delta H_{\mathrm{PBM4}} = \mathcal{O}(\varepsilon h^6 + \varepsilon^2 h^4)$.
   - **Theorem 5.20 & Theorem 6.16 (`Covariant Symplectic Gauge Connection & Čech 1-Cocycle`)**: Across a Jacobi chart switch $\mathcal{T}_A \to \mathcal{T}_B$ and/or dyadic step change $h_A \to h_B$, the non-commutative sandwich $T_{(A,h_A)\to(B,h_B)}^{(\mathrm{PBM4})} = (\mathcal{Z}_4^{\mathrm{BM4},(B)}(h_B))^{-1} \circ \mathbf{M}_{A\to B} \circ \mathcal{Z}_4^{\mathrm{BM4},(A)}(h_A)$ suppresses the $\mathcal{O}(\varepsilon h^4)$ interfacial energy jump down to $\mathcal{O}(\varepsilon h^6)$ during multi-step dynamical evolution and satisfies the exact 3-chart Čech 1-cocycle identity $T_{C\to A} \circ T_{B\to C} \circ T_{A\to B} = \mathrm{Id}$.
4. **Formal Lean 4 Proof Modules Linked (`Appendix A`)**:
   - `03_ConormalCorrectorPBM4_Order6.lean` (`pbm4_universal_obstruction_p0`, `pbm4_effective_order_six_cancellation`, `cotangent_k4_p1_time_reversal_odd`)
   - `06_PalindromicCorrectorZ4_Order6.lean` (`z4_bm4_symmetric_conjugation_weights`, `msym_quadrupole_symmetric`)
   - `07_MultiChartGaugeTransition_Order6.lean` (`pbm4_gauge_operator_conjugation_cocycle`, `pbm4_gauge_transition_order6`)

---

## Section 2 — Paper Showcase Interfaces (`01..06_*.py`) & Dual-Precision Library Bindings

1. **Mandatory Paper Showcase Imports (`Float64`)**:
   - `01_PhaseStateAndSystem_Interface.py`: `ExtendedPhaseState`, `NBodySystem`
   - `02_TreeBuilderAndCorrector_Interface.py`: `build_optimal_tree`, `PalindromicCorrectorZ4`, `CotangentLiftCorrectorP1`, `GlobalCorrectorC4`, `CorrectorC4`, `ProcessedSymplecticIntegrator`, `create_pbm4_integrator`
   - `03_AdaptiveDyadicScheduler_Interface.py`: `DyadicScheduler`
   - `04_SymplecticSubcycle_Interface.py`: `BlanesMoan4Integrator`
   - `05_AtlasTransitions_Interface.py`: `PBM4GaugeTransitionOperator`, `apply_gauge_connection`, `compute_lie_bargmann_invariants`
2. **Quadruple-Precision (`Float128`, 35 Digits `mpmath`) Mirror Bindings**:
   - `normal_cone_nbody_float128.core.state.ExtendedPhaseState`, `core.jacobi.build_tree_from_spec`, `integrators.blanes_moan_4.BlanesMoan4Integrator`, `integrators.correctors.CorrectorC4`, `integrators.processed.apply_gauge_connection`
   - Mandatory runtime assertion on every `Float128` state:
     ```python
     assert state_f128.is_f128 is True
     assert isinstance(state_f128.q[0, 0], mpmath.mpf)
     ```
3. **Explicit Configuration Flags**:
   - `in_shadow_mode = True`
   - `corrector.mode in {"homological", "cotangent_lift", "global"}` (with `"homological"` verified as the default of `CorrectorC4` and `create_pbm4_integrator`)

---

## Section 3 — Initial Conditions, Complete Physical Horizons & Numerical Protocol (`Zero Pre-Pericenter Cutoff`)

1. **Part I (`Volet A` Panels a–b & `ST-2`) — `Float128` Order $6 \to 4$ Crossover & Quadratic Floor on the 3D Eccentric Hierarchical Triple (`[CRIT-7.1.2-02]` Resolved)**:
   - **System**: 3D eccentric inclined hierarchical triple ($m_0 = m_1 = 1.0, m_2 = 0.5, a_{\mathrm{in}} = 1.0, e_{\mathrm{in}} = 0.30, i = 28^\circ, E_0 = 40^\circ$) in `Float128` (`mp.mp.dps = 35`, `is_f128 == True` asserted) across $\varepsilon = (m_2/M_{01})(a_{\mathrm{in}}/R_{\mathrm{out}})^3 \in \{10^{-2}, 10^{-3}, 10^{-4}\}$.
   - **Integration Horizon**: Exact single-chart dyadic horizon $T_{\mathrm{cross}} = 0.32$ (spanning $h \in \{0.32, 0.16, 0.08, 0.04, 0.02, 0.01, 0.005\}$ so $N_{\mathrm{steps}} = T_{\mathrm{cross}}/h \in \{1, 2, 4, 8, 16, 32, 64\}$ is an exact power of two at every step size), plus full-period verification over $T_{\mathrm{bin}} = 2\pi/\sqrt{2} \approx 4.44288$.
   - **Observables**: Relative Hamiltonian error $\Delta E_{\mathrm{rel}}(h, \varepsilon) = |H(\mathbf{z}(T)) - H_0|/|H_0|$ for raw $\mathrm{BM}_4^{(7)}$ and $\mathrm{PBM}_4(\mathcal{Z}_4^{\mathrm{BM4}})$, and **unclipped** empirical slopes $p_{\mathrm{emp}}(h) = \ln(\Delta E(h_k)/\Delta E(h_{k+1}))/\ln(h_k/h_{k+1})$ showing the smooth transition from $p \approx 6.0$ ($\varepsilon \ll h^2$) down to $p \approx 4.0$ ($h \lesssim \sqrt{\varepsilon}$, where the quadratic floor $\mathcal{O}(\varepsilon^2 h^4)$ dominates).
2. **Part II (`ST-1`, `Volet A` Panel c & `Volet B` Panel a) — Direct Experimental Proof of the Universal Obstruction in $\mathcal{P}_0$ (`Theorem 5.14`)**:
   - **System**: 3D eccentric inclined hierarchical triple ($e_{\mathrm{in}} = 0.30, i = 28^\circ$) swept across 8 tidal perturbation scales $\varepsilon \in \{3\times 10^{-2}, 10^{-2}, 3\times 10^{-3}, 10^{-3}, 3\times 10^{-4}, 10^{-4}, 3\times 10^{-5}, 10^{-5}\}$ at fixed step $h = 0.08$ over $T = 1.60$ ($20$ steps).
   - **Two Compared Schemes**:
     1. Raw $\mathrm{BM}_4^{(7)}$ (`BlanesMoan4Integrator`);
     2. Non-perturbative palindromic Lie conjugator $\mathrm{PBM}_4$ ($\mathcal{Z}_4^{\mathrm{BM4}}$, `PalindromicCorrectorZ4`, `mode="homological"`).
   - **Two Asymptotic Scaling Regimes Verified (`[MAJ-7.1.2-02]` Resolved)**:
     - For $\varepsilon \in [10^{-3}, 3\times 10^{-2}] \gtrsim h^2$: quadratic floor $\mathcal{O}(\varepsilon^2 h^4) \propto \varepsilon^2$ dominates;
     - For $\varepsilon \in [10^{-5}, 3\times 10^{-4}] \ll h^2$: pre-asymptotic order-6 channel $\mathcal{O}(\varepsilon h^6) \propto \varepsilon^1$ dominates, yielding a constant plateau gain $\Delta E_{\mathrm{BM4}}/\Delta E_{\mathrm{PBM4}} \propto h^{-2} > 200\times$.
3. **Part III (`Volet A` Panel d, `Volet B` Panel b & `Volet C`) — Shadow Telescoping Amortization & 100-Orbit Secular Endurance (`ST-3`, `[MOD-7.1.2-01]` Resolved)**:
   - **Shadow Telescoping Sweep**: Integration over $M \in \{1, 2, 4, 8, 16, 32, 64, 128, 256\}$ steps comparing un-telescoped step-by-step physical projection ($2M$ corrector calls) vs shadow manifold telescoping $\mathcal{C}_h \circ \Phi_h^M \circ \mathcal{C}_h^{-1}$ ($2$ corrector calls total), verifying bit-exact final state equivalence ($\|\mathbf{z}_{\mathrm{tele}} - \mathbf{z}_{\mathrm{step}}\|_\infty = 0.0$) and quantifying effective force-evaluation amortization $s_{\mathrm{eff}}(M) = 6 + 16/M \to 6.0$.
   - **100-Orbit Secular Endurance (`ST-3`)**: Integration of the 3D eccentric hierarchical triple ($e_{\mathrm{in}}=0.30, i=28^\circ, \varepsilon=4\times 10^{-4}$) over $100$ binary periods ($T_{\mathrm{bin}} = 2\pi/\sqrt{2}$, $40$ steps/orbit, $4{,}000$ total steps) with **fine reference trajectory at $h_{\mathrm{ref}} = h/4$** ($160$ steps/orbit) and dense sampling **`sample_stride = 1`** ($4{,}000$ samples), extracting the fundamental orbital frequency via continuous Hann-windowed Laskar NAFF (`compute_naff_fundamental_freq`).
4. **Part IV (`ST-4`, `Volet B` Panels c–d) — Dynamical Multi-Step Atlas Gauge Transition & Non-Commutative 3-Chart Čech 1-Cocycle (`[CRIT-7.1.2-01]` Resolved)**:
   - **`ST-4a` (Dynamical Multi-Step Chart & Step-Size Transition $\Phi_{h_B}^{(B), N_B} \circ T_{(A,h_A)\to(B,h_B)}^{(\mathrm{PBM4})} \circ \Phi_{h_A}^{(A), N_A}$ in `Float64` & `Float128`)**:
     - Across $h \in \{0.08, 0.04, 0.02, 0.01, 0.005\}$, the system enters the shadow manifold in chart $\mathcal{T}_A = ((0,1),2)$, **executes $N_A = 4$ real symplectic steps $\Phi_h^{(A)}$ in chart $\mathcal{T}_A$**, transitions to chart $\mathcal{T}_B$ (or changes dyadic step $h \to h/2$ and back) via:
       1. `Naive Switch` (no gauge connection, suffering an $\mathcal{O}(\varepsilon h^4)$ interfacial energy jump);
       2. `Obstructed P_0 Switch` (`mode="global"`, leaving the momentum-dependent $\mathcal{O}(\varepsilon h^4)$ jump uncancelled on 3D eccentric orbits);
       3. `Covariant PBM4 Switch` ($T_{(A,h_A)\to(B,h_B)}^{(\mathrm{PBM4})}$ with `mode="homological"` $\mathcal{Z}_4^{\mathrm{BM4}} \in \mathcal{P}_1$),
       and **executes $N_B = 4$ real symplectic steps $\Phi_h^{(B)}$ after the transition** before exiting to physical space!
   - **`ST-4b` (Non-Commutative 3-Chart Čech 1-Cocycle Closure & Lie–Bargmann Conservation on $\mathcal{Z}_4^{\mathrm{BM4}}$)**:
     - Evaluates the cyclic composition across all three 3-body Jacobi charts $\mathcal{T}_A = ((0,1),2) \to \mathcal{T}_B = (0,(1,2)) \to \mathcal{T}_C = ((0,2),1) \to \mathcal{T}_A$ using **`mode="homological"` ($\mathcal{Z}_4^{\mathrm{BM4}}$)** in both `Float64` and `Float128` **without any clamping to `1e-35`**, recording the true floating-point cocycle defect $\|(T_{C\to A} \circ T_{B\to C} \circ T_{A\to B})(\tilde{\mathbf{z}}) - \tilde{\mathbf{z}}\|_\infty$ and total angular momentum conservation $\|\Delta\mathbf{L}\|_2$.

---

## Section 4 — Exact Measurement Protocol & Prohibited Proxies (`Zero Synthetic Proxy`)

1. **Exact Observable Definitions**:
   - Relative energy error: $\Delta E_{\mathrm{rel}} = |H(\mathbf{q}_{\mathrm{phys}}, \mathbf{p}_{\mathrm{phys}}) - H_0| / |H_0|$.
   - Unclipped local empirical slope: $p_{\mathrm{emp}}(h_k) = \ln(\Delta E_{\mathrm{rel}}(h_k)/\Delta E_{\mathrm{rel}}(h_{k+1})) / \ln(h_k/h_{k+1})$.
   - Continuous Laskar NAFF frequency shift: $\delta\nu/\nu_0 = |\nu_{\mathrm{scheme}} - \nu_{\mathrm{ref}(h/4)}| / |\nu_{\mathrm{ref}(h/4)}|$ computed over $4{,}000$ dense samples (`sample_stride = 1`).
2. **Explicit Anti-Proxy Declarations**:
   - Zero `np.clip` on empirical slopes $p_{\mathrm{emp}}(h)$;
   - Zero 0-step static tautology $\mathcal{C}_B \circ \mathcal{C}_B^{-1} = \mathrm{Id}$ in `ST-4a` ($N_A \ge 4$ and $N_B \ge 4$ real integration steps required before and after transition);
   - Zero substitution of `mode="global"` ($\mathcal{P}_0$) or `max(defect, 1e-35)` clamping in `ST-4b`.

---

## Section 5 — Automated Quantitative Assertion Gate (`PASS/FAIL` Contract Verification)

`run_generator_02_universal_obstruction_and_shadow_gauge.py` executes `verify_contract_assertions(results)` prior to persisting its JSON output and enforces the following mandatory `PASS/FAIL` thresholds:

| Assertion ID | Metric / JSON Key | Mathematical Meaning | Mandatory Threshold (`PASS` Condition) |
| :---: | :--- | :--- | :--- |
| `ASSERT-G02-01` | `crossover_3d_f128.float128_verified` & `e_in_3d` | True 35-digit `Float128` crossover sweep on 3D eccentric triple ($e_{\mathrm{in}}=0.30, i=28^\circ$) | `True`, `e_in == 0.30`, `inc_deg == 28.0`, and `max_gain_eps_1e4 >= 500.0` |
| `ASSERT-G02-02` | `crossover_3d_f128.eps_1e-04.local_slopes_pbm4` | Pre-asymptotic order $p \approx 6$ crossing over toward $p \approx 4$ without clipping | `slope_coarse >= 5.5` and `slope_fine <= 4.6` (with `unclipped == True`) |
| `ASSERT-G02-03` | `st1_obstruction_sweep` | Two-regime tidal scaling of $\mathrm{PBM}_4$ ($\mathcal{Z}_4^{\mathrm{BM4}}$) over $\mathrm{BM}_4$ | `max(gain_homological_z4) >= 200.0` |
| `ASSERT-G02-04` | `st3_secular_endurance` | 100-orbit 3D eccentric endurance ($h_{\mathrm{ref}}=h/4$, `sample_stride=1`) | `energy_suppression_factor >= 50.0` and `naff_gain_factor >= 75.0` |
| `ASSERT-G02-05` | `st4_dynamic_gauge_and_cocycle.dynamic_transition` | Multi-step dynamical transition ($N_A \ge 4, N_B \ge 4$) suppressing interfacial jump | `n_steps_before >= 4`, `n_steps_after >= 4`, and `dynamic_jump_suppression_f128 >= 100.0` |
| `ASSERT-G02-06` | `st4_dynamic_gauge_and_cocycle.cech_cocycle_homological` | 3-chart Čech 1-cocycle closure on `mode="homological"` ($\mathcal{Z}_4^{\mathrm{BM4}}$) without clamping | `max_cocycle_defect_f64 <= 1e-14` and `max_cocycle_defect_f128 <= 1e-32` (`clamped == False`) |

---

## Section 6 — Panel-by-Panel Specification of the Segmented Plates & Manuscript Publication Plate (`fig_ch07_pbm4_quadratic_floor_crossover.pdf/.png`)

### 6.0. Canonical $2\times 2$ Manuscript Publication Plate for `§7.1.2` (`fig_ch07_pbm4_quadratic_floor_crossover.pdf/.png`, Iteration 66)
- **Subplot Grid & Dimensions**: `$2\times 2$`, `figsize=(14.0, 10.2)`, `dpi=300`, comparing exclusively the single canonical processed integrator $\mathrm{PBM}_4$ ($\mathcal{Z}_4^{\mathrm{BM4}}$) against unprocessed $\mathrm{BM}_4^{(7)}$ ($s_{\mathrm{eval}} = 6$), with zero obsolete `P0`/`P1` curves and zero internal jargon prefixes.
- **Panel (a)**: *35-Digit `Float128` Energy Error Crossover $|\Delta H(T)|/|H_0|$ vs $h_{\mathrm{eff}} = h/s_{\mathrm{eval}}$ ($e_{\mathrm{in}}=0.30, i=28^\circ$)* across $\varepsilon \in \{10^{-2}, 10^{-3}, 10^{-4}\}$ ($R_{\mathrm{out}} \in \{2.92, 6.30, 13.57\}$), showing up to **$2{,}537\times$** energy gain at $\varepsilon=10^{-4}$ above $h_{\mathrm{cross}}(\varepsilon) \propto \sqrt{\varepsilon}$ and the exact $100\times$/decade vertical separation on the $\mathcal{O}(\varepsilon^2 h^4)$ floor ($1.49\times 10^{-13} \to 2.30\times 10^{-15} \to 2.52\times 10^{-17}$ at $h=0.01$).
- **Panel (b)**: *Unclipped Local Empirical Order $p_{\mathrm{emp}}(h) : 6 \to 4$ vs Geometric Midpoint Macro-Step $\sqrt{h_k h_{k+1}}$* overlaid with the universal analytical two-channel model $p_{\mathrm{th}}(h;\varepsilon) = 4 + 2/[1 + (h_{\mathrm{cross}}(\varepsilon)/h)^2]$ and the flat $\mathrm{BM}_4^{(7)}$ baseline ($p \approx 4.00$).
- **Panel (c)**: *Two-Regime Tidal Perturbation Scaling $\max_t |\Delta H(t)|/|H_0|$ vs $\varepsilon = \eta_S^3 \in [10^{-5}, 3\times 10^{-2}]$ at Fixed $h = 0.08$ (`Float64`, `use_optimized_kernels=True`)* contrasting $\mathrm{BM}_4^{(7)}$ ($\mathcal{O}(\varepsilon h^4) \propto \varepsilon^1$) against $\mathrm{PBM}_4$ transitioning across $\varepsilon_{\mathrm{cross}} \approx h^2 = 6.4\times 10^{-3}$ from the quadratic floor $\mathcal{O}(\varepsilon^2 h^4) \propto \varepsilon^2$ (slope $1.84$) to the pre-asymptotic order-6 channel $\mathcal{O}(\varepsilon h^6) \propto \varepsilon^1$ (slope $0.93$, up to **$236\times$** suppression).
- **Panel (d)**: *35-Digit `Float128` Phase-Space Trajectory Position Error Crossover $\|\mathbf{q}(T) - \mathbf{q}_{\mathrm{ref}}^{\mathrm{f128}}(T)\|_2$ vs $h_{\mathrm{eff}} = h/s_{\mathrm{eval}}$* evaluated against an unoptimized 35-digit `Float128` $\mathrm{BM}_6^{(11)}$ reference trajectory (6th-order Richardson-extrapolated at $h_{\mathrm{ref}} \in \{0.01, 0.005\}$), confirming the exact same $6 \to 4$ crossover (up to **$1{,}335\times$** position gain at $\varepsilon=10^{-4}$) and $100\times$/decade floor spacing ($2.49\times 10^{-13} \to 2.30\times 10^{-15} \to 2.49\times 10^{-17}$ at $h=0.01$).

### 6.1. Volet A — Computational Scaling, Order $6 \to 4$ Crossover & Position Crossover Plate (`gen02_plate_A_computational.pdf/.png`)
- **Subplot Grid & Dimensions**: `$2\times 2$`, `figsize=(14.0, 10.2)`, `dpi=300` (identical to the canonical publication plate `fig_ch07_pbm4_quadratic_floor_crossover.pdf/.png`).

### 6.2. Volet B — Algebraic Obstruction (`ST-1`), Secular NAFF (`ST-3`) & Dynamic Gauge Cocycle (`ST-4`) Plate (`gen02_plate_B_algebraic_invariants.pdf/.png`)
- **Subplot Grid & Dimensions**: `$2\times 2$`, `figsize=(14.0, 10.2)`, `dpi=300`, zero overlapping insets.
- **Panel (a)**: *Protocol `ST-1` — Two-Regime Tidal Suppression Gain of $\mathrm{PBM}_4$ ($\mathcal{Z}_4^{\mathrm{BM4}}$) over $\mathrm{BM}_4^{(7)}$*: suppression gain $\Delta E_{\mathrm{BM4}}/\Delta E_{\mathrm{PBM4}}$ across $\varepsilon \in [10^{-5}, 3\times 10^{-2}]$.
- **Panel (b)**: *Protocol `ST-3` — 100-Orbit 3D Eccentric Secular Endurance & Continuous Laskar NAFF Fidelity*: per-orbit maximum energy error and angular momentum conservation over 100 binary periods ($h_{\mathrm{ref}}=h/4$, `sample_stride=1`).
- **Panel (c)**: *Shadow Manifold Telescoping Amortization*: effective apex evaluations per step $s_{\mathrm{eff}}(M) = 6 + 16/M \to 6.0$ and wall-clock speedup.
- **Panel (d)**: *Protocol `ST-4b` — Non-Commutative 3-Chart Čech 1-Cocycle Closure ($T_{C\to A}\circ T_{B\to C}\circ T_{A\to B} = \mathrm{Id}$) & Lie–Bargmann Conservation on $\mathcal{Z}_4^{\mathrm{BM4}}$* in `Float64` ($\sim 10^{-16}$) and `Float128` ($\sim 10^{-34}$, un-clamped).

### 6.3. Volet C — Physical 3D Orbital Geometry & Shadow-to-Physical Manifold Offset Plate (`gen02_plate_C_physical_dynamics.pdf/.png`)
- **Subplot Grid & Dimensions**: `$1\times 2$`, `figsize=(13.6, 5.6)`, `dpi=300`, zero overlapping insets.
- **Panel (a)**: *3D Inclined Eccentric Hierarchical Triple Geometry ($e_{\mathrm{in}}=0.30, i=28^\circ$) and Non-Zero Radial Momentum $w_{\mathrm{rad}}(t) = (\mathbf{r}\cdot\mathbf{w})/r$* driving the $\mathcal{P}_0$ momentum-parity obstruction at pericenter passages.
- **Panel (b)**: *Canonical Offset $\|\tilde{\mathbf{q}}(t) - \mathbf{q}(t)\|_2 = \mathcal{O}(h^4)$ Between Shadow Manifold and Physical Manifold & Exact Physical Energy Recovery Under `observe_physical_state()`*.

---

## Section 7 — Bijective Traceability Table (`JSON Key <-> Plot Panel <-> Manuscript Table / Text`)

| JSON Key Path (`generator_02_universal_obstruction_and_shadow_gauge.json`) | Volet Plate & Panel ID | Target LaTeX Section / Figure | Certified Numerical Value (`Float64` / `Float128`) |
| :--- | :--- | :--- | :--- |
| `crossover_3d_f128.sweeps` (energy & position) | `Figure 3 — Panels (a, b, d)` | `§7.1.2, Figure 3 (a, b, d)` | `Float128` 3D eccentric ($e_{\mathrm{in}}=0.30, i=28^\circ$): $\varepsilon=10^{-2}$ (`19.1x` energy, `10.1x` pos, slope $5.17\to 4.00$), $\varepsilon=10^{-3}$ (`262.5x` energy, `134.1x` pos, slope $6.12\to 4.05$), $\varepsilon=10^{-4}$ (**`2537.4x` energy**, **`1335.4x` pos**, $\min\Delta E_{\mathrm{rel}} = 2.52\times 10^{-17}$, $\min\|\Delta\mathbf{q}\|_2 = 2.49\times 10^{-17}$, slope $6.23\to 4.53$, unclipped) |
| `st1_obstruction_sweep` | `Figure 3 — Panel (c)` & `Volet B (a)` | `§7.1.2, Figure 3 (c)` | `Float64` (`use_optimized_kernels=True`) at $h=0.08$: $\mathrm{BM}_4^{(7)} \propto \varepsilon^1$ vs $\mathrm{PBM}_4$ up to **`236.1x`** ($\mathcal{O}(\varepsilon^2 h^4) \propto \varepsilon^2$, slope $1.84$ for $\varepsilon\gtrsim h^2$; $\mathcal{O}(\varepsilon h^6) \propto \varepsilon^1$, slope $0.93$ for $\varepsilon\ll h^2$) |
| `shadow_telescoping` | `Volet A (archive)` | `§5.5 / §6.2 (Telescoping Cost)` | $M=256$: $s_{\mathrm{eff}} = 6.06$ (vs $14.03$ un-telescoped), wall-clock speedup `2.53x`, $\|\Delta\mathbf{q}_{\mathrm{tele}}\|_\infty = 0.0$ |
| `st3_secular_endurance` | `Volet B — Panel (b)` | `§7.2.2 (Secular NAFF)` | 100 Orbits ($e_{\mathrm{in}}=0.30, i=28^\circ, h_{\mathrm{ref}}=h/4, N=6000$ samples): $\max|\Delta H|/|H_0|$ reduced `146.7x` (`1.19e-9` $\to$ `8.12e-12`); Laskar NAFF drift reduced `89.0x` (`9.21e-10` $\to$ `1.03e-11`) |
| `st4_dynamic_gauge_and_cocycle.dynamic_transition` | `Volet B — Panel (c)` | `§6.5 / §7.2.1 (Atlas Gauge)` | Multi-step transition ($N_A=4 \to N_B=4$): Naive=`1.85e-14`, Covariant $\mathcal{Z}_4^{\mathrm{BM4}}$=`1.12e-18` (**`16,532x` jump suppression** at $h=0.005$) |
| `st4_dynamic_gauge_and_cocycle.cech_cocycle_homological` | `Volet B — Panel (d)` | `§5.5.4 / §6.5 (Čech 1-Cocycle)` | Unclamped 3-chart Čech 1-cocycle on $\mathcal{Z}_4^{\mathrm{BM4}}$: `Float64` $\le 6.66\times 10^{-16}$ (`1.11e-16` at $h=0.005$), `Float128` $\le 7.98\times 10^{-35}$ ($\|\Delta\mathbf{L}\|_2 \le 2.58\times 10^{-35}$) |
| `physical_3d_and_shadow_offset` | `Volet C — Panels (a, b)` | `§5.5 / §6.2` | Non-zero radial momentum $|w_{\mathrm{rad}}(t)| > 0$ at pericenter; canonical shadow offset $\|\tilde{\mathbf{q}}(t)-\mathbf{q}(t)\|_2 = \mathcal{O}(\varepsilon h^4)$ recovered to $\mathcal{O}(\varepsilon h^6)$ |

