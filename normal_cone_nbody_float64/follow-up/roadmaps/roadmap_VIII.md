# Strategic Roadmap : Remédiation Exhaustive et Rigoureuse des Discrépancies de l'Audit (Roadmap VIII)

## Strategic Overview:
- **Objectif Stratégique Général** : Éradiquer l'intégralité des 14 discrépances (éditoriales, terminologiques, et d'ingénierie numérique IEEE 754) identifiées lors de l'audit contradictoire exhaustif (`AUDIT_EXHAUSTIF_DISCREPANCES_CODE_VS_MONOGRAPHIE.md`), afin d'établir un isomorphisme mathématique, formel et opérationnel absolu entre la monographie de recherche LaTeX, le noyau de théorèmes certifiés en Lean 4, et le moteur de calcul Python `normal_cone_nbody_float64`.
- **Protocole Impératif de Remédiation Pré-Itérative** : Conformément aux directives de l'utilisateur, **avant chaque établissement effectif de plan d'itération (`implementation_plan.md`)**, une ou plusieurs solutions de remédiation technique et formelle (Options A, B, C...) devront être formellement présentées à l'utilisateur pour arbitrage et validation de trajectoire.
- **Règles d'Or & Contraintes Inviolables** :
  1. *Zéro-Adoucissement Strict* : Conservation absolue de la physique non adoucie ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$). Toute singularité est résolue par difféomorphisme analytique de Sundman ou éclatement conormal de Fulton-MacPherson.
  2. *Invariance Symplectique et de Lie-Bargmann* : Préservation intégrale de la 2-forme de Darboux et des 10 invariants de Noether à la précision machine float64 IEEE 754 ($< 10^{-14}$).
  3. *Condition de Frontière de Skeel-Gear* : Maintien strict de la clôture des macro-pas dyadiques $\sum 2^{-k} = 1$ et confinement des transitions d'atlas exclusivement aux frontières de macro-pas.
  4. *Validation Continue et Non-Régression* : 100 % des tests unitaires existants et nouveaux doivent passer avec succès à chaque itération.
- **Documents Maîtres de Référence** :
  - [`AUDIT_EXHAUSTIF_DISCREPANCES_CODE_VS_MONOGRAPHIE.md`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/AUDIT_EXHAUSTIF_DISCREPANCES_CODE_VS_MONOGRAPHIE.md)
  - [`RECETTE_ALGORITHME_INTEGRATION_CONORMALE.md`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/RECETTE_ALGORITHME_INTEGRATION_CONORMALE.md)

---

## Iteration 160 Strategy:
### Iteration 160: Phase 1 — Remédiation et Optimisation du Moteur de Calcul Python (Codebase & Algorithmes)
- **Scope & Objectives**:
  1. **[DISC-NUM-07] Fusion FSAL dans `Yoshida4Integrator` (`integrators/yoshida.py`)** : Remplacer l'appel séquentiel non fusionné à `_strang.step()` par la composition explicite à 4 évaluations effectives de force via les poids d'impulsion FSAL $(0.5 w_1, 0.5(w_1+w_0), 0.5(w_0+w_1), 0.5 w_1)$, réalignant le coût réel sur `evaluations_per_step = 4`.
  2. **[DISC-TERM-01] Formalisation de l'Opérateur de Super-Nœud dans `integrators/common.py`** : Clarifier et exposer l'interface de l'opérateur de super-nœud synchrone pour lever toute ambiguïté de nomenclature entre flot hamiltonien étendu et kick d'impulsion.
  3. **[DISC-TERM-02] Clarification de l'Opérateur de Kick d'Apex dans `potential/screened_apex.py` & `common.py`** : Documenter et expliciter formellement dans le code l'équivalence analytique $V_{\mathrm{kick}} = V_{\mathrm{apex}}^{\mathrm{hier}} + V_{\mathrm{kepler}}^{\mathrm{screened}} = V_{\mathrm{inter}}$ due à l'annulation exacte des monopôles écrantés $\pm \Theta(u_v) \frac{GM_L M_R}{r_v}$.
  4. **[DISC-TERM-03] Gradient dans `compute_screened_composite_potential` (`potential/screened_apex.py`)** : Sécuriser la routine utilitaire locale et annoter explicitement la portée de son gradient partiel par rapport au gradient complet de `compute_apex_forces`.
  5. **[DISC-NUM-01] Rectification des Commentaires de Repli Monotone dans `solvers/kepler_sundman.py`** : Remplacer la mention inexacte `(Brent/bisection)` par la description exacte de la bissection déterministe pure à 60 étapes garantissant $10^{-18}$ de résolution.
  6. **[DISC-NUM-06] Harmonisation du Signe de la Matrice Unimodulaire dans `core/transitions.py`** : Aligner la convention de signe de $\mathbf{M}_{\mathrm{coord}} \in \mathrm{SL}(2, \mathbb{R})$ sur celle du manuscrit tout en garantissant $\det(\mathbf{M}) = +1$.
- **Solutions de Remédiation à Soumettre Préalablement à l'Utilisateur** :
  - *Pour Yoshida 4 [DISC-NUM-07] :*
    - **Option 1A (Recommandée)** : Développer la boucle de fusion FSAL analytique des demi-kicks dans `yoshida.py` (gain immédiat de 33% de performance, conformité exacte avec la Table 2).
    - **Option 1B** : Conserver la composition séquentielle de Strang et corriger la propriété `evaluations_per_step` à 6.
  - *Pour `step_supernode_kick` [DISC-TERM-01] :*
    - **Option 2A (Recommandée)** : Ajouter une fonction wrapper explicite `step_supernode_kick` (ou `step_supernode_flow`) dans `common.py` qui documente et encapsule formellement le flot hamiltonien couplé sous temps de Sundman.
    - **Option 2B** : Conserver l'intégration directe dans `step_drift_and_kepler` et harmoniser les diagrammes d'architecture du papier.
  - *Pour le Signe de $\mathbf{M}_{\mathrm{coord}}$ [DISC-NUM-06] :*
    - **Option 3A (Recommandée)** : Aligner `transitions.py` sur la convention du papier en inversant les signes des coefficients pour que $\mathbf{r}_1' = \mathbf{q}_3 - \mathbf{q}_2$, et vérifier la non-régression de l'atlas cotangent $\operatorname{Sp}(12, \mathbb{R})$.
    - **Option 3B** : Conserver le signe actuel en Python et modifier le texte LaTeX pour refléter la convention $\mathbf{r}_1' = \mathbf{q}_2 - \mathbf{q}_3$.
- **Dependencies & Continuity**: Sécurise l'ensemble du socle Python avant toute modification documentaire ou formelle.
- **Gating Criteria**:
  - `Yoshida4Integrator` validé avec exactement 4 évaluations de force par pas.
  - 100 % de succès sur l'ensemble de la suite de tests unitaires (150/150 tests passants).
  - Aucune dérive de symplecticité ni de conservation d'énergie sur les benchmarks.

---

## Iteration 161 Strategy:
### Iteration 161: Phase 2 — Mise à Niveau du Noyau Formel Lean 4 (Règle de Leibniz Complète pour $\mathcal{C}_6$)
- **Scope & Objectives**:
  1. **[DISC-NUM-08] Règle de Leibniz Complète dans le Certificat Formel Lean 4** :  
     Mettre à jour le fichier `redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/sec05_symplectic_correctors/code/02_CorrectorBM6_Theorems.lean` :
     - Enrichir la définition formelle de `grad_pi` pour inclure le terme variationnel quadratique de Leibniz :
       $$\nabla_{\mathbf{p}_j} K_6 = \frac{\kappa_6}{m_j} \left[ (1 + \|\mathbf{v}_j\|^2) \mathbf{F}_j + 2 (\mathbf{F}_j \cdot \mathbf{v}_j) \mathbf{v}_j \right].$$
     - Démontrer formellement le lemme d'isomorphisme entre le gradient variationnel du champ hamiltonien de Lie et l'implémentation opérationnelle.
- **Solutions de Remédiation à Soumettre Préalablement à l'Utilisateur** :
  - **Option 1 (Recommandée)** : Reformuler directement la structure `CorrectorBM6System` en Lean 4 avec le champ complet à deux termes et ajuster la preuve de non-divergence conormale.
  - **Option 2** : Décomposer le certificat formel en deux théorèmes étanches : (i) la composante de cisaillement linéaire $(1+\|\mathbf{v}\|^2)\mathbf{F}$, et (ii) la composante variationnelle de projection d'impulsion $2(\mathbf{F}\cdot\mathbf{v})\mathbf{v}$.
- **Dependencies & Continuity**: S'appuie sur la validation du gradient exact vérifié à $10^{-11}$ en Python lors de l'Itération 160.
- **Gating Criteria**:
  - Fichier `02_CorrectorBM6_Theorems.lean` compilé avec succès sans `sorry` ni axiome ad-hoc.
  - Isomorphisme parfait certifié entre la définition Lean 4 et les lignes 256-257 de `correctors.py`.

---

## Iteration 162 Strategy:
### Iteration 162: Phase 3 — Harmonisation Textuelle et Précisions Numériques des Chapitres 3, 4 et 5 (LaTeX)
- **Scope & Objectives**:
  1. **[DISC-DOC-02] Scission Éditoriale Ch. 3 / Ch. 5** : Ajouter une note d'articulation épistémique formelle dans `ch03_three_body_paradigm/sec02_classical_regularizations/content.tex` pointant vers la Section 5.2 pour la résolution conormale effective de Sundman par Stumpff-Halley-McGehee.
  2. **[DISC-NUM-02] Formalisation de la T-Invariance des Pas Rétrogrades** : Documenter dans la Proposition 5.2.2 la factorisation algorithmique par involution d'impulsion $T : (\mathbf{r}, \mathbf{w}) \mapsto (\mathbf{r}, -\mathbf{w})$ garantissant la réversibilité bit-exacte sous $\Delta t < 0$.
  3. **[DISC-NUM-03] Garde Matérielle Hyperbolique IEEE 754** : Consigner explicitement la constante de sécurité physique $Z_{\mathrm{MIN}} = -500\,000.0$ dans l'annexe numérique du livre.
  4. **[DISC-NUM-04] Formulation Rationnelle Fermée du Cauchy Mollifier** : Insérer une sous-remarque dans le Chapitre 4 (après l'Éq. 4.69) démontrant la stricte équivalence $\mathcal{C}^\infty$ du quotient d'exponentielles $e^{-1/s}/(e^{-1/s}+e^{-1/(1-s)})$ avec l'intégrale convolutive normalisée.
  5. **[DISC-NUM-05] Stabilisation par Trigger de Schmitt sur le Ratio de Hill** : Ajouter un corollaire à la Proposition 4.5 formalisant la bande d'hystérésis $[\eta_{\mathrm{enter}}, \eta_{\mathrm{exit}}] = [0.60, 0.75]$ comme condition nécessaire d'élimination du broutement topologique.
  6. **[DISC-DOC-03] Inscription Explicite des Constantes Homologiques** : Imprimer explicitement $\kappa_4 = \frac{1}{2880}$ et $\kappa_6 = \frac{1}{120960}$ dans les Équations display de la Section 5.5.
- **Solutions de Remédiation à Soumettre Préalablement à l'Utilisateur** :
  - **Option 1 (Recommandée)** : Intégrer ces précisions au fil du texte sous forme de "Remarks on Computational Implementation & IEEE 754 Stability" dans chaque section respective.
  - **Option 2** : Créer une Annexe dédiée en fin de monographie intitulée "Appendix: Numerical Engineering & Computational Constants of the Conormal Engine".
- **Dependencies & Continuity**: Harmonise la monographie de recherche avec les choix validés en Phase 1 et 2.
- **Gating Criteria**:
  - Manuscrit LaTeX compilable sans erreur (`pdflatex` / `lualatex`).
  - Éradication complète des discrépances documentaires DOC-02, DOC-03, NUM-02, NUM-03, NUM-04, NUM-05.

---

## Iteration 163 Strategy:
### Iteration 163: Phase 4 — Rédaction Intégrale de la Section 5.6 du Manuscrit LaTeX ("Dyadic Multirate Scheduling")
- **Scope & Objectives**:
  1. **[DISC-DOC-01] Rédaction Complète de `sec06_dyadic_multirate_scheduling`** :
     - Créer le fichier `redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/sec06_dyadic_multirate_scheduling/content.tex`.
     - Intégrer formellement les équations (5.106) à (5.132) citées par le code :
       - Formulation de l'échelle conormale locale $\tau_{ij} = \min(\tau_{\mathrm{kin}}, \tau_{\mathrm{dyn}}, \tau_{\mathrm{peri}})$.
       - Condition de partition dyadique de Skeel-Gear : $h_k = h_{\mathrm{macro}} / 2^k$ avec $\sum_{m=1}^{2^k} h_k = h_{\mathrm{macro}}$.
       - Théorème de Découplage Feuille/Apex : démonstration formelle de la neutralisation de $\tau_{\mathrm{peri}}$ pour les paires de feuilles résolues analytiquement par Stumpff-Sundman.
       - Théorème d'Invariance de l'Hamiltonien d'Ombre sous Gel Topologique Intra-Macrobloc.
     - Mettre à jour le fichier maître `main.tex` du Chapitre 5 pour inclure cette nouvelle section.
- **Solutions de Remédiation à Soumettre Préalablement à l'Utilisateur** :
  - **Option 1 (Recommandée)** : Rédiger la Section 5.6 comme une section majeure autonome du Chapitre 5, dotée de ses propres définitions, théorèmes et preuves formelles.
  - **Option 2** : Étendre la sous-section 5.5.4 existante (`sub:dyadic_boundary_telescoping`) pour y fondre l'ensemble des théorèmes et équations dyadiques sans créer de section 5.6 séparée (nécessiterait de mettre à jour les docstrings de `dyadic_scheduler.py`).
- **Dependencies & Continuity**: Comble la discrépance critique documentaire majeure DOC-01.
- **Gating Criteria**:
  - Chapitre 5 compilé avec succès incluant la Section 5.6 complète.
  - Isomorphisme absolu entre les équations (5.106)-(5.132) du texte et les lignes d'exécution de `dyadic_scheduler.py`.

---

## Iteration 164 Strategy:
### Iteration 164: Phase 5 — Certification Globale de Non-Régression, Benchmarks et Synchronisation Démonstrateur
- **Scope & Objectives**:
  1. **Exécution de la Suite Complète des Benchmarks Canoniques** : Figure-8, Burrau, Collision Radiale, Sitnikov, Arenstorf, Jovien, et 4-Corps Hiérarchique sur les 6 intégrateurs (Strang, Yoshida, BM4, BM6, PBM4, PBM6).
  2. **Contrôle Formel de Non-Régression** : Exécution intégrale de la suite de 150+ tests unitaires Python et vérification de la conservation des 10 invariants de Lie-Bargmann à $< 10^{-14}$.
  3. **Régénération des Planches de Publication 300 DPI** : Actualiser les 4 dashboards graphiques (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`, `dashboard_multi_integrator_invariants_300dpi.png`).
  4. **Synchronisation Démonstrateur** : Mettre à jour `showcase_data.json` et resynchroniser `interactive_conormal_showcase.html`.
- **Solutions de Remédiation à Soumettre Préalablement à l'Utilisateur** :
  - **Option 1 (Recommandée)** : Exécution automatisée en batch de tous les benchmarks avec génération d'un rapport de validation comparative finale.
  - **Option 2** : Focus approfondi sur les deux cas extrêmes les plus sensibles (Choc frontal radial $L=0$ et 4-corps hiérarchique au péricentre).
- **Dependencies & Continuity**: Clôture de l'ensemble des remédiations de la Roadmap VIII.
- **Gating Criteria**:
  - 100 % des tests unitaires validés (zéro échec, zéro régression).
  - 100 % des benchmarks validés avec scorecards [PASS].
  - Concordance parfaite et certifiée entre le code Python, le manuscrit LaTeX et les théorèmes Lean 4.

---

## Multi-Iteration Checklist:
- [x] Iteration 160: Phase 1 — Remédiation et Optimisation du Moteur de Calcul Python (Codebase & Algorithmes)
- [x] Iteration 161: Phase 2 — Mise à Niveau du Noyau Formel Lean 4 (Règle de Leibniz Complète pour $\mathcal{C}_6$)
- [x] Iteration 162: Phase 3 — Harmonisation Textuelle et Précisions Numériques des Chapitres 3, 4 et 5 (LaTeX)
- [x] Iteration 163: Phase 4 — Rédaction Intégrale de la Section 5.6 du Manuscrit LaTeX ("Dyadic Multirate Scheduling")
- [x] Iteration 164: Phase 5 — Certification Globale de Non-Régression, Benchmarks et Synchronisation Démonstrateur
