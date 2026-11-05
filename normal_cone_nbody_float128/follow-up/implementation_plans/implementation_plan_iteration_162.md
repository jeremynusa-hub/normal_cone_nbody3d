# Implementation Plan of Iteration 162 : Phase 3 — Harmonisation Textuelle et Précisions Numériques des Chapitres 3, 4 et 5 (LaTeX & Code)
Created at: 2026-10-09T08:50:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur a validé le démarrage de la Phase 3 de la Roadmap VIII avec les arbitrages stricts suivants :
   - `DISC-DOC-02` : Insérer dans le Chapitre 3 (Section 3.2) une note d'articulation conormale vers la Section 5.2 pour la résolution effective et symplectique de Sundman via Stumpff-Halley.
   - `DISC-NUM-02` : Placer le commentaire détaillé sur la $T$-invariance microscopique sous pas rétrogrades $\Delta t < 0$ dans le code Python (`blanes_moan_6.py` / `kepler_sundman.py`), et **non dans la monographie**.
   - `DISC-NUM-03` : Placer le commentaire détaillé sur la garde matérielle hyperbolique IEEE 754 $Z_{\mathrm{MIN}} = -500\,000.0$ dans le code Python (`kepler_sundman.py`), et **non dans la monographie**.
   - `DISC-NUM-04` : Insérer dans le Chapitre 4 (Section 4.2) la formulation analytique rationnelle fermée en exponentielles du Cauchy mollifier et la hiérarchie des splines de Hermite centrées C2 à C8.
   - `DISC-NUM-05` : Formaliser dans le Chapitre 4 (Section 4.3) le protocole couplé à deux étages : (1) l'anticipation prédictive amont de l'entrée dans la sphère de marée de Hill par l'ordonnanceur dyadique ($v_{\mathrm{rad}} < 0, t_{\mathrm{approach}} \le 2\Delta t$), et (2) le verrouillage par bande d'hystérésis de Schmitt $[\eta_{\mathrm{enter}}, \eta_{\mathrm{exit}}] = [0.60, 0.75]$ éliminant tout broutement topologique.
   - `DISC-DOC-03` : Inscrire explicitement en formules display dans le Chapitre 5 (Section 5.5) les constantes homologiques canoniques $\kappa_4 = \frac{1}{2880}$ et $\kappa_6 = \frac{1}{120960}$.
2. **Implicit Strategic Value**: Établir une cohérence absolue et rigoureuse entre le corpus théorique (monographie de recherche) et le moteur de calcul (Python), en respectant la frontière épistémique voulue par l'utilisateur : les concepts physiques et mathématiques généraux dans le texte, et les détails d'ingénierie numérique machine dans les commentaires du code source.
3. **Alignment Criteria**:
   - Manuscrit LaTeX exempt d'erreurs de syntaxe et respectant les environnements typographiques existants (`rembox`, `thmbox`, `propbox`, `corbox`).
   - Code Python documenté de façon exhaustive sans rupture d'implémentation ni altération des invariants.
   - 100 % de succès sur l'ensemble de la suite de tests unitaires (151/151 tests).

## Guidelines:
- **Édition de Code Source & LaTeX** : Interdiction absolue de `write_to_file` avec `overwrite: true` sur tout fichier source ou LaTeX existant. Utilisation exclusive de `replace_file_content` ciblée sur les blocs sémantiques minimaux.
- **Règles d'Or Physiques** : Préservation intégrale du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$) et de l'invariance de Lie-Bargmann.
- **Rigueur Typographique LaTeX** : Respect des balises mathématiques, citations BibTeX existantes et équilibrage strict des environnements.

## Plan:

### Step 1: Articulation Épistémique Ch. 3 -> Ch. 5 (DISC-DOC-02)
- **Context**: [`redaction/unsoftened_conormal_nbody3d/chapters/ch03_three_body_paradigm/sec02_classical_regularizations/content.tex`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch03_three_body_paradigm/sec02_classical_regularizations/content.tex).
- **Objective**: Insérer après le Théorème 3.1 une remarque d'articulation conormale reliant la transformation de Sundman $\mathrm{d}t = r\,\mathrm{d}u$ au solveur symplectique universel de Stumpff-Halley-McGehee explicité en Section 5.2.
- **Deliverables**:
  - `Ch03SundmanBridge`: Section 3.2 enrichie de la remarque d'articulation.
- **Process**: Modifier `content.tex` via `replace_file_content`.

### Step 2: Documentation Python des Invariants Rétrogrades et de la Garde Hyperbolique (DISC-NUM-02 & DISC-NUM-03)
- **Context**: [`normal_cone_nbody_float64/solvers/kepler_sundman.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/kepler_sundman.py) et [`normal_cone_nbody_float64/integrators/blanes_moan_6.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/blanes_moan_6.py).
- **Objective**:
  - Dans `blanes_moan_6.py` et `kepler_sundman.py`, insérer la documentation formelle de la $T$-invariance microscopique et de la factorisation $(\mathbf{r}, \mathbf{w}) \mapsto (\mathbf{r}, -\mathbf{w})$ sous sous-étages rétrogrades $\Delta t < 0$ ($a_i < 0, b_i < 0$).
  - Dans `kepler_sundman.py`, expliciter en détail la justification physique et matérielle IEEE 754 de la constante de sécurité $Z_{\mathrm{MIN}} = -500\,000.0$ prévenant l'overflow dans `sinh` lors de plongées hyperboliques sans altérer la physique unsoftened ($d_0 = 0$).
- **Deliverables**:
  - `PythonRetrogradeDocumentation`: Documentation des pas rétrogrades.
  - `PythonHyperbolicGuardDocumentation`: Documentation de la constante $Z_{\mathrm{MIN}}$.
- **Process**: Modifier les fichiers Python ciblés via `replace_file_content`.

### Step 3: Formulation Fermée du Cauchy Mollifier & Splines de Hermite (DISC-NUM-04)
- **Context**: [`redaction/unsoftened_conormal_nbody3d/chapters/ch04_fulton_macpherson_splitting/sec02_screened_apex_potential/content.tex`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch04_fulton_macpherson_splitting/sec02_screened_apex_potential/content.tex).
- **Objective**: Insérer dans la description de la fonction de screening $\Theta(u)$ (après l'Éq. 4.66) la formulation rationnelle fermée en quotient d'exponentielles $e^{-1/s}/(e^{-1/s}+e^{-1/(1-s)})$ et documenter la hiérarchie des splines de Hermite centrées en $\xi = s - 0.5$ (C2, C4, C6, C8) adaptées aux ordres de composition symplectique.
- **Deliverables**:
  - `Ch04CauchyMollifierClosedForm`: Formules display et remarques de régularité insérées.
- **Process**: Modifier `content.tex` via `replace_file_content`.

### Step 4: Protocole Couplé Anticipation Dyadique & Hystérésis de Hill (DISC-NUM-05)
- **Context**: [`redaction/unsoftened_conormal_nbody3d/chapters/ch04_fulton_macpherson_splitting/sec03_synchronous_supernodes/content.tex`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch04_fulton_macpherson_splitting/sec03_synchronous_supernodes/content.tex).
- **Objective**: Insérer un corollaire formel après la Proposition 4.5 détaillant le protocole couplé à deux barrières : (1) détection prédictive d'approche à $v_{\mathrm{rad}} < 0$ par l'ordonnanceur dyadique pour subdiviser le pas avant l'entrée dans la sphère d'influence, et (2) trigger de Schmitt $[\eta_{\mathrm{enter}}, \eta_{\mathrm{exit}}] = [0.60, 0.75]$ garantissant l'absence de broutement topologique.
- **Deliverables**:
  - `Ch04HillHysteresisCorollary`: Corollaire formel du couplage anticipation-hystérésis.
- **Process**: Modifier `content.tex` via `replace_file_content`.

### Step 5: Inscription Explicite des Constantes Homologiques $\kappa_4, \kappa_6$ (DISC-DOC-03)
- **Context**: [`redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/sec05_symplectic_correctors/content.tex`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/sec05_symplectic_correctors/content.tex).
- **Objective**: Insérer les valeurs numériques exactes $\kappa_4 = \frac{1}{2880}$ et $\kappa_6 = \frac{1}{120960}$ dans les formules display des équations homologiques des correcteurs symplectiques $\mathcal{C}_4$ et $\mathcal{C}_6$.
- **Deliverables**:
  - `Ch05HomologicalConstantsDisplayed`: Constantes imprimées dans les équations.
- **Process**: Modifier `content.tex` via `replace_file_content`.

### Step 6: Contrôle de Non-Régression et Commit Git
- **Context**: Code Python et documents LaTeX modifiés.
- **Objective**:
  - Exécuter la suite de 151 tests unitaires Python (`pytest tests/unit/`).
  - Vérifier l'absence d'erreurs de syntaxe LaTeX.
  - Créer le commit Git scellant la Phase 3 de la Roadmap VIII.
- **Deliverables**:
  - `PytestNonRegressionPassed`: 151/151 tests validés.
  - `GitCommitPhase3`: Commit Git enregistré.
- **Process**: Exécuter `pytest` et `git commit`.

## Acceptance Criteria:
### Step 1:
- La Section 3.2 contient la passerelle formelle vers la Section 5.2 pour le solveur de Stumpff-Halley.
### Step 2:
- `kepler_sundman.py` et `blanes_moan_6.py` contiennent les explications mathématiques de la $T$-invariance sous $\Delta t < 0$ et de la garde $Z_{\mathrm{MIN}}$.
### Step 3:
- La Section 4.2 explicite la formule fermée du quotient d'exponentielles et la hiérarchie C2-C8.
### Step 4:
- La Section 4.3 formalise le corollaire couplant l'anticipation dyadique prédictive amont et la bande d'hystérésis de Hill $[0.60, 0.75]$.
### Step 5:
- La Section 5.5 imprime en display math $\kappa_4 = \frac{1}{2880}$ et $\kappa_6 = \frac{1}{120960}$.
### Step 6:
- 151/151 tests unitaires Python validés à 100 %.
- Commit Git créé sur la branche `main`.

## Observations:
> [!IMPORTANT]
> Conformément à l'arbitrage de l'utilisateur, les détails matériels IEEE 754 (garde hyperbolique) et les symétries discrètes d'implémentation des pas rétrogrades sont consignés exclusivement dans la documentation du code Python, préservant la pureté théorique de la monographie.
