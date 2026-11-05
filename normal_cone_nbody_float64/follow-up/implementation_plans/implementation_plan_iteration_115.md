# Implementation Plan of Iteration 115 : Couche 3 (potential/) - Crible Conormal $\Theta(u_v)$, Potentiel d'Apex $\mathcal{C}^{1,1}$, Détection de Marée de Hill et Newton Direct
Created at: 2026-10-07T17:43:00+02:00

## Context:
1. **Explicit Request**: Implémenter et certifier unitairement la Couche 3 (`potential/`) de la bibliothèque `normal_cone_nbody_float64`, conformément à la feuille de route stratégique Roadmap V. Le développement doit rigoureusement respecter le protocole récurrent à 5 étapes :
   - Étape 1 : Identifier et relire les sections théoriques du papier (Sections 4.1–4.3 de [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/)).
   - Étape 2 : Relire l'ancien code associé (`models/normal_cone_integrator/cas6_nbody3d_async/potential/` et modules associés) et identifier ce qui est conforme, supprimer l'ancien `apex.py` non criblé, éliminer les mentions `cas6`, et remplacer le seuil arbitraire $\delta=0.5$ par le critère intrinsèque de marée de Hill $\eta_S \le 3^{-1/3}$.
   - Étape 3 : Établir un document de spécification claire du code à écrire avec les tests et attendus.
   - Étape 4 : Rédiger le code modulaire optimisé en fonctions pures, avec typage strict et en-tête descriptif citant les sections de la monographie et les théorèmes certifiés Lean 4.
   - Étape 5 : Rédiger et exécuter la suite de tests unitaires miroirs (`tests/unit/potential/`) avec un taux de réussite impératif de 100 %.
2. **Implicit Strategic Value**: Fournir les champs de forces régularisés $\mathcal{C}^{1,1}$ universels ($M_{\mathrm{apex}}, L_{\mathrm{apex}} < \infty$) sur la variété de Fulton--MacPherson $\mathbf{Q}[N]$, garantissant l'annulation rigoureuse des pôles fantômes lors des coïncidences de centres de masse ($u_v \le 1$), tout en préservant l'exactitude asymptotique multipolaire ($u_v \ge 2$) et l'invariance d'échelle sous homothétie.
3. **Alignment Criteria**: Gradient analytique conforme aux différences finies à $< 10^{-11}$, annulation stricte de la force conormale pour $u_v \le 1$ ($\Theta \equiv 0$), invariance d'échelle exacte du détecteur de Hill ($\eta_S(\lambda \mathbf{q}) = \eta_S(\mathbf{q})$), et zéro adoucissement ($\varepsilon = 0, d_0 = 0$).

## Guidelines:
- **Method**: Application stricte du protocole à 5 étapes.
- **Red Lines**:
  - Sanctuarisation absolue en lecture seule du répertoire [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/).
  - Zéro occurrence de `"cas6"` ou `"cas 6"`.
  - Zéro adoucissement ni seuil de distance artificielle ($\varepsilon = 0, d_0 = 0$).
  - Aucun couplage ascendant vers les couches 4, 5, 6 ou 7. Dépendances limitées exclusivement aux Couches 1 (`core/`) et 2 (`solvers/`).
- **Architecture**:
  - `normal_cone_nbody_float64/potential/screening.py` : fonction de coupure $\Theta(u_v)$ sur ratio adimensionné avec modes Cauchy $\mathcal{C}^\infty$ et Hermite $\mathcal{C}^2$ ($\sigma(1)=0, \sigma(2)=1, \sigma'(1)=\sigma'(2)=\sigma''(1)=\sigma''(2)=0$).
  - `normal_cone_nbody_float64/potential/screened_apex.py` : potentiel d'apex criblé $V_{\mathrm{apex}, v}(\mathbf{q}) = V_v(\mathbf{q}) + \Theta(u_v) \frac{G M_L M_R}{\|\mathbf{r}_v\|}$, tenseur quadrupolaire $\mathcal{Q}_v(\mathbf{s})$ et forces analytiques exactes.
  - `normal_cone_nbody_float64/potential/supernode_detector.py` : détection adimensionnée par le ratio de marée de Hill $\eta_S \le 3^{-1/3}$ et le seuil d'énergie naturel de Marchal--Sundman $\delta_{\mathrm{cluster}}^{\mathrm{nat}}$.
  - `normal_cone_nbody_float64/potential/direct.py` : potentiel et forces newtoniennes pures vectorisées $\mathcal{O}(N^2)$ servant de vérité terrain non régularisée.
- **Code Headers**: En-têtes exhaustifs citant les sections du papier de recherche et les théorèmes Lean 4.

## Plan:
### Step 1: Relecture Théorique des Sections 4.1–4.3 de la Monographie
- **Context**: Le Chapitre 4 fonde la décomposition à 5 termes du hamiltonien, la régularité Lipschitz universelle $\mathcal{C}^{1,1}$, la fonction de crible $\Theta(u_v)$ et le critère de compacité de Hill.
- **Objective**: Extraire les formulations mathématiques exactes : polynôme d'Hermite de degré 5, mollificateur de Cauchy, tenseur quadrupolaire $\mathcal{Q}_v$, formule du ratio de Hill $\eta_S$ et seuil naturel $\delta_{\mathrm{cluster}}^{\mathrm{nat}}$.
- **Deliverables**:
  - `Layer3TheorySynthesis`: Synthèse des formules exactes et des labels de théorèmes Lean 4 associés (`apex_universal_lipschitz_bound`, `collision_diagonal_non_transversality`).
- **Process**: Relecture ciblée de `chapters/ch04_fulton_macpherson_splitting/`.

### Step 2: Audit Approfondi du Code Hérité (`cas6_nbody3d_async/potential/`)
- **Context**: L'ancien code comporte `apex.py`, `screened_apex.py`, `screening.py` et `direct.py`.
- **Objective**: Auditer les sources héritées, éliminer l'ancien module `apex.py` redondant non criblé, supprimer les 42 occurrences de `cas6` dans ce répertoire, et remplacer les seuils de distance ad hoc par le ratio de Hill.
- **Deliverables**:
  - `Layer3RefactoringMatrix`: Matrice d'assainissement et de migration des potentiels.
- **Process**: Inspection systématique des sources et confrontation avec `python_codebase_mapping_and_audit_report.md`.

### Step 3: Spécification Formelle et Contrats d'Interface de la Couche 3
- **Context**: Définir les signatures d'API, types et tolérances avant le codage.
- **Objective**: Spécifier les 4 modules et leurs tests unitaires miroirs avec critères d'invariance stricts.
- **Deliverables**:
  - `Layer3ModuleSpecification`: Spécification technique des modules de `potential/`.
- **Process**: Rédaction des spécifications d'invariance et de gradient.

### Step 4: Implémentation Modulaire et Optimisée des 4 Modules de `potential/`
- **Context**: Créer le package `normal_cone_nbody_float64/potential/`.
- **Objective**: Rédiger les 4 modules Python avec en-têtes explicites, vectorisation NumPy et gradients analytiques conformes.
- **Deliverables**:
  - `normal_cone_nbody_float64/potential/__init__.py`
  - `normal_cone_nbody_float64/potential/screening.py`
  - `normal_cone_nbody_float64/potential/screened_apex.py`
  - `normal_cone_nbody_float64/potential/supernode_detector.py`
  - `normal_cone_nbody_float64/potential/direct.py`
- **Process**: Implémentation séquentielle avec typage strict PEP 484/544.

### Step 5: Implémentation et Exécution des Tests Unitaires Miroirs (`tests/unit/potential/`)
- **Context**: Chaque module possède son fichier de test dédié.
- **Objective**: Rédiger et exécuter la suite sous `pytest` avec 100 % de succès.
- **Deliverables**:
  - `tests/unit/potential/test_screening.py`
  - `tests/unit/potential/test_screened_apex.py`
  - `tests/unit/potential/test_supernode_detector.py`
  - `tests/unit/potential/test_direct.py`
- **Process**: Exécution de `pytest tests/unit/potential/` et validation des gradients à $< 10^{-11}$.

## Acceptance Criteria:
### Step 1:
- Recensement complet des références de la monographie : crible conormal (Sec. 4.2), régularité Lipschitz $\mathcal{C}^{1,1}$ (Sec. 4.2), seuil de marée de Hill et Marchal--Sundman (Sec. 4.3).
### Step 2:
- Rapport d'audit de refactorisation attestant la suppression de l'ancien `apex.py` et l'élimination des mentions `cas6`.
### Step 3:
- Spécifications formelles définissant les signatures typées et les tolérances analytiques de gradient.
### Step 4:
- 4 modules de potentiel implémentés avec en-têtes complets citant l'article et Lean 4.
- Zéro adoucissement ($\varepsilon = 0, d_0 = 0$).
### Step 5:
- 100 % de réussite sur `pytest tests/unit/potential/`.
- Annulation exacte du crible conormal : $\Theta(u_v) = 0$ pour tout $u_v \le 1.0$.
- Gradient analytique du potentiel d'apex criblé conforme aux différences finies à $< 10^{-11}$.
- Invariance d'échelle du détecteur de Hill : $\eta_S(\lambda \mathbf{q}) = \eta_S(\mathbf{q})$ à $< 10^{-15}$.
- Énergie potentielle newtonienne directe rigoureusement identique à l'évaluateur de référence.

## Observations:
> [!IMPORTANT]
> Le ratio conormal $u_v = \|\mathbf{r}_v\| / (\bar{\rho}_{L(v)} + \bar{\rho}_{R(v)})$ doit impérativement utiliser l'hyper-rayon normalisé par la masse $\bar{\rho}_C = \rho_C / \sqrt{M_C} \in [L]$ afin que $u_v$ soit rigoureusement adimensionnel et invariant par homothétie spatiale.

> [!NOTE]
> La détection des super-nœuds par le critère de marée de Hill $\eta_S \le 3^{-1/3} \approx 0.69336$ élimine tout seuil arbitraire dépendant de l'échelle spatiale du système, assurant une transférabilité universelle des configurations astronomiques (du système jovien aux amas stellaires denses).
