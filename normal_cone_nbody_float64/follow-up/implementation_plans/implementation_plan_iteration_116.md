# Implementation Plan of Iteration 116 : Couche 4 (integrators/) - Suite Symplectique Modulaire $\mathcal{S}_2, \mathcal{S}_4, \mathrm{BM}_4^{(7)}, \mathrm{BM}_6^{(11)}$, Correcteurs Canoniques $\mathcal{C}_{2p}$ et Intégrateurs Traités $\mathcal{P}\mathrm{BM}$
Created at: 2026-10-07T17:49:00+02:00

## Context:
1. **Explicit Request**: Implémenter et certifier unitairement la Couche 4 (`integrators/`) de la bibliothèque académique `normal_cone_nbody_float64`, conformément à la feuille de route stratégique Roadmap V. Le cycle d'itération doit rigoureusement appliquer le protocole en 5 étapes ordonné par l'utilisateur :
   - Étape 1 : Identifier et relire les sections théoriques du papier (Chapitre 5 de [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/) : Sections 5.1–5.5, Table 1 des coefficients Blanes--Moan, équations homologiques des correcteurs canoniques, télescopage dyadique).
   - Étape 2 : Relire l'ancien code associé (`models/normal_cone_integrator/cas6_nbody3d_async/integrators/`) en s'appuyant sur [`python_codebase_mapping_and_audit_report.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/python_codebase_mapping_and_audit_report.md) ; identifier les composants conformes, éliminer les 38 mentions `cas6`, et concevoir de novo `BM_6_11` (coefficients certifiés à 25 décimales) et `ProcessedSymplecticIntegrator`.
   - Étape 3 : Établir un document de spécification claire du code à écrire avec typage strict PEP 544 (`SymplecticIntegrator` protocol), invariants, ordres de convergence et critères d'acceptation.
   - Étape 4 : Rédiger le code modulaire optimisé en fonctions pures et classes de propagation, avec en-têtes explicites citant l'article de recherche et les théorèmes Lean 4 certifiés (`PaperUnifiedFoundations.lean`, `BM6Order6SymmetricComposition.lean`, `02_LieBargmannExactInvariance_Theorems.lean`).
   - Étape 5 : Rédiger et exécuter la suite de tests unitaires miroirs (`tests/unit/integrators/`) avec un taux de réussite impératif de 100 %.
2. **Implicit Strategic Value**: Fournir le moteur d'intégration symplectique modulaire et haute précision du problème des $N$ corps sans adoucissement, surmontant la barrière d'impossibilité de Sheng--Suzuki ($w_i < 0$) via la régularisation de Stumpff-Sundman et de McGehee des Couches 1–3, et garantissant l'invariance exacte des 10 générateurs de Lie--Bargmann et l'absence de dérive séculaire d'énergie.
3. **Alignment Criteria**:
   - Vérification des ordres formels de convergence numérique : $\mathcal{S}_2$ (ordre 2), $\mathcal{S}_4$ (ordre 4), $\mathrm{BM}_4^{(7)}$ (ordre 4 à erreur résiduelle minimale), $\mathrm{BM}_6^{(11)}$ (ordre 6), $\mathcal{P}\mathrm{BM}_4^{(6)}$ (ordre 6 effectif), $\mathcal{P}\mathrm{BM}_6^{(8)}$ (ordre 8 effectif).
   - Télescopage parfait des correcteurs aux frontières dyadiques : $\mathcal{C}^{-1} \circ \mathcal{C} = \operatorname{Id}$ sans appel supplémentaire au champ de force.
   - Invariance symplectique et réversibilité temporelle exacte $\Phi(-h) \circ \Phi(h) = \operatorname{Id}$.

## Guidelines:
- **Method**: Protocole strict à 5 étapes.
- **Red Lines**:
  - Sanctuarisation absolue en lecture seule du répertoire [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/).
  - Zéro occurrence de `"cas6"` ou `"cas 6"`.
  - Zéro adoucissement ni troncature artificielle ($\varepsilon = 0, d_0 = 0$).
  - Aucun couplage ascendant vers les couches 5 (`system/`), 6 (`diagnostics/`) ou 7 (`benchmarks/`). Dépendances strictes vers `core/`, `solvers/` et `potential/`.
- **Architecture**:
  - `integrators/__init__.py` : exports unifiés des intégrateurs et du protocole.
  - `integrators/protocol.py` : interface formelle `SymplecticIntegrator` (PEP 544), `IntegrationStepResult`, métadonnées de schéma (ordre effectif, nombre d'étages, symétrie palindromique).
  - `integrators/strang.py` : schéma de Strang $\mathcal{S}_2(h) = \Phi_h^B(\frac{h}{2}) \circ \Phi_h^A(h) \circ \Phi_h^B(\frac{h}{2})$ d'ordre 2.
  - `integrators/yoshida.py` : composition triple symétrique de Yoshida $\mathcal{S}_4(h)$ d'ordre 4.
  - `integrators/blanes_moan_4.py` : schéma optimisé $\mathrm{BM}_4^{(7)}$ à 7 étages pour hamiltoniens séparables.
  - `integrators/blanes_moan_6.py` : schéma $\mathrm{BM}_6^{(11)}$ à 11 étages avec coefficients certifiés à 25 décimales (`BM6_DRIFT_WEIGHTS`, `BM6_IMPULSE_WEIGHTS`).
  - `integrators/correctors.py` : correcteurs symplectiques canoniques $\mathcal{C}_{2p}$ (opérateurs de bordure, équation homologique, signe correct $\exp(+h^{2p} X_K)$).
  - `integrators/processed.py` : intégrateur traité $\mathcal{P}\mathrm{BM} = \mathcal{C} \circ \mathcal{K}^m \circ \mathcal{C}^{-1}$ élevant l'ordre effectif de 2 rangs à coût marginal nul.

## Plan:
### Step 1: Relecture Théorique du Chapitre 5 de la Monographie
- **Context**: Le Chapitre 5 formalise les compositions palindromiques, la théorie de Baker--Campbell--Hausdorff (BCH), le dépassement de la barrière de Suzuki et la théorie des intégrateurs traités de Blanes.
- **Objective**: Extraire les coefficients exacts, les équations des correcteurs, les poids d'étage et les théorèmes Lean 4 associés (`BM6Order6SymmetricComposition.lean`, `02_LieBargmannExactInvariance_Theorems.lean`).
- **Deliverables**:
  - `Layer4TheorySynthesis`: Synthèse des coefficients numériques et des formulations d'étage.
- **Process**: Relecture ciblée de `chapters/ch05_higher_order_symplectic/`.

### Step 2: Audit Approfondi du Code Hérité (`cas6_nbody3d_async/integrators/`)
- **Context**: L'ancien code contient les embryons de Strang, Yoshida, BM4 et des correcteurs, mais manque de `BM_6_11` autonome et de `ProcessedSymplecticIntegrator`.
- **Objective**: Recenser les composants, purifier les 38 mentions `cas6`, et concevoir de novo les schémas d'ordre 6 et traités.
- **Deliverables**:
  - `Layer4RefactoringMatrix`: Matrice de migration et de spécification des nouveautés.
- **Process**: Confrontation avec `python_codebase_mapping_and_audit_report.md`.

### Step 3: Spécification Formelle et Contrats d'Interface de la Couche 4
- **Context**: Définir les types, méthodes de pas (`step`), gestion des super-nœuds et de l'arbre de Jacobi.
- **Objective**: Établir le contrat d'interface `SymplecticIntegrator` et spécifier les 7 modules de test unitaire.
- **Deliverables**:
  - `Layer4ModuleSpecification`: Spécifications logicielles et matrices de test.
- **Process**: Définition des signatures `(state, tree, dt, apex, supernode_bodies) -> ExtendedPhaseState`.

### Step 4: Implémentation Modulaire et Optimisée des 7 Modules de `integrators/`
- **Context**: Créer le package `normal_cone_nbody_float64/integrators/`.
- **Objective**: Écrire le code avec vectorisation complète, réutilisation des structures de mémoire et en-têtes certifiés.
- **Deliverables**:
  - `normal_cone_nbody_float64/integrators/__init__.py`
  - `normal_cone_nbody_float64/integrators/protocol.py`
  - `normal_cone_nbody_float64/integrators/strang.py`
  - `normal_cone_nbody_float64/integrators/yoshida.py`
  - `normal_cone_nbody_float64/integrators/blanes_moan_4.py`
  - `normal_cone_nbody_float64/integrators/blanes_moan_6.py`
  - `normal_cone_nbody_float64/integrators/correctors.py`
  - `normal_cone_nbody_float64/integrators/processed.py`
- **Process**: Implémentation séquentielle optimisée.

### Step 5: Implémentation et Exécution des Tests Unitaires Miroirs (`tests/unit/integrators/`)
- **Context**: 7 fichiers de tests unitaires pour valider chaque intégrateur individuellement.
- **Objective**: Exécuter `pytest tests/unit/integrators/` avec 100 % de succès.
- **Deliverables**:
  - `tests/unit/integrators/test_protocol.py`
  - `tests/unit/integrators/test_strang.py`
  - `tests/unit/integrators/test_yoshida.py`
  - `tests/unit/integrators/test_blanes_moan_4.py`
  - `tests/unit/integrators/test_blanes_moan_6.py`
  - `tests/unit/integrators/test_correctors.py`
  - `tests/unit/integrators/test_processed.py`
- **Process**: Validation des taux de convergence par régression log-log de l'erreur d'énergie $\Delta E(h) \propto h^p$.

## Acceptance Criteria:
### Step 1:
- Relecture complète du Chapitre 5 et extraction des coefficients certifiés à 25 décimales.
### Step 2:
- Rapport d'audit attestant la détection et la suppression de toutes les références `cas6`.
### Step 3:
- Protocole formel `SymplecticIntegrator` validant le typage statique mypy/pyright.
### Step 4:
- 7 modules d'intégrateurs opérationnels avec zéro adoucissement ($\varepsilon = 0, d_0 = 0$).
### Step 5:
- 100 % de tests passés sur `pytest tests/unit/integrators/`.
- Pente d'ordre 2 pour Strang, 4 pour Yoshida et BM4, 6 pour BM6, et gain d'ordre sur Processed.
- Télescopage dyadique des correcteurs validé numériquement ($\|\mathcal{C}^{-1}(\mathcal{C}(z)) - z\| < 10^{-14}$).
- Réversibilité temporelle exacte sur tous les schémas palindromiques.

## Observations:
> [!IMPORTANT]
> L'inversion du signe de l'opérateur de correction $\exp(+h^{2p} X_K)$ découle de la parité paire du terme d'erreur d'ordre $2p$ dans l'expansion BCH symétrique, comme identifié lors de l'audit théorique de la monographie LaTeX.

> [!TIP]
> L'architecture des intégrateurs sépare rigoureusement la propagation analytique des paires feuilles (via `kepler_sundman.py`), le kick d'apex (via `screened_apex.py`), et les super-nœuds synchrones (via `supernode.py`), permettant une composition modulaire parfaite dans chaque sous-étape de drift/kick.
