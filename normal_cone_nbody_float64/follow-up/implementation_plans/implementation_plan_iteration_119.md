# Implementation Plan of Iteration 119 : Couche 7 & Certification Globale (`benchmarks/`) - Les 8 Problèmes Canoniques et Courbes Pareto Comparatives
Created at: 2026-10-07T18:10:00+02:00

## Context:
1. **Explicit Request**: Développer, tester unitairement et certifier globalement la Couche 7 (`benchmarks/`) de la bibliothèque `normal_cone_nbody_float64` en suivant strictement le protocole mandaté en 5 étapes.
2. **Implicit Strategic Value**: Livrer le banc d'essai expérimental complet et les cas d'usage canoniques de la monographie de recherche, confrontant les 6 intégrateurs symplectiques ($\mathcal{S}_2, \mathcal{S}_4, \mathrm{BM}_4^{(7)}, \mathrm{BM}_6^{(11)}, \mathcal{P}\mathrm{BM}_4^{(6)}, \mathcal{P}\mathrm{BM}_6^{(8)}$) sur les régimes gravitationnels les plus sévères (collisions radiales, quasi-collisions de Burrau, chorégraphies en huit, orbites d'Arenstorf, système planétaire), et démontrant la supériorité numérique prouvée de l'intégration conormale sans adoucissement.
3. **Alignment Criteria**: Taux de réussite de 100 % sur l'ensemble de la suite de tests unitaires et de benchmarks, absence totale de marqueurs `"cas6"`, zéro adoucissement ($\varepsilon = 0, d_0 = 0$), génération des courbes Pareto travail-précision confirmant les ordres de convergence effectifs (2, 4, 6, 8), et certification globale finale de la bibliothèque.

## Guidelines:
- **Method**: Protocole rigoureux en 5 étapes (Théorie -> Audit hérité -> Spécification formelle -> Implémentation modulaire -> Tests unitaires 1:1).
- **Red Lines**: Interdiction absolue de toute modification dans `redaction/unsoftened_conormal_nbody3d/`. Éradication totale de `"cas6"` ou `"cas 6"`. Zéro adoucissement ni seuillage artificiel de distance ($\varepsilon = 0, d_0 = 0$).
- **Architecture**:
  - `benchmarks/registry.py` : Usine formelle des 8 configurations astrophysiques canoniques :
    1. Collision frontale radiale $L=0$ (franchissement régularisé de singularité).
    2. Problème de Pythagore de Burrau (masses 3-4-5, quasi-collision triple extrême).
    3. Orbite lunaire périodique d'Arenstorf (PCR3BP avec passage rasant).
    4. Orbite chorégraphique en huit de Chenciner--Montgomery (masses égales).
    5. Problème spatial oscillatoire de Sitnikov (symétrie $Z$).
    6. Collision quadruple symétrique ($N=4$).
    7. Système jovien extérieur (Soleil--Jupiter--Saturne, dynamique quasi-séculaire).
    8. Configuration d'évaluation de convergence Pareto travail-précision.
  - `benchmarks/runner.py` : Moteur d'exécution modulaire `BenchmarkRunner` paramétrable par intégrateur, durée, pas et métrologie KPI, avec calcul automatique des courbes travail-précision (évaluations de force vs dérive hamiltonienne).
  - `benchmarks/__init__.py` : Exports du package.
- **Form**: Code, commentaires, docstrings et tests exclusivement en anglais avec citations de la monographie et des théorèmes Lean 4. Plan et liste de tâches en français markdown natif.

## Plan:
### Step 1: Relecture Théorique des Chapitres 6 et 7 de la Monographie
- **Starting Material**: Monographie LaTeX `redaction/unsoftened_conormal_nbody3d/` (Chapitres 6 et 7 sur les cas d'essais numériques et l'analyse comparative des schémas symplectiques).
- **Objective**: Extraire les conditions initiales de haute précision (16 décimales), les périodes analytiques, les intégrales premières théoriques et les invariants de chaque problème canonique.
- **Deliverables**: Spécification théorique consolidée pour la Couche 7 (`Component Theoretical Framework`).
- **Process**: Relecture méthodique des fichiers `.tex` et formalisation des données orbitales exactes.

### Step 2: Audit Approfondi du Code Hérité (`cas6_nbody3d_async/benchmarks/`)
- **Starting Material**: Répertoire hérité `models/normal_cone_integrator/cas6_nbody3d_async/benchmarks/` (notamment `benchmark1_collision.py` à `benchmark7_quadruple_collision.py`, `baselines.py`), et le rapport d'audit `python_codebase_mapping_and_audit_report.md`.
- **Objective**: Recenser les paramètres numériques éprouvés, éliminer les 48 occurrences résiduelles de `cas6`, standardiser les retours d'exécution sous `ExtendedPhaseState`, et modulariser le runner de simulation.
- **Deliverables**: Rapport de refactorisation Couche 7 (`Component Legacy Audit`).
- **Process**: Confrontation ligne à ligne et purification des scripts historiques vers des fonctions réutilisables.

### Step 3: Spécification Formelle et Contrats d'Interface de la Couche 7
- **Starting Material**: Spécification globale `normal_cone_nbody_float64_master_architecture.md` (Section 3.7 & Section 7).
- **Objective**: Définir les types, signatures, exceptions, invariants et critères d'évaluation des modules de `benchmarks/` et de leurs tests unitaires miroirs.
- **Deliverables**: Contrats d'interface du banc d'essai (`Component Interface Contracts`).
- **Process**: Formalisation des classes `BenchmarkRegistry`, `BenchmarkRunner`, `BenchmarkResult`, `ParetoCurve`.

### Step 4: Implémentation Modulaire et Optimisée des Modules de `benchmarks/`
- **Starting Material**: Modules certifiés des Couches 1 à 6 (`core/`, `solvers/`, `potential/`, `integrators/`, `system/`, `diagnostics/`).
- **Objective**: Implémenter les 3 fichiers sous `normal_cone_nbody_float64/benchmarks/` :
  1. `benchmarks/__init__.py` : Exports du package.
  2. `benchmarks/registry.py` : Constructeurs des 8 configurations canoniques avec masses, positions, impulsions et temps caractéristiques.
  3. `benchmarks/runner.py` : Exécuteur `BenchmarkRunner` et fonction `compute_work_precision_pareto`.
- **Deliverables**: Fichiers de production sous `normal_cone_nbody_float64/benchmarks/`.
- **Process**: Écriture modulaire avec typage strict PEP 484/544, zéro allocation superflue et documentation académique complète.

### Step 5: Implémentation et Exécution des Tests Unitaires Miroirs (`tests/unit/benchmarks/`)
- **Starting Material**: Spécifications d'évaluation et modules de `benchmarks/`.
- **Objective**: Développer et exécuter avec 100 % de succès les 2 bancs de test unitaires miroirs :
  1. `tests/unit/benchmarks/test_registry.py` (vérification dimensionnelle, invariances de centre de masse et énergies théoriques des 8 configurations).
  2. `tests/unit/benchmarks/test_runner.py` (exécution de simulations, calcul de Pareto travail-précision, validation de la hiérarchie d'ordre $\mathrm{BM}_6 > \mathrm{BM}_4 > \mathcal{S}_2$).
- **Deliverables**: 2 fichiers de test sous `tests/unit/benchmarks/` avec 100 % de réussite.
- **Process**: Exécution via `pytest -v`, vérification des tolérances et du comportement numérique.

## Acceptance Criteria:
- **Critère 7.1 (Complétude du Registre)** : Les 8 problèmes canoniques sont instanciés avec dimensions exactes, masses positives, et barycentre initial rigoureusement calibré ($\|\mathbf{P}_0\| < 10^{-14}, \|\mathbf{R}_{cm}\| < 10^{-14}$).
- **Critère 7.2 (Robustesse du Runner)** : `BenchmarkRunner.run()` simule avec succès n'importe quelle configuration sur l'ensemble des intégrateurs disponibles et produit un `BenchmarkResult` complet avec rapport KPI.
- **Critère 7.3 (Courbes Pareto Travail-Précision)** : Sur la configuration de test, à nombre d'évaluations de force égal, $\mathrm{BM}_6$ atteint une précision hamiltonienne supérieure d'au moins 2 ordres de grandeur à $\mathcal{S}_2$, confirmant la pente de convergence asymptotique d'ordre 6.
- **Critère 7.4 (Zero-Softening & Zero-Cas6)** : 0 occurrence de `"cas6"` ou `"cas 6"`, 0 softening ($\varepsilon = 0, d_0 = 0$).
- **Critère 7.5 (Taux de Succès Global)** : 100 % de réussite sur l'ensemble de la suite de tests unitaires cumulative de la bibliothèque (115+ tests passés).

## Observations:
> [!IMPORTANT]
> Les benchmarks doivent utiliser des pas de temps adaptés à leurs échelles dynamiques spécifiques pour illustrer la robustesse de l'intégration sans singularité.
> [!NOTE]
> Le runner de benchmark utilise la façade `NBodySystem` certifiée à l'Itération 117 et les métriques de l'Itération 118.
