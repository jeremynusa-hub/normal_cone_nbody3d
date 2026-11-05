# Implementation Plan of Iteration 172 : Génération des Planches de Diagnostic de Coût Profilé (Figure-8, Collision Radiale, 4-Corps Hiérarchique) avec Double Optimisation BM6 (300 DPI)
Created at: 2026-10-09T12:24:10+02:00

## Context:
1. **Explicit Request**: L'utilisateur commande de générer les planches de diagnostic de coût computationnel pour les autres simulations canoniques avec les deux optimisations actives (Kick d'Apex optimisé + Solveur de Kepler optimisé) sous l'intégrateur symplectique d'ordre 6 Blanes-Moan BM6 (« Ok, avec les 2 opti, et BM6, fait les planche des autres simulations »).
2. **Implicit Strategic Value**: Cartographier et caractériser avec précision la signature computationnelle de l'architecture conormale unsoftened sur les différentes classes de dynamiques célestes (régime périodique quasi-circulaire de la Figure-8, singularité de choc frontal à 2 corps, et système hiérarchique multirate à 4 corps), attestant de la robustesse, du débit et de la conservation stricte des invariants physiques sous double optimisation.
3. **Alignment Criteria**: Produire trois planches de diagnostic à 300 DPI distinctes et lisibles dans `brain/`, accompagnées d'un tableau synthétique comparant les temps d'exécution, la répartition des postes de coût, le débit (pas/s), et les erreurs d'énergie, tout en certifiant 100 % de réussite de la suite de tests unitaires.

## Guidelines:
- **Method**: Utiliser l'outil d'instrumentation et de profilage `scratch/run_cost_profiler.py` avec `--apex-variant optimized --kepler-variant optimized --integrator bm6`.
- **Red Lines**: Préservation inconditionnelle du zéro-adoucissement ($\varepsilon = 0, d_0 = 0$), des 10 invariants du groupe de Lie-Bargmann, et de l'intégrité du code cœur.
- **Architecture**: Conservation de la double optimisation modulaire découplée (`OptimizedScreenedApexPotential` et `OptimizedKeplerSolver`).

## Plan:
### Step 1: Validation du paramétrage des trois scénarios dans `scratch/run_cost_profiler.py`
- **Context**: Les scénarios Figure-8, Collision Radiale et 4-Corps Hiérarchique possèdent des durées caractéristiques $t_{\text{span}}$ et pas $\Delta t$ adaptés.
- **Objective**: Vérifier que les paramètres par défaut de chaque scénario sont calibrés pour un profilage représentatif et stable sous BM6.
- **Deliverables**:
  - `ConfigurationAudit`: Validation des paramètres de simulation dans `scratch/run_cost_profiler.py`.
- **Process**:
  - Inspecter les définitions de `figure8` ($t_{\text{span}} \approx 6.3259, \Delta t = 0.01$), `radial_collision` ($t_{\text{span}} \approx 1.1107, \Delta t = 0.005$), et `hierarchical_4body` ($t_{\text{span}} = 20.0, \Delta t = 0.01$).

### Step 2: Exécution et rendu de la planche pour la Chorégraphie Figure-8
- **Context**: Problème à 3 corps de masse égale sur une trajectoire en 8 hautement intriquée sans collision binaire serrée.
- **Objective**: Mesurer le coût de l'apex et de Kepler en régime chorégraphique périodique et générer la planche 300 DPI.
- **Deliverables**:
  - `cost_breakdown_figure8_bm6_fully_optimized_300dpi.png`: Planche 300 DPI dans `brain/`.
- **Process**:
  - Exécuter `python scratch/run_cost_profiler.py --case figure8 --integrator bm6 --apex-variant optimized --kepler-variant optimized`.
  - Vérifier la génération et inspecter la planche.

### Step 3: Exécution et rendu de la planche pour la Collision Radiale
- **Context**: Problème à 2 corps plongeant directement dans la singularité géométrique ($r \to 0, e = 1$).
- **Objective**: Mesurer l'activation du régime de choc frontal et l'isolation conormale de Kepler-Sundman, et générer la planche 300 DPI.
- **Deliverables**:
  - `cost_breakdown_radial_bm6_fully_optimized_300dpi.png`: Planche 300 DPI dans `brain/`.
- **Process**:
  - Exécuter `python scratch/run_cost_profiler.py --case radial_collision --integrator bm6 --apex-variant optimized --kepler-variant optimized`.
  - Vérifier la traversée de singularité et inspecter la planche.

### Step 4: Exécution et rendu de la planche pour le 4-Corps Hiérarchique
- **Context**: Système à 4 corps multi-échelles combinant binaire serrée, hiérarchie composite et super-nœuds synchrones.
- **Objective**: Mesurer l'interaction entre super-nœuds, ordonnanceur dyadique multirate et double optimisation, et générer la planche 300 DPI.
- **Deliverables**:
  - `cost_breakdown_4body_bm6_fully_optimized_300dpi.png`: Planche 300 DPI dans `brain/`.
- **Process**:
  - Exécuter `python scratch/run_cost_profiler.py --case hierarchical_4body --integrator bm6 --apex-variant optimized --kepler-variant optimized`.
  - Vérifier la régularisation synchrone et inspecter la planche.

### Step 5: Synthèse comparative, validation des tests unitaires et enregistrement Git
- **Context**: Synthèse globale des métriques sur les 4 classes de problèmes (avec Burrau).
- **Objective**: Construire le tableau croisé de synthèse des coûts, valider 100 % de la suite de tests unitaires et consigner le commit.
- **Deliverables**:
  - `ComparativeSynthesis`: Tableau croisé des métriques.
  - `GitCommit`: Commit sur `main`.
- **Process**:
  - Exécuter `pytest tests/unit/`.
  - Consigner les modifications de scripts et tests sous commit Git.

## Acceptance Criteria:
### Step 1:
- Les trois scénarios s'exécutent sans divergence ni exception de typage.
### Step 2:
- La planche `cost_breakdown_figure8_bm6_fully_optimized_300dpi.png` est générée à 300 DPI et lisible dans `brain/`.
- L'erreur d'énergie sur la Figure-8 reste confinée à $< 10^{-10}$.
### Step 3:
- La planche `cost_breakdown_radial_bm6_fully_optimized_300dpi.png` est générée à 300 DPI et lisible dans `brain/`.
- Le choc frontal $r \to 0$ est régularisé sans singularité par le solveur conormal.
### Step 4:
- La planche `cost_breakdown_4body_bm6_fully_optimized_300dpi.png` est générée à 300 DPI et lisible dans `brain/`.
- L'ordonnancement dyadique et les super-nœuds synchrones s'activent de façon stable.
### Step 5:
- 100 % des tests unitaires sont validés avec succès (`pytest tests/unit/`).
- Commit consigné sur `main`.

## Observations:
> [!IMPORTANT]
> Chaque simulation opère sous physique conormale strictly unsoftened ($\varepsilon = 0, d_0 = 0$).
> Les artefacts générés fournissent une vision comparative complète du profil de coût à travers toutes les typologies dynamiques du papier de recherche.
