# Implementation Plan of Iteration 169 : Instrumentation du Profiler Fin des Coûts Calculatoires, Visualisation Comparative et Diagnostic sur Cas Référent
Created at: 2026-10-09T11:55:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur demande de concevoir un code permettant de choisir une situation initiale et un intégrateur, de calculer et comptabiliser finement l'ensemble des étapes de calcul et leurs temps d'exécution, puis de générer une planche de diagnostic visuel pour chaque solveur permettant d'identifier immédiatement les étapes les plus consommatrices selon la dynamique physique. Un cas référent simple avec BM6 doit être exécuté immédiatement pour afficher la planche résultante.
2. **Implicit Strategic Value**: Remplacer les métriques de travail théoriques ou partielles par une traçabilité instrumentée rigoureuse du coût calculatoire réel. Cela permet de révéler les goulets d'étranglement algorithmiques (notamment le temps passé dans les sous-pas Leapfrog des super-nœuds en Python pur vs kick d'apex vs solveur de Kepler), de guider de futures accélérations (vectorisation / Numba), et de fonder le diagramme de Pareto sur des mesures empiriques irréfutables.
3. **Alignment Criteria**: L'instrumentation doit avoir un overhead négligeable en production (< 2% de temps supplémentaire), ne pas altérer la conservation des invariants de Lie-Bargmann ni l'ordre de convergence, et fournir une planche analytique haute résolution (300 DPI) immédiatement lisible et didactique.

## Guidelines:
- **Method**: Découper le coût en compteurs physiques et algorithmiques canoniques : évaluations de forces d'apex, itérations képlériennes de Halley/Stumpff, bissections, sous-pas Leapfrog de super-nœuds, calculs d'échelles dyadiques, tests d'atlas, et dérivées directionnelles hessiennes.
- **Red Lines**: Préservation inconditionnelle du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$). Aucune modification intrusive ne doit rompre les propriétés symplectiques ou la réversibilité temporelle microscopique.
- **Architecture**: Profiler modulaire non-intrusif encapsulé dans `normal_cone_nbody_float64/diagnostics/profiler.py`, branché de façon optionnelle sur `NBodySystem`, `HistoryTracker` et les sous-solveurs.

## Plan:
### Step 1: Conception et implémentation du profiler granulaire dans `diagnostics/profiler.py`
- **Context**: Le moteur conormal exécute plusieurs étages à chaque sous-pas dyadique sans instrumentation unifiée de temps CPU par étage ni de compteurs fins d'opérations.
- **Objective**: Construire un modèle de données et un chronomètre de haute précision (`time.perf_counter()`) avec compteurs d'opérations atomiques.
- **Deliverables**:
  - `ComputationalProfiler`: Classe de profilage dotée de context managers légers (`profile_stage("apex")`, `profile_stage("kepler")`, etc.) et de compteurs d'événements.
  - `StageCostMetrics`: Dataclass représentant le temps CPU cumulé et le décompte d'opérations pour chaque composant.
- **Process**: Implémenter le module `profiler.py` avec méthodes d'accumulation rapide et export vers dictionnaire/JSON.

### Step 2: Câblage non-intrusif du profiler dans le moteur conormal (`nbody_system.py`, `solvers/`, `potential/`)
- **Context**: `NBodySystem.step()` orchestre les étapes de calcul sans propager d'objet profiler aux sous-routines.
- **Objective**: Permettre à `NBodySystem` d'accepter un argument `profiler: Optional[ComputationalProfiler] = None` et de le transmettre aux sous-modules (`step_drift_and_kepler`, `compute_apex_forces`, `step_supernode_synchronous`, etc.).
- **Deliverables**:
  - `NBodySystem`: Câblage du profiler optionnel dans `step()` et `evolve()`.
  - `HistoryTracker`: Enregistrement des métriques de coût par macro-pas.
- **Process**: Ajouter les hooks de profilage aux frontières logiques des fonctions cibles sans modifier leur logique mathématique.

### Step 3: Conception de la planche de diagnostic multi-métriques dans `diagnostics/visualizer.py`
- **Context**: L'utilisateur doit pouvoir visualiser d'un coup d'œil où partent les ressources de calcul.
- **Objective**: Concevoir une planche graphique haute résolution (300 DPI) décomposant les coûts sous plusieurs angles complémentaires.
- **Deliverables**:
  - `render_computational_cost_dashboard`: Fonction de visualisation Matplotlib 300 DPI avec 4 sous-panneaux :
    1. Répartition globale du temps CPU (Camembert / barres empilées par étage : Kepler, Apex, Supernœuds, Dyadique, Atlas, Correcteurs).
    2. Dynamique temporelle du coût par pas vs configuration céleste ($r_{\min}(t)$ et niveau dyadique $k(t)$).
    3. Histogramme des opérations élémentaires (itérations Halley, sous-pas Leapfrog, kicks d'apex).
    4. Compromis Coût réel vs Précision de l'énergie.
- **Process**: Implémenter la fonction dans `visualizer.py` en réutilisant la charte graphique scientifique de la monographie.

### Step 4: Création de l'outil CLI interactif/paramétrable `scratch/run_cost_profiler.py`
- **Context**: L'utilisateur souhaite pouvoir choisir une situation initiale et un intégrateur, puis lancer le calcul et obtenir la planche.
- **Objective**: Fournir un script exécutable en ligne de commande ou appelable programmaticement avec sélection simple du cas et du solveur.
- **Deliverables**:
  - `scratch/run_cost_profiler.py`: Script Python avec arguments `--case` (`figure8`, `burrau`, `radial_collision`, `hierarchical_4body`), `--integrator` (`strang`, `yoshida`, `bm4`, `bm6`, `pbm4`, `pbm6`), `--dt`, `--t-span` et `--output`.
- **Process**: Développer le script en s'appuyant sur le registre des benchmarks existant (`benchmarks/registry.py`).

### Step 5: Exécution du diagnostic de référence sur un cas simple avec BM6, validation et clôture git
- **Context**: L'utilisateur demande de tester immédiatement sur un cas initial relativement simple avec BM6 et de lui montrer la planche.
- **Objective**: Exécuter le script sur le problème à 3 corps de Burrau (ou Figure-8) avec BM6, générer la planche 300 DPI, valider l'intégrité globale et consigner les modifications.
- **Deliverables**:
  - `cost_breakdown_burrau_bm6_300dpi.png`: Planche diagnostique générée et copiée dans `brain/`.
  - Suite de tests unitaires : 151/151 tests passants.
  - Commit git sur la branche `main`.
- **Process**: Lancer le script via `run_command`, copier l'image générée dans `brain/`, exécuter `pytest`, et enregistrer le commit git.

## Acceptance Criteria:
### Step 1:
- Le module `profiler.py` fournit des mesures de temps CPU avec une résolution nanoseconde via `perf_counter` et des compteurs cumulatifs exacts.
### Step 2:
- `NBodySystem` s'exécute de manière identique avec ou sans profiler (`profiler=None` n'induit aucun changement de comportement ni dérive numérique).
### Step 3:
- La planche générée est enregistrée à 300 DPI, avec des axes clairs, des pourcentages lisibles et une légende détaillée.
### Step 4:
- Le script CLI prend en charge l'ensemble des 4 cas de référence et des 6 schémas d'intégration.
### Step 5:
- La planche pour le cas initial avec BM6 est générée avec succès dans `brain/`.
- 151/151 tests unitaires passent sans régression.
- Le commit git est consigné sur la branche `main`.

## Observations:
> [!NOTE]
> Le problème à trois corps de Burrau (rencontre serrée avec quasi-chocs et éjection) est particulièrement adapté pour ce benchmark initial car il traverse à la fois des régimes képlériens réguliers, des montées en niveau dyadique et des activations de super-nœuds, permettant d'observer la variation des coûts en fonction de la géométrie de phase.
