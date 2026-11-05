# Implementation Plan of Iteration 173 : Propagation du Flag d'Optimisation Kepler vers les Super-Nœuds Synchrones
Created at: 2026-10-09T12:32:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur commande de propager le flag d'optimisation du solveur Kepler vers les super-nœuds (« Commence par faire en sorte que si le flag d'optimiseur du solveur kepler est actif il se propage aux super noeuds »).
2. **Implicit Strategic Value**: Éliminer le goulet d'étranglement des 68 073 régularisations de paires serrées internes résolues par le solveur de base au sein des super-nœuds, en réutilisant l'inlining de Halley et le micro-noyau scalaire de Stumpff sans modifier la physique exacte ni rompre la rétro-compatibilité.
3. **Alignment Criteria**: Le solveur de super-nœuds (`solvers/supernode.py`) doit détecter automatiquement si un solveur Kepler optimisé est configuré (via le `ContextVar` ou argument explicite), déléguer l'intégration de la paire serrée à `step_kepler_sundman_fast`, préserver une concordance bit-exacte ($< 10^{-14}$) avec l'implémentation de base, et valider 100 % des tests unitaires existants et nouveaux.

## Guidelines:
- **Method**: Utiliser la détection non intrusive via `ContextVar` (`get_current_kepler_solver()`) et le paramètre optionnel `kepler_solver` dans la chaîne d'appel `step_supernode_flow` $\to$ `step_supernode_synchronous` $\to$ `_step_cluster_internal_substep`.
- **Red Lines**: Préservation inconditionnelle du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$), des 10 invariants du groupe de Lie-Bargmann, de l'unimodularité de Battin et de l'équivalence numérique des trajectoires.
- **Architecture**: Découplage propre, tolérance totale aux exécutions autonomes sans optimiseur (fallback baseline garanti).

## Plan:
### Step 1: Implémentation de la propagation du solveur Kepler dans `solvers/supernode.py`
- **Context**: `solvers/supernode.py` importe actuellement en dur `from .kepler_sundman import step_kepler_sundman` pour régulariser les paires serrées ($r_{\min} < 0.5 r_{\text{second}}$).
- **Objective**: Permettre la résolution dynamique de la routine Kepler active via `_resolve_kepler_step` tenant compte du flag / `ContextVar` d'optimisation.
- **Deliverables**:
  - `normal_cone_nbody_float64/solvers/supernode.py`: Prise en charge de `kepler_solver` et propagation à `_step_cluster_internal_substep`, au cas $k=2$ et à la configuration homothétique.
- **Process**:
  - Ajouter la fonction de résolution dynamique `_resolve_kepler_step(kepler_solver=None)`.
  - Mettre à jour `_step_cluster_internal_substep`, `step_supernode_synchronous` pour recevoir et propager le solveur actif.

### Step 2: Raccordement dans `integrators/common.py`
- **Context**: `integrators/common.py` orchestre les étapes de flux (`step_supernode_kick`, `step_supernode_flow`, `step_drift_and_kepler`).
- **Objective**: Transmettre le solveur Kepler actif (issu du `ContextVar` ou passé en argument) aux fonctions de flux de super-nœuds.
- **Deliverables**:
  - `normal_cone_nbody_float64/integrators/common.py`: Transmission du `kepler_solver` dans `step_supernode_kick` et `step_supernode_flow`.
- **Process**:
  - Récupérer `kepler_solver` via `get_current_kepler_solver()` si non fourni explicitement.
  - Passer l'argument à `step_supernode_synchronous`.

### Step 3: Tests unitaires de propagation et non-régression
- **Context**: Garantir que l'activation du solveur optimisé produit un résultat bit-exact tout en accélérant les micro-pas de paires serrées.
- **Objective**: Écrire une suite de tests unitaires dédiée pour attester de la propagation, de la rétro-compatibilité et de la fidélité symplectique.
- **Deliverables**:
  - `tests/unit/solvers/test_supernode_optimized_kepler.py`: Tests unitaires exhaustifs.
- **Process**:
  - Tester l'exécution sous baseline (fallback par défaut).
  - Tester l'exécution sous `OptimizedKeplerSolver` via ContextVar et argument direct.
  - Vérifier l'erreur relative $< 10^{-14}$ et la stricte conservation du moment linéaire.

### Step 4: Évaluation d'impact sur le banc 4-corps hiérarchique et validation complète
- **Context**: Le benchmark 4-corps a révélé 68 073 régularisations de paires serrées au sein des super-nœuds.
- **Objective**: Mesurer l'effet de cette première optimisation sur le temps de calcul des super-nœuds et vérifier l'ensemble de la suite de tests.
- **Deliverables**:
  - `BenchmarkMetrics`: Mesure comparative du temps CPU des super-nœuds.
  - `TestValidation`: Validation de 100 % des tests unitaires (`pytest tests/unit/`).
- **Process**:
  - Lancer un profilage comparatif ciblé sur le 4-corps hiérarchique.
  - Exécuter la suite complète `pytest tests/unit/`.

## Acceptance Criteria:
### Step 1:
- `solvers/supernode.py` résout `step_kepler_sundman_fast` lorsque l'optimiseur Kepler est actif, et conserve `step_kepler_sundman` en mode baseline.
### Step 2:
- `integrators/common.py` propage de manière transparente le solveur Kepler actif vers tous les appels de super-nœuds.
### Step 3:
- Les nouveaux tests unitaires dans `test_supernode_optimized_kepler.py` passent à 100 % avec une dérive $< 10^{-14}$.
### Step 4:
- 100 % des tests de la suite complète (167+ tests) passent sans régression.
- Le temps passé sur les régularisations de paires serrées dans les super-nœuds est réduit.

## Observations:
> [!NOTE]
> Cette étape est la première phase de l'optimisation des super-nœuds, ciblant spécifiquement la délégation aux paires képlériennes serrées (68 073 appels).
> Elle n'introduit aucun nouveau fichier de calcul lourd mais fluidifie le couplage entre les deux solveurs différentiels conormaux.
