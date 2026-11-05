# Implementation Plan of Iteration 174 : Solveur de Super-Nœuds Synchrones Découplé et Micro-Noyau Vectorisé sans Allocation
Created at: 2026-10-09T12:43:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur commande de réaliser l'optimisation des super-nœuds synchrones, tout en conservant une architecture strictement modulaire et optionnelle (« Ok fait le, en modulaire optionnel toujours »).
2. **Implicit Strategic Value**: Résoudre le goulet d'étranglement majeur de 34,4 secondes (88,4 % du temps CPU sur le 4-corps hiérarchique) généré par 306 537 évaluations d'accélérations internes et les allocations répétées de tableaux NumPy 3D dans 204 358 sous-pas de Leapfrog, sans altérer le code cœur historique ni rompre la rétro-compatibilité.
3. **Alignment Criteria**: Produire un module découplé `OptimizedSupernodeSolver` éliminant les allocations en mémoire dans la boucle microscopique, réutilisant les calculs de distances entre l'adaptation temporelle et les sous-pas, inlinant les forces lisses des spectateurs, préservant une concordance bit-exacte ($< 10^{-13}$) avec le solveur de référence, et démontrant une accélération substantielle sur le banc à 4 corps avec 100 % de tests unitaires validés.

## Guidelines:
- **Method**: Concevoir `normal_cone_nbody_float64/solvers/optimized_supernode.py` comme une variante optionnelle interchangeable, activable via `supernode_evaluator` dans `NBodySystem` ou par `ContextVar` dans `integrators/common.py`.
- **Red Lines**: Préservation inconditionnelle du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$), des 10 invariants du groupe de Lie-Bargmann (quantité de mouvement, moment cinétique, centre de masse), et de la réversibilité temporelle microscopique.
- **Architecture**: Tampons pré-alloués (zero-allocation micro-kernel), micro-noyau spécialisé scalaire pour petits clusters ($k \le 4$), pas de fermetures Python recréées dans la boucle, et intégration directe de `step_kepler_sundman_fast`.

## Plan:
### Step 1: Conception et implémentation du module `OptimizedSupernodeSolver`
- **Context**: `solvers/supernode.py` alloue des tenseurs 3D `dr[k,k,3]`, des copies et des fermetures Python à chaque micro-pas de Leapfrog.
- **Objective**: Implémenter `normal_cone_nbody_float64/solvers/optimized_supernode.py` avec tampons pré-alloués, micro-noyau sans allocation pour $k \le 4$, et réutilisation directe des distances calculées pour le pas adaptatif.
- **Deliverables**:
  - `normal_cone_nbody_float64/solvers/optimized_supernode.py`: Classe `OptimizedSupernodeSolver` et micro-noyau `step_supernode_synchronous_fast`.
  - `normal_cone_nbody_float64/solvers/__init__.py`: Export propre du solveur optimisé.
- **Process**:
  - Écrire un micro-noyau d'accélération sans allocation avec tampons pré-alloués et spécialisation $k=3, 4$.
  - Éliminer la fermeture `_compute_smooth_forces` en inlinant la contribution spectatorielle.
  - Fusionner l'évaluation géométrique de l'adaptation temporelle avec le premier demi-kick de Leapfrog.

### Step 2: Branchement modulaire optionnel dans `integrators/common.py` et `system/nbody_system.py`
- **Context**: Le solveur de super-nœuds optimisé doit pouvoir être activé de manière non intrusive, exactement comme `OptimizedScreenedApexPotential` et `OptimizedKeplerSolver`.
- **Objective**: Mettre en place la délégation par `ContextVar` `_CURRENT_SUPERNODE_SOLVER` dans `integrators/common.py` et le support du paramètre `supernode_evaluator` dans `NBodySystem`.
- **Deliverables**:
  - `normal_cone_nbody_float64/integrators/common.py`: Fonctions `set_current_supernode_solver` et `get_current_supernode_solver`.
  - `normal_cone_nbody_float64/system/nbody_system.py`: Prise en charge de `supernode_evaluator` et activation du token.
- **Process**:
  - Définir `_CURRENT_SUPERNODE_SOLVER` ContextVar dans `common.py`.
  - Dériver `step_supernode_synchronous` vers le solveur optimisé s'il est actif.
  - Ajouter `supernode_evaluator` dans `NBodySystem.__init__` et gérer le token dans `step()`.

### Step 3: Tests unitaires d'isomorphisme et de non-régression
- **Context**: Le nouveau solveur optimisé doit être validé de façon exhaustive contre le solveur canonique.
- **Objective**: Rédiger une suite complète de tests unitaires attestant de la concordance numérique, de la conservation du centre de masse et du moment linéaire, et de l'interchangeabilité.
- **Deliverables**:
  - `tests/unit/solvers/test_optimized_supernode.py`: Tests unitaires exhaustifs.
- **Process**:
  - Tester les clusters binaires ($k=2$).
  - Tester les clusters à 3 corps et 4 corps avec et sans paires serrées.
  - Vérifier la conservation de $\sum m_i v_i = 0$ à $< 10^{-15}$ et l'équivalence avec la baseline à $< 10^{-12}$.

### Step 4: Intégration CLI et benchmark comparatif sur le 4-Corps Hiérarchique
- **Context**: Évaluer le gain de performance réel sur le scénario à 4 corps où les super-nœuds représentaient 34,35 s.
- **Objective**: Étendre `scratch/run_cost_profiler.py` avec `--supernode-variant [baseline|optimized]`, exécuter le benchmark comparatif et générer la planche de diagnostic.
- **Deliverables**:
  - `scratch/run_cost_profiler.py`: Support de `--supernode-variant`.
  - `cost_breakdown_4body_bm6_triple_optimized_300dpi.png`: Planche haute résolution dans `brain/`.
  - `BenchmarkReport`: Mesure chiffrée de la réduction du temps de super-nœuds.
- **Process**:
  - Adapter `run_cost_profiler.py`.
  - Exécuter le benchmark sous BM6 avec triple optimisation (Apex + Kepler + Supernode).
  - Inspecter la planche et consigner le commit Git sur `main`.

## Acceptance Criteria:
### Step 1:
- `OptimizedSupernodeSolver` est implémenté sans aucune allocation NumPy dynamique dans la boucle de Leapfrog.
### Step 2:
- La délégation par ContextVar et paramètre `supernode_evaluator` est 100 % rétro-compatible (baseline par défaut).
### Step 3:
- 100 % des tests unitaires de `test_optimized_supernode.py` passent avec une concordance $< 10^{-12}$.
### Step 4:
- La suite complète de tests (175+ tests) passe sans régression.
- Le temps CPU des super-nœuds sur le 4-corps est significativement réduit.
- Commit consigné sur `main`.

## Observations:
> [!IMPORTANT]
> Cette troisième brique d'optimisation complète le triptyque conormal :
> 1. Kick d'Apex (Itération 170) : pré-calcul des couplages et indices non-feuilles.
> 2. Solveur Kepler-Sundman (Itération 171 & 173) : inlining de Halley et micro-noyau Stumpff $c_0..c_3$.
> 3. Solveur de Super-nœuds (Itération 174) : micro-noyau de Leapfrog sans allocation et factorisation géométrique.
> Tous les modules restent des variantes découplées optionnelles n'altérant en rien la physique pure unsoftened ($\varepsilon = 0, d_0 = 0$).
