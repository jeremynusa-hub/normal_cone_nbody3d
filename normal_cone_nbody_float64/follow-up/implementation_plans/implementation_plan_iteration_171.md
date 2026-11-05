# Implementation Plan of Iteration 171 : Module Dérivé du Solveur de Feuilles Képlériennes et Dérive Conormale Optimisée
Created at: 2026-10-09T12:18:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur commande de créer le module dérivé pour optimiser le solveur de feuilles képlériennes et la dérive conormale, à la suite des gisements identifiés lors de l'exploration de l'itération précédente (« fait ce module dérivé »).
2. **Implicit Strategic Value**: Assurer un gain de performance substantiel sur le second goulot d'étranglement de la dynamique conormale unsoftened (qui représente 38 % du temps CPU restant) sans jamais altérer ni régresser les fichiers sources canoniques de référence (`kepler_sundman.py`, `stumpff.py`, `common.py`), tout en garantissant une précision bit-exacte à la précision machine IEEE 754 float64.
3. **Alignment Criteria**: Le nouveau module dérivé doit être 100 % découplé, interchangeable de façon modulaire dans `NBodySystem`, afficher un gain d'environ 2x sur l'étage de résolution képlérienne, préserver les invariants de Lie-Bargmann et la conservation symplectique au niveau du bruit d'arrondi machine ($< 10^{-14}$), et réussir 100 % des tests unitaires existants et nouveaux.

## Guidelines:
- **Method**: Découplage strict via une architecture modulaire variante (comme réalisé pour `OptimizedScreenedApexPotential`), sans altération du code cœur canonique.
- **Red Lines**: Interdiction formelle d'écraser ou de modifier `solvers/kepler_sundman.py`, `solvers/stumpff.py`, ou d'introduire le moindre softening ($\varepsilon = 0, d_0 = 0$).
- **Architecture**:
  - Implémentation du module dérivé `normal_cone_nbody_float64/solvers/optimized_kepler_sundman.py`.
  - Intégration optionnelle non intrusive dans `NBodySystem` via un paramètre optionnel `kepler_evaluator`.
  - Intégration des tests unitaires exhaustifs dans `tests/unit/solvers/test_optimized_kepler_sundman.py`.

## Plan:
### Step 1: Conception et implémentation du module dérivé `solvers/optimized_kepler_sundman.py`
- **Context**: Le solveur canonique `step_kepler_sundman` recalcule inutilement $c_4, c_5$, réévalue Stumpff à la fin du pas et alloue des closures Python.
- **Objective**: Implémenter un solveur de feuilles dérivé avec micro-noyau Stumpff ciblé ($c_0..c_3$), inlining de l'itération de Halley cubique et réutilisation immédiate de la dérivée conormale convergée $f' = r(\tau)$.
- **Deliverables**:
  - `OptimizedKeplerSundmanSolver`: Classe et fonctions dans `normal_cone_nbody_float64/solvers/optimized_kepler_sundman.py`.
- **Process**:
  - Définir `stumpff_c0123_scalar` optimisé et bit-exact.
  - Définir `step_kepler_sundman_fast` avec inlining Halley et réutilisation de $c_0..c_3$ et $\rho_{\text{final}} = f_1$.
  - Intégrer la pré-compilation des projections linéaires fixes Cartésien $\leftrightarrow$ Jacobi par topologie d'arbre pour accélérer le drift conormal.
  - Exporter le module dans `normal_cone_nbody_float64/solvers/__init__.py`.

### Step 2: Intégration modulaire rétrocompatible dans `NBodySystem`
- **Context**: `NBodySystem` orchestre les étapes d'intégration via `step_drift_and_kepler`.
- **Objective**: Permettre l'injection optionnelle du solveur képlérien optimisé via `kepler_evaluator: Optional[Any] = None` dans `NBodySystem.__init__`.
- **Deliverables**:
  - `NBodySystem`: Câblage modulaire rétrocompatible dans `normal_cone_nbody_float64/system/nbody_system.py`.
- **Process**:
  - Ajouter `kepler_evaluator` dans `NBodySystem.__init__`.
  - Si fourni, déléguer le drift et la résolution de feuilles au solveur optimisé; sinon, conserver le comportement canonique actuel.

### Step 3: Tests unitaires et certification numérique bit-exacte
- **Context**: Les invariants physiques (énergie, moment cinétique, déterminant wronskien de Battin) doivent être vérifiés sur des régimes orbitaux variés.
- **Objective**: Valider que le solveur dérivé produit des résultats indiscernables du solveur de référence ($< 10^{-14}$) sur des cas elliptiques, hyperboliques, paraboliques et de quasi-collision.
- **Deliverables**:
  - `TestOptimizedKeplerSundman`: Fichier de test unitaire `tests/unit/solvers/test_optimized_kepler_sundman.py`.
- **Process**:
  - Écrire des tests comparatifs directs face à `step_kepler_sundman`.
  - Tester l'unimodularité du wronskien $|f \dot{g} - \dot{f} g - 1| < 10^{-15}$.
  - Vérifier la réversibilité temporelle microscopique $S(-dt) \circ S(dt) = \text{Id}$.

### Step 4: Benchmark comparatif de coût et génération de planche de diagnostic
- **Context**: L'outil `scratch/run_cost_profiler.py` permet d'instrumenter et de visualiser les gains par étape.
- **Objective**: Mesurer l'impact cumulé de l'optimisation (Apex + Kepler) sur le problème de Burrau avec BM6 et générer la planche de diagnostic haute résolution.
- **Deliverables**:
  - `cost_breakdown_burrau_bm6_fully_optimized_300dpi.png`: Planche de diagnostic 300 DPI dans `brain/`.
- **Process**:
  - Étendre `scratch/run_cost_profiler.py` pour supporter la variante combinée.
  - Exécuter la simulation et comparer les métriques (temps Kepler, débit global, part CPU).

### Step 5: Validation globale de la suite de tests et enregistrement Git
- **Context**: Toute modification doit maintenir l'intégralité des tests unitaires au vert.
- **Objective**: Valider 100 % de la suite de tests unitaires et enregistrer le commit sur `main`.
- **Deliverables**:
  - `GitCommit`: Commit sur la branche `main`.
- **Process**:
  - Exécuter `pytest tests/unit/`.
  - Enregistrer le commit Git avec un message conventionnel.

## Acceptance Criteria:
### Step 1:
- Le module `optimized_kepler_sundman.py` s'importe sans erreur et implémente la résolution sans modification des fichiers existants.
### Step 2:
- `NBodySystem` accepte `kepler_evaluator` sans casser la compatibilité ascendante (100 % des tests existants passent à l'identique).
### Step 3:
- L'écart maximal entre les positions et impulsions résolues par la variante optimisée et la baseline est $< 10^{-14}$.
- La conservation du moment cinétique et l'unimodularité de Battin sont validées à $< 10^{-15}$.
### Step 4:
- Le temps CPU mesuré sur l'étage Kepler présente un facteur d'accélération mesuré $\ge 1,7\mathrm{x}$.
- La planche 300 DPI est générée et consultable dans `brain/`.
### Step 5:
- 100 % des tests unitaires sont validés avec succès.
- Le commit Git est consigné sur `main`.

## Observations:
> [!IMPORTANT]
> Aucune formule d'adoucissement ($\varepsilon > 0$ ou $d_0 > 0$) n'est introduite. La physique conormale unsoftened est strictement préservée.
> Le code canonique `kepler_sundman.py` et `stumpff.py` demeure rigoureusement intact en tant que ground truth académique.
