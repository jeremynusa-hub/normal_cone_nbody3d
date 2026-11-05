# Implementation Plan of Iteration 175 : Micro-Noyau JIT / LLVM Compilé pour Super-Nœuds Synchrones (Bloc 5)
Created at: 2026-10-09T13:00:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur demande de se focaliser sur le « bloc 5 comme un optimiseur modulaire supplémentaire sans modifier le code original ».
2. **Implicit Strategic Value**: Éliminer le goulet d'étranglement restant de 24,73 s (62,7 % du temps CPU sur le 4-corps hiérarchique avec BM6) causé par l'exécution de 204 358 micro-pas de Leapfrog dans l'interpréteur Python, en compilant la boucle microscopique interne directement en code machine natif (via un micro-noyau JIT / LLVM avec Numba ou un micro-moteur natif vectorisé), tout en garantissant une architecture 100 % modulaire, optionnelle et rétro-compatible sans toucher au code d'origine.
3. **Alignment Criteria**: Concevoir `normal_cone_nbody_float64/solvers/jit_supernode.py` avec la classe `JitSupernodeSolver`, interfaçable directement via `supernode_evaluator` dans `NBodySystem` ou la `ContextVar` `_CURRENT_SUPERNODE_SOLVER`. Maintenir une concordance machine ($< 10^{-12}$) avec le solveur de référence, zéro régression sur les 178 tests existants, et mesurer l'accélération sur le benchmark 4-corps avec génération d'une planche de diagnostic 300 DPI.

## Guidelines:
- **Method**: Développer une variante optionnelle découplée `JitSupernodeSolver`. Utiliser Numba (compilation LLVM à la volée des fonctions numériques sans allocation) avec un mécanisme de repli transparent (fallback) sur `OptimizedSupernodeSolver` si le compilateur JIT n'est pas présent dans l'environnement.
- **Red Lines**: Préservation inconditionnelle du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$), des 10 invariants du groupe de Lie-Bargmann (quantité de mouvement, moment cinétique, centre de masse), et de la réversibilité temporelle microscopique. Aucune modification destructive des fichiers sources d'origine (`solvers/supernode.py`, etc.).
- **Architecture**: Ingestion directe des tableaux contigus `float64`, compilation JIT de la boucle microscopique conormale complète (élimination de tout aller-retour Python par sous-pas), déroulement scalaire des paires pour $k=3$ et $k=4$, et inlining direct des forces différentielles spectatoriales.

## Plan:
### Step 1: Mise en place de l'environnement JIT et architecture de `JitSupernodeSolver`
- **Context**: L'environnement Windows dispose de Python 3.12 et NumPy 2.2.6. Numba 0.68.0 avec LLVM précompilé est disponible via pip.
- **Objective**: Installer `numba` dans l'environnement Python et concevoir la structure de classe `JitSupernodeSolver` avec repli automatique gracieux.
- **Deliverables**:
  - `normal_cone_nbody_float64/solvers/jit_supernode.py`: Squelette de classe `JitSupernodeSolver`, détection de disponibilité JIT et mécanisme de repli.
  - `normal_cone_nbody_float64/solvers/__init__.py`: Export propre de `JitSupernodeSolver`.
- **Process**:
  - Installer `numba` via `pip install numba`.
  - Vérifier l'importation et la compilation d'une fonction simple `@njit`.
  - Créer `jit_supernode.py` définissant l'interface et le gestionnaire de compilation.

### Step 2: Implémentation des micro-noyaux compilés JIT pour la boucle conormale
- **Context**: Les fonctions `_substep_cluster_k3`, `_substep_cluster_k4` et l'application des forces spectatoriales de `solvers/optimized_supernode.py` sont écrites en Python pur sans allocation.
- **Objective**: Compiler la boucle microscopique complète en code machine natif via `@njit(fastmath=False)` (ou `fastmath=True` configurable) pour exécuter des milliers de sous-pas sans aucun surcoût d'interpréteur Python.
- **Deliverables**:
  - `normal_cone_nbody_float64/solvers/jit_supernode.py`: Micro-noyaux compilés JIT pour $k=2$, $k=3$, $k=4$ et général $k \ge 5$, ainsi que la boucle globale de sous-échantillonnage conormal.
- **Process**:
  - Adapter les micro-noyaux scalaires pour être compatibles avec le sous-ensemble de types Numba (`numba.njit`).
  - Implémenter la boucle de Leapfrog conormale compilée qui itère les sous-pas directement en mémoire C contiguë.
  - Intégrer la résolution de la paire serrée interne via le solveur de Kepler-Sundman compilé ou régularisé.

### Step 3: Tests unitaires d'isomorphisme et de non-régression
- **Context**: Le solveur JIT doit garantir une équivalence numérique rigoureuse avec `SupernodeSolver` et `OptimizedSupernodeSolver`.
- **Objective**: Écrire une suite de tests unitaires validant la concordance des trajectoires, la conservation du centre de masse et du moment linéaire, et la robustesse du mode fallback.
- **Deliverables**:
  - `tests/unit/solvers/test_jit_supernode.py`: Suite complète de tests unitaires.
- **Process**:
  - Tester les clusters binaires ($k=2$), ternaires ($k=3$), quaternaires ($k=4$), et avec paires serrées.
  - Vérifier la conservation de $\sum m_i v_i = 0$ à $< 10^{-15}$ et la concordance avec la baseline à $< 10^{-12}$.
  - Tester l'activation via `NBodySystem(supernode_evaluator=JitSupernodeSolver())` et via `ContextVar`.
  - Exécuter la suite complète des tests du projet (`pytest tests/unit/`).

### Step 4: Intégration CLI et benchmark comparatif sur le 4-Corps Hiérarchique
- **Context**: Évaluer le gain réel sur la simulation à 4 corps où les super-nœuds représentaient encore 24,73 s sous BM6.
- **Objective**: Étendre `scratch/run_cost_profiler.py` avec `--supernode-variant [baseline|optimized|jit]`, exécuter le benchmark 4-corps et générer la planche 300 DPI.
- **Deliverables**:
  - `scratch/run_cost_profiler.py`: Option `--supernode-variant jit`.
  - `cost_breakdown_hierarchical_4body_bm6_jit_optimized_300dpi.png`: Planche comparative 300 DPI.
  - Commit Git sur `main`.
- **Process**:
  - Adapter le script de profilage pour supporter la variante `jit`.
  - Lancer la simulation 4-corps sous Blanes-Moan BM6 avec la triple optimisation JIT (Apex + Kepler + Supernode JIT).
  - Générer et analyser la planche de diagnostic comparant Baseline (49s), Optimisé (39s) et JIT.
  - Consigner le commit Git sur `main`.

## Acceptance Criteria:
### Step 1:
- `numba` est installé et opérationnel sous Python 3.12 (ou repli fonctionnel sans crash).
- `JitSupernodeSolver` s'instancie proprement et expose la même API que `OptimizedSupernodeSolver`.
### Step 2:
- La boucle microscopique s'exécute en code machine natif sans aucune allocation ni passage dans l'interpréteur Python par sous-pas.
### Step 3:
- 100 % des tests de `test_jit_supernode.py` passent avec une concordance $< 10^{-12}$.
- 100 % des tests de la suite globale du projet passent sans régression.
### Step 4:
- Le temps CPU passé dans les super-nœuds est réduit d'un ordre de grandeur (attendue sous les quelques secondes).
- La planche haute résolution 300 DPI est générée dans le dossier d'artefacts.
- Le commit est consigné sur la branche `main`.

## Observations:
> [!IMPORTANT]
> - Le Bloc 5 constitue l'ultime frontière d'accélération CPU sans changer la physique conormale : compiler la boucle microscopique interne via LLVM/JIT élimine les millions d'opcodes d'interpréteur Python.
> - La modularité est absolue : si `numba` n'est pas installé sur une machine tierce, le système bascule automatiquement sur `OptimizedSupernodeSolver` sans générer d'erreur.
> - Le zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$) et les invariants symplectiques de Lie-Bargmann sont préservés à la précision machine.
