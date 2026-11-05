# Implementation Plan of Iteration 170 : Variante Optimisée du Kick d'Apex et Benchmark Comparatif Multi-Planches (2 Planches 300 DPI)
Created at: 2026-10-09T12:09:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur demande de créer l'optimisation identifiée sur le kick d'apex *sans modifier le code cœur actuel* (proposée comme une variante découplée), puis de benchmarker les différences via le système instrumenté pour produire deux planches de diagnostic complètes (Baseline vs Variante Optimisée).
2. **Implicit Strategic Value**: Valider empiriquement l'hypothèse d'accélération formulée lors de l'Itération 169 (où le kick d'apex consommait 66,2% du temps CPU total). En éliminant les allocations de tenseurs 3D `(N, N, 3)`, en mettant en cache la matrice des couplages de masses $G m_i m_j$ et la liste des paires non-feuilles par topologie d'arbre, on s'attend à diviser par 2 à 3 le temps de l'étage apex sans aucune altération de la physique conormale.
3. **Alignment Criteria**:
   - **Non-intrusion stricte** : Ne pas modifier les fichiers sources existants (`potential/screened_apex.py`, `integrators/common.py`).
   - **Variante modulaire** : Encapsuler l'optimisation dans une classe dédiée `OptimizedScreenedApexPotential` dans un nouveau module `potential/optimized_screened_apex.py`.
   - **Équivalence numérique bit-exacte** : L'écart maximal entre les forces d'apex de référence et optimisées doit être $\le 10^{-15}$ (niveau du bruit d'arrondi IEEE 754).
   - **Livrable visuel** : Générer 2 planches distinctes haute résolution (300 DPI) dans `brain/` : une planche pour la Baseline et une planche pour la Variante Optimisée.

## Guidelines:
- **Architecture**:
  - `normal_cone_nbody_float64/potential/optimized_screened_apex.py`: Nouvelle classe `OptimizedScreenedApexPotential` héritant ou compatible avec l'interface conormale.
  - `NBodySystem`: Accepter un paramètre optionnel `apex_potential: Optional[ScreenedApexPotential] = None` (fallback transparent sur le solveur standard par défaut pour garantir 100% de rétrocompatibilité).
  - `scratch/run_cost_profiler.py`: Support de l'argument `--apex-variant [baseline|optimized]`.
- **Red Lines**: Zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$) inconditionnel. Conservation stricte des 10 invariants du groupe de Lie-Bargmann.

## Plan:

### Step 1: Conception et implémentation de `OptimizedScreenedApexPotential` dans `potential/optimized_screened_apex.py`
- **Context**: `ScreenedApexPotential.compute_apex_forces` alloue des tenseurs 3D `(N, N, 3)` et recalcule la matrice des masses et la liste des paires de feuilles à chaque appel (74 124 fois).
- **Objective**: Développer une variante vectorisée à pré-allocation de mémoire et mise en cache topologique.
- **Deliverables**:
  - `OptimizedScreenedApexPotential`: Classe optimisée implémentant :
    1. Mise en cache de la matrice des couplages de masse $G m_i m_j$.
    2. Cache des indices de paires non-feuilles indexé sur l'identifiant de topologie de l'arbre (`tree.root.id` ou hash de structure), invalidé uniquement lors des 8 transitions d'arbre.
    3. Évaluation vectorisée 2D directe sur les paires non-feuilles sans instanciation de tenseur 3D.
    4. Méthode `step_apex_kick` optimisée réutilisant l'état ou opérant en mémoire tampon contiguë.
- **Process**: Créer le nouveau fichier `potential/optimized_screened_apex.py` et exporter la classe dans `potential/__init__.py`.

### Step 2: Câblage modulaire de la variante dans `NBodySystem` et tests de conformité bit-exacte
- **Context**: `NBodySystem` instancie actuellement en dur `ScreenedApexPotential`.
- **Objective**: Permettre l'injection de dépendance `apex_evaluator: Optional[Any] = None` dans `NBodySystem.__init__`.
- **Deliverables**:
  - `NBodySystem`: Paramètre optionnel `apex_evaluator: Optional[Any] = None` (si `None`, utilise `ScreenedApexPotential` d'origine).
  - `tests/unit/potential/test_optimized_screened_apex.py`: Test unitaire validant l'isomorphisme mathématique bit-exact ($\|F_{\text{opt}} - F_{\text{base}}\| < 10^{-15}$) sur les 4 cas canoniques.
- **Process**: Modifier `NBodySystem` de manière non-intrusive et ajouter le test de validation.

### Step 3: Adaptation de l'outil CLI `scratch/run_cost_profiler.py` pour la comparaison multi-planches
- **Context**: L'utilisateur souhaite comparer baseline et variante optimisée et obtenir deux planches.
- **Objective**: Étendre le script CLI pour supporter `--apex-variant [baseline|optimized]` et permettre l'exécution enchaînée des deux variantes.
- **Deliverables**:
  - `scratch/run_cost_profiler.py`: Support de `--apex-variant` avec génération de noms de fichiers explicites :
    * `cost_breakdown_burrau_bm6_baseline_300dpi.png`
    * `cost_breakdown_burrau_bm6_optimized_300dpi.png`
- **Process**: Mettre à jour `run_cost_profiler.py`.

### Step 4: Exécution du benchmark comparatif sur Burrau avec BM6 et génération des 2 planches 300 DPI
- **Context**: Mesurer les gains réels en temps CPU, débit (pas/s), et part relative de l'étage apex.
- **Objective**: Exécuter successivement la Baseline et la Variante Optimisée sur Burrau avec BM6 ($\Delta t = 0{,}01$, $t \in [0, 16{,}0]$).
- **Deliverables**:
  - `cost_breakdown_burrau_bm6_baseline_300dpi.png` copié dans `brain/`.
  - `cost_breakdown_burrau_bm6_optimized_300dpi.png` copié dans `brain/`.
  - Analyse comparative chiffrée (speedup de l'étage apex, speedup global, débit).
- **Process**: Exécuter via `run_command` et vérifier les sorties graphiques.

### Step 5: Validation globale de la suite de tests unitaires et enregistrement Git
- **Context**: Aucune régression ne doit être introduite.
- **Objective**: Valider l'intégralité de la suite de tests unitaires et consigner les commits sur `main`.
- **Deliverables**:
  - 100% des tests unitaires passants (`pytest tests/unit/`).
  - Commit git sur `main`.
- **Process**: Exécuter `pytest` et valider les modifications.

## Acceptance Criteria:
- Le code cœur historique (`potential/screened_apex.py`, `integrators/common.py`) n'est pas altéré dans sa logique de base.
- L'écart maximal de force entre la baseline et la variante optimisée est $< 10^{-15}$.
- Les deux planches 300 DPI sont générées dans `brain/` et visuellement comparables.
- Le temps CPU de l'étage apex est significativement réduit dans la planche optimisée.
- 100% des tests unitaires sont passants.
