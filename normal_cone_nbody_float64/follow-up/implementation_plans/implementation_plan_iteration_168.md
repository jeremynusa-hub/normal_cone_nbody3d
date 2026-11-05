# Implementation Plan of Iteration 168 : Scellement Architectural Inconditionnel des Super-Nœuds Conormaux et Régénération Certifiée du Benchmark 4-Corps
Created at: 2026-10-09T10:57:30+02:00

## Context:
1. **Explicit Request**: L'utilisateur relève que la désactivation des super-nœuds n'a « aucun sens » et demande pourquoi cette option est pilotable au lieu d'être formellement interdite. Il valide l'intervention pour sceller cette option et relancer le 4-corps avec l'architecture complète.
2. **Implicit Strategic Value**: Assurer l'étanchéité et l'intégrité de la bibliothèque logicielle face au formalisme théorique du Chapitre 4. La régularisation par super-nœuds sous temps conormal déplié $\tau_S$ est le partenaire mathématique indissociable de l'arbre de Jacobi pour toute collision ou approche multiple ($k \ge 3$ corps). Pouvoir désactiver cette régularisation au pilotage constituait une régression architecturale héritée d'anciennes études d'ablation expérimentales.
3. **Alignment Criteria**:
   - `NBodySystem` active systématiquement les super-nœuds par défaut.
   - Toute désactivation éventuelle pour étude d'ablation académique est confinée à un paramètre privé explicite (`_ablation_bypass_supernodes=True`) avec émission d'un warning formel.
   - Le scénario à 4 corps s'exécute avec les super-nœuds conormaux actifs pour tous les schémas.
   - Élimination du pic d'instabilité de Suzuki à $\Delta t = 0.010$ sur BM6 et certification cohérente du tableau 3.B.

## Guidelines:
- **Method**: Verrouiller le comportement par défaut dans `NBodySystem` et nettoyer les scripts d'exposition pour refléter fidèlement le modèle conormal sans adoucissement.
- **Red Lines**: Préservation inconditionnelle du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$). Interdiction absolue d'utiliser `write_to_file` avec `overwrite: true` sur le code source de la bibliothèque.
- **Architecture**: `nbody_system.py` est la façade unifiée, `solvers/supernode.py` est le solveur canonique conormal sous temps déplié $\tau_S$.

## Plan:
### Step 1: Scellement architectural et relégation de l'ablation des super-nœuds dans `nbody_system.py`
1. **Contexte initial**: `NBodySystem.__init__` expose `auto_supernodes: bool = True`, permettant à n'importe quel code appelant de désactiver la régularisation de singularités triples en passant `auto_supernodes=False`.
2. **Objectif**: Rendre la régularisation par super-nœuds intrinsèque et obligatoire, tout en réservant l'ablation à un paramètre privé explicite `_ablation_bypass_supernodes` émettant un avertissement d'intégrité.
3. **Livrables**: `normal_cone_nbody_float64/system/nbody_system.py`.
4. **Procédé**:
   - Mettre à jour `__init__` : `auto_supernodes: Optional[bool] = None`, `_ablation_bypass_supernodes: bool = False`.
   - Si `auto_supernodes is False`, émettre un `UserWarning` indiquant que `auto_supernodes` est déprécié et rediriger vers `_ablation_bypass_supernodes`.
   - Par défaut, `self.auto_supernodes = not _ablation_bypass_supernodes` (donc toujours `True` en usage normal).

### Step 2: Harmonisation et correction de l'orchestration du benchmark 4-corps dans `run_showcase_simulations.py`
1. **Contexte initial**: Dans `scratch/run_showcase_simulations.py`, `run_hierarchical_4body_nested_supernodes()` appelait `NBodySystem` avec `auto_supernodes=False`, sabotant la régularisation du triplet $(1, 2, 3)$.
2. **Objectif**: Activer l'architecture conormale complète (super-nœuds synchrones sous temps déplié $\tau_S$ + ordonnancement dyadique) pour BM6 et l'ensemble des intégrateurs du benchmark 4-corps.
3. **Livrables**: `scratch/run_showcase_simulations.py`.
4. **Procédé**:
   - Supprimer `auto_supernodes=False` de tous les appels `NBodySystem` dans la fonction 4-corps.
   - Vérifier la robustesse de l'exécution sur la grille de pas de test Pareto.

### Step 3: Exécution du benchmark 4-corps et régénération des planches 300 DPI
1. **Contexte initial**: La planche actuelle du 4-corps présente un pic d'erreur à $\Delta t = 0.010$ et tous les statuts à FAIL.
2. **Objectif**: Régénérer la planche comparative 300 DPI et l'application interactive avec des résultats physiques régularisés et certifiés.
3. **Livrables**:
   - `scratch/dashboard_multi_integrator_invariants_4body_300dpi.png`
   - `scratch/dashboard_4body_bm6_300dpi.png`
   - Leurs copies dans `brain/`
   - `interactive_conormal_showcase.html`
4. **Procédé**: Exécution de `run_showcase_simulations.py` et `build_showcase_html.py`.

### Step 4: Validation de non-régression et enregistrement git
1. **Contexte initial**: La suite de tests unitaires doit continuer de passer à 100%.
2. **Objectif**: Exécuter les 151 tests et consigner les modifications par un commit git clair.
3. **Livrables**: Rapport de tests unitaires et commit git sur `main`.

## Acceptance Criteria:
- **Critère 168.1**: `NBodySystem()` instancie systématiquement un système avec super-nœuds actifs par défaut (`self.auto_supernodes == True`).
- **Critère 168.2**: Le paramètre privé `_ablation_bypass_supernodes=True` est requis pour forcer une étude d'ablation sans super-nœuds et émet un avertissement.
- **Critère 168.3**: Sur le graphe de Pareto du 4-corps, le pic aberrant à $\Delta t = 0.010$ de BM6/PBM6 est totalement éradiqué grâce à la régularisation conormale du triplet.
- **Critère 168.4**: 151/151 tests unitaires validés avec succès.

## Observations:
> [!IMPORTANT]
> En scellant architecturalement les super-nœuds, nous éliminons la divergence artificielle entre le formalisme mathématique du Chapitre 4 (qui exige la régularisation simultanée des clusters compacts) et son implémentation logicielle.
