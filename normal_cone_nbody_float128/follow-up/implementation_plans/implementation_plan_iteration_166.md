# Implementation Plan of Iteration 166 : Échelles Bilogarithmiques du Graphe Pareto et Généralisation Multi-Modèles des Planches Comparatives
Created at: 2026-10-09T09:43:00+02:00

## Context:
1. **Explicit Request**: Configurer le graphe de Pareto (3.A) en échelles logarithmiques sur les deux axes ($X$ et $Y$), et décliner la planche de comparaison multi-intégrateurs pour chacun des modèles/scénarios de la suite de simulation.
2. **Implicit Strategic Value**: Fournir un diagnostic d'efficacité travail-précision rigoureux et conforme à l'état de l'art géométrique (Hairer, Lubich, Wanner / Blanes & Moan) où les 9 ordres de grandeur d'erreur relative sont résolus sans écrasement sur l'axe $Y$, et permettre aux pairs relecteurs d'évaluer le comportement comparé des 6 schémas symplectiques sur l'ensemble de la phénoménologie céleste (orbite périodique, rencontre chaotique à 3 corps, singularité de collision radiale avec rebond, et 4 corps hiérarchique multi-échelles avec transitions d'atlas).
3. **Alignment Criteria**: Le panneau (3.A) doit être inconditionnellement configuré en échelles bilogarithmiques $\log_{10}$--$\log_{10}$ avec grille adaptée ; 4 planches comparatives haute résolution (300 DPI) distinctes et équilibrées doivent être générées pour chaque modèle avec le profil KPI adapté dans le tableau 3.B ; 100% des tests unitaires doivent être validés.

## Guidelines:
- **Method**: Édition ciblée de `visualizer.py` via `replace_file_content`, extension du script de simulation comparatif `scratch/run_showcase_simulations.py`, exécution des 4 simulations multi-intégrateurs et synchronisation vers le dossier des artefacts.
- **Red Lines**: Préservation inconditionnelle du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$). Interdiction absolue d'écraser des fichiers sources existants avec `overwrite: true`.
- **Architecture**: `normal_cone_nbody_float64/diagnostics/visualizer.py` pour le moteur de visualisation, `scratch/run_showcase_simulations.py` pour l'orchestration des modèles.

## Plan:

### Step 1: Forçage bilogarithmique inconditionnel du panneau Pareto (3.A) et paramétrage du profil KPI dans `visualizer.py`
- **Context**: Le panneau (3.A) omet `ax.set_xscale("log")` et `ax.set_yscale("log")` lorsque `pareto_data` n'est pas fourni (cas des points uniques par intégrateur), écrasant tous les schémas d'ordre supérieur sur $y=0$. De plus, le profil KPI du tableau 3.B est figé à `"smooth_periodic"`.
- **Objective**:
  - Imposer `ax.set_xscale("log")` et `ax.set_yscale("log")` de manière inconditionnelle sur le panneau (3.A), avec formatage propre des bornes et grille bilogarithmique complète.
  - Ajouter un paramètre `profile: str = "smooth_periodic"` à `render_multi_integrator_invariants_dashboard` pour évaluer les KPI du tableau 3.B selon la nature physique du problème (périodique, chaotique, collision, hiérarchique).
- **Deliverables**:
  - `VisualizerCode`: [`visualizer.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/diagnostics/visualizer.py) enrichi.
- **Process**: Édition ciblée avec `replace_file_content`.

### Step 2: Généralisation de l'orchestration multi-intégrateurs pour chaque modèle céleste
- **Context**: `run_showcase_simulations.py` simule déjà les 6 intégrateurs sur Burrau, Radial Collision et 4-Corps, mais n'archive que des sous-échantillonnages d'énergie sans conserver `all_histories` ni appeler `render_multi_integrator_invariants_dashboard`.
- **Objective**:
  - Adapter `run_burrau_pythagorean`, `run_radial_collision` et `run_hierarchical_4body_nested_supernodes` pour conserver les historiques complets `all_histories`.
  - Appeler `render_multi_integrator_invariants_dashboard` pour chaque scénario en passant le profil KPI adéquat (`chaotic_encounter` pour Burrau, `radial_collision` pour la collision, `hierarchical_multiscale` pour le 4-corps).
  - Définir des chemins de sortie distincts :
    - Figure-8 : `scratch/dashboard_multi_integrator_invariants_figure8_300dpi.png` (et `dashboard_multi_integrator_invariants_300dpi.png`)
    - Burrau : `scratch/dashboard_multi_integrator_invariants_burrau_300dpi.png`
    - Collision radiale : `scratch/dashboard_multi_integrator_invariants_radial_300dpi.png`
    - 4-corps hiérarchique : `scratch/dashboard_multi_integrator_invariants_4body_300dpi.png`
- **Deliverables**:
  - `ShowcaseScript`: `scratch/run_showcase_simulations.py` mis à jour.
- **Process**: Édition ciblée via `replace_file_content`.

### Step 3: Exécution des simulations multi-modèles et synchronisation des planches 300 DPI
- **Context**: Les 4 planches haute résolution doivent être produites et copiées dans l'artifact directory pour mise à disposition de l'utilisateur.
- **Objective**:
  - Exécuter la génération complète des 4 planches comparatives.
  - Copier les 4 fichiers PNG dans `C:\Users\jerem\.gemini\antigravity\brain\59313b1a-0dde-45d6-a476-1bcb513cf82e\`.
- **Deliverables**:
  - `ArtifactBoards`: 4 planches comparatives 300 DPI synchronisées.
- **Process**: Exécution via `run_command` et copie des fichiers.

### Step 4: Validation de non-régression et commit git
- **Context**: Confirmer que l'ensemble des 151 tests unitaires passe et versionner les modifications.
- **Objective**:
  - Exécuter `pytest tests/unit/`.
  - Enregistrer un commit git propre sur la branche `main`.
- **Deliverables**:
  - `GitCommit`: Commit git sur `visualizer.py` et `run_showcase_simulations.py`.
- **Process**: `pytest tests/unit/`, `git add`, `git commit`.

## Acceptance Criteria:
### Step 1:
- Le panneau (3.A) applique inconditionnellement `ax.set_xscale("log")` et `ax.set_yscale("log")`.
- La méthode `render_multi_integrator_invariants_dashboard` accepte un argument `profile`.
### Step 2:
- Les 4 scénarios génèrent chacun leur propre planche comparative multi-intégrateurs.
### Step 3:
- Les 4 planches PNG à 300 DPI existent dans `scratch/` et dans le dossier des artefacts, avec une taille $> 1.5$ Mo chacune.
- Le nuage de points du panneau (3.A) de chaque planche est lisible et distribué sur plusieurs ordres de grandeur verticaux sans superposition à $y=0$.
### Step 4:
- 151/151 tests unitaires validés avec succès en $< 10$s.
- Commit git enregistré sur `main`.

## Observations:
> [!IMPORTANT]
> Pour des scénarios chaotiques ou multi-échelles (Burrau, 4-corps), les schémas d'ordre inférieur (Strang, Yoshida) subissent des erreurs de troncature accumulées plus importantes que sur l'orbite périodique de Chenciner. Le passage en échelle bilogarithmique permettra précisément d'apprécier la séparation nette de performance entre schémas de bas ordre et compositions d'ordre supérieur (BM6/PBM6).
