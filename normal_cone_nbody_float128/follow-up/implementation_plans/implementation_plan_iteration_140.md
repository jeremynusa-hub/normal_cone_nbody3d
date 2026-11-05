# Implementation Plan of Iteration 140 : Aération du Titre et Ajustement des Marges Supérieures sur les Planches Scientifiques (300 DPI)
Created at: 2026-10-08T09:03:00+02:00

## Context:
1. **Explicit Request**: 
   - L'utilisateur signale : "Dans toutes les planche génèrés le titre est collé sur lzs première image du haut, faut aérer".
2. **Implicit Strategic Value**:
   - Rehausser la qualité esthétique et éditoriale des planches scientifiques destinées à la publication (300 DPI).
   - Éliminer tout télescopage ou collision visuelle entre le grand titre principal (`suptitle`) et les en-têtes des panneaux de la première rangée (`axes[0, 0]` et `axes[0, 1]`).
3. **Alignment Criteria**:
   - Configuration d'une marge supérieure explicite (`layout engine rect=[0, 0, 1, 0.96]`) dans `render_10panel_dashboard` et `render_multi_integrator_invariants_dashboard` au sein de [`normal_cone_nbody_float64/diagnostics/visualizer.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/diagnostics/visualizer.py).
   - Positionnement calibré du titre à `y=0.982 - 0.985` garantissant un espacement aéré de plus de 30 points typographiques au-dessus des sous-graphiques.
   - Régénération des 4 planches 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`, `dashboard_multi_integrator_invariants_300dpi.png`) et synchronisation dans les artefacts du brain.
   - 100% de passage de la suite de tests unitaires (132/132 passés).

## Guidelines:
- **Method**: Utiliser l'API officielle de Matplotlib pour `constrained_layout` (`fig.get_layout_engine().set(rect=[0, 0, 1, 0.96])`) pour réserver un bandeau supérieur dédié au titre sans déformer les proportions des subplots.
- **Red Lines**: 
  - Ne pas altérer les données scientifiques, les échelles, ni les calculs physiques sous-jacents.
  - Préserver le mode headless Agg.
- **Architecture**:
  - `diagnostics/visualizer.py` : ajustement du moteur de mise en page et du positionnement du suptitle.
  - `scratch/run_showcase_simulations.py` : réexécution pour regénérer les planches 300 DPI.

## Plan:
### Step 1: Ajustement de l'espacement et des marges dans visualizer.py
- **Context**: Le titre principal était ancré à `y=0.995` avec un `constrained_layout` occupant tout l'espace vertical, provoquant une collision avec la première rangée.
- **Objective**: Ajuster `render_10panel_dashboard` et `render_multi_integrator_invariants_dashboard` avec `layout='constrained'`, `fig.get_layout_engine().set(rect=[0, 0, 1, 0.96])` et `y=0.984`.
- **Deliverables**:
  - `VisualizerLayoutPads`: Code mis à jour dans `diagnostics/visualizer.py`.
- **Process**: Édition via `replace_file_content`.

### Step 2: Régénération des Planches de Publication 300 DPI Aérées
- **Context**: Les 4 planches doivent être recalculées et rendues avec la nouvelle mise en page aérée.
- **Objective**: Exécuter `scratch/run_showcase_simulations.py` et copier les 4 planches PNG mises à jour vers les artefacts du brain.
- **Deliverables**:
  - `AeratedPublicationDashboards`:
    - `dashboard_figure8_bm6_300dpi.png`
    - `dashboard_burrau_bm6_300dpi.png`
    - `dashboard_4body_bm6_300dpi.png`
    - `dashboard_multi_integrator_invariants_300dpi.png`
- **Process**: Exécution du script showcase et copie de fichiers.

### Step 3: Tests Unitaires et Validation Visuelle
- **Context**: S'assurer de la conformité de l'export graphique et de la non-régression.
- **Objective**: Lancer `python -m pytest tests/unit/` et vérifier la validité binaire et les dimensions des fichiers générés.
- **Deliverables**:
  - `LayoutQualityValidation`: Suite de tests validée à 100%.
- **Process**: Exécution de pytest.

## Acceptance Criteria:
### Step 1:
- Les deux méthodes de rendu utilisent `fig.get_layout_engine().set(rect=...)` sans avertissement Matplotlib.
### Step 2:
- Les 4 planches PNG sont générées avec succès avec un espacement net et dégagé entre le grand titre et les graphiques.
### Step 3:
- 132/132 tests unitaires réussis sous `pytest tests/unit/`.

## Observations:
> [!NOTE]
> Cette calibration de mise en page résout un défaut classique d'interaction entre `constrained_layout` et `fig.suptitle` sous Matplotlib en réservant une marge géométrique supérieure explicite pour le titre global.
