# Implementation Plan of Iteration 141 : Agrandissement Majeur de la Marge Supérieure du Titre (> 400px à 300 DPI)
Created at: 2026-10-08T09:13:00+02:00

## Context:
1. **Explicit Request**: 
   - L'utilisateur constate : "Ça change rien", indiquant que l'augmentation de marge précédente (+15 pixels) est restée visuellement imperceptible une fois l'image de 7000 pixels mise à l'échelle sur un écran d'ordinateur standard.
2. **Implicit Strategic Value**:
   - Établir une véritable aération typographique à l'échelle de l'affichage écran (~50 à 70 pixels d'espace blanc franc lors de la visualisation plein écran), ce qui requiert un espacement absolu de **350 à 450 pixels** sur les rendus haute résolution 300 DPI.
3. **Alignment Criteria**:
   - Configuration du rectangle de mise en page à `rect=[0, 0, 1, 0.90]` pour la planche 10-panneaux (au lieu de 0.96) avec `figsize=(16, 25)` et `y=0.965`.
   - Configuration du rectangle à `rect=[0, 0, 1, 0.89]` pour la planche multi-intégrateurs avec `figsize=(16, 20)` et `y=0.955`.
   - Mesure empirique par analyse de pixels validant un intervalle continu de lignes blanches supérieur à **350 pixels** entre le bas du titre et le haut des sous-graphiques.
   - Régénération des 4 planches 300 DPI et copie dans le brain.

## Guidelines:
- **Method**: Mesure pixel-par-pixel via PIL / NumPy pour garantir objectivement le respect de la séparation supérieure avant validation.
- **Red Lines**: 
  - Préserver les résolutions 300 DPI et l'intégrité de tous les graphes scientifiques.
- **Architecture**:
  - `diagnostics/visualizer.py` : ajustement drastique des paramètres de marge haute.
  - `scratch/run_showcase_simulations.py` : recalcul et réexport des planches.

## Plan:
### Step 1: Ajustement Drastique des Marges Supérieures dans visualizer.py
- **Context**: L'analyse de pixels a révélé que la marge antérieure n'était que de 140 pixels à 300 DPI (soit à peine 18 pixels sur écran).
- **Objective**: Passer `rect` à `[0, 0, 1, 0.90]` sur `render_10panel_dashboard` et `[0, 0, 1, 0.89]` sur `render_multi_integrator_invariants_dashboard` avec augmentation de la hauteur totale de figure (`25` et `20` pouces) pour ménager 350-450 pixels de blanc.
- **Deliverables**:
  - `DeepMarginVisualizer`: Code mis à jour dans `diagnostics/visualizer.py`.
- **Process**: Édition via `replace_file_content`.

### Step 2: Régénération et Vérification Pixel des 4 Planches 300 DPI
- **Context**: Les planches doivent être recalculées et certifiées par inspection numérique de l'intervalle de pixels.
- **Objective**: Exécuter `scratch/run_showcase_simulations.py`, mesurer le gap effectif via PIL (qui doit être $\ge 350$ px), et synchroniser les 4 planches PNG dans le brain.
- **Deliverables**:
  - `VerifiedSpaciousDashboards`: Les 4 fichiers PNG 300 DPI dans le brain avec grand dégagement supérieur.
- **Process**: Exécution du script et copie des fichiers.

### Step 3: Tests Unitaires et Validation Visuelle
- **Context**: Valider la conformité de la suite de tests et l'absence de régression.
- **Objective**: Lancer `python -m pytest tests/unit/` et confirmer 100% de réussite.
- **Deliverables**:
  - `TestSuiteValidation`: 132/132 tests validés.
- **Process**: Exécution de pytest.

## Acceptance Criteria:
### Step 1:
- `rect=[0, 0, 1, 0.90]` et `[0, 0, 1, 0.89]` intégrés dans `visualizer.py`.
### Step 2:
- L'analyse d'image sur les planches exportées démontre formellement au moins **350 pixels consécutifs sans texte ni trait** entre le titre et la première rangée.
### Step 3:
- 132/132 tests unitaires réussis sous `pytest tests/unit/`.

## Observations:
> [!IMPORTANT]
> L'analyse quantitative de pixels a démontré qu'un saut de 140 à plus de 400 pixels d'espace blanc est indispensable pour que l'aération soit visible et aérée à l'œil nu lors de la visualisation globale d'une image haute résolution 300 DPI.
