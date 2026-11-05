# Implementation Plan of Iteration 139 : Conception et Génération de la Planche Dédiée de Conservation des Invariants Multi-Intégrateurs (300 DPI)
Created at: 2026-10-08T08:56:00+02:00

## Context:
1. **Explicit Request**: 
   - L'utilisateur demande d'ajouter un graphe prouvant la conservation des invariants sur les planches, interroge sur la permanence des invariants physiques ("sont ils toujours les même ?"), et choisit l'Option A (planche comparative dédiée multi-intégrateurs).
2. **Implicit Strategic Value**:
   - Fournir une preuve graphique irréfutable de la conservation des 10 invariants de Lie-Bargmann ($H, \mathbf{P}, \mathbf{L}, \mathbf{G}$) pour l'ensemble des intégrateurs symplectiques de la bibliothèque (Strang S2, Yoshida S4, BM4, PBM4, BM6, PBM6).
   - Mettre en valeur visuellement la séparation spectrale : prouver que tous les schémas conservent l'impulsion $\mathbf{P}$, le moment cinétique $\mathbf{L}$ et le boost $\mathbf{G}$ au plancher machine IEEE 754 ($\le 10^{-15}$), tandis que l'énergie $H$ discrimine précisément l'ordre effectif de chaque méthode (de $10^{-4}$ à $10^{-13}$).
3. **Alignment Criteria**:
   - Ajout propre et modulaire de la méthode `render_multi_integrator_invariants_dashboard` dans [`normal_cone_nbody_float64/diagnostics/visualizer.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/diagnostics/visualizer.py).
   - Intégration dans `scratch/run_showcase_simulations.py` pour produire [`dashboard_multi_integrator_invariants_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_multi_integrator_invariants_300dpi.png) à 300 DPI.
   - 100% de tests unitaires passés sous `pytest tests/unit/`.

## Guidelines:
- **Method**: Respect strict de la charte graphique et de la palette chromatique établie (Strang=Ambre, Yoshida=Bleu, BM4=Violet, PBM4=Pourpre, BM6=Émeraude, PBM6=Cyan).
- **Red Lines**: 
  - Ne pas dégrader la méthode existante `render_10panel_dashboard`.
  - Pas d'adoucissement ($\varepsilon = 0, d_0 = 0$).
  - Respect du mode Agg headless (`matplotlib.use('Agg')`).
- **Architecture**:
  - `diagnostics/visualizer.py` : implémentation de `render_multi_integrator_invariants_dashboard`.
  - `scratch/run_showcase_simulations.py` : appel et sauvegarde de la planche multi-intégrateurs.
  - `tests/unit/diagnostics/test_visualizer.py` : test unitaire dédié.

## Plan:
### Step 1: Implémentation de render_multi_integrator_invariants_dashboard dans visualizer.py
- **Context**: `PublicationVisualizer` dispose actuellement d'une méthode pour un intégrateur unique. Une planche dédiée comparative à 6 panneaux est requise.
- **Objective**: Implémenter `render_multi_integrator_invariants_dashboard(self, histories: Dict[str, HistoryTracker], title: str, pareto_data: Optional[dict] = None, save_path: Optional[str] = None)` comprenant :
  1. Panneau (1.A) : Dérive relative d'énergie $\log_{10}(\Delta H / |H_0|)(t)$.
  2. Panneau (1.B) : Dérive du moment cinétique total $\log_{10}(\|\Delta \mathbf{L}\|)(t)$.
  3. Panneau (2.A) : Dérive de l'impulsion linéaire totale $\log_{10}(\|\Delta \mathbf{P}\|)(t)$.
  4. Panneau (2.B) : Dérive du centre de masse $\log_{10}(\|\Delta \mathbf{G}\|)(t)$.
  5. Panneau (3.A) : Frontière de Pareto Work-Precision (log-log).
  6. Panneau (3.B) : Tableau comparatif scientifique des KPIs pour les 6 intégrateurs.
- **Deliverables**:
  - `MultiIntegratorVisualizerMethod`: Nouvelle méthode dans `diagnostics/visualizer.py`.
- **Process**: Édition chirurgicale via `replace_file_content`.

### Step 2: Intégration et Génération de la Planche dans le Script Showcase
- **Context**: Les simulations comparatives des 6 schémas sont déjà exécutées sur Figure-8 et Burrau dans `scratch/run_showcase_simulations.py`.
- **Objective**: Collecter les `HistoryTracker` de chaque intégrateur dans `run_figure8_comparison()`, appeler `render_multi_integrator_invariants_dashboard`, enregistrer la planche à 300 DPI (`scratch/dashboard_multi_integrator_invariants_300dpi.png`) et la copier dans les artefacts du brain.
- **Deliverables**:
  - `MultiIntegratorDashboardArtifact`: `dashboard_multi_integrator_invariants_300dpi.png` à 300 DPI dans le brain.
- **Process**: Modification de `scratch/run_showcase_simulations.py` et exécution.

### Step 3: Tests Unitaires et Validation Visuelle
- **Context**: Assurer la couverture et la robustesse de la nouvelle méthode de visualisation.
- **Objective**: Enrichir `tests/unit/diagnostics/test_visualizer.py` pour tester `render_multi_integrator_invariants_dashboard` avec de faux historiques ou des simulations légères, et vérifier le passage à 100% de `pytest tests/unit/`.
- **Deliverables**:
  - `VisualizerUnitTestSuite`: Tests unitaires enrichis et validés.
- **Process**: Lancement de pytest et inspection du fichier PNG généré.

## Acceptance Criteria:
### Step 1:
- `render_multi_integrator_invariants_dashboard` accepte un dictionnaire quelconque d'intégrateurs et gère les échelles semilogarithmiques sans exception mathématique (`np.maximum(..., 1e-18)`).
### Step 2:
- La planche `dashboard_multi_integrator_invariants_300dpi.png` est générée à 300 DPI avec une taille valide (> 1.5 Mo) et des courbes nettes pour les 6 intégrateurs.
### Step 3:
- 100% des tests unitaires réussis sous `pytest tests/unit/`.

## Observations:
> [!IMPORTANT]
> Cette planche prouve empiriquement la préservation géométrique exacte des symétries de translation spatiale ($\mathbf{P}$), rotation ($\mathbf{L}$) et boost galiléen ($\mathbf{G}$) au plancher machine pour tous les intégrateurs, tout en exhibant le gain de 10 ordres de grandeur sur l'énergie $H$ apporté par BM6/PBM6 face à Strang.
