# Implementation Plan of Iteration 165 : Harmonisation des Seuils et Optimisation du Tableau Comparatif Multi-Intégrateurs
Created at: 2026-10-09T09:36:00+02:00

## Context:
1. **Explicit Request**: Retirer la colonne "Order" du tableau 3.B de la planche comparative multi-intégrateurs (`dashboard_multi_integrator_invariants_300dpi.png`), résoudre les faux échecs `FAIL` sur les schémas d'ordre inférieur (Strang, Yoshida, BM4) dus à un seuil inadapté, homogénéiser les seuils de plancher machine entre les graphiques (à $-14.0$ uniforme) et expliciter la convention de lecture logarithmique *(Lower is better)*.
2. **Implicit Strategic Value**: Fournir une planche comparative de publication scientifique rigoureuse, exempte d'ambiguïtés d'interprétation pour les pairs relecteurs, où chaque schéma est évalué de façon équitable selon son ordre de convergence sans pénalisation artificielle, et où les seuils physiques machine sont cohérents.
3. **Alignment Criteria**: Le tableau 3.B ne doit plus comporter de colonne "Order", tous les schémas symplectiques conformes doivent afficher `PASS`, les 4 panneaux d'invariants doivent présenter un seuil homogène à $-14.0$ avec libellé explicite, et la planche régénérée doit préserver la conformité visuelle et typographique.

## Guidelines:
- **Method**: Modification ciblée de `visualizer.py` via `replace_file_content`, tests unitaires de non-régression, et régénération de la planche 300 DPI.
- **Red Lines**: Ne modifier aucune constante physique ($\varepsilon = 0, d_0 = 0$). Ne pas écraser de fichier avec `overwrite: true`.
- **Architecture**: `visualizer.py` dans `normal_cone_nbody_float64/diagnostics/` et script de rendu `scratch/run_showcase_simulations.py`.

## Plan:

### Step 1: Retrait de la colonne "Order" et calibration des tolérances dans le Tableau 3.B
- **Context**: Le tableau 3.B de `render_multi_integrator_invariants_dashboard` contient une colonne "Order" redondante et évalue le statut avec un seuil de pente séculaire fixe $10^{-11}$ propre à l'ordre 6.
- **Objective**:
  - Supprimer la colonne "Order" des en-têtes et des lignes du tableau.
  - Adapter les tolérances d'énergie et de pente séculaire en fonction de l'ordre théorique de chaque intégrateur pour une évaluation équitable et mathématiquement fondée.
  - Mettre à jour les largeurs relatives des colonnes (`colWidths`).
- **Deliverables**:
  - `VisualizerCode`: `normal_cone_nbody_float64/diagnostics/visualizer.py` mis à jour.
- **Process**: Édition ciblée de `visualizer.py` via `replace_file_content`.

### Step 2: Homogénéisation des seuils de plancher machine et mentions explicites sur les axes
- **Context**: Le panneau (2.A) a une ligne à $-15.0$ tandis que les autres sont à $-14.0$, sans mention explicite que la région inférieure représente une précision supérieure.
- **Objective**:
  - Harmoniser la ligne de référence à $-14.0$ sur les 4 panneaux d'invariants (1.A, 1.B, 2.A, 2.B) avec le libellé unifié `Float64 Precision Floor (~10^-14)`.
  - Ajouter la mention `[Lower is better]` sur les étiquettes d'axe des ordonnées pour guider intuitivement le lecteur.
- **Deliverables**:
  - `VisualizerAxes`: Panneaux 1.A, 1.B, 2.A, 2.B homogénéisés dans `visualizer.py`.
- **Process**: Édition ciblée dans `render_multi_integrator_invariants_dashboard`.

### Step 3: Régénération de la planche 300 DPI et synchronisation
- **Context**: La planche `dashboard_multi_integrator_invariants_300dpi.png` doit être actualisée avec le nouveau tableau et les axes homogénéisés.
- **Objective**: Exécuter la génération ciblée du dashboard multi-intégrateurs et synchroniser l'artefact dans le dossier brain.
- **Deliverables**:
  - `UpdatedDashboard`: `scratch/dashboard_multi_integrator_invariants_300dpi.png` et sa copie dans l'artifact directory.
- **Process**: Exécution d'un script scratch dédié ou appel de `PublicationVisualizer`.

### Step 4: Validation de non-régression et commit git
- **Context**: Garantir que la suite de tests unitaires reste à 100% passante.
- **Objective**: Exécuter `pytest tests/unit/` (151 tests) et acter le commit git.
- **Deliverables**:
  - `GitCommit`: Commit git sur `visualizer.py`.
- **Process**: `pytest tests/unit/`, `git add`, `git commit`.

## Acceptance Criteria:
### Step 1:
- Le tableau 3.B comporte exactement 6 colonnes sans la colonne "Order".
- Tous les 6 intégrateurs symplectiques affichent le statut vert [PASS].
### Step 2:
- Les 4 panneaux présentent un seuil identique à $-14.0$ et l'annotation `[Lower is better]`.
### Step 3:
- `dashboard_multi_integrator_invariants_300dpi.png` est régénéré avec succès à 300 DPI sans chevauchement.
### Step 4:
- 151/151 tests unitaires validés avec succès en $< 10$s.
