# Implementation Plan of Iteration 167 : Correction et Forçage Rigoureux des Échelles Bilogarithmiques du Graphe de Pareto sur Tous les Modèles Célestes
Created at: 2026-10-09T09:56:00+02:00

## Context:
1. **Explicit Request**: "T'as oublié de mettre les axe ordonné et absicce en echelle log pour pareto sur certain diagramme".
2. **Implicit Strategic Value**: Assurer une lisibilité académique et éditoriale irréprochable des diagrammes de publication 300 DPI. Même lorsque les données sont groupées sur une plage étroite (par exemple sur la collision radiale où les erreurs oscillent entre $10^{-13}$ et $10^{-11}$), les axes ne doivent jamais se replier vers une numérotation scalaire pseudo-linéaire (`500, 1000, 1500` ou `0.5, 1.0, 1.5 x 10^-11`), mais doivent impérativement conserver la notation exponentielle formelle $10^k$, assortie de véritables courbes de compromis multi-points.
3. **Alignment Criteria**:
   - Absence totale de graduations scalaires linéaires sur le panneau 3.A de n'importe quel modèle.
   - Présence garantie de graduations exponentielles canoniques ($10^k$) sur X et Y via `LogLocator` et `LogFormatterMathtext`.
   - Étendue de l'axe Y harmonisée avec le panneau (1.A) (de l'ordre de $10^{-16}$ à $10^0$).
   - Production de vraies courbes Pareto multi-points pour chacun des 4 modèles célestes (Figure-8, Burrau, Collision radiale, 4-Corps).

## Guidelines:
- **Method**: Utiliser les localisateurs et formateurs standard de `matplotlib.ticker` (`LogLocator(base=10.0)` et `LogFormatterMathtext()`) pour interdire formellement à Matplotlib d'utiliser des étiquettes scalaires.
- **Red Lines**: Préservation inconditionnelle du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$). Interdiction d'utiliser `write_to_file` avec `overwrite: true` sur le code source de la bibliothèque.
- **Architecture**: `visualizer.py` encapsule la logique de rendu graphique, `scratch/run_showcase_simulations.py` orchestre les simulations et la génération des artefacts d'exposition.

## Plan:
### Step 1: Forçage formel du formatage logarithmique et des bornes canoniques dans `visualizer.py`
1. **Contexte initial**: `render_multi_integrator_invariants_dashboard` appelle `ax.set_xscale("log")` et `ax.set_yscale("log")`, mais sans spécifier de `LogFormatterMathtext` ni de bornes canoniques. Lorsque les données s'étendent sur moins de 2 décades, Matplotlib bascule en graduations scalaires.
2. **Objectif**: Verrouiller le formatage mathématique $10^k$ sur les deux axes et imposer une échelle d'ordonnées canonique [$10^{-16}, 10^0$] (ou ajustée si l'erreur dépasse $10^0$) pour que le panneau 3.A affiche toujours des décades explicites.
3. **Livrables**: `normal_cone_nbody_float64/diagnostics/visualizer.py`.
4. **Procédé**:
   - Importer `matplotlib.ticker as ticker`.
   - Configurer `ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))` et `ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext())`.
   - Configurer `ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=12))` et `ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())`.
   - Fixer `ax.set_ylim(bottom=1e-16, top=max(1.0, max_err * 2.0))`.

### Step 2: Calcul des courbes de compromis de Pareto multi-pas pour l'ensemble des 4 modèles dans `run_showcase_simulations.py`
1. **Contexte initial**: Seul le modèle Figure-8 effectuait un balayage sur plusieurs $dt$ pour produire des courbes continues dans `pareto_data`. Les 3 autres modèles n'exécutaient qu'un seul pas $dt$, se résumant à 6 points isolés sur une plage très restreinte.
2. **Objectif**: Calculer de véritables courbes de compromis travail-précision (3 à 4 pas de temps par schéma) pour Burrau, la Collision Radiale et le 4-Corps hiérarchique, sans impacter excessivement le temps de calcul.
3. **Livrables**: `scratch/run_showcase_simulations.py`.
4. **Procédé**:
   - Pour la collision radiale : balayage $dt \in [0.001, 0.002, 0.005, 0.01]$.
   - Pour Burrau : balayage $dt \in [0.002, 0.005, 0.01]$.
   - Pour le 4-corps hiérarchique : balayage $dt \in [0.005, 0.01, 0.02]$.
   - Structuration de `pareto_data` pour chacun et transmission au visualiseur.

### Step 3: Régénération des 4 planches 300 DPI et synchronisation des artefacts
1. **Contexte initial**: Les planches actuelles dans `scratch/` et `brain/` contiennent des panneaux 3.A avec affichage scalaire sur certains modèles.
2. **Objectif**: Régénérer les 4 planches 300 DPI avec les nouvelles courbes et les axes $10^k$ certifiés.
3. **Livrables**:
   - `scratch/dashboard_multi_integrator_invariants_figure8_300dpi.png`
   - `scratch/dashboard_multi_integrator_invariants_burrau_300dpi.png`
   - `scratch/dashboard_multi_integrator_invariants_radial_300dpi.png`
   - `scratch/dashboard_multi_integrator_invariants_4body_300dpi.png`
   - Leurs copies dans `brain/`.
4. **Procédé**: Exécution de `scratch/run_showcase_simulations.py` et copie vers l'annuaire d'artefacts.

### Step 4: Validation de non-régression et enregistrement git
1. **Contexte initial**: Le code doit demeurer 100% stable face à la suite de tests.
2. **Objectif**: Valider les 151 tests unitaires et consigner la modification par un commit git soigné.
3. **Livrables**: Rapport de test et commit git sur `main`.

## Acceptance Criteria:
- **Critère 3A.1**: Aucun des 4 diagrammes de Pareto ne présente d'étiquette scalaire linéaire (`500`, `1000`, `0.5 x 10^-11`, etc.) sur ses axes.
- **Critère 3A.2**: L'axe des ordonnées présente des graduations exponentielles canoniques de $10^{-16}$ à au moins $10^0$ sur chaque planche.
- **Critère 3A.3**: Chacune des 4 planches affiche de véritables courbes continues multi-points de compromis coût/erreur pour les 6 schémas symplectiques.
- **Critère 3A.4**: La suite de tests unitaires (151 tests) s'exécute à 100% avec succès.

## Observations:
> [!NOTE]
> Le repli constaté par l'utilisateur provient d'un comportement par défaut de Matplotlib (`LogLocator` bascule en `ScalarFormatter` dès lors que la plage de données dynamique est inférieure à environ 1.5 à 2 décades). En combinant un balayage multi-dt (qui élargit naturellement la plage dynamique sur plusieurs décades) avec l'imposition explicite de `ticker.LogFormatterMathtext` et d'une borne canonique [$10^{-16}, 10^0$], nous éliminons mathématiquement toute ambiguïté visuelle.
