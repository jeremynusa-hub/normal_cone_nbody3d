# Walkthrough of Iteration 167 : Correction et Forçage Rigoureux des Échelles Bilogarithmiques du Graphe de Pareto sur Tous les Modèles Célestes

## Summary:
L'Itération 167 a remédié de manière définitive à l'apparition de graduations scalaires pseudo-linéaires (`500, 1000, 1500, 2000` et `0.5, 1.0, 1.5 x 10^-11`) sur le panneau de Pareto (3.A). Cette anomalie provenait du comportement par défaut de Matplotlib lorsque la dynamique des données était inférieure à 2 décades (notamment sur la collision radiale et le problème de Burrau à pas unique). Deux actions structurelles ont été menées :
1. Le verrouillage inconditionnel de `ticker.LogLocator(base=10.0)` et `ticker.LogFormatterMathtext()` sur les deux axes dans `visualizer.py`, couplé à une plage canonique d'ordonnées [$10^{-16}, 10^0$] (harmonie directe avec le panneau 1.A).
2. Le calcul de véritables courbes de Pareto multi-points (balayage sur plusieurs pas de temps $\Delta t$) pour chacun des 4 scénarios célestes dans `scratch/run_showcase_simulations.py`, régénérant les 4 planches 300 DPI et l'application interactive HTML.

## Critical Discoveries:
- **Repli scalaire de Matplotlib** : Matplotlib bascule automatiquement sur un formateur scalaire lorsque le ratio max/min d'un axe logarithmique est modeste ($< 1.5$ décade). L'imposition explicite de `LogFormatterMathtext` est indispensable pour forcer l'affichage sous forme $10^k$.
- **Sensibilité du pas $\Delta t$ sur les binaires serrées du 4-corps** : L'inclusion initiale d'un pas trop grossier ($\Delta t = 0.03$) face à la période orbitale de la binaire interne ($T \approx 0.73$) entraînait un sur-découpage dyadique maximal ($2^{17} = 131\,072$ pas). En fixant la liste de test sur des pas réguliers et stables ($\Delta t \in [0.015, 0.01, 0.0075]$), les 6 schémas convergent rapidement et sans saturation.

## Realisation:
- **Timestamp de début** : 2026-10-09T09:56:00+02:00
- **Timestamp de fin** : 2026-10-09T10:41:00+02:00

### Step 1: Forçage formel du formatage logarithmique et des bornes canoniques dans `visualizer.py`
- **Éléments impactés** : [`normal_cone_nbody_float64/diagnostics/visualizer.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/diagnostics/visualizer.py)
- **Description** :
  - Configuration de `ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))` et `ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext())`.
  - Configuration de `ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=12))` et `ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())`.
  - Fixation canonique de l'axe des ordonnées : `y_min = 1e-16`, `y_max = 1e0` (ajustable vers le haut si l'erreur dépasse 1), garantissant l'affichage d'au moins 16 décades formelles ($10^{-16}, 10^{-14}, \dots, 10^0$).
  - Élargissement dynamique des bornes en abscisse si les points sont groupés, assurant l'affichage permanent de plusieurs puissances de dix ($10^2, 10^3, 10^4, 10^5$).

### Step 2: Calcul des courbes de compromis de Pareto multi-pas pour l'ensemble des 4 modèles dans `run_showcase_simulations.py`
- **Éléments impactés** : `scratch/run_showcase_simulations.py`
- **Description** :
  - Intégration du balayage multi-$\Delta t$ pour Burrau ($\Delta t \in [0.015, 0.01, 0.0075, 0.005]$).
  - Intégration du balayage multi-$\Delta t$ pour la Collision Radiale ($\Delta t \in [0.02, 0.01, 0.005, 0.002, 0.001]$).
  - Intégration du balayage multi-$\Delta t$ pour le 4-Corps Hiérarchique ($\Delta t \in [0.015, 0.01, 0.0075]$).
  - Transmission des dictionnaires `pareto_data` complets aux 4 appels de `vis.render_multi_integrator_invariants_dashboard`.
  - Mise à jour du template interactif `scratch/template_showcase.html` et de `scratch/build_showcase_html.py` pour afficher dynamiquement les courbes et les graduations logarithmiques selon le scénario actif.

### Step 3: Régénération des 4 planches 300 DPI et synchronisation des artefacts
- **Éléments impactés** :
  - [`dashboard_multi_integrator_invariants_figure8_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_multi_integrator_invariants_figure8_300dpi.png)
  - [`dashboard_multi_integrator_invariants_burrau_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_multi_integrator_invariants_burrau_300dpi.png)
  - [`dashboard_multi_integrator_invariants_radial_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_multi_integrator_invariants_radial_300dpi.png)
  - [`dashboard_multi_integrator_invariants_4body_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_multi_integrator_invariants_4body_300dpi.png)
  - [`interactive_conormal_showcase.html`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/interactive_conormal_showcase.html)
- **Description** :
  - Exécution complète du script `run_showcase_simulations.py` (code retour 0).
  - Copie automatique des planches 300 DPI vers le dossier d'artefacts `brain/`.
  - Recompilation de l'application HTML autonome interactive via `build_showcase_html.py`.

### Step 4: Validation de non-régression et enregistrement git
- **Éléments impactés** : Suite de tests unitaires, dépôt git local sur branche `main`.
- **Description** :
  - Exécution réussie des 151 tests unitaires en 6.68 secondes (`pytest tests/unit/`).
  - Enregistrement du commit git `e0e952a` sur `main`.

## Commit Details:
- **Hash** : `e0e952a`
- **Message** : `fix(visualizer): enforce strict LogFormatterMathtext and multi-point Pareto curves across all celestial models`

## Verification:
### Acceptance Verification:
- **Critère 3A.1 [Conforme]** : Les axes X et Y du panneau 3.A affichent exclusivement des étiquettes mathématiques $10^k$, éliminant tout repli scalaire linéaire.
- **Critère 3A.2 [Conforme]** : L'axe des ordonnées s'étend de $10^{-16}$ à $10^0$, assurant une correspondance visuelle immédiate avec le panneau 1.A.
- **Critère 3A.3 [Conforme]** : Chacune des 4 planches dispose de véritables courbes de Pareto continues multi-points illustrant la convergence des 6 schémas.
- **Critère 3A.4 [Conforme]** : 151/151 tests unitaires validés avec succès sans la moindre régression.

## Next Steps:
- Explorer l'ajout d'une ligne directrice asymptotique théorique $\mathcal{O}(h^2)$, $\mathcal{O}(h^4)$, $\mathcal{O}(h^6)$ sur le graphe de Pareto pour matérialiser visuellement la pente d'ordre attendue de chaque famille de schémas.
