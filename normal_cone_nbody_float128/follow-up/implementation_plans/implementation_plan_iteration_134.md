# Implementation Plan of Iteration 134 : Régénération des Planches Scientifiques 300 DPI et des Données de la Web App Interactive
Created at: 2026-10-07T22:00:00+02:00

## Context:
1. **Explicit Request**:
   - « relance la génération des planche et des données de la webapp »
2. **Implicit Strategic Value**:
   - Mettre à jour l'ensemble des artefacts visuels scientifiques (tableaux de bord 300 DPI 10 panneaux) et le jeu de données `showcase_data.json` avec le solveur conormal désormais 100% conforme (flot hamiltonien complet $(\mathbf{q}, \mathbf{p})$ des correcteurs $\mathcal{C}_4/\mathcal{C}_6$, régularisation conormale des supernœuds avec $d_{\text{safe}} \ge 10^{-14}$, élimination du verrou $N=3$ et intégration du dwell-time dynamique $\tau_{\min}$).
   - Répercuter le nouveau record de précision énergétique sur Burrau ($4.49 \times 10^{-11}$) et sur l'ensemble des scénarios dans l'application web interactive `interactive_conormal_showcase.html`.
3. **Alignment Criteria**:
   - Exécution complète du script de benchmark et génération de données [`scratch/run_showcase_simulations.py`](file:///C:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scratch/run_showcase_simulations.py) sur les 4 scénarios canoniques (Figure-8, Burrau, Collision Radiale, 4-Corps hiérarchique).
   - Production des 3 planches scientifiques 300 DPI sans erreur et synchronisation dans le dossier des artefacts du brain.
   - Synchronisation de `showcase_data.json` et injection dans `interactive_conormal_showcase.html`.
   - Validation de l'intégrité de la webapp et des visualisations interactives (arbres de Jacobi, courbes d'énergie).

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des scénarios.
- **Règle d'Édition** : Utiliser exclusivement `replace_file_content` pour modifier le code source existant.

## Plan:
### Step 1: Exécution de la Suite de Simulation Multi-Scénarios
- **Contexte & Point de départ** : Script [`scratch/run_showcase_simulations.py`](file:///C:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scratch/run_showcase_simulations.py) configuré pour exécuter les 4 scénarios clés sur les 6 intégrateurs symplectiques.
- **Objectif** : Générer les données de trajectoires, diagnostics célestes et planches 300 DPI avec le solveur conormal pleinement harmonisé.
- **Livrable** : `Component Name`: `ShowcaseSimulationsExecution` produisant `scratch/showcase_data.json` et les 3 fichiers PNG `dashboard_*_300dpi.png`.
- **Processus** :
  - Lancer `python scratch/run_showcase_simulations.py` via `run_command`.
  - Suivre l'exécution asynchrone jusqu'à complétion.

### Step 2: Synchronisation des Planches 300 DPI dans le Brain Artifacts
- **Contexte & Point de départ** : Les planches 300 DPI sont générées dans `scratch/`.
- **Objectif** : Copier les 3 planches haute résolution dans le répertoire des artefacts du brain (`<appDataDir>\brain\<conversation-id>/`) pour affichage utilisateur.
- **Livrables** :
  - `dashboard_figure8_bm6_300dpi.png`
  - `dashboard_burrau_bm6_300dpi.png`
  - `dashboard_4body_bm6_300dpi.png`

### Step 3: Injection des Données Harmoniques dans la Web App Interactive
- **Contexte & Point de départ** : `interactive_conormal_showcase.html` héberge `SHOWCASE_DATA` inline pour une portabilité autonome sans serveur local.
- **Objectif** : Synchroniser `SHOWCASE_DATA` avec le contenu fraîchement généré de `showcase_data.json`.
- **Livrable** : `Component Name`: `InteractiveShowcaseUpdate` dans [`interactive_conormal_showcase.html`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/interactive_conormal_showcase.html).
- **Processus** :
  - Injecter le JSON dans la variable `SHOWCASE_DATA`.
  - Vérifier la validité de la syntaxe JavaScript et la persistance de l'arbre et du graphe d'énergie.

### Step 4: Validation Globale, Notes et Rapport Final
- **Contexte & Point de départ** : Artefacts régénérés.
- **Objectif** : Vérifier que 100% des indicateurs sont conformes, mettre à jour `notes.md` et produire `walkthrough.md`.
- **Processus** :
  - Vérifier l'accès aux planches et le bon chargement des données.
  - Consigner l'action dans `notes.md`.
  - Rédiger `walkthrough.md`.

## Acceptance Criteria:
- **Critère 1** : Les 3 planches 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`) sont régénérées et présentes dans le dossier des artefacts.
- **Critère 2** : `showcase_data.json` est régénéré avec les données du solveur conormal harmonisé.
- **Critère 3** : `interactive_conormal_showcase.html` est mis à jour et opérationnel.
