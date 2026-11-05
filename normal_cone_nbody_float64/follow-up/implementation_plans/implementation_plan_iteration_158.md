# Plan d'Implémentation de l'Itération 158 : Re-génération Intégrale des Planches de Rendu 300 DPI et Resynchronisation du Démonstrateur Interactif
Created at: 2026-10-08T15:31:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur demande : « relance la création des planche de rendu ».
2. **Implicit Strategic Value**: Suite aux remédiations théoriques et numériques chirurgicales de l'Itération 156 (commit `f7e0514`) et à la certification unanime par le grand audit contradictoire à 6 sous-agents de l'Itération 157, les planches graphiques de publication haute résolution 300 DPI et le jeu de données du démonstrateur interactif doivent être régénérés pour refléter fidèlement le comportement du solveur immunisé contre l'overflow hyperbolique et doté du correcteur $K_6$ à parité impaire exacte.
3. **Alignment Criteria**:
   - Production des 4 planches haute résolution 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`, `dashboard_multi_integrator_invariants_300dpi.png`) sans collision de texte, avec aération typographique conforme (> 1.9x taille de police).
   - Export de `scratch/showcase_data.json` et resynchronisation de `interactive_conormal_showcase.html` dans le répertoire des artefacts du brain.
   - Synchronisation des 4 fichiers PNG dans le répertoire des artefacts du brain pour consultation directe par l'utilisateur.
   - Préservation stricte du zéro-adoucissement ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$).

## Guidelines:
- **Physique** : Zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$) sur l'ensemble des scénarios.
- **Rendu** : Résolution 300 DPI headless, respect strict des marges typographiques et des profils de tolérance KPI contextuels.
- **Intégrité** : Exécution via `scratch/run_showcase_simulations.py` et synchronisation dans les artefacts du brain.

## Plan:

### Step 1: Exécution des Simulations Comparatives Multi-Scénarios & Rendu des 4 Planches 300 DPI
- **Context**: `scratch/run_showcase_simulations.py` orchestre la simulation des 4 scénarios canoniques (Figure-8, Burrau, Collision radiale, 4-corps hiérarchique) sur l'ensemble des 6 intégrateurs symplectiques (Strang, Yoshida, BM4, PBM4, BM6, PBM6).
- **Objective**: Exécuter le pipeline de simulation avec le code certifié (commit `f7e0514`), calculer les KPI célestes et générer les 4 planches de publication 300 DPI dans `scratch/`.
- **Deliverables**:
  - `SimulationRun`: Exécution complète sans erreur des 4 scénarios.
  - `PublicationDashboards`: Génération des 4 fichiers PNG 300 DPI dans `scratch/`.
- **Process**: Lancer l'exécution du script `scratch/run_showcase_simulations.py` via Python dans l'environnement du projet.

### Step 2: Synchronisation des Planches 300 DPI et du Démonstrateur Interactif dans les Artefacts du Brain
- **Context**: Les artefacts visuels et les jeux de données doivent être directement accessibles et affichables dans le dossier des artefacts du brain (`<appDataDir>\brain\<conversation-id>/`).
- **Objective**: Copier les 4 planches PNG 300 DPI vers le dossier des artefacts du brain, synchroniser `showcase_data.json` et mettre à jour `interactive_conormal_showcase.html`.
- **Deliverables**:
  - `BrainDashboards`: Les 4 images PNG 300 DPI dans le dossier des artefacts.
  - `InteractiveShowcase`: `interactive_conormal_showcase.html` resynchronisé avec les nouvelles trajectoires.
- **Process**: Exécuter `scratch/sync_showcase_html.py` et copier les planches générées dans le brain.

### Step 3: Vérification de l'Intégrité Visuelle et Certification Finale
- **Context**: Contrôler les métriques de rendu, les tailles de fichiers et la conformité des KPIs.
- **Objective**: S'assurer que les 4 planches sont générées sans troncature, avec des tailles de fichiers valides (> 1 Mo chacune) et que les KPIs célestes sont au statut [PASS].
- **Deliverables**:
  - `VisualCertification`: Validation objective des planches et des données.
- **Process**: Vérifier l'existence et les propriétés des artefacts générés et certifier l'Itération 158.

## Acceptance Criteria:
### Step 1:
- `scratch/run_showcase_simulations.py` s'exécute à 100% sans exception ni avertissement d'overflow.
- Les 4 fichiers PNG 300 DPI sont créés dans `scratch/` avec des tailles > 1 Mo.
### Step 2:
- Les 4 fichiers PNG 300 DPI sont copiés dans le dossier des artefacts du brain.
- `interactive_conormal_showcase.html` est mis à jour avec le nouveau bloc `const SHOWCASE_DATA = ...`.
### Step 3:
- Les 4 planches sont vérifiées intègres et accessibles via les liens `file:///`.
- 100% des tests unitaires restent au statut PASS.

## Observations:
> [!NOTE]
> Le temps d'exécution de la suite complète des 4 scénarios (incluant les comparaisons multi-intégrateurs et le rendu graphique 300 DPI sous Matplotlib) est typiquement de l'ordre de 45 à 60 secondes.
