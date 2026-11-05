# Implementation Plan of Iteration 138 : Régénération des Planches de Publication 300 DPI et des Données de la Webapp avec la Hiérarchie C6/C8
Created at: 2026-10-08T08:43:00+02:00

## Context:
1. **Explicit Request**: 
   - L'utilisateur demande : "Régénérer plan et app", commandant la réexécution des simulations complètes, la mise à jour des planches scientifiques 300 DPI et la resynchronisation des données de l'application interactive `interactive_conormal_showcase.html`.
2. **Implicit Strategic Value**:
   - Aligner l'ensemble des artefacts visuels, métriques célestes et télémétries interactives sur l'état le plus abouti de la bibliothèque :
     - Éradication des seuils dyadiques arbitraires via le scheduler dynamique continu invariant d'échelle $\tau_{\text{enc}} = \min(\tau_{\text{kin}}, \tau_{\text{dyn}}, \tau_{\text{peri}})$.
     - Hiérarchie d'écrantage Hermite symétrique $\mathcal{C}^6$ (degré 13) pour PBM4/BM6 et $\mathcal{C}^8$ (degré 17) pour PBM6, éliminant l'extinction exponentielle IEEE 754.
     - Documenter de manière probante le gain d'exactitude sur le système 4-corps hiérarchique ($1.14 \times 10^{-8}$ sur BM6, réduction d'un facteur 17.5 du saut d'énergie péricentrique).
3. **Alignment Criteria**:
   - Exécution intégrale et sans régression de `scratch/run_showcase_simulations.py` sur les 4 scénarios canoniques (Figure-8, Burrau, Collision radiale avec rebond McGehee, 4-corps hiérarchique avec 15 transitions d'atlas).
   - Génération des 3 planches haute résolution 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`) avec 100% de certifications KPI [PASS].
   - Injection automatique des nouvelles données dans `interactive_conormal_showcase.html` et synchronisation dans le brain.

## Guidelines:
- **Method**: Utiliser `scratch/run_showcase_simulations.py` pour orchestrer les calculs comparatifs des 6 intégrateurs (Strang S2, Yoshida S4, BM4, PBM4, BM6, PBM6), puis synchroniser l'application HTML via `scratch/sync_showcase_html.py`.
- **Red Lines**: 
  - Physique sans adoucissement stricte ($\varepsilon = 0, d_0 = 0$) sur toutes les simulations.
  - Aucune modification du manuscrit LaTeX (`redaction/`).
  - Vérification de la non-régression de l'ensemble des 131 tests unitaires.
- **Architecture**:
  - `scratch/run_showcase_simulations.py` $\to$ calcul et rendu matplotlib 300 DPI.
  - `scratch/sync_showcase_html.py` $\to$ injection JSON dans l'artefact HTML du brain.
  - Synchronisation des fichiers PNG et JSON vers l'espace de stockage des artefacts du brain.

## Plan:
### Step 1: Exécution des Simulations Comparatives Multi-Scénarios
- **Context**: La bibliothèque intègre les nouveaux profils $\mathcal{C}^6$ et $\mathcal{C}^8$ ainsi que le scheduler dyadique continu.
- **Objective**: Exécuter `scratch/run_showcase_simulations.py` pour simuler les 4 scénarios canoniques avec les 6 schémas symplectiques et générer `showcase_data.json` ainsi que les 3 planches 300 DPI.
- **Deliverables**:
  - `ShowcaseSimulationPayload`: `scratch/showcase_data.json` recalculé.
  - `PublicationDashboards`: `scratch/dashboard_figure8_bm6_300dpi.png`, `scratch/dashboard_burrau_bm6_300dpi.png`, `scratch/dashboard_4body_bm6_300dpi.png`.
- **Process**: Lancement en tâche de commande PowerShell supervisée.

### Step 2: Synchronisation de la Webapp et des Artefacts du Brain
- **Context**: Les données numériques et planches générées doivent être injectées dans l'application web interactive et copiées dans le brain.
- **Objective**: Exécuter `scratch/sync_showcase_html.py` pour injecter le JSON dans `interactive_conormal_showcase.html`, et copier les planches PNG et `showcase_data.json` vers le répertoire des artefacts (`C:\Users\jerem\.gemini\antigravity\brain\59313b1a-0dde-45d6-a476-1bcb513cf82e\`).
- **Deliverables**:
  - `InteractiveShowcaseApp`: `interactive_conormal_showcase.html` resynchronisé.
  - `BrainArtifactCopies`: `showcase_data.json` et les 3 planches PNG mis à jour dans le brain.
- **Process**: Exécution du script de synchronisation et copie de fichiers.

### Step 3: Contrôle Qualité et Validation d'Intégrité
- **Context**: S'assurer de la validité de toutes les planches et de la cohérence de l'interface interactive.
- **Objective**: Vérifier les métriques clés (erreur max $\sim 10^{-8}$ sur 4-corps, $\sim 10^{-11}$ sur Burrau, $\sim 10^{-12}$ sur Figure-8), vérifier la présence et la taille non nulle des planches 300 DPI, et valider la suite de tests unitaires (131/131 passés).
- **Deliverables**:
  - `QualityValidationReport`: Certification des KPIs et validation des tests.
- **Process**: Inspection de sortie et exécution de pytest.

## Acceptance Criteria:
### Step 1:
- Simulation des 4 scénarios achevée avec succès sans division par zéro ni crash.
- `scratch/showcase_data.json` généré avec une taille > 5 Mo contenant les 6 intégrateurs.
### Step 2:
- `interactive_conormal_showcase.html` mis à jour avec `const SHOWCASE_DATA` synchronisé.
- Les 3 planches PNG 300 DPI copiées dans le brain avec des tailles valides (> 1.5 Mo chacune).
### Step 3:
- BM6 sur 4-corps hiérarchique certifié avec une erreur maximale $\le 1.2 \times 10^{-8}$.
- 131/131 tests unitaires validés avec succès sous `pytest tests/unit/`.

## Observations:
> [!IMPORTANT]
> Cette régénération consacre la transition complète vers la théorie d'écrantage Hermite symétrique $\mathcal{C}^6 / \mathcal{C}^8$ et l'ordonnancement dynamique continu, éliminant définitivement les artefacts de calcul et les plateaux d'underflow IEEE 754.
