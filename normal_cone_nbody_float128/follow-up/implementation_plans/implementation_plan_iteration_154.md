# Implementation Plan of Iteration 154 : Re-Certification Complète des 8 Benchmarks, Planches 300 DPI et Appli Démonstrateur
Created at: 2026-10-08T12:45:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur poursuit la Roadmap VII vers sa dernière étape, la Phase 5 (Itération 154), consacrée à la re-certification empirique complète du moteur après l'ensemble des remédiations théoriques et numériques des Phases 1 à 4.
2. **Implicit Strategic Value**: Attester formellement que l'éradication des 4 divergences critiques (solveur Kepler-Sundman universel robuste à l'inversion et aux pas rétrogrades, régularisation canonique par isolation pairwise dans les super-nœuds, gradient complet de Leibniz pour $\mathcal{C}_6$ avec télescopage d'ombre et connexion de jauge, et requalification épistémique de $\mathrm{PBM}_4$) maintient 100% des benchmarks au vert [PASS], préserve les invariants de Lie-Bargmann à la précision machine float64, et produit des artefacts de publication (tableaux, planches 300 DPI, webapp démonstrateur) parfaitement resynchronisés.
3. **Alignment Criteria**:
   - 100% des 8 benchmarks canoniques validés sur les 6 intégrateurs (Strang, Yoshida, BM4, PBM4, BM6, PBM6).
   - Invariants de Lie-Bargmann préservés ($\|\Delta \mathbf{P}\| < 10^{-14}, \|\Delta \mathbf{L}\| < 10^{-14}$ sur orbites régulières).
   - Régénération des 4 tableaux de bord haute résolution 300 DPI sans chevauchement typographique.
   - Synchronisation exacte de `scratch/showcase_data.json` et de l'artefact `interactive_conormal_showcase.html`.
   - 100% des 148 tests unitaires au vert.

## Guidelines:
- **Method**:
  - Exécuter le script de certification globale `scratch/run_full_benchmarks_certification.py`.
  - Exécuter `scratch/run_showcase_simulations.py` pour régénérer les 4 planches 300 DPI et exporter `scratch/showcase_data.json`.
  - Exécuter `scratch/sync_showcase_html.py` pour mettre à jour `interactive_conormal_showcase.html`.
  - Lancer la suite complète `pytest tests/unit/`.
- **Red Lines**:
  - Zéro adoucissement artificiel ($\varepsilon = 0, d_0 = 0$).
  - Aucun échec ou dégradation de précision sur les benchmarks.
  - Typographie et espacements des planches 300 DPI rigoureusement préservés.
- **Architecture**:
  - Respect de l'architecture modulaire en 7 couches de `normal_cone_nbody_float64`.

## Plan:
### Step 1: Exécution de la Suite Complète des 8 Benchmarks sur les 6 Intégrateurs
- **Context**: Le script `scratch/run_full_benchmarks_certification.py` évalue les 8 scénarios astrophysiques canoniques sur l'ensemble des 6 intégrateurs symplectiques.
- **Objective**: Quantifier les erreurs relatives d'énergie et la conservation des invariants de Lie-Bargmann sur les 48 configurations.
- **Deliverables**:
  - `GlobalBenchmarkCertification`: Rapport d'exécution et tableau de bord comparatif des 48 simulations.
- **Process**: Lancer `python scratch/run_full_benchmarks_certification.py` avec `$env:PYTHONPATH="."`.

### Step 2: Régénération des 4 Planches de Publication 300 DPI
- **Context**: Les planches haute résolution 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`, `dashboard_multi_integrator_invariants_300dpi.png`) doivent refléter le solveur sécurisé avec le gradient de Leibniz $\mathcal{C}_6$ et l'isolation pairwise des super-nœuds.
- **Objective**: Générer les 4 figures de publication avec une résolution de 300 DPI et des KPIs à 100% [PASS].
- **Deliverables**:
  - `PublicationDashboards300DPI`: 4 fichiers PNG dans le répertoire de travail et le brain.
- **Process**: Lancer `python scratch/run_showcase_simulations.py` avec `$env:PYTHONPATH="."`.

### Step 3: Resynchronisation des Données du Démonstrateur Interactif
- **Context**: `scratch/showcase_data.json` contient l'ensemble des données multi-intégrateurs et la topologie des arbres de Jacobi. `interactive_conormal_showcase.html` doit intégrer ces nouvelles données.
- **Objective**: Synchroniser `SHOWCASE_DATA` dans le démonstrateur web standalone HTML5/SVG/Canvas.
- **Deliverables**:
  - `SynchronizedInteractiveShowcase`: `interactive_conormal_showcase.html` mis à jour.
- **Process**: Exécuter `python scratch/sync_showcase_html.py`.

### Step 4: Validation de l'Ensemble de la Suite des Tests Unitaires (148/148)
- **Context**: La suite de tests unitaires garantit la non-régression de l'ensemble de la bibliothèque.
- **Objective**: Vérifier 100% de succès sur 148 tests.
- **Deliverables**:
  - `UnitTestsCertification`: 148/148 tests réussis.
- **Process**: Exécuter `$env:PYTHONPATH="."; pytest tests/unit/`.

### Step 5: Git Commit de Phase 5 et Clôture de la Roadmap VII
- **Context**: Phase 5 et Roadmap VII achevées.
- **Objective**: Enregistrer le commit atomique de clôture.
- **Deliverables**:
  - `GitCommitPhase5`: Commit git documenté sur `main`.
- **Process**: `git add` et `git commit`.

## Acceptance Criteria:
### Step 1:
- Les 8 benchmarks s'exécutent sans divergence ni NaN sur les 6 intégrateurs.
- Figure-8 atteint une précision d'énergie $< 10^{-12}$ sur BM6.
- Les collisions radiale et quadruple maintiennent une erreur finie sous régularisation sans division par zéro.
### Step 2:
- Les 4 fichiers PNG 300 DPI sont générés, valides, et de taille $> 1$ Mo.
- Tous les voyants de conformité affichent [PASS].
### Step 3:
- `scratch/sync_showcase_html.py` s'exécute avec la mention `SUCCESS`.
- `interactive_conormal_showcase.html` contient le jeu de données complet incluant les 6 intégrateurs et les arbres hiérarchiques.
### Step 4:
- 148 tests sur 148 passent avec succès en $< 10$s.
### Step 5:
- Commit git atomique enregistré sur `main`.

## Observations:
> [!IMPORTANT]
> Cette Phase 5 scelle l'intégralité de la Roadmap VII en démontrant que le solveur unifié conormal atteint des performances inégalées sans aucun compromis d'adoucissement ($\varepsilon = 0, d_0 = 0$) sur l'ensemble du corpus de test astrophysique.
