# Implementation Plan of Iteration 164 : Phase 5 — Certification Globale, Benchmarks et Synchronisation Démonstrateur
Created at: 2026-10-09T09:15:00+02:00

## Context:
1. **Explicit Request**: Exécution de la Phase 5 (Option 5A) de la Roadmap VIII : certification globale de l'ensemble des 8 benchmarks canoniques sur les 6 intégrateurs symplectiques, régénération des 4 dashboards de publication 300 DPI, synchronisation du démonstrateur interactif HTML5 et validation de non-régression à 100%.
2. **Implicit Strategic Value**: Parachever la Roadmap VIII en apportant la preuve empirique, visuelle et numérique irréfutable que les remédiations apportées aux discrépances de l'audit (enrichissement de Leibniz pour $\mathcal{C}_6$ en Lean 4 et Python, triggers spectraux de Hill $[2^{-1/3}, 2^{-1/6}]$, synchronisation de Skeel-Gear, découplage feuille/apex, section 5.6) fonctionnent en parfaite harmonie sans la moindre régression, sous zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$).
3. **Alignment Criteria**: Les 8 benchmarks doivent passer avec succès (statut [PASS] sur toutes les métriques célestes et invariants de Lie-Bargmann), les 4 planches 300 DPI doivent être générées sans collision de titre et avec les envelopes de Brouwer exactes, le showcase web doit être synchronisé avec les données bit-exactes, et les 151 tests unitaires doivent être validés à 100%.

## Guidelines:
- **Method**: Exécution séquentielle et contrôlée des scripts de benchmark et de visualisation, capture des métriques, et comparaison avec les baselines certifiées.
- **Red Lines**: Préservation absolue de la physique non adoucie ($d_0 = 0, \varepsilon = 0, d_{\mathrm{safe}} = 0$). Aucune modification ad-hoc de seuils pour masquer une anomalie.
- **Architecture**: Utilisation exclusive du pipeline canonique `normal_cone_nbody_float64` et des scripts de certification dédiés dans `scratch/`.

## Plan:

### Step 1: Exécution de la suite complète de certification des 8 benchmarks canoniques
- **Context**: Les 8 scénarios (Figure-8, Burrau, Collision radiale, Sitnikov, Arenstorf, Jovien, Quadruple collision, 4-corps hiérarchique) doivent être évalués sur les 6 intégrateurs (Strang S2, Yoshida S4, BM4, BM6, PBM4, PBM6) avec les constantes spectrales de Hill $[2^{-1/3}, 2^{-1/6}]$.
- **Objective**: Exécuter `scratch/run_full_benchmarks_certification.py` et certifier que l'ensemble des 48 configurations conserve l'énergie et les invariants de Lie-Bargmann ($< 10^{-14}$ sur orbites lisses, précision bornée sur quasi-collisions).
- **Deliverables**:
  - `CertificationReport`: Rapport d'exécution console et métriques d'erreur relative d'énergie et de moment cinétique.
- **Process**: Lancer le runner de benchmark global et enregistrer les scorecards de conformité.

### Step 2: Régénération haute résolution 300 DPI des 4 planches de publication
- **Context**: Les 4 dashboards de publication (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`, `dashboard_multi_integrator_invariants_300dpi.png`) doivent intégrer les dernières évolutions algorithmiques.
- **Objective**: Exécuter `scratch/run_showcase_simulations.py` pour régénérer les 4 planches PNG à 300 DPI avec typographie aérée (gap $\approx 2.0 \times F$), cône de bruit de Brouwer $\pm 5\sigma$ et table de certification à 5 colonnes.
- **Deliverables**:
  - `PublicationDashboards`: Les 4 fichiers PNG dans le répertoire de travail et copiés dans l'artifact directory.
- **Process**: Exécution du script de simulation et de rendu graphique, vérification des dimensions et de la lisibilité.

### Step 3: Exportation et synchronisation du démonstrateur interactif
- **Context**: Le démonstrateur web `interactive_conormal_showcase.html` et son fichier de données `showcase_data.json` doivent refléter fidèlement les trajectoires et dynamiques d'arbre de Jacobi.
- **Objective**: Exporter `showcase_data.json` et resynchroniser la constante `SHOWCASE_DATA` dans `interactive_conormal_showcase.html`.
- **Deliverables**:
  - `ShowcaseData`: `showcase_data.json` actualisé.
  - `InteractiveShowcase`: `interactive_conormal_showcase.html` synchronisé et autonome.
- **Process**: Injection des données JSON sérialisées dans la balise script du fichier HTML.

### Step 4: Validation intégrale de non-régression de la suite de tests unitaires
- **Context**: Toutes les modifications apportées durant la Roadmap VIII doivent garantir un taux de succès de 100% sur les tests.
- **Objective**: Exécuter `pytest tests/unit/` (151 tests) et s'assurer que zéro avertissement critique n'est émis.
- **Deliverables**:
  - `UnitTestSuite`: 151/151 tests unitaires validés avec succès.
- **Process**: Exécution de `python -m pytest tests/unit/`.

### Step 5: Bilan de clôture de la Roadmap VIII et commit git
- **Context**: Toutes les 14 discrépances de l'audit sont résolues.
- **Objective**: Dresser la synthèse finale de conformité, acter le commit git de clôture de la Phase 5 et mettre à jour la documentation maître.
- **Deliverables**:
  - `GitCommit`: Commit git final marquant l'achèvement de la Roadmap VIII.
  - `RoadmapStatus`: Checklist complète de la Roadmap VIII validée à 100%.
- **Process**: `git add`, `git commit`, mise à jour des artefacts de planning.

## Acceptance Criteria:
### Step 1:
- 48 configurations (8 benchmarks x 6 intégrateurs) exécutées sans crash, sans NaN et avec statut [PASS].
- Conservation des invariants de Lie-Bargmann $\Delta P, \Delta L < 10^{-14}$ sur orbites régulières.
### Step 2:
- 4 planches PNG 300 DPI générées avec taille de fichier > 1.5 Mo chacune.
- Marge d'en-tête aérée ($> 140$ px réels à 300 DPI) sans collision avec les titres de sous-graphes.
### Step 3:
- `showcase_data.json` et `interactive_conormal_showcase.html` synchronisés avec les 4 scénarios et les 6 intégrateurs.
### Step 4:
- 151/151 tests unitaires validés avec succès en $< 10$s.
### Step 5:
- Commit git enregistré et statut Roadmap VIII complété à 100%.

## Observations:
> [!IMPORTANT]
> Cette Phase 5 constitue le jalon de certification ultime de la Roadmap VIII. L'exécution de la suite complète de benchmarks (Option 5A) valide de façon définitive la robustesse de l'intégrateur conormal unifié face à l'ensemble des scénarios astrophysiques canoniques.
