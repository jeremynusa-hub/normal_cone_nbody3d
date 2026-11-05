# Plan d'Implémentation de l'Itération 148 : Phase 5 — Découplage Conormal Feuille/Apex, Certification des 8 Benchmarks et Planches 300 DPI
Created at: 2026-10-08T11:28:00+02:00

## Context:
1. **Explicit Request**: Planifier rigoureusement l'Option B (découplage conormal feuille / apex dans l'ordonnanceur dyadique `DyadicScheduler`) afin d'aligner l'évaluation des échelles de temps et de péricentres sur le domaine réel du potentiel d'apex $V_{\mathrm{apex}}^{\mathrm{hier}}$, éliminer le sur-raffinement superflu à $k=17$ ($131\,072$ sous-pas) sur les paires de feuilles déjà résolues analytiquement par Stumpff-Sundman, exécuter la certification complète des 8 benchmarks canoniques sur les 6 schémas symplectiques, régénérer les 4 planches 300 DPI et resynchroniser le démonstrateur interactif.
2. **Implicit Strategic Value**: Rétablir l'isomorphisme mathématique absolu entre la décomposition hamiltonienne de Fulton--MacPherson ($H = H_{\mathrm{base}} + H_{\mathrm{leaf}} + V_{\mathrm{apex}}^{\mathrm{hier}}$) et la discrétisation multi-taux. Puisque l'attraction mutuelle des paires de feuilles est totalement exclue du potentiel d'apex et intégrée sans erreur de troncature de pas dans $H_{\mathrm{leaf}}$, la fréquence singulière interne d'une feuille au péricentre ($r_{\mathrm{peri}} \to 0$) ne doit en aucun cas dicter la subdivision temporelle des impulsions d'apex $\mathcal{K}^{\mathrm{apex}}$, dont le champ d'action est minoré par la distance d'isolement aux corps tiers ($d_0 > 0$). Ce découpage rigoureux divise par $\sim 100$ le temps de calcul des approches serrées tout en préservant l'invariance d'échelle théorique, le zéro adoucissement ($\varepsilon = 0$) et les invariants de Lie--Bargmann.
3. **Alignment Criteria**:
   - `DyadicScheduler.compute_dynamical_level` extrait les paires de feuilles actives de l'arbre de Jacobi $\mathcal{T}$ et restreint l'anticipation singulière du péricentre ($r_{\mathrm{peri}} = a(1-e)$) exclusivement aux paires non-feuilles qui contribuent effectivement au gradient d'apex.
   - Les paires de feuilles actives sont régulées par leur échelle dynamique orbitale naturelle ($\tau_{\mathrm{kin}}, \tau_{\mathrm{dyn}}$ à l'échelle de l'orbite), évitant l'écroulement artificiel vers $k=17$ lors des frôlements internes.
   - Suite à cette correction, les 8 benchmarks s'exécutent de façon fluide et complète sur les 6 intégrateurs en moins de 3 minutes au total.
   - Les 4 planches haute résolution 300 DPI et le showcase interactif sont régénérés avec succès.

## Guidelines:
- **Rigueur Mathématique & Isomorphisme** : Respecter scrupuleusement la partition du Hamiltonien $H = T_0 + T_{\mathrm{comp}} + H_{\mathrm{leaf}} + V_{\mathrm{apex}}^{\mathrm{hier}}$ (Chapitres 4 et 5). Aucun terme gravitationnel interne de feuille ne doit contaminer l'ordonnanceur d'apex.
- **Physique Non Adoucie Stricte** : Interdiction absolue de réintroduire un adoucissement $\varepsilon > 0$, un seuil de distance rigide $d_0 > 0$, ou un clamping arbitraire sur $k_{\max}$. La profondeur dyadique maximale théorique reste $k_{\max}=17$, accessible naturellement si une vraie collision non-feuille (d'apex) le requiert.
- **Invariance de Lie-Bargmann** : Préservation machine ($< 10^{-14}$) de la quantité de mouvement totale $\mathbf{P}$, du moment cinétique $\mathbf{L}$, et de l'énergie conjuguée $E$.
- **Non-Régression & Couverture** : 100 % de succès sur l'ensemble de la suite de tests unitaires (actuellement 139 tests).

## Plan:
### Step 1: Implémentation du Découplage Feuille/Apex dans DyadicScheduler
- **Context**: `normal_cone_nbody_float64/system/dyadic_scheduler.py` calcule actuellement $\tau_{\min}$ et anticipe $r_{\mathrm{peri}}$ indifféremment sur toutes les paires $i < j$, ignorant que les paires de feuilles sont intégrées analytiquement par Stumpff-Sundman et absentes de $V_{\mathrm{apex}}$.
- **Objective**: Mettre à jour `compute_dynamical_level` pour extraire les paires de feuilles de l'arbre $\mathcal{T}$ actif et n'appliquer l'anticipation singulière de péricentre qu'aux interactions non-feuilles (apex).
- **Deliverables**:
  - `normal_cone_nbody_float64/system/dyadic_scheduler.py`
- **Process**:
  1. Inspecter `tree.leaf_pairs` dans `compute_dynamical_level`.
  2. Construire l'ensemble des paires de feuilles canoniques `leaf_pair_set = set(tuple(sorted(n.body_indices)) for n in tree.leaf_pairs)`.
  3. Pour chaque paire $(i, j)$ :
     - Si $(i, j) \in \text{leaf\_pair\_set}$ : calculer $\tau_{\mathrm{kin}} = r_{ij} / v_{\mathrm{rel}}$ et $\tau_{\mathrm{dyn}} = \sqrt{r_{ij}^3 / GM}$, mais ne pas calculer la contraction vers $r_{\mathrm{peri}}$ singulier car le canal est exact.
     - Si $(i, j) \notin \text{leaf\_pair\_set}$ : appliquer l'anticipation képlérienne complète du péricentre $r_{\mathrm{peri}} = a(1-e)$ pour $v_{\mathrm{rad}} < 0$, assurant un raffinement dyadique maximal si un corps tiers pénètre l'apex avant bascule de Jacobi.
  4. Valider la cohérence pour le cas où `tree is None` (croisière globale).

### Step 2: Validation Unitaire et Tests de Régression Ciblés sur DyadicScheduler
- **Context**: Le comportement du planificateur dyadique doit être formellement vérifié sur des cas limites (binaire isolée serrée, passage au péricentre d'une feuille, approche d'un tiers perturbateur).
- **Objective**: Écrire des tests unitaires dédiés certifiant que le niveau dyadique reste modéré lors du passage au péricentre d'une feuille active et s'élève rigoureusement lors de l'intrusion d'un tiers non-feuille.
- **Deliverables**:
  - `tests/unit/system/test_dyadic_scheduler_leaf_decoupling.py`
- **Process**:
  1. Tester une binaire serrée en feuille : vérifier que $k$ reste borné ($k \le 6$) au péricentre $r \sim 10^{-4}$.
  2. Tester un intrus tiers non-feuille approchant : vérifier que $k$ augmente bien dynamiquement ($k \ge 8$) pour protéger l'apex kick.
  3. Exécuter `pytest tests/unit/system/`.

### Step 3: Exécution et Certification Complète des 8 Benchmarks Canoniques
- **Context**: Avec l'ordonnanceur dyadique corrigé, les 8 benchmarks canoniques peuvent être simulés sans blocage temporel artificiel.
- **Objective**: Évaluer et certifier l'ensemble des 8 scénarios (Figure-8, Burrau, Radial collision, 4-corps hiérarchique, Arenstorf, Sitnikov, Jovien, Collision quadruple) sur les 6 schémas symplectiques (Strang, Yoshida, BM4, PBM4, BM6, PBM6).
- **Deliverables**:
  - `scratch/run_full_benchmarks_certification.py` (actualisé et instrumenté)
  - Tableau de certification consolidé avec erreurs relatives maximales, erreurs finales et conservations $\Delta P, \Delta L$.
- **Process**:
  1. Exécuter le script de certification en mode unbuffered (`python -u scratch/run_full_benchmarks_certification.py`).
  2. Vérifier que chaque benchmark se termine en quelques secondes, y compris le 4-corps hiérarchique.
  3. Consigner les résultats de convergence et de conservation d'énergie.

### Step 4: Régénération des 4 Planches Haute Résolution 300 DPI
- **Context**: Les planches de publication doivent être produites avec les données réelles issues du code remédié et une calibration typographique irréprochable (aération des titres et diagrammes de Pareto complets).
- **Objective**: Générer les 4 tableaux de bord haute résolution 300 DPI et les copier dans le répertoire de restitution.
- **Deliverables**:
  - `scratch/dashboard_figure8_bm6_300dpi.png`
  - `scratch/dashboard_burrau_bm6_300dpi.png`
  - `scratch/dashboard_4body_bm6_300dpi.png`
  - `scratch/dashboard_multi_integrator_invariants_300dpi.png`
- **Process**:
  1. Exécuter `python scratch/run_showcase_simulations.py`.
  2. Vérifier visuellement et programmatiquement la présence des 4 fichiers PNG 300 DPI.
  3. Copier les planches dans l'espace d'artefacts `<appDataDir>\brain\<conversation-id>/`.

### Step 5: Resynchronisation du Démonstrateur Interactif et de son Dataset
- **Context**: L'interface HTML interactive [`interactive_conormal_showcase.html`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/interactive_conormal_showcase.html) et le fichier de données `showcase_data.json` doivent refléter les trajectoires et métriques réelles.
- **Objective**: Exporter les données complètes et vérifier la cohérence du visualiseur dynamique.
- **Deliverables**:
  - `scratch/showcase_data.json`
  - `interactive_conormal_showcase.html`
- **Process**:
  1. Valider l'export de `showcase_data.json` généré par `run_showcase_simulations.py`.
  2. Mettre à jour l'injection des trajectoires et des états d'arbres si nécessaire.

### Step 6: Validation Globale des Tests et Commit Git de Clôture de la Roadmap VI
- **Context**: Clôture officielle de la Phase 5 (Itération 148) et achèvement complet de la Roadmap VI.
- **Objective**: Valider 100 % des tests unitaires et enregistrer un commit git retraçant l'ensemble des remédiations.
- **Deliverables**:
  - Rapport pytest `100% passed` sans avertissement critique.
  - Git commit certifiant la Phase 5.
- **Process**:
  1. Exécuter `pytest tests/unit/`.
  2. Préparer et valider le commit git.

## Acceptance Criteria:
- **Critère Step 1 (Découplage Feuille/Apex)** : `DyadicScheduler.compute_dynamical_level` ne réduit pas $\tau_{\mathrm{pair}}$ via $r_{\mathrm{peri}}$ pour les paires appartenant à `tree.leaf_pairs`. Le péricentre d'une feuille active à $r \sim 10^{-4}$ produit un niveau dyadique $k \le 6$ (au lieu de $k=17$).
- **Critère Step 2 (Tests Unitaires Dédiés)** : Les nouveaux tests unitaires dans `tests/unit/system/test_dyadic_scheduler_leaf_decoupling.py` passent à 100 %, confirmant la différenciation feuille vs non-feuille.
- **Critère Step 3 (Exécution des 8 Benchmarks)** : L'ensemble des 48 simulations (8 benchmarks $\times$ 6 intégrateurs) s'exécute sans divergence, sans NaN, avec $\Delta \|\mathbf{P}\| < 10^{-14}$ et $\Delta \|\mathbf{L}\| < 10^{-14}$. L'erreur d'énergie sur le 4-corps hiérarchique sous BM6 est $< 10^{-7}$ en moins de 10 secondes d'exécution.
- **Critère Step 4 (Planches 300 DPI)** : Les 4 planches PNG à 300 DPI sont produites, vérifiées non vides, avec espacement typographique conforme et courbes de Pareto complètes pour tous les intégrateurs.
- **Critère Step 5 (Démonstrateur Interactif)** : `showcase_data.json` est valide, formaté et synchronisé avec `interactive_conormal_showcase.html`.
- **Critère Step 6 (Non-régression & Git)** : La suite complète `pytest tests/unit/` réussit à 100 % (au moins 141 tests). Le commit git de clôture est enregistré.

## Observations:
> [!IMPORTANT]
> Le découplage feuille / apex dans `DyadicScheduler` est le point d'ancrage exact de la théorie de Fulton--MacPherson : l'erreur de pas de temps de l'opérateur splitting provient exclusivement des crochets de Poisson avec le potentiel d'apex $V_{\mathrm{apex}}^{\mathrm{hier}}$. L'interaction interne d'une feuille n'étant pas dans $V_{\mathrm{apex}}$, elle ne doit pas piloter la subdivision des kicks d'apex.
>
> [!NOTE]
> La profondeur maximale théorique $k_{\max}=17$ reste strictement préservée dans les constantes et dans le scheduler pour les vraies collisions d'apex non résolues par une feuille de Jacobi.
