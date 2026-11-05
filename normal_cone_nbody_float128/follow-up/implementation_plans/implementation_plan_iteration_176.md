# Implementation Plan of Iteration 176 : Génération des Planches de Diagnostic du 4-Corps avec Optimisations Triples/JIT pour les 6 Intégrateurs Symplectiques
Created at: 2026-10-09T13:07:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur demande : « lance la génération des planches pour 4-body avec toutes les optimisation précédentes mais en faisant varier l'intégrateur en faisant strang, yoshida, bm4, bm6, pbm4, pbm6 ».
2. **Implicit Strategic Value**: Fournir une cartographie exhaustive et comparative de la répartition du temps de calcul (Screened Apex, Kepler Leaf, Synchronous Supernodes, Dyadic Sub-stepping) et de la précision d'intégration hamiltonienne sur le problème à 4 corps hiérarchique hautement collisionnel sous les 6 schémas symplectiques de la monographie, tous équipés de la triple optimisation conormale (Apex vectorisé sans allocation, Kepler-Sundman à dérivée réutilisée et Super-nœuds compilés JIT / LLVM).
3. **Alignment Criteria**: Générer les 6 planches de publication 300 DPI à 4 panneaux, vérifier la préservation du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$) et des invariants de Lie-Bargmann pour chaque schéma, enregistrer les planches dans le répertoire d'artefacts du brain, et extraire un tableau comparatif rigoureux (débit en pas/s, temps d'étape, erreur relative maximale en énergie $\Delta H / |H_0|$).

## Guidelines:
- **Method**: Utiliser le moteur de profilage fine-grained [`scratch/run_cost_profiler.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scratch/run_cost_profiler.py) avec les drapeaux `--case hierarchical_4body`, `--apex-variant optimized`, `--kepler-variant optimized`, `--supernode-variant jit` pour chaque intégrateur dans (`strang`, `yoshida`, `bm4`, `bm6`, `pbm4`, `pbm6`).
- **Red Lines**: Aucun adoucissement physique ($\varepsilon = 0, d_0 = 0$). Les paramètres physiques et dynamiques ($dt=0.01$, $t \in [0, 8.0]$, masses, conditions initiales) sont rigoureusement identiques pour les 6 exécutions.
- **Architecture**: Conservation de la modularité ; chaque intégrateur bénéficie automatiquement des variantes optimisées via `NBodySystem(apex_evaluator=..., kepler_solver=..., supernode_evaluator=...)`.

## Plan:
### Step 1: Configuration et automatisation de la série de profilage 4-corps multi-intégrateurs
- **Context**: `scratch/run_cost_profiler.py` implémente déjà la commande unitaire `--integrator [strang|yoshida|bm4|bm6|pbm4|pbm6]` et la prise en charge de la variante `jit`.
- **Objective**: Mettre au point une fonction / boucle d'exécution séquencée ou parallélisée pour profiler les 6 schémas consécutivement sans interruption.
- **Deliverables**:
  - `scratch/run_cost_profiler.py` (ou script d'orchestration dédié `scratch/run_all_4body_profilers.py`): Prise en charge de `--integrator all` pour exécuter la batterie complète en une seule commande.
- **Process**:
  - Ajouter l'option `--integrator all` dans `scratch/run_cost_profiler.py` qui itère sur `["strang", "yoshida", "bm4", "bm6", "pbm4", "pbm6"]`.
  - Assurer la traçabilité des métriques et le nommage systématique des sorties 300 DPI.

### Step 2: Exécution des simulations profilées sous les 6 schémas (Strang, Yoshida, BM4, BM6, PBM4, PBM6)
- **Context**: Chaque simulation à 4 corps sur $t \in [0, 8.0]$ avec $dt=0.01$ compte 802 pas macro et jusqu'à 102 000 activations de sous-amas conormaux.
- **Objective**: Exécuter les 6 simulations complètes avec la triple optimisation conormale JIT (`apex: optimized`, `kepler: optimized`, `supernode: jit`).
- **Deliverables**:
  - Télémétrie complète et profils d'exécution pour les 6 intégrateurs.
- **Process**:
  - Lancer l'exécution des 6 intégrateurs avec chronométrage précis et extraction des profils temporels $r_{\min}(t)$ et $k(t)$.
  - Capturer les erreurs d'énergie $\Delta H / |H_0|$ et les temps CPU par étape (Apex, Kepler, Supernode, Dyadic, Atlas).

### Step 3: Génération des 6 planches haute résolution 300 DPI et miroitement dans les artefacts
- **Context**: Le visualiseur `PublicationVisualizer` produit des tableaux de bord 300 DPI à 4 panneaux (camembert de décomposition wall-clock, profil temporel dynamique de coût, histogramme log des opérations élémentaires, tableau récapitulatif).
- **Objective**: Générer et sauvegarder les 6 planches 300 DPI et les copier dans le dossier d'artefacts `brain/`.
- **Deliverables**:
  - `cost_breakdown_hierarchical_4body_strang_jit_optimized_300dpi.png`
  - `cost_breakdown_hierarchical_4body_yoshida_jit_optimized_300dpi.png`
  - `cost_breakdown_hierarchical_4body_bm4_jit_optimized_300dpi.png`
  - `cost_breakdown_hierarchical_4body_bm6_jit_optimized_300dpi.png`
  - `cost_breakdown_hierarchical_4body_pbm4_jit_optimized_300dpi.png`
  - `cost_breakdown_hierarchical_4body_pbm6_jit_optimized_300dpi.png`
- **Process**:
  - Sauvegarder les 6 fichiers PNG 300 DPI dans `scratch/` et les dupliquer dans l'artifact directory.

### Step 4: Analyse comparative des coûts computationnels et synthèse d'efficacité
- **Context**: Disposer des métriques comparatives permet d'évaluer le rapport coût/précision de chaque schéma sous la pile optimisée conormale.
- **Objective**: Rédiger une synthèse chiffrée comparative dans le walkthrough et les notes opérationnelles.
- **Deliverables**:
  - Tableau de synthèse comparatif des 6 schémas (temps total, temps super-nœuds, temps apex, débit pas/s, conservation $\Delta H / |H_0|$).

## Acceptance Criteria:
### Step 1:
- `scratch/run_cost_profiler.py` supporte `--integrator all` sans régression sur les intégrateurs individuels.
### Step 2:
- Les 6 simulations (`strang`, `yoshida`, `bm4`, `bm6`, `pbm4`, `pbm6`) terminent avec succès sans divergence, sans division par zéro et sans NaN.
### Step 3:
- Les 6 planches haute résolution 300 DPI sont générées et présentes dans le répertoire des artefacts du brain.
### Step 4:
- Les métriques comparatives (temps CPU, répartition par étape, précision machine) sont compilées et documentées.

## Observations:
> [!IMPORTANT]
> - Les intégrateurs traités PBM4 et PBM6 appliquent des correcteurs symplectiques canoniques aux frontières de macro-pas. Sous le mode d'intégration d'ombre (`in_shadow_mode`), le coût des correcteurs est télescopé entre macro-pas consécutifs.
> - La variante JIT compile les sous-pas microscopiques en code machine natif, réduisant drastiquement le temps total pour l'ensemble des 6 schémas.
