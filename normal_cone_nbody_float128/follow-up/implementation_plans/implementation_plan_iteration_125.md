# Implementation Plan of Iteration 125 : Implémentation Intégrale des Correctifs d'Audit P0 & P1 et Certification de Résilience
Created at: 2026-10-07T19:12:45+02:00

## Context:
1. **Explicit Request**: « evidemment, corrige tout »
2. **Implicit Strategic Value**: Résoudre définitivement l'ensemble des 7 anomalies structurelles et failles diagnostiques mises au jour par le triple audit sous-agentique contradictoire de l'Itération 124, éliminer le blowup numérique à $t \approx 3.3$ sur les super-nœuds, rétablir l'écrantage d'apex $\mathcal{C}^{1,1}$ dans les kicks symplectiques, calibrer les profils de benchmarks et consolider la matrice de tests unitaires avec des scénarios dynamiques réels.
3. **Alignment Criteria**: Suppression de toutes les erreurs de divergence, rétablissement de l'horloge conormale dépliée $\mathrm{d}t_S = \chi_S^2\,\mathrm{d}\tau_S$ dans `supernode.py`, dérive cinématique exacte du centre de masse des super-nœuds dans `common.py`, activation de `include_screened_monopole=True`, généralisation du détecteur de Hill à tout $N \ge 4$, et passage à 100 % [PASS] sur tous les tests unitaires et benchmarks célestes.

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des modules.
- **Granularité Sémantique** : Utiliser exclusivement `replace_file_content` pour modifier les fichiers sources existants.
- **Isomorphisme Mathématique** : Rétablir la stricte conformité aux équations des Chapitres 3, 4 et 5 de la monographie.

## Plan:
### Step 1: Correctifs P0 du Solveur de Super-Nœud et Dérive Cinématique (`solvers/supernode.py`, `integrators/common.py`, `system/nbody_system.py`)
- **Matériel de départ** : [`solvers/supernode.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/supernode.py), [`integrators/common.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/common.py), [`system/nbody_system.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/system/nbody_system.py).
- **Objectif** :
  1. Dans `solvers/supernode.py`, implémenter la véritable horloge conormale de Sundman dépliée $\mathrm{d}t_S = \chi_S^2\,\mathrm{d}\tau_S$ adaptant le sous-pas $\Delta t_k = \rho_S(t_k)\,\Delta \tau_S$, et préserver la résolution analytique exacte pour la binaire interne la plus serrée.
  2. Dans `integrators/common.py`, faire translater le barycentre $\mathbf{q}_{\mathrm{cdm}} \leftarrow \mathbf{q}_{\mathrm{cdm}} + \mathrm{d}t \cdot (\mathbf{p}_{\mathrm{cdm}} / M_S)$ et ne pas geler les nœuds composites parents.
  3. Dans `system/nbody_system.py` et `common.py`, gérer les super-nœuds sous forme de liste de clusters disjoints `List[Tuple[int, ...]]` au lieu de les fusionner dans un unique `Set[int]`.
- **Livrable** : Code source corrigé et validé sur le cas 4-corps à $t=3.3$.
- **Processus** : Éditions ciblées via `replace_file_content` et vérification interactive.

### Step 2: Réactivation de l'Apex $\mathcal{C}^{1,1}$ et Extension du Détecteur de Hill (`potential/screened_apex.py`, `potential/supernode_detector.py`)
- **Matériel de départ** : [`potential/screened_apex.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screened_apex.py), [`potential/supernode_detector.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/supernode_detector.py).
- **Objectif** :
  1. Dans `screened_apex.py`, activer `include_screened_monopole=True` dans `step_apex_kick` pour déployer le potentiel d'apex lisse $\mathcal{C}^{1,1}$ $W_v = \Theta(u_v) \frac{G M_L M_R}{\|\mathbf{r}_v\|}$ et ses gradients analytiques complets.
  2. Dans `supernode_detector.py`, supprimer la restriction `if state.N <= 3` et appliquer universellement le critère de compacité de Marchal--Sundman $\delta_{\mathrm{cluster}}^{\mathrm{nat}}(S)$ pour tout sous-amas $|S| \ge 3$ jusqu'à $N$ corps.
- **Livrable** : Régularité $\mathcal{C}^{1,1}$ opérationnelle et détection de super-nœud validée sur collision quadruple.
- **Processus** : Éditions chirurgicales via `replace_file_content`.

### Step 3: Calibrage Contextuel des Benchmarks et Correcteurs Symplectiques (`benchmarks/runner.py`, `benchmarks/registry.py`, `integrators/correctors.py`)
- **Matériel de départ** : [`benchmarks/runner.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/benchmarks/runner.py), [`benchmarks/registry.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/benchmarks/registry.py), [`integrators/correctors.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/correctors.py).
- **Objectif** :
  1. Dans `runner.py`, passer le profil contextuel de chaque benchmark (`"chaotic_encounter"`, `"hierarchical_multiscale"` ou `"smooth_periodic"`) lors de l'appel à `compute_celestial_kpis()`.
  2. Dans `registry.py`, annoter chaque benchmark avec son profil canonique.
  3. Dans `correctors.py`, ajouter la dérive d'impulsion $\mathbf{p} \leftarrow \mathbf{p} - \delta \nabla_{\mathbf{q}} K$ pour préserver rigoureusement la forme symplectique de Darboux.
- **Livrable** : Validation 100 % [PASS] sur `BenchmarkRunner.run()` pour tous les problèmes canoniques.
- **Processus** : Mises à jour ciblées via `replace_file_content`.

### Step 4: Renforcement de la Suite de Tests Unitaires et Validation Globale (`tests/unit/`)
- **Matériel de départ** : [`tests/unit/system/test_atlas_manager.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/system/test_atlas_manager.py), [`tests/unit/solvers/test_supernode.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/solvers/test_supernode.py), [`tests/unit/solvers/test_mcgehee.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/solvers/test_mcgehee.py).
- **Objectif** :
  1. Ajouter un test d'intégration dynamique subissant un basculement d'atlas effectif avec conservation d'énergie post-transition.
  2. Ajouter un test unitaire validant la traversée d'un super-nœud sans blowup sur le cas 4-corps à $t=3.3$.
  3. Valider l'intégralité de la suite `pytest tests/unit/` (100 % passants).
- **Livrable** : Suite de tests enrichie, robuste et certifiée à 100 %.
- **Processus** : Ajout de tests ciblés et exécution via `run_command`.

## Acceptance Criteria:
- **Critère 1** : L'instabilité numérique à $t \approx 3.3$ sur le 4-corps est éliminée ($\Delta H / |H_0| < 10^{-4}$ au lieu de $+78\,000\,\%$).
- **Critère 2** : Les super-nœuds multiples sont traités indépendamment sans mélange barycentrique.
- **Critère 3** : Le potentiel d'apex lisse $\mathcal{C}^{1,1}$ $W_v$ est pleinement actif lors des kicks.
- **Critère 4** : `BenchmarkRunner.run()` certifie les benchmarks chaotiques et multi-échelles avec leurs profils adaptés.
- **Critère 5** : 100 % de succès sur l'ensemble des tests unitaires (`pytest tests/unit/`).

## Observations:
> [!IMPORTANT]
> L'intégration du temps déplié de Sundman $\mathrm{d}t_S = \chi_S^2\,\mathrm{d}\tau_S$ au sein du super-nœud doit préserver l'invariance exacte du moment linéaire total $\mathbf{P}_S$ et du moment cinétique interne $\mathbf{L}_S$ au seuil float64 ($< 10^{-14}$).
> [!WARNING]
> La dérive du centre de masse du super-nœud dans `common.py` doit être rigoureusement synchronisée avec les sous-pas d'intégration pour garantir la conservation globale de l'invariance galiléenne $\mathbf{G}(t)$.
