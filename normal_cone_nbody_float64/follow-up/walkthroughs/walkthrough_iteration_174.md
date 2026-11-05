# Walkthrough - Itération 174 : Solveur de Super-Nœuds Synchrones Découplé et Micro-Noyau Vectorisé sans Allocation

L'objectif de cette itération était de concevoir, implémenter et valider une variante optimisée découplée du solveur de super-nœuds synchrones ([`solvers/optimized_supernode.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/optimized_supernode.py)), dotée d'un micro-noyau sans allocation dynamique et avec déroulement scalaire pour petits amas ($k \le 4$), tout en maintenant une architecture strictement modulaire, optionnelle et 100 % rétro-compatible avec le solveur canonique de référence.

---

## 1. Réalisations Techniques

### A. Conception du Micro-Noyau sans Allocation (`solvers/optimized_supernode.py`)
- **Élimination des allocations mémoires dans la boucle microscopique** :
  - Dans la version de référence, chaque micro-pas de Leapfrog allouait des tenseurs 3D `dr[k, k, 3]`, des matrices de distances `d2[k, k]`, des masques diagonaux booléens et des matrices inverses `inv_d3`, soit plus de 4,5 millions d'appels à l'API C de NumPy sur la simulation à 4 corps.
  - Le nouveau module introduit des micro-noyaux spécialisés sans aucune allocation :
    - `_substep_cluster_k3` pour $k=3$ (3 paires scalaires explicites).
    - `_substep_cluster_k4` pour $k=4$ (6 paires scalaires explicites).
    - `_substep_cluster_general` pour $k \ge 5$ avec tampons de paires pré-alloués et calculs vectorisés 2D en place.
- **Factorisation et réutilisation géométrique directe** :
  - Les distances $d_{ij}$ et carrés $d^2_{ij}$ calculés pour le choix adaptatif du pas conormal sous horloge de Sundman ($dt_S = \rho_S d\tau_S$) sont directement réutilisés pour le premier demi-kick de Leapfrog sans ré-évaluation.
- **Inlining des interactions avec les spectateurs** :
  - Élimination de la fermeture Python `_compute_smooth_forces` recréée à chaque sous-pas ; les forces différentielles spectatoriales sont appliquées directement par paires en place sur le vecteur vitesse $v$.
- **Intégration directe du micro-noyau Kepler** :
  - Les paires serrées régularisées ($r_{\min} < 0.5 r_{\text{second}}$) et les super-nœuds binaires ($k=2$) sont résolus via `step_kepler_sundman_fast`.

### B. Branchement Modulaire Non Intrusif (`integrators/common.py` et `system/nbody_system.py`)
- **Gestion par ContextVar** :
  - Ajout de `_CURRENT_SUPERNODE_SOLVER` et des fonctions `set_current_supernode_solver` / `get_current_supernode_solver` dans [`integrators/common.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/common.py).
  - `step_supernode_flow` vérifie automatiquement si un solveur de super-nœuds optimisé est actif dans le contexte courant ; si oui, il délègue vers celui-ci, sinon il exécute la baseline canonique sans modification de signature ni rupture de compatibilité.
- **Intégration dans `NBodySystem`** :
  - Paramètre `supernode_evaluator: Optional[Any] = None` dans `NBodySystem.__init__`.
  - Activation automatique du jeton de contexte dans `NBodySystem.step()` avec réinitialisation garantie dans le bloc `finally:`.

### C. Validation Unitaire Exhaustive (`tests/unit/solvers/test_optimized_supernode.py`)
- 7 tests unitaires dédiés validés à 100 % :
  1. `test_optimized_supernode_binary_concordance` : concordance bit-exacte ($< 10^{-14}$) sur les clusters $k=2$.
  2. `test_optimized_supernode_k3_concordance_and_momentum` : concordance haute précision ($< 10^{-12}$) et conservation du moment linéaire ($\sum m_i v_i = 0$ à $< 10^{-15}$) et cinétique ($< 10^{-13}$) sur 3 corps.
  3. `test_optimized_supernode_k3_close_pair_concordance` : stabilité et équivalence avec paires serrées régularisées.
  4. `test_optimized_supernode_k4_concordance` : validation du déroulement à 6 paires sur 4 corps ($< 10^{-12}$).
  5. `test_optimized_supernode_homothetic_collapse` : symétrie d'effondrement radial en configuration centrale.
  6. `test_optimized_supernode_reversibility` : réversibilité temporelle microscopique $S(-dt) \circ S(dt) = \mathrm{Id}$.
  7. `test_optimized_supernode_contextvar_and_nbody_system` : cycle de vie du jeton de contexte et isolation dans `NBodySystem`.
- **Suite globale de tests** : 178/178 tests réussis en 7.62 s.

---

## 2. Benchmark Comparatif sur le 4-Corps Hiérarchique (Blanes-Moan BM6)

Le banc d'essai complet sur le 4-corps hiérarchique avec super-nœuds imbriqués a été exécuté avec la triple optimisation active :
- **Apex Kick** : `OptimizedScreenedApexPotential` (paires non-feuilles pré-mises en cache).
- **Kepler Leaf** : `OptimizedKeplerSolver` (Halley inliné, Stumpff $c_0..c_3$).
- **Supernode** : `OptimizedSupernodeSolver` (micro-noyau scalaire sans allocation).

### Synthèse Comparative des Performances :

| Métrique de Calcul | Baseline (It. 172) | Triple Optimisé (It. 174) | Gain Constaté |
| :--- | :---: | :---: | :---: |
| **Temps CPU Total Simulation** | 49.10 s | **39.43 s** | **-9.67 s (-19.7%)** |
| **Temps CPU Étape Super-nœuds** | 34.41 s | **24.73 s** | **-9.68 s (-28.1%)** |
| **Temps CPU Étape Apex Kick** | 4.17 s | **4.17 s** | Stable (déjà optimisé 2.6x) |
| **Temps CPU Étape Kepler** | 0.11 s | **0.11 s** | Stable (0.4% du temps total) |
| **Débit Global (steps/s)** | 16.3 steps/s | **20.3 steps/s** | **+24.5% de débit** |
| **Erreur d'Énergie $\Delta H / \|H_0\|$** | $8.221 \times 10^{-3}$ | **$8.221 \times 10^{-3}$** | **Concordance bit-exacte** |
| **Sous-pas Leapfrog Super-nœuds** | 204 358 | 204 358 | 100% exécutés sans allocation |
| **Évaluations Accélérations Internes** | 306 537 | 306 537 | 100% exécutées en micro-noyau |

### Analyse d'Efficacité :
- La suppression des allocations dynamiques et des fermetures Python dans les 204 358 sous-pas de Leapfrog a éliminé **près de 10 secondes pleines** de pure surcharge d'allocateur et de dispatch NumPy.
- La planche haute résolution à 300 DPI a été générée et synchronisée dans le dossier d'artefacts :
  [`cost_breakdown_hierarchical_4body_bm6_triple_optimized_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/cost_breakdown_hierarchical_4body_bm6_triple_optimized_300dpi.png).

---

## 3. Conformité et Enregistrement

1. **Intégrité physique** : Préservation inconditionnelle du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$) et des invariants de Lie-Bargmann.
2. **Modularité** : Le solveur de référence reste le solveur par défaut ; l'optimisation s'active de manière optionnelle via le paramètre ou le ContextVar.
3. **Commit Git** : Consigné sur `main` sous le hash `a7d883f`.
