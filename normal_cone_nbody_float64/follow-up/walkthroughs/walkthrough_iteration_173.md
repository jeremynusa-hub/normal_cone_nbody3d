# Walkthrough - Itération 173 : Propagation du Flag d'Optimisation Kepler vers les Super-Nœuds Synchrones

L'objectif de cette itération était de propager de manière transparente et non intrusive le flag d'optimisation du solveur Kepler-Sundman (`OptimizedKeplerSolver`) vers le solveur de super-nœuds synchrones ([`solvers/supernode.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/supernode.py)), afin que les paires serrées isolées internes ($r_{\min} < 0.5 r_{\text{second}}$), les super-nœuds binaires ($k=2$) et les effondrements homothétiques bénéficient du micro-noyau inliné `step_kepler_sundman_fast` (micro-noyau Stumpff $c_0..c_3$, inlining de Halley sans fermeture, réutilisation directe de $\rho_{\text{final}}$ sans second appel de Stumpff).

---

## 1. Réalisations Techniques

### A. Aiguillage Dynamique dans `solvers/supernode.py`
- Implémentation de la fonction helper `_resolve_kepler_step(kepler_solver=None)` :
  - Interroge le paramètre `kepler_solver` ou, à défaut, le `ContextVar` d'exécution actif `get_current_kepler_solver()`.
  - Si un solveur optimisé est actif, retourne `step_kepler_sundman_fast`.
  - Sinon, retourne le solveur canonique de référence `step_kepler_sundman`.
- Propagation du paramètre optionnel `kepler_solver` et de la routine résolue dans :
  1. `_step_cluster_internal_substep` : régularisation analytique des paires serrées (ligne 202).
  2. `step_supernode_synchronous` (branche binaire $k=2$, ligne 296).
  3. `step_supernode_synchronous` (branche homothétique centrale, ligne 343).
  4. La boucle adaptative conormale interne sous temps déplié.

### B. Transmission via `integrators/common.py` et `solvers/optimized_kepler_sundman.py`
- Dans [`integrators/common.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/common.py) :
  - `step_supernode_flow` (et son alias `step_supernode_kick`) accepte désormais `kepler_solver: Optional[Any] = None` et résout automatiquement le `ContextVar` actif pour le transmettre à `step_supernode_synchronous`.
  - `step_drift_and_kepler` propage son `active_kepler` vers `step_supernode_flow`.
- Dans [`solvers/optimized_kepler_sundman.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/optimized_kepler_sundman.py) :
  - `step_supernode_flow` est appelé en lui transmettant explicitement `kepler_solver=self`.

### C. Validation Unitaire Dédiée
- Création du fichier de test [`tests/unit/solvers/test_supernode_optimized_kepler.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/solvers/test_supernode_optimized_kepler.py) couvrant 4 aspects critiques :
  1. `test_resolve_kepler_step_routing` : routage dynamique par défaut, explicite, via ContextVar, et réversion post-reset.
  2. `test_supernode_binary_concordance` : concordance bit-exacte ($< 10^{-14}$) sur les clusters binaires ($k=2$).
  3. `test_supernode_close_pair_regularization_concordance` : concordance ($< 10^{-13}$) et conservation stricte du moment linéaire ($\Delta P < 10^{-14}$) lors de la régularisation de paires serrées.
  4. `test_step_supernode_flow_propagation` : propagation transparente à travers la fonction de flux `step_supernode_flow`.

---

## 2. Résultats de Profilage et Analyse d'Impact sur le 4-Corps Hiérarchique

L'exécution du banc d'essai à 4 corps hiérarchique avec super-nœuds sous BM6 a révélé l'impact direct de cette propagation :

| Métrique Opérationnelle | Avant Propagation (It. 172) | Après Propagation (It. 173) | Gain Constaté |
| :--- | :---: | :---: | :---: |
| **Itérations de Halley** | 289 013 | **153 344** | **-135 669 itérations (-47%)** |
| **Évaluations de Stumpff** | 289 013 | **153 344** | **-135 669 évaluations (-47%)** |
| **Régularisations paires serrées** | 68 073 | 68 073 | 100% régularisées via micro-noyau |
| **Temps CPU Kepler Leaf** | 0.109 s | 0.113 s | Équivalent |
| **Temps CPU Super-nœuds** | 34.41 s | 34.35 s | $\approx -0.07$ s |
| **Erreur d'énergie $\Delta H/\|H_0\|$** | $8.221 \times 10^{-3}$ | $8.221 \times 10^{-3}$ | **Identique bit-exacte** |

### Enseignement Clé sur le Goulet d'Étranglement Restant :
La division par deux du nombre d'itérations de Halley/Stumpff démontre que la propagation vers les 68 073 paires régularisées fonctionne parfaitement et a éliminé tout le surcoût de ré-évaluation de Stumpff à convergence.
Cependant, le temps global des super-nœuds demeure à **$\approx 34.3$ secondes**. Cela confirme de façon rigoureuse notre diagnostic précédent :
> **Le véritable goulet d'étranglement des 34 secondes ne vient pas de la résolution de la paire képlérienne (qui ne consomme que $\approx 0.1$ s), mais des 306 537 évaluations d'accélérations internes (`_compute_internal_cluster_acc`) et des allocations répétées de tenseurs 3D NumPy (`dr[k,k,3]`, `np.zeros_like`, `np.copy`) dans les 204 358 sous-pas de Leapfrog.**

---

## 3. Conformité et Enregistrement

1. **Tests unitaires** : 171/171 tests réussis (100% de succès en 7.48 s).
2. **Physique conormale** : Zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$) et invariants de Lie-Bargmann préservés.
3. **Commit Git** : Consigné sur `main` sous le hash `521be28`.
