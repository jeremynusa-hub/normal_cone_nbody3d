# Implementation Plan of Iteration 133 : Mise en Conformité Théorique Complète du Code Source (Supernœuds Conormaux, Correcteurs Symplectiques et Atlas Dynamique)
Created at: 2026-10-07T21:50:00+02:00

## Context:
1. **Explicit Request**:
   - « Met en conformité le code »
   - Consécutif à l'audit contradictoire de l'Itération 132 et au relevé des 8 défauts et écarts théoriques répertoriés dans [`audit_defauts_et_non_conformites.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/audit_defauts_et_non_conformites.md) et [`confrontation_couples_theoriques_vs_implementation.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/confrontation_couples_theoriques_vs_implementation.md).
2. **Implicit Strategic Value**:
   - Élever la bibliothèque `normal_cone_nbody_float64/` au statut d'isomorphisme mathématique absolu avec les 9 couples théoriques formels (*Problem Formulation* / *Proof*) de Cas 6 et le manuscrit LaTeX.
   - Éliminer les approximations ad hoc et les masquages de bugs (Leapfrog brut dans les supernœuds, verrouillage artificiel $N \ge 4$ sur les collisions triples, omission de la hessienne dans les correcteurs symplectiques $\mathcal{C}_{2p}$).
   - Préserver la parfaite stabilité numérique sous zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$) sur l'ensemble des scénarios astrophysiques.
3. **Alignment Criteria**:
   - Les correcteurs $\mathcal{C}_4$ et $\mathcal{C}_6$ appliquent le flot hamiltonien complet de $K$ (dérive spatiale $\dot{q}$ et impulsion de force $\dot{p}$ via la hessienne directionnelle).
   - Le solveur de supernœud `solvers/supernode.py` désingularise les collisions internes de sous-paires sous l'horloge dépliée de Sundman $\mathrm{d}t_S = \rho_S \mathrm{d}\tau_S$ sans singularité newtonienne brute.
   - Le détecteur de supernœuds `potential/supernode_detector.py` prend en charge sans restriction artificielle les collisions triples à $N=3$ et préserve la symétrie des amas tétraédriques à $N=4$.
   - L'`AtlasManager` active un `tau_min` dynamique basé sur l'échelle de temps de traversée physique $\frac{16}{15} \frac{r_{\min}}{v_{\max}}$.
   - 100% de la suite de tests unitaires (125 tests) passe sans avertissement.

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des scénarios.
- **Règle d'Édition** : Utiliser exclusivement `replace_file_content` pour modifier le code Python existant, ne jamais utiliser `write_to_file` avec `overwrite: true`.
- **Non-Régression Strict** : Préserver la précision record sur Burrau ($9.5 \times 10^{-11}$), Figure-8 ($3.46 \times 10^{-12}$) et Radial Collision ($7.34 \times 10^{-12}$).

## Plan:
### Step 1: Mise en Conformité des Correcteurs Symplectiques C4 et C6 (Couples 7 & 9)
- **Contexte & Point de départ** : `integrators/correctors.py` applique uniquement la dérive spatiale $q \leftarrow q + \delta \nabla_p K$ et pose $p' = p$, omettant le kick en impulsion $\dot{p} = -\nabla_q K$ prescrit par les théorèmes 7.1 et 9.1.
- **Objectif** : Implémenter l'action complète du champ hamiltonien $X_K$ en intégrant la variation d'impulsion $\delta \mathbf{p}_i = \kappa \sum_j \mathbf{H}_{ij} \mathbf{v}_j$ via dérivation directionnelle exacte de la force d'apex le long de la vitesse $\mathbf{v} = \mathbf{p} / \mathbf{m}$.
- **Livrable** : `Component Name`: `SymplecticCorrectorsCompleteFlow` dans [`integrators/correctors.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/correctors.py).
- **Processus** :
  - Dans `CorrectorC4.apply` et `CorrectorC6.apply`, calculer le gradient directionnel de force $D_{\mathbf{v}} \mathbf{F}_i \approx \frac{\mathbf{F}_i(\mathbf{q} + \eta \mathbf{v}) - \mathbf{F}_i(\mathbf{q} - \eta \mathbf{v})}{2\eta}$ avec $\eta = 10^{-7} \max(1, \|\mathbf{q}\|)$.
  - Mettre à jour l'impulsion $\mathbf{p}_i' = \mathbf{p}_i + \text{step\_coeff} \cdot D_{\mathbf{v}} \mathbf{F}_i \cdot m_i$ pour former une composition de Lie-Trotter symétrique d'ordre 2 pour le flot de correction.
  - Vérifier la réversibilité analytique et l'annulation exacte $\mathcal{C}^{-1} \circ \mathcal{C} = \mathrm{Id}$.

### Step 2: Régularisation Conormale du Solveur Interne des Super-Nœuds (Couple 6)
- **Contexte & Point de départ** : `solvers/supernode.py` utilise un Leapfrog cartésien direct avec `inv_d3 = G / (d2 * d)`, non régularisé en cas de collision interne.
- **Objectif** : Implémenter une intégration désingularisée sous l'horloge dépliée conormale de Sundman $\mathrm{d}t_S = \chi_S^2 \mathrm{d}\tau_S = \rho_S \mathrm{d}\tau_S$ avec régularisation conormale des collisions internes ou isolement de la paire interne dominante sous variables de Stumpff.
- **Livrable** : `Component Name`: `ConormalSundmanSupernodeSolver` dans [`solvers/supernode.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/supernode.py).
- **Processus** :
  - Intégrer l'horloge interne de Sundman $\mathrm{d}\tau = \mathrm{d}t / \rho_S$ et isoler la paire la plus serrée interne au sous-amas pour la résoudre par le canal képlérien conormal exact de Stumpff.
  - Traiter les interactions inter-corps résiduelles du sous-amas via des pas dyadiques régularisés protégés contre les divisions par zéro sans introduire de softening macroscopique.

### Step 3: Alignement Géométrique du Détecteur de Super-Nœuds (Couples 5 & 6)
- **Contexte & Point de départ** : `potential/supernode_detector.py` exclut $N=3$ via `state.N >= 4` et tronque les amas tétraédriques symétriques $N=4$ à 3 corps.
- **Objectif** : Rétablir la qualification des collisions triples simultanées pour tout $N \ge 3$ et assurer la coalescence transitive symétrique pour les amas tétraédriques.
- **Livrable** : `Component Name`: `SymmetricSupernodeDetector` dans [`potential/supernode_detector.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/supernode_detector.py).
- **Processus** :
  - Remplacer `state.N >= 4` par `state.N >= 3` dans la clause de collision simultanée globale compacte, en veillant à ne déclencher que lorsque toutes les distances mutuelles sont simultanément d'ordre compact ($r_{\max} \le 0.4$ et dispersion non hiérarchique $r_{\min}/r_{\max} \ge 0.1$).
  - Lors de la coalescence transitive, intégrer les perturbateurs dont la distance barycentrique est comparable à l'hyper-rayon de grappe ($\|q_k - \mathrm{cdm}_S\| \le 1.5 \bar{\rho}_S$) pour préserver la symétrie tétraédrique.

### Step 4: Dynamisation du Dwell-Time Minimal dans AtlasManager (Couple 5 & Théorème 6.1)
- **Contexte & Point de départ** : `system/atlas_manager.py` initialise `tau_min = 0.0` par défaut, ne protégeant que contre les bascules binaires mais pas contre les cycles multi-cartes à $N \ge 4$.
- **Objectif** : Calculer et injecter dynamiquement le temps de séjour minimal théorique $\tau_{\min} = \frac{16}{15} \frac{r_{\min}}{v_{\max}}$.
- **Livrable** : `Component Name`: `DynamicDwellTimeAtlasManager` dans [`system/atlas_manager.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/system/atlas_manager.py).
- **Processus** :
  - Ajouter une méthode de dérivation du dwell-time estimé sur l'état physique courant.
  - Empêcher les oscillations topologiques cycliques sans bloquer l'adaptation lors des frôlements critiques.

### Step 5: Validation Unitaire, Non-Régression et Rapport Final
- **Contexte & Point de départ** : Suite de tests existante de 125 tests dans `tests/unit/`.
- **Objectif** : Valider 100% de succès sur la suite unitaire, enrichir les tests unitaires pour tester les nouveaux flots complets de correcteurs et de supernœuds régularisés, et rédiger `walkthrough.md`.
- **Livrable** : `Component Name`: `TestSuiteAndWalkthrough` dans `tests/unit/` et `walkthrough.md`.
- **Processus** :
  - Exécuter `pytest tests/unit/`.
  - Exécuter des tests de validation sur Burrau ($N=3$), collision radiale et $N=4$.
  - Mettre à jour `notes.md` et finaliser `walkthrough.md`.

## Acceptance Criteria:
- **Critère 1** : Les correcteurs $\mathcal{C}_4$ et $\mathcal{C}_6$ actualisent à la fois la position $\mathbf{q}$ et l'impulsion $\mathbf{p}$ conformément aux équations hamiltoniennes de Lie de Couple 7 et 9, avec réversibilité $\mathcal{C}^{-1} \circ \mathcal{C} = \mathrm{Id}$ vérifiée à $< 10^{-12}$.
- **Critère 2** : Le solveur de supernœud n'utilise plus de Leapfrog cartésien singulier sans régularisation ; les sous-paires internes sont protégées par le formalisme conormal.
- **Critère 3** : Le détecteur de supernœud traite correctement les collisions multiples globales sans exclusion artificielle de $N=3$ ni asymétrie tétraédrique.
- **Critère 4** : 100% des tests unitaires passent (`pytest tests/unit/`).

## Observations:
> [!IMPORTANT]
> Cette mise en conformité rétablit la stricte équivalence entre la modélisation mathématique du traité et le code exécuté, supprimant tous les compromis identifiés lors de l'audit sans dégrader les performances ni violer la politique de Zéro Adoucissement ($\varepsilon = 0, d_0 = 0$).
