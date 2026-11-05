# Implementation Plan of Iteration 160 : Phase 1 — Remédiation et Optimisation du Moteur de Calcul Python (Codebase & Algorithmes)
Created at: 2026-10-09T08:28:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur a validé l'arbitrage pré-itératif pour la Phase 1 en sélectionnant :
   - `1 : A` (Fusion FSAL analytique des demi-kicks dans `Yoshida4Integrator` pour réduire le nombre d'évaluations de force de 6 à 4 par pas).
   - `2 : A` (Formalisation explicite de l'opérateur de super-nœud synchrone `step_supernode_flow` avec alias `step_supernode_kick` dans `integrators/common.py`).
   - `3 : A` (Harmonisation du signe de la matrice unimodulaire $\mathbf{M}_{\mathrm{coord}} \in \mathrm{SL}(2, \mathbb{R})$ dans `core/transitions.py` sur l'Équation 4.116 du manuscrit de recherche).
2. **Implicit Strategic Value**: Éliminer immédiatement les 6 discrépances logicielles et algorithmiques identifiées par l'audit contradictoire dans le code Python avant toute modification du manuscrit LaTeX ou du noyau Lean 4. Optimiser l'efficacité de calcul de Yoshida 4 de 33%, clarifier les interfaces de splitting pour les super-nœuds, assurer une conformité géométrique absolue sur les matrices de transition de Jacobi, et documenter sans équivoque l'annulation analytique des monopôles d'apex.
3. **Alignment Criteria**:
   - `Yoshida4Integrator.step()` exécute exactement 4 évaluations de force par pas avec la chaîne d'impulsions fusionnées $(0.5 w_1, 0.5(w_1+w_0), 0.5(w_0+w_1), 0.5 w_1)$, confirmant la convergence d'ordre 4.
   - `step_supernode_flow` (et son alias `step_supernode_kick`) est formellement défini et exporté dans [`common.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/common.py), documentant le flot hamiltonien couplé sous temps de Sundman.
   - Les coefficients de `SL2TransitionMatrix.from_masses` dans [`transitions.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/core/transitions.py) sont strictement conformes à l'Équation (4.116) avec $\det(\mathbf{M}) = +1$ et relèvement symplectique $\operatorname{Sp}(12, \mathbb{R})$ invariant.
   - La démonstration de l'annulation analytique des monopôles écrantés $V_{\mathrm{kick}} = V_{\mathrm{apex}}^{\mathrm{hier}} + V_{\mathrm{kepler}}^{\mathrm{screened}} = V_{\mathrm{inter}}$ est formellement documentée dans [`screened_apex.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screened_apex.py).
   - Le commentaire de repli de Halley dans [`kepler_sundman.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/kepler_sundman.py) décrit fidèlement la bissection déterministe à 60 étapes.
   - 100 % de succès sur l'ensemble de la suite de 150 tests unitaires.

## Guidelines:
- **Physique & Géométrie** : Zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$). Invariance absolue sous translations $\mathbb{R}^3$, rotations $\mathrm{SO}(3)$ et boosts galiléens.
- **Symplecticité** : Conservation de la 2-forme de Darboux à la précision machine float64 IEEE 754 ($< 10^{-14}$).
- **Édition de Code Source** : Interdiction absolue de `write_to_file` avec `overwrite: true` sur tout fichier source existant. Utilisation exclusive de `replace_file_content` ciblée sur les blocs sémantiques minimaux.
- **Non-Régression** : Tous les tests unitaires et benchmarks existants doivent demeurer strictement validés.

## Plan:

### Step 1: Implémentation de la Fusion FSAL dans `Yoshida4Integrator`
- **Context**: [`integrators/yoshida.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/yoshida.py) boucle actuellement sur `YOSHIDA_WEIGHTS` en appelant `_strang.step()`, ce qui effectue $2 \times 3 = 6$ évaluations de force au lieu des 4 évaluations permises par la fusion des demi-kicks adjacents.
- **Objective**: Définir les poids d'impulsion fusionnés FSAL et réécrire la méthode `step()` pour enchaîner 4 kicks et 3 dérives, divisant par 1.5 le nombre de calculs de force.
- **Deliverables**:
  - `FSALYoshidaWeights`: Poids d'impulsion explicites `YOSHIDA_IMPULSE_WEIGHTS = (0.5 * YOSHIDA_W1, 0.5 * (YOSHIDA_W1 + YOSHIDA_W0), 0.5 * (YOSHIDA_W0 + YOSHIDA_W1), 0.5 * YOSHIDA_W1)`.
  - `OptimizedYoshidaStep`: Méthode `step()` restructurée avec la séquence `Kick -> Drift -> Kick -> Drift -> Kick -> Drift -> Kick`.
- **Process**: Modifier [`yoshida.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/yoshida.py) par `replace_file_content`.

### Step 2: Formalisation de l'Opérateur de Super-Nœud dans `integrators/common.py`
- **Context**: Le super-nœud synchrone est résolu par `step_supernode_synchronous` directement dans `step_drift_and_kepler`, sans fonction dédiée exposant l'interface attendue par les spécifications de kick/flot.
- **Objective**: Exposer formellement la fonction `step_supernode_flow` (avec alias `step_supernode_kick`) dans [`common.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/common.py) et documenter la nature du flot hamiltonien étendu sans duplication barycentrique.
- **Deliverables**:
  - `StepSupernodeFlowFunction`: Fonction `step_supernode_flow(state, tree, dt, ...)` et alias `step_supernode_kick`.
  - `UpdatedDriftAndKepler`: Branchement explicite dans `step_drift_and_kepler`.
- **Process**: Modifier [`common.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/common.py) par `replace_file_content`.

### Step 3: Harmonisation de l'Orientation de Jacobi dans `core/transitions.py`
- **Context**: Dans [`core/transitions.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/core/transitions.py), la matrice $\mathbf{M}_{\mathrm{coord}}$ porte un signe global $-1$ par rapport à l'Équation (4.116) du papier.
- **Objective**: Aligner les coefficients $(a, b, c, d)$ sur l'Équation (4.116) en définissant $(-m_A/m_{AB}, 1, -m_B m_{ABC}/(m_{AB} m_{BC}), -m_C/m_{BC})$, vérifier $\det(\mathbf{M}) = +1$ et adapter les tests unitaires associés dans [`test_transitions.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/core/test_transitions.py).
- **Deliverables**:
  - `AlignedSL2TransitionMatrix`: Coefficients conformes à l'Éq. 4.116.
  - `SynchronizedTransitionTests`: Tests unitaires vérifiant $\det=+1$, la composition unimodulaire et le relèvement symplectique $\operatorname{Sp}(12, \mathbb{R})$.
- **Process**: Modifier [`transitions.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/core/transitions.py) et [`test_transitions.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/core/test_transitions.py) par `replace_file_content`.

### Step 4: Documentation Analytique et Sécurisation des Gradients (`screened_apex.py` & `kepler_sundman.py`)
- **Context**: L'audit a pointé le besoin d'expliciter l'annulation analytique des monopôles d'apex dans le kick, de préciser la portée locale de `compute_screened_composite_potential`, et de rectifier la mention "Brent/bisection" dans le solveur de Sundman.
- **Objective**:
  1. Ajouter dans [`screened_apex.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screened_apex.py) et [`common.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/common.py) l'explication mathématique de $V_{\mathrm{kick}} = V_{\mathrm{inter}}$ due à l'annulation mutuelle des monopôles écrantés.
  2. Documenter la nature de gradient partiel de `compute_screened_composite_potential`.
  3. Mettre à jour le commentaire dans [`kepler_sundman.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/kepler_sundman.py) pour indiquer la bissection déterministe pure à 60 étapes.
- **Deliverables**:
  - `ScreenedApexDocumentation`: Docstrings enrichies et rigoureuses.
  - `KeplerSundmanCommentUpdate`: Commentaire technique exact.
- **Process**: Modifier [`screened_apex.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screened_apex.py) et [`kepler_sundman.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/kepler_sundman.py) par `replace_file_content`.

### Step 5: Validation Unitaire Intégrale, Non-Régression et Clôture Git
- **Context**: S'assurer que l'intégralité des 150 tests unitaires passe sans avertissement et sceller les modifications de la Phase 1 par un commit Git.
- **Objective**: Exécuter `pytest tests/unit/`, vérifier la convergence de Yoshida 4, et enregistrer le commit Git propre.
- **Deliverables**:
  - `PytestValidation`: 150/150 tests unitaires passants.
  - `GitCommitPhase1`: Commit Git référençant la Phase 1 de la Roadmap VIII.
- **Process**: Exécuter `pytest` et `git commit` via `run_command`.

## Acceptance Criteria:
### Step 1:
- `Yoshida4Integrator.step()` effectue exactement 4 évaluations de `compute_apex_forces` par pas (au lieu de 6).
- La convergence d'ordre 4 de Yoshida 4 est certifiée avec un taux d'erreur rigoureusement identique ou supérieur.
### Step 2:
- `step_supernode_flow` et `step_supernode_kick` sont importables depuis `normal_cone_nbody_float64.integrators.common`.
- L'intégration des super-nœuds préserve la quantité de mouvement totale à $< 10^{-15}$.
### Step 3:
- `SL2TransitionMatrix.from_masses` génère des matrices avec $b = +1.0$ (conforme à l'Éq. 4.116).
- $\det(\mathbf{M}) = +1.0$ à l'erreur machine ($< 10^{-15}$).
- 100% des tests de `tests/unit/core/test_transitions.py` passent avec succès.
### Step 4:
- Les docstrings de `screened_apex.py` détaillent formellement l'annulation analytique des monopôles.
- Le commentaire de la ligne 124 de `kepler_sundman.py` mentionne la bissection pure déterministe à 60 étapes.
### Step 5:
- La commande `python -m pytest tests/unit/` valide 150/150 tests sans échec.
- Le commit Git scellant la Phase 1 est créé avec succès.

## Observations:
> [!IMPORTANT]
> La modification du signe de $\mathbf{M}_{\mathrm{coord}}$ dans `transitions.py` correspond à la permutation de l'orientation du vecteur relatif $\mathbf{r}_1' = \mathbf{q}_3 - \mathbf{q}_2$ (conforme à la convention de Delves (1960) adoptée dans le Chapitre 4). Puisque la transformation de phase cotangente relève $\mathbf{M}$ et $(\mathbf{M}^{-1})^T$, cette inversion globale conserve rigoureusement la symplecticité $\operatorname{Sp}(12, \mathbb{R})$.
