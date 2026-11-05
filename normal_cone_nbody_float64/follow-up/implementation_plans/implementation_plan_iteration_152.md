# Implementation Plan of Iteration 152 : Complétion de Leibniz sur le Gradient de $\mathcal{C}_6$ et Télescopage d'Ombre Inter-Pas
Created at: 2026-10-08T12:25:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur poursuit la Roadmap VII vers la Phase 3 (Itération 152) pour éradiquer les Divergences Critiques 3 et 4 identifiées par l'audit contradictoire à 9 sous-agents : compléter rigoureusement le gradient d'impulsion $\nabla_p K_6$ dans `normal_cone_nbody_float64/integrators/correctors.py` avec le terme de Leibniz direct $+ \frac{2}{m_j^2} (\mathbf{F}_j \cdot \mathbf{v}_j) \mathbf{F}_j$, brancher l'opérateur de connexion de jauge de Lie `apply_gauge_connection` lors des transitions d'atlas dans `system/nbody_system.py`, et optimiser la chaîne télescopique $\mathcal{C} \circ \mathcal{S}^N \circ \mathcal{C}^{-1}$ pour supprimer les évaluations redondantes de correcteur.
2. **Implicit Strategic Value**: Les intégrateurs traités (Processed Symplectic Integrators, PBM) reposent sur la conjugaison hamiltonienne $\Psi(h) = \mathcal{C}(h) \circ \mathcal{S}(h) \circ \mathcal{C}^{-1}(h)$. Dans une trajectoire à plusieurs pas, la propriété de télescopage $\mathcal{C}^{-1} \circ \mathcal{C} = \mathrm{Id}$ doit annuler toutes les évaluations de correcteurs internes, réduisant le coût de $2M$ évaluations à seulement 2 évaluations aux extrémités. Cependant, si l'atlas de Jacobi bascule de $\mathcal{T}_A$ à $\mathcal{T}_B$, le non-commutateur $[\mathcal{X}_{K_A}, \mathcal{X}_{K_B}] \ne 0$ brise le télescopage et injecte un saut interfacial dans l'hamiltonien d'ombre ; l'application de la connexion de jauge $\mathcal{J}_{A \to B} = \mathcal{C}_B^{-1} \circ \mathcal{C}_A$ rétablit la continuité symplectique. De surcroît, l'exactitude de $\nabla_p K_6$ garantit que le flux de Verlet pour $\mathcal{C}_6$ dérive rigoureusement du générateur scalaire $K_6$.
3. **Alignment Criteria**:
   - Sanctuarisation absolue de `redaction/unsoftened_conormal_nbody3d/` (lecture seule).
   - Zéro adoucissement résiduel ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$).
   - Exactitude analytique de $\nabla_p K_6$ et de la force conjuguée vérifiée contre différences finies centrales d'ordre 4 ($< 10^{-10}$).
   - Préservation stricte de l'inversion algébrique exacte $\mathcal{C}_6^{-1} \circ \mathcal{C}_6 = \mathrm{Id}$ à la précision machine float64 ($< 10^{-14}$).
   - Préservation de tous les 10 invariants de Lie-Bargmann à travers la connexion de jauge $\mathcal{J}_{A \to B}$ ($< 10^{-14}$).
   - 100% des tests unitaires validés (au moins 146/146 tests existants et nouveaux tests passés).

## Guidelines:
- **Méthode** :
  - Dans `integrators/correctors.py`, implémenter le gradient canonique complet de $K_6(\mathbf{q}, \mathbf{p}) = \kappa_6 \sum_i \frac{1}{m_i} (\mathbf{F}_i \cdot \mathbf{v}_i)^2$ :
    $$\dot{\mathbf{q}}_j = \nabla_{\mathbf{p}_j} K_6 = \kappa_6 \frac{2}{m_j^2} (\mathbf{F}_j \cdot \mathbf{v}_j) \mathbf{F}_j$$
    $$\dot{\mathbf{p}}_j = -\nabla_{\mathbf{q}_j} K_6 = -\kappa_6 D_{\mathbf{u}} \mathbf{F}_j \quad \text{avec} \quad \mathbf{u}_i = \frac{2}{m_i} (\mathbf{F}_i \cdot \mathbf{v}_i) \mathbf{v}_i.$$
  - Dans `system/nbody_system.py`, interconnecter `apply_gauge_connection` lorsque `was_flipped` est True lors d'une transition d'arbre, préservant la continuité de l'hamiltonien d'ombre.
  - Dans `system/nbody_system.py`, concevoir le mode d'intégration d'ombre pour `evolve()` permettant le télescopage macro-pas $\mathcal{C} \circ \mathcal{S}^M \circ \mathcal{C}^{-1}$ tout en assurant l'exactitude de l'enregistrement de l'historique physique.
- **Lignes Rouges** :
  - Aucun adoucissement ($\varepsilon = 0$).
  - Aucun bruit de différences finies résiduel dans le gradient de $K_6$.
  - Aucune dégradation de la réversibilité temporelle microscopique ($< 10^{-14}$).
- **Architecture** :
  - Maintenir les signatures de `SymplecticCorrector`, `CorrectorC6`, `ProcessedSymplecticIntegrator` et `NBodySystem`.

## Plan:
### Step 1: Formulation Mathématique et Validation Finie du Gradient Canonique de $K_6$
- **Context**: L'audit contradictoire a révélé que `CorrectorC6.apply()` utilisait une formulation advective heuristique au lieu du gradient variationnel exact de $K_6 = \kappa_6 \sum_i \frac{1}{m_i} (\mathbf{F}_i \cdot \mathbf{v}_i)^2$.
- **Objective**: Écrire un script de validation dans `scratch/verify_k6_gradient.py` confrontant les dérivées analytiques $\nabla_p K_6$ et $-\nabla_q K_6$ aux différences finies centrales d'ordre 4.
- **Deliverables**:
  - `K6GradientVerificationScript`: Script dans `scratch/verify_k6_gradient.py`.
- **Process**: Calculer $K_6(q, p)$ scalaire et comparer $\nabla_p K_6$ et $-\nabla_q K_6$ contre $\frac{-f(x+2\eta) + 8f(x+\eta) - 8f(x-\eta) + f(x-2\eta)}{12\eta}$.

### Step 2: Implémentation du Gradient de Leibniz Canonique dans `CorrectorC6`
- **Context**: `CorrectorC6.apply()` doit exécuter le flot de Verlet pour le champ hamiltonien $\mathcal{X}_{K_6}$.
- **Objective**: Remplacer l'échelle advective dans `CorrectorC6.apply()` par le terme direct de Leibniz $\dot{\mathbf{q}}_j = \kappa_6 \frac{2}{m_j^2} (\mathbf{F}_j \cdot \mathbf{v}_j) \mathbf{F}_j$ et la dérivée directionnelle $\dot{\mathbf{p}}_j = -\kappa_6 D_{\mathbf{u}} \mathbf{F}_j$ avec $\mathbf{u}_i = \frac{2}{m_i} (\mathbf{F}_i \cdot \mathbf{v}_i) \mathbf{v}_i$.
- **Deliverables**:
  - `CanonicalCorrectorC6`: Méthode `apply()` mise à jour dans `normal_cone_nbody_float64/integrators/correctors.py`.
- **Process**: Mettre à jour `CorrectorC6.apply()` et vérifier l'inversion exacte $\mathcal{C}_6^{-1} \circ \mathcal{C}_6 = \mathrm{Id}$ ($< 10^{-14}$).

### Step 3: Branchement de la Connexion de Jauge Inter-Atlas dans `NBodySystem`
- **Context**: Lors d'un basculement de carte de Jacobi ($\mathcal{T}_A \to \mathcal{T}_B$), le champ de correction de $\mathcal{T}_A$ ne commute pas avec celui de $\mathcal{T}_B$. L'opérateur de jauge `apply_gauge_connection` ($\mathcal{J}_{A \to B} = \mathcal{C}_B^{-1} \circ \mathcal{C}_A$) doit être appliqué.
- **Objective**: Dans `NBodySystem.step()`, lorsque `was_flipped` est True et que l'intégrateur est traité (ou sous télescopage), appliquer `apply_gauge_connection` pour transférer continûment l'état sur l'hamiltonien d'ombre de la nouvelle carte.
- **Deliverables**:
  - `AtlasGaugeConnectionBridge`: Branchement dans `normal_cone_nbody_float64/system/nbody_system.py`.
- **Process**: Intégrer l'appel à `apply_gauge_connection` lors des transitions détectées par `evaluate_tree_adaptation`.

### Step 4: Optimisation du Télescopage d'Ombre Inter-Pas dans `NBodySystem`
- **Context**: L'application de $\mathcal{C}^{-1}$ à l'entrée et $\mathcal{C}$ à la sortie de chaque macro-pas coûte $2M$ évaluations de correction, détruisant l'avantage du télescopage théorique où les correcteurs intermédiaires s'annulent ($\mathcal{C} \circ \mathcal{C}^{-1} = \mathrm{Id}$).
- **Objective**: Permettre à `NBodySystem.evolve()` d'exécuter la chaîne télescopée globale $\mathcal{C} \circ (\mathcal{S}_{\mathrm{macro}})^M \circ \mathcal{C}^{-1}$ tout en préservant l'enregistrement régulier de l'historique physique sans multiplier les évaluations superflues.
- **Deliverables**:
  - `OptimizedShadowTelescoping`: Méthode d'évolution optimisée dans `normal_cone_nbody_float64/system/nbody_system.py`.
- **Process**: Structurer l'état d'ombre persistant ou le télescopage dans `evolve()`.

### Step 5: Tests Unitaires et Certification de Non-Régression
- **Context**: Les tests unitaires dans `tests/unit/integrators/test_correctors.py` et sur l'ensemble de la suite doivent confirmer l'exactitude de $\mathcal{C}_6$ et la stabilité globale.
- **Objective**: Ajouter des tests unitaires vérifiant $\nabla_p K_6$, la connexion de jauge inter-atlas, et exécuter la suite complète de 146+ tests.
- **Deliverables**:
  - `TestCorrectorsLeibniz`: Nouveaux tests unitaires dans `tests/unit/integrators/test_correctors.py`.
  - `FullTestSuitePass`: Validation complète via `pytest tests/unit/`.
- **Process**: Écrire et lancer les tests unitaires.

### Step 6: Git Commit de Phase 3
- **Context**: Phase 3 de la Roadmap VII achevée.
- **Objective**: Enregistrer un commit atomique documenté.
- **Deliverables**:
  - `GitCommit`: Commit git avec message standardisé.
- **Process**: `git add` et `git commit`.

## Acceptance Criteria:
### Step 1:
- Les composantes de $\nabla_p K_6$ et $-\nabla_q K_6$ concordent avec les différences finies centrales d'ordre 4 avec une erreur relative maximale $< 10^{-9}$.
### Step 2:
- L'opérateur $\mathcal{C}_6$ conserve la propriété d'inversion algébrique exacte $\mathcal{C}_6^{-1} \circ \mathcal{C}_6 = \mathrm{Id}$ avec une erreur sur $q$ et $p$ inférieure à $10^{-14}$.
- L'impulsion totale $\sum \mathbf{p}_i$ est conservée à $< 10^{-15}$ lors de l'application de $\mathcal{C}_6$.
### Step 3:
- Lors d'une transition d'arbre provoquée, `apply_gauge_connection` préserve les 10 invariants de Lie-Bargmann à $< 10^{-14}$.
### Step 4:
- Le nombre d'évaluations de correcteurs dans une simulation de $M$ macro-pas sans transition d'arbre est réduit de $2M$ à 2.
### Step 5:
- 100% des tests unitaires (au moins 148+ tests) passent avec succès sans avertissement.
### Step 6:
- Commit git propre enregistré sur `main`.

## Observations:
> [!IMPORTANT]
> Le potentiel $K_6 = \kappa_6 \sum_i \frac{1}{m_i} (\mathbf{F}_i \cdot \mathbf{v}_i)^2$ dépend quadratiquement des vitesses relatives par le biais de l'impulsion. Sa dérivée par rapport à l'impulsion $\nabla_p K_6$ comporte à la fois un terme proportionnel à la force d'apex locale $\mathbf{F}_i$ (terme direct de Leibniz) et un terme projeté sur la hessienne directionnelle $D_v \mathbf{F}$. L'inclusion du terme direct est indispensable pour garantir que le flot de Verlet intègre rigoureusement le champ hamiltonien $\mathcal{X}_{K_6}$.
