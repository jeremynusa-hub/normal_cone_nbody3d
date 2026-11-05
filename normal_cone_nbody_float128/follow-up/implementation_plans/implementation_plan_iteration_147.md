# Implementation Plan of Iteration 147 : Phase 4 — Remédiation de la Jauge de Lie et des Correcteurs Symplectiques Traités
Created at: 2026-10-08T10:50:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur a confirmé la continuation de la feuille de route (« Continue la roadmap donc » / « go »). Conformément à la Roadmap VI, cela active la **Phase 4 (Itération 147)** consacrée à la remédiation de la géométrie de jauge de Lie et des correcteurs symplectiques traités $\mathcal{C}_4$ et $\mathcal{C}_6$ pour PBM4 et PBM6.
2. **Implicit Strategic Value**: Éradiquer les deux dernières sources d'imperfection théorique et de bruit numérique dans la chaîne symplectique : le plancher de bruit à $10^{-9}$ induit par le stencil par différences finies ($\eta = 10^{-7}$) dans l'évaluation de $D_v \mathbf{F}$, les constantes homologiques génériques de Leapfrog ($1/2880, 1/120960$) inadaptées aux intégrateurs de Blanes--Moan, et la rupture de jauge de Lie $[X_{K_A}, X_{K_B}] \neq 0$ lors des transitions d'atlas qui injectait un saut d'énergie interfacial $\Delta H = h^{2p}(H_{2p+1}^{(B)} - H_{2p+1}^{(A)})$.
3. **Alignment Criteria**: Produit hessien-vecteur analytique $D_v \mathbf{F}$ calculé à la précision machine float64 IEEE 754 ($< 10^{-16}$) sans paramètre $\eta$ ; constantes homologiques $\kappa_4 \approx 0.03678699950623901$ et $\kappa_6$ calibrées analytiquement sur les quadratures BCH du manuscrit ; opérateur de connexion de jauge $\mathcal{J}_{A \to B} = \mathcal{C}_B^{-1} \circ \mathcal{C}_A$ implémenté et certifié ; 100 % de succès sur l'ensemble de la suite de tests unitaires.

---

## Guidelines:
- **Method**: Respect strict des équations homologiques $\{H, K_{2p}\} = H_{2p+1}$ du Chapitre 5 (Section 5.5, Théorèmes 5.16 et 5.17) et des théorèmes certifiés en Lean 4 (`01_CorrectorBM4_Theorems.lean` et `02_CorrectorBM6_Theorems.lean`).
- **Red Lines**: Zéro adoucissement ni garde-fou artificiel ($\varepsilon = 0, d_0 = 0$). Interdiction de réintroduire des différences finies ou des approximations numériques dans les dérivées de force. Sanctuarisation absolue du répertoire du manuscrit `redaction/unsoftened_conormal_nbody3d/`.
- **Architecture**: Maintien de l'isolation modulaire : `correctors.py` encapsule les opérateurs de correction et leurs dérivées analytiques, `processed.py` gère le décorateur PBM et la connexion de jauge, et `nbody_system.py` orchestre les transitions aux frontières macro.

---

## Plan:

### Step 1: Dérivation et Implémentation du Produit Hessien-Vecteur Analytique dans `correctors.py`
- **Context**: Actuellement, `_directional_force_derivative` dans `correctors.py` évalue $D_v \mathbf{F}_i$ via un stencil centré par différences finies avec un pas $\eta = 10^{-7}$. En float64, cette approximation perd la moitié des chiffres significatifs et injecte un plancher d'erreur à $\approx 10^{-9}$, brisant l'ordre 8 asymptotique de PBM6.
- **Objective**: Remplacer le stencil par différences finies par l'évaluation analytique exacte du produit hessien-vecteur pour les forces newtoniennes d'apex :
  $$D_v \mathbf{F}_{ij} = G m_i m_j \left[ \frac{\Delta \mathbf{v}_{ij}}{r_{ij}^3} - 3 \frac{(\mathbf{r}_{ij} \cdot \Delta \mathbf{v}_{ij})}{r_{ij}^5} \mathbf{r}_{ij} \right], \quad D_v \mathbf{F}_i = \sum_{j \neq i, (i,j) \notin \text{masked}} D_v \mathbf{F}_{ij}$$
  où $\mathbf{r}_{ij} = \mathbf{q}_j - \mathbf{q}_i$ et $\Delta \mathbf{v}_{ij} = \mathbf{v}_j - \mathbf{v}_i$.
- **Deliverables**:
  - `integrators/correctors.py`: Fonction vectorisée `compute_analytical_directional_force_derivative(eval_apex, state, tree, supernode_bodies, v_dir)` calculant analytiquement le gradient directionnel avec masquage strict des paires de feuilles et des sous-amas de super-nœuds, remplaçant l'ancienne fonction à différences finies.
- **Process**:
  1. Dériver et coder la formule tensorielle vectorisée de la dérivée directionnelle sans division par zéro ($r_{ij} < 10^{-15}$).
  2. Masquer les paires de feuilles de `tree.leaf_pairs` et les corps appartenant à `supernode_bodies`.
  3. Vérifier que $\sum_i D_v \mathbf{F}_i = \mathbf{0}$ à la précision machine float64 ($< 10^{-16}$).
  4. Mettre à jour `CorrectorC4.apply` et `CorrectorC6.apply` pour invoquer exclusivement cette dérivée analytique.

### Step 2: Calibrage Analytique des Constantes Homologiques $\kappa_4$ et $\kappa_6$ dans `correctors.py`
- **Context**: Actuellement, `CorrectorC4` utilise $\kappa_4 = 1/2880 \approx 0.000347$ et `CorrectorC6` utilise $\kappa_6 = 1/120960 \approx 8.26 \times 10^{-6}$, qui sont les constantes du schéma Leapfrog à 2 étages. Or, la composition Blanes--Moan BM4 à 7 étages présente un résidu homologique distinct donné par l'Équation (5.51) du manuscrit :
  $$\kappa_4 = \sum_{j=1}^7 b_j \left( \sum_{k=1}^j a_k - \frac{1}{2} \right)^3 - \frac{1}{24} \sum_{j=1}^7 b_j a_j^2 \approx 0.03678699950623901$$
- **Objective**: Définir les constantes canoniques $\kappa_4$ et $\kappa_6$ adaptées à BM4 et BM6 et les instancier par défaut dans les correcteurs.
- **Deliverables**:
  - `core/constants.py`: Déclaration des constantes `KAPPA4_BM4 = 0.03678699950623901` et `KAPPA6_BM6 = 0.0004128525`.
  - `integrators/correctors.py`: Mise à jour des valeurs par défaut de `kappa4` et `kappa6` dans `CorrectorC4.__init__` et `CorrectorC6.__init__`.
- **Process**:
  1. Ajouter `KAPPA4_BM4` et `KAPPA6_BM6` dans `core/constants.py`.
  2. Importer et assigner ces valeurs par défaut dans `CorrectorC4` et `CorrectorC6`.
  3. Mettre à jour la documentation des classes avec les références du manuscrit.

### Step 3: Opérateur de Connexion de Jauge Symplectique $\mathcal{J}_{A \to B}$ dans `processed.py` et `system/nbody_system.py`
- **Context**: Lors d'un changement d'arbre de Jacobi $T_A \to T_B$ à la frontière de macro-pas, l'hamiltonien d'ombre subit un saut interfacial $\Delta H = h^{2p}(H_{2p+1}^{(B)} - H_{2p+1}^{(A)})$ dû à la non-commutativité $[X_{K_A}, X_{K_B}] \neq 0$.
- **Objective**: Implémenter formellement l'opérateur de connexion de jauge symplectique $\mathcal{J}_{A \to B} \coloneqq \mathcal{C}_B^{-1} \circ \mathcal{C}_A$ réconciliant les coordonnées d'ombre de deux cartes d'atlas différentes, et s'assurer que la séquence de macro-pas dans `NBodySystem` opère cette transformation sans discontinuité.
- **Deliverables**:
  - `integrators/processed.py`: Fonction canonique `apply_gauge_connection(state, old_tree, new_tree, dt, corrector, apex, supernode_bodies) -> ExtendedPhaseState`.
  - `system/nbody_system.py`: Intégration explicite de la connexion de jauge lors des transitions d'atlas à la frontière macro pour `ProcessedSymplecticIntegrator`.
- **Process**:
  1. Implémenter `apply_gauge_connection` dans `processed.py` avec court-circuit identité si `old_tree is new_tree`.
  2. Vérifier la compatibilité avec la boucle de macro-pas de `NBodySystem.step`.
  3. Valider que la connexion préserve la 2-forme de Darboux et les 10 invariants de Lie--Bargmann à la précision machine.

### Step 4: Tests Unitaires Formels et Validation de Non-Régression
- **Context**: La suite de tests unitaires (135 tests) doit valider la dérivée analytique, l'inversion machine, le calibrage homologique et la connexion de jauge.
- **Objective**: Attester de l'éradication du bruit numérique et du fonctionnement parfait des correcteurs.
- **Deliverables**:
  - `tests/unit/integrators/test_correctors.py`: Nouveaux tests pour le produit hessien-vecteur analytique (vérifiant l'accord avec les différences finies à $< 10^{-11}$ et l'annulation stricte de la force totale), et pour les constantes $\kappa_4, \kappa_6$.
  - `tests/unit/integrators/test_processed.py`: Test unitaire de l'opérateur de connexion de jauge $\mathcal{J}_{A \to B}$ sur une transition d'atlas réelle à 3 et 4 corps.
- **Process**:
  1. Enrichir les modules de tests unitaires existants.
  2. Exécuter `pytest tests/unit/` et s'assurer du taux de succès 100 %.

### Step 5: Validation Git Commit
- **Context**: Clôture de la Phase 4 de la Roadmap VI.
- **Objective**: Consigner les remédiations de la jauge de Lie et des correcteurs dans un commit git traçable.
- **Deliverables**:
  - Commit Git structuré : `feat(integrators,system): remediate Lie gauge connection, analytical Hessian-vector product, and homological constants - Phase 4 of Roadmap VI (Iteration 147)`.
- **Process**:
  1. Vérifier l'état git via `git status`.
  2. Indexer les fichiers modifiés et valider le commit.

---

## Acceptance Criteria:

### Step 1:
- `compute_analytical_directional_force_derivative` évalue exactement la dérivée directionnelle sans utiliser de pas $\eta$ fini.
- L'invariance par translation $\sum_i D_v \mathbf{F}_i = \mathbf{0}$ est respectée à $< 10^{-15}$.
- L'accord entre le calcul analytique et les différences finies centrées ($\eta \in [10^{-5}, 10^{-6}]$) est vérifié à $< 10^{-10}$.

### Step 2:
- `KAPPA4_BM4` vaut exactement `0.03678699950623901` et est la valeur par défaut de `CorrectorC4`.
- `KAPPA6_BM6` est défini et assigné par défaut dans `CorrectorC6`.
- Les correcteurs préservent l'inversion algébrique exacte $\mathcal{C}^{-1} \circ \mathcal{C} = \operatorname{Id}$ à $< 10^{-14}$.

### Step 3:
- L'opérateur de connexion de jauge $\mathcal{J}_{A \to B} = \mathcal{C}_B^{-1} \circ \mathcal{C}_A$ est implémenté et satisfait $\mathcal{J}_{A \to A} = \operatorname{Id}$ à $< 10^{-15}$.
- Lors d'une transition d'arbre $T_A \to T_B$, l'énergie de l'état physique raccordé ne présente aucun saut non-physique.

### Step 4:
- 100 % de succès sur l'ensemble de la suite de tests unitaires `tests/unit/` (>= 138 tests).

### Step 5:
- Commit Git tracé et archivé.

---

## Observations:
> [!IMPORTANT]
> L'élimination du paramètre de différences finies $\eta = 10^{-7}$ au profit du produit hessien-vecteur analytique ferme l'un des compromis numériques les plus pénalisants du solveur, restaurant l'exactitude machine float64 dans l'évaluation des commutateurs de Lie.

> [!NOTE]
> Le calibrage de $\kappa_4 \approx 0.03679$ (au lieu de $1/2880 \approx 0.000347$, soit un facteur 106x !) permet enfin à PBM4 d'annuler rigoureusement le terme d'erreur $H_5$ du schéma Blanes--Moan, assurant une véritable pente de convergence d'ordre 6.
