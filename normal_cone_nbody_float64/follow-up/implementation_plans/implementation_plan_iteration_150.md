# Implementation Plan of Iteration 150 : Sécurisation Absolue du Solveur Kepler-Sundman (Amorce Asymptotique, Fallback Monotone Garanti, Suppression du Clamping z)
Created at: 2026-10-08T12:09:40+02:00

## Context:
1. **Explicit Request**: L'utilisateur a validé l'ouverture de la Roadmap VII suite à la présentation des conclusions de l'audit contradictoire à 9 sous-agents. La Phase 1 cible l'éradication immédiate des défaillances critiques identifiées dans `normal_cone_nbody_float64/solvers/kepler_sundman.py` (cycle limite d'ordre 3 d'Halley et stagnation silencieuse à $\tau = 10^{-15}$, saturation du pas sur plongeon hyperbolique, et clamping dur sur $z$).
2. **Implicit Strategic Value**: Le solveur de Kepler-Sundman constitue le noyau analytique 2-corps universel sur lequel reposent toutes les paires de feuilles et les compositions symplectiques d'ordre élevé. Une défaillance de convergence non détectée injecte une erreur de position de 45% à 98% tout en masquant l'erreur via un Wronskien unimodulaire trompeur. Sécuriser ce solveur avec une garantie absolue de convergence monotone (100% sur 1000 orbites aléatoires extrêmes) élimine la principale source potentielle de gel dynamique et de sauts d'énergie lors des pas rétrogrades.
3. **Alignment Criteria**:
   - Sanctuarisation absolue de `redaction/unsoftened_conormal_nbody3d/` (lecture seule).
   - Zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$).
   - Réversibilité temporelle microscopique $\mathcal{S}(-\Delta t) \circ \mathcal{S}(\Delta t) = \operatorname{Id}$ validée à $< 10^{-14}$ sur 100% de 1000 orbites aléatoires incluant $e \ge 0.999$, $\Delta t > P$ et $v_0 \ge 500$.
   - Aucun gel silencieux de particule : interdiction formelle de renvoyer $\tau = 10^{-15}$ si $|F(\tau) - \Delta t| > \mathrm{tol} \cdot \Delta t$.
   - 100% des tests unitaires validés (au moins 142/142 tests passés).

## Guidelines:
- **Méthode** :
  - Conserver la méthode cubique d'Halley comme accélérateur d'ordre 3 dans son bassin d'attraction régulier.
  - Implémenter une zone de confiance sur le pas d'Halley : rejeter le saut si le dénominateur $2(F')^2 - FF'' \le 0.1(F')^2$ ou si $|\delta\tau| > 0.5\tau$, et basculer sur un pas de Newton amorti.
  - Adjoindre un solveur de secours bracketé monotone garanti (méthode de Brent / Dekker ou bissection) exploitant la stricte positivité de la dérivée temporelle $F'(\tau) = r(\tau) > 0$.
- **Lignes Rouges** :
  - Ne jamais introduire de softening ni de seuil de distance artificiel.
  - Ne jamais clamper $z$ à une valeur arbitraire déconnectée de la fonction de Stumpff (utiliser les formes trigonométriques/hyperboliques normalisées stables).
- **Architecture** :
  - Isoler les fonctions d'amorce, d'itération et de fallback dans `kepler_sundman.py` avec typage strict et docstrings mathématiques.

## Plan:
### Step 1: Reproduction et Isolation des Régimes de Défaillance
- **Context**: Les sous-agents 2 et 3 ont documenté des défaillances sur orbites inward excentriques ($e \ge 0.9, \sigma_0 < 0$) et hyperboliques ($v_0 \ge 50$).
- **Objective**: Écrire un script de test reproduisant les 12 cas d'échec sur 1000 orbites et mesurant précisément le comportement du dénominateur d'Halley et de l'erreur d'énergie.
- **Deliverables**:
  - `ReproduceFailureScript`: Script de diagnostic reproduisant la défaillance à $\tau = 10^{-15}$ et le blocage hyperbolique.
- **Process**: Évaluer `solve_sundman_time_equation` sur les configurations problématiques identifiées par l'audit.

### Step 2: Amorce Asymptotique Hyperbolique et Suppression du Clamping z
- **Context**: Pour $v_0 \ge 50$ ($\alpha \ll 0$), l'amorce cubique $\tau_0 = (6\Delta t/GM)^{1/3}$ est surévaluée d'un facteur 20, provoquant un plafonnement du pas d'Halley à $2/v_0$ et une non-convergence après 25 itérations. De plus, $z$ est clampé à $-500\,000$.
- **Objective**: Formuler une amorce adaptée au régime dominé par l'énergie cinétique $\tau_0 \approx \Delta t / (r_0 + \frac{1}{2} v_0 \Delta t)$ et supprimer le clamping dur sur $z$ en utilisant les formes stables de $\sinh(v)/v$.
- **Deliverables**:
  - `HyperbolicInitialGuess`: Calcul de $\tau_0$ adapté au signe de $\alpha$ et au ratio $v_0^2 / (2GM/r_0)$.
  - `UnclampedStumpff`: Gestion sans clamping arbitraire dans `kepler_sundman.py`.
- **Process**: Modifier `_initial_guess_sundman` et l'évaluation de $z$ dans `kepler_sundman.py`.

### Step 3: Contrôle de Pas de Halley et Fallback Monotone Garanti (Brent / Bissection)
- **Context**: Lorsque le dénominateur d'Halley dégénère ou que $\tau$ oscille en cycle limite, la boucle s'arrête à 25 itérations et renvoie silencieusement $\tau = 10^{-15}$.
- **Objective**: Sécuriser chaque pas d'Halley par une région de confiance et, en cas de non-convergence après 12 itérations, basculer sur un encadrement monotone de Brent/bissection garanti sur $[0, \tau_{\max}]$ car $F(\tau)$ est strictement croissante ($F'(\tau) = r(\tau) > 0$).
- **Deliverables**:
  - `SafeHalleyIterator`: Boucle itérative avec région de confiance et rejet de pas dégénérés.
  - `MonotoneBracketFallback`: Résolveur garanti par encadrement monotone éliminant tout cycle limite.
- **Process**: Implémenter l'algorithme robuste dans `solve_sundman_time_equation`.

### Step 4: Seuil Invariant de Collision Radiale
- **Context**: La ligne 190 utilise `is_head_on = bool(L_sq < 1e-24 * (GM * rho0 + 1e-15))` qui dépend de $\rho_0(t)$, brisant la symétrie temporelle $\mathcal{S}(-\Delta t) \circ \mathcal{S}(\Delta t) = \operatorname{Id}$.
- **Objective**: Remplacer ce critère dépendant de la position instantanée par un critère invariant d'échelle adimensionnel $L^2 / (G^2 M^2 \mu^2) < \epsilon_{\mathrm{rad}}$.
- **Deliverables**:
  - `InvariantHeadOnThreshold`: Critère indépendant de $\rho_0$.
- **Process**: Modifier la condition de collision radiale dans `step_kepler_sundman`.

### Step 5: Banc d'Essai Massif et Non-Régression
- **Context**: La suite existante compte 142 tests unitaires.
- **Objective**: Valider 1000 orbites aléatoires sans le moindre échec de réversibilité ($\Delta \le 10^{-14}$) et vérifier la non-régression intégrale.
- **Deliverables**:
  - `MassiveKeplerTest`: Test unitaire enrichi validant 1000 orbites extrêmes.
- **Process**: Exécuter `pytest` sur l'ensemble de la suite de tests unitaires.

### Step 6: Git Commit de Phase 1
- **Context**: Phase 1 de la Roadmap VII complétée.
- **Objective**: Enregistrer un commit atomique documenté.
- **Deliverables**:
  - `GitCommit`: Commit git avec message standardisé.
- **Process**: `git add` et `git commit`.

## Acceptance Criteria:
### Step 1:
- Les cas d'échec d'Halley identifiés par l'audit sont reproduits et mesurés avec précision.
### Step 2:
- Pour $v_0 \in [50, 1000]$, l'amorce hyperbolique converge en moins de 8 itérations sans stagnation ni clamping dur.
### Step 3:
- L'algorithme résout $F(\tau) = \Delta t$ avec une erreur résiduelle $|F(\tau) - \Delta t| \le 10^{-12} \cdot \Delta t$ sans jamais renvoyer silencieusement $\tau = 10^{-15}$.
### Step 4:
- Le seuil de collision radiale est strictement réversible par inversion du temps.
### Step 5:
- 1000/1000 orbites aléatoires (elliptiques excentriques, paraboliques, hyperboliques rapides) passent le test de réversibilité avec une erreur $\|\mathbf{r}(-\Delta t) - \mathbf{r}_0\| < 10^{-13}$.
- 142+/142+ tests unitaires passent avec succès (100%).
### Step 6:
- Commit git enregistré sur la branche `main`.

## Observations:
> [!IMPORTANT]
> La stricte positivité de la dérivée $F'(\tau) = r(\tau) > 0$ garantit que la fonction $F(\tau) = \int_0^\tau r(s)\mathrm{d}s$ est un difféomorphisme strictement croissant de $[0, \infty[$ vers $[0, \infty[$. L'équation $F(\tau) = \Delta t$ admet donc **une et une seule solution réelle unique**. L'échec d'Halley n'était qu'un artefact d'algorithme local, qu'un solveur de secours bracketé élimine de manière mathématiquement prouvée.
