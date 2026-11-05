# Implementation Plan of Iteration 151 : Régularisation Canonique de l'Intégrateur Interne de Super-Nœuds sous Temps Déplié
Created at: 2026-10-08T12:17:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur poursuit la Roadmap VII vers la Phase 2 (Itération 151) pour éradiquer la Divergence Critique 2 identifiée par l'audit contradictoire à 9 sous-agents : la défaillance du solveur interne de super-nœuds `normal_cone_nbody_float64/solvers/supernode.py` en cas de quasi-collision binaire interne ($b \le 0.02$).
2. **Implicit Strategic Value**: Les super-nœuds synchrones (Chapitre 4, Sections 4.3 & 4.4) constituent le rempart mathématique contre la divergence du commutateur en $\mathcal{O}(\rho_S^{-5/2})$. Si une quasi-collision binaire survient au sein d'un super-nœud ($k \ge 3$), l'hyper-rayon global $\rho_S$ reste borné inférieurement par la distance aux spectateurs, de sorte que l'horloge globale $dt_S = \rho_S d\tau_S$ ne ralentit pas assez pour résoudre la singularité locale $r_{ij} \to 0$. L'ancien Leapfrog physique subissait alors une explosion d'erreur en $\Delta E/E \propto b^{-2}$ ou un épuisement silencieux du compteur de pas (`t_accum < abs_dt`). La régularisation canonique par isolation analytique sous Stumpff-Sundman de la sous-paire compacte ou par régularisation temporelle locale garantit l'intégrité symplectique et l'absence totale de divergence en $b^{-2}$.
3. **Alignment Criteria**:
   - Sanctuarisation absolue de `redaction/unsoftened_conormal_nbody3d/` (lecture seule).
   - Zéro adoucissement résiduel ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$).
   - Protection mathématique stricte contre la division par zéro ($d_{ij} = 0$) dans `_compute_internal_cluster_acc`.
   - Interdiction formelle de retour silencieux incomplet : levée explicite d'une exception `RuntimeError` si `t_accum < abs_dt`.
   - Conservation de l'énergie et stabilité de trajectoire vérifiées sur un banc d'essai de quasi-collisions à paramètre d'impact arbitrairement petit ($b \in [10^{-1}, 10^{-6}, 0.0]$).
   - 100% des tests unitaires validés (au moins 143/143 tests passés).

## Guidelines:
- **Méthode** :
  - Dans `_compute_internal_cluster_acc`, substituer au calcul direct `inv_d3 = G / (d ** 3)` un filtrage vectoriel strict `inv_d3 = np.where(d2 > 0.0, G / (d2 * d), 0.0)` pour garantir l'absence totale de `NaN` même en cas de choc ponctuel direct.
  - Décomposer l'interaction interne du super-nœud en séparant la sous-paire la plus serrée $\min_{i<j} r_{ij}$ intégrée analytiquement par l'opérateur universel de Stumpff-Sundman `step_kepler_sundman`, et les interactions résiduelles intégrées par split symplectique de Strang, ou appliquer un schéma temps-transformé canonique réversible.
  - Sécuriser la boucle temporelle : si `t_accum < abs_dt` lorsque `steps == max_steps`, lever immédiatement une `RuntimeError` circonstanciée plutôt que de renvoyer un état transitoire corrompu.
- **Lignes Rouges** :
  - Aucun adoucissement de distance ($d_{\mathrm{safe}} = 0$).
  - Aucun contournement par division masquée qui altère la conservation du moment linéaire ($\sum \mathbf{p}_i = \text{cste}$) ou angulaire ($\sum \mathbf{q}_i \times \mathbf{p}_i = \text{cste}$).
- **Architecture** :
  - Préserver l'interface publique `step_supernode_synchronous` et la structure `SupernodeCluster`.

## Plan:
### Step 1: Diagnostic et Isolation de l'Explosion d'Énergie et de l'Épuisement de Pas
- **Context**: L'audit et le test préliminaire ont révélé qu'une quasi-collision binaire interne à $b = 10^{-3}$ épuise les 100 000 pas en atteignant seulement $t \approx 0.099$ sur $0.2$, éjectant un état non convergé à $v \approx 250$.
- **Objective**: Écrire un script de diagnostic quantifiant précisément la dérive d'énergie et l'arrêt prématuré de la boucle temporelle pour $b \in [10^{-1}, 10^{-3}, 10^{-6}, 0.0]$.
- **Deliverables**:
  - `SupernodeDiagnosticScript`: Script de test diagnostique dans `scratch/diagnose_supernode_collision.py`.
- **Process**: Évaluer `step_supernode_synchronous` sur des triplets avec paramètres d'impact décroissants.

### Step 2: Sécurisation Anti-Division par Zéro et Exception d'Incomplétude Temporelle
- **Context**: En cas de choc géométrique strict ($d_{ij} = 0$), `inv_d3 = G / (d**3)` produit une division par zéro et injecte des `NaN`. De plus, la sortie de boucle à `steps == max_steps` n'est pas signalée.
- **Objective**: Protéger `_compute_internal_cluster_acc` avec `np.where(d2 > 0.0, G / (d2 * d), 0.0)` et ajouter une exception explicite si `t_accum < abs_dt`.
- **Deliverables**:
  - `ZeroDivisionGuard`: Filtrage vectoriel dans `_compute_internal_cluster_acc`.
  - `IncompleteStepException`: Contrôle d'incomplétude temporelle avec levée de `RuntimeError`.
- **Process**: Mettre à jour `_compute_internal_cluster_acc` et la terminaison de boucle dans `solvers/supernode.py`.

### Step 3: Régularisation Canonique par Découplage de Sous-Paire Compacte (Stumpff-Sundman)
- **Context**: Dans un amas de $k \ge 3$ corps, si une paire $(i, j)$ plonge vers une quasi-collision ($r_{ij} \to 0$), l'hyper-rayon global $\rho_S$ reste dominé par les autres corps distants. L'intégration globale en temps physique conduit à une instabilité en $b^{-2}$.
- **Objective**: Implémenter une décomposition symplectique locale à l'intérieur du super-nœud : isoler dynamiquement la sous-paire la plus serrée $(i^*, j^*)$ et propager son mouvement relatif par l'opérateur universel désingularisé `step_kepler_sundman`, tandis que le barycentre de la paire et les corps tiers évoluent sous les forces de marée lisses via Strang splitting.
- **Deliverables**:
  - `SupernodePairwiseRegularizer`: Moteur de propagation régularisé dans `step_supernode_synchronous`.
- **Process**: Refondre la boucle interne de `step_supernode_synchronous` en exploitant `step_kepler_sundman`.

### Step 4: Banc d'Essai de Collision Binaire Interne ($b \in [10^{-1}, 10^{-6}, 0.0]$)
- **Context**: La suite existante `tests/unit/solvers/test_supernode.py` ne teste que des configurations avec vitesses divergentes ou homothéties pures.
- **Objective**: Ajouter des tests unitaires rigoureux avec quasi-collisions frontales pour $b \in [0.1, 10^{-3}, 10^{-6}, 0.0]$, prouvant la conservation d'énergie, l'absence de NaN et la continuité temporelle.
- **Deliverables**:
  - `TestSupernodeBinaryCollision`: Nouvelles fonctions de test dans `tests/unit/solvers/test_supernode.py`.
- **Process**: Écrire et exécuter les tests via `pytest`.

### Step 5: Certification de Non-Régression Globale
- **Context**: 143 tests unitaires existants doivent tous passer avec succès.
- **Objective**: S'assurer qu'aucune régression n'est introduite sur les benchmarks 4-corps, Burrau, et les propriétés de symétrie/invariance de Lie-Bargmann.
- **Deliverables**:
  - `FullTestSuitePass`: Exécution complète de `pytest tests/unit/`.
- **Process**: Lancer la suite de tests et vérifier le taux de succès à 100%.

### Step 6: Git Commit de Phase 2
- **Context**: Phase 2 de la Roadmap VII complétée.
- **Objective**: Enregistrer un commit atomique documenté.
- **Deliverables**:
  - `GitCommit`: Commit git avec message standardisé.
- **Process**: `git add` et `git commit`.

## Acceptance Criteria:
### Step 1:
- Le comportement de divergence en $b^{-2}$ et l'épuisement de `max_steps` sont quantifiés et reproductibles dans `scratch/diagnose_supernode_collision.py`.
### Step 2:
- Aucun `RuntimeWarning` ni `NaN` n'est généré même si deux corps sont placés à la même coordonnée exacte ($d_{ij} = 0.0$).
- Si `steps == max_steps` avec `t_accum < abs_dt`, une exception explicite `RuntimeError` est levée.
### Step 3:
- Lors d'une quasi-collision binaire interne à $b \le 10^{-4}$, l'erreur relative d'énergie reste bornée ($< 10^{-4}$) au lieu d'exploser à $10^{15}$.
### Step 4:
- Les nouveaux tests pour $b \in [0.1, 10^{-3}, 10^{-6}, 0.0]$ passent avec succès dans `test_supernode.py`.
### Step 5:
- 100% des tests unitaires (144+ tests) passent avec succès.
### Step 6:
- Commit git propre enregistré sur `main`.

## Observations:
> [!IMPORTANT]
> Dans un amas à $k \ge 3$ corps, la singularité d'une quasi-collision binaire interne ($r_{ij} \to 0$) est une sous-variété de codimension 3 du cône normal où $\rho_S > 0$. L'isolation analytique de cette sous-paire par `step_kepler_sundman` (qui résout exactement le difféomorphisme $dt = r_{ij} ds$) transpose la désingularisation locale de Sundman au cœur même de l'amas synchrone sans compromettre la synchronisation globale avec les spectateurs.
