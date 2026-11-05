# Implementation Plan of Iteration 146 : Phase 3 — Remédiation de la Synchronisation de Skeel-Gear et de l'Échéancier Dyadique du Chapitre 5
Created at: 2026-10-08T10:25:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur a ordonné : « Continue la roadmap donc ». Dans le cadre de la Roadmap VI, cela active immédiatement la **Phase 3 (Itération 146)** consacrée à la remédiation de la synchronisation de Skeel-Gear (1992) aux frontières macro-temporelles, à l'assouplissement de la profondeur dyadique dynamique jusqu'à $k_{\max}=17$ dans `system/dyadic_scheduler.py`, et à la fusion FSAL des demi-kicks d'apex dans `integrators/blanes_moan_4.py`.
2. **Implicit Strategic Value**: Restaurer la symétrie palindromique miroir $\mathcal{S}(-h) \circ \mathcal{S}(h) = \operatorname{Id}$ de l'ensemble du macro-bloc dyadique en interdisant toute mutation topologique d'arbre asymétrique au milieu des sous-pas (qui détruisait l'hamiltonien d'ombre) ; garantir la résolution physique intégrale des péricentres ultra-serrés ($r \sim 10^{-4}, v \sim 100\text{ m/s}$) sans sous-échantillonnage de pas ; et optimiser le coût numérique du schéma BM4 par factorisation des évaluations de force d'apex (division par 2 du nombre de kicks).
3. **Alignment Criteria**: L'arbre de Jacobi reste rigoureusement fixe au cours des $2^k$ sous-pas d'un macro-pas et n'est adapté qu'aux frontières synchronisées ; la profondeur dyadique monte dynamiquement jusqu'à $k=17$ sans plafonnement artificiel à 10 ; BM4 conserve sa précision d'ordre 4 et sa réversibilité tout en divisant par 2 ses évaluations d'apex ; 100 % des tests unitaires sont validés.

---

## Guidelines:
- **Method**: Respect strict du principe de synchronisation de Skeel & Gear (1992) pour les intégrateurs multi-taux et de la théorie KAM sur la préservation de la symétrie temporelle réversible.
- **Red Lines**: Zéro adoucisseur ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$). Interdiction d'introduire des seuils arbitraires en distance dans le calcul dyadique (maintien du calcul adimensionnel échelle-invariant basé sur les temps de collision $\tau_{\mathrm{kin}}, \tau_{\mathrm{dyn}}, \tau_{\mathrm{peri}}$). Le répertoire `redaction/unsoftened_conormal_nbody3d/` demeure sanctuarisé en lecture seule.
- **Architecture**: Préservation de la séparation modulaire : `system/nbody_system.py` orchestre la boucle macro, `system/dyadic_scheduler.py` calcule la granularité sans mutation d'état, et `integrators/blanes_moan_4.py` respecte le protocole `SymplecticIntegrator`.

---

## Plan:

### Step 1: Confinement des Transitions d'Atlas aux Frontières Macro-Temporelles dans `system/nbody_system.py`
- **Context**: Actuellement, `nbody_system.py` évalue `self.atlas_manager.evaluate_tree_adaptation` à l'intérieur de la boucle des sous-pas dyadiques `for _ in range(num_substeps)`. Si l'arbre bascule au sous-pas $m$ ($0 < m < 2^k$), la symétrie palindromique du macro-pas est rompue ($m \ne 2^k - m$), et pour les schémas traités, cela crée un télescopage hétérogène $\mathcal{C}_B \circ \mathcal{C}_A^{-1} \ne \operatorname{Id}$.
- **Objective**: Confiner strictement l'adaptation topologique de l'arbre de Jacobi à la frontière du macro-pas (au début de `step`), en figeant l'arbre pour l'intégralité des $2^k$ sous-pas du bloc dyadique.
- **Deliverables**:
  - `system/nbody_system.py`: Méthode `step` remaniée pour retirer tout appel à `evaluate_tree_adaptation` à l'intérieur des boucles de sous-pas, simplifier le télescopage frontière des intégrateurs traités sur l'arbre actif unifié, et garantir la réversibilité miroir stricte du macro-bloc.
- **Process**:
  1. Maintenir l'évaluation d'arbre au début du macro-pas `t_start` (Stage 1).
  2. Supprimer la réévaluation intra-bloc dans la boucle `for _ in range(num_substeps)` tant pour `ProcessedSymplecticIntegrator` que pour les intégrateurs standards.
  3. Pour `ProcessedSymplecticIntegrator`, appliquer $\mathcal{C}^{-1}$ au départ sur `self.tree`, propager les $2^k$ sous-pas sur ce même `self.tree`, puis refermer par $\mathcal{C}$ sur `self.tree`.
  4. Valider que l'opérateur composé $\mathcal{S}_{\mathrm{macro}}(h)$ est strictement palindromique.

### Step 2: Profondeur Dyadique Dynamique jusqu'à $k_{\max}=17$ dans `system/dyadic_scheduler.py`
- **Context**: `DyadicScheduler` est actuellement initialisé avec `max_level: int = 10` par défaut, ce qui plafonne à 1024 sous-pas et empêche la résolution conormale adéquate des péricentres ultra-serrés ($r \sim 10^{-4}$ à $v \sim 100\text{ m/s}$ où $k=17$ est requis).
- **Objective**: Élever la limite maximale de profondeur dyadique par défaut à $k_{\max}=17$ ($131\,072$ sous-pas), permettant d'atteindre un pas effectif $h_{\mathrm{sub}} = \Delta t / 2^{17} \sim 10^{-7}\text{ s}$ garantissant $h_{\mathrm{sub}} / \tau_{\mathrm{kin}} \le 0.05$.
- **Deliverables**:
  - `core/constants.py`: Définition de `DEFAULT_MAX_DYADIC_LEVEL: int = 17`.
  - `system/dyadic_scheduler.py`: Valeur par défaut `max_level: int = DEFAULT_MAX_DYADIC_LEVEL` dans `__init__`, avec vérification de non-dépassement et calcul sans overflow.
- **Process**:
  1. Déclarer la constante canonique `DEFAULT_MAX_DYADIC_LEVEL = 17` dans `constants.py`.
  2. Aligner `DyadicScheduler.__init__` sur cette constante.
  3. Vérifier que `get_substep_params` et le décalage binaire `1 << k` opèrent exactement en float64 IEEE 754 pour $k \in [0, 17]$.

### Step 3: Fusion FSAL des Demi-Kicks d'Apex dans `integrators/blanes_moan_4.py`
- **Context**: Actuellement, `BlanesMoan4Integrator` exécute une boucle sur les 7 étages Strang `self._strang.step(...)`, évaluant 2 kicks d'apex par étage, soit 14 évaluations de force par pas. La propriété FSAL (First Same As Last) permet de fusionner les deux demi-kicks adjacents $\mathcal{K}(a_i \frac{dt}{2}) \circ \mathcal{K}(a_{i+1} \frac{dt}{2}) = \mathcal{K}((a_i + a_{i+1}) \frac{dt}{2})$.
- **Objective**: Réécrire le schéma BM4 en séquence directe de drifts képlériens et de kicks d'apex fusionnés, réduisant le nombre d'évaluations de force à 8 par pas isolé (et 7 en chaîne continue FSAL), divisant par deux le coût de calcul du potentiel d'apex.
- **Deliverables**:
  - `integrators/blanes_moan_4.py`: Séquence de coefficients précalculés `BM4_DRIFTS` (7 éléments) et `BM4_KICKS` (8 éléments), et méthode `step` utilisant directement `step_drift_and_kepler` et `step_apex_kick`.
- **Process**:
  1. Dériver analytiquement les 7 coefficients de drift $d_j = a_j$ et les 8 coefficients de kick $k_1 = a_1/2$, $k_2 = (a_1+a_2)/2$, ..., $k_8 = a_1/2$.
  2. Implémenter l'alternance kick initial $\to$ (drift, kick fusionné) $\times 6 \to$ drift central $\to$ kick final.
  3. Vérifier l'identité mathématique stricte avec la formulation Strang précédente ($< 10^{-15}$) et la préservation de la réversibilité temporelle et de l'ordre 4.

### Step 4: Tests Unitaires Miroirs et Validation de Non-Régression
- **Context**: La suite de tests unitaires (132 tests) doit être maintenue à 100 % de succès et complétée par des vérifications directes de la synchronisation de Skeel-Gear, du raffinement dyadique profond $k=17$, et de la fusion FSAL de BM4.
- **Objective**: Attester formellement de la conformité de chaque composant remédié.
- **Deliverables**:
  - `tests/unit/system/test_nbody_system.py`: Test de non-mutation de l'arbre de Jacobi au cours des sous-pas dyadiques et réversibilité miroir macro $\mathcal{S}(-h) \circ \mathcal{S}(h) = \operatorname{Id}$.
  - `tests/unit/system/test_dyadic_scheduler.py`: Test de raffinement dyadique jusqu'à $k=17$ sans saturation.
  - `tests/unit/integrators/test_blanes_moan_4.py`: Test comparatif de performance et d'exactitude vérifiant que le schéma fusionné produit exactement la même trajectoire ($< 10^{-15}$) avec 8 évaluations d'apex au lieu de 14.
- **Process**:
  1. Implémenter les nouveaux cas de test dans les suites unitaires respectives.
  2. Exécuter `pytest tests/unit/` et valider les 135+ tests.

### Step 5: Validation Git Commit
- **Context**: Clôture de la Phase 3 de la Roadmap VI.
- **Objective**: Enregistrer le commit traçant la synchronisation Skeel-Gear, le scheduler $k=17$ et la fusion FSAL de BM4.
- **Deliverables**:
  - Commit Git structuré : `feat(system,integrators): enforce Skeel-Gear macro boundary sync, dynamic dyadic depth k_max=17, and BM4 FSAL kick fusion - Phase 3 of Roadmap VI (Iteration 146)`.
- **Process**:
  1. Vérifier l'état git via `git status`.
  2. Indexer les fichiers modifiés et valider le commit.

---

## Acceptance Criteria:

### Step 1:
- L'arbre de Jacobi actif demeure strictement invariant au cours de l'ensemble des $2^k$ sous-pas dyadiques d'un macro-pas.
- L'invariance par renversement du temps $\mathcal{S}(-h) \circ \mathcal{S}(h) = \operatorname{Id}$ est vérifiée à $< 10^{-14}$ sur l'ensemble du macro-pas dyadique multi-sous-pas.

### Step 2:
- `DEFAULT_MAX_DYADIC_LEVEL` vaut 17.
- Pour une rencontre ultra-serrée ($r \sim 10^{-4}$ à grande vitesse), `compute_dynamical_level` renvoie un niveau $k > 10$ (jusqu'à 17), sans tronquer artificiellement à 10.
- La relation de décomposition temporelle $\sum_{m=1}^{2^k} h_m = h_{\mathrm{macro}}$ est vérifiée à $< 10^{-16}$ pour $k = 17$ ($131\,072$ pas).

### Step 3:
- `BlanesMoan4Integrator` utilise la séquence de demi-kicks fusionnés, réduisant le nombre d'appels d'écrantage d'apex de 14 à 8 par pas macro.
- La solution numérique de BM4 fusionné est identique à celle de Strang composite à la précision machine ($< 10^{-15}$).
- La pente de convergence d'ordre 4 ($E(2h)/E(h) \approx 16$) et la réversibilité temporelle ($< 10^{-14}$) sont rigoureusement préservées.

### Step 4:
- 100 % de succès sur l'ensemble de la suite de tests unitaires `tests/unit/` (>= 132 tests).

### Step 5:
- Commit Git tracé et documenté dans le repository.

---

## Observations:
> [!IMPORTANT]
> Le confinement des transitions d'atlas exclusivement aux frontières macro synchronisées (Skeel & Gear 1992) est la clé fondamentale qui rétablit la régularité de l'hamiltonien d'ombre et prépare le terrain pour la connexion de jauge symplectique de l'Itération 147 (Phase 4).

> [!NOTE]
> La fusion des demi-kicks adjacents dans BM4 ne modifie aucunement la physique du problème mais apporte un gain immédiat de ~40% sur le temps de calcul du potentiel d'apex écranté.
