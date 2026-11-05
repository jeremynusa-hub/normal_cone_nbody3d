# Implementation Plan of Iteration 117 : Couche 5 (`system/`) - Orchestrateur `NBodySystem`, Gestionnaire d'Atlas `AtlasManager` à Hystérésis $\gamma = 0.6$, Ordonnanceur Dyadique Multi-Taux `DyadicScheduler` et Traçage Temporel `HistoryTracker`
Created at: 2026-10-07T18:01:00+02:00

## Context:
1. **Explicit Request**: Développer, tester unitairement de manière exhaustive et certifier formellement la Couche 5 (`system/`) de la bibliothèque `normal_cone_nbody_float64` en suivant strictement le protocole mandaté en 5 étapes.
2. **Implicit Strategic Value**: Fournir la couche d'orchestration macroscopique et d'intégration système unifiée reliant les coordonnées géométriques (Couche 1), les solveurs analytiques (Couche 2), les potentiels d'apex régularisés (Couche 3) et la suite d'intégrateurs symplectiques (Couche 4) au sein d'une façade cohérente, robuste et modulaire, tout en garantissant la dynamique d'atlas sans bégaiement (théorème d'hystérésis $\gamma = 0.6$) et l'ordonnancement temporel asynchrone exact.
3. **Alignment Criteria**: Taux de réussite de 100 % sur l'ensemble de la suite de tests unitaires de la couche `tests/unit/system/` et cumulatif (100+ tests passés), absence totale de marqueurs `"cas6"`, zéro adoucissement ($\varepsilon = 0, d_0 = 0$), exactitude machine des transitions d'atlas dans $\mathrm{SL}(2, \mathbb{R})$, télescopage dyadique parfait et conservation stricte des 10 invariants de Lie--Bargmann.

## Guidelines:
- **Method**: Protocole rigoureux en 5 étapes (Théorie -> Audit hérité -> Spécification formelle -> Implémentation modulaire -> Tests unitaires 1:1).
- **Red Lines**: Interdiction absolue de toute modification dans `redaction/unsoftened_conormal_nbody3d/`. Éradication totale de `"cas6"` ou `"cas 6"`. Zéro adoucissement ni seuillage artificiel de distance ($\varepsilon = 0, d_0 = 0$).
- **Architecture**:
  - `system/nbody_system.py` : Façade d'orchestration de simulation encapsulant l'état étendu $\mathcal{P}_{\mathrm{ext}}$, l'arbre d'atlas courant, l'intégrateur choisi, le gestionnaire d'atlas et l'historique.
  - `system/atlas_manager.py` : Machine à états finis régissant la validité topologique de l'arbre de Jacobi, la détection des inversions de hiérarchie et l'application des transitions $\mathrm{SL}(2, \mathbb{R})$ avec bande d'hystérésis stricte $\gamma = 0.6$.
  - `system/dyadic_scheduler.py` : Ordonnanceur multi-taux asynchrone calculant les pas de temps dyadiques $h_v = 2^{-k_v} h_{\mathrm{macro}}$ selon la vitesse et la courbure locale, synchronisant les sous-blocs et garantissant le télescopage parfait des opérateurs.
  - `system/history.py` : Collecteur chronologique optimisé pour l'enregistrement économe des trajectoires, des états étendus et de la métrologie des 10 invariants de Lie--Bargmann.
- **Form**: Code, commentaires, docstrings et tests exclusivement en anglais avec citations de la monographie et des théorèmes Lean 4. Plan et liste de tâches en français markdown natif.

## Plan:
### Step 1: Relecture Théorique des Chapitres 4, 5 et 6 de la Monographie
- **Starting Material**: Monographie LaTeX `redaction/unsoftened_conormal_nbody3d/` (Chapitre 4 Section 4.3, Chapitre 5 Section 5.1 & 5.5, Chapitre 6 Section 6.1–6.3).
- **Objective**: Extraire les règles mathématiques exactes de transition d'atlas sous hystérésis $\gamma = 0.6$, le critère anti-bégaiement (Théorème 6.1), la règle de dérivation du pas dyadique $k_v = \lceil \log_2(h_{\mathrm{macro}} / \tau_v) \rceil$, et les théorèmes Lean 4 (`PaperUnifiedFoundations.lean`, `02_LieBargmannExactInvariance_Theorems.lean`).
- **Deliverables**: Spécification théorique consolidée pour la Couche 5 (`Component Theoretical Framework`).
- **Process**: Relecture ciblée des fichiers `.tex` et formalisation des invariants opérationnels.

### Step 2: Audit Approfondi du Code Hérité (`cas6_nbody3d_async/tree_atlas.py`, `cas6_nbody3d_async/hierarchical_cluster.py`)
- **Starting Material**: Code legacy dans `models/normal_cone_integrator/cas6_nbody3d_async/tree_atlas.py`, `strang_propagator.py`, et le rapport d'audit `python_codebase_mapping_and_audit_report.md`.
- **Objective**: Recenser les algorithmes validés de flip topologique, identifier les 32 occurrences résiduelles de `cas6`, éliminer les couplages directs avec l'état dynamique mutable, et concevoir de novo l'ordonnanceur dyadique formel.
- **Deliverables**: Rapport de refactorisation Couche 5 (`Component Legacy Audit`).
- **Process**: Confrontation ligne à ligne et élimination des limitations du prototype historique.

### Step 3: Spécification Formelle et Contrats d'Interface de la Couche 5
- **Starting Material**: Spécification globale `normal_cone_nbody_float64_master_architecture.md` (Section 3.5 & 6.5).
- **Objective**: Définir les types, signatures, exceptions, invariants et critères d'évaluation des 4 modules de `system/` et de leurs 4 bancs de tests miroirs.
- **Deliverables**: Contrats d'interface d'orchestration (`Component Interface Contracts`).
- **Process**: Formalisation des classes `NBodySystem`, `AtlasManager`, `DyadicScheduler`, `HistoryTracker`.

### Step 4: Implémentation Modulaire et Optimisée des 4 Modules de `system/`
- **Starting Material**: Modules certifiés des Couches 1 à 4 (`core/`, `solvers/`, `potential/`, `integrators/`).
- **Objective**: Implémenter les 4 fichiers de production sous `normal_cone_nbody_float64/system/` :
  1. `system/__init__.py` : Exports du package.
  2. `system/atlas_manager.py` : Détection des inversions de ratio de distance $r_v / r_{\mathrm{parent}} > \gamma = 0.6$, flips d'arbres unimodulaires $\mathrm{SL}(2, \mathbb{R})$, maintien de l'arbre adapté.
  3. `system/dyadic_scheduler.py` : Calcul des niveaux de raffinement dyadique $k \in [0, k_{\max}]$, planification des sous-pas, synchronisation palindromique et télescopage.
  4. `system/history.py` : Journalisation d'états compacts, calcul vectorisé des séries d'invariants et export de diagnostics.
  5. `system/nbody_system.py` : Façade intégratrice avec méthodes `step(dt)`, `evolve(t_span, dt)`, `get_state()`, `get_invariants()`, permutation d'intégrateurs sans perte de contexte.
- **Deliverables**: Fichiers de production sous `normal_cone_nbody_float64/system/`.
- **Process**: Écriture modulaire avec typage strict PEP 484/544, zéro allocation superflue et documentation académique complète.

### Step 5: Implémentation et Exécution des Tests Unitaires Miroirs (`tests/unit/system/`)
- **Starting Material**: Spécifications d'évaluation et modules de `system/`.
- **Objective**: Développer et exécuter avec 100 % de succès les 4 bancs de test unitaires miroirs :
  1. `tests/unit/system/test_atlas_manager.py` (détection d'hystérésis, anti-bégaiement, invariance symplectique des transitions).
  2. `tests/unit/system/test_dyadic_scheduler.py` (calcul des niveaux dyadiques, décomposition récursive de pas, synchronisation temporelle).
  3. `tests/unit/system/test_history.py` (enregistrement d'échantillons, invariance des dimensions, consistance des métriques).
  4. `tests/unit/system/test_nbody_system.py` (façade unifiée, boucle `evolve()`, conservation de Lie--Bargmann, interchangeabilité des intégrateurs).
- **Deliverables**: 4 fichiers de test sous `tests/unit/system/` avec 100 % de réussite.
- **Process**: Exécution via `pytest -v`, vérification des tolérances et du comportement asymptotique.

## Acceptance Criteria:
- **Critère 5.1 (Anti-Bégaiement & Hystérésis)** : `AtlasManager` n'effectue aucun flip oscillatoire répété sur des configurations oscillant autour du seuil grâce à la marge $\gamma = 0.6$ (Théorème 6.1).
- **Critère 5.2 (Symplecticité des Transitions)** : La transition d'atlas préserve les coordonnées canoniques et les 10 générateurs de Lie--Bargmann avec une erreur $< 10^{-14}$.
- **Critère 5.3 (Intégrité Temporelle Dyadique)** : `DyadicScheduler` garantit que la somme des sous-pas dyadiques recompose exactement le macro-pas $\sum 2^{-k} h = h_{\mathrm{macro}}$ sans résidu temporel ($|t_{\mathrm{final}} - T| < 10^{-15}$).
- **Critère 5.4 (Façade NBodySystem & Invariance)** : La méthode `system.evolve(t_span=1.0)` maintient la conservation de l'impulsion $P_{\mathrm{tot}}$ et du moment cinétique $L_{\mathrm{tot}}$ à $< 10^{-13}$ quel que soit l'intégrateur utilisé.
- **Critère 5.5 (Interchangeabilité des Intégrateurs)** : Le système permet de commuter à chaud entre $\mathcal{S}_2, \mathcal{S}_4, \mathrm{BM}_4, \mathrm{BM}_6$ et $\mathcal{P}\mathrm{BM}$ sans réinitialisation destructive.
- **Critère 5.6 (Zero-Softening & Zero-Cas6)** : 0 occurrence de `"cas6"` ou `"cas 6"`, 0 softening ($\varepsilon = 0, d_0 = 0$).
- **Critère 5.7 (Taux de Succès Global)** : 100 % de réussite sur `tests/unit/system/` et sur la suite cumulative (Core + Solvers + Potential + Integrators + System, soit plus de 105 tests).

## Observations:
> [!IMPORTANT]
> L'ordonnanceur dyadique doit impérativement respecter la parité palindromique des pas pour préserver la réversibilité temporelle microscopique du schéma global.
> [!NOTE]
> Le gestionnaire d'atlas `AtlasManager` opère de manière purement fonctionnelle sur la structure immuable `JacobiTree`, séparant l'évaluation topologique de l'état dynamique $\mathcal{P}_{\mathrm{ext}}$.
