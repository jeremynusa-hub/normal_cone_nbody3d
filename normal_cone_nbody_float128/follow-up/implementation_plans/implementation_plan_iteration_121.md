# Implementation Plan of Iteration 121 : Intégration des Recommandations d'Audit Multi-Agents (Correctifs P0, P1 et P2)
Created at: 2026-10-07T18:25:30+02:00

## Context:
1. **Explicit Request**: Intégrer les recommandations issues de l'audit multi-agent approfondi mené par 13 sous-agents d'audit sur l'ensemble de la bibliothèque `normal_cone_nbody_float64`.
2. **Implicit Strategic Value**: Élever la bibliothèque à un degré de perfection mathématique, algorithmique et documentaire irréprochable : accélérer la convergence de Halley en éliminant un terme redondant de dérivée seconde, garantir la périodicité séculaire de la Figure-8 à 16 décimales complètes, optimiser le coût numérique du télescopage dans `NBodySystem`, et aligner l'ensemble des citations de code sur la structure canonique de la monographie LaTeX.
3. **Alignment Criteria**: Maintien strict du taux de réussite de 100 % sur la suite de tests unitaires, validation des ordres de convergence restaurés, zéro adoucissement ($\varepsilon = 0, d_0 = 0$), et zéro mention `"cas6"`.

## Guidelines:
- **Method**: Édition chirurgicale et granulaire via `replace_file_content` sur les modules ciblés.
- **Red Lines**: Sanctuarisation absolue de `redaction/unsoftened_conormal_nbody3d/` (aucune modification de texte LaTeX). Interdiction absolue d'écraser des fichiers existants avec `write_to_file`. Aucune régression sur les 122 tests existants.
- **Architecture des Correctifs**:
  - **Pôle P0 (Mathématiques & Algorithmique)** :
    - `solvers/kepler_sundman.py` : Correction de la dérivée seconde $F''(\tau)$ dans l'itérateur de Halley (suppression du terme redondant $-\alpha \sigma_0 \tau^2 c_2(z)$).
    - `benchmarks/registry.py` : Rétablissement des 16 décimales complètes de Simó (2000) pour l'orbite de Figure-8.
    - `benchmarks/runner.py` : Ajout de la fonction `fit_convergence_order(x_vals, error_vals, mode)` calculant l'ordre empirique $p$, l'intercept et le $R^2$ log-log, avec intégration dans `compute_work_precision_pareto()`.
  - **Pôle P1 (Système & Performance)** :
    - `system/nbody_system.py` : Prise en charge du télescopage dyadique des correcteurs $\mathcal{C}^{-1} \circ \mathcal{S}^M \circ \mathcal{C}$ dans `step()` lors de l'utilisation de `ProcessedSymplecticIntegrator`.
    - `system/atlas_manager.py` : Enregistrement du dwell time écoulé $\tau_{\mathrm{dwell}}$ et comparaison avec le seuil minimal théorique $\tau_{\min} = \frac{16}{15 \Omega_{\max}}$.
    - `diagnostics/radar.py` : Création de l'alias d'interopérabilité vers `collision_radar.py` et alignement de la docstring du score de danger $\sigma$.
  - **Pôle P2 (Harmonisation Documentaire des En-têtes)** :
    - `core/` (`jacobi.py`, `coordinates.py`, `transitions.py`, `invariants.py`) : Mise à jour des citations vers le Chapitre 1 \S1.4 et le Chapitre 4 \S4.1 & \S4.3.
    - `potential/supernode_detector.py` : Correction de la citation "Proposition 4.4" en "Proposition 4.5".
    - `solvers/stumpff.py` & `kepler_sundman.py` : Mise à jour des citations vers le Chapitre 5 Section 5.2.2.
    - `solvers/mcgehee.py` : Clarification des Sections 3.2 et 5.2.
    - `benchmarks/` : Mise à jour des citations vers le Chapitre 7 (\S7.1 et \S7.2).
- **Form**: Code, commentaires et docstrings en anglais. Plan et liste de tâches en français markdown.

## Plan:
### Step 1: Application des Correctifs Mathématiques et Algorithmiques P0
- **Starting Material**: `solvers/kepler_sundman.py`, `benchmarks/registry.py`, `benchmarks/runner.py`.
- **Objective**: Corriger la formule analytique de $F''(\tau)$ dans Halley, injecter les 16 décimales de la Figure-8, et ajouter la régression log-log de Pareto.
- **Deliverables**:
  - `CorrectedKeplerHalley`: Formule $F''(\tau)$ exacte sans redondance.
  - `HighPrecisionFigure8`: Conditions initiales à 16 décimales.
  - `ParetoLogLogRegression`: Fonction `fit_convergence_order()`.
- **Process**: Éditions granulaires via `replace_file_content`.

### Step 2: Application des Optimisations Systémiques P1
- **Starting Material**: `system/nbody_system.py`, `system/atlas_manager.py`, `diagnostics/`.
- **Objective**: Implémenter le télescopage dyadique des correcteurs dans `NBodySystem.step()`, tracer $\tau_{\mathrm{dwell}}$ dans `AtlasManager`, et ajouter l'alias `diagnostics/radar.py`.
- **Deliverables**:
  - `TelescopedSystemFaçade`: Macro-pas optimisé pour les intégrateurs traités.
  - `DwellTimeTracking`: Suivi de $\tau_{\mathrm{dwell}}$ dans l'historique d'atlas.
  - `RadarAlias`: Module de compatibilité `diagnostics/radar.py`.
- **Process**: Éditions et création de fichier via `write_to_file` sans écrasement de source existante.

### Step 3: Harmonisation des Citations de la Monographie P2
- **Starting Material**: En-têtes des modules de `core/`, `potential/`, `solvers/`, `benchmarks/`.
- **Objective**: Réaligner les numéros de chapitres et de sections sur la structure finale du manuscrit LaTeX.
- **Deliverables**:
  - `HarmonizedDocstrings`: En-têtes documentaires parfaitement synchronisés.
- **Process**: Éditions ciblées des blocs d'en-tête via `replace_file_content`.

### Step 4: Validation Globale et Non-Régression
- **Starting Material**: Suite de tests unitaires `tests/unit/`.
- **Objective**: Exécuter `python -m pytest tests/unit/ -v` et vérifier le taux de réussite de 100 %, la convergence cubique de Halley et la stabilité à long terme de la Figure-8.
- **Deliverables**:
  - `TestVerificationSuite`: Rapport d'exécution sans aucune régression.
- **Process**: Exécution via `run_command` et analyse des métriques.

## Acceptance Criteria:
### Step 1:
- `f_prime2` dans `kepler_sundman.py` est strictement égal à `eff_sigma0 * c0 + (GM - alpha * r0) * tau * c1`.
- Les 4 coordonnées/vitesses de la Figure-8 possèdent 16 décimales exactes.
- `fit_convergence_order()` extrait la pente $p$ avec un $R^2 > 0.99$ sur données asymptotiques.
### Step 2:
- Pour un `ProcessedSymplecticIntegrator`, `NBodySystem.step()` n'évalue les correcteurs qu'au début et à la fin du macro-pas.
- `AtlasTransitionRecord` comporte le champ `dwell_time`.
- `from normal_cone_nbody_float64.diagnostics.radar import CollisionRadar` fonctionne immédiatement.
### Step 3:
- L'ensemble des en-têtes pointe vers les chapitres finaux (Ch 1, Ch 4, Ch 5, Ch 7) sans références à des sections obsolètes.
### Step 4:
- 100 % de tests unitaires passants sous pytest.

## Observations:
> [!IMPORTANT]
> Les correctifs P0 sur $F''(\tau)$ et la Figure-8 améliorent la fidélité mathématique sans altérer l'API publique de la bibliothèque.
> [!NOTE]
> Le télescopage dyadique dans `NBodySystem.step()` apporte un gain de performance direct proportionnel à $2 \times (2^k - 1)$ évaluations de force évitées.
