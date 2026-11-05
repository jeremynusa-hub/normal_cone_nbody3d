# Implementation Plan of Iteration 113 : Implémentation et Certification Unitaire de la Couche 1 (core/) - Algèbre, État Étendu, Invariants et Géométrie de Jacobi
Created at: 2026-10-07T17:32:00+02:00

## Context:
1. **Explicit Request**: Développer et tester exhaustivement la Couche 1 (`core/`) de la bibliothèque `normal_cone_nbody_float64`, conformément à la feuille de route stratégique Roadmap V validée. Le développement doit STRICTEMENT respecter le protocole récurrent à 5 étapes :
   - Étape 1 : Identifier et relire les sections théoriques du papier (Chapitres 1, 2 et Sections 3.1–3.4 de [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/)).
   - Étape 2 : Relire l'ancien code associé (`models/normal_cone_integrator/cas6_nbody3d_async/core/`) et identifier ce qui est conforme, les 26 occurrences de `cas6` à éradiquer, et le découplage de `JacobiNode`.
   - Étape 3 : Établir un document de spécification claire du code à écrire avec les tests et attendus.
   - Étape 4 : Rédiger le code modulaire optimisé, décomposé en fonctions pures, avec abstractions ciblées (`dataclasses`, `typing.Protocol`) et en-tête descriptif citant les sections du papier de recherche.
   - Étape 5 : Rédiger et exécuter la suite de tests unitaires miroirs (`tests/unit/core/`) avec un taux de réussite de 100 %.
2. **Implicit Strategic Value**: Poser le socle algébrique et géométrique invariant sur lequel reposeront tous les solveurs analytiques (Couche 2), champs de potentiel (Couche 3) et intégrateurs symplectiques (Couche 4). Garantir que l'état étendu $\mathcal{P}_{\mathrm{ext}}$ et les transformations canoniques sont exempts de toute dérive numérique, sans allocation mémoire redondante lors des calculs cinématiques.
3. **Alignment Criteria**: Les 6 modules de `core/` et leurs 6 fichiers de tests unitaires miroirs doivent vérifier exactement les théorèmes formels Lean 4 associés (fermeture de Lie-Bargmann à $< 10^{-14}$, symplecticity Darboux, dimension adimensionnelle de $u_v$ via $\bar{\rho}_C \in [L]$, unimodularité $\mathrm{SL}(2, \mathbb{R})$ $\det = 1.0$), bannir tout vestige de `"cas6"` et respecter l'interdiction absolue d'adoucissement ($\varepsilon = 0, d_0 = 0$).

## Guidelines:
- **Method**: Application stricte du protocole à 5 étapes (Théorie $\to$ Audit du code hérité $\to$ Spécification $\to$ Implémentation modulaire $\to$ Tests unitaires).
- **Red Lines**:
  - Sanctuarisation absolue en lecture seule du répertoire [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/).
  - Zéro occurrence de `"cas6"`, `"cas 6"`, `"conormal deformation"`, `"conormal cone"`.
  - Zéro adoucissement ni troncature ($\varepsilon = 0, d_0 = 0$).
  - Aucun couplage ascendant : la Couche 1 ne peut dépendre d'aucune couche supérieure.
- **Architecture**:
  - `normal_cone_nbody_float64/core/constants.py` : constantes physiques, seuil de Hill $\eta_{\mathrm{crit}} = 3^{-1/3}$, hystérésis $\gamma = 0.6$, coefficients Blanes--Moan $\mathrm{BM}_6^{(11)}$ à 25 décimales.
  - `normal_cone_nbody_float64/core/state.py` : dataclass immuable `ExtendedPhaseState` sur $\mathcal{P}_{\mathrm{ext}} = T^*\mathbf{Q}[N] \times \mathbb{R}_t \times \mathbb{R}_E$.
  - `normal_cone_nbody_float64/core/invariants.py` : calcul des 10 générateurs de Lie-Bargmann ($\mathbf{P}, \mathbf{L}, \mathbf{G}, E$) et vérification d'algèbre de Poisson.
  - `normal_cone_nbody_float64/core/jacobi.py` : topologie immuable d'arbre binaire de Jacobi, masses réduites $\mu_v$, hyper-rayon normalisé $\bar{\rho}_C = \rho_C / \sqrt{M_C} \in [L]$, génération combinatoire des $(2N-3)!!$ arbres.
  - `normal_cone_nbody_float64/core/coordinates.py` : bijections symplectiques exactes $(\mathbf{q}, \mathbf{p}) \leftrightarrow (\mathbf{r}_v, \mathbf{w}_v)$ et rétraction de Rees à $t \to 0$.
  - `normal_cone_nbody_float64/core/transitions.py` : rotations d'arbres $\Psi_{A \to B} \in \operatorname{Sp}(6N-6, \mathbb{R})$ via transformations unimodulaires $\mathrm{SL}(2, \mathbb{R})$.
- **Code Headers**: Chaque fichier source débute par un bloc de documentation indiquant son rôle, ses équations de référence dans la monographie et ses ancres de théorèmes Lean 4 certifiés.

## Plan:
### Step 1: Relecture Théorique des Chapitres 1, 2 et Sections 3.1–3.4
- **Context**: La monographie dans `redaction/unsoftened_conormal_nbody3d/` formalise les espaces géométriques, l'algèbre de Lie-Bargmann et la décomposition de Fulton--MacPherson.
- **Objective**: Extraire les définitions, formules canoniques et identités exactes requises pour le socle `core/`.
- **Deliverables**:
  - `Layer1TheorySynthesis`: Synthèse des formules mathématiques exactes et des labels LaTeX correspondants.
- **Process**: Relecture ciblée de `chapters/ch01/content.tex`, `chapters/ch02/content.tex` et `chapters/ch03/content.tex`.

### Step 2: Audit Approfondi du Code Hérité (`cas6_nbody3d_async/core/`)
- **Context**: L'ancien code contient 6 modules fonctionnels dans `models/normal_cone_integrator/cas6_nbody3d_async/core/`.
- **Objective**: Recenser exhaustivement le code réutilisable, les 26 mentions de `cas6` à éliminer, le découplage nécessaire de `JacobiNode`, et la correction dimensionnelle de l'hyper-rayon $\bar{\rho}_C$.
- **Deliverables**:
  - `Layer1RefactoringMatrix`: Matrice détaillée des fonctions à migrer, assainir ou réécrire.
- **Process**: Inspection fichier par fichier des sources de l'ancien répertoire `core/`.

### Step 3: Spécification Formelle et Contrats d'Interface de la Couche 1
- **Context**: La conception doit être arrêtée avant l'écriture du code.
- **Objective**: Définir les signatures d'API précises, types NumPy stricts, structures de données immuables et critères d'acceptation de test unitaire pour chaque module.
- **Deliverables**:
  - `Layer1ModuleSpecification`: Spécifications des classes, fonctions pures, en-têtes et signatures.
- **Process**: Rédaction de la spécification technique conformément au document maître.

### Step 4: Implémentation Modulaire et Optimisée des 6 Modules de `core/`
- **Context**: Le répertoire cible `normal_cone_nbody_float64/core/` doit être créé.
- **Objective**: Rédiger les 6 modules Python avec leurs en-têtes descriptifs et leurs optimisations vectorielles pures.
- **Deliverables**:
  - `normal_cone_nbody_float64/__init__.py`
  - `normal_cone_nbody_float64/core/__init__.py`
  - `normal_cone_nbody_float64/core/constants.py`
  - `normal_cone_nbody_float64/core/state.py`
  - `normal_cone_nbody_float64/core/invariants.py`
  - `normal_cone_nbody_float64/core/jacobi.py`
  - `normal_cone_nbody_float64/core/coordinates.py`
  - `normal_cone_nbody_float64/core/transitions.py`
- **Process**: Implémentation modulaire séquentielle sans écrasement destructif, avec typage rigoureux.

### Step 5: Implémentation et Exécution des Tests Unitaires Miroirs (`tests/unit/core/`)
- **Context**: Chaque module Python possède son fichier de test unitaire miroir dédié.
- **Objective**: Rédiger et exécuter la suite de tests complète sous `pytest`, vérifiant 100 % des cas passants et la conservation des invariants mathématiques.
- **Deliverables**:
  - `tests/unit/core/test_constants.py`
  - `tests/unit/core/test_state.py`
  - `tests/unit/core/test_invariants.py`
  - `tests/unit/core/test_jacobi.py`
  - `tests/unit/core/test_coordinates.py`
  - `tests/unit/core/test_transitions.py`
- **Process**: Exécution de `pytest tests/unit/core/` via `run_command` et vérification des tolérances $< 10^{-14}$.

## Acceptance Criteria:
### Step 1:
- Identification complète des références de l'article : $\mathcal{P}_{\mathrm{ext}}$ (Sec. 2.1), 10 invariants (Sec. 2.3), décomposition de Jacobi (Sec. 3.1), transformations symplectiques (Sec. 3.2), atlas $\mathrm{SL}(2, \mathbb{R})$ (Sec. 3.4), constantes $\mathrm{BM}_6^{(11)}$ (Table 1).
### Step 2:
- Rapport d'audit de refactorisation confirmant la suppression des 26 occurrences de `cas6` et le découplage de `JacobiNode`.
### Step 3:
- Spécifications typées complètes pour les 6 modules.
### Step 4:
- 6 modules implémentés avec en-têtes conformes citant les sections du papier et les théorèmes Lean 4.
- Zéro adoucissement ($\varepsilon = 0, d_0 = 0$).
### Step 5:
- 100 % des tests unitaires de `tests/unit/core/` réussis avec succès.
- Conservation de l'algèbre de Lie-Bargmann prouvée à $< 10^{-14}$.
- Symplecticité de transition vérifiée à $< 10^{-15}$.

## Observations:
> [!IMPORTANT]
> Le découplage entre la topologie de l'arbre de Jacobi (`JacobiNode`) et les coordonnées dynamiques instantanées constitue un impératif architectural majeur issu du triple audit. `JacobiNode` devient une structure immuable purement combinatoire et géométrique, tandis que les positions et impulsions restent hébergées dans `ExtendedPhaseState`.
