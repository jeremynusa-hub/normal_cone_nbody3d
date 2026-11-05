# Implementation Plan of Iteration 110 : Analyse Systématique & Indexation Exhaustive du Codebase Existant pour la Refactorisation en Bibliothèque Scientifique Autonome
Created at: 2026-10-07T16:26:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur demande de planifier l'analyse du code existant (`models/normal_cone_integrator/cas6_nbody3d_async/`) et son indexation exhaustive dans un fichier Markdown du brain, en vue de restructurer le code en une version bibliothèque ("lib") publiable dans l'article académique, scientifiquement exacte, au code et commentaires en anglais, modulaire, exempte de référence à "cas 6", sans tricks numériques, permettant de tester facilement différents intégrateurs sur une même simulation et d'exécuter des benchmarks préconfigurés. L'utilisateur demande ensuite de produire un document de rapport dans la vue.
2. **Implicit Strategic Value**: Transformer un prototype de recherche historique hautement performant mais organiquement intriqué (marqué par les étapes successives des roadmaps I à IV) en une bibliothèque logicielle autonome d'excellence académique, dotée d'une architecture propre, d'abstractions rigoureuses conformes aux théorèmes certifiés en Lean 4, et directement citable/reproductible par la communauté de physique mathématique et de mécanique céleste.
3. **Alignment Criteria**: L'indexation et l'analyse devront couvrir 100 % des composants actuels sans angle mort, isoler précisément chaque écart (terminologie, tricks, langue, couplages rigides, absence d'interface d'échange des intégrateurs), et poser les spécifications complètes de l'API cible avant toute modification de code.

## Guidelines:
- **Method**: Audit d'ingénierie logicielle systématique couplé à une analyse de conformité mathématique stricte (confrontation directe du code face aux monographies théoriques et aux modules formels certifiés).
- **Red Lines**:
  - Sanctuarisation absolue du répertoire [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/) : interdiction formelle d'y écrire ou de le modifier.
  - Zéro régression sur le principe d'airain du Zéro Clamping ($\varepsilon = 0, d_0 = 0$).
  - Aucun contournement par des forces ad-hoc ou des pseudo-régularisations empiriques.
  - Aucune modification destructive prématurée du code de production existant dans cette itération d'analyse et de planification.
- **Architecture**: Conception orientée objet et fonctionnelle pure pour la bibliothèque cible, privilégiant le *Strategy pattern* pour les schémas d'intégration, le *Registry pattern* pour les scénarios célestes, et des protocoles typés stricts (`Protocol`, `dataclasses`, `numpy.typing`).

## Plan:
### Step 1: Inventaire Exhaustif et Cartographie Structurelle du Codebase Existant
- **Context**: Le dossier `models/normal_cone_integrator/cas6_nbody3d_async/` comporte 11 modules centraux, ainsi que les sous-dossiers `solar/`, `benchmarks/` et `tests/`.
- **Objective**: Recenser et cartographier exhaustivement chaque module Python : localisation, taille, classes, fonctions publiques, entrées/sorties et matrice de dépendances mutuelles.
- **Deliverables**:
  - `CodebaseInventory`: Matrice complète des modules et graphe des flux de données.
- **Process**: Inspection systématique des fichiers via `view_file` et extraction des signatures publiques et imports croisés.

### Step 2: Audit Critique et Rédaction de l'Indexation Exhaustive dans le Brain (`codebase_analysis_and_index.md`)
- **Context**: L'indexation doit servir de cartographie de vérité terrain pour guider la refactorisation de la bibliothèque.
- **Objective**: Analyser chaque composant sous le prisme des 6 exigences de l'utilisateur et compiler les résultats dans `codebase_analysis_and_index.md` dans le brain.
- **Deliverables**:
  - `CodebaseAnalysisIndex`: Fichier Markdown persistant `C:\Users\jerem\.gemini\antigravity\brain\59313b1a-0dde-45d6-a476-1bcb513cf82e\codebase_analysis_and_index.md`.
- **Process**:
  1. *Critère 1 (Purger "cas 6")* : Dénombrer et localiser toutes les occurrences dans les noms de modules, fonctions, variables, docstrings et chemins.
  2. *Critère 2 (Conformité théorique & Zéro trick)* : Identifier les zones de conformité exacte aux Piliers I à V et débusquer les éventuels artefacts historiques (seuils arbitraires, branches mortes, contournements temporaires).
  3. *Critère 3 (Langue anglaise académique)* : Évaluer la qualité linguistique (docstrings en français à traduire, nomenclature hybride).
  4. *Critère 4 (Modularité)* : Identifier les couplages trop forts (ex: apex liant coordonnées et évaluation d'arbres, propagateurs monolithiques).
  5. *Critère 5 (Interchangeabilité des intégrateurs)* : Évaluer la faisabilité d'une signature standard unique pour `Strang S2`, `Blanes-Moan BM4-7`, `Blanes-Moan BM6-11`, et les correcteurs $\mathcal{C}_4 / \mathcal{C}_6$.
  6. *Critère 6 (Jeux de tests et benchmarks préconfigurés)* : Classifier les scénarios physiques existants (Burrau, Sitnikov, Figure-8, Choc radial, Quadruple, Hill, Laplace, JPL DE441) pour les transformer en fixtures réutilisables.

### Step 3: Spécification Architecturale de la Bibliothèque Scientifique Cible ("lib")
- **Context**: La bibliothèque doit être structurée pour publication académique et utilisable comme package autonome.
- **Objective**: Formaliser le blueprint complet de la future bibliothèque (nom de package suggéré : `conormal_nbody` ou `pyconormal`), son arborescence modulaire, ses contrats d'interface et ses protocoles de simulation.
- **Deliverables**:
  - `TargetLibrarySpecification`: Architecture cible, abstractions et signatures de classes documentées dans l'index du brain.
- **Process**: Définition des modules de la bibliothèque (`core/`, `integrators/`, `clusters/`, `scheduling/`, `benchmarks/`, `data/`).

### Step 4: Élaboration et Restitution du Rapport de Synthèse dans la Vue
- **Context**: L'utilisateur requiert un document de rapport accessible dans sa vue de travail.
- **Objective**: Rédiger une synthèse exécutive de l'analyse, exposant le diagnostic, la matrice des écarts, la structure de la future bibliothèque et la feuille de route de mise en œuvre.
- **Deliverables**:
  - `ExecutiveReport`: Rapport de synthèse didactique et structuré restitué dans la réponse conversationnelle.
- **Process**: Consolidation des résultats, confrontation des compromis architecturaux et proposition des étapes suivantes.

## Acceptance Criteria:
### Step 1:
- Recensement documenté de l'ensemble des modules Python de `cas6_nbody3d_async` sans exception.
- Graphe de dépendance des modules complètement explicité.
### Step 2:
- Document `codebase_analysis_and_index.md` rédigé dans le brain, exhaustif, couvrant chaque module selon les 6 axes d'évaluation.
- Zéro écriture dans `redaction/unsoftened_conormal_nbody3d/`.
### Step 3:
- Spécification architecturale claire proposant des interfaces d'intégrateurs unifiées (`IntegratorProtocol`) et un chargeur de cas tests prédéfinis.
### Step 4:
- Rapport exécutif clair, rigoureux, sans jargon vague, permettant à l'utilisateur de valider la trajectoire de refactorisation.

## Observations:
> [!IMPORTANT]
> Cette itération est une itération d'analyse approfondie, d'indexation exhaustive et de cadrage architectural. Aucune modification de code destructive ne sera effectuée sur le moteur opérationnel actuel pendant cette phase, garantissant l'intégrité absolue des 159 tests unitaires passants et du sanctuaire de rédaction.
