# Implementation Plan of Iteration 111 : Conception de l'Architecture Logicielle Académique Conforme à l'Article et Spécification de l'Objet d'État Dynamique Unifié
Created at: 2026-10-07T16:38:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur demande de se positionner comme l'ingénieur chercheur chargé de réaliser le code à intégrer et présenter dans l'article [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/). Pour cela, il est impératif de respecter les mathématiques rédigées dans l'article et les conventions de nommage, en s'appuyant sur `AGENT_GUIDE.md` et `conormal_vs_normal_cone_foundations.md`. L'architecture doit rigoureusement séparer : (1) les éléments classiques élémentaires potentiellement optimisables, (2) les composants modulaires (solveurs/intégrateurs), (3) les tests unitaires et scénarios, (4) les outils de diagnostique, et (5) définir un objet central permettant de stocker l'état instantané du système, d'exécuter les solveurs dessus, tout en conservant l'historique au besoin.
2. **Implicit Strategic Value**: Assurer une parfaite cohérence épistémique et formelle entre le texte mathématique de l'article de recherche et la bibliothèque logicielle d'accompagnement. Cet alignement est indispensable pour garantir l'acceptation par les pairs dans une revue académique de premier plan et offrir à la communauté une bibliothèque de mécanique céleste géométrique propre, modulaire et extensible.
3. **Alignment Criteria**: L'architecture proposée devra refléter fidèlement la taxonomie de Fulton (1984), Verdier (1976) et Calaque-Safronov (2024), bannir formellement les néologismes proscrits ("conormal deformation", "conormal cone"), isoler le calcul pur dans des noyaux optimisables sans overhead d'allocation, et doter l'utilisateur d'une classe d'état `ConormalSystem` claire, élégante et unifiée.

## Guidelines:
- **Method**: Démarche d'ingénierie mathématique descendante (*Theory-to-Code Pipeline*) : chaque classe et fonction Python découle directement d'une définition, équation ou théorème de l'article de recherche.
- **Red Lines**:
  - Sanctuarisation absolue du répertoire [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/) : interdiction formelle d'y écrire ou de le modifier.
  - Zéro Clamping inviolable ($\varepsilon = 0, d_0 = 0$).
  - Respect strict des proscriptions sémantiques de `AGENT_GUIDE.md` (interdiction absolue de "conormal deformation", de "conormal cone", et de localisation erronée du kick).
  - Aucune modification destructive prématurée du code avant validation complète du design.
- **Architecture**: Architecture en 4 couches concentriques étanches (Noyau élémentaire $\to$ Algorithmes & Solveurs $\to$ Système & Orchestration $\to$ Diagnostiques & Scénarios), orientée objet moderne avec protocoles typés stricts (`typing.Protocol`, `dataclasses`, `numpy.typing`).

## Plan:
### Step 1: Cadrage Épistémique & Alignement Terminologique avec l'Article
- **Context**: L'article `unsoftened_conormal_nbody3d` et le document de doctrine `conormal_vs_normal_cone_foundations.md` fixent la vérité mathématique du projet.
- **Objective**: Établir la nomenclature exacte, les symboles et les définitions formelles reliant le texte LaTeX de l'article à l'implémentation Python.
- **Deliverables**:
  - `EpistemicGlossary`: Matrice de correspondance concepts mathématiques $\leftrightarrow$ identifiants de code Python.
- **Process**: Extraction méthodique des symboles des chapitres 2 à 6 de l'article et validation de leur conformité aux règles d'`AGENT_GUIDE.md`.

### Step 2: Découpage Modulaire en 4 Couches Étanches
- **Context**: Le code doit séparer nettement l'élémentaire optimisable, le modulaire interchangeable, les scénarios et les outils d'instrumentation.
- **Objective**: Structurer le package cible en 4 couches rigoureuses :
  - *Couche 1 (Socle Classique & Optimisable)* : algèbre vectorielle 3D, états cartésiens, cométrique cinétique de Jacobi, séries de Stumpff, solveur Sundman-Halley, rebond de McGehee.
  - *Couche 2 (Moteurs Modulaires / Solveurs)* : coupure conormale lisse $\Theta(u_v)$, maintenance d'atlas $\mathrm{SL}(2, \mathbb{R})$, détection adimensionnée des super-nœuds, famille d'intégrateurs symplectiques ($\mathcal{S}_2, \mathrm{BM}_4, \mathrm{BM}_6, \mathcal{C}_4, \mathcal{C}_6$).
  - *Couche 3 (Diagnostique & Métrologie)* : invariants de Lie-Bargmann, radar d'anticipation d'événements, projection sur les écrans $\mathbf{Q}[N]$, visualisations 10 panneaux 300 DPI.
  - *Couche 4 (Scénarios & Bancs d'Essais)* : fixtures déclaratives pour les 11 cas célèbres (Burrau, Figure-8, Sitnikov, Quadruple, Hill, Laplace, DE441).
- **Deliverables**:
  - `LayeredArchitectureBlueprint`: Schéma fonctionnel et découpage modulaire validé.
- **Process**: Spécification des dépendances autorisées (unidirectionnelles) interdisant tout couplage ascendant ou circulaire.

### Step 3: Spécification Formelle de l'Objet Maître d'État et de Simulation (`ConormalSystem`)
- **Context**: L'utilisateur demande de définir un "objet" central stockant l'état instantané, permettant d'exécuter les solveurs dessus et de conserver l'historique si nécessaire.
- **Objective**: Concevoir la classe unifiée `ConormalSystem` regroupant :
  - État instantané canonique : $t, \mathbf{q}, \mathbf{p}, \mathbf{m}$, topologie active $\mathcal{T}$, super-nœuds actifs, énergie initiale $E_0$.
  - Interface d'exécution fluide : `system.step(dt, integrator="bm4_7")`, `system.evolve(t_span, dt, ...)`.
  - Gestionnaire d'historique détachable (`HistoryTracker`) : enregistrement optionnel (`record_history=True`) sans pénalité de performance en mode croisière (zéro copie si inactif).
  - Capacités de rejeu et d'inversion temporelle microscopique ($\mathcal{R}_{\text{ext}}$).
- **Deliverables**:
  - `ConormalSystemSpecification`: Définition complète de l'API de `ConormalSystem` et de ses protocoles associés.
- **Process**: Conception d'une interface élégante, intuitive et hautement performante.

### Step 4: Élaboration du Rapport d'Ingénierie Chercheur et Restitution dans la Vue
- **Context**: L'utilisateur attend la posture et la restitution détaillée de l'ingénieur chercheur.
- **Objective**: Rédiger une synthèse exhaustive exposant les choix de conception, la justification mathématique au regard de l'article, les signatures d'API et les diagrammes d'architecture.
- **Deliverables**:
  - `ResearchEngineerReport`: Restitution complète et structurée dans la réponse conversationnelle.
- **Process**: Synthèse didactique, rigoureuse et argumentée pour validation par l'utilisateur.

## Acceptance Criteria:
### Step 1:
- Zéro occurrence de termes proscrits ("conormal deformation", "conormal cone", "cas 6") dans la conception.
- Alignement strict avec les théorèmes et notations des chapitres 2 à 6 de l'article.
### Step 2:
- Séparation explicite et étanche des 4 couches avec graphe de dépendance unidirectionnel.
- Noyau élémentaire (Couche 1) conçu pour être optimisable (vectorisable NumPy / C-ABI / Numba sans refonte d'interface).
### Step 3:
- Spécification formelle de `ConormalSystem` intégrant l'état instantané, l'exécution polymorphe de solveurs et la gestion optionnelle d'historique.
### Step 4:
- Rapport d'ingénierie chercheur limpide, exhaustif et directement actionnable pour l'implémentation.

## Observations:
> [!IMPORTANT]
> Conformément à la règle de sanctuarisation, le répertoire `redaction/unsoftened_conormal_nbody3d/` est consulté exclusivement en lecture seule. Tous les développements logiciels seront logés dans le nouveau package `conormal_nbody/`.
