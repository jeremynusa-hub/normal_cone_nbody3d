# Implementation Plan of Iteration 124 : Multi-Audit Sous-Agentique du Flot d'Intégration, de la Complétude des Tests et de l'Isomorphisme Théorique
Created at: 2026-10-07T19:04:30+02:00

## Context:
1. **Explicit Request**: « Y a encore des fails, et le solveur semble pas tout tester correctement, lance des audit sous agentique pour comprendre le flot du code et le comparer au papier de recherche et au code d'avant »
2. **Implicit Strategic Value**: Assurer que la bibliothèque [`normal_cone_nbody_float64`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/) ne masque aucune erreur de conception ou de discrétisation, que le pipeline de calcul (coordonnées, intégrateurs, super-nœuds, transitions d'atlas, ordonnanceur dyadique) reflète fidèlement la théorie de la déformation conormale des Chapitres 3, 4 et 5, et que la matrice de tests et de benchmarks évalue rigoureusement toutes les composantes sans faux positifs ni tests tronqués.
3. **Alignment Criteria**: Production d'un audit sous-agentique contradictoire à 3 volets étanches fournissant une cartographie précise du flot d'exécution, une comparaison ligne à ligne avec le code hérité, une validation théorique avec le manuscrit LaTeX, et un diagnostic sans complaisance des zones de tests manquantes ou des points de défaillance.

## Guidelines:
- **Méthode** : Déploiement de trois sous-agents spécialisés et étanches via `invoke_subagent` (Flot & Architecture, Alignement Monographie LaTeX, Héritage & Matrice de Vérification).
- **Lignes Rouges** : Ne modifier aucun fichier source du manuscrit `redaction/unsoftened_conormal_nbody3d/` ; conserver le zéro adoucissement strict ($\varepsilon = 0, d_0 = 0$) partout.
- **Rigueur Analytique** : Analyser les mécanismes réels d'exécution et non les déclarations nominales (traçage des appels, activation effective des super-nœuds, switches d'arbres pendant les sous-pas, préservation symplectique).

## Plan:
### Step 1: Définition et Déploiement du Triple Audit Sous-Agentique
- **Matériel de départ** : Code source complet de [`normal_cone_nbody_float64`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/), manuscrit LaTeX dans `redaction/unsoftened_conormal_nbody3d/`, et code hérité dans `models/normal_cone_integrator/cas6_nbody3d_async/`.
- **Objectif** : Lancer trois sous-agents avec des mandats exhaustifs, précis et complémentaires.
- **Livrable** : `Subagent Deployments` (Flot & Architecture, Théorie LaTeX, Code Hérité & Tests).
- **Processus** : Configuration des prompts d'initialisation et lancement concurrent via `invoke_subagent`.

### Step 2: Synthèse Contradictoire et Cartographie du Flot d'Exécution Réel
- **Matériel de départ** : Rapports intermédiaires des 3 sous-agents.
- **Objectif** : Confronter les conclusions des auditeurs sur le flot effectif (`NBodySystem` -> `DyadicScheduler` -> `Integrator` -> `Jacobi` / `Supernodes` / `Apex`).
- **Livrable** : `Audit Synthesis Report` documentant les déviations structurelles, les simplifications abusives et les divergences d'implémentation.
- **Processus** : Compilation critique et mise en exergue des points névralgiques où des défaillances ou des absences de test ont été identifiées.

### Step 3: Plan d'Action Correctif et Validation des Solveurs
- **Matériel de départ** : Rapport de synthèse de l'audit.
- **Objectif** : Établir la feuille de route des correctifs nécessaires pour corriger les failles du solveur et renforcer la couverture des tests.
- **Livrable** : `Actionable Remediation Roadmap`.
- **Processus** : Priorisation P0/P1 des anomalies de flot et des tests à implémenter.

## Acceptance Criteria:
- **Critère 1** : Trois sous-agents déployés avec des périmètres étanches couvrant l'intégralité du cycle de vie du code.
- **Critère 2** : Élucidation complète du chemin d'exécution réel (flot des données, conversions Jacobi/Cartésien, gestion des super-nœuds, switches d'atlas).
- **Critère 3** : Identification sans ambiguïté des causes des « fails » et des lacunes de couverture dans les tests unitaires et les benchmarks.
- **Critère 4** : Rapport de synthèse consolidé accessible dans le brain.

## Observations:
> [!IMPORTANT]
> L'audit doit porter une attention critique sur le couplage entre l'adaptation d'atlas sous hystérésis $\gamma=0.6$, les intégrateurs multi-étages ($\mathrm{BM}_4$, $\mathrm{BM}_6$) et l'ordonnancement dyadique : le changement d'arbre de Jacobi est-il autorisé au sein des sous-pas ou uniquement à la frontière macroscopique ?
> [!WARNING]
> La régularisation des super-nœuds doit être auscultée : le solveur synchrone $k$-corps est-il véritablement appelé et résolu lors des collisions multiples, ou s'agit-il d'un passage en drift découplé ?
