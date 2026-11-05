# Implementation Plan of Iteration 112 : Planification Stratégique de la Roadmap V et Découpage Granulaire de la Bibliothèque normal_cone_nbody_float64
Created at: 2026-10-07T17:30:30+02:00

## Context:
1. **Explicit Request**: L'utilisateur demande de planifier la roadmap de création de la bibliothèque en découpant le projet en plusieurs itérations granulaires. Pour chaque itération, le processus doit STRICTEMENT suivre ce protocole :
   - D'abord identifier et relire les sections théoriques du papier.
   - Puis relire l'ancien code associé et identifier la conformité ou les éléments manquants.
   - Puis établir un document de spécification claire du code à écrire avec les tests et attendus.
   - Puis rédiger ce code (modulaire, décomposé en fonctions pures au rôle limpide, utilisant des abstractions adaptées tout en conservant une vision d'optimisation) en commençant chaque fichier par un commentaire descriptif référençant les sections de l'article de recherche.
   - Puis tester ce code.
2. **Implicit Strategic Value**: Établir une trajectoire d'ingénierie logicielle irréprochable qui garantit l'alignement absolu entre la théorie de l'article, les preuves formelles en Lean 4 et l'implémentation Python, tout en prévenant les dérives architecturales, les régressions et les surcoûts d'allocation mémoire.
3. **Alignment Criteria**: La feuille de route stratégique `roadmap.md` doit partitionner la bibliothèque en 7 itérations opérationnelles calquées sur les 7 couches de l'architecture maîtresse, imposer le protocole à 5 étapes comme contrat d'itération inviolable, et recueillir la validation formelle de l'utilisateur avant d'engager le code de la Couche 1.

## Guidelines:
- **Method**: Démarche ascendante par couches étanches (*Layered Theory-to-Code Pipeline*) régie par le protocole strict à 5 étapes : (1) Relecture théorique $\to$ (2) Audit du code hérité $\to$ (3) Spécification formelle $\to$ (4) Implémentation modulaire et optimisée $\to$ (5) Tests unitaires miroirs (1:1).
- **Red Lines**:
  - Sanctuarisation absolue du répertoire [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/) : interdiction formelle d'y écrire.
  - Zéro adoucissement ni troncature ($\varepsilon = 0, d_0 = 0$).
  - Broutement topologique strictement interdit via le bandeau d'hystérésis asymétrique $\gamma = 0.6$.
  - Interdiction d'utiliser le terme `"cas6"` dans le nouveau package.
- **Architecture**: 7 couches étanches (Core, Solvers, Potential, Integrators, System, Diagnostics, Benchmarks) avec protocoles typés stricts (`typing.Protocol`), dataclasses immuables et calculs vectorisés purs.
- **Form**: Rédaction des artefacts de planification en français, en Markdown natif sans balises HTML ni styles inline.

## Plan:
### Step 1: Formalisation de la Roadmap V Stratégique et Découpage des 7 Itérations
- **Context**: Le document maître finalisé [`normal_cone_nbody_float64_master_architecture.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/normal_cone_nbody_float64_master_architecture.md) spécifie 28 modules répartis en 7 couches, 28 modules de tests et 8 bancs d'essais canoniques.
- **Objective**: Établir dans [`roadmap.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/roadmap.md) le cadencement séquentiel des itérations 113 à 119 en y intégrant le protocole obligatoire à 5 étapes.
- **Deliverables**:
  - `StrategicRoadmapV`: Fichier [`roadmap.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/roadmap.md) structuré avec objectifs macro, continuités et critères de passage (gating) pour chaque itération.
- **Process**: Découpage biunivoque couche par couche, avec identification systématique des chapitres du papier et des sous-répertoires de l'ancien code.

### Step 2: Cadrage Opérationnel de l'Itération 112 et Alignement des Contrats de Tâche
- **Context**: L'Itération 112 formalise la trajectoire globale et verrouille le protocole d'ingénierie avant la première ligne de code.
- **Objective**: Synchroniser les artefacts de l'itération en cours (`task.md` et `implementation_plan.md`) en conformité stricte avec les règles du skill `planning-contracts`.
- **Deliverables**:
  - `Iteration112ImplementationPlan`: Fichier [`implementation_plan.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/implementation_plan.md).
  - `Iteration112TaskList`: Fichier [`task.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/task.md) reflétant rigoureusement les étapes d'exploration, de réalisation et de vérification.
- **Process**: Archivage préalable des artefacts d'Itération 111 (`*_iteration_111.md`) et instanciation des nouveaux fichiers.

### Step 3: Préparation du Protocole d'Exécution de l'Itération 113 (Couche 1 `core/`)
- **Context**: L'Itération 113 constituera la première étape de développement concret sur le socle `core/`.
- **Objective**: Identifier avec précision les 6 modules de la Couche 1, leurs correspondances théoriques (Chapitres 1 à 3), leur code hérité à assainir et le plan de test unitaire miroir associé.
- **Deliverables**:
  - `Layer1ExecutionBrief`: Fiche de cadrage opérationnel pour la Couche 1 intégrée dans la restitution utilisateur.
- **Process**: Synthèse croisée entre le rapport d'audit du code hérité et les théorèmes Lean 4 associés aux invariants et à la géométrie de Jacobi.

## Acceptance Criteria:
### Step 1:
- Fichier [`roadmap.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/roadmap.md) créé et contenant la stratégie macro, les 7 itérations d'implémentation (113 à 119) et la checklist multi-itérations.
- Chaque bloc d'itération détaille explicitement le protocole à 5 étapes (théorie, code hérité, spécification, code modulaire optimisé avec en-têtes, tests unitaires).
### Step 2:
- Parfaite équivalence structurelle 1-to-1 entre `implementation_plan.md` et la section `Realisation Phase` de `task.md`.
- Statuts de progression respectant exclusivement la syntaxe contractuelle (`[ ]`, `[/]`, `[x]`).
### Step 3:
- Restitution claire à l'utilisateur présentant la feuille de route et sollicitant son approbation formelle avant le déclenchement de l'Itération 113.

## Observations:
> [!IMPORTANT]
> Conformément à l'exigence de l'utilisateur, chaque module Python créé au cours des itérations 113 à 119 débutera obligatoirement par un en-tête descriptif exhaustif indiquant :
> 1. Son rôle fonctionnel précis et son périmètre mathématique.
> 2. Les références explicites aux chapitres, sections, propositions et théorèmes du papier de recherche [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/).
> 3. Les théorèmes certifiés correspondants en Lean 4.
> 4. Les garanties d'invariance et de performance (zéro adoucissement $\varepsilon = 0$, conservation symplectique $< 10^{-14}$, zéro allocation superflue dans les boucles intérieures).
