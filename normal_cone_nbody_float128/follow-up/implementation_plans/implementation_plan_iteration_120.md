# Implementation Plan of Iteration 120 : Audit Sous-Agentique Multi-Blocs de la Bibliothèque (Conformité aux Citations et Exactitude Mathématique)
Created at: 2026-10-07T18:21:00+02:00

## Context:
1. **Explicit Request**: Lancer un audit sous-agentique avec un sous-agent par bloc fonctionnel (au moins un par solveur, un par intégrateur et un par dossier conteneur). Chaque sous-agent doit vérifier que la section du papier de recherche est bien mentionnée dans l'en-tête du code, correspond au code développé, et que les formulations mathématiques sont rigoureusement exactes et cohérentes avec la monographie LaTeX.
2. **Implicit Strategic Value**: Établir une garantie académique et éditoriale absolue de premier ordre avant soumission ou diffusion publique : attester que chaque module de `normal_cone_nbody_float64` est mathématiquement isomorphe à sa théorie source, sans simplification abusive, sans glissement de signe, et avec une traçabilité documentaire parfaite vers les chapitres de la monographie.
3. **Alignment Criteria**: Cartographie exhaustive des 13 blocs fonctionnels audités par des sous-agents autonomes, rapports d'audit individuels détaillés, identification d'éventuelles divergences ou manques de citation, et synthèse consolidée avec matrice de conformité mathématique.

## Guidelines:
- **Method**: Déploiement concurrent de sous-agents de type `research` équipés d'outils de lecture en lecture seule (`view_file`), garantissant le sanctuarisation absolue du manuscrit LaTeX `redaction/unsoftened_conormal_nbody3d/` et du code source pendant l'audit.
- **Red Lines**: Interdiction stricte de modifier les fichiers sources ou les fichiers `.tex` pendant l'audit. Aucun usage de commandes shell pour altérer l'espace de travail.
- **Architecture de l'Audit**:
  - **Pôle Dossiers Conteneurs (5 agents)** :
    1. `core/` (Algèbre symplectique, invariants de Lie--Bargmann, état étendu, coordonnées de Jacobi, transitions $\mathrm{SL}(2, \mathbb{R})$).
    2. `potential/` (Crible conormal $\Theta(u_v)$, apex $\mathcal{C}^{1,1}$, détecteur de marée de Hill, Newton direct).
    3. `system/` (Orchestrateur `NBodySystem`, `AtlasManager` à hystérésis $\gamma=0.6$, `DyadicScheduler`, `HistoryTracker`).
    4. `diagnostics/` (Radar de collision, les 6 KPIs célestes, visualiseur 10 panneaux).
    5. `benchmarks/` (Registre des 8 cas canoniques, runner de simulation, analyseur Pareto travail-précision).
  - **Pôle Solveurs Individuels (4 agents)** :
    6. `solvers/stumpff.py` (Séries de Stumpff $c_0-c_5$, récurrences différentielles, formes analytiques).
    7. `solvers/kepler_sundman.py` (Solveur universel Kepler--Sundman, Halley cubique, régularisation sans adoucissement).
    8. `solvers/mcgehee.py` (Involution de McGehee $\mathcal{B}_C^2 = \operatorname{Id}_{\mathbf{M}_C}$, isométrie cinétique).
    9. `solvers/supernode.py` (Super-nœuds synchrones, horloge dépliée $\mathrm{d}t_S = \chi_S^2 \mathrm{d}\tau_S$, anti-Zénon).
  - **Pôle Intégrateurs Individuels (4 agents)** :
    10. `integrators/strang.py` (Schéma palindromique $\mathcal{S}_2$ d'ordre 2, réversibilité, préservation de Casimir).
    11. `integrators/yoshida.py` (Schéma de composition $\mathcal{S}_4$ d'ordre 4, coefficients de Suzuki-Yoshida).
    12. `integrators/blanes_moan_4.py` & `integrators/blanes_moan_6.py` (Compositions optimales $\mathrm{BM}_4^{(7)}$ et $\mathrm{BM}_6^{(11)}$, forme FSAL, sous-pas rétrogrades bornés, constantes d'ombre).
    13. `integrators/correctors.py` & `integrators/processed.py` (Correcteurs $\mathcal{C}_{2p}$, signe inversé $\exp(+h^{2p} X_K)$, intégrateurs traités $\mathcal{P}\mathrm{BM}$, télescopage dyadique).
- **Form**: Prompts d'audit unifiés avec protocole d'évaluation en 3 volets (Citations, Isomorphisme mathématique, Recommandations). Rapport de synthèse rédigé en français avec liens markdown cliquables vers les fichiers audités.

## Plan:
### Step 1: Cadrage et Définition des Missions d'Audit
- **Starting Material**: Structure de la bibliothèque `normal_cone_nbody_float64/` et manuscrit LaTeX dans `redaction/unsoftened_conormal_nbody3d/`.
- **Objective**: Structurer les mandats précis et les assignations de fichiers/sections pour les 13 sous-agents d'audit.
- **Deliverables**:
  - `AuditMissionMatrix`: Grille d'assignation explicite associant chaque sous-agent à ses fichiers cibles de code et ses chapitres/sections LaTeX de référence.
- **Process**: Cartographier les dépendances documentaires et valider les périmètres d'intervention.

### Step 2: Lancement Parallèle des 13 Sous-Agents d'Audit
- **Starting Material**: Outil `invoke_subagent` avec TypeName `research` et modèle `inherit`.
- **Objective**: Déployer simultanément les 13 agents d'audit indépendants pour examiner le code et le texte théorique.
- **Deliverables**:
  - `SubagentsDeployment`: 13 instances de sous-agents en exécution d'audit.
- **Process**: Invoquer les agents avec des instructions strictes (lecture via `view_file`, zéro modification, analyse critique comparative).

### Step 3: Collecte et Analyse des Rapports d'Audit Individuels
- **Starting Material**: Messages de restitution retournés par les 13 sous-agents.
- **Objective**: Recueillir les conclusions de chaque agent, extraire les verdicts de citation et d'exactitude mathématique, et recenser les écarts éventuels.
- **Deliverables**:
  - `IndividualAuditReports`: Corpus de 13 évaluations détaillées.
- **Process**: Synthétiser les constats, confronter les appréciations sur les formules mathématiques et relever les citations manquantes ou imprécises.

### Step 4: Rédaction du Rapport Global de Certification et Matrice de Synthèse
- **Starting Material**: Corpus des rapports d'audit individuels.
- **Objective**: Rédiger un artefact complet de synthèse d'audit `multi_agent_mathematical_and_citation_audit.md` avec matrice de conformité globale.
- **Deliverables**:
  - `AuditSynthesisArtifact`: Document de synthèse d'audit certifiant la conformité du code au papier.
- **Process**: Consolidation structurée avec tableau récapitulatif des 13 blocs, citations validées et plan de retouches documentaires éventuelles.

## Acceptance Criteria:
### Step 1:
- Les 13 sous-agents couvrent 100 % des fichiers de la bibliothèque logicielle (dossiers conteneurs, solveurs et intégrateurs).
### Step 2:
- Chaque sous-agent reçoit un mandat précis avec les chemins absolus du code à auditer et des chapitres LaTeX associés.
### Step 3:
- Les 13 sous-agents retournent un verdict explicite sur (1) la présence et pertinence des citations, et (2) l'exactitude des équations mathématiques.
### Step 4:
- L'artefact de synthèse consolidé offre une traçabilité intégrale et met en évidence toute discordance théorique ou lacune de documentation.

## Observations:
> [!IMPORTANT]
> Les sous-agents opèrent en mode strictement lecture seule (`research`) afin de préserver l'intégrité absolue de la monographie LaTeX et du code source.
> [!NOTE]
> La décomposition en 13 sous-agents garantit une granularité maximale : 5 dossiers conteneurs globaux, 4 solveurs élémentaires et 4 briques d'intégration symplectique.
