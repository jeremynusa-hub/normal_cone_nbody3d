# Implementation Plan of Iteration 135 : Audit Contradictoire Multi-Agent et Télémétrie Fine du Système 4-Corps Hiérarchique
Created at: 2026-10-07T22:20:00+02:00

## Context:
1. **Explicit Request**:
   - « Sur le Système 4-Corps Hiérarchique on voit encore des saut d'erreur d'energie, identifie comment découper l'analyse du code pour lancer un audit de sous agents en réfléchissant a la granularité et les thématique au vu de la source of truth qu'est l'article, mais qu'on peut aussi lire les anciens couple de document ou l'ancien code »
   - Décision validée via consultation : Lancement de l'audit complet à 4 sous-agents spécialisés avec script de télémétrie pas-à-pas pour isoler l'instant exact du saut.
2. **Implicit Strategic Value**:
   - Élucider la cause racine du résidu d'énergie ($1.96 \times 10^{-7}$) et des sauts de courbure observés sur le 4-corps hiérarchique, pour élever ce scénario multi-échelles complexe au même niveau de pureté machine que Burrau ($4.49 \times 10^{-11}$) et Figure-8 ($2.12 \times 10^{-12}$).
   - Confronter de manière étanche et contradictoire la géométrie de Jacobi à 4 corps, les boîtes de Rees, les correcteurs symplectiques sous switches d'atlas, et les leçons tirées des anciens "kicks virtuels" du solveur CAS6.
3. **Alignment Criteria**:
   - Source of truth primaire : Manuscrit LaTeX (`redaction/unsoftened_conormal_nbody3d/`).
   - Sources secondaires : Les 9 couples formels Markdown (`formal/Formal/DeformationToTheNormalCone/Cas6_NBody3D/docs/`) et le code historique (`models/normal_cone_integrator/cas6_nbody3d_async/`).
   - Télémétrie pas-à-pas avec horodatage, distances mutuelles, topologie d'arbre active, et variation symplectique locale.

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des scénarios.
- **Règle d'Édition** : Utiliser exclusivement `replace_file_content` pour modifier le code existant.

## Plan:
### Step 1: Déploiement du Script de Télémétrie Fine & Profilage Pas-à-Pas
- **Contexte & Point de départ** : Le système à 4 corps est simulé sur $t \in [0, 8.0]$ avec $dt = 0.01$ dans `scratch/run_showcase_simulations.py`.
- **Objectif** : Isoler l'instant précis $t_{\text{jump}}$, la topologie de l'arbre $\mathcal{T}$, les distances mutuelles $r_{ij}$ des 6 paires, le niveau dyadique $k$ et la composante d'opérateur incriminée (Drift, Kepler, Apex Kick, ou Corrector).
- **Livrable** : `Component Name`: `TelemetryProbe` dans `scratch/diagnose_4body_jump.py` produisant un rapport télémétrique détaillé.

### Step 2: Lancement Concurrent des 4 Sous-Agents d'Audit
- **Contexte & Point de départ** : 4 axes d'analyse orthogonaux définis par le cadrage systémique.
- **Objectif** : Déployer 4 sous-agents spécialisés via `invoke_subagent` :
  - **Sous-Agent 1 (Géométrie de Jacobi & Transitions d'Atlas à $N=4$)** :
    - *Périmètre :* `core/jacobi.py`, `core/coordinates.py`, `core/transitions.py`, `system/atlas_manager.py`.
    - *Mandat :* Audit de la topologie peigne vs double-paire binaire $((0,1),(2,3))$, des projecteurs symplectiques dans $\mathbb{R}^{18}$, de la complétude de $\mathrm{SL}(2, \mathbb{R})$, et du comportement lors des 15 transitions d'atlas.
  - **Sous-Agent 2 (Potentiel d'Apex Écranté, Boîtes de Rees & Kicks Virtuels)** :
    - *Périmètre :* `potential/screened_apex.py`, `potential/screening.py`, `potential/direct.py`, `integrators/common.py`.
    - *Mandat :* Vérification de la partition exacte de la force $F_{\text{Newton}} = F_{\text{leaf}} + F_{\text{super}} + F_{\text{apex}}$, de la distribution des gradients sur les 4 corps, et de l'absence de force virtuelle ou de fuite monopolaire lors des frôlements entre corps de sous-nœuds composites.
  - **Sous-Agent 3 (Correcteurs Symplectiques, Algèbre de Lie & Ordonnancement Dyadique)** :
    - *Périmètre :* `integrators/correctors.py`, `integrators/processed.py`, `integrators/blanes_moan_6.py`, `system/dyadic_scheduler.py`, `system/nbody_system.py`.
    - *Mandat :* Audit du télescopage frontière $\mathcal{C}_A^{-1} \circ \mathcal{C}_B$ sous transitions de cartes, de la consistance de Skeel-Gear lors d'approches binaires simultanées, et de la réversibilité temporelle sous pas négatifs.
  - **Sous-Agent 4 (Archéologie de l'Ancien Code & Précédents Historiques CAS6)** :
    - *Périmètre :* `models/normal_cone_integrator/cas6_nbody3d_async/` (code et dossiers `follow-up/walkthroughs/`).
    - *Mandat :* Exhumation des défauts historiques identiques résolus aux itérations 63-65 (kicks virtuels d'apex) et 108 (synchronisation multirate), et confrontation avec la structure actuelle.

### Step 3: Synthèse Contradictoire et Diagnostic Formel
- **Contexte & Point de départ** : Rapports des 4 sous-agents et données de la télémétrie.
- **Objectif** : Réunir les conclusions, croiser les données empiriques avec les analyses de code et la théorie du manuscrit, et isoler formellement les causes exactes des sauts d'énergie.
- **Livrable** : Document de synthèse dans le brain [`diagnostic_sauts_energie_4corps.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/diagnostic_sauts_energie_4corps.md).

### Step 4: Élaboration de la Feuille de Route de Remédiation & Walkthrough
- **Contexte & Point de départ** : Diagnostic formel certifié.
- **Objectif** : Formaliser le plan d'action correctif et consigner les résultats dans `notes.md` et `walkthrough.md`.

## Acceptance Criteria:
- **Critère 1** : Le script de télémétrie localise le pas précis et les grandeurs physiques associées au premier saut d'énergie significatif.
- **Critère 2** : Les 4 sous-agents rendent des rapports d'audit complets, documentés avec citations d'équations et de code.
- **Critère 3** : Le document de synthèse identifie formellement la cause racine du saut d'énergie sur le 4-corps.
