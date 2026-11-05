# Implementation Plan of Iteration 143 : Audit Contradictoire Multi-Agentique Exhaustif des Chapitres 3, 4 et 5 (9 Auditeurs)
Created at: 2026-10-08T09:40:00+02:00

## Context:
1. **Explicit Request**:
   - L'utilisateur constate la persistance de sauts d'énergie et exige : "Lance un audit sous agentique multiple, trois agent par chapitre du papier de recherche : 3,4,5. Sur chaque chapitre, chaque agent devra : relire exhaustivement le chapitre adressé, identifier les étapes du code qui dépendent de ce chapitre, faire un audit impitoyable de ce code. Il y aura trois agents par chapitre car un doit être expert ingénieur de recherche de la traduction des maths en code, l'un doit être expert numérique des schéma d'intégration d'équation, le dernier est expert en algèbre de lie et homologie avec application numérique".
2. **Implicit Strategic Value**:
   - Révéler les éventuelles approximations masquées, les ruptures d'isomorphisme mathématique, les pertes d'ordre symplectique ou les discontinuités de jauge au sein du solveur conormal, en croisant trois perspectives épistémiques indépendantes et étanches pour chaque chapitre théorique.
3. **Alignment Criteria**:
   - Déploiement de 9 sous-agents d'audit autonomes répartis sur 3 triades disciplinaires (Chapitres 3, 4, 5).
   - Lecture exhaustive des sources LaTeX originales (`redaction/unsoftened_conormal_nbody3d/chapters/ch03`, `ch04`, `ch05`).
   - Cartographie exhaustive du code Python dépendant de chaque chapitre.
   - Restitution d'un rapport contradictoire consolidé (`audit_theorique_chapitres_3_4_5.md`) classant chaque point par degré de criticité (Bloquant / Déviation / Conforme).

## Guidelines:
- **Method**:
  * Découpage en 3 triades étanches :
    - **Axe 1 (Ingénieur Recherche Math-to-Code)** : Traduction des équations, lemmes et théorèmes du papier en structures logiques du code.
    - **Axe 2 (Analyste Numérique)** : Conditionnement, stabilité sous pas négatifs/hyperboliques, annulation soustractive, convergence effective.
    - **Axe 3 (Expert Lie & Homologie)** : Préservation symplectique de Darboux, crochets de Poisson, équations homologiques, commutateurs et télescopage.
- **Red Lines**:
  * Aucun code source Python ne doit être modifié durant cette itération d'audit pur.
  * Interdiction d'éditer le manuscrit LaTeX de recherche.
  * Autonomie totale des sous-agents avec restitution de preuves textuelles et de citations d'équations.
- **Architecture**:
  * Manuscrit LaTeX : `redaction/unsoftened_conormal_nbody3d/chapters/`
  * Bibliothèque cible : `normal_cone_nbody_float64/`
  * Rapport de synthèse : `audit_theorique_chapitres_3_4_5.md` dans le brain.

## Plan:
### Step 1: Déploiement Triadique sur le Chapitre 3 (Paradigme 3-Corps & Régularisations)
- **Context**: Le Chapitre 3 traite du problème à 3 corps, des coordonnées conormales, du difféomorphisme de Sundman, des fonctions universelles de Stumpff et de l'opérateur de rebond de McGehee.
- **Objective**: Auditer `solvers/stumpff.py`, `solvers/kepler_sundman.py`, `solvers/mcgehee.py` et `integrators/common.py` (`step_drift_and_kepler`).
- **Deliverables**:
  - `Chapter3AuditReports`: 3 rapports d'audit spécialisés (Math-to-Code, Numérique, Algèbre de Lie).
- **Process**: Lancement concurrent de 3 sous-agents d'audit via `invoke_subagent`.

### Step 2: Déploiement Triadique sur le Chapitre 4 (Splitting FM, Apex Écranté & Super-nœuds)
- **Context**: Le Chapitre 4 formalise la déformation de Rees, le potentiel d'apex écranté $\Theta(u_v)$, les forces régularisées et l'horloge dépliée anti-Zénon des super-nœuds synchrones.
- **Objective**: Auditer `potential/screening.py`, `potential/screened_apex.py`, `potential/supernode_detector.py`, `solvers/supernode.py` et `integrators/common.py` (`step_apex_kick`).
- **Deliverables**:
  - `Chapter4AuditReports`: 3 rapports d'audit spécialisés (Math-to-Code, Numérique, Algèbre de Lie).
- **Process**: Lancement concurrent de 3 sous-agents d'audit via `invoke_subagent`.

### Step 3: Déploiement Triadique sur le Chapitre 5 (Intégrateurs d'Ordre Élevé, PBM & Dyadique)
- **Context**: Le Chapitre 5 définit les schémas de Blanes-Moan (BM4, BM6), les équations homologiques des correcteurs symplectiques $C_{2p}$, les intégrateurs traités PBM et l'ordonnancement asynchrone dyadique.
- **Objective**: Auditer `integrators/` (`strang.py`, `yoshida.py`, `blanes_moan_4.py`, `blanes_moan_6.py`, `correctors.py`, `processed.py`), `system/dyadic_scheduler.py` et `system/nbody_system.py`.
- **Deliverables**:
  - `Chapter5AuditReports`: 3 rapports d'audit spécialisés (Math-to-Code, Numérique, Algèbre de Lie).
- **Process**: Lancement concurrent de 3 sous-agents d'audit via `invoke_subagent`.

### Step 4: Synthèse Contradictoire Globale et Cartographie des Anomalies
- **Context**: Les 9 rapports doivent être confrontés pour isoler l'origine exacte des sauts d'énergie résiduels.
- **Objective**: Produire le document consolidé `audit_theorique_chapitres_3_4_5.md` avec matrice de conformité et plan d'action d'ingénierie.
- **Deliverables**:
  - `MasterAuditReport`: Rapport consolidé final dans le brain.
- **Process**: Compilation, recoupement contradictoire et validation de non-régression.

## Acceptance Criteria:
### Step 1:
- Confronter rigoureusement l'inversion cubique de Halley, le Wronskien unimodulaire et le rebond de McGehee avec les équations du Chapitre 3.
### Step 2:
- Démontrer la conformité de la partition de forces $F = F_{\text{leaf}} + F_{\text{super}} + F_{\text{apex}}$ et vérifier l'absence d'erreurs de troncature sur les super-nœuds par rapport au Chapitre 4.
### Step 3:
- Évaluer la validité des équations homologiques des correcteurs $C_4/C_6$, le télescopage dyadique sous transitions d'atlas et les coefficients de Blanes-Moan face au Chapitre 5.
### Step 4:
- Cartographie sans complaisance des causes exactes des sauts d'énergie avec propositions de correction formelles.

## Observations:
> [!IMPORTANT]
> Cet audit à 9 agents constitue le passage au crible le plus exhaustif jamais réalisé sur la base de code, confrontant le code ligne par ligne avec les 3 chapitres fondateurs du manuscrit de recherche sous trois regards d'experts étanches.
