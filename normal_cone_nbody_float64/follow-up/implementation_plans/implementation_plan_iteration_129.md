# Implementation Plan of Iteration 129 : Réanalyse Exhaustive Multi-Agents du Calcul Hamiltonien au Sein des Nœuds (Kicks Virtuels, Splitting Conormal et Transitions d'Atlas) sur le Cas Burrau
Created at: 2026-10-07T20:20:00+02:00

## Context:
1. **Explicit Request**:
   - « Quand je regarde le cas bureau 3, l'erreur sur l'energie jump quand les 2 feuile d'un noeuds intermediaire interagisse fortement (particule m2 et m3), pourtant la même interaction entre m1 et m2 quelques secondes aprés ne créé aucune erreur. »
   - « Ca ressemble à un défaut qu'on avait déjà eu dans l'ancien solveur je pense. »
   - « Réanalyse exhaustivement le code notament sur les bons calculs de hamiltonien (de memoire on avait des kick virtuels) au sein d'un noeuds et compare ca avec le papier théorique et l'ancien code via l'appels d'agents diversifiés »
2. **Implicit Strategic Value**:
   - Résoudre le saut d'erreur d'énergie lors des interactions péricentriques fortes entre corps non-analytiques dans l'arbre actif (e.g. $m_2$ et $m_3$ vs $m_1$ et $m_2$ dans le problème de Burrau).
   - Élucider le mécanisme exact des « kicks virtuels » / splitting conormal interne au sein d'un nœud composite tel que développé dans l'ancien solveur ou prescrit par la décomposition de Fulton-MacPherson / Rees du papier théorique.
   - Confronter rigoureusement le code actuel (`normal_cone_nbody_float64/`), l'ancien code (`models/normal_cone_integrator/cas6_nbody3d_async/`), et le manuscrit théorique (`redaction/unsoftened_conormal_nbody3d/`).

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des scénarios.
- **Règle d'Édition** : Utiliser exclusivement `replace_file_content` pour les modifications du code Python, ne jamais écraser avec `overwrite: true`.

## Plan:
### Step 1: Déploiement d'un Audit Multi-Agents Spécialisé
- **Agent 1 (Legacy Codebase Archeologist)** :
  - Fouiller l'ancien code (`models/normal_cone_integrator/cas6_nbody3d_async/` et ses walkthroughs historiques) pour identifier le concept et l'implémentation des « kicks virtuels », du splitting hamiltonien par nœud, et de la gestion des collisions binaires secondaires dans Burrau.
- **Agent 2 (Monograph LaTeX Theorist)** :
  - Examiner dans `redaction/unsoftened_conormal_nbody3d/` (Chapitres 3, 4 et 5) les formules exactes du splitting d'apex, de l'interaction résiduelle au sein d'un cluster, du potentiel régularisé dans les coordonnées conormales de Rees et des conditions de transition d'atlas.
- **Agent 3 (Dynamic Diagnostics & Execution Profiler)** :
  - Profiler pas à pas la trajectoire de Burrau sur `normal_cone_nbody_float64/` pour localiser le timestamp exact $t_{\text{jump}}$, la topologie de l'arbre à cet instant, les distances relatives $r_{12}, r_{23}, r_{13}$, les composantes du gradient d'apex, et le comportement de l'`AtlasManager` et du `DyadicScheduler`.

### Step 2: Synthèse Contradictoire et Diagnostic de l'Écart Hamiltonien
- Réunir les conclusions des 3 sous-agents dans un rapport d'audit exhaustif.
- Déterminer si le saut d'erreur provient :
  1. D'un retard ou blocage de transition d'atlas (hystérésis $\gamma=0.6$, dwell time, ou seuil de transition empêchant $(m_2, m_3)$ de devenir la paire de feuilles analytique).
  2. D'un terme d'interaction interne au nœud composite (« kick virtuel » / potentiel résiduel au barycentre) manquant ou mal pondéré dans l'apex conormal.
  3. D'une activation inadéquate du super-nœud synchrone sur les corps $\{m_2, m_3\}$.

### Step 3: Conception et Implémentation du Correctif Algorithmique
- Corriger le module cible identifié (`system/atlas_manager.py`, `potential/screened_apex.py`, `integrators/common.py` ou `solvers/supernode.py`).
- Garantir la stricte conservation symplectique, l'invariance de Lie-Bargmann et le zéro-adoucissement ($\varepsilon = 0$).

### Step 4: Validation Numérique, Régénération des Données et Walkthrough
- Valider la disparition du jump d'énergie sur le problème de Burrau.
- Régénérer `showcase_data.json` et les dashboards de publication 300 DPI.
- Valider 100% de la suite de tests unitaires (`pytest tests/unit/`).

## Acceptance Criteria:
- **Critère 1** : L'interaction forte entre $m_2$ et $m_3$ dans le cas de Burrau ne produit plus de saut anormal d'erreur d'énergie par rapport à l'interaction entre $m_1$ et $m_2$.
- **Critère 2** : L'isomorphisme mathématique avec le papier théorique et la parité avec les mécanismes validés de l'ancien code sont formellement démontrés.
- **Critère 3** : 100% des tests unitaires passent (`pytest tests/unit/`).
