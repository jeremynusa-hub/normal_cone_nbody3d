# Implementation Plan of Iteration 127 : Visualisation Dynamique des Sous-Nœuds d'Arbres de Jacobi et Comparaison Multi-Intégrateurs Avec vs Sans Correcteurs Symplectiques
Created at: 2026-10-07T19:57:00+02:00

## Context:
1. **Explicit Request**:
   - « Dans les visuels on voit les super noeuds mais pas les sous noeuds, et les arbre semblent pas evoluer. »
   - « L'évolution des energie semble pas documenté pour tout les cas d'intégrateur. »
   - « De plus dans tout les cas, ou sont les intégrateurs blanes moan avec vs sans correcteur symplectic ? »
2. **Implicit Strategic Value**:
   Assurer une transparence épistémique et pédagogique totale sur les mécanismes de la géométrie conormale :
   - Démontrer visuellement la vie dynamique de l'arbre de Jacobi lors des basculements d'atlas (17 transitions sur le 4-corps), en distinguant clairement la racine, les sous-nœuds composites intermédiaires ($\sigma_v$ avec leur masse réduite $\mu_v$), les paires de feuilles analytiques et les super-nœuds compacts de Hill.
   - Fournir l'historique complet de conservation d'énergie pour tous les intégrateurs de la hiérarchie conormale (`strang`, `yoshida`, `blanes_moan_4`, `processed_bm4`, `blanes_moan_6`, `processed_bm6`) sur chaque scénario physique.
   - Mettre en lumière la théorie du Chapitre 5 (Section 5.5) démontrant le saut d'ordre effectif conféré par la conjugaison aux correcteurs symplectiques homologiques ($\mathcal{C}_4$ et $\mathcal{C}_6$) à zéro évaluation de force supplémentaire par télescopage frontière.
3. **Alignment Criteria**:
   - `showcase_data.json` enrichi avec la structure arborescente complète à chaque pas échantillonné et les courbes d'énergie des 6 intégrateurs sur les 4 scénarios.
   - Refonte du composant SVG dans `interactive_conormal_showcase.html` pour tracer un arbre binaire dynamique récursif (nœuds, arêtes, sous-nœuds composites étiquetés, paires de feuilles et surbrillance des super-nœuds).
   - Ajout dans la WebApp d'un sélecteur comparatif dédié « Avec vs Sans Correcteur Symplectique » (BM4 vs PBM4, BM6 vs PBM6).
   - Régénération des dashboards scientifiques 300 DPI et déploiement dans le répertoire des artefacts.

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des scénarios.
- **Rigueur Mathématique** : Respecter scrupuleusement la formulation de Rees et de Darboux pour les sous-nœuds de Jacobi et les opérateurs de conjugaison $\Psi_{2p+2} = \mathcal{C}_{2p} \circ S_{2p} \circ \mathcal{C}_{2p}^{-1}$.

## Plan:
### Step 1: Enrichissement du Générateur de Données (`scratch/run_showcase_simulations.py`)
- **Matériel de départ** : [`scratch/run_showcase_simulations.py`](file:///C:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scratch/run_showcase_simulations.py).
- **Objectif** :
  1. Implémenter une fonction de sérialisation récursive `serialize_tree(node)` capturant la topologie exacte (racine, sous-nœuds composites avec $\mu_v$, et feuilles $m_i$) et l'enregistrer dans `tree_history` pour chaque pas échantillonné.
  2. Intégrer les 6 solveurs (`strang`, `yoshida`, `blanes_moan_4`, `processed_bm4`, `blanes_moan_6`, `processed_bm6`) sur chaque scénario (Figure-8, Burrau, Radial, 4-corps) et exporter leurs séries d'erreurs d'énergie relatives $\log_{10} |\Delta H / |H_0||$.
- **Livrable** : Script enrichi et prêt pour le calcul multi-intégrateurs.
- **Processus** : Édition via `replace_file_content`.

### Step 2: Refonte Dynamique de l'Arbre et de la Comparaison d'Énergie dans la WebApp (`interactive_conormal_showcase.html`)
- **Matériel de départ** : [`interactive_conormal_showcase.html`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/interactive_conormal_showcase.html).
- **Objectif** :
  1. Remplacer le rendu SVG statique par un visualiseur d'arbre hiérarchique récursif calculant dynamiquement les coordonnées $(x, y)$ de chaque nœud : racine, sous-nœuds intermédiaires $\sigma_v$ avec badge de masse réduite, et nœuds terminaux $m_i$ avec leurs couleurs respectives.
  2. Matérialiser visuellement les transitions d'atlas en animant la reconfiguration de l'arbre lors des basculements de cartes de Jacobi.
  3. Mettre à jour le graphique d'énergie pour supporter les 6 intégrateurs avec un mode de filtrage dédié « Avec vs Sans Correcteur Symplectique » mettant en évidence le saut spectaculaire de précision entre BM4/PBM4 et BM6/PBM6.
- **Livrable** : Code HTML/JS/CSS de l'application web interactive actualisé.
- **Processus** : Édition ciblée via `replace_file_content`.

### Step 3: Exécution des Simulations, Synchronisation et Déploiement
- **Matériel de départ** : Scripts actualisés dans `scratch/`.
- **Objectif** :
  1. Lancer l'exécution de `run_showcase_simulations.py` pour générer `showcase_data.json` et les 3 planches 300 DPI.
  2. Injecter les données dans `interactive_conormal_showcase.html` via `sync_showcase_html.py`.
  3. Déployer l'ensemble des fichiers dans le dossier brain.
- **Livrables** : `showcase_data.json`, `interactive_conormal_showcase.html` et les 3 images 300 DPI actualisées.
- **Processus** : Commandes asynchrones et synchronisation.

### Step 4: Contrôle Qualité et Validation Scientifique
- **Matériel de départ** : Artefacts déployés.
- **Objectif** :
  1. Valider que l'arbre de Jacobi évolue pas à pas lors des 17 transitions et met en valeur les sous-nœuds composites.
  2. Vérifier que les courbes d'énergie des 6 intégrateurs sont visibles et cohérentes sur tous les scénarios.
  3. Documenter les résultats dans `notes.md` et finaliser `walkthrough.md`.
- **Livrable** : Walkthrough de validation complète de l'Itération 127.
- **Processus** : Inspection visuelle et métrique.

## Acceptance Criteria:
- **Critère 1** : L'arbre de Jacobi dans la WebApp est généré dynamiquement à partir de la topologie réelle, affichant distinctement la racine, les sous-nœuds composites avec masse réduite $\mu_v$, et les feuilles $m_i$.
- **Critère 2** : L'évolution de l'arbre reflète fidèlement les basculements d'atlas (changement des paires de feuilles lors des transitions).
- **Critère 3** : Le graphique d'énergie documente les 6 intégrateurs (`strang`, `yoshida`, `blanes_moan_4`, `processed_bm4`, `blanes_moan_6`, `processed_bm6`) sur chaque scénario.
- **Critère 4** : La comparaison explicite BM4 vs PBM4 et BM6 vs PBM6 montre clairement le gain d'ordre apporté par les correcteurs symplectiques $\mathcal{C}_4$ et $\mathcal{C}_6$.
- **Critère 5** : Les 3 planches 300 DPI sont régénérées et certifiées sans régression.

## Observations:
> [!NOTE]
> La décomposition récursive de l'arbre de Jacobi permet de visualiser en direct la structure de Whitney et la désingularisation de Fulton-MacPherson à chaque pas physique.
> [!IMPORTANT]
> Les intégrateurs traités (PBM) appliquent les correcteurs $\mathcal{C}_{2p}$ uniquement aux frontières de blocs grâce au télescopage dyadique, conservant un coût de calcul identique aux intégrateurs de base tout en gagnant +2 ordres effectifs.
