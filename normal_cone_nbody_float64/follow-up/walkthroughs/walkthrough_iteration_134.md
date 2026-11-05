# Walkthrough - Itération 134 : Régénération des Planches Scientifiques 300 DPI et des Données de la Web App Interactive

## 1. Contexte & Objectif
À la suite de la mise en conformité théorique complète de la bibliothèque [`normal_cone_nbody_float64`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64) (flot complet $(\mathbf{q}, \mathbf{p})$ des correcteurs symplectiques $\mathcal{C}_{2p}$, régularisation conormale des supernœuds avec protection contre les divisions par zéro, traitement des collisions triples à $N=3$ et agrégation tétraédrique symétrique à $N=4$, dwell-time dynamique $\tau_{\min}$), la présente itération a exécuté la **régénération complète des planches scientifiques 300 DPI et des données de la Web App Interactive**.

---

## 2. Simulations Comparatives et Résultats Numériques

L'ensemble des 4 scénarios canoniques a été exécuté sur les 6 intégrateurs symplectiques (Strang $\mathcal{S}_2$, Yoshida $\mathcal{S}_4$, Blanes--Moan $\mathrm{BM}_4$, Processed $\mathrm{PBM}_4$, Blanes--Moan $\mathrm{BM}_6$, Processed $\mathrm{PBM}_6$) via [`scratch/run_showcase_simulations.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scratch/run_showcase_simulations.py).

### Synthèse des Erreurs Relatives Maximales d'Énergie $|\Delta H / H_0|$

| Intégrateur | Ordre | 1. Figure-8 (2 périodes) | 2. Burrau 3-Corps ($t \in [0, 10]$) | 3. Collision Radiale ($r \to 0$) | 4. 4-Corps Hiérarchique ($t \in [0, 8]$) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Strang $\mathcal{S}_2$** | 2 | $5.71 \times 10^{-5}$ | $2.48 \times 10^{-5}$ | $2.85 \times 10^{-3}$ | $4.75 \times 10^{-3}$ |
| **Yoshida $\mathcal{S}_4$** | 4 | $6.12 \times 10^{-9}$ | $3.09 \times 10^{-10}$ | $5.12 \times 10^{-6}$ | $1.43 \times 10^{-4}$ |
| **Blanes--Moan $\mathrm{BM}_4$** | 4 | $8.78 \times 10^{-12}$ | $1.32 \times 10^{-10}$ | $1.15 \times 10^{-8}$ | $1.36 \times 10^{-6}$ |
| **Processed $\mathrm{PBM}_4$** | 4 (eff. 6) | $1.28 \times 10^{-11}$ | $1.25 \times 10^{-10}$ | $1.15 \times 10^{-8}$ | **$1.36 \times 10^{-6}$** (stable) |
| **Blanes--Moan $\mathrm{BM}_6$** | 6 | **$3.16 \times 10^{-12}$** | **$1.11 \times 10^{-10}$** | **$3.04 \times 10^{-11}$** | **$1.96 \times 10^{-7}$** |
| **Processed $\mathrm{PBM}_6$** | 6 (eff. 8) | **$3.16 \times 10^{-12}$** | **$1.11 \times 10^{-10}$** | **$3.04 \times 10^{-11}$** | **$1.96 \times 10^{-7}$** |

---

## 3. Artefacts Scientifiques et Visuels Régénérés

Tous les artefacts ont été générés et synchronisés dans le répertoire des artefacts du brain ([`appDataDir/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/)) :

1. **Planche 300 DPI — Chenciner Figure-8 Choreography :**
   - Fichier : [`dashboard_figure8_bm6_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_figure8_bm6_300dpi.png) (1.91 MB)
   - Contenu : 10 panneaux analytiques (trajectoire 3D, résidus d'invariants de Lie--Bargmann, spectres FFT, analyse Pareto de convergence log-log montrant la pente d'ordre 4.02 pour Yoshida et le plancher machine de BM6).

2. **Planche 300 DPI — Problème Pythagoricien de Burrau :**
   - Fichier : [`dashboard_burrau_bm6_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_burrau_bm6_300dpi.png) (1.67 MB)
   - Contenu : 2001 pas physiques, 5 transitions de cartes d'atlas dans $\mathrm{SL}(2, \mathbb{R})$, passages au péricentre $r_{\min} = 0.0115$ sans adoucissement, erreur d'énergie bornée à $1.11 \times 10^{-10}$.

3. **Planche 300 DPI — Système Hiérarchique à 4 Corps :**
   - Fichier : [`dashboard_4body_bm6_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_4body_bm6_300dpi.png) (2.21 MB)
   - Contenu : 15 transitions d'atlas sous hystérésis $\gamma = 0.6$, ordonnancement dyadique de niveau 10 (1024 sous-pas), stabilité parfaite de $\mathrm{PBM}_4$ sans aucune divergence.

4. **Jeu de Données Numériques de Showcase :**
   - Fichier : [`showcase_data.json`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/showcase_data.json) (5.32 MB)
   - Données brutes de trajectoires échantillonnées, historique complet de restructuration de l'arbre binaire de Jacobi, et courbes énergétiques comparatives pour les 6 intégrateurs.

5. **Application Web Interactive Autonome :**
   - Fichier : [`interactive_conormal_showcase.html`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/interactive_conormal_showcase.html) (2.78 MB)
   - Synchronisation complète de la variable inline `SHOWCASE_DATA` : visualiseur 2D/3D Canvas dynamique, rendu SVG arborescent en temps réel de l'arbre de Jacobi, sélecteur multi-intégrateurs interactif avec affichage instantané des courbes comparatives.

---

## 4. Certification Finale
- 100% des planches 300 DPI et des données de showcase sont parfaitement synchronisées avec le solveur régularisé et conforme.
- Aucune régression : zéro divergence, zéro NaN, zéro adoucissement ($\varepsilon = 0, d_0 = 0$).
