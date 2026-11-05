# Walkthrough of Iteration 130 : Éradication Définitive du Saut d'Énergie à 10⁻⁸ par Anticipation Képlérienne du Péricentre dans l'Ordonnanceur Dyadique

## Summary:
L'Itération 130 a éradiqué définitivement le saut d'erreur d'énergie résiduel à $10^{-8}$ observé sur le graphique de comparaison multi-intégrateurs du problème de Burrau ($m_1=3, m_2=4, m_3=5$) dans l'application web interactive `interactive_conormal_showcase.html`.
Ce saut touchait spécifiquement les schémas d'ordre 4 : Blanes-Moan 4 (BM4, qui atteignait $9.28 \times 10^{-9}$) et Yoshida 4 ($6.70 \times 10^{-8}$) lors du passage au péricentre à $t \approx 8.76$ ($r_{\text{peri}} \approx 0.0085$).
Une analyse fine a révélé que l'ordonnanceur dyadique (`DyadicScheduler.compute_dynamical_level`) évaluait la proximité des corps de manière purement **réactive** à partir de leur distance statique au début du pas macro ($r(t_n) = 0.0890 > 0.05$), n'attribuant que le niveau dyadique 4 ($2^4 = 16$ sous-pas). En raison d'une forte vitesse radiale convergente ($v_{\text{rad}} \approx -13.4\,\mathrm{m/s}$), les deux corps plongeaient à $r_{\text{peri}} = 0.0085$ au sein de ce même intervalle macro, provoquant un saut d'erreur de troncature sur le kick d'apex $\Delta H \approx 1.18 \times 10^{-7}$.
Grâce à l'implémentation d'une **anticipation képlérienne et balistique à deux corps** dans l'ordonnanceur dyadique (calcul prédictif de $r_{\text{peri}} = a(1-e)$ pour les paires en approche $v_{\text{rad}} < 0$), l'ordonnanceur anticipe le péricentre et affine préventivement le maillage dyadique jusqu'au niveau 9 ($2^9 = 512$ sous-pas).
Les résultats numériques confirment la suppression totale de toute discontinuité :
- **Blanes-Moan 4 (BM4)** : l'erreur relative chute de $9.28 \times 10^{-9}$ à **$1.324 \times 10^{-10}$** (gain d'un facteur 70, profil parfaitement lisse).
- **Yoshida 4** : l'erreur relative chute de $6.70 \times 10^{-8}$ à **$3.086 \times 10^{-10}$** (gain d'un facteur 217).
- **Processed BM4 (PBM4)** : atteint **$1.197 \times 10^{-10}$**.
- **Blanes-Moan 6 (BM6)** : préserve son exactitude optimale à **$1.112 \times 10^{-10}$**.
- **Processed BM6 (PBM6)** : atteint **$1.311 \times 10^{-10}$**.
- **Strang (S2)** : $2.484 \times 10^{-5}$ (pente d'ordre 2 pure).

Tous les 125 tests unitaires passent à 100% en 8.02s. L'ensemble des 4 simulations astrophysiques, les 3 dashboards 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`), `showcase_data.json` et `interactive_conormal_showcase.html` ont été régénérés et synchronisés.

---

## Critical Discoveries:
1. **Origine du saut à $10^{-8}$ : Déficit Prédictif de l'Ordonnanceur Dyadique** :
   - À $t = 8.750$, la distance mutuelle entre $m_2$ et $m_3$ est $r_{23} = 0.1477$.
   - À $t = 8.755$, $r_{23} = 0.0890$.
   - L'algorithme évaluait `min_sep = 0.0890`. Le barème par paliers attribuait `prox_k = 4` ($2^4 = 16$ sous-pas) car $0.0890 > 0.05$.
   - Cependant, la vitesse d'approche relative était $v_{\text{rad}} = \frac{\mathbf{r} \cdot \mathbf{v}}{r} \approx -13.4\,\mathrm{m/s}$.
   - Au cours du pas macro $\Delta t = 0.005$, les deux corps franchissaient le péricentre à $r_{\text{peri}} = 0.008516$. Traverser une distance aussi infime avec seulement 16 sous-pas de kick d'apex entraînait une erreur de troncature locale d'apex de $\sim 10^{-7}$, se manifestant par une marche d'escalier à $10^{-8}$ pour les schémas d'ordre 4.
2. **Anticipation Képlérienne Balistique à Deux Corps** :
   - Pour toute paire présentant une vitesse radiale convergente ($v_{\text{rad}} < 0$), l'énergie orbitale locale à deux corps $E = \frac{1}{2} v^2 - \frac{\mu}{r}$ et le moment cinétique spécifique $h = \|\mathbf{r} \times \mathbf{v}\|$ permettent de déduire le demi-grand axe $a$ et l'excentricité $e = \sqrt{1 + 2 E h^2 / \mu^2}$.
   - Le péricentre théorique $r_{\text{peri}} = a(1-e)$ ainsi que la projection linéaire balistique $r_{\text{proj}} = r + v_{\text{rad}} \Delta t$ fournissent la distance minimale attendue sur le bloc à venir : $r_{\text{projected}} = \min(r_{ij}, r_{\text{proj}}, r_{\text{peri}})$.
   - Lorsque $r_{\text{projected}} \le 0.015$, le niveau dyadique est instantanément porté à $k = 9$ ($2^9 = 512$ sous-pas) *avant* l'engagement dans le puits gravitationnel, éliminant tout résidu de troncature.
3. **Harmonie Totale des Ordres 4 et 6** :
   - Sur l'ensemble du profil $t \in [0, 10]$ de Burrau, aucune courbe d'énergie d'aucun intégrateur d'ordre supérieur ne dépasse désormais $3.1 \times 10^{-10}$.
   - La hiérarchie asymptotique Strang ($2.5 \times 10^{-5}$) $\gg$ Yoshida ($3.1 \times 10^{-10}$) $\approx$ BM4 ($1.3 \times 10^{-10}$) $\approx$ BM6 ($1.1 \times 10^{-10}$) est certifiée.

---

## Realization:
- **Started At**: 2026-10-07T20:35:00+02:00
- **Completed At**: 2026-10-07T20:46:30+02:00

### Step 1: Diagnostic Numérique Fin du Saut à 10⁻⁸
- **Impacted Elements**: `normal_cone_nbody_float64/system/dyadic_scheduler.py`, `scratch/`
- **Actions & Implementation**:
  - Traçage précis de l'intégration pas à pas pour localiser l'instant du saut d'énergie ($t \in [8.750, 8.760]$).
  - Évaluation des grandeurs cinématiques : $r_{23} = 0.0890 \to 0.0085$, $v_{\text{rad}} \approx -13.4$.
  - Constat de la sous-estimation du niveau dyadique (niveau 4 au lieu du niveau 8-9 requis).

### Step 2: Conception et Implémentation de l'Anticipation Képlérienne du Péricentre
- **Impacted Elements**: [`normal_cone_nbody_float64/system/dyadic_scheduler.py`](file:///C:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/system/dyadic_scheduler.py#L140-L180)
- **Actions & Implementation**:
  - Édition ciblée de `compute_dynamical_level` via `replace_file_content` :
    - Détection des paires en approche ($v_{\text{rad}} < 0$).
    - Calcul du péricentre képlérien local à deux corps $r_{\text{peri}} = a(1-e)$ et de l'extrapolation linéaire.
    - Raffinement prédictif : $r_{\text{projected}} \le 0.015 \implies k=9$, $r_{\text{projected}} \le 0.03 \implies k=7$.
  - Maintien du zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$).

### Step 3: Contrôle Qualité, Régénération des Simulations et Synchronisation
- **Impacted Elements**: `scratch/run_showcase_simulations.py`, `scratch/sync_showcase_html.py`, `scratch/showcase_data.json`, `interactive_conormal_showcase.html`, `brain/`
- **Actions & Implementation**:
  - Exécution complète de `pytest tests/unit/` : 125/125 tests réussis à 100% en 8.02s.
  - Exécution de `run_showcase_simulations.py` :
    - Burrau : BM4 à $1.324 \times 10^{-10}$, Yoshida à $3.086 \times 10^{-10}$, BM6 à $1.112 \times 10^{-10}$, PBM4 à $1.197 \times 10^{-10}$, PBM6 à $1.311 \times 10^{-10}$.
    - Figure-8 : BM6 à $3.161 \times 10^{-12}$, scaling Pareto certifié.
    - Collision radiale : rebond McGehee avec erreur $3.039 \times 10^{-11}$.
    - 4-corps hiérarchique : 15 transitions d'atlas, niveau dyadique 10, PBM6 à $1.914 \times 10^{-7}$.
  - Régénération des 3 dashboards 300 DPI et copie dans l'artifact directory.
  - Exécution de `sync_showcase_html.py` mettant à jour `interactive_conormal_showcase.html` avec les nouvelles données.

---

## Design Choices, Trade-offs and Warnings:
- **Anticipation Képlérienne Locale sans Coût Global** :
  L'anticipation képlérienne $r_{\text{peri}} = a(1-e)$ est évaluée analytiquement sur les paires présentant $v_{\text{rad}} < 0$ et $r < 0.2$. Son coût d'évaluation est négligeable ($\mathcal{O}(N^2)$ opérations scalaires rapides) par rapport à l'évaluation des forces d'apex, tout en évitant les sauts de troncature catastrophiques.
- **Zéro Adoucissement Intact** :
  Aucun paramètre de softening $\varepsilon$ ni de seuil de lissage $d_0$ n'a été introduit. La régularité découle exclusivement de la géométrie de Fulton-MacPherson et de la résolution dyadique prédictive.

---

## Verification:
### Acceptance Verification:
- [x] L'erreur d'énergie sur le problème de Burrau ne présente plus aucun saut à $10^{-8}$ : BM4 passe de $9.28 \times 10^{-9}$ à **$1.324 \times 10^{-10}$**, Yoshida de $6.70 \times 10^{-8}$ à **$3.086 \times 10^{-10}$**.
- [x] Aucun saut de troncature d'apex n'est constaté au passage au péricentre critique $r_{\text{peri}} \approx 0.0085$ ($t \approx 8.76$).
- [x] 100% des tests unitaires passent sans régression : **125/125 réussis** en 8.02s.
- [x] `showcase_data.json`, les 3 planches de publication 300 DPI et `interactive_conormal_showcase.html` sont synchronisés.
- [x] `notes.md` est mis à jour selon la règle append-only.

---

## Next Steps:
- Explorer l'introduction d'un correcteur d'ordre supérieur adapté aux transitions d'atlas multiples pour les intégrateurs traités PBM sur le 4-corps hiérarchique.
- Étudier la parallélisation GPU / multi-cœurs du calcul dyadique pour des systèmes à grand nombre de corps ($N \ge 100$).
