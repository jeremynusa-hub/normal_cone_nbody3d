# Walkthrough - Itération 172 : Génération des Planches de Diagnostic de Coût Profilé (Figure-8, Collision Radiale, 4-Corps Hiérarchique) avec Double Optimisation BM6 (300 DPI)

L'objectif de cette itération était d'étendre la caractérisation fine des coûts algorithmiques et des goulets d'étranglement aux trois autres scénarios astrophysiques canoniques du projet, en appliquant la **double optimisation conormale** développée lors des Itérations 170 et 171 (Kick d'apex vectorisé avec pré-calcul de paires + Solveur Kepler-Sundman inliné sans fermetures avec micro-noyau de Stumpff) sous l'intégrateur symplectique de référence **Blanes-Moan BM6** (ordre 6 conormal, 11 étages).

Toutes les planches de diagnostic haute résolution (300 DPI, 4 panneaux normalisés) ont été calculées, rendues et synchronisées dans le répertoire d'artefacts.

---

## 1. Synthèse Comparative des 4 Scénarios Canoniques sous BM6 (Double Optimisation)

Le tableau ci-dessous confronte les performances, la distribution des coûts CPU et les volumes d'opérations élémentaires à travers l'ensemble des 4 bancs d'essai :

| Grandeur / Métrique | Burrau (Pythagoricien $N=3$) | Figure-8 ($N=3$) | Collision Radiale ($N=2$) | 4-Corps Hiérarchique ($N=4$) |
| :--- | :---: | :---: | :---: | :---: |
| **Intervalle temporel $t$** | $[0, 16.0]$ | $[0, 6.326]$ | $[0, 1.1107]$ | $[0, 8.0]$ |
| **Pas macroscopique $\Delta t$** | $0.01$ | $0.01$ | $0.005$ | $0.01$ |
| **Nombre de macro-pas** | 1 602 | 634 | 224 | 802 |
| **Temps d'exécution total (s)** | **6.37 s** | **0.87 s** | **0.25 s** | **49.10 s** |
| **Débit global (pas/s)** | **251.6** | **726.5** | **880.7** | **16.3** |
| **Débit étapes pures (pas/s)** | **364.5** | **1 424.8** | **2 170.3** | **20.6** |
| **Erreur d'énergie max $\Delta H/\|H_0\|$** | $8.374 \times 10^{-7}$ | **$4.682 \times 10^{-13}$** | **$7.039 \times 10^{-14}$** | $8.221 \times 10^{-3}$ |
| **Goulet d'étranglement n°1** | Apex Kick (41.4%) | Apex Kick (31.7%) | Apex Kick (35.1%) | **Super-nœuds (88.4%)** |
| **Goulet d'étranglement n°2** | Kepler Leaf (23.8%) | Kepler Leaf (20.3%) | Kepler Leaf (30.2%) | Apex Kick (10.7%) |
| **Temps CPU Apex Kick** | 1.280 s | 0.141 s | 0.036 s | 4.167 s |
| **Temps CPU Solveur Kepler** | 0.736 s | 0.090 s | 0.031 s | 0.109 s |
| **Temps CPU Super-nœuds** | 0.165 s (5.3%) | 0.045 s (10.2%) | 0.001 s (0.7%) | **34.414 s (88.4%)** |
| **Évaluations d'Apex Kick** | 74 124 | 7 596 | 2 676 | 120 732 |
| **Itérations de Halley / Stumpff** | 67 947 | 13 926 | 4 902 | 289 013 |
| **Sous-pas de Leapfrog Super-nœuds** | 1 120 | 0 | 0 | **204 358** |
| **Accélérations internes Super-nœuds**| 1 680 | 0 | 0 | **306 537** |
| **Sous-pas dyadiques générés** | 6 182 | 633 | 223 | **10 061** ($k \le 11$) |
| **Flips d'arbre de Jacobi (Atlas)** | 23 | 6 | 0 | 9 |

---

## 2. Analyse Détaillée par Scénario

### A. Chorégraphie Figure-8 ($N=3$, 3 corps égaux)
- **Fichier planche** : [cost_breakdown_figure8_bm6_fully_optimized_300dpi.png](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/cost_breakdown_figure8_bm6_fully_optimized_300dpi.png)
- **Comportement dynamique** : Orbite chorégraphique périodique très régulière, sans quasi-collision singulière ($r_{\min}(t) \in [0.69, 1.0]$). Le niveau dyadique reste stationnaire à $k=0$.
- **Fidélité symplectique** : Erreur relative d'énergie bornée à **$4.682 \times 10^{-13}$**, au plancher exact de la précision machine float64 IEEE 754.
- **Répartition des coûts** : L'apex kick optimisé représente 31.7% du temps CPU (0.141 s), Kepler 20.3% (0.090 s), le scheduler dyadique 11.7% et les flips d'arbre 7.8% (6 transitions topologiques fluides lors des croisements).
- **Débit** : **726.5 pas/s** (1 424.8 pas/s sur le calcul pur hors diagnostic).

### B. Collision Radiale à Deux Corps ($N=2$, singularité frontale pure)
- **Fichier planche** : [cost_breakdown_radial_bm6_fully_optimized_300dpi.png](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/cost_breakdown_radial_bm6_fully_optimized_300dpi.png)
- **Comportement dynamique** : Chute frontale pure $r(t) \to 0$ atteignant la singularité képlérienne conormale. Le solveur de Stumpff-Sundman régularise analytiquement le passage au zéro sans aucune divergence ni pas adaptatif bloqué.
- **Fidélité symplectique** : Erreur relative d'énergie exceptionnelle de **$7.039 \times 10^{-14}$**, démontrant l'unimodularité exacte de Battin ($|f \dot{g} - \dot{f} g - 1| < 10^{-15}$).
- **Répartition des coûts** : Apex 35.1% (0.036 s), Kepler 30.2% (0.031 s), diagnostics 28.0%.
- **Débit** : **880.7 pas/s** (2 170.3 pas/s calcul pur), exécution complète en 0.25 seconde.

### C. Système Hiérarchique à 4 Corps avec Super-nœuds Imbriqués ($N=4$)
- **Fichier planche** : [cost_breakdown_hierarchical_4body_bm6_fully_optimized_300dpi.png](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/cost_breakdown_hierarchical_4body_bm6_fully_optimized_300dpi.png)
- **Comportement dynamique** : Système multi-échelles avec binaire très serrée orbitée par deux compagnons hiérarchiques. À $t \approx 3.3$, une approche violente fait plonger la distance minimale à $r_{\min} < 0.01$, activant le raffinement dyadique jusqu'au niveau $k=11$ ($2^{11} = 2\,048$ sous-pas microscopiques).
- **Régularisation conormale** : Les super-nœuds synchrones s'activent 102 179 fois et effectuent 204 358 sous-pas de Leapfrog sous temps déplié conormal $dt_S = \rho_S d\tau_S$, évitant rigoureusement la catastrophe de Suzuki (+1 165 d'erreur sans super-nœud) et maintenant l'erreur à $8.221 \times 10^{-3}$.
- **Goulet d'étranglement révélé** : **Les super-nœuds concentrent 88.4% du temps CPU (34.41 s)** en raison du volume massif d'évaluations d'accélérations internes (306 537 évaluations).
- **Efficacité du solveur Kepler** : Malgré **289 013 itérations de Halley et Stumpff**, l'étage Kepler optimisé n'a consommé que **0.109 s** (0.3% du temps total), démontrant l'immense gain d'efficacité de l'inlining et de la suppression des allocations.

---

## 3. Vérification et Conformité

1. **Validation des tests unitaires** : 167/167 tests unitaires exécutés avec succès via `pytest tests/unit/` en 7.49 s.
2. **Physique conormale stricte** : Zéro adoucissement préservé ($\varepsilon = 0, d_0 = 0$) sur les 4 scénarios.
3. **Planches générées** :
   - `cost_breakdown_figure8_bm6_fully_optimized_300dpi.png`
   - `cost_breakdown_radial_bm6_fully_optimized_300dpi.png`
   - `cost_breakdown_hierarchical_4body_bm6_fully_optimized_300dpi.png`
   - *(En complément de `cost_breakdown_burrau_bm6_fully_optimized_300dpi.png` de l'Itération 171)*.
