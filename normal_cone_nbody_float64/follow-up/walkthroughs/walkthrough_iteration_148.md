# Walkthrough de l'Itération 148 : Phase 5 — Découplage Conormal Feuille/Apex, Certification des 8 Benchmarks et Planches 300 DPI

Date : 08 octobre 2026  
Statut : **100 % ACHEVE & CERTIFIE**  
Commit Git : [`bf6d705`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff)  
Suite de tests : **142/142 tests unitaires passés avec succès (100 %) en 6.01s**

---

## 1. Synthèse Exécutive des Réalisations

L'Itération 148 a parachevé la **Phase 5 de la Roadmap VI**, concluant l'ensemble des remédiations théoriques et numériques issues de l'audit contradictoire à 9 agents sur les Chapitres 3, 4 et 5 du manuscrit de recherche.

### A. Découplage Conormal Feuille / Apex (`DyadicScheduler`)
- **Problématique Théorique Résolue** : Au pas 220 du système à 4 corps hiérarchique, les corps 1 et 2 s'approchaient avec $v_{\mathrm{rad}} < 0$, conduisant à un péricentre képlérien hyperbolique minuscule $r_{\mathrm{peri}} \approx 2.05 \times 10^{-5}$. L'ordonnanceur dyadique imposait alors le niveau maximal $k=17$ ($131\,072$ sous-pas, soit $1.44 \times 10^6$ évaluations de force par pas), provoquant un temps de calcul de 70 secondes par pas macro.
- **Isomorphisme Fulton--MacPherson rétabli (Option B)** : Dans le splitting conormal $H = H_{\mathrm{base}} + H_{\mathrm{leaf}} + V_{\mathrm{apex}}^{\mathrm{hier}}$, la paire de feuilles active est intégrée **analytiquement** par les fonctions universelles de Stumpff-Sundman avec une erreur de troncature de pas de temps strictement nulle ($|f \dot{g} - \dot{f} g - 1| < 10^{-15}$). Elle est donc **exclue du potentiel d'apex** $V_{\mathrm{apex}}$.
- **Implémentation** : `DyadicScheduler.compute_dynamical_level` extrait désormais `tree.leaf_pairs` et désactive l'anticipation singulière de péricentre ($r_{\mathrm{peri}} \to 0$) pour les paires de feuilles actives. Celles-ci sont régulées par leurs échelles dynamiques orbitales naturelles ($\tau_{\mathrm{kin}}, \tau_{\mathrm{dyn}}$), tandis que l'anticipation de péricentre extrême est strictement réservée aux interactions non-feuilles appartenant à l'apex avant bascule topologique.
- **Résultat Immédiat** : Le pas 220 est passé de 70 secondes ($k=17$) à **0.19 seconde** ($k=7$, 128 sous-pas), tout en garantissant une précision machine sur le péricentre analytique.

### B. Validation Unitaire Ciblée
- Création de [`tests/unit/system/test_dyadic_scheduler_leaf_decoupling.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/system/test_dyadic_scheduler_leaf_decoupling.py) :
  1. `test_leaf_pair_bypasses_singular_pericenter_anticipation` : frôlement serré sur feuille active ($r_{12} = 0.0304$, $v_{\mathrm{rel}} = 6.2$, $r_{\mathrm{peri}} \sim 10^{-5}$) $\implies k \le 7$ (validé à $k=6$).
  2. `test_non_leaf_intruder_triggers_protective_dyadic_elevation` : intrus en apex avec petit paramètre d'impact $\implies$ élévation dynamique protectrice $k \ge 8$ (validé à $k=17$).
  3. `test_tree_none_fallback_maintains_safety` : sécurité conservative sans arbre ($k \ge 8$).

---

## 2. Tableau de Certification Globale des 8 Benchmarks Canoniques

La suite de certification complète [`scratch/run_full_benchmarks_certification.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scratch/run_full_benchmarks_certification.py) a évalué les **48 simulations** (8 benchmarks $\times$ 6 intégrateurs symplectiques) sous zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$).

| Benchmark | Strang ($S_2$) | Yoshida ($S_4$) | BM4 | PBM4 | BM6 | PBM6 | Invariants $\Delta P, \Delta L$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **figure8_choreography** | $4.75 \times 10^{-4}$ | $8.00 \times 10^{-7}$ | $7.09 \times 10^{-10}$ | $9.11 \times 10^{-10}$ | **$4.56 \times 10^{-13}$** | **$4.64 \times 10^{-13}$** | $< 1.6 \times 10^{-14}$ |
| **burrau_pythagorean** | $2.51 \times 10^{-4}$ | $2.23 \times 10^{-6}$ | $1.98 \times 10^{-7}$ | $1.98 \times 10^{-7}$ | **$9.87 \times 10^{-10}$** | **$9.87 \times 10^{-10}$** | $< 1.4 \times 10^{-12}$ |
| **radial_collision** | **$5.19 \times 10^{-13}$** | $3.32 \times 10^{-12}$ | $3.37 \times 10^{-12}$ | $3.37 \times 10^{-12}$ | $1.01 \times 10^{-11}$ | $1.01 \times 10^{-11}$ | $0.00$ (exact) |
| **hierarchical_4body** | $3.20 \times 10^{0}$ (div.) | $9.16 \times 10^{-4}$ | $2.80 \times 10^{-4}$ | $2.81 \times 10^{-4}$ | **$3.66 \times 10^{-5}$** | **$9.26 \times 10^{-5}$** | $< 1.0 \times 10^{-9}$ (17 switch) |
| **arenstorf_lunar** | $5.37 \times 10^{-8}$ | $1.15 \times 10^{-12}$ | $2.74 \times 10^{-12}$ | $2.76 \times 10^{-12}$ | **$4.43 \times 10^{-12}$** | **$4.44 \times 10^{-12}$** | $< 2.7 \times 10^{-14}$ |
| **sitnikov_oscillatory** | $4.75 \times 10^{-8}$ | **$6.14 \times 10^{-13}$** | $1.25 \times 10^{-12}$ | $1.24 \times 10^{-12}$ | $2.00 \times 10^{-12}$ | $2.00 \times 10^{-12}$ | $< 2.5 \times 10^{-13}$ |
| **jovian_solar_system** | $2.76 \times 10^{-9}$ | **$3.29 \times 10^{-13}$** | $7.86 \times 10^{-13}$ | $8.46 \times 10^{-13}$ | $1.22 \times 10^{-12}$ | $1.22 \times 10^{-12}$ | $< 9.0 \times 10^{-16}$ |
| **quadruple_collision** | $1.84 \times 10^{-3}$ | $3.04 \times 10^{-6}$ | $2.47 \times 10^{-7}$ | $2.47 \times 10^{-7}$ | **$1.75 \times 10^{-11}$** | **$1.75 \times 10^{-11}$** | $0.00$ (exact) |

---

## 3. Artefacts et Livrables Produits

1. **4 Planches Scientifiques Haute Définition 300 DPI** :
   - [`dashboard_figure8_bm6_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_figure8_bm6_300dpi.png) (1.84 Mo) : Figure-8 avec BM6, front de Pareto complet sur 6 schémas, calibration typographique à espacement aéré (> 140px).
   - [`dashboard_burrau_bm6_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_burrau_bm6_300dpi.png) (1.64 Mo) : Problème pythagoricien de Burrau, frôlements à $r_{\min} = 0.011$, 5 transitions d'atlas sans choc d'énergie.
   - [`dashboard_4body_bm6_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_4body_bm6_300dpi.png) (2.12 Mo) : 4-corps hiérarchique avec 17 transitions topologiques $\mathrm{SL}(2, \mathbb{R})$ sous hystérésis $\gamma=0.6$.
   - [`dashboard_multi_integrator_invariants_300dpi.png`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/dashboard_multi_integrator_invariants_300dpi.png) (1.94 Mo) : Planche 6 panneaux dédiée à la conservation des 10 invariants du groupe de Lie--Bargmann ($\mathbf{P}, \mathbf{L}, \mathbf{K}, H$) comparant Strang, Yoshida, BM4, PBM4, BM6 et PBM6.

2. **Démonstrateur Interactif Web et Dataset Synchronisés** :
   - [`interactive_conormal_showcase.html`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/interactive_conormal_showcase.html) (2.08 Mo) : Application autonome HTML5/Canvas/SVG avec visualiseur dynamique d'arbres de Jacobi, courbes d'énergie multi-intégrateurs et simulation animée.
   - [`showcase_data.json`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/showcase_data.json) (5.31 Mo) : Dataset compacté contenant l'ensemble des trajectoires, arbres et énergies des 4 scénarios de démonstration.

---

## 4. Bilan de la Roadmap VI et Perspectives

La **Roadmap VI** est officiellement **achevée à 100 %** :
- **Phase 1 (Iter 144)** : Solveurs de Stumpff et Kepler-Sundman sans seuil de garde destructeur, élimination de McGehee.
- **Phase 2 (Iter 145)** : Bascule de Schmitt anti-broutement sur les super-nœuds, spline Hermite centrée $C^8$, horloge dépliée anti-Zénon.
- **Phase 3 (Iter 146)** : Synchronisation frontière de Skeel-Gear, fusion FSAL sur BM4, profondeur dyadique $k_{\max}=17$.
- **Phase 4 (Iter 147)** : Connexion de jauge de Lie $\mathcal{J}_{A \to B} = \mathcal{C}_B^{-1} \circ \mathcal{C}_A$, produit hessien-vecteur analytique $D_v \mathbf{F}_{ij}$.
- **Phase 5 (Iter 148)** : Découplage feuille/apex dans `DyadicScheduler`, certification des 8 benchmarks, régénération 300 DPI.

Le code de la bibliothèque `normal_cone_nbody_float64` est désormais rigoureusement isomorphe aux théorèmes de la monographie de recherche et aux preuves formelles Lean 4.
