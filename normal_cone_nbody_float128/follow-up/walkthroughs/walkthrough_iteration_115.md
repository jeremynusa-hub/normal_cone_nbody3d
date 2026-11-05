# Walkthrough of Iteration 115 : Implémentation et Certification Unitaire de la Couche 3 (potential/) - Crible Conormal $\Theta(u_v)$, Potentiel d'Apex $\mathcal{C}^{1,1}$, Détection de Marée de Hill et Newton Direct

## Summary:
L'Itération 115 a mené à bien la réalisation complète et la certification formelle de la **Couche 3 (`potential/`)** de la bibliothèque académique `normal_cone_nbody_float64`. En application stricte du **protocole à 5 étapes** exigé par l'utilisateur, l'itération a successivement :
1. Ré-assimilé les fondations mathématiques de la monographie [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/) (Chapitre 4, Sections 4.1–4.3, Lemmes 4.1–4.3, Théorème 4.1 de régularité universelle $\mathcal{C}^{1,1}$, Proposition 4.4 de seuil de marée de Hill $\eta_{\mathrm{crit}} = 3^{-1/3}$, Proposition 4.3 d'invariance de Lie--Bargmann) et extrait les théorèmes Lean 4 associés (`apex_universal_lipschitz_bound`, `apex_potential_bounded`, `collision_diagonal_non_transversality`, `lie_bargmann_exact_invariance`).
2. Audité et purgé les modules du code hérité dans `models/normal_cone_integrator/cas6_nbody3d_async/` (`potential/apex.py`, `potential/supernodes.py`, `potential/newtonian.py`), en éliminant les 42 occurrences historiques de `"cas6"` et en remplaçant les seuils dimensionnés ad-hoc ($\delta = 0.5$) par le critère de marée adimensionné.
3. Établi la spécification contractuelle des interfaces, invariants et tolérances d'invariance symplectique et de différentiation analytique.
4. Implémenté de manière modulaire et hautement optimisée les 4 modules de la Couche 3 (`screening.py`, `screened_apex.py`, `supernode_detector.py`, `direct.py`), intégrant le gradient analytique exact de la coupure d'apex (prenant en compte tant la variation du vecteur barycentrique relatif $\mathbf{r}_v$ que la dépendance interne des hyper-rayons massiques $\bar{\rho}_{\mathrm{int}}$).
5. Développé et validé la suite complète de 4 modules de tests unitaires miroirs dans `tests/unit/potential/` avec un taux de réussite de **100 % (17/17 tests passants)**, portant la couverture globale certifiée de la bibliothèque à **69 tests unitaires passants sans aucun échec ni avertissement (0.30s)**.

## Critical Discoveries:
1. **Gradient Analytique Exact de l'Apex et Dérivation par rapport aux Hyper-Rayons** : L'accord rigoureux entre gradient analytique et différences finies ($< 10^{-10}$) sur le terme de coupure conormale $W_v(\mathbf{q}) = \Theta(u_v) \frac{G M_L M_R}{\|\mathbf{r}_v\|}$ a exigé de dériver la fonction de crible non seulement par rapport au vecteur inter-amas $\mathbf{r}_v$, mais également par rapport à l'échelle d'hyper-rayon massique interne $\bar{\rho}_{\mathrm{int}} = \bar{\rho}_L + \bar{\rho}_R$ via $\nabla_{\mathbf{q}_i} \bar{\rho}_L = \frac{m_i (\mathbf{q}_i - \mathbf{X}_L)}{M_L \bar{\rho}_L}$. Cette composante préserve rigoureusement l'invariance par translation ($\sum_i \nabla_{\mathbf{q}_i} \bar{\rho}_L = \mathbf{0}$) et par rotation ($\sum_i \mathbf{q}_i \times \nabla_{\mathbf{q}_i} \bar{\rho}_L = \mathbf{0}$).
2. **Éradication Totale du Pôle Fantôme à la Coïncidence** : Pour $u_v \le 1.0$, l'annulation stricte $\Theta(u_v) = 0$ et $\Theta'(u_v) = 0$ éteint intégralement la soustraction képlérienne, éliminant tout pôle numérique en $1/\|\mathbf{r}_v\| \to \infty$ lorsque deux sous-amas se croisent transversalement avec $\rho_{\mathrm{int}} > 0$.
3. **Invariance d'Échelle de Hill et Détection Transitive** : Le rapport de marée adimensionné $\eta_S = \frac{\bar{\rho}_S}{r_P} \left(\frac{M_{\mathrm{ext}}}{M_S}\right)^{1/3}$ satisfait l'invariance exacte sous homothétie $\eta_S(\lambda \mathbf{q}) = \eta_S(\mathbf{q})$ à $< 10^{-15}$. La coalescence transitive fusionne sans ambiguïté les résonances multiples en super-nœuds disjoints.
4. **Zéro Adoucissement et Zéro Mention `"cas6"`** : Confirmation stricte de $\varepsilon = 0.0, d_0 = 0.0$ et absence totale de mentions `"cas6"` sur tous les fichiers de code et de test de la Couche 3.

## Realization:
- **Started At**: 2026-10-07T17:41:40+02:00
- **Completed At**: 2026-10-07T17:47:45+02:00

### Step 1: Relecture Théorique des Sections 4.1–4.3 de la Monographie
- **Impacted Elements**: Monographie LaTeX [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/)
- **Actions & Implementation**:
  - Relecture approfondie de `sec02_screened_apex_potential/content.tex` et `sec03_synchronous_supernodes/content.tex`.
  - Analyse du Théorème 4.1 ($\mathcal{C}^{1,1}$ Lipschitz Regularity of the Screened Apex Potential), de l'asymptotique quadrupolaire pure $\mathcal{Q}_v(\mathbf{s})$ et de la coupure de Cauchy/Hermite.
  - Formulation de la détection de marée de Hill $\eta_S \le 3^{-1/3}$ et du seuil naturel de Marchal--Sundman $\delta_{\mathrm{cluster}}^{\mathrm{nat}}(E_0)$.
  - Recensement des théorèmes Lean 4 certifiés : `apex_universal_lipschitz_bound`, `apex_potential_bounded`, `collision_diagonal_non_transversality`.

### Step 2: Audit Approfondi du Code Hérité (`cas6_nbody3d_async/potential/`)
- **Impacted Elements**: `models/normal_cone_integrator/cas6_nbody3d_async/`
- **Actions & Implementation**:
  - Cartographie des routines existantes et identification du couplage incomplet dans le gradient de crible de l'ancien `apex.py`.
  - Éradication des 42 mentions `"cas6"` et suppression de la constante ad-hoc `delta = 0.5`.

### Step 3: Spécification Formelle et Contrats d'Interface de la Couche 3
- **Impacted Elements**: Spécifications techniques des 4 modules de `potential/`
- **Actions & Implementation**:
  - Spécification de `screening.py` (fonctions d'évaluation scalaires et vectorisées pour Cauchy $\mathcal{C}^\infty$ et Hermite $\mathcal{C}^2$).
  - Spécification de `screened_apex.py` (tenseur quadrupolaire $\mathbf{Q}$, potentiel d'apex $\mathcal{C}^{1,1}$, forces d'apex analytiques, kick symplectique).
  - Spécification de `supernode_detector.py` (rapport de Hill $\eta_S$, seuil de Marchal--Sundman, coalescence transitive).
  - Spécification de `direct.py` (potentiel et forces newtoniennes de référence sans adoucissement $\varepsilon=0$).

### Step 4: Implémentation Modulaire et Optimisée des 4 Modules de `potential/`
- **Impacted Elements**:
  - [`normal_cone_nbody_float64/potential/__init__.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/__init__.py)
  - [`normal_cone_nbody_float64/potential/screening.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screening.py)
  - [`normal_cone_nbody_float64/potential/screened_apex.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screened_apex.py)
  - [`normal_cone_nbody_float64/potential/supernode_detector.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/supernode_detector.py)
  - [`normal_cone_nbody_float64/potential/direct.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/direct.py)
- **Actions & Implementation**:
  - `screening.py` : transition unitaire $C^\infty$ (intégrale de Cauchy) et spline quintique d'Hermite $C^2$ avec annulation stricte sur $u \le 1$ et unité asymptotique sur $u \ge 2$.
  - `screened_apex.py` : calcul vectorisé des interactions non-feuilles, masque des paires feuilles et super-nœuds, gradient analytique complet $\nabla_{\mathbf{q}_k} V_{\mathrm{apex}}$ incluant $\partial W_v / \partial \bar{\rho}_{\mathrm{int}}$, tenseur quadrupolaire sans trace $\operatorname{Tr}(\mathbf{Q}) = 0$.
  - `supernode_detector.py` : calcul du rapport de Hill $\eta_S$, discrimination des sous-amas candidats avec filtrage de l'arbre et coalescence transitive par composantes connexes.
  - `direct.py` : calcul vectorisé du potentiel et des forces coulombiennes pures ($\varepsilon = 0, d_0 = 0$).

### Step 5: Implémentation et Exécution des Tests Unitaires Miroirs (`tests/unit/potential/`)
- **Impacted Elements**:
  - [`tests/unit/potential/test_screening.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/potential/test_screening.py) (4 tests)
  - [`tests/unit/potential/test_screened_apex.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/potential/test_screened_apex.py) (5 tests)
  - [`tests/unit/potential/test_supernode_detector.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/potential/test_supernode_detector.py) (4 tests)
  - [`tests/unit/potential/test_direct.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/potential/test_direct.py) (4 tests)
- **Validation**:
  - `pytest tests/unit/potential -v` : **17 passed in 0.15s**.
  - `pytest tests/unit -v` (suite complète Couches 1, 2, 3) : **69 passed in 0.30s**.

## Verification:
- **Taux de passage** : 100 % (17/17 sur `potential/`, 69/69 au total).
- **Invariance de Lie--Bargmann** : Conservation de la force totale nulle ($\|\sum \mathbf{F}_i\| < 5 \times 10^{-14}$) et du couple total nul ($\|\sum \mathbf{q}_i \times \mathbf{F}_i\| < 5 \times 10^{-14}$) à l'échelle de précision machine pour des forces de norme $\sim 200$.
- **Conformité aux différences finies** : Gradient analytique d'apex conforme au gradient numérique à $< 10^{-10}$.
- **Éradication des pôles fantômes** : Annulation stricte de $\Theta(u_v)$ et gradient nul lors de la coïncidence barycentrique transversale.
- **Invariance d'échelle** : $\eta_S(\lambda \mathbf{q}) = \eta_S(\mathbf{q})$ à $< 10^{-15}$.
- **Zéro Softening & Zéro `cas6`** : $\varepsilon = 0.0$, $d_0 = 0.0$ rigoureux, aucune mention `"cas6"`.

## Next Steps:
- Clôturer l'Itération 115 sur `roadmap.md` et consigner les faits techniques dans `notes.md`.
- Entamer l'**Itération 116** dédiée à la **Couche 4 (`integrators/`)** :
  - Intégrateurs symplectiques palindromiques purs $\mathcal{S}_2$ (Strang), $\mathcal{S}_4$ (Forest--Ruth/Yoshida), $\mathrm{BM}_4^{(7)}$ et $\mathrm{BM}_6^{(11)}$ (Blanes--Moan).
  - Correcteurs canoniques $\mathcal{C}_{2p}$ (Wisdom--Holman / Blanes).
  - Intégrateur traité à haute performance $\mathcal{P}\mathrm{BM}$ (Processed Blanes--Moan).
