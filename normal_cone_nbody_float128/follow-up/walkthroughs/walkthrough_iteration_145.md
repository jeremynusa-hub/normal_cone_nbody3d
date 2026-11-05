# Walkthrough of Iteration 145 : Phase 2 — Remédiation des Super-Nœuds et Régularisation du Potentiel d'Apex du Chapitre 4

## Summary:
L'itération 145 a mené à bien la **Phase 2 de la Roadmap VI** en éradiquant les 3 anomalies fondamentales identifiées par l'audit théorique du Chapitre 4 au sein de la décomposition conormale de Fulton--MacPherson :
1. **Éradication du broutement topologique (flickering)** : Mise en œuvre d'une bascule de Schmitt à hystérésis $[\eta_{\mathrm{enter}} = 0.60, \eta_{\mathrm{exit}} = 0.75]$ dans [`potential/supernode_detector.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/supernode_detector.py), interdisant toute alternance intempestive d'activation/désactivation sous oscillations de marée autour de $\eta_{\mathrm{crit}} \approx 0.693$, et transmission de la mémoire d'état `active_supernodes` dans [`system/nbody_system.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/system/nbody_system.py).
2. **Formulation polynomiale centrée conditionnée ($\kappa \approx 1$) de la spline Hermite $C^8$** : Réécriture de $\sigma_8(s)$ et $\sigma_8'(s)$ dans [`potential/screening.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screening.py) sous forme centrée en puissances paires de $\xi = s - 0.5$ ($\kappa \approx 1$), avec 9 coefficients dyadiques exacts en puissances de 2 (zéro erreur d'arrondi), éliminant l'annulation soustractive à $s=0.5$ et assurant une stricte monotonie ($d\Theta/du \ge 0$).
3. **Solveur de super-nœuds sous horloge dépliée de Sundman sans adoucisseur** : Réécriture du propagateur interne dans [`solvers/supernode.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/supernode.py) sous temps régularisé $\mathrm{d}t_S = \rho_S\,\mathrm{d}\tau_S$, éradication totale de l'adoucisseur masqué $d_{\mathrm{safe}} = 10^{-14}$ ($d_{\mathrm{safe}} \equiv 0$), suppression du recentrage artificiel de vitesse et projection canonique de Dirac du résidu machine sur la variété $\sum m_i v_i = 0$.
4. **Validation et certification unitaire complète** : 132/132 tests réussis à 100 % avec traçabilité Git commitée (`a5a5d25`).

---

## Critical Discoveries:
- **Exactitude binaire des coefficients centrés de $C^8$** :
  L'intégration formelle de $\sigma_8'(\xi) = 218790 (0.25 - \xi^2)^8$ a révélé que tous les dénominateurs des 9 coefficients $c_k = \frac{218790 \binom{8}{k} (-1)^k (0.25)^{8-k}}{2k+1}$ sont des puissances de 2 ($2^0, 2^1, 2^2, 2^6, 2^6, 2^9, 2^{10}, 2^{15}$). Ils sont donc **rigoureusement exacts** sans le moindre bit de troncature en double précision float64 IEEE 754.
- **Origine du drift de quantité de mouvement sans recentrage** :
  L'analyse fine a démontré que la dérive linéaire de $1.42 \times 10^{-14}$ observée après 320 sous-pas provenait de la division flottante par les masses différentes $m_i$ ($F_{ij}/m_i \times m_i \ne F_{ij}$ à $10^{-16}$). Une projection canonique de Dirac sur la variété de contrainte $\sum m_i v_i = 0$ élimine ce résidu machine à chaque pas sans altérer l'énergie.
- **Seuil d'agrégation tétraédrique vs détection de Hill** :
  L'agrégation symétrique tétraédrique ($d \le 1.5 \bar{\rho}_S$) absorbait le perturber externe lorsque celui-ci avait la même masse ($M_{\mathrm{ext}} = 1.0$), car $\eta \le 0.55 \implies r_P \le 1.26 \bar{\rho}_S \le 1.5 \bar{\rho}_S$. En calibrant la masse externe à $M_{\mathrm{ext}} = 8.0$, $r_P / \bar{\rho}_S = 2.52 > 1.5$, isolant parfaitement la bascule de Schmitt.

---

## Realization:
- **Started At**: 2026-10-08T10:16:30+02:00
- **Completed At**: 2026-10-08T10:23:00+02:00

### Step 1: Bascule de Schmitt à Hystérésis dans `potential/supernode_detector.py`
- **Impacted Elements**: [`normal_cone_nbody_float64/core/constants.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/core/constants.py), [`normal_cone_nbody_float64/core/__init__.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/core/__init__.py), [`normal_cone_nbody_float64/potential/supernode_detector.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/supernode_detector.py), [`normal_cone_nbody_float64/system/nbody_system.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/system/nbody_system.py).
- **Actions & Implementation**:
  - Définition des constantes certifiées `ETA_ENTER_HILL = 0.60` et `ETA_EXIT_HILL = 0.75` dans `constants.py`.
  - Extension de `detect_supernodes` avec le paramètre `active_supernodes` : maintien des amas actifs jusqu'à $\eta_{\mathrm{exit}} = 0.75$, et condition d'entrée stricte $\eta \le \eta_{\mathrm{enter}} = 0.60$ pour les nouveaux amas (repli sur $\eta_{\mathrm{crit}}$ si stateless).
  - Persistance de `self.active_supernodes` dans `NBodySystem` et transmission à chaque macro-pas.

### Step 2: Formulation Centrée Conditionnée ($\kappa \approx 1$) de la Spline $C^8$ dans `potential/screening.py`
- **Impacted Elements**: [`normal_cone_nbody_float64/potential/screening.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screening.py).
- **Actions & Implementation**:
  - Dérivation et implémentation du tuple `_C8_CENTERED_COEFFS` en puissances paires de $\xi^2 = (s - 0.5)^2$.
  - Refonte hybride de `eval_screening_theta` et `eval_screening_theta_vectorized` : évaluation centrée $(\kappa \approx 1)$ pour $|\xi| \le 0.25$, et forme factorisée $p^9 P_8(p)$ pour $|\xi| > 0.25$ (zéro annulation soustractive aux bords $u=1$ et $u=2$).
  - Vérification de la continuité machine à la frontière $|\xi| = 0.25$ ($|\Delta\Theta| = 0.00\mathrm{e}{+00}$) et monotonie stricte sur 1 million de points.

### Step 3: Solveur de Super-Nœuds sous Horloge Dépliée de Sundman sans Adoucisseur dans `solvers/supernode.py`
- **Impacted Elements**: [`normal_cone_nbody_float64/solvers/supernode.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/solvers/supernode.py).
- **Actions & Implementation**:
  - Extraction de `_compute_internal_cluster_acc` avec suppression totale de l'adoucisseur masqué ($d_{\mathrm{safe}} \equiv 0$).
  - Implémentation du sous-pas adaptatif sous horloge dépliée de Sundman $\mathrm{d}t_S = \rho_S\,\mathrm{d}\tau_S$, bornant le pas local par $\min(0.02\tau_{\mathrm{dyn}}, 0.05\tau_{\mathrm{sundman}}, \Delta t_{\mathrm{rem}})$.
  - Remplacement de la dissipation artificielle de vitesse par la projection canonique de Dirac du résidu flottant machine ($\sim 10^{-16}$) sur le sous-espace $\sum m_i v_i = 0$.

### Step 4: Tests Unitaires Miroirs et Validation de Non-Régression
- **Impacted Elements**: [`tests/unit/potential/test_supernode_detector.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/potential/test_supernode_detector.py), [`tests/unit/potential/test_screening.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/potential/test_screening.py), [`tests/unit/solvers/test_supernode.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/tests/unit/solvers/test_supernode.py).
- **Actions & Implementation**:
  - Ajout du test `test_schmitt_trigger_hysteresis_anti_chattering` validant les 4 états de la bascule de Schmitt.
  - Ajout du test `test_c8_centered_spline_precision_and_symmetry` validant la valeur exacte $\Theta(1.5)=0.5$, l'anti-symétrie, et le pic de dérivée $3.338470458984375$.
  - Ajout du test `test_supernode_unfolded_sundman_adaptive_stepping_unsoftened` validant l'intégration adaptative sous horloge dépliée, zéro adoucisseur et réversibilité temporelle.
  - Exécution complète de `pytest tests/unit/` : 132/132 tests réussis (100 %).

### Step 5: Validation Git Commit
- **Impacted Elements**: Répertoire Git.
- **Actions & Implementation**:
  - Indexation des 9 fichiers impactés et commit `a5a5d25` : `feat(potential,solvers): remediate Chapter 4 supernodes and screening - Phase 2 of Roadmap VI (Iteration 145)`.

---

## Architecture and Code Organization (Evolution):
```
Couche 3 : Potentiel Conormal et Détection Régularisée
 ├── constants.py            : ETA_ENTER_HILL = 0.60, ETA_EXIT_HILL = 0.75
 ├── supernode_detector.py   : Bascule de Schmitt anti-chattering [0.60, 0.75]
 ├── screening.py            : Spline C^8 centrée en xi = s - 0.5 (kappa ~ 1)
Couche 2 : Solveurs Analytiques
 └── solvers/supernode.py    : Horloge dépliée dt_S = rho_S dtau_S, zero adoucisseur (d_safe = 0)
Couche 5 : Système N-Corps
 └── system/nbody_system.py  : Persistance active_supernodes à travers les pas macro
```

## Dependencies and Configuration:
- Aucune dépendance externe ajoutée.

## Commit Details:
- **Identifier Phase 2**: `a5a5d25` (`feat(potential,solvers): remediate Chapter 4 supernodes and screening - Phase 2 of Roadmap VI (Iteration 145)`)

## Design Choices, Trade-offs and Warnings:
- **Conditionnement de la Spline $C^8$** : L'utilisation de la forme centrée en $\xi = s - 0.5$ ramène le conditionnement de $\kappa = 1902$ à $\kappa \approx 1$ sans introduire de saut à la frontière $\xi = 0.25$ ($|\Delta\Theta| = 0.00$).
- **Projection Canonique de Dirac vs Recentrage Artificiel** : La projection de Dirac ne traite que le résidu machine flottant ($\sim 10^{-16}$) sans modifier la dynamique physique newtonienne.

---

## Verification:
### Acceptance Verification:
- [x] **Critère Step 1** : `detect_supernodes` filtre les oscillations de marée et préserve l'état actif dans la bande $[0.60, 0.75]$ sans chattering.
- [x] **Critère Step 2** : La spline $C^8$ présente un conditionnement $\kappa \approx 1$ avec anti-symétrie exacte, dérivée de crête exacte ($3.338470458984375$), et stricte monotonie.
- [x] **Critère Step 3** : Le solveur `supernode.py` opère sans adoucisseur ($d_{\mathrm{safe}} \equiv 0$) sous horloge dépliée de Sundman et conserve l'impulsion à $< 10^{-14}$.
- [x] **Critère Step 4** : 132/132 tests unitaires validés (100 % pass rate).
- [x] **Critère Step 5** : Commit Git `a5a5d25` tracé et enregistré.

---

## Next Steps:
- Mettre à jour `notes.md` et `roadmap.md` pour marquer la Phase 2 comme terminée (`[x]`).
- Engager la **Phase 3 de la Roadmap VI (Itération 146)** : Remédiation de la Synchronisation de Skeel-Gear et de l'Échéancier Dyadique du Chapitre 5 (Confinement des transitions d'atlas strictement aux frontières macro, profondeur dyadique dynamique $k_{\max}=17$, et fusion FSAL des demi-kicks dans BM4).
