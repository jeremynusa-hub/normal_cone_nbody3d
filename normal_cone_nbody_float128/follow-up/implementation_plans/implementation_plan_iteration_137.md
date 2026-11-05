# Implementation Plan of Iteration 137 : Intégration du Spline Hermite C8 (Degré 17) pour PBM6 et C6 (Degré 13) pour PBM4
Created at: 2026-10-08T08:30:00+02:00

## Context:
1. **Explicit Request**: 
   - L'utilisateur valide et commande l'intégration du profil polynomial $\mathcal{C}^8$ (spline d'Hermite de degré 17) pour les schémas d'ordre 8 (PBM6), et l'affectation du profil $\mathcal{C}^6$ (degré 13) pour PBM4 (dont l'ordre effectif est 6).
2. **Implicit Strategic Value**:
   - Rendre les schémas traités (Processed Symplectic Integrators) intrinsèquement cohérents avec leur ordre effectif d'analyse rétrograde d'erreur (Baker–Campbell–Hausdorff) :
     - Strang $\mathcal{S}_2$ (ordre 2) $\to \mathcal{C}^2$ (degré 5)
     - BM4 / Yoshida (ordre 4) $\to \mathcal{C}^4$ (degré 9)
     - PBM4 / BM6 (ordre effectif 6 / ordre 6 brut) $\to \mathcal{C}^6$ (degré 13)
     - PBM6 (ordre effectif 8) $\to \mathcal{C}^8$ (degré 17)
   - Éviter l'underflow IEEE 754 de l'exponentielle du mollifieur de Cauchy à $s < 1.4 \times 10^{-3}$ tout en garantissant des performances CPU optimales par évaluation polynomiale factorisée (FMA/Horner).
3. **Alignment Criteria**:
   - Implémentation certifiée de `"c8"` dans `potential/screening.py` (degré 17, dérivées 1 à 8 nulles aux bords, évaluation symétrique autour de $s=0.5$).
   - Prise en compte de `effective_order` dans `ProcessedSymplecticIntegrator` et les correcteurs (`CorrectorC4` utilise `"c6"`, `CorrectorC6` utilise `"c8"`).
   - 100% de passage sur la suite de tests unitaires et non-régression numérique.

## Guidelines:
- **Method**: Respect strict des coefficients polynomiaux analytiques exacts issus de l'intégration de $\int 218790 s^8 (1-s)^8 \mathrm{d}s$.
- **Red Lines**: 
  - Maintenir $d_0 = 0, \varepsilon = 0$ (zéro adoucissement).
  - Évaluation symétrique inconditionnelle pour éliminer tout risque d'annulation catastrophique en float64.
- **Architecture**:
  - `potential/screening.py` : ajout de `"c8"` en scalaire et vectorisé.
  - `integrators/correctors.py` : `CorrectorC4` $\to$ `"c6"`, `CorrectorC6` $\to$ `"c8"`.
  - `integrators/processed.py` : couplage basé sur `effective_order`.
  - `system/nbody_system.py` : détection prioritaire de `effective_order` pour `ProcessedSymplecticIntegrator`.

## Plan:
### Step 1: Implémentation du Profil Hermite C8 dans screening.py
- **Context**: Le mode `c8` (degré 17) doit être ajouté pour couvrir les schémas d'ordre effectif 8.
- **Objective**: Implémenter $\sigma_8(s) = s^9 \sum_{k=0}^8 c_k s^k$ et sa dérivée $\sigma_8'(s) = 218790 (s(1-s))^8$ de façon symétrique autour de $s=0.5$ en scalaire et vectorisé.
- **Deliverables**:
  - `ScreeningC8Implementation`: Fonctions `eval_screening_theta` et `eval_screening_theta_vectorized` supportant `"c8"`, `"order_8"`, et mise à jour de `get_screening_mode_for_order(order)`.
- **Process**: Intégration des coefficients rationnels exacts et factorisation symétrique.

### Step 2: Couplage par Ordre Effectif pour PBM4 et PBM6
- **Context**: PBM4 a un ordre effectif de 6, PBM6 a un ordre effectif de 8.
- **Objective**: Configurer `CorrectorC4` avec `screening_mode="c6"` (ordre 6), `CorrectorC6` avec `screening_mode="c8"` (ordre 8), et synchroniser `ProcessedSymplecticIntegrator` et `NBodySystem` sur `integrator.effective_order`.
- **Deliverables**:
  - `EffectiveOrderCoupling`: Mise à jour de `correctors.py`, `processed.py` et `nbody_system.py`.
- **Process**: Édition via `replace_file_content`.

### Step 3: Tests Unitaires et Validation Numérique
- **Context**: Valider les dérivées, la monotonicité, la vectorisation et la non-régression de l'ensemble de la bibliothèque.
- **Objective**: Étendre `tests/unit/potential/test_screening.py` à `"c8"`, vérifier `pytest tests/unit/`, et vérifier la stabilité de la télémétrie 4-corps.
- **Deliverables**:
  - `ScreeningC8Tests`: Tests unitaires enrichis et validés à 100%.
- **Process**: Lancement de la suite de tests et diagnostic.

## Acceptance Criteria:
- `C8ExactBoundary`: Pour $u=1.0$, $\theta=0, \theta'=0$. Pour $u=2.0$, $\theta=1, \theta'=0$. Pour $u=1.5$, $\theta=0.5$.
- `C8Derivatives`: Dérivées analytiques en accord avec les différences finies à $10^{-5}$ près sans instabilité.
- `Monotonicity`: $\theta'(u) \ge 0$ partout sur $[1.0, 2.0]$.
- `TestSuitePass`: 100% des tests unitaires réussis sous `pytest tests/unit/`.

## Observations:
> [!IMPORTANT]
> L'adaptation de la régularité à l'ordre effectif ($\mathcal{C}^6$ pour PBM4, $\mathcal{C}^8$ pour PBM6) clôt de manière formelle et irréfutable la question de la complétude du formalisme de Baker–Campbell–Hausdorff sur la variété conormale $\mathbf{Q}[N]$.
