# Walkthrough of Iteration 136 : Éradication des Seuils Hardcodés de Paliers Dyadiques, Régularité de Découpage Adaptée au Schéma et Conformité Monographique

## Summary:
L'itération 136 a restauré la conformité théorique et l'invariance d'échelle de la bibliothèque conormale conformément à la Section 4.2 et au Chapitre 5 de la monographie :
1. **Rétractation formelle et clarification épistémique** : Reconnaissance de l'inexistence du soi-disant « Lemme 4.3 » (hallucination antérieure) et réalignement sur les vrais théorèmes du manuscrit (Théorème 4.1 de régularité Lipschitz $\mathcal{C}^{1,1}$, Lemme 4.1 d'éclatement conormal, Théorème 4.2 d'annulation multipolaire).
2. **Éradication totale des seuils métriques hardcodés** : Suppression de tous les paliers arbitraires `prox_k` (`min_r <= 0.015`, etc.) dans `dyadic_scheduler.py`. Le niveau de raffinement $k$ est désormais déterminé de façon continue et invariante d'échelle par le ratio entre le pas macro et les échelles de temps caractéristiques physiques locales ($\tau_{\text{kin}} = r / v_{\text{rel}}$, $\tau_{\text{dyn}} = \sqrt{r^3 / GM}$, et anticipation péricentrique $r_{\text{peri}} = a(1-e)$).
3. **Hiérarchie de fonctions de découpage adaptées à l'ordre du schéma** : Implémentation dans `screening.py` de splines Hermite de degrés adaptés ($\mathcal{C}^2$ quintique pour Strang 2, $\mathcal{C}^4$ nonique pour BM4/Yoshida 4 avec dérivées 1 à 4 nulles aux bords, $\mathcal{C}^6$ tridécylique pour BM6 avec dérivées 1 à 6 nulles aux bords, et mollifieur de Cauchy $\mathcal{C}^\infty$), avec évaluation symétrique stable autour de $s=0.5$ éliminant l'annulation catastrophique en float64.
4. **Couplage automatique intégrateur $\leftrightarrow$ apex** : `NBodySystem` et les classes d'intégrateurs configurent automatiquement `ScreenedApexPotential` selon l'ordre du schéma.
5. **Harmonisation monopolaire** : Définition de `include_screened_monopole=False` par défaut dans `compute_apex_forces`, unifiant le calcul de la force de kick $V_{\text{inter}}$ entre `step_apex_kick` et `correctors.py`.
6. **Gain numérique certifié** : 131 tests unitaires validés (100%), et sur le système 4-corps hiérarchique, le saut d'erreur d'énergie de BM6 au passage péricentrique ($t=3.3000$) chute de $1.96 \times 10^{-7}$ à $1.12 \times 10^{-8}$, et l'erreur relative maximale sur l'ensemble de la trajectoire est réduite à **$1.14 \times 10^{-8}$** (gain supérieur à un facteur 17) à travers 15 transitions d'atlas sans dérive séculaire.

## Critical Discoveries:
- **Annulation catastrophique sur les splines de haut degré** : Lors de l'évaluation directe du polynôme de degré 13 pour $\mathcal{C}^6$ au voisinage de $s=1$, les puissances de $s$ et les grands coefficients alternés ($24024, 20020, 16380, \dots$) généraient une perte de précision en float64 conduisant à $\theta(u) > 1 + 1.88 \times 10^{-12}$. La résolution par symétrie exacte autour de $s=0.5$ ($\sigma(s) = 1 - \sigma(1-s)$ pour $s > 0.5$) a ramené la précision à $10^{-16}$ et garanti $0 \le \theta \le 1$ à précision machine.

## Realization:
- **Start Time**: 2026-10-08T08:18:00+02:00
- **Writing Time**: 2026-10-08T08:26:00+02:00

### Step 1: Rétractation Formelle et Rectification Épistémique
- **Impacted Elements**: `implementation_plan.md`, `task.md`, `notes.md`.
- **Description**: Retrait définitif de la citation hallucinée du « Lemme 4.3 ». Clarification de la taxonomie des résultats formels du Chapitre 4 de la monographie.

### Step 2: Refonte Physique Continue du Planificateur Dyadique
- **Impacted Elements**: [`system/dyadic_scheduler.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/system/dyadic_scheduler.py).
- **Description**: Élimination de `prox_k` et des seuils en distance `min_r <= 0.015`. Implémentation du calcul vectorisé de $\tau_{\text{kin}} = r / v_{\text{rel}}$, $\tau_{\text{dyn}} = \sqrt{r^3 / GM}$, et pour les paires en rapprochement ($v_{\text{rad}} < 0$), anticipation de $\tau_{\text{peri}} = \min(r_{\text{peri}}/v_{\text{peri}}, \sqrt{r_{\text{peri}}^3/GM})$. Le niveau dyadique découle continûment de $k = \lceil \log_2(\Delta t_{\text{macro}} / (\eta \tau_{\text{enc}})) \rceil$, respectant strictement l'invariance d'échelle sous homothéties.

### Step 3: Fonctions de Découpage Conormales Adaptées à l'Ordre du Schéma
- **Impacted Elements**: [`potential/screening.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screening.py), [`potential/__init__.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/__init__.py).
- **Description**: Implémentation de la fonction `get_screening_mode_for_order(order)` et enrichissement d'`eval_screening_theta` et `eval_screening_theta_vectorized` avec les profils :
  - `"c2"` (quintique, ordre 2) : $\sigma_2(s) = 10s^3 - 15s^4 + 6s^5$, $\sigma' = 30 s^2 (1-s)^2$.
  - `"c4"` (nonique, ordre 4) : $\sigma_4(s) = s^5 (126 - 420s + 540s^2 - 315s^3 + 70s^4)$, $\sigma' = 630 s^4 (1-s)^4$.
  - `"c6"` (tridécylique, ordre 6) : $\sigma_6(s) = s^7 (1716 - 9009s + 20020s^2 - 24024s^3 + 16380s^4 - 6006s^5 + 924s^6)$, $\sigma' = 12012 s^6 (1-s)^6$.
  - `"cauchy"` (mollifieur $\mathcal{C}^\infty$ Gevrey-2).
  Évaluation symétrique stable pour tous les splines garantissant la non-discontinuité de l'analyse rétrograde d'erreur jusqu'à l'ordre du schéma.

### Step 4: Couplage Automatique Intégrateur-Apex et Harmonisation Monopolaire
- **Impacted Elements**: 
  - [`potential/screened_apex.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/potential/screened_apex.py)
  - [`integrators/strang.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/strang.py)
  - [`integrators/yoshida.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/yoshida.py)
  - [`integrators/blanes_moan_4.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/blanes_moan_4.py)
  - [`integrators/blanes_moan_6.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/blanes_moan_6.py)
  - [`integrators/correctors.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/integrators/correctors.py)
  - [`system/nbody_system.py`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/system/nbody_system.py)
- **Description**: `ScreenedApexPotential` accepte désormais `order: Optional[int]` et résout son profil optimal. `NBodySystem` extrait l'ordre et le mode de son intégrateur pour initialiser son évaluateur d'apex. Harmonisation du défaut `include_screened_monopole=False` assurant que les dérivées directionnelles des correcteurs et les kicks d'apex calculent exactement le même champ $V_{\text{inter}}$.

### Step 5: Vérification Numérique et Benchmarking
- **Impacted Elements**: `tests/unit/`, `scratch/diagnose_4body_jump.py`.
- **Description**: Exécution complète de la suite unitaire avec 131 tests réussis sur 131 (100%). Exécution de la télémétrie fine 4-corps montrant l'abaissement de l'erreur d'énergie de BM6 à $1.14 \times 10^{-8}$.

## Architecture and Code Organization (Evolution):
```
[NBodySystem]
     │
     ├── integrator (order = 2, 4, 6)
     │        │
     │        ▼
     ├── apex = ScreenedApexPotential(order=integrator.order)
     │        │
     │        ▼
     │   eval_screening_theta(u, mode="c2" | "c4" | "c6" | "cauchy")
     │   (Hermite splines C^{2p} with vanished boundary derivatives)
     │
     └── scheduler = DyadicScheduler(dt_macro, G)
              │
              ▼
         compute_dynamical_level(state, dt)
         tau_enc = min(tau_kin, tau_dyn, tau_peri)
         k = ceil(log2(dt / (eta * tau_enc)))
         (Strictly continuous, scale-invariant, no hardcoded distances)
```

## Design Choices, Trade-offs and Warnings:
- **Évaluation symétrique des splines** : Pour $s > 0.5$, l'utilisation de $1.0 - \sigma(1 - s)$ supprime rigoureusement les risques de débordement $\theta > 1.0$ dus aux erreurs d'arrondi sur les grands monômes de degré 9 ou 13.
- **Préservation de l'hamiltonien modifié** : L'utilisation d'un spline $\mathcal{C}^4$ pour BM4 et $\mathcal{C}^6$ pour BM6 garantit que l'hamiltonien effectif de l'intégrateur ne subit aucune singularité différentielle lors des transitions d'apex.
- **Point d'attention pour les prochains agents** : Ne jamais réintroduire de table empirique de distances dans `dyadic_scheduler.py`. Tout raffinement dyadique doit demeurer une fonctionnelle purement cinématique et dynamique des ratios sans dimension $\Delta t / \tau$.

## Verification:
- **Suite de tests unitaires** : `pytest tests/unit/` $\to$ **131 tests passés avec succès** (100% en 6.15s).
- **Télémétrie 4-corps hiérarchique (`blanes_moan_6`)** :
  - Erreur relative maximale : **$1.1415 \times 10^{-8}$** (contre $1.9602 \times 10^{-7}$ avant l'itération).
  - Saut d'énergie au pas péricentrique critique (pas 220, $t=3.3000$) : **$1.1230 \times 10^{-8}$** (contre $1.9597 \times 10^{-7}$ avant l'itération, gain de $17.5\times$).
  - Niveau dyadique sélectionné automatiquement : niveau 10 (1024 sous-pas, $\Delta t_{\text{sub}} = 1.46 \times 10^{-5}$ s).
  - Transitions d'atlas $\mathrm{SL}(2, \mathbb{R})$ : 15 transitions parfaitement stables et continues sans chattering.

## Next Steps:
- Mettre à jour les planches graphiques haute résolution (300 DPI) et le fichier de données de la vitrine interactive (`showcase_data.json`) pour refléter cette précision de $1.14 \times 10^{-8}$ sur le 4-corps hiérarchique.
- Explorer la formulation des correcteurs symplectiques $C_4$ en coordonnées de Jacobi relatives afin d'éviter l'amortissement d'amplitude lors des sous-pas dyadiques très fins ($h \to 0$).
