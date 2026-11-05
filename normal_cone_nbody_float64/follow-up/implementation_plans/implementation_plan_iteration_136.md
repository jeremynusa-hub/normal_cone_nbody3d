# Implementation Plan of Iteration 136 : Éradication des Seuils Hardcodés de Paliers Dyadiques, Régularité de Découpage Adaptée au Schéma et Conformité Monographique
Created at: 2026-10-08T08:18:00+02:00

## Context:
1. **Explicit Request**: 
   - L'utilisateur exige une relecture rigoureuse de la monographie pour identifier et éradiquer les erreurs conceptuelles critiques.
   - Suppression catégorique de tout seuil hardcodé sur les paliers dyadiques (ex: `prox_k = 9 if min_r <= 0.015`, etc.).
   - Rétractation formelle de la citation hallucinée d'un soi-disant "Lemme 4.3" (inexistant dans le manuscrit).
   - Audit et traitement de la fonction de régularité conormale $\Theta(u_v)$ : la théorie démontre une régularité $\mathcal{C}^\infty$ (mollifieur de Cauchy Gevrey-2), mais sur le plan computationnel, elle doit être adaptée à l'ordre de chaque schéma numérique ($\mathcal{C}^2$ pour Strang 2, $\mathcal{C}^4$ pour BM4, $\mathcal{C}^6$ ou lisse pour BM6) afin de préserver l'analyse rétrograde d'erreur sans tronquer l'ordre symplectique.
2. **Implicit Strategic Value**:
   - Rétablir l'invariance d'échelle absolue de la mécanique céleste sans régulateur métrique arbitraire ($d_0 = 0, \varepsilon = 0$).
   - Préserver l'exactitude des ordres symplectiques effectifs (ordre 4 pour BM4, ordre 6 pour BM6) en éliminant les sauts de dérivées dans l'anneau de transition $1 < u < 2$.
   - Éliminer le saut d'erreur d'énergie constaté lors du péricentre serré du système 4-corps hiérarchique ($t=3.3000$).
3. **Alignment Criteria**:
   - Absence totale de seuils métriques absolus dans `dyadic_scheduler.py` : le niveau dyadique $k$ découle continûment et dynamiquement du ratio entre le pas macro et les temps caractéristiques physiques ($\tau_{\text{kin}}, \tau_{\text{dyn}}$).
   - Disponibilité dans `potential/screening.py` de profils de transition polynomiaux et lisses certifiés pour chaque ordre de schéma ($\mathcal{C}^2$ quintique, $\mathcal{C}^4$ nonique, $\mathcal{C}^6$ tridécylique, et $\mathcal{C}^\infty$ Cauchy), automatiquement sélectionnés par l'intégrateur associé.
   - Suppression du résidu de saut d'énergie artificiel lors des frôlements péricentriques dans le benchmark 4-corps.

## Guidelines:
- **Method**: Respect strict de la géométrie de Fulton-MacPherson et des formulations de la Section 4.2 de la monographie.
- **Red Lines**: 
  - Interdiction absolue de réintroduire tout seuil de distance hardcodé (`if r < 0.015 ...`).
  - Interdiction absolue d'adoucir la physique ($d_0 = 0, \varepsilon = 0$).
  - Ne jamais modifier les fichiers LaTeX sous `redaction/unsoftened_conormal_nbody3d/`.
- **Architecture**:
  - `system/dyadic_scheduler.py` : calcul continu $k = \lceil \log_2(\Delta t / \Delta t_{\text{dyn}}) \rceil$.
  - `potential/screening.py` : profils de régularité `order_2` ($\mathcal{C}^2$), `order_4` ($\mathcal{C}^4$), `order_6` ($\mathcal{C}^6$), et `cauchy` ($\mathcal{C}^\infty$).
  - `potential/screened_apex.py` : couplage transparent avec le mode de régularité requis.
  - `integrators/correctors.py` : cohérence stricte de la partition d'énergie avec `screened_apex.py`.

## Plan:
### Step 1: Rétractation Formelle et Rectification Épistémique
- **Context**: Mention erronée d'un "Lemme 4.3" dans les échanges passés.
- **Objective**: Clarifier la structure exacte du Chapitre 4 de la monographie : Théorème 4.1 (Régularité Lipschitz $\mathcal{C}^{1,1}$ universelle), Lemme 4.1 (Identité d'éclatement conormale), Définition 4.1 & Théorème 4.2 (Annulation multipolaire).
- **Deliverables**:
  - `MonographCitationCorrection`: Rectification formelle consignée dans le plan et les notes.
- **Process**: Vérification directe de l'arborescence et des labels LaTeX du Chapitre 4.

### Step 2: Refonte Physique Continue du Planificateur Dyadique
- **Context**: `normal_cone_nbody_float64/system/dyadic_scheduler.py` contenait des paliers arbitraires `prox_k` hardcodés par plages métriques `min_r <= 0.015`, violant l'invariance d'échelle.
- **Objective**: Éradiquer ces paliers et calculer le niveau dyadique $k$ exclusivement à partir des temps caractéristiques dynamiques et cinématiques invariants d'échelle :
  $\tau_{\text{kin}} = r_{ij} / v_{\text{rel}, ij}$, $\tau_{\text{dyn}} = \sqrt{r_{ij}^3 / (G (m_i + m_j))}$, et anticipation képlérienne du temps jusqu'au péricentre.
- **Deliverables**:
  - `ScaleInvariantDyadicScheduler`: Méthode `compute_dynamical_level` refondue sans aucun seuil métrique.
- **Process**: Remplacement des blocs `if min_r <= 0.015...` par la formule universelle continue $k = \lceil \log_2(\Delta t_{\text{macro}} / (\eta \tau_{\text{enc}})) \rceil$.

### Step 3: Fonctions de Découpage Conormales Adaptées à l'Ordre du Schéma
- **Context**: La monographie stipule que le mollifieur $\mathcal{C}^\infty$ garantit la théorie, mais qu'un proxy polynomial doit avoir une régularité matching l'ordre du schéma pour que l'analyse rétrograde d'erreur (modified Hamiltonian) ne soit pas tronquée. L'implémentation actuelle ne proposait qu'un spline quintique $\mathcal{C}^2$ (saut de dérivée troisième $\sigma''' = 60 \ne 0$) et un mollifieur transcendental.
- **Objective**: Équiper `potential/screening.py` de splines Hermite de degrés adaptés :
  - Ordre 2 ($\mathcal{C}^2$, degré 5) : pour Strang $\mathcal{S}_2$.
  - Ordre 4 ($\mathcal{C}^4$, degré 9) : $\sigma_4(s) = s^5 (126 - 420 s + 540 s^2 - 315 s^3 + 70 s^4)$, dérivées 1 à 4 nulles aux bords, pour BM4 et Yoshida 4.
  - Ordre 6 ($\mathcal{C}^6$, degré 13) : $\sigma_6(s)$ avec dérivées 1 à 6 nulles aux bords, ou mollifieur $\mathcal{C}^\infty$, pour BM6.
- **Deliverables**:
  - `OrderAdaptedScreening`: Fonctions `eval_screening_theta` et version vectorisée supportant `"c2"`, `"c4"`, `"c6"` et `"cauchy"`.
- **Process**: Implémentation analytique exacte, dérivées premières continues, et tests unitaires de régularité aux frontières $u=1$ et $u=2$.

### Step 4: Couplage Automatique Intégrateur-Apex et Harmonisation Monopolaire
- **Context**: Chaque intégrateur doit automatiquement transmettre son exigence de régularité à son évaluateur d'apex, et la gestion du monopôle écranté doit être parfaitement cohérente entre `screened_apex.py` et les correcteurs.
- **Objective**: Assurer que BM4 active le spline $\mathcal{C}^4$, BM6 active $\mathcal{C}^6$ ou $\mathcal{C}^\infty$, et Strang active $\mathcal{C}^2$. Harmoniser les valeurs par défaut de `include_screened_monopole`.
- **Deliverables**:
  - `IntegratorScreeningCoupling`: Initialisation harmonisée dans `NBodySystem`, `StrangIntegrator`, `BlanesMoan4Integrator`, `BlanesMoan6Integrator`.
- **Process**: Édition via `replace_file_content` des classes d'intégrateurs et de `NBodySystem`.

### Step 5: Vérification Numérique et Benchmarking
- **Context**: Le système 4-corps hiérarchique présentait un saut d'énergie de $+1.96 \times 10^{-7}$ au pas 220 sous BM6.
- **Objective**: Vérifier que l'intégration du système 4-corps avec le planificateur physique continu et le découpage $\mathcal{C}^6$/Cauchy élimine le saut artificiel et maintient la stabilité symplectique.
- **Deliverables**:
  - `DiagnosticVerificationRun`: Exécution de `scratch/diagnose_4body_jump.py` et validation de la suite unitaire.
- **Process**: Lancement des tests unitaires et de la télémétrie pas-à-pas.

## Acceptance Criteria:
- `ScaleInvariance`: `grep -rn "prox_k =" normal_cone_nbody_float64/` ne renvoie aucun résultat. Aucun seuil hardcodé en distance n'existe dans le scheduler.
- `SmoothC4C6Splines`: Les fonctions de découpage $\mathcal{C}^4$ et $\mathcal{C}^6$ satisfont $|\sigma(1)| < 10^{-15}$, $|\sigma(2) - 1| < 10^{-15}$, $|\sigma'(1)| < 10^{-15}$, $|\sigma'(2)| < 10^{-15}$, et les dérivées successives requises s'annulent rigoureusement aux frontières.
- `SymplecticOrderIntegrity`: BM4 utilise au minimum un profil $\mathcal{C}^4$, BM6 utilise un profil $\mathcal{C}^6$ ou $\mathcal{C}^\infty$.
- `TestSuitePass`: Tous les tests unitaires sous `tests/unit/` s'exécutent avec succès (`pytest tests/unit/`).

## Observations:
> [!IMPORTANT]
> L'invariance d'échelle est le pilier fondamental de la géométrie de Fulton-MacPherson : toute dépendance à une unité de longueur externe (comme un seuil à $0.015$) rompt le groupe d'homothétie de la mécanique céleste.
> [!NOTE]
> La transition d'un spline $\mathcal{C}^2$ à un spline $\mathcal{C}^4$ ou $\mathcal{C}^6$ rétablit la validité de l'hamiltonien modifié dans le théorème de Baker-Campbell-Hausdorff jusqu'à l'ordre du schéma considéré.
