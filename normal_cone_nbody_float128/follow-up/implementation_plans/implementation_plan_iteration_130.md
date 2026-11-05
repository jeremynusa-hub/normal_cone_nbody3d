# Implementation Plan of Iteration 130 : Éradication Définitive du Saut d'Énergie à 10⁻⁸ par Anticipation Képlérienne du Péricentre dans l'Ordonnanceur Dyadique
Created at: 2026-10-07T20:35:00+02:00

## Context:
1. **Explicit Request**:
   - « Je vois toujours des sauts dans la web app a 10-8 »
2. **Implicit Strategic Value**:
   - Atteindre une exactitude et une régularité numérique absolue sur toutes les familles d'intégrateurs symplectiques (Strang S2, Yoshida S4, Blanes-Moan BM4, PBM4, BM6, PBM6) sans discontinuité ni saut d'énergie artificiel lors des frôlements péricentriques extrêmes ($r < 0.01$).
   - Éliminer le défaut de réactivité de l'ordonnanceur dyadique en le dotant d'une capacité prédictive balistique et képlérienne (anticipation de la distance au péricentre $r_{\text{peri}} = a(1-e)$ et projection linéaire $r(t) + v_{\text{rad}} \Delta t$).
3. **Alignment Criteria**:
   - L'erreur relative d'énergie maximale pour BM4 et Yoshida sur le scénario de Burrau chute strictement sous $10^{-9}$ (suppression intégrale du saut à $10^{-8}$).
   - BM6 et PBM6 maintiennent une conservation d'énergie optimale sous $2 \times 10^{-10}$.
   - 100% des tests unitaires passent (`pytest tests/unit/`).
   - `showcase_data.json`, les planches 300 DPI et l'application interactive `interactive_conormal_showcase.html` sont régénérés et synchronisés.

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des scénarios.
- **Règle d'Édition** : Utiliser exclusivement `replace_file_content` pour les modifications du code Python, ne jamais écraser avec `overwrite: true`.

## Plan:
### Step 1: Diagnostic Numérique Fin du Saut à 10⁻⁸
- Identifier la cause précise du saut à $10^{-8}$ résiduel sur BM4 ($9.28 \times 10^{-9}$) et Yoshida ($6.70 \times 10^{-8}$) à $t \approx 8.76$ dans le problème de Burrau.
- Quantifier le délai de réponse de l'ordonnanceur dyadique (`compute_dynamical_level`) lors d'un effondrement radial rapide ($v_{\text{rad}} < 0$).

### Step 2: Conception et Implémentation de l'Anticipation Képlérienne du Péricentre
- Enrichir `normal_cone_nbody_float64/system/dyadic_scheduler.py` avec la projection prédictive balistique et képlérienne à deux corps.
- Adapter la profondeur dyadique $k$ dès l'entrée dans le cône d'approche ($r_{\text{projected}} \le 0.015 \implies k=9$, $r_{\text{projected}} \le 0.03 \implies k=7$) avant le franchissement du péricentre.

### Step 3: Contrôle Qualité, Régénération des Simulations et Synchronisation
- Vérifier la régression numérique sur Burrau : chute de BM4 à $\sim 1.3 \times 10^{-10}$ et de Yoshida à $\sim 3.0 \times 10^{-10}$.
- Exécuter la suite unitaire complète (`pytest tests/unit/`).
- Régénérer `scratch/showcase_data.json` et les planches 300 DPI, et synchroniser `interactive_conormal_showcase.html`.

## Acceptance Criteria:
- **Critère 1** : L'erreur d'énergie sur le problème de Burrau ne présente plus aucun saut à $10^{-8}$ ; BM4 et Yoshida atteignent une précision inférieure à $5 \times 10^{-10}$.
- **Critère 2** : Aucun saut de troncature d'apex n'est constaté au passage au péricentre $r_{\text{peri}} \approx 0.0085$.
- **Critère 3** : 100% des tests unitaires passent sans régression (`pytest tests/unit/`).

## Observations:
> [!NOTE]
> Le solveur analytique de Stumpff-Sundman résout exactement la dynamique képlérienne à deux corps. Le saut à $10^{-8}$ n'était pas une erreur du solveur analytique lui-même, mais une sous-résolution du macro-pas dans le kick d'apex d'arrière-plan ($F_{\text{apex}}$) lorsque le couple plongeait de $r=0.089$ à $r=0.0085$ au sein d'un pas macro sans raffinement prédictif préalable.
