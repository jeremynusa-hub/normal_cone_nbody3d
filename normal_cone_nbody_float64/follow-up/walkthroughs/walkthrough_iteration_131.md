# Walkthrough of Iteration 131 : Audit Contradictoire Approfondi et Remédiation des Correcteurs Symplectiques Traités (PBM) sur les Systèmes Multi-Corps à Transitions d'Atlas

## Summary:
L'Itération 131 a élucidé et résolu l'anomalie de saut d'énergie et de divergence (`nan` / $+14\,333$) observée sur l'intégrateur traité Blanes-Moan 4 (PBM4) dans le scénario à 4 corps hiérarchique (`run_hierarchical_4body_nested_supernodes`).
Deux audits contradictoires indépendants (Auditeur Théorique du Manuscrit LaTeX et Profileur Numérique de Flot d'Exécution) ont été déployés pour diagnostiquer la mécanique des correcteurs symplectiques $C_{2p}(h)$ et du télescopage dyadique de frontière $\Psi = C \circ S^M \circ C^{-1}$.

L'audit contradictoire a révélé le mécanisme exact du décrochage :
1. Dans `NBodySystem.step()`, la branche dédiée aux schémas traités déléguait l'intégration dyadique à `ProcessedSymplecticIntegrator.integrate()`, qui **gelait la topologie de l'arbre de Jacobi (`tree`) sur l'ensemble des sous-pas dyadiques** (jusqu'à 512 sous-pas). À l'inverse, les intégrateurs de base (BM4, BM6) réévaluaient l'adaptation de l'arbre à chaque sous-pas.
2. Au pas 412 ($t \approx 6.1838$), les corps 1 et 2 plongeaient dans un passage péricentrique rasant ($r_{12} \approx 8.7 \times 10^{-4}$). L'arbre de Jacobi n'ayant pas muté, l'interaction $(1, 2)$ était calculée via la force d'apex non régularisée $1/r^2 \sim 1.3 \times 10^6$ au lieu d'être résolue analytiquement par le solveur universel de Stumpff-Sundman. Cette omission injectait un kick gigantesque propulsant l'énergie à $+45\,815$ et causant un débordement `nan`.
3. De plus, lors d'une transition d'arbre $\mathcal{T}_A \to \mathcal{T}_B$, le gradient d'apex $F_{\text{apex}}$ subit un saut discontinu fini correspondant au basculement des paires feuilles masquées ($F_{\text{apex}}(\mathcal{T}_B) - F_{\text{apex}}(\mathcal{T}_A) = F_{\text{Newton}, 01} - F_{\text{Newton}, 12}$). Appliquer $C(\mathcal{T}_B)$ en sortie d'un bloc entré avec $C^{-1}(\mathcal{T}_A)$ brisait le télescopage frontière $C \circ C^{-1} = \mathrm{Id}$.

Après implémentation d'une **adaptation topologique continue à chaque sous-pas dyadique** et d'une **fermeture de télescopage frontière sécurisée sur la carte d'entrée** (`corr_tree = tree_enter if had_transition else self.tree`), PBM4 s'exécute avec une régularité absolue :
- **PBM4 sur le cas à 4 corps** : l'erreur relative chute de `nan` / $14\,333$ à **$1.370 \times 10^{-6}$**, rejoignant exactement la précision de BM4 ($1.370 \times 10^{-6}$).
- **PBM6 sur le cas à 4 corps** : atteint **$1.968 \times 10^{-7}$**, rejoignant exactement BM6 ($1.968 \times 10^{-7}$).
- L'ensemble des 6 schémas sur le 4-corps est parfaitement hiérarchisé : Strang ($4.8 \times 10^{-3}$) $\gg$ Yoshida ($4.0 \times 10^{-4}$) $\gg$ BM4/PBM4 ($1.37 \times 10^{-6}$) $\gg$ BM6/PBM6 ($1.97 \times 10^{-7}$).
- Les 125 tests unitaires passent à 100% en 11.95s. La vitrine interactive et les planches 300 DPI ont été intégralement resynchronisées.

---

## Critical Discoveries:
1. **Désynchronisation de l'Atlas dans `integrate()`** :
   - Pour les schémas non traités (Strang, Yoshida, BM4, BM6), `NBodySystem.step()` exécutait :
     ```python
     for _ in range(num_substeps):
         self.tree, _ = self.atlas_manager.evaluate_tree_adaptation(curr_state, self.tree)
         curr_state = self.integrator.step(curr_state, self.tree, dt=dt_substep, ...)
     ```
   - Pour `ProcessedSymplecticIntegrator`, `integrate()` exécutait une boucle interne `for _ in range(num_steps): curr_state = self.base_integrator.step(curr_state, tree, ...)` où `tree` demeurait figé sur la topologie du début du pas macro.
   - En conséquence, lors d'approches serrées rapides ($r_{12} \to 8.7 \times 10^{-4}$ au pas 412), la paire en collision n'était pas promue feuille analytique et subissait l'explosion en $1/r^2$.
2. **Discontinuité de l'Opérateur de Correction à Travers les Transitions $\mathrm{SL}(2, \mathbb{R})$** :
   - L'opérateur de correction $C_{2p}(h) = \exp(-h^{2p} X_K)$ dépend des forces d'apex $F_{\text{apex}}(q; \mathcal{T})$.
   - Lorsque l'arbre de Jacobi bascule de $\mathcal{T}_A$ (paire feuille $(0, 1)$) à $\mathcal{T}_B$ (paire feuille $(1, 2)$), la force d'apex change de manière discontinue car les paires masquées diffèrent :
     $$F_{\text{apex}}(q; \mathcal{T}_B) - F_{\text{apex}}(q; \mathcal{T}_A) = F_{\text{Newton}, 01} - F_{\text{Newton}, 12} \sim \mathcal{O}\left(\frac{G}{r^2}\right).$$
   - Si le macro-pas commence avec $C^{-1}(\mathcal{T}_A)$ et se termine avec $C(\mathcal{T}_B)$, l'opérateur net traversé est :
     $$C(\mathcal{T}_B) \circ C^{-1}(\mathcal{T}_A) \neq \mathrm{Id},$$
     injectant un déplacement spatial discontinu non-physique $\Delta q \propto \kappa_4 h^4 G / r^2$.
3. **Principe de Clôture Sécurisée sur la Carte d'Entrée** :
   - Lorsqu'une transition d'atlas survient au sein d'un pas macro, clôturer le bloc avec l'arbre d'entrée $\mathcal{T}_A$ (`corr_tree = tree_enter`) garantit que $C(\mathcal{T}_A) \circ C^{-1}(\mathcal{T}_A) = \mathrm{Id}$ s'annule identiquement à précision machine ($< 10^{-15}$).
   - La trajectoire physique est préservée sans perturbation parasite de jauge, tout en bénéficiant de l'adaptation sous-pas dynamique de l'arbre.

---

## Realization:
- **Started At**: 2026-10-07T20:50:00+02:00
- **Completed At**: 2026-10-07T21:24:00+02:00

### Step 1: Déploiement de 2 Audits Contradictoires Spécialisés
- **Impacted Elements**: `redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/`, `normal_cone_nbody_float64/integrators/correctors.py`, `system/nbody_system.py`
- **Actions & Implementation**:
  - Déploiement de 2 sous-agents :
    - `c0a80636` (Théoricien du Manuscrit) : confrontation avec la Section 5.5, analyse des équations homologiques $\{H, K_{2p}\} = H_{2p+1}$, vérification des certificats Lean 4 (`01_CorrectorBM4_Theorems.lean`).
    - `ffb731e8` (Profileur Numérique) : traçage chronologique du crash au pas 412 ($t \approx 6.1838$) et démonstration empirique de la résolution par adaptation sous-pas.

### Step 2: Synthèse Contradictoire et Diagnostic Formel de la Pathologie
- **Actions & Implementation**:
  - Confrontation des deux rapports établissant la cause double : gel de l'arbre dans la boucle interne de sous-pas et saut de force d'apex aux transitions de cartes.

### Step 3: Conception et Implémentation du Correctif
- **Impacted Elements**: [`normal_cone_nbody_float64/system/nbody_system.py`](file:///C:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/system/nbody_system.py#L132-L165)
- **Actions & Implementation**:
  - Modification de `NBodySystem.step()` via `replace_file_content` :
    - Évaluation active de `self.atlas_manager.evaluate_tree_adaptation` à **chaque sous-pas dyadique** pour les schémas traités PBM.
    - Suivi du drapeau `had_transition`.
    - Sélection de l'arbre de sortie `corr_tree = tree_enter if had_transition else self.tree` pour garantir $C \circ C^{-1} = \mathrm{Id}$ en cas de bifurcation topologique.

### Step 4: Contrôle Qualité, Validation Numérique et Walkthrough
- **Impacted Elements**: `scratch/run_showcase_simulations.py`, `scratch/showcase_data.json`, `scratch/sync_showcase_html.py`, `interactive_conormal_showcase.html`, `brain/`
- **Actions & Implementation**:
  - Exécution complète des 125 tests unitaires : 100% de succès en 11.95s.
  - Régénération complète des 4 scénarios de simulation astrophysique (`run_showcase_simulations.py`) :
    - Figure-8 : PBM4 à $6.786 \times 10^{-12}$, BM6/PBM6 à $3.161 \times 10^{-12}$.
    - Burrau : BM4 à $1.324 \times 10^{-10}$, PBM4 à $1.527 \times 10^{-10}$, BM6/PBM6 à $1.112 \times 10^{-10}$.
    - Collision radiale avec rebond McGehee : erreur $3.039 \times 10^{-11}$.
    - 4-corps hiérarchique : PBM4 à **$1.370 \times 10^{-6}$**, PBM6 à **$1.968 \times 10^{-7}$**.
  - Synchronisation de `showcase_data.json` (5.32 Mo) et régénération des planches 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`).

---

## Design Choices, Trade-offs and Warnings:
- **Clôture Symplectique Neutre lors des Changements d'Arbre** :
  Lorsqu'une transition d'atlas $\mathrm{SL}(2, \mathbb{R})$ a lieu, utiliser `tree_enter` pour le correcteur de sortie annule exactement la transformation de jauge initiale ($C_{\text{enter}} \circ C_{\text{enter}}^{-1} = \mathrm{Id}$). Le pas se réduit alors à une intégration symplectique d'ordre 4 (ou 6) pur, sans discontinuité spatiale.
- **Zéro Adoucissement Strict** :
  $\varepsilon = 0$ et $d_0 = 0$ restent strictement observés sur l'ensemble de la bibliothèque et des scénarios.

---

## Verification:
### Acceptance Verification:
- [x] Déploiement et synthèse de deux audits contradictoires indépendants.
- [x] Diagnostic de l'origine du saut de PBM4 (gel de l'arbre et saut discontinu de $F_{\text{apex}}$).
- [x] Implémentation du correctif dans `system/nbody_system.py`.
- [x] Validation numérique : PBM4 chute de `nan` / $14\,333$ à **$1.370 \times 10^{-6}$**, encadré par BM4 et BM6.
- [x] 125/125 tests unitaires passent à 100% en 11.95s.
- [x] `showcase_data.json`, planches 300 DPI et `interactive_conormal_showcase.html` resynchronisés.
- [x] `notes.md` mis à jour de manière persistante et append-only.

---

## Next Steps:
- Formaliser la théorie de la transformation de jauge cotangente continue pour les correcteurs d'ordre supérieur à travers les bifurcations d'arbres $\mathrm{SL}(2, \mathbb{R})$.
- Évaluer le bénéfice d'un correcteur d'ordre 8 complet sur les orbites quasi-périodiques à $N=4$ sans transition d'arbre.
