# Implementation Plan of Iteration 156 : Remédiation Ciblée Post-Audit (Sécurisation Anti-Overflow Hyperbolique & Alignement de Parité de $K_6$)
Created at: 2026-10-08T13:55:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur a validé expressément ("Oui") la proposition d'exécuter immédiatement les remédiations chirurgicales formulées par le collège d'audit contradictoire de l'Itération 155.
2. **Implicit Strategic Value**: L'éradication du risque d'`OverflowError` dans le solveur de Kepler-Sundman immunise le code contre tout crash numérique lors de survols hyperboliques extrêmes à grand pas. L'alignement de la parité de $K_6$ rétablit la solvabilité rigoureuse de l'équation homologique $\{H, K_6\} = H_7$ face à l'erreur BCH d'ordre 7, conférant à la bibliothèque une conformité théorique absolue avec la monographie et les théorèmes Lean 4.
3. **Alignment Criteria**:
   - Zéro adoucissement strict ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$).
   - Préservation machine des 10 invariants de Lie-Bargmann ($< 2 \times 10^{-15}$).
   - 100% des tests unitaires passants avec succès.
   - Non-régression certifiée sur les 8 benchmarks canoniques.

## Guidelines:
- Modifications chirurgicales et granulaires via `replace_file_content` (interdiction stricte de réécriture complète par overwrite sur les sources).
- Validation empirique outillée systématique avant et après chaque modification.
- Conservation de l'invariance symplectique et de la réversibilité microscopique.

## Plan:

### Phase 1: Sécurisation Anti-Overflow du Solveur de Kepler-Sundman
- **Context**: `solvers/kepler_sundman.py` crashait par `OverflowError: math range error` lorsque $\alpha \ll 0$ combiné à $dt \ge 1$ car $z_t = \alpha \tau^2$ n'était pas borné dans `_eval_sundman` et l'amorce de Simpson sur-estimait $\tau_0$.
- **Objective**: Clamper $z_t = \max(\alpha \tau^2, Z_{\mathrm{MIN\_HYPERBOLIC}})$ dans `_eval_sundman`, borner l'amorce de Simpson pour $\alpha < 0$, et étendre la projection de Battin.
- **Deliverables**:
  - `kepler_sundman.py` sécurisé contre tout dépassement d'exposant float64.
  - Test unitaire dédié dans `tests/unit/solvers/test_kepler_sundman.py` éprouvant le cas d'overflow critique ($r_0=1, v_0=100, \alpha=-10\,000, dt=100$).
- **Process**:
  1. Éditer `kepler_sundman.py` avec `replace_file_content`.
  2. Ajouter le cas de test hyperbolique extrême.
  3. Exécuter `pytest` sur `test_kepler_sundman.py`.

### Phase 2: Alignement de la Parité en Impulsion du Générateur $K_6$
- **Context**: Dans `integrators/correctors.py`, $K_6$ était implémenté comme $\sum (F_i \cdot p_i)^2$ (degré 2 pair en impulsion), ce qui rendait $\{H, K_6\}$ impair, incapable d'annuler le terme d'erreur BCH $H_7$ (pair).
- **Objective**: Réaligner $K_6$ sur le manuscrit LaTeX (Éq. 130) et le certificat Lean 4 (`02_CorrectorBM6_Theorems.lean`) sous forme impaire en impulsion :
  $$K_6 = \kappa_6 \sum_i \left(1 + \frac{\|\mathbf{p}_i\|^2}{m_i^2}\right) \left( \mathbf{F}_i \cdot \frac{\mathbf{p}_i}{m_i} \right)$$
  Dériver analytiquement le gradient variationnel $\nabla_p K_6$ et la force conjuguée $-\nabla_q K_6$.
- **Deliverables**:
  - `integrators/correctors.py` mis à jour avec le générateur impair et ses dérivées analytiques exactes.
  - Tests unitaires dans `tests/unit/integrators/test_correctors.py` validant l'inversion exacte $\mathcal{C}_6^{-1} \circ \mathcal{C}_6 = \mathrm{Id}$ et la concordance avec les différences finies.
- **Process**:
  1. Éditer `correctors.py`.
  2. Mettre à jour `test_correctors.py`.
  3. Valider l'accord numérique et l'inversion.

### Phase 3: Harmonisation Épistémique de $\mathrm{PBM}_4$ & Re-Certification Globale
- **Context**: Harmoniser la sémantique de $\mathrm{PBM}_4$ entre code et théorie (ordre 4 optimisé à -35% d'erreur) et attester de la non-régression intégrale.
- **Objective**: Ajuster `effective_order` dans `processed.py` et lancer la suite complète des tests unitaires et benchmarks.
- **Deliverables**:
  - `integrators/processed.py` harmonisé.
  - 100% de la suite de tests unitaires passante.
- **Process**:
  1. Éditer `processed.py`.
  2. Exécuter la suite de tests complète via pytest.

## Acceptance Criteria:
### Phase 1:
- Le cas $r_0=1, v_0=100, \alpha=-10\,000, dt=100$ s'exécute sans `OverflowError` et converge à $< 10^{-15}$.
- 100% des tests de `test_kepler_sundman.py` passent.
### Phase 2:
- $K_6$ est strictement de parité impaire en $\mathbf{p}$ ($K_6(q, -p) = -K_6(q, p)$).
- L'inversion $\mathcal{C}_6^{-1} \circ \mathcal{C}_6 = \mathrm{Id}$ est vérifiée à $< 10^{-14}$.
### Phase 3:
- 100% des tests unitaires passent sans régression.
