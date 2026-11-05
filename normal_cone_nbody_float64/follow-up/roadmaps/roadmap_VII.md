# Strategic Roadmap : Éradication des Défaillances Critiques Découvertes par l'Audit et Harmonisation Théorique Intégrale (Roadmap VII)

## Strategic Overview:
- **Objectif Stratégique Général** : Éradiquer définitivement les 4 divergences critiques découvertes lors de l'audit contradictoire à 9 sous-agents (Itération 149), à savoir : le piège de Halley / stall hyperbolique dans le solveur de Kepler-Sundman, l'explosion d'énergie en collision interne dans le solveur de super-nœuds non régularisé, l'omission du terme direct de Leibniz dans le gradient du correcteur $\mathcal{C}_6$, et la clarification épistémique de l'ordre de convergence de $\mathrm{PBM}_4$.
- **Règles d'Or & Contraintes Inviolables** :
  1. *Sanctuarisation Documentaire* : Le répertoire du manuscrit de recherche [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/) reste strictement en lecture seule, sauf en Phase 4 où des rectifications ciblées expressément validées pourront être apportées.
  2. *Zéro-Adoucissement Strict* : Aucun adoucissement artificiel ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$). Toute singularité doit être résolue par difféomorphisme analytique de Sundman ou régularisation canonique.
  3. *Invariance Symplectique et de Lie-Bargmann* : Préservation absolue de la forme de Darboux et des 10 invariants de Noether à la précision machine float64 IEEE 754 ($< 10^{-14}$).
  4. *Convergence et Réversibilité Monotones* : Aucun solveur ne doit échouer silencieusement sans lever d'exception ni basculer sur un solveur monotone bracketé garanti.
  5. *Certification Continue* : 100 % des tests unitaires existants et nouveaux doivent passer avec succès à chaque itération.
- **Document Maître de Référence** : [`audit_theorique_post_remediation_chapitres_3_4_5.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/audit_theorique_post_remediation_chapitres_3_4_5.md).

---

## Iteration 150 Strategy:
### Iteration 150: Phase 1 — Sécurisation Absolue du Solveur Kepler-Sundman (Amorce Asymptotique, Fallback Monotone Garanti, Suppression du Clamping z)
- **Scope & Objectives**:
  - `solvers/kepler_sundman.py` :
    1. Résoudre le piège d'Halley en rejetant les pas où $2(F')^2 - FF'' \le 0.1(F')^2$ ou $|\delta\tau| > 0.5\tau$, et substituer un pas de Newton amorti ou une bissection.
    2. Implémenter un fallback de sécurité garanti basé sur la méthode de Brent/Dekker monotone sur $F(\tau) = \Delta t$ dès que 15 itérations d'Halley échouent à converger, interdisant tout retour silencieux à $\tau = 10^{-15}$.
    3. Formuler l'amorce asymptotique hyperbolique $\tau_0 \approx \Delta t / (r_0 + \frac{1}{2} v_0 \Delta t)$ pour les régimes à haute énergie ($\alpha \ll 0, v_0 \ge 50$), éliminant le "Halley stall".
    4. Supprimer le clamping arbitraire $z = -500\,000$ en exploitant la forme fermée stable $\sinh(v)/v$ avec protection contre l'overflow exponentiel via division analytique normalisée.
    5. Substituer au seuil de collision radiale variable `rho0` un critère invariant d'échelle $\frac{L^2}{G^2 M^2 \mu^2} < \epsilon_{\mathrm{rad}}$.
  - Enrichir `tests/unit/solvers/test_kepler_sundman.py` avec 1000 orbites aléatoires à haute excentricité ($e \ge 0.999$) et plonges hyperboliques extrêmes ($v_0 \ge 500$), garantissant 100% de réversibilité microscopique ($\mathcal{S}(-\Delta t) \circ \mathcal{S}(\Delta t) = \mathrm{Id}$ à $< 10^{-14}$).
- **Dependencies & Continuity**: Sécurise le cœur analytique universel 2-corps de l'ensemble de la bibliothèque.
- **Gating Criteria**:
  - Zéro défaillance sur 1000 orbites aléatoires (réversibilité $< 10^{-14}$ sur 100% des cas).
  - Absence totale de clamping dur sur $z$ et de retour silencieux à $\tau = 10^{-15}$.
  - 100% des tests unitaires validés.

---

## Iteration 151 Strategy:
### Iteration 151: Phase 2 — Régularisation Canonique de l'Intégrateur Interne de Super-Nœuds sous Temps Déplié
- **Scope & Objectives**:
  - `solvers/supernode.py` :
    1. Remplacer le Leapfrog à pas physique variable par une formulation symplectique canonique sous le temps déplié de Poincaré-Sundman $H \mapsto \rho_S (H - E)$ sur l'espace des phases étendu, ou par une décomposition hiérarchique locale isolant les sous-binaires compactes formées au sein du super-nœud ($k \ge 3$) sous Stumpff-Sundman.
    2. Ajouter un filtre vectoriel anti-zéro `np.where(d2 > 0, G / (d2 * d), 0.0)` dans `_compute_internal_cluster_acc` pour prévenir tout `NaN` lors de chocs géométriques stricts ($d_{ij} = 0$).
    3. Ajouter une exception explicite si `t_accum < abs_dt` lorsque `steps == max_steps`.
  - Enrichir `tests/unit/solvers/test_supernode.py` avec des tests de quasi-collision binaire interne pour des paramètres d'impact $b \in [10^{-1}, 10^{-6}, 0.0]$, démontrant la stabilité de l'énergie sans explosion en $b^{-2}$.
- **Dependencies & Continuity**: Nécessite la solidité du solveur de base (Itération 150).
- **Gating Criteria**:
  - Conservation de l'énergie vérifiée lors de quasi-collisions serrées dans un super-nœud sans divergence exponentielle.
  - Zéro `NaN` et zéro division par zéro non masquée.

---

## Iteration 152 Strategy:
### Iteration 152: Phase 3 — Complétion de Leibniz sur le Gradient de $\mathcal{C}_6$ et Télescopage d'Ombre Inter-Pas
- **Scope & Objectives**:
  - `integrators/correctors.py` :
    1. Ajouter le terme direct de Leibniz $+ \frac{2}{m_j^2} (F_j \cdot v_j) F_j$ dans `CorrectorC6.apply()` pour compléter rigoureusement le gradient d'impulsion $\nabla_p K_6$.
  - `system/nbody_system.py` :
    2. Optimiser la chaîne télescopique $\mathcal{C} \circ \mathcal{S}^N \circ \mathcal{C}^{-1}$ pour éviter d'évaluer le correcteur physique à chaque macro-pas pour le simple logging, réduisant le coût computationnel de $2M$ à 2 évaluations de correction.
    3. Brancher l'opérateur de connexion de jauge `apply_gauge_connection` lors des transitions d'atlas en cours d'intégration.
  - Valider sur `tests/unit/integrators/test_correctors.py` et les benchmarks à 4 corps.
- **Dependencies & Continuity**: S'appuie sur la régularité des forces d'apex et de Kepler.
- **Gating Criteria**:
  - Exactitude de $\nabla_p K_6$ vérifiée contre différences finies à $< 10^{-11}$.
  - Préservation de la chaîne télescopique sans surcoût d'évaluation.

---

## Iteration 153 Strategy:
### Iteration 153: Phase 4 — Harmonisation Mathématique et Typographique de la Monographie
- **Scope & Objectives**:
  - `redaction/unsoftened_conormal_nbody3d/chapters/ch03_three_body_paradigm/sec02_suzuki_barrier_mcgehee/content.tex` :
    1. Rectifier la coquille typographique à la ligne 82 où le dénominateur $z$ manquait dans l'expression de la dérivée de Stumpff $c_k'(z) = \frac{c_{k-1}(z) - k c_k(z)}{2z}$.
  - `redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/` :
    2. Clarifier le statut épistémique de $\mathrm{PBM}_4$ : expliciter qu'il s'agit d'un schéma d'ordre 4 à constante d'erreur réduite (de 35%) et non d'ordre effectif 6, en raison de l'impossibilité d'annuler les 6 crochets de Hall de $\mathcal{L}_5$ avec un générateur scalaire unique.
- **Dependencies & Continuity**: Formalisation théorique des découvertes de l'audit.
- **Gating Criteria**:
  - Manuscrit LaTeX compilable sans erreur.
  - Parfaite cohérence entre le texte de recherche et le comportement numérique certifié du code.

---

## Iteration 154 Strategy:
### Iteration 154: Phase 5 — Re-Certification Complète des 8 Benchmarks, Planches 300 DPI et Appli Démonstrateur
- **Scope & Objectives**:
  - Exécuter la suite complète des 8 benchmarks canoniques (Figure-8, Burrau, Radial collision, 4-corps hiérarchique, Arenstorf, Sitnikov, Jovien, Collision quadruple) sur l'ensemble des 6 intégrateurs.
  - Régénérer les 4 planches de publication haute résolution 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`, `dashboard_multi_integrator_invariants_300dpi.png`).
  - Mettre à jour `showcase_data.json` et resynchroniser `interactive_conormal_showcase.html`.
- **Dependencies & Continuity**: Aboutissement de l'ensemble des remédiations des Itérations 150 à 153.
- **Gating Criteria**:
  - 100% des benchmarks certifiés [PASS] sur tous les intégrateurs.
  - 100% des tests unitaires validés (aucun test ignoré ou en échec).
  - Planches 300 DPI et application interactive parfaitement à jour.

---

## Multi-Iteration Checklist:
- [x] Iteration 150: Phase 1 — Sécurisation Absolue du Solveur Kepler-Sundman (Amorce Asymptotique, Fallback Monotone Garanti, Suppression du Clamping z)
- [x] Iteration 151: Phase 2 — Régularisation Canonique de l'Intégrateur Interne de Super-Nœuds sous Temps Déplié
- [x] Iteration 152: Phase 3 — Complétion de Leibniz sur le Gradient de $\mathcal{C}_6$ et Télescopage d'Ombre Inter-Pas
- [x] Iteration 153: Phase 4 — Harmonisation Mathématique et Typographique de la Monographie
- [x] Iteration 154: Phase 5 — Re-Certification Complète des 8 Benchmarks, Planches 300 DPI et Appli Démonstrateur
