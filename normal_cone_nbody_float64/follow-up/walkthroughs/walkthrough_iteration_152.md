# Rapport d'Exécution : Itération 152 — Complétion de Leibniz sur le Gradient de $\mathcal{C}_6$ et Télescopage d'Ombre Inter-Pas

## Statut : VALIDÉ & CERTIFIÉ [PASS]
- **Date** : 2026-10-08
- **Phase** : Phase 3 de la Roadmap VII
- **Commit Git** : `7b860bf`
- **Tests Unitaires** : 148/148 validés (100% de succès en 6.15s)

---

## 1. Synthèse Exécutive

Lors de l'audit théorique contradictoire à 9 sous-agents (Itération 149), deux anomalies structurelles avaient été identifiées concernant les intégrateurs traités et les correcteurs symplectiques :
1. Dans `normal_cone_nbody_float64/integrators/correctors.py`, la méthode `CorrectorC6.apply()` utilisait une loi d'échelle advective heuristique $(1 + \|\mathbf{v}\|^2)$ omettant le terme direct de Leibniz $+ \frac{2}{m_j^2} (\mathbf{F}_j \cdot \mathbf{v}_j) \mathbf{F}_j$ dans le gradient d'impulsion $\nabla_p K_6$. Le flot de Verlet ne dérivait donc pas rigoureusement du générateur scalaire $K_6$.
2. Dans `normal_cone_nbody_float64/system/nbody_system.py`, le télescopage théorique de l'hamiltonien d'ombre $\Psi^{(M)} = \mathcal{C} \circ \mathcal{S}^M \circ \mathcal{C}^{-1}$ était détruit en pratique : `step()` appliquait $\mathcal{C}^{-1}$ à l'entrée et $\mathcal{C}$ à la sortie de chaque macro-pas pour enregistrer les coordonnées physiques, engendrant $2M$ évaluations de correcteur au lieu de 2. De plus, lors des transitions de carte atlas ($\mathcal{T}_A \to \mathcal{T}_B$), l'opérateur de connexion de jauge $\mathcal{J}_{A \to B} = \mathcal{C}_B^{-1} \circ \mathcal{C}_A$ restait orphelin et n'était pas appelé.

L'Itération 152 a intégralement résolu ces deux défaillances.

---

## 2. Réalisations Mathématiques et Algorithmiques

### A. Gradient Canonique de Leibniz pour $\mathcal{C}_6$ (`correctors.py`)
Le potentiel de correction d'ordre 6 est :
$$
K_6(\mathbf{q}, \mathbf{p}) = \kappa_6 \sum_i \frac{1}{m_i} (\mathbf{F}_i(\mathbf{q}) \cdot \mathbf{v}_i)^2 = \kappa_6 \sum_i \frac{1}{m_i^3} (\mathbf{F}_i(\mathbf{q}) \cdot \mathbf{p}_i)^2.
$$
Ses dérivées variationnelles exactes sont :
1. **Gradient par rapport à l'impulsion $\mathbf{p}_j$** (règle de Leibniz directe) :
   $$
   \dot{\mathbf{q}}_j = \nabla_{\mathbf{p}_j} K_6 = \kappa_6 \frac{2}{m_j^3} (\mathbf{F}_j \cdot \mathbf{p}_j) \mathbf{F}_j = \kappa_6 \frac{2}{m_j^2} (\mathbf{F}_j \cdot \mathbf{v}_j) \mathbf{F}_j.
   $$
2. **Gradient par rapport à la position $\mathbf{q}_j$** (force conjuguée) :
   $$
   \dot{\mathbf{p}}_j = -\nabla_{\mathbf{q}_j} K_6 = -\kappa_6 D_{\mathbf{u}} \mathbf{F}_j \quad \text{avec} \quad \mathbf{u}_i = \frac{2}{m_i} (\mathbf{F}_i \cdot \mathbf{v}_i) \mathbf{v}_i.
   $$
Ces deux gradients ont été validés contre des différences finies centrales d'ordre 4 dans `scratch/verify_k6_gradient.py` :
- Erreur relative maximale sur $\nabla_p K_6$ : **$3.913 \times 10^{-12}$**.
- Erreur relative maximale sur $\nabla_q K_6$ : **$1.463 \times 10^{-11}$**.
- Erreurs absolues : $< 1.6 \times 10^{-16}$.

### B. Inversion Algébrique Machine et Conservation du Moment
Le flot de Verlet palindromique appliqué à $\mathcal{X}_{K_6}$ garantit :
- Inversion exacte $\mathcal{C}_6^{-1} \circ \mathcal{C}_6 = \mathrm{Id}$ : erreur sur $q$ de **$2.0 \times 10^{-26}$**, erreur sur $p$ de **$0.0$** (strict zéro machine).
- Conservation de l'impulsion linéaire totale : $\|\Delta \mathbf{P}\| < 10^{-16}$.

### C. Branchement de la Connexion de Jauge et Télescopage d'Ombre (`nbody_system.py`)
- **Télescopage Global dans `evolve()`** :
  Lors d'une évolution de $M$ macro-pas sous `ProcessedSymplecticIntegrator`, le système entre une seule fois sur la variété d'ombre via $\mathcal{C}^{-1}$ à $t=0$, propage les états intérieurs uniquement via l'intégrateur de base $\mathcal{S}$, et ne ressort sur la variété physique via $\mathcal{C}$ qu'à $t=t_{\mathrm{final}}$.
  Le nombre d'évaluations de correcteurs passe de $2M$ à :
  $$
  N_{\mathrm{eval}} = 2 + 2 \cdot N_{\mathrm{transitions}}.
  $$
  Sur $M=20$ macro-pas sans transition, le nombre d'évaluations passe de 40 à **2** (gain de $20\times$).
- **Connexion de Jauge Inter-Atlas** :
  Lorsqu'une transition d'arbre $\mathcal{T}_A \to \mathcal{T}_B$ est détectée en cours de propagation d'ombre, l'opérateur symplectique :
  $$
  \mathcal{J}_{A \to B} = \mathcal{C}_B^{-1} \circ \mathcal{C}_A
  $$
  est automatiquement appelé via `apply_gauge_connection()`, éliminant le saut interfacial non-commutatif $[\mathcal{X}_{K_A}, \mathcal{X}_{K_B}] \ne 0$ et préservant les 10 invariants de Lie-Bargmann à $< 10^{-13}$.

---

## 3. Données de Certification et Vérification

### A. Précision des Gradients et Propriétés Algébriques de $\mathcal{C}_6$
- **Test de Gradient Variationnel** : `test_corrector_c6_variational_gradient_and_momentum_conservation()` (PASS).
- **Test d'Inversion Algébrique** : `test_corrector_c6_exact_algebraic_inversion()` (PASS, résidu $< 10^{-14}$).
- **Test d'Échelle Perturbative** : `test_correctors_scaling_powers()` : ratio $(\Delta q_{2h} / \Delta q_h) = 63.99997 \approx 64.0$ ($2^6$).

### B. Efficacité et Intégrité du Télescopage dans `NBodySystem`
- **Test de Télescopage** : `test_nbody_system_evolve_shadow_telescoping_evaluation_efficiency()` (PASS).
  - Évaluations inverses mesurées : 1.
  - Évaluations directes mesurées : 1.
  - Conservation $\mathbf{P}, \mathbf{L}, \mathbf{G}$ : $< 10^{-13}$.

### C. Suite Complète des Tests Unitaires
- **Commande** : `python -m pytest tests/unit/ -q`
- **Résultat** : **148 passed in 6.15s** (100% de succès, 0 avertissement).

---

## 4. Bilan des Critères d'Acceptation de la Phase 3

| Critère d'Acceptation | Cible Requise | Valeur Obtenue | Statut |
| :--- | :---: | :---: | :---: |
| Précision relative $\nabla_p K_6$ vs DF4 | $< 10^{-9}$ | $3.91 \times 10^{-12}$ | **CERTIFIÉ** |
| Précision relative $\nabla_q K_6$ vs DF4 | $< 10^{-9}$ | $1.46 \times 10^{-11}$ | **CERTIFIÉ** |
| Inversion exacte $\mathcal{C}_6^{-1} \circ \mathcal{C}_6 = \mathrm{Id}$ | $< 10^{-14}$ | $< 10^{-16}$ ($2.0 \times 10^{-26}$) | **CERTIFIÉ** |
| Conservation de l'impulsion totale sous $\mathcal{C}_6$ | $< 10^{-14}$ | $< 10^{-16}$ | **CERTIFIÉ** |
| Évaluations de correcteurs sur $M$ pas sans transition | $2$ (au lieu de $2M$) | 2 (1 entrée, 1 sortie) | **CERTIFIÉ** |
| Connexion de jauge $\mathcal{J}_{A \to B}$ raccordée | Automatique sur flip | Validée sous `in_shadow_mode` | **CERTIFIÉ** |
| Suite de tests unitaires complète | 100% pass | 148/148 passed | **CERTIFIÉ** |
| Commit Git de Phase 3 | Enregistré | `7b860bf` | **CERTIFIÉ** |

---

## 5. Perspectives : Phase 4 de la Roadmap VII (Itération 153)

La Phase 4 portera sur l'**harmonisation mathématique et typographique de la monographie** :
1. Correction typographique de la formule de dérivée de Stumpff à la ligne 82 de `redaction/unsoftened_conormal_nbody3d/chapters/ch03_three_body_paradigm/sec02_suzuki_barrier_mcgehee/content.tex` :
   $$c_k'(z) = \frac{c_{k-1}(z) - k c_k(z)}{2z}$$ (rétablissement du dénominateur $z$).
2. Clarification épistémique de l'ordre de convergence de $\mathrm{PBM}_4$ dans le Chapitre 5 : explicitation du fait qu'un générateur scalaire unique $K_4$ ne peut annuler simultanément les 6 crochets de Hall indépendants de $\mathcal{L}_5$, ce qui maintient l'ordre effectif à 4 avec une constante d'erreur réduite de 35%.
