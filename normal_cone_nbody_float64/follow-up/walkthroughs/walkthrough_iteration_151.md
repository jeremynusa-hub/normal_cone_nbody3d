# Rapport d'Exécution : Itération 151 — Régularisation Canonique de l'Intégrateur Interne de Super-Nœuds sous Temps Déplié

## Statut : VALIDÉ & CERTIFIÉ [PASS]
- **Date** : 2026-10-08
- **Phase** : Phase 2 de la Roadmap VII
- **Commit Git** : `f884c54`
- **Tests Unitaires** : 146/146 validés (100% de succès en 6.31s)

---

## 1. Synthèse Exécutive

Lors de l'audit théorique contradictoire à 9 sous-agents (Itération 149), une divergence critique avait été mise au jour dans `solvers/supernode.py` :
1. Dans un amas synchrone à $k \ge 3$ corps, lorsqu'une sous-paire $(i^*, j^*)$ subit une quasi-collision binaire interne ($r_{ij} \to 0$), l'hyper-rayon global $\rho_S$ reste borné inférieurement par la présence des spectateurs tiers. L'horloge globale $dt_S = \rho_S d\tau_S$ ne se dilate donc pas suffisamment pour désingulariser le choc binaire local.
2. L'ancien Leapfrog adaptatif en temps physique subissait une explosion de l'erreur d'énergie inversement proportionnelle au carré du paramètre d'impact :
   $$
   \frac{\Delta E}{E} \propto b^{-2}
   $$
   atteignant $6.68 \times 10^3$ pour $b=10^{-3}$, $6.60 \times 10^{10}$ pour $b=10^{-6}$, et $1.92 \times 10^{11}$ pour $b=0.0$.
3. La boucle temporelle s'arrêtait prématurément par épuisement de `max_steps = 100_000` à $t \approx 0.099$ (sur un pas macro de $0.2$), retournant silencieusement un état non convergé avec des corps à $v \approx 250$ m/s au péricentre.
4. En cas de collision géométrique stricte ($d_{ij} = 0$), `inv_d3 = G / (d**3)` provoquait une division par zéro et injectait des `NaN`.

L'Itération 151 a totalement résolu ces défaillances en intégrant une isolation analytique de sous-paire compacte par l'opérateur universel de Stumpff-Sundman au sein du super-nœud.

---

## 2. Réalisations Mathématiques et Algorithmiques

### A. Filtrage Anti-Division par Zéro et Exception d'Incomplétude
- **Protection Vectorielle** : Dans `_compute_internal_cluster_acc`, le calcul de la force de gravité interne a été immunisé contre $d_{ij} = 0$ :
  ```python
  inv_d3 = np.zeros_like(d2)
  valid = d2 > 0.0
  inv_d3[valid] = G / (d2[valid] * d[valid])
  ```
  garantissant l'absence absolue de `NaN` ou de `RuntimeWarning`.
- **Exception d'Incomplétude Temporelle** : Si la boucle de sous-pas atteint `max_steps` sans avoir couvert l'intervalle temporel (`t_accum < abs_dt`), une `RuntimeError` circonstanciée est levée, interdisant tout retour silencieux d'état tronqué.

### B. Régularisation Canonique par Découplage de Stumpff-Sundman
- **Critère d'Isolation de Sous-Paire** : Au sein du super-nœud, lorsque la distance minimale entre deux corps vérifie :
  $$
  r_{\min} < 0.5 \cdot r_{\mathrm{second}} \quad \text{et} \quad r_{\min} < 0.2 \cdot \rho_S
  $$
  la paire la plus serrée $(i^*, j^*)$ est isolée dynamiquement.
- **Propagation Symplectique Hiérarchisée (Strang Splitting)** :
  1. Demi-kick de marée gravitationnelle appliqué par les spectateurs tiers sur $i^*$ et $j^*$.
  2. Résolution analytique exacte du mouvement képlérien relatif $\mathbf{r}^* = \mathbf{q}_{j^*} - \mathbf{q}_{i^*}$ et $\mathbf{w}^*$ par `step_kepler_sundman` sous le difféomorphisme universel $dt = r^* ds$.
  3. Dérive barycentrique de la paire et dérive des spectateurs tiers.
  4. Second demi-kick de marée gravitationnelle.
- **Formule Unimodulaire de Battin** : Dans `kepler_sundman.py`, la relation unimodulaire exacte de Battin :
  $$
  \dot{g} = \frac{1 + \dot{f} g}{f}
  $$
  a été imposée, garantissant la conservation du Wronskien et du moment cinétique à la précision machine float64 ($2.04 \times 10^{-14}$).

---

## 3. Données de Certification et Vérification

### A. Banc d'Essai de Quasi-Collision Binaire Interne ($b \in [0.5, 0.0]$)
Le script de diagnostic `scratch/diagnose_supernode_collision.py` a été exécuté sur un amas de 3 corps ($m = [1.0, 1.0, 0.5]$) avec plonge vers le péricentre :

| Paramètre d'Impact $b$ | Ancienne Erreur Rel. $\Delta E/E$ | Nouvelle Erreur Rel. $\Delta E/E$ | Nombre de Sous-Pas | Statut |
| :---: | :---: | :---: | :---: | :---: |
| $0.5$ (large) | $6.92 \times 10^{-4}$ | $6.92 \times 10^{-4}$ | 5 | **CONFORME** |
| $0.1$ | $1.09 \times 10^{-2}$ | $1.09 \times 10^{-2}$ | 11 | **CONFORME** |
| $10^{-2}$ | $7.08 \times 10^{-1}$ | $5.99 \times 10^{-3}$ | 33 | **RÉGULARISÉ** |
| $10^{-3}$ | $6.68 \times 10^{3}$ (échec) | $6.75 \times 10^{-3}$ | 45 | **RÉGULARISÉ** |
| $10^{-4}$ | $6.60 \times 10^{6}$ (échec) | $6.87 \times 10^{-3}$ | 49 | **RÉGULARISÉ** |
| $10^{-6}$ | $6.60 \times 10^{10}$ (échec) | $6.91 \times 10^{-3}$ | 52 | **RÉGULARISÉ** |
| $0.0$ (collision frontale) | $1.92 \times 10^{11}$ (échec) | $6.91 \times 10^{-3}$ | 52 | **RÉGULARISÉ** |

- **Invariance de Lie-Bargmann** :
  - Impulsion linéaire : $|\Delta \mathbf{P}| < 10^{-16}$.
  - Moment cinétique : $|\Delta \mathbf{L}| < 2.04 \times 10^{-14}$.
- **Éradication du Gel Zénon** : Le nombre de sous-pas reste strictement borné entre 5 et 52, éliminant totalement l'épuisement des pas.

### B. Suite Complète des Tests Unitaires
- **Commande** : `python -m pytest tests/unit/ -q`
- **Résultat** : **146 passed in 6.31s** (100% de succès, 3 nouveaux tests ajoutés).

---

## 4. Bilan des Critères d'Acceptation de la Phase 2

| Critère d'Acceptation | Cible Requise | Valeur Obtenue | Statut |
| :--- | :---: | :---: | :---: |
| Erreur relative d'énergie à $b \le 10^{-4}$ | Bornée ($< 1\%$) | $6.87 \times 10^{-3}$ ($0.69\%$) | **CERTIFIÉ** |
| Erreur relative d'énergie en choc frontal ($b = 0.0$) | Bornée ($< 1\%$) | $6.91 \times 10^{-3}$ ($0.69\%$) | **CERTIFIÉ** |
| Nombre de sous-pas requis pour $b \to 0$ | Borné ($N \ll 100\,000$) | 52 sous-pas | **CERTIFIÉ** |
| Protection anti-division par zéro ($d_{ij} = 0$) | 0 `NaN`, 0 `RuntimeWarning` | 0 `NaN`, 0 avertissement | **CERTIFIÉ** |
| Exception explicite si `t_accum < abs_dt` | Levée de `RuntimeError` | Validée par test unitaire | **CERTIFIÉ** |
| Conservation du moment cinétique SO(3) | $|\Delta \mathbf{L}| < 10^{-13}$ | $2.04 \times 10^{-14}$ | **CERTIFIÉ** |
| Suite de tests unitaires complète | 100% pass | 146/146 passed | **CERTIFIÉ** |
| Commit Git de Phase 2 | Enregistré | `f884c54` | **CERTIFIÉ** |

---

## 5. Perspectives : Phase 3 de la Roadmap VII (Itération 152)

La Phase 3 s'attaque à la **complétion de Leibniz sur le gradient de $\mathcal{C}_6$ et au télescopage d'ombre inter-pas** :
1. Ajout dans `integrators/correctors.py` (`CorrectorC6.apply()`) du terme direct de la règle de Leibniz $+ \frac{2}{m_j^2} (\mathbf{F}_j \cdot \mathbf{v}_j) \mathbf{F}_j$ dans $\nabla_p K_6$.
2. Optimisation de la chaîne télescopique $\mathcal{C} \circ \mathcal{S}^N \circ \mathcal{C}^{-1}$ dans `system/nbody_system.py` pour éviter d'évaluer le correcteur à chaque macro-pas uniquement pour la télémétrie.
3. Branchement de l'opérateur de connexion de jauge `apply_gauge_connection` lors des transitions de cartes atlas.
4. Validation contre différences finies centrales d'ordre 4 ($< 10^{-11}$).
