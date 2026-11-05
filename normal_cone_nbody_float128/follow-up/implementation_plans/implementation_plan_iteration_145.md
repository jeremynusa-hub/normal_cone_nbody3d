# Implementation Plan of Iteration 145 : Phase 2 — Remédiation des Super-Nœuds et Régularisation du Potentiel d'Apex du Chapitre 4
Created at: 2026-10-08T10:16:30+02:00

## Context:
1. **Explicit Request**: L'utilisateur a ordonné : « Continue la roadmap donc ». Dans le cadre de la Roadmap VI, cela active immédiatement la **Phase 2 (Itération 145)** consacrée à la remédiation théorique et numérique des composants du Chapitre 4 du manuscrit (Super-nœuds, Écrantage Hermite $C^8$, et Horloge dépliée de Sundman).
2. **Implicit Strategic Value**: Assurer que la décomposition conormale de Fulton--MacPherson ne souffre d'aucune discontinuité numérique : éliminer le broutement d'activation des super-nœuds (flickering), garantir un gradient d'écrantage $C^8$ strictement monotone et insensible au bruit machine à $s=0.5$ ($\kappa \approx 1$), et implémenter l'intégration interne des super-nœuds sous temps déplié de Sundman sans adoucisseur résiduel ($d_{\mathrm{safe}} \equiv 0$).
3. **Alignment Criteria**: L'état du super-nœud reste stable sous oscillations de marée autour de $\eta_{\mathrm{crit}}$ grâce à l'hystérésis $[0.60, 0.75]$ ; la spline $C^8$ évaluée en forme centrée $\xi = s - 0.5$ s'annule aux dérivées $1 \dots 8$ sans saut numérique ; le solveur de super-nœuds conserve l'impulsion et l'énergie sans adoucissement masqué ; 100 % des tests unitaires sont validés.

---

## Guidelines:
- **Method**: Respect strict de la déformation vers le cône normal (Chapitre 4). Formulation polynomiale centrée en $\xi = s - 0.5$ pour les splines d'ordre élevé afin de préserver un conditionnement $\kappa \approx 1$.
- **Red Lines**: Interdiction absolue de tout adoucisseur résiduel ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$). Interdiction de tout recentrage artificiel non-canonique de vitesse ou de quantité de mouvement. Le répertoire `redaction/unsoftened_conormal_nbody3d/` demeure sanctuarisé en lecture seule.
- **Architecture**: Maintien de l'indépendance modulaire entre `potential/supernode_detector.py` (détection), `potential/screening.py` (écrantage conormal), `solvers/supernode.py` (flot interne régularisé), et `system/nbody_system.py` (orchestration de l'état).

---

## Plan:

### Step 1: Bascule de Schmitt à Hystérésis dans `potential/supernode_detector.py`
- **Context**: Actuellement, `detect_supernodes` évalue le seuil critique brut $\eta_S \le \eta_{\mathrm{crit}} \approx 0.693361$ sans mémoire d'état, causant des basculements intempestifs (chattering) à chaque pas lors des passages frôlants.
- **Objective**: Mettre en place une hystérésis de Schmitt à deux seuils $[\eta_{\mathrm{enter}} = 0.60, \eta_{\mathrm{exit}} = 0.75]$ conditionnée par l'historique des amas actifs.
- **Deliverables**:
  - `potential/supernode_detector.py`: Fonction `detect_supernodes` enrichie avec le paramètre `active_clusters: Optional[List[Tuple[int, ...]]] = None` et constantes canoniques `ETA_ENTER_HILL = 0.60`, `ETA_EXIT_HILL = 0.75`.
  - `system/nbody_system.py`: Mémorisation d'état `self.active_supernodes` transmise à chaque pas macro pour stabiliser la topologie des super-nœuds.
- **Process**:
  1. Définir les constantes `ETA_ENTER_HILL = 0.60` et `ETA_EXIT_HILL = 0.75` conformes au Chapitre 4.
  2. Ajuster le critère de détection : un amas candidat actif se maintient si $\eta_S \le \eta_{\mathrm{exit}}$, tandis qu'un nouvel amas nécessite $\eta_S \le \eta_{\mathrm{enter}}$.
  3. Mettre à jour `NBodySystem` pour persister `self.active_supernodes` entre les pas.

### Step 2: Formulation Centrée Conditionnée ($\kappa \approx 1$) de la Spline $C^8$ dans `potential/screening.py`
- **Context**: L'évaluation monomiale non centrée de la spline Hermite de degré 17 pour $C^8$ induit un conditionnement $\kappa \approx 1902$ et un bruit de $-6.6 \times 10^{-14}$ à $s=0.5$ avec possible inversion de signe de $d\Theta/du$.
- **Objective**: Réécrire l'évaluation de $\sigma_8(s)$ et $\sigma_8'(s)$ sous forme strictement centrée en $\xi = s - 0.5 \in [-0.5, 0.5]$.
- **Deliverables**:
  - `potential/screening.py`: Fonction `eval_screening_theta` et version vectorisée `eval_screening_theta_vectorized` avec polynôme centré en puissances paires de $\xi^2 = (s - 0.5)^2$.
- **Process**:
  1. Exprimer $\sigma_8'(\xi) = 218790 (0.25 - \xi^2)^8$ et son intégrale exacte $\sigma_8(\xi) = 0.5 + \xi \sum_{k=0}^8 c_k \xi^{2k}$.
  2. Dériver les coefficients exacts $c_k = \frac{218790 \binom{8}{k} (-1)^k (0.25)^{8-k}}{2k+1}$.
  3. Remplacer les formes monomiales instables dans `eval_screening_theta` et `eval_screening_theta_vectorized`.
  4. Vérifier la stricte monotonie, la parité impaire par rapport à $\xi$, et l'absence de saut à $\xi = 0$.

### Step 3: Solveur de Super-Nœuds sous Horloge Dépliée de Sundman sans Adoucisseur dans `solvers/supernode.py`
- **Context**: `solvers/supernode.py` utilise un Leapfrog cartésien à pas physique fixe avec un garde-fou $d_{\mathrm{safe}} = 10^{-14}$ (adoucissement masqué) et un recentrage non canonique de vitesse $v_{\mathrm{cur}} -= P_{\mathrm{rel}} / M_{\mathrm{tot}}$.
- **Objective**: Régulariser l'évolution interne par difféomorphisme temporel de Sundman $\mathrm{d}t_S = \rho_S\,\mathrm{d}\tau_S$, éradiquer $d_{\mathrm{safe}}$, et éliminer la dissipation artificielle de quantité de mouvement.
- **Deliverables**:
  - `solvers/supernode.py`: Implémentation conforme au Théorème 4.2 sous pas de temps Sundmanisé $\Delta \tau_S = \Delta t / \rho_S$, avec conservation exacte et analytique de l'impulsion linéaire ($\sum F_i \equiv 0$) sans soustraction artificielle, et zéro adoucisseur ($d_{\mathrm{safe}} = 0$).
- **Process**:
  1. Remplacer le pas cartésien fixe par une intégration sous variable temporelle locale régularisée $\tau_S$ gouvernée par l'hyper-rayon $\rho_S$.
  2. Supprimer la clause $d_{\mathrm{safe}} = \max(d, 10^{-14})$ pour opérer en gravité newtonienne pure non adoucie.
  3. Supprimer le recentrage artificiel `v_cur -= (P_rel_cdm / m_tot)` et vérifier que l'accélération interne newtonienne préserve naturellement $\sum m_i a_i = 0$ à la précision machine.

### Step 4: Tests Unitaires Miroirs et Validation de Non-Régression
- **Context**: L'intégrité de la suite de tests unitaires (129 tests) doit être strictement préservée, avec des ajouts spécifiques testant l'hystérésis, la stabilité de la spline $C^8$ et le solveur de super-nœuds.
- **Objective**: Écrire des tests unitaires ciblés pour chaque composant remédié et attester du 100 % de succès.
- **Deliverables**:
  - `tests/unit/potential/test_supernode_detector.py`: Test de l'hystérésis de Schmitt (non-chattering sur oscillations $\eta \in [0.55, 0.72]$).
  - `tests/unit/potential/test_screening.py`: Test de continuité, dérivées nulles d'ordre 1 à 8, et conditionnement $\kappa \approx 1$ à $s=0.5$.
  - `tests/unit/solvers/test_supernode.py`: Test de conservation de l'impulsion et de l'énergie sans $d_{\mathrm{safe}}$ sous approche serrée.
- **Process**:
  1. Écrire les nouveaux cas de test dans les suites unitaires respectives.
  2. Exécuter `pytest tests/unit/` et corriger tout écart éventuel.

### Step 5: Validation Git Commit
- **Context**: Clôture de la Phase 2 de la Roadmap VI.
- **Objective**: Tracer de manière atomique et exhaustive les modifications dans le journal Git.
- **Deliverables**:
  - Commit Git structuré : `feat(potential,solvers): remediate Chapter 4 supernodes and screening - Phase 2 of Roadmap VI (Iteration 145)`.
- **Process**:
  1. Vérifier l'état git via `git status`.
  2. Indexer les fichiers modifiés et valider le commit.

---

## Acceptance Criteria:

### Step 1:
- `detect_supernodes` filtre les oscillations périodiques de marée sans changer d'état tant que $\eta_S$ reste dans la bande $[0.60, 0.75]$.
- `NBodySystem` conserve ses super-nœuds actifs de manière cohérente à travers les pas de temps.

### Step 2:
- La fonction de screening $C^8$ présente un conditionnement $\kappa \approx 1$ à $s=0.5$ sans saut détectable ($|\Theta(0.5^+) - \Theta(0.5^-)| < 10^{-16}$) et avec $d\Theta/du \ge 0$ strictement monotone sur $[0, 1]$.
- Les dérivées d'ordre 1 à 8 s'annulent identiquement à $u=1.0$ et $u=2.0$ à la précision machine.

### Step 3:
- Zéro adoucissement résiduel ($d_{\mathrm{safe}} \equiv 0$) dans `solvers/supernode.py`.
- L'impulsion interne $\mathbf{P}_S$ est conservée à $< 10^{-14}$ par construction symplectique sans soustraction artificielle.
- Pas d'explosion d'énergie sur les approches serrées multi-corps.

### Step 4:
- 100 % de succès sur l'ensemble de la suite de tests unitaires `tests/unit/` (>= 129 tests).

### Step 5:
- Commit Git tracé et documenté dans le repository.

---

## Observations:
> [!IMPORTANT]
> La bascule de Schmitt à hystérésis $[\eta_{\mathrm{enter}}, \eta_{\mathrm{exit}}] = [0.60, 0.75]$ garantit la conformité avec le Théorème 6.1 et empêche les chocs de splitting artificiels observés dans les régimes chaotiques.

> [!NOTE]
> La forme centrée en $\xi = s - 0.5$ pour la spline $C^8$ exploite la parité impaire de $\sigma_8(\xi) - 0.5$, réduisant de 17 à 9 le nombre de monômes à évaluer et éliminant tout risque de dépassement ou de soustraction catastrophique.
