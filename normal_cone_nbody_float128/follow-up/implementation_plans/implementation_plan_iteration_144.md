# Implementation Plan of Iteration 144 : Phase 1 — Remédiation et Certification des Solveurs Universels du Chapitre 3 (Stumpff, Sundman-Halley, McGehee)
Created at: 2026-10-08T09:48:50+02:00

## Context:
1. **Explicit Request**:
   - L'utilisateur a ordonné : "Lance une nouvelle roadmap pour corriger et tester tout ça", faisant suite aux conclusions de l'audit contradictoire à 9 agents consigné dans [`audit_theorique_chapitres_3_4_5.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/audit_theorique_chapitres_3_4_5.md).
2. **Implicit Strategic Value**:
   - Assurer que les briques analytiques et universelles de base (solveurs de Stumpff, Sundman et McGehee) fonctionnent sans clamping arbitraire, sans perte de chiffres significatifs sur les dérivées et sans distorsion de volume d'espace des phases, garantissant que les intégrateurs symplectiques d'ordre élevé (notamment Blanes-Moan BM6 aux pas rétrogrades) s'appuient sur un sous-espace képlérien exact à la précision machine float64.
3. **Alignment Criteria**:
   - Clamping artificiel $z < -100$ éradiqué dans `solvers/kepler_sundman.py`, remplacé par un garde-fou conforme à la capacité réelle float64 ($z \approx -500\,000$).
   - Annulation soustractive des dérivées $c_k'(z)$ dans `solvers/stumpff.py` résolue par expansion de Taylor/Horner jusqu'à $|z| < 0.2$, ramenant l'erreur relative sous $10^{-14}$.
   - Non-canonical velocity flip de `solvers/mcgehee.py` sécurisé pour préserver la structure canonique sans saut de phase discontinu.
   - 100 % de succès sur les tests unitaires existants et enrichis.

## Guidelines:
- **Method**:
  - Analyse mathématique des domaines de convergence et de précision machine float64.
  - Remplacement direct des formules défaillantes par des développements analytiques stables.
  - Vérification continue du déterminant unimodulaire du Wronskien $|f \dot{g} - \dot{f} g - 1| < 10^{-15}$.
- **Red Lines**:
  - Zéro modification du manuscrit LaTeX [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/).
  - Interdiction formelle d'introduire tout paramètre ad hoc d'adoucissement ($\varepsilon = 0, d_0 = 0$).
  - Ne jamais utiliser `write_to_file` avec `overwrite: true` sur des fichiers source existants ; utiliser `replace_file_content` de façon chirurgicale.
- **Architecture**:
  - `normal_cone_nbody_float64/solvers/kepler_sundman.py`
  - `normal_cone_nbody_float64/solvers/stumpff.py`
  - `normal_cone_nbody_float64/solvers/mcgehee.py`
  - `tests/unit/solvers/test_kepler_sundman.py`, `tests/unit/solvers/test_stumpff.py`, `tests/unit/solvers/test_mcgehee.py`

## Plan:
### Step 1: Élimination du Clamping Destructeur dans `solvers/kepler_sundman.py`
- **Context**: L'audit numérique (Agent 3.2, `9443d7c2`) a révélé que la ligne 90 de `kepler_sundman.py` impose `z = max(z, -100.0)`. Lors d'étapes rétrogrades de BM6 ($\Delta t < 0$) à grande vitesse, $z = \alpha s^2$ peut atteindre $-500\,000$. En clampant à $-100$, Halley évalue une équation universelle fausse, menant à une divergence immédiate ($\Delta r = 52.6$).
- **Objective**: Déplacer le seuil de garde à la limite réelle de float64 ($\cosh\sqrt{-z} \le 10^{308} \implies -z \le 710^2 \approx 500\,000$), soit $z \ge -500\,000$, et s'assurer que le dénominateur de Halley ne subit pas d'overflow.
- **Deliverables**:
  - `KeplerSundmanSolverRemediated`: `normal_cone_nbody_float64/solvers/kepler_sundman.py` mis à jour sans clamping destructeur.
- **Process**: Édition chirurgicale via `replace_file_content`.

### Step 2: Résolution de l'Annulation Soustractive dans `solvers/stumpff.py`
- **Context**: L'audit (Agent 3.2, `9443d7c2`) a identifié que le calcul de $c_k'(z) = (c_{k-1}(z) - k c_k(z))/(2z)$ pour $z \in [10^{-4}, 0.1]$ perd 4 chiffres significatifs en raison d'une annulation entre $c_{k-1}$ et $k c_k$.
- **Objective**: Étendre le seuil de série de Taylor/Horner pour les dérivées $c_k'(z)$ de $|z| < 10^{-4}$ à $|z| < 0.2$ avec 8 termes, garantissant une précision relative $\le 10^{-14}$ sur l'ensemble du domaine de transition.
- **Deliverables**:
  - `StumpffDerivativesRemediated`: `normal_cone_nbody_float64/solvers/stumpff.py` corrigé.
- **Process**: Édition ciblée via `replace_file_content`.

### Step 3: Sécurisation Canonique du Rebond de McGehee dans `solvers/mcgehee.py`
- **Context**: L'audit (Agent 3.3, `4a314cc6`) a prouvé algébriquement que $\det \mathrm{d}\mathcal{B}_C = -1$. De plus, la clause $r < 10^{-15}$ applique une inversion non canonique $w \mapsto -w$ qui flippe le moment cinétique $L \mapsto -L$.
- **Objective**: Harmoniser le traitement limite dans `mcgehee.py` pour préserver la symétrie radiale et documenter formellement dans le code le statut théorique du rebond de McGehee (reflétant le Théorème 3.4 du Chapitre 3 qui établit la nécessité de la déformation vers le cône normal).
- **Deliverables**:
  - `McGeheeSolverSecured`: `normal_cone_nbody_float64/solvers/mcgehee.py` sécurisé et commenté.
- **Process**: Édition chirurgicale via `replace_file_content`.

### Step 4: Tests Unitaires Miroirs et Validation de Non-Régression
- **Context**: Les modifications apportées aux solveurs de base doivent être validées unitairement sur des cas limites (étapes rétrogrades hyperboliques extrêmes, transition $z \in [-0.2, 0.2]$, précision machine du Wronskien).
- **Objective**: Enrichir `tests/unit/solvers/test_kepler_sundman.py`, `tests/unit/solvers/test_stumpff.py` et `tests/unit/solvers/test_mcgehee.py`, et exécuter la suite de tests complète (132 tests).
- **Deliverables**:
  - `EnhancedSolversTests`: Tests enrichis dans `tests/unit/solvers/`.
  - `TestSuiteValidation`: 100 % de succès sur l'ensemble de la suite.
- **Process**: Enrichissement des tests et exécution via `execution-engine` / pytest.

### Step 5: Validation Git Commit
- **Context**: La fin de l'itération requiert la traçabilité git des modifications.
- **Objective**: Créer un commit git propre et documenté pour la Phase 1 de la Roadmap VI.
- **Deliverables**:
  - `GitCommitPhase1`: Commit git des fichiers modifiés.
- **Process**: `git add` et `git commit`.

## Acceptance Criteria:
### Step 1:
- `kepler_sundman.solve_kepler_universal` résout sans exception ni troncature les pas rétrogrades avec $z$ descendant jusqu'à $-100\,000$.
- Wronskien unimodulaire $|f \dot{g} - \dot{f} g - 1| < 10^{-15}$ vérifié sur ces configurations hyperboliques rétrogrades.
### Step 2:
- Erreur relative sur $c_k'(z)$ inférieure à $10^{-14}$ sur l'ensemble de l'intervalle $z \in [-0.2, 0.2]$ pour $k=0, 1, 2, 3$.
### Step 3:
- Élimination de toute inversion artificielle de moment cinétique dans `mcgehee.py` ; invariance de la parité vérifiée.
### Step 4:
- 100 % des tests unitaires de `tests/unit/solvers/` passent avec succès.
- Taux de réussite global de 100 % sur l'ensemble de la bibliothèque (zéro régression).
### Step 5:
- Commit git créé avec message explicatif référençant l'Itération 144 et les anomalies du Chapitre 3.

## Observations:
> [!IMPORTANT]
> Cette Phase 1 de la Roadmap VI assainit le fondement numérique de la bibliothèque. L'élimination du clamping $z < -100$ débloque la robustesse de Blanes-Moan BM6 lors de frôlements rapides à pas rétrogrades.
