# Implementation Plan of Iteration 114 : Couche 2 (solvers/) - Séries de Stumpff, Flot de Sundman-Halley, Rebond de McGehee et Super-Nœuds Synchrones
Created at: 2026-10-07T17:38:00+02:00

## Context:
1. **Explicit Request**: Implémenter et certifier unitairement la Couche 2 (`solvers/`) de la bibliothèque `normal_cone_nbody_float64`, conformément à la feuille de route stratégique Roadmap V. Le développement doit rigoureusement respecter le protocole récurrent à 5 étapes :
   - Étape 1 : Identifier et relire les sections théoriques du papier (Sections 3.3–3.5, 4.4 et Appendice A de [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/)).
   - Étape 2 : Relire l'ancien code associé (`models/normal_cone_integrator/cas6_nbody3d_async/solvers/` et modules associés) et identifier ce qui est conforme, les mentions de `cas6` à éradiquer, et les composants manquants (Halley complet avec dérivées secondes, gestion rétrograde $\Delta t < 0$).
   - Étape 3 : Établir un document de spécification claire du code à écrire avec les tests et attendus.
   - Étape 4 : Rédiger le code modulaire optimisé en fonctions pures, avec typage strict et en-tête descriptif citant les théorèmes du papier de recherche et les preuves Lean 4.
   - Étape 5 : Rédiger et exécuter la suite de tests unitaires miroirs (`tests/unit/solvers/`) avec un taux de réussite impératif de 100 %.
2. **Implicit Strategic Value**: Fournir les moteurs analytiques fondamentaux pour la résolution exacte du mouvement à deux corps non adouci et la régularisation de collision sans paramètre de ramollissement ($\varepsilon = 0, d_0 = 0$). Ces solveurs constituent le cœur cinématique sur lequel reposeront les étapes de dérive (Drift) des intégrateurs de Blanes--Moan et le franchissement déterministe des singularités.
3. **Alignment Criteria**: Respect de l'identité du Wronskien $|f \dot{g} - \dot{f} g - 1| < 10^{-15}$, réversibilité temporelle exacte $\Phi(-\Delta t) \circ \Phi(\Delta t) = \operatorname{Id}$, involution stricte du rebond $\mathcal{B}_C^2 = \operatorname{Id}_{\mathbf{M}_C}$ à la précision machine, et absence totale de blocage zénonien sur les super-nœuds.

## Guidelines:
- **Method**: Application stricte du protocole à 5 étapes sans court-circuit.
- **Red Lines**:
  - Sanctuarisation absolue en lecture seule du répertoire [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/).
  - Zéro occurrence de `"cas6"` ou `"cas 6"`.
  - Zéro adoucissement ni seuil de distance artificielle ($\varepsilon = 0, d_0 = 0$).
  - Aucun couplage ascendant vers les couches 3, 4, 5, 6 ou 7. Dépendances limitées exclusivement à la Couche 1 (`core/`).
- **Architecture**:
  - `normal_cone_nbody_float64/solvers/stumpff.py` : fonctions universelles $c_0(z)$ à $c_5(z)$ avec développement de Horner d'ordre 6 pour $|z| < 0.1$ et relations différentielles récursives.
  - `normal_cone_nbody_float64/solvers/kepler_sundman.py` : solveur universel de Kepler-Sundman via itération cubique de Halley avec première et seconde dérivées exactes, propagation par $f, g, \dot{f}, \dot{g}$ et support rétrograde $\Delta t < 0$.
  - `normal_cone_nbody_float64/solvers/mcgehee.py` : coordonnées de blow-up sur la variété de collision $\mathbf{M}_C$ et involution analytique $\mathcal{B}_C$ préservant la mesure symplectique.
  - `normal_cone_nbody_float64/solvers/supernode.py` : solveur régularisé de super-nœud synchrone sous horloge dépliée $\mathrm{d}t_S = \chi_S^2 \mathrm{d}\tau_S = \rho_S \mathrm{d}\tau_S$ garantissant $N_{\mathrm{sub}} < \infty$.
- **Code Headers**: En-têtes conformes au standard académique avec citations explicites des sections de la monographie et des théorèmes formels Lean 4.

## Plan:
### Step 1: Relecture Théorique des Sections 3.3–3.5, 4.4 et Appendice A
- **Context**: La monographie détaille la régularisation de Sundman, les fonctions de Stumpff, la variété de McGehee et le franchissement de la barrière de Sheng--Suzuki.
- **Objective**: Extraire l'ensemble des équations analytiques exactes : coefficients des séries de Stumpff, formule d'itération cubique de Halley avec dérivée seconde complète $F''(\chi)$, expressions des fonctions $f, g, \dot{f}, \dot{g}$, involution de McGehee $\mathcal{B}_C$ et horloge dépliée de Sundman.
- **Deliverables**:
  - `Layer2TheorySynthesis`: Synthèse des équations analytiques et des théorèmes certifiés Lean 4 associés.
- **Process**: Consultation ciblée des chapitres 3 et 4 et de l'annexe mathématique de la monographie.

### Step 2: Audit Approfondi du Code Hérité (`cas6_nbody3d_async/solvers/`)
- **Context**: L'ancien code comporte des implémentations de Stumpff, Kepler conormal et amas hiérarchiques.
- **Objective**: Auditer `stumpff.py`, `conormal_kepler.py` et `hierarchical_cluster.py` pour cartographier le code réutilisable, éliminer le marquage `cas6`, corriger les lacunes (dérivées d'ordre 2 pour Halley, pas rétrogrades stables).
- **Deliverables**:
  - `Layer2RefactoringMatrix`: Matrice de migration et de purification des algorithmes de résolution.
- **Process**: Inspection systématique des sources et confrontation avec le rapport d'audit.

### Step 3: Spécification Formelle et Contrats d'Interface de la Couche 2
- **Context**: Fixer les signatures d'API, structures de données et critères d'invariance avant tout codage.
- **Objective**: Définir les interfaces de `stumpff.py`, `kepler_sundman.py`, `mcgehee.py` et `supernode.py`, ainsi que les tolérances des tests unitaires miroirs.
- **Deliverables**:
  - `Layer2ModuleSpecification`: Spécification technique complète des 4 modules de solveurs.
- **Process**: Rédaction des spécifications d'invariance et des contrats d'interface.

### Step 4: Implémentation Modulaire et Optimisée des 4 Modules de `solvers/`
- **Context**: Le répertoire cible `normal_cone_nbody_float64/solvers/` doit être créé.
- **Objective**: Écrire les 4 modules Python avec en-têtes académiques exhaustifs, vectorisation NumPy et algorithmes de résolution optimisés.
- **Deliverables**:
  - `normal_cone_nbody_float64/solvers/__init__.py`
  - `normal_cone_nbody_float64/solvers/stumpff.py`
  - `normal_cone_nbody_float64/solvers/kepler_sundman.py`
  - `normal_cone_nbody_float64/solvers/mcgehee.py`
  - `normal_cone_nbody_float64/solvers/supernode.py`
- **Process**: Implémentation modulaire séquentielle sans écrasement destructif.

### Step 5: Implémentation et Exécution des Tests Unitaires Miroirs (`tests/unit/solvers/`)
- **Context**: Chaque module possède son fichier de test unitaire miroir dédié.
- **Objective**: Implémenter et exécuter la suite de tests sous `pytest` avec 100 % de succès.
- **Deliverables**:
  - `tests/unit/solvers/test_stumpff.py`
  - `tests/unit/solvers/test_kepler_sundman.py`
  - `tests/unit/solvers/test_mcgehee.py`
  - `tests/unit/solvers/test_supernode.py`
- **Process**: Exécution via `pytest` et vérification des critères de précision machine.

## Acceptance Criteria:
### Step 1:
- Recensement complet des références théoriques : séries de Stumpff (Sec. 3.3), Halley cubique (Sec. 3.3.2), McGehee (Sec. 3.5), super-nœuds anti-Zénon (Sec. 4.4).
### Step 2:
- Rapport d'audit de refactorisation attestant la purification de tout terme `cas6` et l'intégration des dérivées secondes complètes de Halley.
### Step 3:
- Spécifications formelles précisant les types d'entrée/sortie et les bornes d'erreur attendues.
### Step 4:
- 4 modules de solveurs implémentés avec en-têtes complets citant l'article et Lean 4.
- Zéro adoucissement ni découpage artificiel ($\varepsilon = 0, d_0 = 0$).
### Step 5:
- 100 % de réussite sur `pytest tests/unit/solvers/`.
- Unimodularité du Wronskien vérifiée à $|f \dot{g} - \dot{f} g - 1| < 10^{-15}$.
- Réversibilité temporelle validée : $\|\Phi(-\Delta t) \circ \Phi(\Delta t)(\mathbf{r}, \mathbf{w}) - (\mathbf{r}, \mathbf{w})\| < 10^{-14}$.
- Involution de McGehee $\mathcal{B}_C^2 \equiv \operatorname{Id}$ à l'égalité exacte.
- Régularisation anti-Zénon validée sur collision frontale à moment cinétique nul ($L=0$).

## Observations:
> [!IMPORTANT]
> L'itération de Halley doit obligatoirement inclure le terme de courbure complet issu de la dérivée seconde $F''(\chi) = \sigma_0 c_0(\alpha\chi^2) + (\mu - \alpha r_0)\chi c_1(\alpha\chi^2) - \alpha \sigma_0 \chi^2 c_2(\alpha\chi^2)$ afin de garantir la vitesse de convergence cubique $\mathcal{O}(|\chi - \chi^*|^3)$ même en régime de quasi-collision frontale où $r_0 \to 0$.

> [!NOTE]
> La résolution des sous-pas rétrogrades ($\Delta t < 0$) est indispensable pour le bon fonctionnement des intégrateurs de Blanes--Moan d'ordres élevés ($\mathrm{BM}_4^{(7)}, \mathrm{BM}_6^{(11)}$), dont certains coefficients $w_i$ sont négatifs. Le solveur de Kepler-Sundman doit être rigoureusement symétrique par rapport au renversement du temps.
