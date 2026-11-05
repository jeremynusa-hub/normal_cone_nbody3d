# Implementation Plan of Iteration 161 : Phase 2 — Mise à Niveau du Noyau Formel Lean 4 (Règle de Leibniz Complète pour $\mathcal{C}_6$)
Created at: 2026-10-09T08:37:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur a validé l'arbitrage pré-itératif pour la Phase 2 en retenant :
   - `Option A` : Isomorphisme intégral monolithique dans Lean 4 — injection directe du terme variationnel quadratique complet de Leibniz dans `makeConormalGenerator6.grad_pi` au sein de `02_CorrectorBM6_Theorems.lean` :
     $$\nabla_{\boldsymbol{\pi}} K_6 = \kappa (1 + \|\boldsymbol{\pi}\|^2)\mathbf{F} + 2\kappa(\mathbf{F}\cdot\boldsymbol{\pi})\boldsymbol{\pi},$$
     et adaptation de l'ensemble des théorèmes et preuves constructives (`conormal_generator6_is_odd`, majoration lipschitzienne conormale et certificat maître `cone_normal_corrector_bm6_certified`) sans aucun axiome ad-hoc ni `sorry`.
2. **Implicit Strategic Value**: Établir un isomorphisme mathématique et formel absolu entre le noyau de théorèmes certifiés en Lean 4 et les lignes 256-257 de `normal_cone_nbody_float64/integrators/correctors.py`. Garantir que la dérivation formelle du champ hamiltonien de Lie associé au correcteur d'ordre 6 $\mathcal{C}_6$ de Blanes--Moan reflète exactement la règle de dérivation du produit (Leibniz) sur la composante d'impulsion sans approximation ni omission de terme.
3. **Alignment Criteria**:
   - `makeConormalGenerator6.grad_pi` dans [`02_CorrectorBM6_Theorems.lean`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/sec05_symplectic_correctors/code/02_CorrectorBM6_Theorems.lean) intègre les deux composantes vectorielles : le cisaillement $\kappa (1 + \|\boldsymbol{\pi}\|^2)\mathbf{F}$ et la projection quadratique $2\kappa(\mathbf{F}\cdot\boldsymbol{\pi})\boldsymbol{\pi}$.
   - Le théorème de parité temporelle `conormal_generator6_is_odd` est formellement démontré en prouvant que le terme quadratique vérifie $(-1) \times (-1) = +1$ sous inversion d'impulsion $\boldsymbol{\pi} \mapsto -\boldsymbol{\pi}$.
   - Les bornes lipschitziennes de `ConeNormalCorrectorBM6ProblemData` et le certificat maître `cone_normal_corrector_bm6_certified` sont rigoureusement maintenus avec zéro `sorry` et zéro axiome custom.
   - Les 151 tests unitaires Python demeurent à 100 % de succès.

## Guidelines:
- **Standards Formels** : Respect strict du skill `lean-proof-craft` : zéro `sorry`, zéro axiome ad-hoc, preuves constructives modulaires $\le 30$ lignes, utilisation des tactiques canoniques de Mathlib 4 (`ring`, `linarith`, `nlinarith`, `norm_num`, `ext`, `dsimp`).
- **Édition de Code Source** : Interdiction absolue de `write_to_file` avec `overwrite: true` sur tout fichier source existant. Utilisation exclusive de `replace_file_content` ciblée sur les blocs sémantiques minimaux.
- **Isomorphisme Code-Théorie** : Concordance 1:1 bit-exacte entre la définition formelle Lean 4 et l'implémentation opérationnelle Python.

## Plan:

### Step 1: Formalisation Algébrique du Vecteur de Leibniz dans Vec3
- **Context**: [`02_CorrectorBM6_Theorems.lean`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/sec05_symplectic_correctors/code/02_CorrectorBM6_Theorems.lean) définit la structure `Vec3` et ses lemmes algébriques (`dot_comm`, `dot_neg_right`, `smul_neg`, etc.).
- **Objective**: Fournir les lemmes élémentaires dans `Vec3` sur le comportement du terme bilinéaire $(\mathbf{F} \cdot (-\mathbf{v})) (-\mathbf{v}) = (\mathbf{F} \cdot \mathbf{v}) \mathbf{v}$ et l'addition vectorielle sous symétrie temporelle pour faciliter la preuve de parité.
- **Deliverables**:
  - `Vec3LeibnizLemmas`: Lemmes `smul_dot_neg_neg` ou simplification directe dans `Vec3`.
- **Process**: Modifier `02_CorrectorBM6_Theorems.lean` par `replace_file_content`.

### Step 2: Enrichissement de `makeConormalGenerator6` et Démonstration de Parité
- **Context**: `makeConormalGenerator6` (lignes 189-202) définit `grad_pi` avec seulement le premier terme.
- **Objective**:
  - Remplacer `grad_pi` par la somme vectorielle exacte des deux termes de Leibniz :
    `Vec3.add (Vec3.smul (kappa * (1 + pi_sq)) (apex.grad_V p)) (Vec3.smul (kappa * 2 * Vec3.dot (apex.grad_V p) p.pi) p.pi)`.
  - Adapter la preuve de `conormal_generator6_is_odd` pour certifier l'invariance temporelle de `grad_pi`.
- **Deliverables**:
  - `LeibnizConormalGenerator6`: Structure `makeConormalGenerator6` complète.
  - `ParityTheoremAdapted`: Preuve constructive de `conormal_generator6_is_odd` validée.
- **Process**: Modifier `02_CorrectorBM6_Theorems.lean` par `replace_file_content`.

### Step 3: Majoration Lipschitzienne et Adaptation du Certificat Maître
- **Context**: `makeCertifiedBM6Data` et `cone_normal_corrector_bm6_certified` établissent la majoration lipschitzienne de `grad_pi` (`hK_grad_pi_bound`).
- **Objective**:
  - Ajuster la constante de majoration $M_{\pi}$ dans `makeCertifiedBM6Data` pour borner analytiquement la norme de `grad_pi` avec les deux termes (terme direct + terme projectif).
  - Certifier la preuve constructive de `hK_grad_pi_bound` et le certificat maître `cone_normal_corrector_bm6_certified`.
- **Deliverables**:
  - `CertifiedBM6DataAdapted`: Données certifiées conformes avec Leibniz.
  - `MasterCertificateBM6Certified`: Théorème `cone_normal_corrector_bm6_certified` prouvé sans `sorry`.
- **Process**: Modifier `02_CorrectorBM6_Theorems.lean` par `replace_file_content`.

### Step 4: Validation de Compilation Formelle et Contrôle de Non-Régression
- **Context**: Le fichier Lean 4 doit compiler sans erreur ni avertissement.
- **Objective**:
  - Vérifier la validité formelle du fichier Lean 4.
  - Exécuter la suite de tests unitaires Python `pytest tests/unit/` pour certifier la non-régression.
  - Créer le commit Git de clôture de la Phase 2.
- **Deliverables**:
  - `Lean4Verification`: Fichier compilé avec succès sans `sorry`.
  - `PytestNonRegression`: 151/151 tests unitaires validés.
  - `GitCommitPhase2`: Commit Git de la Phase 2.
- **Process**: Exécuter les commandes de validation et `git commit`.

## Acceptance Criteria:
### Step 1:
- `Vec3` dispose des identités algébriques requises pour factoriser $(\mathbf{F} \cdot (-\mathbf{v})) (-\mathbf{v}) = (\mathbf{F} \cdot \mathbf{v}) \mathbf{v}$.
### Step 2:
- `makeConormalGenerator6.grad_pi` intègre rigoureusement les deux termes $\kappa (1 + \|\boldsymbol{\pi}\|^2)\mathbf{F} + 2\kappa(\mathbf{F}\cdot\boldsymbol{\pi})\boldsymbol{\pi}$.
- `conormal_generator6_is_odd` est prouvé de manière constructive sans `sorry`.
### Step 3:
- `makeCertifiedBM6Data` majore `grad_pi` par une constante analytique $M_{\pi} > 0$ rigoureusement démontrée.
- Le certificat maître `cone_normal_corrector_bm6_certified` est clos sans `sorry` ni axiome custom.
### Step 4:
- Zéro `sorry` et zéro axiome custom dans `02_CorrectorBM6_Theorems.lean`.
- 151 / 151 tests unitaires Python validés.
- Commit Git créé avec succès.

## Observations:
> [!IMPORTANT]
> L'ajout du terme variationnel quadratique de Leibniz $2\kappa(\mathbf{F}\cdot\boldsymbol{\pi})\boldsymbol{\pi}$ dans $\nabla_{\boldsymbol{\pi}} K_6$ respecte rigoureusement la parité temporelle impaire de $K_6$ car $(-1) \times (-1) = +1$. Par suite, le crochet de Lie-Poisson $\{H, K_6\}$ demeure de parité strictement paire en impulsion, garantissant la résolubilité exacte de l'équation homologique $\{H, K_6\} = H_7$ face au résidu d'ordre 6 de Blanes--Moan.
