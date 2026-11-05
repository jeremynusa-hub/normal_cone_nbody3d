# Implementation Plan of Iteration 153 : Harmonisation Mathématique et Typographique de la Monographie
Created at: 2026-10-08T12:35:00+02:00

## Context:
1. **Explicit Request**: L'utilisateur poursuit la Roadmap VII vers la Phase 4 (Itération 153) consacrée à l'harmonisation mathématique et typographique du manuscrit de recherche LaTeX (`redaction/unsoftened_conormal_nbody3d/`), suite aux conclusions sans appel de l'audit contradictoire à 9 sous-agents.
2. **Implicit Strategic Value**: Assurer un isomorphisme parfait et une vérité épistémique irréprochable entre le formalisme théorique du papier de recherche et le comportement effectif du moteur de calcul certifié. La correction de la coquille typographique sur la dérivée de Stumpff garantit la rigueur mathématique formelle, tandis que la clarification épistémique sur $\mathrm{PBM}_4$ dissipe toute ambiguïté théorique concernant l'élimination des crochets de Lie de $\mathcal{L}_5$ et justifie la supériorité de $\mathrm{PBM}_6$.
3. **Alignment Criteria**:
   - Levée autorisée et strictement circonscrite de la sanctuarisation de `redaction/unsoftened_conormal_nbody3d/` pour les deux sections ciblées (`sec02_suzuki_barrier_mcgehee` et `sec05_symplectic_correctors` du Chapitre 5).
   - Zéro modification sur les autres fichiers ou chapitres du manuscrit.
   - Préservation stricte de la syntaxe LaTeX, des labels de référence croisée et des badges Lean 4.
   - 100% des tests unitaires de la bibliothèque Python (148/148) maintenus au vert.

## Guidelines:
- **Method**:
  - Utiliser `replace_file_content` pour des modifications chirurgicales et ciblées dans les fichiers `.tex`.
  - Dans `redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/sec02_suzuki_barrier_mcgehee/content.tex`, rétablir le dénominateur $z$ dans le premier terme de la dérivée de Stumpff : $c_k'(z) = \frac{c_{k-1}(z) - k c_k(z)}{2z}$.
  - Dans `redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/sec05_symplectic_correctors/content.tex`, expliciter dans une remarque académique dédiée et ajuster l'énoncé du Théorème 5.7 pour stipuler que $\mathrm{PBM}_4$ est un schéma d'ordre 4 à constante d'erreur réduite de 35% (projection scalaire optimale de l'équation homologique), l'élimination simultanée des 6 crochets de Hall indépendants de $\mathcal{L}_5$ n'étant accessible qu'aux noyaux de composition d'ordre supérieur comme $\mathrm{BM}_6^{(11)}$ conjugué à $\mathcal{C}_6$.
- **Red Lines**:
  - Ne jamais introduire d'adoucissement ($\varepsilon = 0, d_0 = 0$).
  - Ne pas altérer la structure globale des chapitres ni rompre les références `\ref` ou `\eqref`.
  - Aucune régression dans la suite de tests unitaires.
- **Architecture**:
  - Respecter scrupuleusement la typographie et le style LaTeX du document (environnements `thmbox`, `rembox`, `align`, notations conormales).

## Plan:
### Step 1: Rectification de la Formule Différentielle de Stumpff dans le Chapitre 5 (sec02)
- **Context**: À la ligne 82 de `redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/sec02_suzuki_barrier_mcgehee/content.tex`, l'expression imprimée $c_k'(z) = \frac{1}{2} c_{k-1}(z) - \frac{k}{2z} c_k(z)$ omet par inadvertance typographique le facteur $1/z$ sur le premier terme.
- **Objective**: Corriger la formule en $c_k'(z) = \frac{c_{k-1}(z) - k c_k(z)}{2z}$.
- **Deliverables**:
  - `StumpffDerivativeFix`: Édition chirurgicale dans `sec02_suzuki_barrier_mcgehee/content.tex`.
- **Process**: Appliquer `replace_file_content` sur la ligne 82 pour insérer l'expression rigoureuse.

### Step 2: Clarification Épistémique de l'Ordre de Convergence de $\mathrm{PBM}_4$ dans le Chapitre 5 (sec05)
- **Context**: L'audit théorique et numérique a prouvé qu'un unique champ hamiltonien scalaire $X_{K_4}$ ne possède qu'un degré de liberté ($\kappa_4$) et ne peut annuler simultanément les 6 crochets de Hall indépendants de $\mathcal{L}_5$. $\mathrm{PBM}_4$ conserve une convergence d'ordre 4 avec une constante d'erreur réduite de 35%, tandis que l'ordre effectif élevé est conféré par $\mathrm{PBM}_6$.
- **Objective**: Clarifier le statut de $\mathrm{PBM}_4$ dans `sec05_symplectic_correctors/content.tex` via une remarque théorique rigoureuse explicitant la géométrie de Lie-Poisson et l'obstruction dimensionnelle de $\mathcal{L}_5$.
- **Deliverables**:
  - `PBM4EpistemicClarification`: Précision textuelle et remarque théorique dans `sec05_symplectic_correctors/content.tex`.
- **Process**: Appliquer `replace_file_content` pour intégrer la clarification homologique et la distinction entre $\mathrm{PBM}_4$ (ordre 4 avec $-35\%$ d'erreur) et $\mathrm{PBM}_6$ (ordre effectif supérieur).

### Step 3: Contrôle Syntaxique et Cohérence Globale du Manuscrit LaTeX
- **Context**: Les fichiers LaTeX modifiés doivent être exempts d'erreurs de balises, de délimiteurs de formules ou d'environnements mal fermés.
- **Objective**: Vérifier l'intégrité syntaxique des deux fichiers modifiés.
- **Deliverables**:
  - `LaTeXSyntaxAuditReport`: Validation de la balance des délimiteurs et des environnements.
- **Process**: Script d'analyse syntaxique Python vérifiant la fermeture des environnements et accolades.

### Step 4: Tests Unitaires de Non-Régression
- **Context**: La bibliothèque Python doit rester en parfait état de fonctionnement (148/148 tests passés).
- **Objective**: Exécuter `pytest tests/unit/`.
- **Deliverables**:
  - `UnitTestsCertification`: 148/148 tests validés.
- **Process**: Lancement de la suite de tests.

### Step 5: Git Commit de Phase 4
- **Context**: Phase 4 de la Roadmap VII achevée.
- **Objective**: Enregistrer un commit atomique documenté.
- **Deliverables**:
  - `GitCommitPhase4`: Commit git documenté sur `main`.
- **Process**: `git add` et `git commit`.

## Acceptance Criteria:
### Step 1:
- L'équation dans `sec02_suzuki_barrier_mcgehee/content.tex` affiche explicitement $c_k'(z) = \frac{c_{k-1}(z) - k c_k(z)}{2z}$ sans modifier le reste du texte de la preuve.
### Step 2:
- Le texte de `sec05_symplectic_correctors/content.tex` explicite que $\mathrm{PBM}_4$ réduit la constante d'erreur de 35% tout en maintenant un ordre 4 asymptotically strict en raison de la dimension 6 de l'espace des crochets de $\mathcal{L}_5$.
### Step 3:
- Les deux fichiers `.tex` ont 100% de concordance syntaxique (balises, accolades, environnements).
### Step 4:
- 100% des tests unitaires (148/148) passent avec succès.
### Step 5:
- Commit git atomique enregistré sur `main`.

## Observations:
> [!IMPORTANT]
> L'honnêteté épistémique est un principe directeur fondamental. Requalifier $\mathrm{PBM}_4$ comme un intégrateur d'ordre 4 optimisé (avec réduction de 35% de la constante d'erreur de Blanes-Moan) renforce considérablement la crédibilité académique de la monographie face aux auditeurs et pairs, tout en mettant en valeur la nécessité de $\mathrm{PBM}_6$ pour franchir les ordres supérieurs.
