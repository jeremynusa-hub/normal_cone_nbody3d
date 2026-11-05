# Implementation Plan of Iteration 131 : Audit Contradictoire Approfondi et Remédiation des Correcteurs Symplectiques Traités (PBM) sur les Systèmes Multi-Corps à Transitions d'Atlas
Created at: 2026-10-07T20:50:00+02:00

## Context:
1. **Explicit Request**:
   - « Dans le cas à 4 corpsn je vois un saut d'nergie pour PBM4 et pas pour BM4 ou BM6 qui 'lencadrent, c'est étrange non ? Lance 2 audit pour tester la bonne implémentation des correcteur symplectic, elle me semble douteuse »
2. **Implicit Strategic Value**:
   - Élucider rigoureusement le comportement des intégrateurs symplectiques traités (Processed Blanes-Moan : PBM4 et PBM6) lors de simulations astrophysiques multi-corps subissant des transitions d'atlas de Jacobi $\mathrm{SL}(2, \mathbb{R})$ répétées.
   - Comprendre et éradiquer la cause mathématique ou algorithmique pour laquelle PBM4 produit un saut d'énergie ou une divergence (`nan`) sur le 4-corps hiérarchique alors que BM4 ($1.37 \times 10^{-6}$) et BM6 ($3.65 \times 10^{-5}$) demeurent stables.
   - Aligner strictement l'implémentation des opérateurs de correction $C_{2p}(h)$ et du télescopage frontière avec les spécifications théoriques de la Section 5.5 du manuscrit de recherche.
3. **Alignment Criteria**:
   - Deux rapports d'audit indépendants et étanches (Théorie & Manuscrit vs Flot d'Exécution & Mécanique Numérique) analysent en profondeur les équations homologiques, l'action des correcteurs sur les coordonnées, et la compatibilité entre la conjugaison $\Psi(h) = C(h) \circ S(h) \circ C(h)^{-1}$ et les changements de topologie d'arbre.
   - La pathologie de PBM4 sur le scénario à 4 corps est formellement diagnostiquée et éliminée.
   - 100% des tests unitaires passent (`pytest tests/unit/`).

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des scénarios.
- **Règle d'Édition** : Utiliser exclusivement `replace_file_content` pour les modifications du code Python, ne jamais écraser avec `overwrite: true`.

## Plan:
### Step 1: Déploiement de 2 Audits Contradictoires Spécialisés
- **Audit 1 (Théoricien du Manuscrit & Fondements Symplectiques)** :
  - Examiner la Section 5.5 de `redaction/unsoftened_conormal_nbody3d/chapters/ch05_higher_order_symplectic/content.tex`.
  - Auditer la formulation du générateur hamiltonien $K_{2p}$, les crochets $\{H, K_{2p}\} = H_{2p+1}$, le signe et la parité de l'opérateur inverse $C^{-1}$, ainsi que l'hypothèse de base sous-jacente du télescopage dyadique ($C^{-1} \circ C = \mathrm{Id}$).
  - Analyser l'impact d'une permutation ou d'une transition d'atlas de Jacobi survenant au cours d'un pas macro ou entre l'application de $C$ et de $C^{-1}$.
- **Audit 2 (Profileur Numérique & Flot d'Exécution des Correcteurs)** :
  - Profiler pas à pas l'exécution de PBM4 et PBM6 sur le problème à 4 corps (`scratch/run_showcase_simulations.py` et `system/nbody_system.py`).
  - Examiner comment `ProcessedSymplecticIntegrator` dans `integrators/processed.py` et `integrators/correctors.py` calcule les déplacements de position $\delta q$ et d'impulsion $\delta p$.
  - Identifier le pas précis où PBM4 diverge ou produit `nan`, inspecter les valeurs des forces de correction, les dénominateurs des distances mutuelles et l'interaction avec le `DyadicScheduler`.

### Step 2: Synthèse Contradictoire et Diagnostic Formel de la Pathologie
- Confronter les conclusions théoriques et numériques des deux sous-agents.
- Établir sans équivoque si le défaut provient :
  1. D'une formule erronée ou d'une singularité dans le générateur de correction $K_4$ de `correctors.py` lors d'approches multi-corps complexes.
  2. D'une rupture de conjugaison lors d'un switch d'atlas (si l'état entre dans le pas avec l'atlas $A$, subit un switch vers $B$ à l'intérieur du pas macro, et sort avec $C_B^{-1}$ au lieu de $C_A^{-1}$, brisant la symplecticité).
  3. D'une incompatibilité entre le télescopage dyadique des sous-pas ($h / 2^k$) et l'application du correcteur sur le macro-pas global.

### Step 3: Conception et Implémentation du Correctif
- Corriger les modules incriminés (`integrators/correctors.py`, `integrators/processed.py`, ou `system/nbody_system.py`) via `replace_file_content`.
- Assurer la stabilité absolue de PBM4 et PBM6 sans dégrader leur ordre effectif.

### Step 4: Contrôle Qualité, Validation Numérique et Walkthrough
- Valider la suite de tests unitaires (`pytest tests/unit/`).
- Vérifier la trajectoire et la conservation d'énergie de PBM4 sur le cas à 4 corps.
- Régénérer `showcase_data.json` et synchroniser `interactive_conormal_showcase.html`.

## Acceptance Criteria:
- **Critère 1** : L'anomalie de PBM4 sur le problème à 4 corps est résolue ; PBM4 présente une erreur stable et bornée comparable ou supérieure à BM4.
- **Critère 2** : L'implémentation des correcteurs symplectiques est formellement certifiée conforme à la Section 5.5 du manuscrit de recherche et robuste aux transitions d'atlas.
- **Critère 3** : 100% des tests unitaires passent sans régression (`pytest tests/unit/`).

## Observations:
> [!IMPORTANT]
> Un intégrateur traité (Processed Integrator) repose sur la conjugaison $\Psi(h) = C(h) \circ S(h) \circ C(h)^{-1}$. Lorsque le système d'évolution change dynamiquement de carte d'atlas ou de coordonnées de Jacobi au cours de la trajectoire, la commutativité du télescopage frontière $C^{-1} \circ C = \mathrm{Id}$ peut être rompue si la métrique de transition n'est pas strictement symplectique par rapport à l'opérateur de Lie de $K_{2p}$.
