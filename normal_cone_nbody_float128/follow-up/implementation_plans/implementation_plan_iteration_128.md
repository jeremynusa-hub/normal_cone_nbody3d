# Implementation Plan of Iteration 128 : Correction de la Persistance Dynamique de l'Arbre de Jacobi et du Graphe d'Énergie en Lecture Continue
Created at: 2026-10-07T20:15:00+02:00

## Context:
1. **Explicit Request**:
   - « les arbre de jacobi et l'evolution comparé d'nergtie disparaissent quand on met play. »
2. **Implicit Strategic Value**:
   - Garantir une expérience visuelle et pédagogique fluide et sans accroc dans la vitrine conormale `interactive_conormal_showcase.html`.
   - Permettre à l'utilisateur d'observer en direct le comportement dynamique complet lors de l'animation : trajectoires orbitales, restructuration de l'arbre de Jacobi aux 17 basculements d'atlas, activation des super-nœuds compacts de Hill, et défilement synchronisé du curseur temporel sur le comparateur d'énergie multi-intégrateurs sans effacement des courbes ni du graphe SVG.
3. **Root Cause Analysis**:
   - Dans `interactive_conormal_showcase.html`, la fonction `runAnimationLoop()` invoque à chaque tick `renderAll()`, laquelle fait appel à `getScenarioData()`.
   - La fonction `getScenarioData()` reconstruisait un objet synthétique pour chaque scénario en **omettant** les propriétés critiques `tree_history` et `integrator_energies`.
   - En conséquence, lors du playback :
     1. `scData.tree_history` devenait `undefined`, assignant `treeSnap = null`. `renderJacobiTree` exécutait alors `svg.innerHTML = ''; return;`, vidant instantanément l'arbre de Jacobi de la page.
     2. `scData.integrator_energies` devenait `undefined`, de sorte que `integData` valait `{}`. `renderEnergyComparisonChart` effaçait le canvas via `ctx.clearRect(0, 0, w, h)` et ne traçait aucune courbe.
   - De surcroît, l'absence de garde-fous sur l'échantillonnage temporel (`tree_history.length` vs `positions.length`), l'absence de mémorisation du dernier arbre valide (`lastValidTreeSnap`), et l'absence de blocs `try/catch` rendaient l'animation vulnérable à toute exception non interceptée.

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des scénarios.
- **Règle d'Édition** : Utiliser exclusivement `replace_file_content` pour les modifications du code HTML/JS, ne jamais écraser avec `overwrite: true`.

## Plan:
### Step 1: Correction Fondamentale de `getScenarioData()` dans `interactive_conormal_showcase.html`
- **Fichier** : [`interactive_conormal_showcase.html`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/interactive_conormal_showcase.html)
- **Objectif** :
  - Transmettre systématiquement `tree_history` et `integrator_energies` dans le retour de `getScenarioData()` pour les 4 scénarios (`figure8`, `burrau`, `radial`, `hierarchical_4body`).
  - Lier correctement `transitions` et `supernodes` pour chaque scénario.
- **Processus** : Édition ciblée via `replace_file_content`.

### Step 2: Sécurisation et Robustesse de `renderJacobiTree` et `renderEnergyComparisonChart`
- **Fichier** : [`interactive_conormal_showcase.html`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/interactive_conormal_showcase.html)
- **Objectif** :
  - Dans `renderAll` : borner l'index de snapshot `treeSnap` via `Math.min(Math.max(0, currentStep), scData.tree_history.length - 1)`.
  - Dans `renderJacobiTree` : introduire une variable `lastValidTreeSnap` de repli, ne jamais effacer le SVG si aucun arbre n'est disponible, sécuriser les calculs de coordonnées de nœuds enfants (`leftCoord`, `rightCoord`), et encapsuler le rendu dans un bloc `try/catch`.
  - Dans `renderEnergyComparisonChart` : vérifier que les dimensions du canvas sont strictement positives (`w > 0 && h > 0`), borner le curseur temporel à `[0, numSteps - 1]`, et encapsuler dans un bloc `try/catch`.
- **Processus** : Édition ciblée via `replace_file_content`.

### Step 3: Validation Automatisée et Contrôle de Non-Régression
- **Objectif** :
  - Exécuter la suite de tests unitaires complète (`pytest tests/unit/`).
  - Vérifier par un script de test headless/statique que la structure HTML et les scripts JS s'exécutent sans erreur de syntaxe ni référence indéfinie.
  - Documenter la résolution dans `notes.md` et produire le `walkthrough.md` de clôture.
- **Processus** : Tests unitaires, inspection de code, validation.

## Acceptance Criteria:
- **Critère 1** : L'arbre de Jacobi reste continuellement affiché et se met à jour en temps réel à chaque pas de simulation pendant la lecture continue (`isPlaying = true`).
- **Critère 2** : Le graphe comparatif d'énergie multi-intégrateurs conserve toutes ses courbes visibles pendant la lecture continue, et le curseur vertical se déplace harmonieusement avec le pas de temps sans effacer les tracés.
- **Critère 3** : Le comportement est parfaitement stable et testé sur les 4 scénarios (Figure-8, Burrau, Collision radiale, 4-corps hiérarchique).
- **Critère 4** : 100% des tests unitaires (`pytest tests/unit/`) passent avec succès.
