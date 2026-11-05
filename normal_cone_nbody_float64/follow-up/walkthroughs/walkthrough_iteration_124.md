# Walkthrough de l'Itération 124 : Multi-Audit Sous-Agentique du Flot d'Intégration, de la Complétude des Tests et de l'Isomorphisme Théorique

## Summary:

Dans cette Itération 124, en réponse directe à la demande de l'utilisateur :  
*« Y a encore des fails, et le solveur semble pas tout tester correctement, lance des audit sous agentique pour comprendre le flot du code et le comparer au papier de recherche et au code d'avant »*,  
nous avons déployé un **collège d'audit contradictoire indépendant composé de 3 sous-agents experts** :

1. **Sous-Agent 1 (`76eec236`) — Flot d'Exécution & Mécanique d'Intégration Dynamique** : Traçage chronologique du pipeline `NBodySystem.step()` $\to$ `DyadicScheduler` $\to$ `SymplecticIntegrator` $\to$ `JacobiCoordinates` $\to$ `Supernodes` $\to$ `ScreenedApex`.
2. **Sous-Agent 2 (`9320740e`) — Isomorphisme Mathématique avec la Monographie LaTeX** : Confrontation terme à terme avec les Chapitres 3, 4 et 5 de `redaction/unsoftened_conormal_nbody3d/` et les lemmes formellement certifiés en Lean 4.
3. **Sous-Agent 3 (`53abf383`) — Parité avec le Code Hérité, Couverture de Tests & Diagnostic des Échecs** : Rapprochement ligne à ligne avec `models/normal_cone_integrator/cas6_nbody3d_async/`, audit critique des 123 tests unitaires et reproduction pas à pas du mécanisme de divergence numérique.

Les trois rapports ont été consolidés dans l'artéfact de synthèse majeur :  
[`multi_agent_flow_and_testing_audit.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/multi_agent_flow_and_testing_audit.md).

---

## Critical Discoveries (Les 7 Découvertes Clés de l'Audit) :

1. **Isolation Expérimentale du Mécanisme de Blowup ($+78\,000\,\%$) à $t \approx 3.3$** :
   - À $t = 3.285$, l'arbre de Jacobi pivote et qualifie le triplet $(1, 2, 3)$ en super-nœud synchrone.
   - Dans `step_drift_and_kepler` (ligne 67), comme la paire de feuilles $(1, 2)$ appartient au super-nœud, **le solveur universel de Kepler-Sundman est sauté**.
   - La paire est confiée à `step_supernode_synchronous` (`solvers/supernode.py`), qui applique un Leapfrog cartésien direct non régularisé à pas fixe $h = \mathrm{d}t / n_{\mathrm{steps}}$, **sans aucune horloge conormale de Sundman dépliée** ($\mathrm{d}t_S = \chi_S^2\,\mathrm{d}\tau_S$).
   - Au sous-pas 13, $d_{12}$ chute à $2.37 \times 10^{-3}$ ; la force non régularisée $1/d_{12}^2 \approx 1.8 \times 10^5$ explose, faisant sauter l'énergie totale de $H = -3.20$ à $+2489.97$ !
2. **Gel Cinématique du Centre de Masse des Super-Nœuds** :
   - Le nœud composite reliant le super-nœud au reste de la topologie est sauté dans le drift cinétique.
   - Dans `step_supernode_synchronous`, les positions sont actualisées via `q_new = x_cur + q_cdm`, mais $\mathbf{q}_{\mathrm{cdm}}$ reste statiquement figé. Le super-nœud ne translate pas avec sa vitesse orbitale globale autour des autres corps.
3. **Fusion Indiscriminée de Tous les Super-Nœuds** :
   - Dans `nbody_system.py`, tous les clusters détectés sont fusionnés en un unique ensemble `sn_set: Set[int]`. Deux sous-amas physiquement disjoints sont ainsi amalgamés en un seul objet barycentrique.
4. **Désactivation Silencieuse du Potentiel d'Apex Écranté $\Theta(u_v)$** :
   - Dans `potential/screened_apex.py` (ligne 297), `step_apex_kick` impose `include_screened_monopole=False`. Le terme monopole $\mathcal{C}^{1,1}$ $W_v = \Theta(u_v) \frac{G M_L M_R}{\|\mathbf{r}_v\|}$ (Théorème 4.2 du papier) est entièrement neutralisé à l'exécution.
5. **Violation de l'Algèbre RKN dans $\mathrm{BM}_6^{(11)}$** :
   - L'ordre 6 de Blanes--Moan repose sur l'hypothèse $[B, [B, A]] = 0$. Comme l'opérateur $\Phi$ inclut les résolutions képlériennes non-linéaires dépendant de $q$ et $p$, $[B, [B, A]] \ne 0$. L'ordre 6 effectif dégénère en régime de marée intense, alors que $\mathrm{BM}_4^{(7)}$ conserve son ordre 4 inconditionnel (composition de Strang symétrique).
6. **Miscalibration Systématique des KPIs dans `BenchmarkRunner.run()`** :
   - `runner.py` (ligne 120) appelle `compute_celestial_kpis()` sans spécifier de profil contextuel, retombant sur `"smooth_periodic"`. Ce profil rejette toute pente séculaire $> 10^{-11}$, marquant `symplecticity_passed: False` (FAIL) sur Burrau alors que son erreur d'énergie est excellente ($6.9 \times 10^{-9}$).
7. **Angles Morts de la Suite de Tests Unitaires** :
   - Les 123 tests unitaires n'exécutent que des cas quasi-statiques. Aucun test n'intègre dynamiquement une transition d'atlas sous mouvement réel (l'unique test d'atlas asserte explicitement `assert len(transition_history) == 0`).
   - `apply_mcgehee_bounce` dans `solvers/mcgehee.py` est du code mort non connecté au solveur.

---

## Actionable Remediation Roadmap :

| Priorité | Cible | Action Corrective Précise |
| :---: | :--- | :--- |
| **P0.1** | `solvers/supernode.py` | Intégrer la véritable horloge conormale de Sundman dépliée $\mathrm{d}t_S = \chi_S^2\,\mathrm{d}\tau_S$ dans la boucle de sous-pas et préserver le solveur képlérien de feuille interne. |
| **P0.2** | `integrators/common.py` | Dériver le centre de masse du super-nœud $\mathbf{q}_{\mathrm{cdm}} \leftarrow \mathbf{q}_{\mathrm{cdm}} + \mathrm{d}t \cdot (\mathbf{p}_{\mathrm{cdm}} / M_S)$ et itérer sur `List[Tuple[int, ...]]` sans aplatissement global. |
| **P0.3** | `potential/screened_apex.py` | Activer `include_screened_monopole=True` dans `step_apex_kick` avec le gradient analytique complet du Théorème 4.2. |
| **P0.4** | `potential/supernode_detector.py` | Supprimer la condition `if state.N <= 3` pour permettre la qualification de super-nœud sur les collisions centrales à $N \ge 4$. |
| **P0.5** | `benchmarks/runner.py` | Associer à chaque benchmark son profil de tolérance adapté (`"chaotic_encounter"`, `"hierarchical_multiscale"` ou `"smooth_periodic"`). |
| **P1.1** | `integrators/correctors.py` | Corriger la symplecticité de $\mathcal{C}_4, \mathcal{C}_6$ en appliquant la dérive des impulsions $\mathbf{p} \leftarrow \mathbf{p} - \delta \nabla_{\mathbf{q}} K$. |
| **P1.2** | `tests/unit/` | Ajouter les tests unitaires dynamiques d'atlas flip sous intégration et de collision quadruple sous super-nœud. |
