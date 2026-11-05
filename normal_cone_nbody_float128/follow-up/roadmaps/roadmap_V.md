# Strategic Roadmap : Création et Certification de la Bibliothèque Symplectique Académique `normal_cone_nbody_float64` (Roadmap V)

## Strategic Overview:
- **Objectif Stratégique Général** : Développer, vérifier et certifier de bout en bout la bibliothèque de mécanique céleste symplectique `normal_cone_nbody_float64`, accompagnant la monographie de recherche dans [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/).
- **Protocole Opérationnel Strict à 5 Étapes pour Chaque Itération** :
  1. *Étape Théorique* : Identification et relecture méthodique des chapitres, sections, définitions et théorèmes du papier LaTeX [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/).
  2. *Étape d'Audit du Code Hérité* : Relecture approfondie des modules correspondants dans [`models/normal_cone_integrator/cas6_nbody3d_async/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/models/normal_cone_integrator/cas6_nbody3d_async/) en s'appuyant sur [`python_codebase_mapping_and_audit_report.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/python_codebase_mapping_and_audit_report.md) ; identification des composants conformes, des reliquats de nommage (`cas6`) à bannir et des parties manquantes.
  3. *Étape de Spécification Précise* : Rédaction d'un document de spécification claire du code à écrire, de l'API modulaire, des invariants de référence et des tests unitaires avec assertions et tolérances strictes.
  4. *Étape d'Implémentation Modulaire & Optimisée* : Rédaction du code source avec en-tête descriptif exhaustif référençant les sections de l'article de recherche, décomposition en fonctions pures au rôle limpide, abstractions ciblées (`Protocol`, `dataclasses`) et conception prête pour l'optimisation bas-niveau (vectorisation NumPy, zéro allocation superflue).
  5. *Étape de Test & Validation* : Renseignement et exécution des tests unitaires miroirs (`tests/unit/`) avec taux de succès impératif de 100 % avant tout passage à l'itération suivante.
- **Contraintes Techniques & Inviolables** :
  - Sanctuarisation absolue en lecture seule du répertoire [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/).
  - Zéro adoucissement ni troncature artificielle ($\varepsilon = 0, d_0 = 0$).
  - Éradication totale du terme `"cas6"` et respect des dénominations de Fulton--MacPherson et Calaque--Safronov.
  - Conservation discrète exacte des 10 générateurs de Lie-Bargmann au seuil IEEE 754 ($< 10^{-14}$).
- **Document Maître de Référence** : [`normal_cone_nbody_float64_master_architecture.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/normal_cone_nbody_float64_master_architecture.md).

---

## Iteration 112 Strategy:
### Iteration 112: Planification Stratégique de la Roadmap V et Structuration du Cycle de Vie d'Implémentation
- **Scope & Objectives**:
  - Définir la structure globale de Roadmap V (itérations 112 à 119) selon le protocole granulaire strict à 5 étapes.
  - Mettre en place les contrats de planification formels (`roadmap.md`, `implementation_plan.md`, `task.md`) conformément aux directives du skill `planning-contracts`.
  - Recueillir la validation formelle de l'utilisateur sur la feuille de route avant le premier pas d'implémentation.
- **Dependencies & Continuity**: S'appuie sur le document maître finalisé et les 3 audits multi-agents de l'Itération 111.
- **Gating Criteria**:
  - `roadmap.md` créé et aligné sur les 7 couches de la bibliothèque.
  - `implementation_plan.md` et `task.md` d'Itération 112 rédigés en français sans balises HTML.
  - Validation explicite par l'utilisateur du découpage en 7 itérations de développement.

---

## Iteration 113 Strategy:
### Iteration 113: Couche 1 (`core/`) - Algèbre Fondamentale, État Étendu $\mathcal{P}_{\mathrm{ext}}$, Invariants et Arbres de Jacobi
- **Scope & Objectives**:
  - Implémenter les modules du socle fondamental : `core/constants.py`, `core/state.py`, `core/invariants.py`, `core/jacobi.py`, `core/coordinates.py`, `core/transitions.py`.
  - Relecture théorique : Chapitres 1, 2 et Sections 3.1–3.4 (variété $\mathcal{P}_{\mathrm{ext}}$, groupe de Galilée étendu, treillis des partitions $\Pi(N)$ et compactification de Fulton--MacPherson).
  - Relecture du code hérité : `cas6_nbody3d_async/core/` (nettoyage de 26 occurrences de `cas6`, découplage `JacobiNode`/état dynamique, adoption de $\bar{\rho}_C \in [L]$).
  - Implémentation des tests unitaires miroirs : `tests/unit/core/test_*.py` (6 fichiers de test).
- **Dependencies & Continuity**: Première brique d'exécution logicielle, pose les fondations de données pour toutes les couches supérieures.
- **Gating Criteria**:
  - 100 % de succès sur les tests de `tests/unit/core/`.
  - Exactitude symplectique Darboux $\mathbf{J} \mathbf{\Omega} \mathbf{J}^T = \mathbf{\Omega}$ à $< 10^{-15}$.
  - Déterminant unimodulaire $\det \mathbf{M}_{\mathrm{coord}} \equiv +1.0$ sur les rotations d'arbres.

---

## Iteration 114 Strategy:
### Iteration 114: Couche 2 (`solvers/`) - Séries de Stumpff, Flot de Sundman-Halley, Rebond de McGehee et Super-Nœuds
- **Scope & Objectives**:
  - Implémenter les solveurs singuliers : `solvers/stumpff.py`, `solvers/kepler_sundman.py`, `solvers/mcgehee.py`, `solvers/supernode.py`.
  - Relecture théorique : Sections 3.3–3.5, Section 4.4, Appendice A (dissémination de la barrière de Sheng--Suzuki via les fonctions de Stumpff, collision frontale sans régularisation ad hoc, involution unitaire $\mathcal{B}_C^2 = \operatorname{Id}_{\mathbf{M}_C}$, horloge dépliée de Sundman $\mathrm{d}t_S = \chi_S^2\,\mathrm{d}\tau_S$).
  - Relecture du code hérité : `cas6_nbody3d_async/solvers/` (intégration de la convergence cubique de Halley avec dérivées secondes complètes, réversibilité temporelle exacte $\Delta t < 0$).
  - Implémentation des tests unitaires miroirs : `tests/unit/solvers/test_*.py` (4 fichiers de test).
- **Dependencies & Continuity**: Requiert les structures d'état et les constantes de la Couche 1 (`core/`).
- **Gating Criteria**:
  - Involution stricte $\mathcal{B}_C^2 \equiv \operatorname{Id}_{\mathbf{M}_C}$ à l'égalité binaire près.
  - Unimodularité $|f \dot{g} - \dot{f} g - 1| < 10^{-15}$ du flot universel de Kepler.
  - Régularisation anti-Zénon prouvée sur collision frontale à moment cinétique nul ($L=0$).

---

## Iteration 115 Strategy:
### Iteration 115: Couche 3 (`potential/`) - Crible Conormal $\Theta(u_v)$, Potentiel d'Apex $\mathcal{C}^{1,1}$, Détection de Marée de Hill et Newton Direct
- **Scope & Objectives**:
  - Implémenter les potentiels et détecteurs : `potential/screening.py`, `potential/screened_apex.py`, `potential/supernode_detector.py`, `potential/direct.py`.
  - Relecture théorique : Section 4.2–4.3, Lemme 4.1, Théorème 4.2 (crible conormal lisse $\Theta(u_v)$ sur ratio adimensionné, bornes de Lipschitz universelles $M_{\mathrm{apex}}, L_{\mathrm{apex}} < \infty$, seuil naturel de marée de Hill $\eta_S \le 3^{-1/3}$).
  - Relecture du code hérité : `cas6_nbody3d_async/potential/` (suppression de l'ancien `apex.py`, transition de $\delta=0.5$ au critère de Hill formel).
  - Implémentation des tests unitaires miroirs : `tests/unit/potential/test_*.py` (4 fichiers de test).
- **Dependencies & Continuity**: Requiert la géométrie de Jacobi de Couche 1 et alimentera les champs de force des intégrateurs.
- **Gating Criteria**:
  - Gradient analytique conforme aux différences finies à $< 10^{-11}$.
  - Annulation rigoureuse des pôles fantômes multipolaires en coïncidence barycentrique ($u_v \le 1 \implies \Theta = 0$).
  - Invariance d'échelle exacte du détecteur de Hill : $\eta_S(\lambda \mathbf{q}) = \eta_S(\mathbf{q})$.

---

## Iteration 116 Strategy:
### Iteration 116: Couche 4 (`integrators/`) - Suite Symplectique Modulaire $\mathcal{S}_2, \mathcal{S}_4, \mathrm{BM}_4^{(7)}, \mathrm{BM}_6^{(11)}$, Correcteurs $\mathcal{C}_{2p}$ et Intégrateurs Traités $\mathcal{P}\mathrm{BM}$
- **Scope & Objectives**:
  - Implémenter les schémas symplectiques : `integrators/protocol.py`, `integrators/strang.py`, `integrators/yoshida.py`, `integrators/blanes_moan_4.py`, `integrators/blanes_moan_6.py`, `integrators/correctors.py`, `integrators/processed.py`.
  - Relecture théorique : Chapitre 5 (compositions palindromiques de Blanes--Moan, forme FSAL, équations homologiques des correcteurs $\{H, K\} = H_{2p+1}$, signe de l'opérateur inverse $\exp(+h^{2p} X_K)$, télescopage dyadique).
  - Relecture du code hérité : `cas6_nbody3d_async/integrators/` (développement *de novo* de `BlanesMoan6Integrator` avec les coefficients à 25 décimales et de `ProcessedSymplecticIntegrator`).
  - Implémentation des tests unitaires miroirs : `tests/unit/integrators/test_*.py` (7 fichiers de test).
- **Dependencies & Continuity**: Assemble les solveurs de Couche 2 et les champs de potentiel de Couche 3 sous l'interface typée `SymplecticIntegrator`.
- **Gating Criteria**:
  - Ordres de convergence effectifs vérifiés numériquement : 2, 4, 6 et 8.
  - Télescopage dyadique parfait des correcteurs : $\mathcal{C}^{-1} \circ \mathcal{C} = \operatorname{Id}$ avec zéro évaluation de force supplémentaire dans l'intérieur du pas.
  - Respect du critère de Suzuki sur les sous-pas rétrogrades ($w_i < 0$).

---

## Iteration 117 Strategy:
### Iteration 117: Couche 5 (`system/`) - Orchestrateur `NBodySystem`, Atlas à Hystérésis, Ordonnanceur Dyadique Multi-Taux et Historique
- **Scope & Objectives**:
  - Implémenter l'orchestration système : `system/nbody_system.py`, `system/atlas_manager.py`, `system/dyadic_scheduler.py`, `system/history.py`.
  - Relecture théorique : Section 4.5, Section 5.6 (transitions d'atlas symplectiques selon l'hystérésis $\gamma = 0.6$, temps de séjour minimal $\tau_{\min} > 0$, pas dyadiques synchronisés $h / 2^m$, backward error analysis de Skeel & Gear).
  - Relecture du code hérité : `cas6_nbody3d_async/system/` (conception du planificateur dyadique formel et de la classe unifiée `NBodySystem`).
  - Implémentation des tests unitaires miroirs : `tests/unit/system/test_*.py` (4 fichiers de test).
- **Dependencies & Continuity**: Intègre toutes les couches précédentes au sein de l'objet d'état unifié.
- **Gating Criteria**:
  - Absence absolue de broutement topologique (impossibilité géométrique $1 < 0.36$).
  - Conservation symplectique asynchrone $\mathrm{d}\Omega_{\mathrm{async}} = 0$ sur les subdivisions dyadiques.
  - Zéro fuite mémoire sur $10^6$ pas avec `HistoryTracker`.

---

## Iteration 118 Strategy:
### Iteration 118: Couche 6 (`diagnostics/`) - Métrologie des 6 KPIs Célestes, Radar d'Anticipation et Visualiseur 10 Panneaux
- **Scope & Objectives**:
  - Implémenter les outils métrologiques et graphiques : `diagnostics/metrics.py`, `diagnostics/radar.py`, `diagnostics/plotting.py`.
  - Relecture théorique : Section 5.3, Chapitre 6 (les 6 critères célestes canoniques, régression de dérive de phase $\mathcal{O}(h^p t)$, pente séculaire d'énergie $\mathrm{Slope}_H \equiv 0$, distance aux diviseurs exceptionnels $E_I \subset \partial\mathbf{Q}[N]$).
  - Relecture du code hérité : `cas6_nbody3d_async/diagnostics/` (harmonisation des visualisations sans serveur graphique, rendu 300 DPI publication).
  - Implémentation des tests unitaires miroirs : `tests/unit/diagnostics/test_*.py` (3 fichiers de test).
- **Dependencies & Continuity**: Fournit les instruments d'évaluation pour la validation de la Couche 7.
- **Gating Criteria**:
  - Détection prédictive des collisions multiples par le radar à l'approche de $\partial\mathbf{Q}[N]$.
  - Tracé autonome 10 panneaux conforme au standard des figures du papier.
  - Calcul certifié du résidu de troncature de Lie-Bargmann $\Delta_{\mathrm{trunc}}$.

---

## Iteration 119 Strategy:
### Iteration 119: Couche 7 & Certification Globale (`benchmarks/`) - Les 8 Problèmes Canoniques et Courbes Pareto Comparatives
- **Scope & Objectives**:
  - Implémenter le banc de simulation : `benchmarks/registry.py`, `benchmarks/runner.py`.
  - Relecture théorique : Chapitres 6 et 7 de l'article (les 8 configurations de test : collision frontale, Burrau, Arenstorf, Figure-8, Sitnikov, Système jovien extérieur, etc.).
  - Relecture du code hérité : `cas6_nbody3d_async/benchmarks/`.
  - Implémenter et exécuter la suite de benchmarks comparatifs (`tests/simulation_benchmarks/`) confrontant les 6 intégrateurs ($\mathcal{S}_2, \mathcal{S}_4, \mathrm{BM}_4^{(7)}, \mathrm{BM}_6^{(11)}, \mathcal{P}\mathrm{BM}_4^{(6)}, \mathcal{P}\mathrm{BM}_6^{(8)}$).
  - Génération des courbes Pareto travail-précision pour publication.
- **Dependencies & Continuity**: Étape terminale de certification de l'ensemble de la bibliothèque logicielle.
- **Gating Criteria**:
  - Les 8 benchmarks exécutés avec succès sur les 6 schémas symplectiques.
  - Supériorité de $\mathrm{BM}_6^{(11)}$ et $\mathcal{P}\mathrm{BM}_6^{(8)}$ démontrée expérimentalement sur la conservation d'énergie et la dérive de phase.
  - Respect à 100 % de l'ensemble des 28 tests unitaires et des 8 bancs d'essais.

---

## Multi-Iteration Checklist:
- [x] Iteration 112: Planification Stratégique de la Roadmap V et Structuration du Cycle de Vie d'Implémentation
- [x] Iteration 113: Couche 1 (`core/`) - Algèbre Fondamentale, État Étendu $\mathcal{P}_{\mathrm{ext}}$, Invariants et Arbres de Jacobi
- [x] Iteration 114: Couche 2 (`solvers/`) - Séries de Stumpff, Flot de Sundman-Halley, Rebond de McGehee et Super-Nœuds
- [x] Iteration 115: Couche 3 (`potential/`) - Crible Conormal $\Theta(u_v)$, Potentiel d'Apex $\mathcal{C}^{1,1}$, Détection de Marée de Hill et Newton Direct
- [x] Iteration 116: Couche 4 (`integrators/`) - Suite Symplectique Modulaire $\mathcal{S}_2, \mathcal{S}_4, \mathrm{BM}_4^{(7)}, \mathrm{BM}_6^{(11)}$, Correcteurs $\mathcal{C}_{2p}$ et Intégrateurs Traités $\mathcal{P}\mathrm{BM}$
- [x] Iteration 117: Couche 5 (`system/`) - Orchestrateur `NBodySystem`, Atlas à Hystérésis, Ordonnanceur Dyadique Multi-Taux et Historique
- [x] Iteration 118: Couche 6 (`diagnostics/`) - Métrologie des 6 KPIs Célestes, Radar d'Anticipation et Visualiseur 10 Panneaux
- [x] Iteration 119: Couche 7 & Certification Globale (`benchmarks/`) - Les 8 Problèmes Canoniques et Courbes Pareto Comparatives
- [x] Iteration 120: Audit Multi-Agents Contradictoire de la Bibliothèque Symplectique (13 sous-agents)
- [x] Iteration 121: Intégration des Recommandations d'Audit Multi-Agents (Correctifs P0, P1 et P2)
- [x] Iteration 122: Démonstration Comparative de Surperformance & Visualisations Dynamiques Multi-Échelles
- [x] Iteration 123: Audit et Re-étalonnage des Profils de Tolérance Céleste KPI
- [x] Iteration 124: Audit Multi-Agent Épistémique du Flux d'Exécution et Couverture de Tests (3 sous-agents)
- [x] Iteration 125: Remédiation Complète des 7 Points d'Audit et Éradication de la Rupture à 4 Corps
- [x] Iteration 126: Régénération Complète des Planches 300 DPI et Synchronisation Définitive des Données Webapp
- [x] Iteration 127: Enrichissement Webapp Démonstrateur, Topologie SVG et Comparaison Traités PBM4/6
- [x] Iteration 128: Résolution du Bug d'Affichage Dynamique Webapp
- [x] Iteration 129: Audit Archéologique et Remédiation du Saut d'Énergie Burrau (Contradiction 3-Agents)
- [x] Iteration 130: Remédiation du Plateau 10^-8 sur Burrau par Anticipation Dyadique Képlérienne
- [x] Iteration 131: Remédiation du Saut d'Énergie PBM4 sur le Système Hiérarchique à 4 Corps
- [x] Iteration 132: Audit Contradictoire Approfondi de la Conformité Théorique sur 5 Axes (3 sous-agents)
- [x] Iteration 133: Harmonisation Théorique Complète des 9 Couples et Élimination des 8 Compromis d'Audit
- [x] Iteration 134: Régénération des Planches de Publication 300 DPI et du Dataset Interactif Conormal
- [x] Iteration 135: Diagnostic Télémetrique et Archéologique du Saut à 4 Corps (4 sous-agents)
- [x] Iteration 136: Éradication des Seuils Arbitraires, Écrantages Hermite C2-C6 et Élimination du Lemme 4.3
- [x] Iteration 137: Écrantage Heptadécique C8 pour Schémas d'Ordre 8 (PBM6) et Tridécique C6 pour PBM4
- [x] Iteration 138: Régénération des Planches et Appli Interactive sous Hiérarchie C6/C8 et Dyadique Dynamique
- [x] Iteration 139: Conception de la Planche Dédiée aux Invariants de Lie-Bargmann Multi-Intégrateurs
- [x] Iteration 140: Aération des Titres et Marges Hautes des Planches 300 DPI
- [x] Iteration 141: Aération Profonde des Marges Supérieures (> 400px à 300 DPI)
- [x] Iteration 142: Étalonnage Typographique Analytique et Complétude Pareto des Schémas Traités
- [x] Iteration 143: Audit Contradictoire Multi-Agentique Exhaustif des Chapitres 3, 4 et 5 (9 Auditeurs)
