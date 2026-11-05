# Implementation Plan of Iteration 118 : Couche 6 (`diagnostics/`) - Métrologie des 6 KPIs Célestes, Radar d'Anticipation et Visualiseur 10 Panneaux
Created at: 2026-10-07T18:06:00+02:00

## Context:
1. **Explicit Request**: Développer, tester unitairement de façon exhaustive et certifier formellement la Couche 6 (`diagnostics/`) de la bibliothèque `normal_cone_nbody_float64` en suivant scrupuleusement le protocole mandaté en 5 étapes.
2. **Implicit Strategic Value**: Fournir la suite instrumentale d'analyse quantitative, de métrologie symplectique et de visualisation publication-ready indispensable pour certifier les benchmarks de l'Itération 119, en mesurant rigoureusement la conservation des premiers intégrales, l'anticipation prédictive des collisions multiples aux diviseurs exceptionnels $E_I \subset \partial\mathbf{Q}[N]$, et la restitution graphique automatisée sans dépendance graphique interactive.
3. **Alignment Criteria**: Taux de réussite de 100 % sur l'ensemble de la suite de tests unitaires de la couche `tests/unit/diagnostics/` et cumulatif (115+ tests passés), absence totale de marqueurs `"cas6"`, zéro adoucissement ($\varepsilon = 0, d_0 = 0$), extraction autonome des 6 KPIs célestes canoniques, et génération sans faille de planches 10 panneaux vectorielles / matricielles à 300 DPI sous backend `Agg`.

## Guidelines:
- **Method**: Protocole rigoureux en 5 étapes (Théorie -> Audit hérité -> Spécification formelle -> Implémentation modulaire -> Tests unitaires 1:1).
- **Red Lines**: Interdiction absolue de toute modification dans `redaction/unsoftened_conormal_nbody3d/`. Éradication totale de `"cas6"` ou `"cas 6"`. Zéro adoucissement ni seuillage artificiel de distance ($\varepsilon = 0, d_0 = 0$).
- **Architecture**:
  - `diagnostics/collision_radar.py` : Détecteur prédictif de proximité des diviseurs exceptionnels de collision de Fulton--MacPherson $E_I$, estimant le temps de collision balistique $t_{\mathrm{coll}} = \min_{i<j} r_{ij} / |\dot{r}_{ij}|$ et signalant les configurations à haut risque de singularité.
  - `diagnostics/metrics.py` : Calculateur unifié des 6 KPIs célestes canoniques (dérive d'énergie absolue et relative, pente séculaire $\mathrm{Slope}_H$, conservation vectorielle de $\mathbf{P}$ et $\mathbf{L}$, dérive du centre de masse $\mathbf{G}$, erreur du Wronskien de Kepler, et temps de séjour d'atlas $\tau_{\mathrm{dwell}}$).
  - `diagnostics/visualizer.py` : Générateur de planches académiques multi-panneaux (10 panneaux d'analyse orbitale et symplectique) fonctionnant en mode non-interactif (`matplotlib.use('Agg')`).
  - `diagnostics/__init__.py` : Façade et exportations du module.
- **Form**: Code, commentaires, docstrings et tests exclusivement en anglais avec citations de la monographie et des théorèmes Lean 4. Plan et liste de tâches en français markdown natif.

## Plan:
### Step 1: Relecture Théorique des Chapitres 5 et 6 de la Monographie
- **Starting Material**: Monographie LaTeX `redaction/unsoftened_conormal_nbody3d/` (Section 5.3 sur la backward error analysis et l'absence de dérive séculaire séculaire $\mathrm{Slope}_H = 0$, Chapitre 6 sur les diviseurs exceptionnels $E_I$ et la métrologie céleste).
- **Objective**: Extraire les définitions formelles des 6 KPIs célestes, la formulation de la régression linéaire de dérive d'énergie, la distance riemannienne aux diviseurs exceptionnels de collision et le théorème de non-divergence séculaire.
- **Deliverables**: Spécification théorique consolidée pour la Couche 6 (`Component Theoretical Framework`).
- **Process**: Relecture méthodique des fichiers `.tex` et formalisation des métriques analytiques.

### Step 2: Audit Approfondi du Code Hérité (`cas6_nbody3d_async/diagnostics/`)
- **Starting Material**: Répertoire hérité `models/normal_cone_integrator/cas6_nbody3d_async/diagnostics/` (notamment `collision_radar.py`, `energy_monitor.py`, `dashboard.py`), et le rapport d'audit `python_codebase_mapping_and_audit_report.md`.
- **Objective**: Recenser les algorithmes validés de calcul de temps à la collision, identifier et éliminer les 18 occurrences résiduelles de `cas6`, supprimer tout appel bloquant `plt.show()` au profit du backend `Agg`, et encapsuler les calculs métriques dans des fonctions pures et structures `dataclass`.
- **Deliverables**: Rapport de refactorisation Couche 6 (`Component Legacy Audit`).
- **Process**: Confrontation ligne à ligne et élimination des dépendances d'affichage GUI directes.

### Step 3: Spécification Formelle et Contrats d'Interface de la Couche 6
- **Starting Material**: Spécification globale `normal_cone_nbody_float64_master_architecture.md` (Section 3.6 & 6.6).
- **Objective**: Définir les types, signatures, exceptions, invariants et critères d'évaluation des modules de `diagnostics/` et de leurs tests miroirs.
- **Deliverables**: Contrats d'interface diagnostique (`Component Interface Contracts`).
- **Process**: Formalisation des classes `CollisionRadar`, `CelestialKPIs`, `SimulationDiagnostics`, `PublicationVisualizer`.

### Step 4: Implémentation Modulaire et Optimisée des Modules de `diagnostics/`
- **Starting Material**: Modules certifiés des Couches 1 à 5 (`core/`, `solvers/`, `potential/`, `integrators/`, `system/`).
- **Objective**: Implémenter les 4 fichiers sous `normal_cone_nbody_float64/diagnostics/` :
  1. `diagnostics/__init__.py` : Exports du package.
  2. `diagnostics/collision_radar.py` : Calcul vectorisé de la matrice des distances mutuelles $r_{ij}$, vitesse relative radiale $v_{r, ij} = (\mathbf{q}_i - \mathbf{q}_j) \cdot (\mathbf{v}_i - \mathbf{v}_j) / r_{ij}$, estimation balistique du temps de première approche $t_{\mathrm{approach}}$, et classification des régimes de proximité singulière.
  3. `diagnostics/metrics.py` : Calculateur des 6 KPIs fondamentaux, extraction de séries temporelles depuis `HistoryTracker`, régression linéaire des pentes séculaires ($\mathrm{d}(\Delta H)/\mathrm{d}t$), et diagnostic global `KPIReport`.
  4. `diagnostics/visualizer.py` : Visualiseur 10 panneaux exportant des planches publiables (PNG / PDF 300 DPI) illustrant l'évolution spatiale, les invariants de Lie--Bargmann, les transitions d'atlas, les pas dyadiques et la distance minimale.
- **Deliverables**: Fichiers de production sous `normal_cone_nbody_float64/diagnostics/`.
- **Process**: Écriture modulaire avec typage strict PEP 484/544, vectorisation NumPy intégrale et rendu graphique headless garanti.

### Step 5: Implémentation et Exécution des Tests Unitaires Miroirs (`tests/unit/diagnostics/`)
- **Starting Material**: Spécifications d'évaluation et modules de `diagnostics/`.
- **Objective**: Développer et exécuter avec 100 % de succès les 3 bancs de test unitaires miroirs :
  1. `tests/unit/diagnostics/test_collision_radar.py` (détection exacte de paires en collision, estimation du temps d'impact, invariance par translation).
  2. `tests/unit/diagnostics/test_metrics.py` (calcul des 6 KPIs sur trajectoire canonique, nullité prouvée de la pente séculaire $\mathrm{Slope}_H < 10^{-12}$, intégrité du rapport).
  3. `tests/unit/diagnostics/test_visualizer.py` (génération headless sans erreur X11, validation de la présence des 10 panneaux, sauvegarde en mémoire tampon ou fichier temporaire).
- **Deliverables**: 3 fichiers de test sous `tests/unit/diagnostics/` avec 100 % de réussite.
- **Process**: Exécution via `pytest -v`, vérification des tolérances et validation graphique.

## Acceptance Criteria:
- **Critère 6.1 (Radar de Collision Prédictif)** : Pour une paire de corps en chute libre radiale l'un vers l'autre avec vitesse relative constante, le temps d'impact calculé $t_{\mathrm{approach}}$ est exact à $< 10^{-14}$.
- **Critère 6.2 (Calcul des 6 KPIs Célestes)** : La structure `KPIReport` calcule fidèlement les 6 grandeurs : dérive relative d'énergie, pente séculaire $\mathrm{Slope}_H$, dérive d'impulsion $\Delta P$, dérive de moment cinétique $\Delta L$, Wronskien moyen et temps de séjour d'atlas $\tau_{\mathrm{dwell}}$.
- **Critère 6.3 (Pente Séculaire Symplectique)** : Sur une orbite quasi-périodique de Kepler intégrée avec $\mathrm{BM}_4$ ou $\mathrm{BM}_6$, la pente séculaire ajustée $\mathrm{Slope}_H = |\mathrm{d}(\Delta H)/\mathrm{d}t|$ est inférieure à $10^{-11}$, prouvant l'absence de dérive d'énergie séculaire artificielle.
- **Critère 6.4 (Rendu Visuel Headless Garanti)** : `PublicationVisualizer` génère l'ensemble des 10 figures en mémoire sans nécessiter de serveur X11 ni d'interface graphique fenêtrée, produisant un fichier d'image valide vérifié par lecture binaire.
- **Critère 6.5 (Zero-Softening & Zero-Cas6)** : 0 occurrence de `"cas6"` ou `"cas 6"`, 0 softening ($\varepsilon = 0, d_0 = 0$).
- **Critère 6.6 (Taux de Succès Global)** : 100 % de réussite sur `tests/unit/diagnostics/` et sur la suite cumulative (> 115 tests).

## Observations:
> [!IMPORTANT]
> Le module de visualisation doit forcer le backend non-interactif `matplotlib.use('Agg')` dès l'import pour garantir une exécution parfaite dans les environnements CI/CD ou sans écran.
> [!NOTE]
> Le radar d'anticipation de collision ne modifie en rien la dynamique physique : c'est un instrument métrologique pur observant la géométrie de configuration.
