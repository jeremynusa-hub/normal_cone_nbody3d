# Implementation Plan of Iteration 177 : Forecast à 1 an du Système Solaire Étendu (46 Corps) avec PBM6 et Triple Optimisation & Confrontation NASA JPL Horizons
Created at: 2026-10-09T13:17:30+02:00

## Context:
1. **Explicit Request**: L'utilisateur demande : « Reprends la configuration du système solaire avec les 46 corps de la dernière fois. Et forecast 1 an des positions de ces corps avec pbm6 et toute nos optimisation. Compare les résultats avec la nasa. »
2. **Implicit Strategic Value**: Valider la pertinence céleste et l'exactitude physique de notre moteur conormal Fulton-MacPherson à haute précision d'ordre 8 (PBM6) couplé à la triple optimisation (Apex vectorisé sans allocation, solveur Kepler-Sundman universel inliné et super-nœuds synchrones JIT / LLVM) sur une échelle astrophysique réelle multi-échelles (des comètes plongeantes et lunes intérieures jusqu'aux confins transneptuniens comme Sedna), et attester de la fidélité d'intégration en confrontant directement les positions prédites à $+1\text{ an}$ aux éphémérides officielles de référence de la NASA (JPL Horizons DE440).
3. **Alignment Criteria**: Intégrer les 46 corps sur 365,25 jours (époque initiale 2026-09-29 -> époque finale 2027-09-29) sans aucun adoucissement ($\varepsilon = 0, d_0 = 0$), préserver strictement les 10 invariants du groupe de Lie-Bargmann, interroger l'API REST officielle NASA JPL Horizons pour récupérer les vecteurs d'état exacts des 46 corps à l'époque finale, calculer les écarts euclidiens en AU et en kilomètres, et produire un tableau récapitulatif détaillé et certifié ainsi qu'une planche de diagnostic 300 DPI.

## Guidelines:
- **Method**: Charger l'état étendu à 46 corps depuis [`models/solar_system/data/solar_system_extended_state.json`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/models/solar_system/data/solar_system_extended_state.json) ($JD = 2461312.5$). Utiliser `GAUSSIAN_G = 0.00029591220828559115\text{ AU}^3 / (M_\odot \cdot d^2)$. Configurer l'arbre de Jacobi hiérarchique optimal via `AtlasManager` et exécuter sous `ProcessedSymplecticIntegrator(BlanesMoan6Integrator(G=GAUSSIAN_G))` avec `apex_evaluator=OptimizedScreenedApexPotential`, `kepler_solver=OptimizedKeplerSolver`, et `supernode_evaluator=JitSupernodeSolver`.
- **Red Lines**: Aucun adoucissement physique ($\varepsilon = 0, d_0 = 0$). Conservation rigoureuse de la masse et des constantes astronomiques officielles IAU/JPL DE440.
- **Comparison Baseline**: Vecteurs de position barycentriques ICRF/J2000 (`@0`) issus de l'API REST officielle NASA JPL Horizons à l'époque 2027-09-29 00:00:00 TDB.

## Plan:
### Step 1: Préparation du module d'orchestration et pont Système Solaire pour la bibliothèque conormale
- **Context**: L'état initial des 46 corps est sérialisé dans `models/solar_system/data/solar_system_extended_state.json`. La bibliothèque `normal_cone_nbody_float64` utilise `ExtendedPhaseState` et `NBodySystem`.
- **Objective**: Concevoir un script dédié `scratch/forecast_solar_system_46bodies_1yr.py` initialisant `NBodySystem` avec les 46 corps, `G=GAUSSIAN_G`, l'intégrateur PBM6 et la triple optimisation complète.
- **Deliverables**:
  - `scratch/forecast_solar_system_46bodies_1yr.py`: Pipeline complet de simulation et de collecte des trajectoires.
- **Process**:
  - Charger l'état 46 corps en unités canoniques ($AU$, $M_\odot$, jours).
  - Configurer `AtlasManager` pour générer l'arbre de Jacobi multi-échelles initial (couplage Soleil-Planètes et Planètes-Satellites).
  - Instancier PBM6 avec la triple optimisation (`apex: optimized`, `kepler: optimized`, `supernode: jit`).

### Step 2: Intégration prospective conormale sur 1 an (2026-09-29 -> 2027-09-29)
- **Context**: 365,0 jours de simulation sous pas macro calibré $\Delta t_{\text{macro}} = 0{,}5\text{ jour}$ ($12\text{ heures}$, soit exactement 730 pas macro). Ce pas alloue 176 pas par période orbitale de Mercure (phase $\Delta \theta \approx 2{,}0^\circ$) et laisse le découplage conormal de Stumpff-Sundman intégrer analytiquement les feuilles serrées (Phobos, Mimas, Io) complété par le découpage dyadique pour les perturbations de marée.
- **Objective**: Exécuter la propagation conormale sans aucune divergence ni perte de régularité sous temps déplié et super-nœuds JIT.
- **Deliverables**:
  - Positions et vitesses finales simulées $q_{\text{sim}}, v_{\text{sim}}$ des 46 corps à $t = 365{,}0$ jours (2027-09-29).
  - Télémétrie de conservation des invariants (quantité de mouvement, moment cinétique, énergie totale).
- **Process**:
  - Exécuter `nbody.evolve(t_span=365.0, dt=dt)`.
  - Extraire l'état physique final via le dé-processing canonique $C_6(h)$.

### Step 3: Acquisition des éphémérides officielles NASA JPL Horizons à $t = +1\text{ an}$ (2027-09-29)
- **Context**: L'API REST NASA JPL Horizons fournit les vecteurs d'état officiels au format JSON/VECTORS en coordonnées ICRF/J2000 barycentriques (`@0`).
- **Objective**: Télécharger les positions officielles des 46 corps célestes pour la date du 2027-09-29.
- **Deliverables**:
  - `models/solar_system/data/jpl_horizons_1yr_2027_09_29.json`: Cache des éphémérides officielles NASA.
- **Process**:
  - Interroger l'API JPL pour chacun des 46 identifiants officiels (`10` pour le Soleil, `399` pour la Terre, `301` pour la Lune, etc.).
  - Parser les blocs de vecteurs cartésiens et convertir en coordonnées AU et AU/jour.

### Step 4: Confrontation fine, analyse des écarts et rendu du tableau de bord 300 DPI
- **Context**: Comparer les positions simulées conormales $q_{\text{sim}}$ aux positions NASA $q_{\text{NASA}}$ pour chaque catégorie de corps (Soleil, planètes telluriques, planètes géantes, lunes, astéroïdes, comètes).
- **Objective**: Calculer l'erreur spatiale absolue $\|q_{\text{sim}} - q_{\text{NASA}}\|$ en AU et en kilomètres, l'erreur relative, analyser les contributions physiques (perturbations relativistes post-newtoniennes, marées), et générer une planche de diagnostic scientifique 300 DPI.
- **Deliverables**:
  - `brain/solar_system_46bodies_1yr_pbm6_vs_nasa_300dpi.png`: Planche de diagnostic 300 DPI.
  - `models/solar_system/data/solar_system_1yr_forecast_metrics.json`: Données numériques complètes.
  - Synthèse comparative détaillée dans le walkthrough.

## Acceptance Criteria:
### Step 1:
- Le chargeur de données instancie un `ExtendedPhaseState` valide à 46 corps avec masses solaires exactes.
### Step 2:
- L'intégration de 365 jours se déroule sans interruption, sans NaN et avec conservation de l'énergie et des invariants de Lie-Bargmann.
### Step 3:
- Les 46 corps disposent de leurs éphémérides officielles NASA JPL Horizons téléchargées et validées.
### Step 4:
- Les écarts orbitaux Terre, Mars, Jupiter, etc. sont rigoureusement quantifiés (attendus au millième d'AU ou mieux pour un modèle purement gravitationnel Newtonien N-corps non relativiste sans forces non-gravitationnelles de dégazage cométaire).
- La planche 300 DPI est générée et disponible dans les artéfacts.

## Observations:
> [!IMPORTANT]
> - Dans un modèle gravitationnel purement newtonien ponctuel non adouci ($\varepsilon = 0$), les effets relativistes généraux de Schwarzschild (avance du périhélie de Mercure de 43'' d'arc par siècle, $\sim 10^{-7}\text{ AU}$ sur 1 an), la forme quadripolaire $J_2$ des planètes et le dégazage non-gravitationnel des comètes actives introduisent un écart physique naturel prévisible par rapport aux éphémérides numériques intégrées DE440 de la NASA.
> - Pour les planètes majeures (Terre, Vénus, Mars, Jupiter, Saturne), la concordance newtonienne sur 1 an est attendue avec une précision de l'ordre de $10^{-4}$ à $10^{-5}\text{ AU}$ ($10^3$ à $10^4\text{ km}$).
