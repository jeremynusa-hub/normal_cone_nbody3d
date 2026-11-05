# Implementation Plan of Iteration 178 : Prévision sur 1 an du Système Solaire Étendu (46 corps) sous pas macroscopique de 3h et PBM6
Created at: 2026-10-09T14:17:00+02:00

## Context:
1. **Demande explicite** : Lancer la prévision complète sur 1 an (365.0 jours, 2 920 macro-pas) du système solaire étendu à 46 corps avec l'intégrateur traité de Blanes-Moan d'ordre 6 (PBM6), un pas macroscopique de 3 heures ($\Delta t = 0.125\text{ j}$), l'ordonnanceur dyadique dynamique actif (résolvant l'orbite de Phobos à $k=6$, soit 64 sous-pas de 2.81 min par macro-pas, totalisant 186 880 sous-pas de PBM6), et confronter les positions finales aux données officielles de la NASA (JPL Horizons au 29 septembre 2027).
2. **Valeur stratégique implicite** : Valider la robustesse, la conservation symplectique et la haute précision de la géométrie conormale dérégularisée sous ordonnancement dyadique continu et super-nœuds synchrones JIT sur une révolution planétaire complète, sans gel computationnel, en capturant fidèlement les satellites ultra-rapides et les petits corps.
3. **Critères d'alignement** : Conservation de l'énergie hamiltonienne $\Delta H / |H_0| \le 10^{-9}$, conservation du moment cinétique $\Delta L / |L_0| \le 10^{-9}$, dérive du moment linéaire $\le 10^{-16}$, écarts orbitaux planétaires kilométriques conformes à l'état de l'art, et génération d'un tableau comparatif et d'une planche graphique de publication 300 DPI.

## Guidelines:
- **Méthode** : Intégration en mode ombre (`in_shadow_mode=True`) pour télescopage des correcteurs symplectiques, pas macro $\Delta t = 0.125\text{ j}$ (3h), ordonnanceur dyadique dynamique sous le nouveau plancher numérique corrigé ($10^{-50}$ pour $GM$, $10^{-100}$ pour $(GM)^2$), triple optimisation active (Screened Apex optimisé C8 d'ordre 8, Kepler-Sundman vectorisé optimisé, Super-nœuds JIT LLVM).
- **Lignes rouges** : Zéro-adoucissement strict ($\varepsilon = 0, d_0 = 0$), interdiction formelle de désactiver les super-nœuds synchrones (`auto_supernodes = True` impératif), pas de simplification empirique non conforme à la variété conormale $Q[N]$.
- **Architecture** : Script de production dédié `scratch/forecast_solar_system_46bodies_1yr_3h.py` avec points d'étape de progression (tous les 250 pas / ~31 jours), sauvegarde de l'état simulé dans `models/solar_system/data/solar_system_simulated_1yr_3h_pbm6.json`, confrontation avec `models/solar_system/data/jpl_horizons_1yr_2027_09_29.json` déjà disponible, et génération d'une planche 300 DPI dédiée.

## Plan:
### Step 1: Configuration et écriture du script de forecast de production sur 1 an sous pas de 3h
- **Contexte** : Données initiales à 46 corps disponibles au 29-09-2026, éphémérides JPL Horizons de référence au 29-09-2027 déjà acquises.
- **Objectif** : Écrire le script complet `scratch/forecast_solar_system_46bodies_1yr_3h.py` orchestrant l'intégration sur 2 920 pas de $\Delta t = 0.125\text{ j}$ avec journalisation des temps d'exécution, conservation des invariants et sauvegarde d'états intermédiaires.
- **Livrables** :
  - `ScriptProductionForecast`: `scratch/forecast_solar_system_46bodies_1yr_3h.py`.
- **Processus** : Initialisation de l'état `ExtendedPhaseState`, arbre de Jacobi canonique, intégrateur PBM6 avec télescopage d'ombre, `DyadicScheduler` avec $\Delta t = 0.125\text{ j}$ et niveau dynamique actif, boucle de 2 920 pas avec télémétrie par tranche de 250 pas.

### Step 2: Exécution de l'intégration symplectique sur 1 an (365.0 jours)
- **Contexte** : Script de production prêt, environnement Python calibré.
- **Objectif** : Exécuter l'intégration complète des 2 920 pas en arrière-plan, surveiller sa progression et enregistrer l'état final physique.
- **Livrables** :
  - `SimulatedExtendedState1Yr`: `models/solar_system/data/solar_system_simulated_1yr_3h_pbm6.json`.
- **Processus** : Lancement asynchrone via `run_command`, suivi automatique de progression jusqu'à terminaison, extraction de l'état physique via application du correcteur symplectique final $C(h)$.

### Step 3: Confrontation analytique contre les éphémérides officielles NASA JPL Horizons
- **Contexte** : Données simulées au 29-09-2027 et vérité terrain JPL Horizons barycentrique `@0`.
- **Objectif** : Calculer les résidus de position $\|\mathbf{q}_{\text{sim}} - \mathbf{q}_{\text{NASA}}\|$ (en km et en AU) et de vitesse $\|\mathbf{v}_{\text{sim}} - \mathbf{v}_{\text{NASA}}\|$ pour les 46 corps célestes, et comparer aux résultats obtenus précédemment avec le pas de 12h.
- **Livrables** :
  - `ComparisonReport`: Rapport d'erreurs et statistiques agrégées (médiane, percentiles 90%, planètes, lunes, astéroïdes).
- **Processus** : Chargement des deux fichiers JSON, appariement par nom de corps céleste, calcul vectoriel des normes d'écart, tri hiérarchique et bilan analytique.

### Step 4: Génération de la planche de rendu 300 DPI et rapport de synthèse
- **Contexte** : Trajectoires et résidus d'erreurs calculés sur 1 an.
- **Objectif** : Générer une planche graphique haute résolution 300 DPI à 4 panneaux illustrant les trajectoires 3D, les profils de conservation des invariants (énergie, moments), l'histogramme des niveaux dyadiques et la distribution des erreurs de position par corps.
- **Livrables** :
  - `PublicationDashboard300DPI`: `C:\Users\jerem\.gemini\antigravity\brain\59313b1a-0dde-45d6-a476-1bcb513cf82e\solar_system_46bodies_1yr_3h_pbm6_vs_nasa_300dpi.png`.
  - `WalkthroughDocument`: `walkthrough.md`.
- **Processus** : Rendu matplotlib haute fidélité sous `Agg`, sauvegarde en PNG 300 DPI, mise à jour des artifacts `task.md` et `walkthrough.md`.

## Acceptance Criteria:
### Step 1:
- Le script `scratch/forecast_solar_system_46bodies_1yr_3h.py` est syntaxiquement valide, importe tous les modules optimisés et configure $\Delta t_{\text{macro}} = 0.125\text{ j}$ avec `default_level=0, max_level=16`.
### Step 2:
- Les 2 920 macro-pas sont exécutés avec succès sans NaN, sans exception et sans décrochage de phase.
- Conservation de l'énergie hamiltonienne $\Delta H / |H_0| \le 1.0 \times 10^{-9}$ sur l'année entière.
- Conservation du moment cinétique relatif $\Delta L / |L_0| \le 1.0 \times 10^{-9}$.
- Dérive de quantité de mouvement totale $\le 1.0 \times 10^{-16}$.
### Step 3:
- Les 46 corps sont confrontés un à un aux éphémérides officielles JPL Horizons au 29 septembre 2027.
- L'erreur médiane de position reste confinée dans l'ordre de grandeur kilométrique/sub-millimétrique pour les planètes majeures.
### Step 4:
- La planche 300 DPI est générée sans collision typographique, lisible et enregistrée dans le répertoire des artifacts.
- Le fichier `walkthrough.md` documente de manière exhaustive l'ensemble des résultats et métriques.

## Observations:
> [!NOTE]
> L'intégration sur 1 an complet à $\Delta t_{\text{macro}} = 3\text{ h}$ ($0.125\text{ j}$) représente 2 920 macro-pas et 186 880 sous-pas de PBM6. Le temps d'exécution estimé est d'environ 60 à 70 minutes sur le processeur hôte.
> [!IMPORTANT]
> Conformément à la règle de zéro-adoucissement strict et à la géométrie conormale, les super-nœuds synchrones restent intégralement actifs (`auto_supernodes = True`), garantissant la désingularisation sans artefact.
