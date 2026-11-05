# Implementation Plan of Iteration 126 : Régénération des Planches Scientifiques 300 DPI et Synchronisation des Données de la WebApp
Created at: 2026-10-07T19:40:00+02:00

## Context:
1. **Explicit Request**: « relance la création des planche et des données de la webapp du coup »
2. **Implicit Strategic Value**: Actualiser tous les visuels scientifiques haute résolution (300 DPI) et les données interactives de la WebApp pour refléter l'élimination définitive du blowup 4-corps (erreur relative passée de $+78\,000\,\%$ à $4.73 \times 10^{-13}$), la cinématique galiléenne exacte sans dérive double, et la certification formelle 100 % [PASS] sur le benchmark chaotique de Burrau et les transitions d'atlas.
3. **Alignment Criteria**: Génération complète sans avertissement ni instabilité des 3 dashboards scientifiques 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`), mise à jour du fichier `showcase_data.json`, synchronisation de l'application web [`interactive_conormal_showcase.html`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/interactive_conormal_showcase.html) et copie de tous les artefacts dans le répertoire brain.

## Guidelines:
- **Sanctuarisation Absolue** : Ne jamais modifier aucun fichier dans `redaction/unsoftened_conormal_nbody3d/`.
- **Zéro Adoucissement Strict** : Maintenir rigoureusement $\varepsilon = 0, d_0 = 0$ sur l'ensemble des scénarios.
- **Intégrité Numérique** : Utiliser les solveurs et intégrateurs conormaux corrigés (`blanes_moan_6`, `blanes_moan_4`, `yoshida`, `strang`, `processed_bm6`) avec l'adaptation de l'arbre de Jacobi aux sous-pas dyadiques.
- **Reproductibilité Scientifique** : Garantir la fidélité des métriques célestes affichées sur les 10 panneaux de chaque planche.

## Plan:
### Step 1: Exécution de la Suite de Simulation et Régénération des Planches Scientifiques (`scratch/run_showcase_simulations.py`)
- **Matériel de départ** : [`scratch/run_showcase_simulations.py`](file:///C:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scratch/run_showcase_simulations.py), bibliothèque [`normal_cone_nbody_float64/`](file:///C:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/normal_cone_nbody_float64/).
- **Objectif** : Exécuter les 4 scénarios canoniques (Figure-8 avec comparatif multi-intégrateurs et Pareto, Burrau pythagoricien à $t=10$, choc radial avec rebond de McGehee, et système 4-corps hiérarchique avec 5 basculements d'atlas).
- **Livrables** : `dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`, et `showcase_data.json` générés dans `scratch/`.
- **Processus** : Lancement du script via `run_command` et suivi de l'exécution asynchrone.

### Step 2: Synchronisation des Données JSON et Actualisation de la WebApp Interactive (`scratch/sync_showcase_html.py`, `interactive_conormal_showcase.html`)
- **Matériel de départ** : [`scratch/sync_showcase_html.py`](file:///C:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/scratch/sync_showcase_html.py), `scratch/showcase_data.json`, [`interactive_conormal_showcase.html`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/interactive_conormal_showcase.html).
- **Objectif** : Injecter le nouveau payload JSON consolidé dans la constante `SHOWCASE_DATA` de la WebApp interactive.
- **Livrable** : `interactive_conormal_showcase.html` synchronisé avec les trajectoires résolues sans explosion numérique.
- **Processus** : Exécution du script de synchronisation.

### Step 3: Copie et Déploiement des Artefacts Visuels dans le Répertoire Brain (`brain/`)
- **Matériel de départ** : Fichiers images et JSON générés dans `scratch/`.
- **Objectif** : Déployer les nouvelles planches 300 DPI et le fichier `showcase_data.json` directement dans le répertoire des artefacts brain (`C:\Users\jerem\.gemini\antigravity\brain\59313b1a-0dde-45d6-a476-1bcb513cf82e\`).
- **Livrables** : Artefacts finaux prêts pour consultation et affichage par l'utilisateur.
- **Processus** : Copie automatisée via script Python.

### Step 4: Contrôle Qualité et Certification Visuelle des Indicateurs Célestes
- **Matériel de départ** : Fichiers PNG générés et métriques des rapports de certification.
- **Objectif** : Valider l'intégrité des 10 panneaux graphiques (trajectoires 3D, spectres d'énergie, dérive séculaire, invariants de Lie-Bargmann, radar conormal, transitions d'atlas, et tableau récapitulatif 100 % [PASS]).
- **Livrable** : Rapport de vérification et validation complète de l'Itération 126.
- **Processus** : Inspection des tailles de fichiers, des métriques et compte-rendu exhaustif.

## Acceptance Criteria:
- **Critère 1** : Les 3 planches scientifiques 300 DPI sont générées avec succès (taille > 1 Mo chacune) et exemptes de tout crash ou explosion d'échelle.
- **Critère 2** : Le panneau 10 de chaque planche certifie l'ensemble des KPIs au statut vert `[PASS]`.
- **Critère 3** : Le fichier `showcase_data.json` contient l'ensemble des trajectoires recalculées avec l'arbre de Jacobi adaptatif aux sous-pas dyadiques.
- **Critère 4** : L'application web interactive `interactive_conormal_showcase.html` est synchronisée et fonctionnelle.
- **Critère 5** : L'invariance stricte sans adoucissement ($\varepsilon = 0, d_0 = 0$) est préservée sur tous les graphiques.

## Observations:
> [!NOTE]
> La simulation Figure-8 inclut un échantillonnage Pareto sur 5 pas de temps distincts ($dt \in [0.0075, 0.03]$) pour les 4 intégrateurs, démontrant visuellement la supériorité asymptotique de BM6.
> [!IMPORTANT]
> Les planches sont rendues à 300 DPI sans interface graphique interactive (`Agg` backend), assurant un rendu haute précision pour publication scientifique.
