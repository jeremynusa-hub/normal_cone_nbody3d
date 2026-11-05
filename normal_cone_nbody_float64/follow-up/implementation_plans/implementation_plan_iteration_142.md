# Implementation Plan of Iteration 142 : Calibrage Typographique Analytique des Marges de Titres
Created at: 2026-10-08T09:25:00+02:00

## Context:
1. **Explicit Request**:
   - L'utilisateur recadre fermement : "Là c'est beaucoup trop, choisis une hauteur adapté a la police et fait de vrais calcul pas du doigt mouillé".
2. **Implicit Strategic Value**:
   - Établir une modélisation mathématique et typographique rigoureuse de la mise en page, fondée sur les canons de l'édition scientifique (Bringhurst, Chicago Manual of Style) reliant la taille de police du titre ($F = 18\text{ pt}$) aux marges d'aération ($M_{\text{top}} \approx 1.25\,F$, $M_{\text{gap}} \approx 1.8 \text{ à } 2.0\,F$), au lieu d'ajustements pifométriques ou disproportionnés.
3. **Alignment Criteria**:
   - Pour un titre en corps 18 pt ($0.25\text{ in} = 75\text{ px}$ à 300 DPI), le blanc résiduel sous le titre doit s'établir rigoureusement entre $1.8\times$ et $2.0\times$ la taille de la police, soit $34\text{ pt}$ à $36\text{ pt}$ ($140\text{ px}$ à $150\text{ px}$ réels à 300 DPI), avec une marge supérieure de $1.25\times F \approx 22.5\text{ pt}$ ($94\text{ px}$).
   - Hauteurs de canevas redimensionnées à leur juste échelle ($22.0\text{ in}$ pour 10 panneaux, $17.0\text{ in}$ pour multi-intégrateurs) avec $H_{\text{banner}} \approx 1.0\text{ in}$ ($72\text{ pt} = 4.0\,F$).
   - Régénération des 4 planches 300 DPI et certification formelle par analyse PIL/NumPy.
   - Validation 100% de la suite de tests unitaires (132/132 tests).

## Guidelines:
- **Method**: Dérivation analytique exacte des coordonnées normalisées de figure Matplotlib :
  * $H_{\text{header}} = M_{\text{top}} + F + M_{\text{gap}} = (1.25 + 1.00 + 1.85) \times \frac{18}{72} = 1.025\text{ in}$.
  * Bounding box subplots : $\mathrm{rect}_{\text{top}} = 1.0 - \frac{H_{\text{header}}}{H_{\text{fig}}}$.
  * Position verticale du titre : $y_{\text{center}} = 1.0 - \frac{M_{\text{top}} + 0.5\,F}{H_{\text{fig}}}$.
- **Red Lines**:
  * Zéro arbitraire : chaque pixel de marge doit découler de la taille de police $F = 18\text{ pt}$.
  * Ne pas déformer ni comprimer les sous-panneaux scientifiques.
- **Architecture**:
  * `diagnostics/visualizer.py` : intégration des formules typographiques analytiques dans `render_10panel_dashboard` et `render_multi_integrator_invariants_dashboard`.
  * `scratch/run_showcase_simulations.py` : régénération et synchronisation vers le brain.

## Plan:
### Step 1: Modélisation Typographique et Intégration dans visualizer.py
- **Context**: Le titre en corps 18 pt nécessite une bande d'en-tête de $1.0\text{ in}$ ($72\text{ pt}$), répartie en $22.5\text{ pt}$ de marge haute ($1.25\,F$), $18\text{ pt}$ de texte ($1.0\,F$), et $34\text{ pt}$ à $36\text{ pt}$ de dégagement bas ($1.9\text{ à } 2.0\,F$).
- **Objective**: Paramétrer `PublicationVisualizer` avec les valeurs analytiques exactes :
  * 10 panneaux : `figsize=(16, 22.0)`, `rect=[0, 0, 1, 0.955]`, `y=0.980`, `va='center'`.
  * Multi-intégrateurs : `figsize=(16, 17.0)`, `rect=[0, 0, 1, 0.940]`, `y=0.974`, `va='center'`.
- **Deliverables**:
  * `AnalyticalTypographicVisualizer`: Code mis à jour dans `normal_cone_nbody_float64/diagnostics/visualizer.py`.
- **Process**: Édition ciblée via `replace_file_content`.

### Step 2: Régénération des 4 Planches 300 DPI et Validation Pixel
- **Context**: Les planches doivent être recalculées et certifiées matriciellement.
- **Objective**: Exécuter `scratch/run_showcase_simulations.py`, mesurer avec exactitude les intervalles de pixels de blanc sur les images 300 DPI exportées, et vérifier que le ratio d'espacement est exactement $\approx 1.85 \text{ à } 2.0 \times \text{font size}$ ($140-150\text{ px}$). Copier dans le brain.
- **Deliverables**:
  * `HarmoniousDashboards`: Les 4 fichiers PNG 300 DPI dans le brain avec espacement typographique certifié.
- **Process**: Exécution du script et synchronisation.

### Step 3: Tests de Non-Régression et Documentation
- **Context**: Assurer la solidité globale de l'écosystème de tests et consigner les équations.
- **Objective**: Lancer `pytest tests/unit/`, mettre à jour `walkthrough.md`, `task.md`, et consigner la note atomique dans `notes.md`.
- **Deliverables**:
  * `CertifiedVisualizerSuite`: 132/132 tests réussis.
- **Process**: Exécution des tests et consolidation documentaire.

## Acceptance Criteria:
### Step 1:
- Paramètres `rect` et `y` dans `visualizer.py` rigoureusement déduits de la formule analytique $F = 18\text{ pt}$.
### Step 2:
- Mesure matricielle démontrant formellement un blanc résiduel compris entre $140\text{ px}$ et $155\text{ px}$ à 300 DPI ($1.8\times$ à $2.0\times$ le corps 18 pt), sans aucun vide disproportionné.
### Step 3:
- 132/132 tests unitaires validés sous `pytest tests/unit/`.

## Observations:
> [!IMPORTANT]
> Pour un corps de $18\text{ pt}$ à 300 DPI ($75\text{ px}$), la règle d'art typographique impose une respiration de $1.5\times$ à $2.0\times$ la hauteur d'œil, soit $112\text{ px}$ à $150\text{ px}$. Les 428 pixels précédents représentaient près de $6\times$ la hauteur de police, créant une impression de vide béant désormais résolue par ce calcul analytique.
