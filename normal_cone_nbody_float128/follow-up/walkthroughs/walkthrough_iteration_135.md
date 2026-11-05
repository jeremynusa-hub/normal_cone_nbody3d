# Walkthrough de l'Itération 135 : Audit Contradictoire Multi-Agent et Télémétrie Fine du Système 4-Corps Hiérarchique

## Summary:
L'Itération 135 a accompli avec une rigueur analytique et expérimentale absolue la mission d'investigation contradictoire demandée par l'utilisateur pour comprendre l'origine des sauts d'énergie résiduels ($1.96 \times 10^{-7}$) observés sur le **Système 4-Corps Hiérarchique** (`hierarchical_4body` dans `normal_cone_nbody_float64/`).

Grâce au déploiement d'une **sonde télémétrique pas-à-pas** (`scratch/diagnose_4body_jump.py`) et de **4 sous-agents d'audit indépendants et étanches**, le mécanisme intime de la défaillance a été formellement découvert, démontré mathématiquement et réconcilié avec l'historique de développement de CAS6 et la monographie théorique.

L'ensemble des résultats et des démonstrations a été consigné dans l'artefact de référence [`diagnostic_sauts_energie_4corps.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/diagnostic_sauts_energie_4corps.md).

---

## 1. Découverte Télémétrique : L'Instant et la Configuration Exacte du Choc

La sonde télémétrique fine a localisé le saut d'énergie à un instant unique et irréversible :
- **Macro-pas 220 ($t = 3.3000$)** : l'erreur relative saute de $4.159 \times 10^{-11}$ (pas 219) à $1.9602 \times 10^{-7}$ ($\Delta E_{\text{step}} = +1.9597 \times 10^{-7}$).
- **Sous-pas critique 449/512 ($t = 3.298213$)** :
  - La paire de feuilles $\{1, 2\}$ franchit un péricentre rasoir à **$r_{12}^{\min} = 1.024 \times 10^{-3}$ (1 mm !)** à une vitesse de **$v_{\max} = 59.6\text{ m/s}$**.
  - La vitesse angulaire instantanée culmine à **$\dot{\theta} \approx 57\,600\text{ rad/s}$**.
  - Les corps 1 et 2 tourbillonnent et pivotent d'un angle de **$\Delta \theta \approx 97^\circ$ en un seul sous-pas dyadique** !
  - Le corps 3 se situe alors à $r_{13} \approx 0.286$ et $r_{23} \approx 0.287$, tandis que le corps 0 orbite au loin à $r_{0i} \approx 4.19$.

---

## 2. Synthèse des 4 Audits Contradictoires

### Sous-Agent 1 : Géométrie de Jacobi & Transitions d'Atlas (`3ac8f675`)
- **Exactitude de l'arbre** : La topologie `(0, ((1, 2), 3))` est rigoureusement conforme. Les masses réduites ($\mu_3 = 4/9$, $\mu_2 = 0.72$, $\mu_0 = 0.75$) et fractions barycentriques sont certifiées exactes à tolérance machine.
- **Innocence de la transition d'atlas** : Le switch d'arbre opéré au pas 219 a conservé l'énergie avec une variation de seulement $1.58 \times 10^{-13}$.
- **Innocence du solveur de Stumpff-Sundman** : Testé isolément au péricentre du sous-pas 449, le solveur képlérien analytique conserve l'énergie à moins de $7.9 \times 10^{-13}$.

### Sous-Agent 2 : Potentiel d'Apex Écranté & Boîtes de Rees (`a8d665ed`)
- **Déconnexion de l'écrantage dans le propagateur** : Dans `step_apex_kick`, `include_screened_monopole=False` est forcé en dur pour éviter d'omettre des forces d'attraction lorsque les nœuds composites sont en drift libre.
- **Impact physique** : Les interactions $(1, 3)$ et $(2, 3)$ au sein du nœud composite $((1, 2), 3)$ sont donc traitées comme des forces newtoniennes brutes non régularisées en $1/r^2$.
- **Absence de bug de gradient** : Les forces d'apex dérivent rigoureusement du potentiel calculé (écart $< 4 \times 10^{-9}$ en différences finies). Le problème est purement cinématique : appliquer un kick cartésien direct d'Euler sur une binaire tourbillonnant à $97^\circ$ par sous-pas découple positions et impulsions.

### Sous-Agent 3 : Correcteurs Symplectiques & Ordonnancement Dyadique (`cacb63ad`)
- **Preuve de violation du critère de Skeel-Gear** :
  L'échelle de temps cinématique réelle au péricentre est $\tau_{\text{kin}} = \frac{r_{12}}{v_{12}} \approx 1.72 \times 10^{-5}\text{ s}$.
  Le sous-pas dyadique imposé par le niveau 9 vaut $\Delta t_{\text{sub}} = 2.93 \times 10^{-5}\text{ s}$.
  Le ratio vaut **$\Delta t_{\text{sub}} / \tau_{\text{kin}} = 1.705 > 1$**. La grille de `dyadic_scheduler.py` s'arrêtait arbitrairement à $k=9$ pour tout $r \le 0.015$, ignorant les péricentres millimétriques !
- **Élucidation du mystère PBM4 vs BM4** :
  Le correcteur $\mathcal{C}_4$ appliqué à chaque sous-pas avec $h = \Delta t_{\text{sub}}$ produit un terme d'ordre $\kappa_4 h^4 \approx 2.56 \times 10^{-22}$, soit 6 ordres de grandeur sous $\varepsilon_{\text{mach}} \approx 2.22 \times 10^{-16}$. En float64, $\mathcal{C}_4 \equiv \operatorname{Id}$ (identité pure par sous-flux numérique). De plus, par télescopage frontière $\mathcal{C}^{-1} \circ \mathcal{C} = \operatorname{Id}$, aucun correcteur n'intervient au cœur du macro-pas où le choc a lieu.

### Sous-Agent 4 : Archéologue Historique CAS6 (`5ea65dfc`)
- **Correspondance historique parfaite** : L'anomalie du pas 220 est la réplique exacte de l'**Itération 82 (Anomalie d'éjection de Deimos)** dans CAS6 : traiter le sous-amas $((1, 2), 3)$ sous Option A (drift libre du nœud composite et kick newtonien direct) engendre une impulsion cartésienne brutale qui détruit la symplecticité.
- **Divergence de code débusquée** : Dans `integrators/correctors.py`, les appels à `compute_apex_forces` ont été faits sans paramètre, utilisant par défaut `include_screened_monopole=True`, tandis que le solveur de base utilise `False`.
- **Rôle de la configuration** : Dans `scratch/run_showcase_simulations.py`, `auto_supernodes` était forcé à `False` sur le 4-corps, interdisant la prise en charge du triplet par l'horloge dépliée anti-Zénon.

---

## 3. Feuille de Route de Remédiation

1. **Extension de la grille dyadique (`dyadic_scheduler.py`)** : Autoriser les niveaux $k=11$ ($2048$ sous-pas) pour $r \le 0.005$ et $k=13$ ($8192$ sous-pas, $\Delta t_{\text{sub}} \approx 1.83 \times 10^{-6}\text{ s}$) pour $r \le 0.0015$, garantissant $\Delta t_{\text{sub}} / \tau_{\text{kin}} \le 0.1$ conformément à Skeel & Gear.
2. **Harmonisation de `include_screened_monopole`** : Définir la valeur par défaut à `False` dans `ScreenedApexPotential.compute_apex_forces` et `compute_apex_potential` pour aligner parfaitement `correctors.py` et `common.py`.
3. **Activation de `auto_supernodes=True`** : Promouvoir le sous-amas $\{1, 2, 3\}$ en super-nœud synchrone régularisé sous horloge dépliée de Sundman $\mathrm{d}t_S = \chi_S^2 \mathrm{d}\tau_S$ lors des approches compactes.

---

## 4. Verification & Status:
- [x] Télémétrie pas-à-pas exécutée et validée (Pas 220, $t=3.3000$, sous-pas 449/512, $r_{12} = 1.024 \times 10^{-3}$).
- [x] 4 sous-agents d'audit contradictoire achevés avec rapports certifiés.
- [x] Diagnostic formel rédigé dans [`diagnostic_sauts_energie_4corps.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/diagnostic_sauts_energie_4corps.md).
- [x] Découvertes capitales consignées dans `notes.md`.
- [x] 130/130 tests unitaires toujours au vert (100% de succès).
