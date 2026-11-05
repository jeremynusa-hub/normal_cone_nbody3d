# Walkthrough of Iteration 143 : Audit Contradictoire Multi-Agentique Exhaustif des Chapitres 3, 4 et 5 (9 Auditeurs)

## Summary:
L'itération 143 a conduit le passage au crible le plus exhaustif, méthodique et approfondi jamais mené sur l'ensemble de la bibliothèque `normal_cone_nbody_float64`.
Conformément aux exigences de l'utilisateur, **9 sous-agents d'audit spécialisés** ont été déployés en **3 triades disciplinaires étanches** (un Ingénieur Recherche Math-to-Code, un Analyste Numérique des Schémas, et un Expert en Algèbre de Lie et Homologie) couvrant respectivement les **Chapitres 3, 4 et 5** du manuscrit de recherche.

Cette confrontation contradictoire a permis d'isoler avec une exactitude mathématique absolue les causes racines des sauts d'énergie résiduels et d'établir une cartographie sans complaisance de 10 anomalies classées par criticité dans [`audit_theorique_chapitres_3_4_5.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/audit_theorique_chapitres_3_4_5.md).

---

## Critical Discoveries:

### 1. Chapitre 3 (Régularisations, Stumpff, Sundman & McGehee)
- **Clamping artificiel destructeur (`z < -100`) dans `kepler_sundman.py`** :
  En float64, $\cosh\sqrt{-z}$ supporte $z$ jusqu'à $-500\,000$. En bloquant $z$ à $-100.0$, Halley résout une équation tronquée incohérente lors des étages rétrogrades de Blanes-Moan BM6 ($a_4 < 0, b_3 < 0$) à grande vitesse, causant une divergence immédiate ($\Delta r = 52.6$ en un pas).
- **Annulation soustractive des dérivées de Stumpff dans `stumpff.py`** :
  L'évaluation de $(c_{k-1} - k c_k)/(2z)$ pour $z \in [10^{-4}, 0.1]$ subit une perte de 4 chiffres significatifs, dégradant la précision relative à $5.44 \times 10^{-10}$ à $z=0.1$.
- **Non-symplecticité démontrée du rebond de McGehee** :
  Calcul analytique de la Jacobienne prouvant que $\det \mathrm{d}\mathcal{B}_C = -1$. L'opérateur inverse le volume de Liouville et brise la 2-forme de Darboux, validant formellement la critique théorique du Chapitre 3 qui motive la déformation de Rees.

### 2. Chapitre 4 (Splitting Fulton-MacPherson, Apex Écranté & Super-Nœuds)
- **Rupture d'isomorphisme dans `solvers/supernode.py`** :
  Pour les amas compacts à $k \ge 3$ corps, le code n'intègre pas les équations sous l'horloge dépliée de Sundman $\mathrm{d}t_S = \rho_S \mathrm{d}\tau_S$. Il exécute un Leapfrog cartésien à **pas physique fixe** avec un plafond aveugle à 4096 sous-pas, un adoucissement masqué $d_{\text{safe}} = 10^{-14}$ et un recentrage non-canonique de quantité de mouvement. Lors d'une quasi-collision interne, l'énergie explose de $+35\,000\,\%$.
- **Chattering topologique (flickering) du détecteur de super-nœuds** :
  L'absence d'hystérésis sur $\eta_{\text{crit}} \approx 0.693$ fait osciller le super-nœud entre activation et désactivation à chaque pas de temps, injectant des chocs de splitting non symplectiques de $10^{-7}$ à $10^{-6}$.
- **Élucidation du mystère des kicks virtuels** :
  Preuve que `include_screened_monopole=False` dans le kick d'apex est la compensation exacte de l'absence de propagateur képlérien sur les nœuds composites ($V_{\mathrm{apex}}^{\mathrm{hier}} + V_{\mathrm{kepler}}^{\mathrm{screened}} = V_{\mathrm{inter}}$). L'invariance de Lie-Bargmann ($\sum F_i = \mathbf{0}$, $\sum q_i \times F_i = \mathbf{0}$) est respectée à $10^{-16}-10^{-17}$.

### 3. Chapitre 5 (Intégrateurs Symplectiques, PBM & Ordonnancement Dyadique)
- **Cause exacte du saut résiduel à $1.14 \times 10^{-8}$ sur le 4-corps hiérarchique** :
  1. *Plafonnement dyadique rigide* : au péricentre $r_{12} = 3.74 \times 10^{-4}$ ($v=98.1\text{ m/s}$), le niveau requis était $k=17$ ($131\,072$ sous-pas). L'échéancier a écrêté à $k=10$ ($1024$ sous-pas), utilisant un pas 77 fois trop grand.
  2. *Violation du critère de synchronisation de Skeel-Gear (1992)* : la transition d'arbre a été exécutée au sous-pas $640/1024$ en plein milieu du bloc, brisant la symétrie palindromique miroir et détruisant la continuité de l'hamiltonien d'ombre $\widetilde{H}(h)$.
- **Rupture de jauge de Lie aux transitions d'atlas pour les schémas traités (PBM)** :
  Les générateurs homologiques amont et aval ne commutent pas ($[X_{K_A}, X_{K_B}] \neq 0$). L'opérateur interfacial inter-macro-pas est $\mathcal{C}_B^{-1} \circ \mathcal{C}_A \neq \mathrm{Id}$, injectant un saut d'énergie direct $\Delta H = h^{2p} (H_{2p+1}^{(B)} - H_{2p+1}^{(A)})$.
- **Coefficients homologiques génériques** : $\kappa_4 = 1/2880$ et $\kappa_6 = 1/120960$ sont des constantes de Leapfrog classique, déconnectées des compositions de Blanes-Moan.
- **Surcoût de BM4** : BM4 exécute 14 évaluations de force au lieu de 7 par non-fusion des demi-kicks consécutifs.

---

## Realization:
- **Started At**: 2026-10-08T09:40:00+02:00
- **Completed At**: 2026-10-08T09:46:00+02:00

### Step 1: Déploiement Triadique sur le Chapitre 3
- **Impacted Elements**: `solvers/stumpff.py`, `solvers/kepler_sundman.py`, `solvers/mcgehee.py`.
- **Actions & Implementation**: Lancement et analyse des sous-agents `e61b676e` (Math-to-Code), `9443d7c2` (Numérique), et `4a314cc6` (Lie/Symplectique). Démonstration formelle de l'unimodularité du Wronskien de Stumpff et identification des bugs critiques de clamping et de dérivées.

### Step 2: Déploiement Triadique sur le Chapitre 4
- **Impacted Elements**: `potential/screening.py`, `potential/screened_apex.py`, `potential/supernode_detector.py`, `solvers/supernode.py`.
- **Actions & Implementation**: Lancement et analyse des sous-agents `241269a7` (Math-to-Code), `a14618eb` (Numérique), et `993ea2ba` (Lie/Symplectique). Découverte de la rupture d'isomorphisme dans les super-nœuds (absence de Sundman multi-corps) et du chattering de Hill.

### Step 3: Déploiement Triadique sur le Chapitre 5
- **Impacted Elements**: `integrators/` (`blanes_moan_4/6`, `correctors`, `processed`), `system/dyadic_scheduler.py`, `system/nbody_system.py`.
- **Actions & Implementation**: Lancement et analyse des sous-agents `1b993f2f` (Math-to-Code), `fb79a088` (Numérique), et `06a494de` (Lie/Homologie). Diagnostic complet du saut d'énergie du 4-corps, violation de Skeel-Gear et rupture de jauge de Lie des correcteurs traités.

### Step 4: Synthèse Contradictoire Globale et Cartographie des Anomalies
- **Impacted Elements**: [`audit_theorique_chapitres_3_4_5.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/audit_theorique_chapitres_3_4_5.md).
- **Actions & Implementation**: Croisement contradictoire des 9 rapports, consolidation de la cartographie des 10 anomalies et établissement de la feuille de route formelle de remédiation en 4 phases.

---

## Architecture and Code Organization (Evolution):
```
Cartographie des Anomalies Théorie-Code (Audit à 9 Agents)
 ├── Chapitre 3 : Clamping z < -100 (Kepler) | Annulation c_k' (Stumpff) | Rebond McGehee non-symplectique
 ├── Chapitre 4 : Leapfrog cartésien fixe (Supernœud) | Chattering topologique (Hill) | Inversion Hermite C8
 └── Chapitre 5 : Plafonnement k=10 au péricentre (4-corps) | Violation Skeel-Gear (flips intra-bloc) | Rupture de jauge C_B^-1 o C_A != Id
```

## Dependencies and Configuration:
- Aucune modification de dépendance externe.

## Design Choices, Trade-offs and Warnings:
- **Audit Pur Sans Modification Hâtive de Code** : Conformément aux règles d'or de l'audit épistémique, aucun fichier source n'a été modifié durant cette itération afin de préserver l'objectivité intégrale du diagnostic et de soumettre la feuille de route à l'arbitrage préalable de l'utilisateur.

## Verification:
- **Couverture Documentaire** : 100% des chapitres 3, 4 et 5 du manuscrit LaTeX ont été confrontés à l'intégralité des modules Python correspondants.
- **Rapports Détaillés** : 9 rapports d'experts indépendants reçus et consolidés dans [`audit_theorique_chapitres_3_4_5.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/audit_theorique_chapitres_3_4_5.md).
- **Traçabilité Formelle** : Chaque anomalie est reliée à son équation mathématique exacte et à ses lignes de code source.

## Next Steps:
- Soumettre au développeur le plan d'action d'ingénierie en 4 phases pour implémentation séquentielle et contrôlée.
