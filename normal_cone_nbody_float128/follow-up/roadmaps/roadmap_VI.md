# Strategic Roadmap : Remédiation Théorique Complète et Certification Formelle des Schémas Symplectiques (Roadmap VI)

## Strategic Overview:
- **Objectif Stratégique Général** : Éradiquer l'ensemble des 10 anomalies théoriques et numériques identifiées par l'audit contradictoire à 9 agents sur les Chapitres 3, 4 et 5 du manuscrit de recherche, afin d'éliminer définitivement les sauts résiduels d'énergie, de restaurer l'isomorphisme mathématique absolu entre la théorie et le code Python `normal_cone_nbody_float64`, et de certifier les performances ultimes sur l'ensemble des benchmarks astrophysiques canoniques.
- **Règles d'Or & Contraintes Inviolables** :
  1. *Sanctuarisation Documentaire* : Le répertoire du manuscrit de recherche [`redaction/unsoftened_conormal_nbody3d/`](file:///c:/Users/jerem/Antigravity%20Workspaces/antigravity-projects/mathsStuff/redaction/unsoftened_conormal_nbody3d/) est strictement en lecture seule.
  2. *Physique Non Adoucie Stricte* : Zéro adoucissement ni seuil de sécurité masqué ($\varepsilon = 0, d_0 = 0, d_{\mathrm{safe}} = 0$). Toute singularité doit être résolue par difféomorphisme analytique de Sundman ou éclatement conormal de Fulton--MacPherson.
  3. *Invariance Symplectique et de Lie-Bargmann* : Conservation stricte des 10 invariants de Noether ($\mathbf{P}, \mathbf{L}, \mathbf{K}, H$) à la précision machine float64 IEEE 754 ($< 10^{-14}$).
  4. *Règle de Synchronisation de Skeel-Gear (1992)* : Préservation stricte de la symétrie palindromique miroir des macro-blocs dyadiques et synchronisation des transitions d'atlas exclusivement aux frontières de macro-pas.
  5. *Certification Continue* : 100 % des tests unitaires existants et nouveaux doivent passer avec succès à chaque itération.
- **Document Maître de Référence** : [`audit_theorique_chapitres_3_4_5.md`](file:///C:/Users/jerem/.gemini/antigravity/brain/59313b1a-0dde-45d6-a476-1bcb513cf82e/audit_theorique_chapitres_3_4_5.md).

---

## Iteration 144 Strategy:
### Iteration 144: Phase 1 — Remédiation et Certification des Solveurs Universels du Chapitre 3 (Stumpff, Sundman-Halley, McGehee)
- **Scope & Objectives**:
  - `solvers/kepler_sundman.py` : Supprimer le clamping destructeur `z < -100` qui fait diverger l'itérateur cubique de Halley lors des pas rétrogrades ($\Delta t < 0$) de Blanes-Moan ; calibrer le garde-fou numérique à la limite réelle de $\cosh\sqrt{-z}$ en float64 ($z \approx -500\,000$).
  - `solvers/stumpff.py` : Résoudre l'annulation soustractive sur $(c_{k-1} - k c_k)/(2z)$ pour $z \in [10^{-4}, 0.1]$ en étendant le développement de Taylor de Horner jusqu'à $|z| < 0.2$ pour les dérivées universelles $c_k'(z)$, garantissant une précision relative machine ($< 10^{-15}$).
  - `solvers/mcgehee.py` : Purge intégrale du module et de ses tests, le modèle physique étant strictement celui de masses ponctuelles pures sans rebond matériel ni artefact de collision.
  - Enrichir et valider la suite de tests unitaires `tests/unit/solvers/`.
- **Dependencies & Continuity**: Établit le socle d'exactitude képlérienne indispensable pour toutes les compositions d'ordre élevé sous pas négatifs.
- **Gating Criteria**:
  - Wronskien universel $|f \dot{g} - \dot{f} g - 1| < 10^{-15}$ validé sur orbites hyperboliques rétrogrades avec $z$ jusqu'à $-100\,000$.
  - Dérivées de Stumpff $c_k'(z)$ vérifiées à $< 10^{-14}$ par rapport aux formes fermées analytiques sur tout le domaine $[-0.2, 0.2]$.
  - 100 % de succès sur `tests/unit/solvers/` et non-régression globale (129/129 tests).

---

## Iteration 145 Strategy:
### Iteration 145: Phase 2 — Remédiation des Super-Nœuds et Régularisation du Potentiel d'Apex du Chapitre 4 (Screening, Schmitt Hysteresis, Horloge Dépliée)
- **Scope & Objectives**:
  - `potential/supernode_detector.py` : Implémenter une bascule de Schmitt à hystérésis $[\eta_{\min}, \eta_{\max}] = [0.60, 0.75]$ autour du ratio de marée critique $\eta_{\mathrm{crit}} \approx 3^{-1/3} \approx 0.6933$ pour éliminer le broutement topologique (flickering).
  - `potential/screening.py` : Reformuler l'évaluation du polynôme d'écrantage Hermite $C^8$ en forme centrée $\sigma_8(s)$ autour de $s=0.5$, supprimant le conditionnement monomial défaillant ($\kappa \approx 1902$) et garantissant l'annulation exacte des dérivées d'ordre 1 à 8 aux bords $u=1$ et $u=2$.
  - `solvers/supernode.py` : Réécrire l'intégrateur interne sous temps déplié de Sundman $\mathrm{d}t_S = \chi_S^2\,\mathrm{d}\tau_S$ pour les sous-amas denses $k \ge 3$, éradiquer définitivement l'adoucissement masqué $d_{\mathrm{safe}} = 10^{-14}$ et la soustraction de quantité de mouvement non canonique.
  - Enrichir et valider les tests unitaires `tests/unit/potential/` et `tests/unit/solvers/test_supernode.py`.
- **Dependencies & Continuity**: S'appuie sur la régularité des solveurs de base (Itération 144) et sécurise l'évaluation des forces d'apex.
- **Gating Criteria**:
  - Absence totale d'oscillations d'activation/désactivation du super-nœud sur trajectoire oscillante autour de $\eta_{\mathrm{crit}}$.
  - Gradient analytique d'apex conforme aux différences finies à $< 10^{-11}$ pour l'écrantage $C^8$ sans anomalie en $s=0.5$.
  - Zéro adoucissement résiduel ($d_{\mathrm{safe}} \equiv 0$) dans `supernode.py`.

---

## Iteration 146 Strategy:
### Iteration 146: Phase 3 — Remédiation de la Synchronisation de Skeel-Gear et de l'Échéancier Dyadique du Chapitre 5
- **Scope & Objectives**:
  - `system/nbody_system.py` : Conformer rigoureusement l'évolution dyadique à la règle de Skeel & Gear (1992) en confinant les transitions topologiques d'arbres de Jacobi strictement aux frontières des macro-pas synchronisés, interdisant tout basculement intra-bloc pour garantir la symétrie palindromique du propagateur composite.
  - `system/dyadic_scheduler.py` : Remplacer l'écrêtage arbitraire $k_{\max}=10$ par une profondeur dyadique dynamique allant jusqu'à $k_{\max}=17$ ($131\,072$ sous-pas) en présence de péricentres ultra-rapides ($v \sim 100\text{ m/s}, r \sim 10^{-4}$), garantissant $h_{\mathrm{sub}} / \tau_{\mathrm{kin}} < 0.1$.
  - `integrators/blanes_moan_4.py` : Implémenter la fusion FSAL (First Same As Last) sur les demi-kicks consécutifs d'apex, réduisant le coût computationnel de 14 à 7 évaluations de force par pas macro.
  - Enrichir et valider les tests unitaires `tests/unit/system/` et `tests/unit/integrators/`.
- **Dependencies & Continuity**: Requiert la stabilité de l'apex et des super-nœuds (Itération 145).
- **Gating Criteria**:
  - Preuve formelle de symétrie palindromique miroir $\mathcal{S}(-h) \circ \mathcal{S}(h) = \operatorname{Id}$ sur l'ensemble du macro-pas dyadique.
  - Résolution complète du péricentre à $r_{12} = 3.74 \times 10^{-4}$ sur le 4-corps avec pas dyadique adéquat ($k \ge 14$).
  - Gain de performance mesuré sur BM4 (7 évaluations au lieu de 14).

---

## Iteration 147 Strategy:
### Iteration 147: Phase 4 — Remédiation de la Jauge de Lie et des Correcteurs Symplectiques Traités
- **Scope & Objectives**:
  - `integrators/processed.py` & `system/nbody_system.py` : Développer l'opérateur de connexion de jauge symplectique $\mathcal{J}_{A \to B} = \mathcal{C}_B \circ \mathcal{C}_A^{-1}$ appliqué lors des transitions d'atlas de Jacobi, annulant le saut interfacial $\Delta H = h^{2p}(H_{2p+1}^{(B)} - H_{2p+1}^{(A)})$ dû à la non-commutativité $[X_{K_A}, X_{K_B}] \neq 0$.
  - `integrators/correctors.py` : Remplacer les approximations par différences finies ($\eta = 10^{-7}$) par l'évaluation analytique exacte de la dérivée directionnelle de force (produit hessien-vecteur), éliminant le plancher de bruit à $10^{-9}$.
  - `integrators/correctors.py` : Calibrer analytiquement les constantes homologiques $\kappa_4$ et $\kappa_6$ à partir des coefficients réels de Blanes--Moan plutôt que des constantes génériques de Leapfrog.
  - Valider l'ordre effectif 6 de PBM4 et l'ordre effectif 8 de PBM6.
- **Dependencies & Continuity**: Nécessite l'alignement des macro-pas de l'Itération 146 pour opérer les connexions de jauge aux frontières.
- **Gating Criteria**:
  - Inversion exacte $\mathcal{C}^{-1} \circ \mathcal{C} = \operatorname{Id}$ à $< 10^{-15}$.
  - Continuité de l'énergie de l'hamiltonien d'ombre à travers les transitions d'atlas pour PBM4 et PBM6.
  - Convergence effective vérifiée d'ordre 6 sur PBM4 et d'ordre 8 sur PBM6.

---

## Iteration 148 Strategy:
### Iteration 148: Phase 5 — Certification Globale des 8 Benchmarks, Planches 300 DPI et Appli Démonstrateur
- **Scope & Objectives**:
  - Exécuter la suite complète des 8 benchmarks canoniques (Figure-8, Burrau, Radial collision, 4-corps hiérarchique, Arenstorf, Sitnikov, Jovien, Collision quadruple) sur l'ensemble des 6 intégrateurs.
  - Constater et documenter la chute de l'erreur d'énergie résiduelle sur le système à 4 corps bien en-deçà du plateau de $1.14 \times 10^{-8}$.
  - Régénérer les 4 planches de publication haute résolution 300 DPI (`dashboard_figure8_bm6_300dpi.png`, `dashboard_burrau_bm6_300dpi.png`, `dashboard_4body_bm6_300dpi.png`, `dashboard_multi_integrator_invariants_300dpi.png`).
  - Mettre à jour `showcase_data.json` et resynchroniser `interactive_conormal_showcase.html`.
- **Dependencies & Continuity**: Aboutissement de l'ensemble des remédiations des Itérations 144 à 147.
- **Gating Criteria**:
  - Erreur d'énergie sur le 4-corps hiérarchique réduite d'au moins 2 ordres de grandeur ($< 10^{-10}$).
  - Les 4 planches 300 DPI générées avec calibration typographique et scorecards 100% conformes.
  - 100 % des tests unitaires et d'intégration validés.

---

## Multi-Iteration Checklist:
- [x] Iteration 144: Phase 1 — Remédiation et Certification des Solveurs Universels du Chapitre 3 (Stumpff, Sundman-Halley, McGehee)
- [x] Iteration 145: Phase 2 — Remédiation des Super-Nœuds et Régularisation du Potentiel d'Apex du Chapitre 4 (Screening, Schmitt Hysteresis, Horloge Dépliée)
- [x] Iteration 146: Phase 3 — Remédiation de la Synchronisation de Skeel-Gear et de l'Échéancier Dyadique du Chapitre 5
- [x] Iteration 147: Phase 4 — Remédiation de la Jauge de Lie et des Correcteurs Symplectiques Traités
- [x] Iteration 148: Phase 5 — Certification Globale des 8 Benchmarks, Planches 300 DPI et Appli Démonstrateur
