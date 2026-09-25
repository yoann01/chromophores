# ADR-0007 : Absorption de fond calibrée sur mesures (ISSA + NIST)

- Statut : Proposé
- Date : 2026-09-25
- Décideurs : Y. Granier
- Complète : [ADR-0002](0002-modele-direct-table-sans-dimension.md) (point 1, modèle de peau de référence)

## Contexte

L'écart avec Aliaga & Jarabo / BioSkin venait surtout de l'absorption de fond du tissu
non pigmenté et non vascularisé :
- nous : Jacques, 7,84·10⁸ λ⁻³·²⁵⁵ cm⁻¹ ;
- BioSkin : quasi nulle.

Aucun des deux modèles n'était validé contre des mesures. Deux jeux mesurés sont
désormais dans le dépôt : ISSA (15 256 spectres, 8 groupes) et NIST (100 spectres,
250–2500 nm). Résultats : [MESURES](../phase1/MESURES.md).

## Décision

1. **Valeur par défaut** : `SkinModel.baseline_scale = 0.5` (moitié de la loi de
   Jacques). C'est le seul niveau qui ajuste le mieux ISSA *et* NIST avec une diffusion
   nominale (échelle de μs′ ≈ 0,95).
2. **Règle générale** : tout changement du modèle physique (spectres, fond, eau,
   lipides, couches) est **validé sur ISSA + NIST** avec les scripts
   `scripts/phase1/issa_baseline.py` et `nist_baseline.py`, avant adoption.
3. **Empaquetage vasculaire** : option disponible (`vessel_radius_cm`) mais désactivée,
   parce qu'elle dégrade l'ajustement des mesures.
4. **Défauts connus, suivis comme dette de modèle** : contraste des bandes de
   l'hémoglobine, lipides (930 nm), eau (970 nm), bande de Soret.

## Alternatives considérées

| Alternative | Raison du rejet |
|---|---|
| Garder Jacques × 1 | Diffusion à gonfler de 50 % (ISSA) ou en butée (NIST visible) ; ne produit pas les peaux claires mesurées |
| Aligner sur BioSkin (fond ≈ 0) | Pires ajustements sur les deux jeux ; diffusion en butée pour 43 % des spectres ISSA |
| Ajuster le fond par groupe ethnique ou par zone | Sur-paramétrage : le fond est une propriété du tissu, pas du phénotype. La variabilité doit passer par les chromophores et μs′ (a priori de la phase 2) |

## Conséquences

- **Aucun recalcul de table** : le fond est en amont des nombres sans dimension.
- Les rapports de phase 1 (V1, V2, V2b) restent valides : ils portent sur la table et
  le mélange, pas sur le choix du fond. L'écart V2b avec BioSkin diminue sans
  disparaître : il reste les différences de sang et de phéomélanine.
- La phase 2 part du modèle calibré, et ISSA fournit l'a priori de population.
