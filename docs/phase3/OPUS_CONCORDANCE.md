# Concordance avec le moteur Opus (modèle direct)

Date : 2026-09-30. Mesures faites côté Opus par l'IA du projet moteur (outil `dump_strat`, rendus de
dalle stratifiée) ; ajustements sur mesures faits ici. Aucun code de production d'Opus modifié.

## Étape 1 : coefficients, couche par couche

- **Derme** : concordance à 1,00–1,01 une fois le fond et l'eau alignés. Hémoglobine, bilirubine et
  carotène concordent. Aux réglages par défaut, le derme d'Opus absorbe ×1,66 à 650 nm (fond ×1,
  sans eau).
- **Épiderme** : l'écart de 21–28 % venait de la **constante de phéomélanine de Python, 10 fois trop
  faible** (Donner & Jensen 2006, éq. 4, en mm⁻¹). Opus est fidèle au papier. Corrigé ici
  ([PHEO_CHECK](../phase2/PHEO_CHECK.md)).
- **Diffusion** : μs′ d'Opus ≈ 0,76 × le nôtre à 550–650 nm (0,86 à 450 nm), pente plus forte.
- Autres différences : couche cornée de 18 µm (Opus seulement), derme isotrope dans Opus, épiderme
  par défaut 0,25 mm (a priori visage ici : 0,21 mm).

## Étape 2 : transport (dalle plate stratifiée)

| Comparaison | ΔE00 (5 tons) |
|---|---|
| Opus contre MC Python, **mêmes coefficients**, rebonds illimités | **0,2 – 0,35** (Opus ≈ 1,5 % plus sombre) |
| Idem avec `max_depth = 64` (défaut) | jusqu'à 0,9 sur la peau la plus claire (≈ 2 %) |
| Modèle Python par défaut (avant correction) contre Opus | 4,3 – 4,9 |
| Python, phéomélanine corrigée | 2,6 – 4,4 |
| Python, phéomélanine corrigée, fond ×1, sans eau | 2,1 – 3,2 |

**Le transport stratifié d'Opus concorde avec la référence Python** (critère albédo de l'ADR-0008
rév. 2, ΔE00 < 0,5, satisfait sur dalle). Les écarts restants viennent des constantes et du modèle.
Contrôle hors peau : dalle grise homogène conforme à Chandrasekhar à < 0,5 %.

Pièges du montage (côté test, pas d'Opus) : `SPECTRAL_SSS_CMP_SCALE` divise σ (1 unité = 1/scale mm ;
dalle semi-infinie de 40 mm → `SCALE=0.1`) ; le blanc de référence demande
`SPECTRAL_MAX_DEPTH=100000`.

## Le fond ×0,5 et l'eau survivent-ils à la correction ?

Voir [MODEL_VARIANTS](../phase2/MODEL_VARIANTS.md) (298 ISSA 400–700 nm, 100 NIST 400–1000 nm) :

| Modèle | RMS NIST médian | Sang ajusté (NIST) | Échelle μs′ ajustée (NIST) |
|---|---|---|---|
| Fond ×0,5, eau (actuel) | **0,0132** | 1,7 % | 0,96 |
| Fond ×1, eau | 0,0179 | 3,0 % | 1,62 |
| Fond ×0,5, sans eau | 0,0275 | 10,9 % | 2,17 |
| Fond ×1, sans eau (Opus) | 0,0284 | 13,2 % | 2,99 |

- **Oui, les deux survivent.** Sans eau, l'ajustement du proche infrarouge devient deux fois moins
  bon et absurde (sang 11–13 %, μs′ × 2–3) : le sang sert à imiter les bandes de l'eau.
- Dans le visible seul (ISSA), fond et eau ne se distinguent pas : le fond ×1 est compensé par une
  diffusion × 1,6 et un épiderme plus mince. Or **Opus combine fond ×1 et μs′ × 0,76** : les deux
  écarts assombrissent le rouge dans le même sens.
- Conséquence : l'alignement doit se faire **côté Opus** (fond, eau, puis μs′), par un ADR. Il
  reste à sourcer μs′ et le fond de Jacques 2013 (données fournies par l'utilisateur, en attente)
  avant de le rédiger.

## Suite

1. Isoler un par un les écarts restants (μs′, couche cornée, g du derme, carotène) avec un
   `SkinModel` « constantes Opus » en Python : le transport étant concordant, le MC Python peut
   servir de substitut au moteur. Il faut pour cela les formules exactes d'Opus.
2. ADR d'alignement des constantes (fond, eau, μs′), une fois Jacques 2013 intégré.
3. `max_depth = 64` : biais d'environ 2 % sur peau claire, à surveiller (roulette russe plutôt
   qu'une coupure franche ?).
