# Topo : Masse & Walster 2023, *Artist Friendly Physically Plausible Attenuation*

Framestore, DigiPro '23, [doi:10.1145/3603521.3604298](https://doi.org/10.1145/3603521.3604298). PDF :
[`pubs/3603521.3604298.pdf`](pubs/3603521.3604298.pdf).

## 1. Résumé

### Problème visé
Dans le paramétrage artiste albédo / atténuation (RGB), les deux cartes sont indépendantes. On peut
alors créer des coefficients de diffusion impossibles, qui donnent une diffusion verdâtre, des
dégradés artificiels et des inversions de couleur en transmission.

### Idée
Le coefficient de diffusion de la nature est contraint : il se situe entre **Mie** (grosses
particules, blanc) et **Rayleigh** (petites particules, en λ⁻⁴). Un seul scalaire *w* règle le
mélange : s = w·r + (1 − w)·m, avec r ≈ (0,553 ; 0,808 ; 1,64) en sRGB et m = (1, 1, 1).
L'absorption, elle, n'est pas contrainte.

### Trois méthodes
| Méthode | Hypothèse | Usage |
|---|---|---|
| Diffusion constante | s fixé par *w* ; l'albédo pilote l'absorption (albédo sombre = plus opaque) | Peau (w ≈ 0,64, tiré de Jacques 2015) |
| Absorption constante | absorption blanche ; σs projeté sur la direction de s | Glace, eau, cristaux |
| Injection de Rayleigh | atténuation normalisée (l'opacité ne vient que de la densité) | FX (fumée, explosions) |

Chaîne de la méthode « peau » : albédo multiple A → albédo simple α (polynôme de Chiang 2016) →
atténuation = α / s.

## 2. Mise en regard avec le projet

### Convergences
- **Diffusion fixe, couleur portée par l'absorption** : c'est notre conclusion de phase 2 (la
  diffusion n'est pas identifiable depuis le RGB, il faut une constante).
- Leur Figure 2 (albédo et atténuation pour des plages de mélanine et d'hémoglobine) part du même
  point que nous. Ils y ont renoncé faute de pouvoir extraire les concentrations des images de
  référence : c'est ce que fait notre pipeline de phase 2.
- **L'artiste agit sur l'albédo, la transmission suit.** C'est notre contrainte de production.

### Limites
1. **Milieu unique.** La mélanine devient une absorption répartie dans tout le volume. Mesuré
   ([FILTER_FEASIBILITY](phase3/FILTER_FEASIBILITY.md)) : à albédo identique, le rayon moyen du
   profil tombe à 0,36–0,67 fois la référence bicouche sur peau mate ou foncée. La peau foncée
   devient trop opaque, alors que le transport rouge dans le derme reste en réalité en grande
   partie intact. Sur peau très claire, en revanche, le milieu unique est excellent.
2. **Spectre à 3 échantillons.** Rayleigh est évalué aux longueurs d'onde « dominantes » des
   primaires RGB. Notre modèle utilise μs′ = 46·(λ/500)^−1,421 (Jacques), qui est déjà un mélange
   Mie + Rayleigh condensé en une loi de puissance.
3. **Conventions.** Chez Jacques, *f*Ray est une fraction de μs′ (diffusion réduite) à 500 nm ;
   leur *w* mélange des σs normalisés par la moyenne RGB. Ce ne sont pas les mêmes grandeurs.
   Rayleigh (g ≈ 0) et Mie (g ≈ 0,9) n'ont pas la même anisotropie, alors que la fonction de phase
   reste réglée à la main (ils le notent au §8.2).
4. **Albédo multiple → simple** par un polynôme ajusté pour un réglage précis. Notre LUT (α, g)
   (ADR-0006) tient compte de g et de Fresnel.

## 3. Ce qu'on en retient

- **Base de comparaison** pour les rendus de phase 3 (même albédo, milieu unique contre bicouche).
- **Repli pour les texels hors peau** (masque de phase 2) : la version spectrale de la « diffusion
  constante » donne une réponse plausible sans inventer de chromophores.
- **Argument de communication** : les studios corrigent ce problème par des contraintes ad hoc ;
  une représentation par chromophores l'évite par construction.
