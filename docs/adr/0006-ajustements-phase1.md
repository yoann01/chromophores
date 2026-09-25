# ADR-0006 : Ajustements issus de la phase 1 (mélange à g partagé, interpolation préservant l'albédo, critère V2 relatif au bruit)

- Statut : Proposé
- Date : 2026-09-25
- Décideurs : Y. Granier
- Amende : [ADR-0002](0002-modele-direct-table-sans-dimension.md) (points 4 et 6), [ADR-0005](0005-validation-et-criteres.md) (niveau V2)

## Contexte

La mise en œuvre de la phase 1 a révélé quatre problèmes :

1. **Interpolation instable des paramètres du mélange.** Chaque cellule est ajustée
   correctement : erreur sur A médiane de 0,02 %, profil au niveau du bruit Monte Carlo.
   Mais les paramètres K=2 libres (α₁, g₁, s₁, α₂, g₂, s₂, w) admettent plusieurs
   solutions quasi équivalentes, et le fit passe de l'une à l'autre entre cellules
   voisines. Interpolés en multilinéaire, ils produisent aux points milieux une erreur
   d'albédo médiane de 3 à 11 %, et jusqu'à 37 % sur une peau moyenne.
2. **Portée de la LUT homogène.** Avec α ≥ 0,1, les cellules très sombres (A < 10⁻³,
   épiderme très pigmenté en UV/bleu) ne sont pas atteignables.
3. **Critère V2 inatteignable tel quel.** À 3·10⁴ photons par cellule, le **plancher de
   bruit** du Monte Carlo (écart entre deux tirages indépendants d'une même cellule)
   vaut déjà 0,15–0,28 en RMSE log non pondérée et 0,04–0,23 pondérée par l'énergie.
   Le seuil de 0,05 de l'ADR-0005 n'est donc pas mesurable sans référence au bruit.
4. **Convention de réflectance à expliciter.** Les tables donnent A par photon
   *transmis* dans le tissu. L'albédo vu par une caméra en polarisation croisée vaut
   (1 − R_sp)·A, avec R_sp ≈ 2,8 % pour n = 1,4.

## Décision

1. **Mélange K=2 à g partagé.** Les deux milieux prennent l'anisotropie physique du
   tissu, g(λ). Il reste 5 paramètres par cellule (α₁, s₁, α₂, s₂, w). Cela réduit les
   dégénérescences, pour une perte de qualité par cellule négligeable (RMSE énergie
   médiane de 0,088 contre 0,082).
2. **Lissage régularisé entre cellules.** Après le fit initial, on fait 3 passes
   rouge-noir : chaque cellule est ré-ajustée avec une pénalité √(μ·n)·(x − x̄_voisins),
   avec μ = 0,02.
3. **Interpolation qui préserve l'albédo** (moteur et référence Python). Après
   interpolation, on décale tous les u_k = log(1 − α_k) d'un même δ pour que
   Σ w_k A_h(α_k, g) = A_table. L'albédo, donc la couleur, est alors exact par
   construction, et seule la forme du profil dépend de l'interpolation. Dans le
   moteur, il suffit de quelques itérations de Newton ou de bisection sur la LUT 1D.
4. **LUT homogène v2** : α de 0,001 à 1 − 2·10⁻⁵ (48 valeurs).
5. **Critère V2 reformulé** :
   - métrique principale : RMSE du log-profil **pondérée par l'énergie** (part de
     l'énergie réfléchie par bin), qui reflète ce qui se voit au rendu ;
   - critère : écart au plancher de bruit de la référence < 0,05 (p95), **et**
     ΔE00 d'albédo < 0,5 ;
   - la RMSE non pondérée et le plancher de bruit sont toujours publiés à côté.
6. **Convention** : tables et LUT en réflectance par photon transmis ; les albédos
   colorimétriques sont multipliés par (1 − R_sp).

## Alternatives considérées

| Alternative | Raison du rejet |
|---|---|
| Garder g libre par milieu (ADR-0002 initial) | Interpolation instable (pire cas : 132 % d'erreur d'albédo aux points milieux contre 24 % avec g partagé, avant lissage et correction) |
| Plus de photons pour baisser le plancher de bruit | Utile, mais relève de la farm (×10 photons = ×10 en temps) ; ne règle pas la multimodalité |
| Grille plus fine | Réduit l'interpolation, mais pas la multimodalité ; prévu pour la table de production (farm) |
| Réseau pour lisser le mapping cellule → paramètres | Contraire à l'ADR-0001 tant qu'une solution explicite suffit |

## Conséquences

- La couleur ne dépend plus de l'interpolation du mélange. L'erreur résiduelle vient
  seulement de la table multicouche (V1).
- Le moteur doit implémenter la correction d'albédo : une petite boucle par longueur
  d'onde héroïque, sur une LUT 1D.
- On perd la liberté d'un g par milieu. Si la qualité de profil devient limitante sur
  peau foncée, on pourra revenir à un g libre pour le 2ᵉ milieu seulement, avec lissage.

## Questions ouvertes

- Budget photons et résolution de la table de production (farm) : viser un plancher de
  bruit < 0,02 ?
