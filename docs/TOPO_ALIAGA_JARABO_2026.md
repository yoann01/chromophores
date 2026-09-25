# Topo : Aliaga & Jarabo 2026, *Spectral Subsurface Scattering from RGB via Biophysical Skin Inversion*

Meta Reality Labs Research, arXiv 2606.27604v1 (25 juin 2026). PDF :
[`pubs/aliaga_2606.27604v1.pdf`](pubs/aliaga_2606.27604v1.pdf).

## 1. Résumé

### Problème visé
En production, le SSS path-tracé utilise un milieu homogène unique (α, σt, g) :
- α est obtenu par *albedo inversion* (Chiang et al. 2016, Wrenninge et al. 2017) ;
- la distance de diffusion et l'anisotropie sont réglées à la main.

Le papier y voit deux défauts :
1. **Distance et couleur sont décorrélées**, alors que chez l'humain elles sont liées par
   les chromophores.
2. **Un seul milieu ne représente pas la peau multicouche.**

### Pipeline
1. **Modèle biophysique** (celui d'Aliaga et al. 2023) :
   - deux couches : épiderme de 50 à 350 µm sur un derme semi-infini ;
   - diffusion de Jacques, avec g(λ) = 0,62 + 0,00029·λ ;
   - bilirubine, β-carotène et eau fixés à des valeurs basales ;
   - 5 paramètres libres : p = (f_mel, f_blood, t_epi, SO₂, b_mel).
2. **Dataset Monte Carlo GPU** :
   - 25 000 tons de peau pour l'entraînement (Halton) et 5 000 pour la validation (uniforme) ;
   - 376 longueurs d'onde de 250 à 1000 nm, 10⁶ photons par λ, sur 2 A100 ;
   - on garde R_d(λ) et le profil R_d(r, λ) sur 512 bins logarithmiques jusqu'à 5 cm.
3. **Ajustement d'un mélange de K milieux** :
   - LUT homogène P(r | α, g) de 800×800 et R(α, g), invariantes d'échelle en σt :
     P(r | α, σt) = σt² · P(σt·r | α, 1) ;
   - fit par ton et par λ avec Adam, perte = RMSE du log-profil + 5 × erreur de réflectance ;
   - progression K=1 → 2 → 3, en 46 h sur un A100.
4. **Distillation en un réseau chaîné**, entraîné par étapes gelées :
   - `Enc(RGB) → p` (5D) : entraîné sans supervision des paramètres, uniquement sur la
     reconstruction de couleur. L'exposition γ = Y(A_RGB)/Y(Â) est analytique, donc le
     réseau n'exploite de fait que la **chromaticité** ;
   - `Dec_R(p) → R(λ)` sur 376 longueurs d'onde (architecture d'Aliaga 2023) ;
   - `Dec_T(p, R̂) → θ(λ)` : un MLP à 4 couches qui produit 376 × 12 valeurs, soit
     (α, g, σt, w) pour chacun des 3 milieux ;
   - pertes : Huber sur les paramètres, réflectance via la LUT (ΔE de 11 → 0,062) et
     lissage spectral ;
   - un facteur de « profondeur de slab » appris, c(λ), sert à obtenir un gradient sur σt.
5. **Rendu** : dans Mitsuba 3, chaque random walk tire au départ un milieu selon les
   poids w_k.

### Résultats annoncés
- **Rendu d'un slab** (ΔE76 dans le visible) : K=1 donne 51,4 ; K=2 donne 3,80 ; K=3 donne
  3,79 ; le réseau donne 3,19, donc mieux que le fit par ton.
- **K=2 suffit.** Le passage de K=1 à K=2 réduit l'angle spectral (SAM) de 80 % dans le
  visible (bandes de l'hémoglobine) ; K=3 n'apporte qu'un gain « structurel ».
- **L'anisotropie libre est indispensable** : ΔE 0,5 en isotrope contre 0,019 avec g libre.
- **Le 3ᵉ milieu joue le rôle de l'épiderme** : α ≈ 0,34, concentré sur l'UV et le bleu,
  avec une distribution de g bimodale.

### Limites reconnues
- Hypothèse d'incidence normale : la translucidité est trop « sèche » en incidence rasante.
- Les lèvres sortent de la variété des peaux d'entraînement.
- Aucune validation sur des mesures in vivo.
- La LUT ignore σt.

## 2. Points à clarifier ou critiquer

- **Incohérences internes :**
  - l'encart de la Fig. 6 donne des poids de 0,96 / 0,04 / 0,01, alors que le texte donne
    un poids de 0,35 au 3ᵉ milieu ;
  - l'erreur à K=1 vaut 34,7 dans la légende de la Fig. 5 mais 51,4 dans la Table 1 ;
  - l'étape 3 est étiquetée « DecR » au lieu de « DecT » ;
  - trois ΔE différents (0,062 ; 3,19 ; 0,019) sans que la métrique de chacun soit claire.
- **Validation 100 % in silico**, contre leur propre Monte Carlo, sur un slab en incidence
  normale.
- **« Le réseau bat le plafond d'optimisation »** : cela montre surtout que la perte du fit
  (profil + R) ne correspond pas à la métrique d'évaluation. Ce n'est pas un gain réel.
- **Pas d'incertitude.** L'encodeur choisit un point de l'ensemble des métamères, fixé par
  son biais inductif et par la distribution d'entraînement.
- **La « base de tons de peau » est synthétique.** Ce sont des échantillons de Halton de
  l'espace 5D passés dans leur Monte Carlo, pas des mesures.

## 3. Mise en parallèle avec notre travail

| Aspect | Aliaga & Jarabo | Nous |
|---|---|---|
| Modèle de peau | 2 couches, 5 paramètres, 250–1000 nm, g(λ), eau | 2 couches, 8 paramètres (dont échelle de μs′, carotène, bilirubine), 380–780 nm, g = 0,8 fixe. **Le leur est plus complet spectralement** |
| Diffusion (μs′) | **Fixe (Jacques) pour tout le monde** | Paramètre libre |
| Ambiguïté RGB | Non analysée ; l'encodeur n'utilise que la chromaticité (2 degrés de liberté pour 5 paramètres) | Quantifiée : Cramér-Rao et métamères |
| Homogénéisation | LUT homogène (α, g) + fit d'un mélange par ton et par λ, puis distillation en réseau | Table 3D sans dimension **de la peau multicouche elle-même**, puis moment matching pour 1 milieu |
| Échec de K=1 | ΔE 51, « ne découple pas forme et amplitude » | Écart de profil de 5 à 38 %, pire sur peau foncée dans le vert, mais A et le rayon moyen sont exacts par construction |
| Chaîne de capture | RGB supposé colorimétrique (D65 + CMF) | Calibration caméra : ΔE00 de 2,4 avec une ColorChecker |
| Ombrage résiduel | γ de luminance : robuste, mais perd l'information d'amplitude | Gain nuisible : cohérent avec leur choix. Sang et mélanine restent estimables, la diffusion non |
| Moteur | Mitsuba 3, spectral par λ | Hero wavelength, shade-before-hit |
| Stockage | 376 × 12 valeurs par texel (non discuté) | 5–6 scalaires par vertex + une table globale d'environ 8 Ko |

### Ce qui converge
- **Un seul milieu ne suffit pas**, et l'écart vient de l'épiderme absorbant. Leur 3ᵉ
  milieu « épidermique » correspond à notre pire cas : peau foncée à 550 nm, 38 % d'écart
  sur la queue du profil.
- **Carotène et bilirubine sont inobservables** : ils les fixent, nous les trouvons non
  identifiables.
- **La translucidité doit dépendre des chromophores**, pas d'un réglage artistique
  indépendant.
- **L'exposition doit être découplée**, ce qui rejoint notre analyse du gain nuisible.

### Ce qu'on apporte qu'ils n'ont pas
1. **Leur « translucidité depuis le RGB » est en partie une hypothèse.** Comme la diffusion
   est la même pour tous, la translucidité ne varie qu'avec l'absorption. Nos métamères
   montrent que μs′ est justement le paramètre le moins contraint : 97 % de sa plage à RGB
   identique. La variabilité de la diffusion entre individus, âges et zones du visage est
   absente de leur modèle.
2. **Leur partie directe peut devenir une table.** Ils notent que l'invariance d'échelle ne
   tient pas en multicouche. Notre regroupement en nombres sans dimension (τ_e,
   μa_d/μs′, d_e·μs′) la rétablit. Avec leur g(λ), la table aurait 4 dimensions, peut-être
   5 si la diffusion diffère entre couches. Leurs 9,4 millions de simulations Monte Carlo
   (25 000 tons × 376 λ) et les 46 h de fit deviendraient un fit par cellule de table,
   indépendant des spectres de chromophores. `Dec_T` deviendrait une lookup exacte. Seul
   `Enc` (RGB → paramètres), la partie mal posée, resterait à apprendre ou à régulariser,
   idéalement avec une incertitude.
3. **Les questions de production** qu'ils n'abordent pas :
   - calibration caméra ;
   - stockage par vertex en shade-before-hit ;
   - a priori par acteur et par région (lèvres) ;
   - combinaison du mélange avec le hero wavelength. Comme w_k(λ) et σt,k(λ) varient avec
     λ, choisir le milieu à la longueur d'onde héroïque impose de pondérer les 3 autres
     (MIS spectral).

### Ce qu'on devrait leur reprendre
- **K=2 plutôt que 1** : le 2ᵉ milieu est nécessaire, et un 3ᵉ n'apporte presque rien.
- **g libre** dans le fit.
- **Perte en log-profil**, qui donne autant de poids au pic qu'à la queue.
- **Plage 250–1000 nm, g(λ) et eau**, pour être comparables.
- **Bins radiaux logarithmiques** dans le Monte Carlo.
- **Leurs 10 tons représentatifs (Fig. 2, paramètres et RGB publiés)** comme banc de test
  commun.

### Positionnement proposé
« Séparer le bien-posé du mal-posé ».
- **Partie directe, exacte et tabulée** : chromophores → mélange de 2 milieux, pour chaque
  λ héroïque.
- **Partie inverse, probabiliste** : RGB calibré → chromophores avec incertitude, plus un
  a priori sur la diffusion par acteur ou par région, calibré par quelques mesures
  ponctuelles.
