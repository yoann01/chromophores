# ADR-0004 : Représentation moteur — chromophores par vertex, évaluation au hit, mélange K=2 en hero wavelength

- Statut : Accepté
- Date : 2026-09-25
- Décideurs : Y. Granier

## Contexte

Le moteur maison est spectral, en hero wavelength (4 longueurs d'onde par chemin, Wilkie
et al. 2014), avec une architecture shade-before-hit inspirée de Manuka (Fascione et al.
2018) :
- **au dicing** : les paramètres de matériau sont évalués puis stockés par vertex de
  micropolygone ;
- **au hit** : seule l'évaluation spectrale a lieu, aux longueurs d'onde du chemin.

Il faut donc une représentation compacte par vertex et évaluable à n'importe quel λ.

Aliaga & Jarabo produisent 376 × 12 valeurs par texel. Ils sélectionnent un milieu par
random walk, sans traiter le cas où poids et extinction varient entre les longueurs
d'onde d'un même chemin.

## Décision

1. **Stockage par vertex** : 5 à 6 scalaires en half-float :
   - charge de mélanine ;
   - β mélanine ;
   - fraction sanguine ;
   - SO₂ ;
   - échelle de μs′ ;
   - épaisseur de l'épiderme (ou charge de mélanine et Cm séparés, selon l'ADR-0003).

   Carotène et bilirubine sont des constantes de matériau (ou des cartes optionnelles).
2. **Évaluation au hit, pour chacune des 4 longueurs d'onde** :
   μa,épi(λ), μa,derme(λ), μs′(λ) (analytique + LUT 1D pour l'hémoglobine et l'eau) →
   coordonnées sans dimension → **table globale de l'ADR-0002** → (α_k, g_k, σt,k, w_k)
   pour k = 1, 2.
3. **Sélection du milieu et MIS spectral** :
   - on tire le milieu k avec la probabilité w_k(λ_hero) ;
   - les longueurs d'onde secondaires λ_j sont pondérées par la combinaison MIS (balance
     heuristic) des pdf de sélection et de marche sur les 4 longueurs d'onde.
4. **Échantillonnage des distances** : on utilise le σt,k de λ_hero, et les autres
   longueurs d'onde reçoivent des poids par rapport des pdf (*spectral tracking* /
   one-sample MIS). Sur peau moyenne, σt varie d'un facteur ≈ 2 entre 450 et 750 nm : le
   problème est bien conditionné, mais c'est à mesurer.
5. **Filtrage** (mip / dicing) : on filtre dans un **espace linéaire en absorption**
   (épaisseurs optiques : charge de mélanine, charge de sang), pas dans les paramètres
   bruts, pour limiter le biais non linéaire (taches, pores, barbe).
6. **Contrôles artistiques** : ce sont des multiplicateurs ou décalages sur les
   paramètres stockés (rougeur, bronzage, translucidité), évalués au dicing. Ils ne
   touchent jamais la table.

## Alternatives considérées

| Alternative | Raison du rejet |
|---|---|
| Stocker des paramètres spectraux précalculés par texel (Aliaga & Jarabo) | 376 × 12 valeurs, ou une base spectrale : stockage et bande passante incompatibles avec le shade-before-hit |
| *Uplifting* RGB → spectre (Jakob & Hanika 2019) des paramètres SSS | Perd la physique : les bandes de l'hémoglobine ne sont pas reproduites |
| Évaluer un réseau au hit | Coût par hit et dépendance d'inférence dans le moteur |
| Un seul milieu (K=1) | Repli bon marché autorisé (LOD, rendus lointains), pas le défaut |

## Conséquences

- **Positives**
  - L'empreinte mémoire par vertex est celle d'environ 2 textures RGB.
  - Le coût au hit est une douzaine d'opérations analytiques plus une interpolation 4D
    par longueur d'onde.
  - Les paramètres ont un sens physique, ce qui en fait des contrôles artistiques
    intuitifs.
- **Négatives / risques**
  - Le moteur doit gérer le MIS spectral pour la sélection du mélange et le tracking
    chromatique. C'est un chantier en soi s'il n'existe pas déjà.
  - La variance peut augmenter quand w_k(λ) varie fortement entre les longueurs d'onde
    d'un même chemin (bleu contre rouge sur peau foncée). À mesurer.
  - Le modèle impose une géométrie « semi-infinie » au SSS : les parties minces (oreilles,
    narines) restent traitées par le random walk géométrique habituel du moteur. Il faut
    vérifier que la table (derme semi-infini) reste acceptable dans ces zones.

## Questions ouvertes

- Le random walk actuel gère-t-il déjà un σt chromatique par longueur d'onde héroïque ?
- Budget mémoire par vertex acceptable (5 ou 6 half) ?
- Le moteur échantillonne-t-il hors du visible (UV/NIR), ou s'en tient-on à 380–780 nm ?
