# ADR-0001 : Séparer le modèle direct (bien posé, tabulé) de l'inversion (mal posée, probabiliste)

- Statut : Proposé
- Date : 2026-09-25
- Décideurs : Y. Granier (à valider)

## Contexte

Objectif : rendu de production de digital doubles. On part d'albédos RGB calibrés
(TexturingXYZ VFace), rendus par un moteur spectral maison en hero wavelength et en
shade-before-hit (inspiré de Manuka). Le SSS doit être piloté par des paramètres
biophysiques (chromophores), et non par un rayon de diffusion RGB réglé à la main.

L'état de l'art le plus proche (Aliaga & Jarabo 2026, voir
[topo](../TOPO_ALIAGA_JARABO_2026.md)) traite toute la chaîne RGB → paramètres de rendu
spectraux avec un **réseau chaîné**. Ce réseau est entraîné sur 25 000 tons simulés par un
**Monte Carlo GPU** (≈ 10¹³ photons, plus 46 h de fit sur A100). Nous n'avons ni cette
infrastructure GPU ni ce dataset.

Nos analyses ([AXE_DE_RECHERCHE](../AXE_DE_RECHERCHE.md), [PRODUCTION](../PRODUCTION.md))
montrent deux choses :
- **La chaîne se sépare en deux problèmes de nature différente :**
  - *chromophores → paramètres de random walk* est **bien posé**. C'est un calcul direct,
    que la réduction sans dimension rend tabulable indépendamment des spectres et de λ ;
  - *RGB → chromophores* est **mal posé**. Il y a 3 mesures pour 5 à 8 paramètres. Le RGB
    contraint environ 2 à 3 directions (charge de mélanine, sang, un peu SO₂) ; les
    métamères RGB couvrent 97 % de la plage de μs′.
- **Un réseau de bout en bout mélange les deux problèmes.** L'a priori y est implicite
  (biais inductif, distribution d'entraînement), sans incertitude, et on ne peut ni
  l'auditer ni l'ajuster par acteur.

## Décision

On architecture le système en **deux blocs indépendants, reliés par l'espace des
chromophores** :

1. **Bloc direct (déterministe, exact à la discrétisation près).**
   Chromophores + λ → coefficients optiques → nombres sans dimension → table précalculée
   → paramètres d'un mélange de milieux pour le random walk. Aucun apprentissage.
   Détails : [ADR-0002](0002-modele-direct-table-sans-dimension.md).
2. **Bloc inverse (probabiliste, a priori explicite).**
   RGB calibré → chromophores (estimation MAP) **plus une incertitude**. Il n'estime que
   le sous-espace identifiable ; le reste vient d'a priori hiérarchiques (population,
   région, acteur). Détails : [ADR-0003](0003-inversion-rgb-lut3d-a-priori.md).

L'interface entre les deux blocs est un **vecteur de chromophores**. C'est aussi ce que
stocke le moteur ([ADR-0004](0004-representation-moteur-shade-before-hit.md)) et ce que
les artistes peuvent éditer.

## Alternatives considérées

| Alternative | Raison du rejet (pour l'instant) |
|---|---|
| Réseau chaîné de bout en bout (Aliaga & Jarabo) | Demande GPU et dataset massifs ; a priori implicite ; pas d'incertitude ; sortie de 376 × 12 valeurs par texel, incompatible avec un stockage compact par vertex |
| Statu quo : *albedo inversion* RGB (Chiang 2016 / Wrenninge 2017) + rayon artistique | Ne résout ni la corrélation couleur/translucidité ni le spectral |
| Inversion Monte Carlo par texel | Coût prohibitif (10⁷–10⁸ texels) |
| Modèle direct appris par réseau (surrogate) | Inutile puisque la table est exacte et petite ; peut revenir si la dimension de la table explose (> 5D) |

## Conséquences

- **Positives**
  - Chaque bloc se valide séparément : le direct contre le Monte Carlo, l'inverse contre
    des mesures ou des métamères.
  - L'a priori devient un **objet explicite**, versionnable et ajustable par acteur ou par
    région.
  - Pas de dépendance à un GPU d'entraînement : tout le calcul lourd est parallélisable
    sur la render farm (CPU).
  - Le lissage spectral est garanti par construction, puisque les paramètres dérivent de
    μa(λ) physique.
- **Négatives / risques**
  - La table directe suppose un modèle de peau figé (2 couches, même diffusion dans les
    deux couches, incidence normale). Changer le modèle impose de recalculer la table.
  - L'inversion avec a priori explicite demande de construire ces a priori : c'est du
    travail de données (VFace, mesures ponctuelles).
  - On perd l'éventuel bénéfice « régularisant » d'un réseau partagé entre tons, noté par
    Aliaga & Jarabo. À compenser par la régularisation spatiale et les a priori.

## Questions ouvertes

- Faut-il exposer aux artistes les chromophores bruts, ou un espace dérivé plus intuitif
  (teinte de mélanine, rougeur, translucidité) ?
- Où placer la frontière si une 3ᵉ couche (derme papillaire / plexus) devient nécessaire ?
