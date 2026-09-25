# ADR-0002 : Modèle direct — peau bicouche → table sans dimension → mélange K=2 par cellule

- Statut : Accepté
- Date : 2026-09-25
- Décideurs : Y. Granier

## Contexte

Le random walk du moteur attend, pour chaque longueur d'onde, les paramètres (α, σt, g)
d'un ou plusieurs milieux homogènes. Il faut donc un moyen de passer des chromophores à
ces paramètres.

- **Aliaga & Jarabo** ajustent un mélange de K milieux **par ton et par λ** (9,4 millions
  de couples), puis distillent le résultat dans un réseau. Résultats : K=1 est
  insuffisant, K=2 atteint le plafond de qualité, K=3 n'apporte rien de mesurable.
- **Notre prototype** montre qu'à λ fixé, pour une peau bicouche (couches de même indice,
  même diffusion μs et g dans les deux couches, derme semi-infini), la réponse diffuse ne
  dépend que de nombres sans dimension. Une table 3D de 12×14×6 cellules (≈ 5 min CPU)
  reproduit le Monte Carlo direct avec une erreur médiane de 2 % (max 13 %).
- **Un seul milieu homogène** ajusté par moments donne un écart de profil de 5 à 38 %.
  Le pire cas est une peau foncée à 550 nm, où la queue du profil est sous-estimée.

## Décision

1. **Modèle de peau de référence** : deux couches, paramétrisation Donner & Jensen
   (2006) étendue et alignée sur Aliaga & Jarabo pour pouvoir comparer :
   - plage 250–1000 nm ;
   - g(λ) = 0,62 + 0,00029·λ ;
   - eau dans le derme ;
   - carotène et bilirubine tabulés (PhotochemCAD), fixés à des valeurs basales mais
     modifiables ;
   - μs′ de Jacques avec un **facteur d'échelle libre** (voir ADR-0003 pour son a priori).
2. **Table sans dimension en 4D** : τ_e = μa,épi·d_épi, μa,derme/μs′, d_épi·μs′ et g.
   - Axes logarithmiques pour les trois premiers ; g sur [0,6 ; 0,95].
   - Taille cible d'environ 16×20×8×6 cellules (≈ 15 000).
   - Le rapport de diffusion épiderme/derme est **fixé à 1**. S'il devient libre, il
     ajoute une 5ᵉ dimension (question ouverte).
3. **Contenu d'une cellule** : réflectance totale A, moments radiaux ⟨r⟩ et ⟨r²⟩, et un
   profil R(r) sur ~64 bins logarithmiques, le tout en unités de 1/μs′.
4. **Mélange K=2 ajusté une fois par cellule**, pas par ton.
   - Paramètres : (α₁, g₁, σt,₁·/μs′, α₂, g₂, σt,₂·/μs′, w).
   - Perte : RMSE du log-profil + pénalité sur A (comme Aliaga & Jarabo), initialisée par
     moment matching.
   - Pour ce fit, on utilise une LUT homogène (α, g) → (R, P(r)), invariante d'échelle.
5. **Remplissage** :
   - phase 1 : Monte Carlo CPU (numba), distribué sur la farm ;
   - vérification croisée : *adding-doubling* multicouche (Prahl, `iadpython`) pour A, et
     diffusion bicouche dans le domaine fréquentiel pour les moments ;
   - phase 3 : **re-remplissage par le moteur maison** (rendu d'un slab sous éclairage
     ponctuel), pour que la table et le rendu partagent exactement les mêmes conditions
     aux limites.
6. **Évaluation au hit** : μa(λ) analytique → coordonnées sans dimension → interpolation
   multilinéaire (en log) → paramètres du mélange, remis à l'échelle par μs′(λ).

## Alternatives considérées

| Alternative | Raison du rejet |
|---|---|
| Dataset par ton × λ (Aliaga & Jarabo) | ~600× plus de simulations pour la même information |
| Un seul milieu homogène (K=1, moment matching) | Queue du profil fausse sur peau foncée (38 %) ; conservé comme repli bon marché et comme initialisation |
| K=3 | Aucun gain mesurable chez Aliaga & Jarabo ; +50 % de coût de stockage et de sélection |
| Modèles analytiques seuls (dipôle / multipôle) | Biais de 12 à 51 % mesuré contre le Monte Carlo ; conservés pour la vérification |
| Random walk réellement multicouche dans le moteur | Plus exact, mais demande une géométrie « coque » épidermique et coûte plus cher. Gardé comme **référence de validation** (ADR-0005) |

## Conséquences

- **Positives**
  - La table ne dépend ni des spectres des chromophores ni de la plage spectrale :
    changer un spectre ou ajouter un chromophore ne demande pas de recalcul.
  - Petite (quelques centaines de Ko à quelques Mo) : un seul asset global.
  - Coût de production : ~15 000 × 10⁵–10⁶ photons, soit des heures sur quelques
    machines, ou des minutes sur la farm.
- **Négatives / risques**
  - Hypothèses figées : incidence normale et bord de Fresnel lisse (même limite que
    Aliaga & Jarabo en incidence rasante), n = 1,4, derme semi-infini.
  - Interpolation 4D : des zones à forte courbure (τ_e élevé) peuvent demander un
    raffinement local de la grille.
  - Les paramètres du mélange peuvent être multimodaux entre cellules voisines
    (dégénérescences). Il faut imposer une continuité d'une cellule à l'autre (fit
    initialisé depuis les voisines, ordre canonique des composantes).

## Questions ouvertes

- La diffusion est-elle identique dans l'épiderme et le derme ? Sinon, faut-il une 5ᵉ
  dimension ou une hypothèse de rapport fixe ?
- Faut-il une dimension d'angle d'incidence plus tard, pour la translucidité en rasant ?
- Où s'arrêter côté UV et NIR pour la production ? (Le moteur échantillonne-t-il au-delà
  de 380–780 nm ?)
