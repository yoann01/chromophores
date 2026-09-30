# ADR-0008 : Représentation moteur de la phase 3 — chromophores vers le milieu stratifié du moteur ; ni filtre épidermique, ni mélange K=2 par défaut

- Statut : Proposé (révision 2)
- Date : 2026-09-30
- Décideurs : Y. Granier
- Amende (si accepté) : [ADR-0004](0004-representation-moteur-shade-before-hit.md) (points 2 à 4), [ADR-0005](0005-validation-et-criteres.md) (critère bloquant de la phase 3)

## Historique

La **révision 1** proposait deux candidats départagés par un test d'image : A, le mélange K=2
(ADR-0004) ; B, l'épiderme en filtre d'interface avec une épaisseur optique effective τ_eff sur un
random walk dermique. Elle reposait sur une hypothèse fausse : que le moteur ne sait parcourir
que des milieux homogènes. La révision 2 la remplace.

## Contexte

1. **Le moteur sait déjà faire le bicouche exact.** Vérifié dans le code du moteur :
   - le chemin de production `skin_stratified_medium` parcourt épiderme et derme en
     *delta-tracking* stratifié, sans approximation de couche ;
   - le « voile » (`boundary_transmittance`, medium.hpp:337, appliqué à path_integrator.cpp:1466)
     multiplie le débit spectral à chaque **franchissement** de la surface, mais **pas aux
     réflexions internes** (path_integrator.cpp:1440). Il n'est plus actif que sous
     `SPECTRAL_SKIN_LEGACY_VEIL=1`, et l'ADR 0014 rév. du moteur l'a abandonné pour cette raison.
2. **L'option B de la révision 1 approche ce que le moteur fait exactement.** Dans l'étude de
   faisabilité ([FILTER_FEASIBILITY](../phase3/FILTER_FEASIBILITY.md)), le filtre était une couche
   sous l'interface de Fresnel : il s'appliquait donc *aussi* aux réflexions internes. Le voile
   du moteur, qui ne s'y applique pas, ferait encore moins bien que les chiffres mesurés, et un
   τ_eff ajusté sur l'albédo compenserait ce défaut seulement en moyenne, avec une dépendance à la
   géométrie et à l'éclairage.
3. **Le mélange K=2 de l'ADR-0004** servait à émuler le bicouche avec des milieux homogènes.
   Ses deux risques ouverts (MIS spectral sur la sélection du milieu, géométrie mince) et le
   blocage V2-forme n'ont plus de raison d'être si le moteur parcourt directement les deux couches.
4. **L'étude de faisabilité reste utile** comme mesure :
   - le milieu unique (Masse & Walster 2023) comprime le profil de 0,36 à 0,67 fois sur peau
     pigmentée ;
   - un filtre ne reproduit pas le pic de diffusion simple de l'épiderme pigmenté.

   Ces deux résultats justifient le bicouche exact.

## Décision

1. **Pas de filtre épidermique** : l'option B est abandonnée, et le voile historique n'est pas
   réactivé.
2. **Représentation moteur :**
   - par vertex : les chromophores de l'ADR-0004 (point 1, inchangé) ;
   - au hit, pour chaque longueur d'onde du chemin : on calcule analytiquement les coefficients
     de chaque couche avec le modèle `SkinModel` (ADR-0007), soit μa,épi(λ), μa,derme(λ), μs(λ),
     g(λ) et l'épaisseur de l'épiderme ;
   - ces coefficients sont transmis à `skin_stratified_medium`.

   Il n'y a ni table 4D ni mélange au hit.
3. **Le mélange K=2 n'est plus le défaut.** Il reste disponible comme LOD bon marché (rendus
   lointains), comme l'ADR-0004 le prévoyait déjà pour K=1. La table 4D reste l'outil de
   l'**inversion** (phase 2), où seul l'albédo compte, et l'albédo passe V1 et V2.
4. **V2-forme n'est plus bloquant pour la phase 3**, et la table « farm » n'est pas nécessaire
   au rendu.
5. **Nouveau critère bloquant, V4 (validation croisée) :** la référence Python et le milieu
   stratifié du moteur reçoivent les mêmes coefficients, sur une dalle. Deux tests :
   - albédo spectral : ΔE00 < 0,5 ;
   - test de bord (demi-plan éclairé et spot de 1 mm, 400–700 nm, D65) : ΔE00 p95 < 1.

   Le panel de peaux est tiré de l'a priori de population, peaux foncées comprises. Ensuite, sur
   la géométrie courbe ou mince (oreille, narine), c'est le moteur qui fait référence : la
   référence Python ne connaît que la dalle.
6. **Masse & Walster** sert de base de comparaison pour les rendus et de repli pour les texels
   hors peau (masque de la phase 2).

## Alternatives considérées

| Alternative | Raison du rejet |
|---|---|
| Filtre épidermique à τ_eff (révision 1, option B) | Approche ce que le moteur fait exactement ; le voile ne s'applique pas aux réflexions internes ; τ_eff dépendrait de la géométrie et de l'éclairage |
| Mélange K=2 par défaut (ADR-0004) | Émule un bicouche que le moteur sait parcourir ; forme du profil hors critère ; MIS de sélection à écrire ; conservé comme LOD |
| Milieu unique (Masse & Walster) | Profil comprimé sur peau pigmentée (mesuré) |

## Conséquences

- **Positives**
  - Le modèle direct du rendu est exactement celui de la référence : il n'y a plus d'erreur de
    représentation, seulement le bruit Monte Carlo.
  - La phase 3 est débloquée sans calcul sur la farm.
  - Le stockage par vertex est inchangé, et le coût au hit se réduit à quelques évaluations
    analytiques et de LUT 1D (hémoglobine, eau) par longueur d'onde.
- **Négatives / risques**
  - Le coût de rendu est celui du delta-tracking stratifié, déjà payé en production.
  - L'inversion (table 4D, dalle semi-infinie) et le rendu (géométrie réelle) peuvent diverger
    sur les zones minces : l'albédo mesuré y contient de la transmission. Il faudra le mesurer
    sur des rendus, et c'est un sujet pour la phase 4 (mesures).
  - Les conventions doivent être alignées entre la référence Python et le moteur : n = 1,4,
    interface interne sans saut d'indice, g(λ), albédo par photon transmis.

## Questions ouvertes

- Comment `skin_stratified_medium` définit-il la profondeur de la couche : le long de la normale
  d'entrée, ou par distance à la surface ? Sur un épiderme mince et une géométrie courbe, la
  différence peut compter.
- Quelles entrées accepte-t-il (coefficients spectraux par couche, fonction de λ, g par couche) ?
  Accepte-t-il g(λ) ?
- Peut-on l'exécuter hors rendu, sur une dalle, pour V4 ?
- La diaphonie mélanine / sang (clôture de la phase 2) reste un sujet séparé.
