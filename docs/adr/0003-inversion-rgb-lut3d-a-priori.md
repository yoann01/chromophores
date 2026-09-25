# ADR-0003 : Inversion RGB → chromophores par LUT 3D, avec a priori explicite et incertitude

- Statut : Proposé
- Date : 2026-09-25
- Décideurs : Y. Granier (à valider)

## Contexte

L'entrée de production est un albédo VFace calibré, présumé en polarisation croisée (à
confirmer). Nos analyses montrent :

- **Le RGB contraint** la charge de mélanine (Cm × épaisseur, avec un ratio eu/phéo mal
  séparé) et le volume sanguin. Sur peau claire, il contraint aussi un peu la SO₂ et la
  diffusion.
- **Sur peau foncée, le RGB ne dit presque rien du derme.**
- **La couleur diffuse est bien contrainte**, même sous d'autres illuminants
  (métamères ≤ 1,6 ΔE00). **La translucidité ne l'est pas** (μs′ libre à 97 %).
- **Un gain achromatique inconnu** (occlusion ou exposition résiduelle) se confond avec
  μs′.
- **Une calibration caméra sur ColorChecker** laisse 2,4 ΔE00 médian sur la peau ; une
  calibration sur des spectres de peau descend à 0,4.

Les textures font 4K à 16K : l'inversion doit être rapide, déterministe et inspectable.

## Décision

1. **Espace d'entrée.** L'entrée est un RGB **linéaire, calibré, avec primaires et blanc
   connus** (espace à fixer : sRGB linéaire ou ACEScg). Le modèle direct utilisé par
   l'inversion intègre la **sensibilité spectrale de la caméra** quand elle est connue ;
   sinon, on utilise une matrice de calibration ajustée sur des spectres de peau.
2. **Paramètres estimés.** On estime le **sous-espace identifiable** :
   - charge de mélanine (Cm·d_épi) ;
   - teinte de mélanine (β), avec un a priori fort ;
   - fraction sanguine ;
   - gain achromatique (nuisance, marginalisé ou fixé si on dispose d'une AO).

   SO₂, l'épaisseur de l'épiderme, l'échelle de μs′, le carotène et la bilirubine
   viennent de l'**a priori**. Ils sont mis à jour par des mesures quand elles existent.
3. **A priori hiérarchiques explicites**, sous forme de fichiers versionnés :
   - *population* : obtenu par bayésien empirique sur la bibliothèque VFace (on inverse
     tous les visages avec un a priori large, puis on réajuste l'a priori sur la
     distribution obtenue) ;
   - *région* : via les masques VFace (lèvres, paupières, zone de barbe, oreilles) ;
   - *acteur* : ajusté sur 5 à 10 mesures ponctuelles (spectre + profil radial ; voir
     ADR-0005).
4. **Mise en œuvre sous forme de LUT 3D (type OCIO)** :
   - une grille RGB de 64³ ; chaque nœud contient la solution MAP (Gauss-Newton ou LM
     sur le modèle direct de l'ADR-0002) et la covariance a posteriori (Laplace) ;
   - application par texel en interpolation trilinéaire ;
   - une LUT par combinaison (a priori population × région), plus une variante par
     acteur si besoin.
5. **Sorties** : cartes de chromophores et **cartes d'incertitude** (écart-type par
   paramètre). Les cartes d'incertitude servent au contrôle qualité et à l'art direction.
6. **Post-traitement spatial optionnel** : régularisation qui tient compte des
   statistiques spatiales propres à la mélanine (taches nettes) et au sang (rougeurs
   diffuses). Elle opère sur les cartes, pas dans la LUT.
7. **Ombrage résiduel** : on divise l'albédo par une AO calculée depuis le déplacement
   VFace quand c'est possible ; sinon, on marginalise le gain.

## Alternatives considérées

| Alternative | Raison du rejet |
|---|---|
| Encodeur neuronal (Aliaga & Jarabo `Enc`) | A priori implicite, pas d'incertitude, demande un dataset ; peut revenir comme accélérateur si la LUT devient trop grande (entrées supplémentaires) |
| Optimisation par texel à chaque fois | Juste, mais lente sur 16K ; conservée comme **référence** pour valider la LUT |
| ICA mélanine/hémoglobine (Tsumura) | Pas de modèle physique de diffusion, pas de spectres, pas d'incertitude |
| Estimer tous les paramètres (8D) sans a priori | Mal posé ; produit du bruit dans les directions non identifiables |

## Conséquences

- **Positives**
  - Rapide (une lookup trilinéaire par texel), déterministe, intégrable dans un pipeline
    couleur existant.
  - L'a priori est inspectable et ajustable : on sait *pourquoi* un texel a telle
    translucidité.
  - L'incertitude signale les zones où une mesure ou un réglage artistique est
    nécessaire.
- **Négatives / risques**
  - La LUT dépend de l'a priori : chaque changement d'a priori demande de la recalculer
    (bon marché : 262 000 petites optimisations).
  - Une LUT ne voit pas le contexte spatial : les ambiguïtés locales ne sont levées que
    par le post-traitement spatial.
  - Qualité dépendante de la calibration amont (caméra, polarisation, exposition).

## Questions ouvertes

- VFace : espace colorimétrique de livraison, chaîne de calibration, polarisation, AO
  résiduelle ?
- Disponibilité des sensibilités spectrales de la caméra de capture ?
- Format de la LUT : `.cube` / OCIO (3 sorties par LUT, donc plusieurs LUT) ou format
  maison multicanal ?
