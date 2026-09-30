# ADR-0008 : Représentation moteur de la phase 3 — filtre épidermique sur random walk dermique, ou mélange K=2, départagés par un test d'image

- Statut : Proposé
- Date : 2026-09-30
- Décideurs : Y. Granier
- Amende (si accepté) : [ADR-0004](0004-representation-moteur-shade-before-hit.md) (points 2 à 4), [ADR-0005](0005-validation-et-criteres.md) (critère bloquant de la phase 3)

## Contexte

1. **La phase 3 est bloquée par V2-forme.** Le mélange K=2 interpolé reproduit l'albédo
   (ADR-0006) mais pas la forme du profil radial : écart médian de 0,25 au-dessus du plancher de
   bruit, contre un critère de 0,05. La piste prévue était une table « farm » (plus de photons,
   grille plus dense).
2. **ADR-0004 laisse deux risques ouverts** : le MIS spectral sur la sélection du milieu (w_k varie
   avec λ), et la géométrie mince (oreilles, narines), le mélange étant ajusté sur un derme
   semi-infini.
3. **Une autre représentation est possible** : l'épiderme (50–150 µm) comme **filtre d'interface**
   (transmittance appliquée à l'entrée et à la sortie du random walk), et le derme comme **milieu
   homogène physique** parcouru par le random walk natif du moteur. Plus de table 4D ni de
   mélange au hit, et la géométrie réelle est respectée.
4. **Étude de faisabilité** ([FILTER_FEASIBILITY](../phase3/FILTER_FEASIBILITY.md), 4 peaux ×
   3 longueurs d'onde, Monte Carlo bicouche de référence). Écart de forme au-dessus du plancher de
   bruit, en médiane / maximum :

   | Représentation | Forme | Albédo |
   |---|---|---|
   | Filtre d'interface, τ_eff ajusté sur l'albédo | **0,14** / 1,88 | exact sauf peau très claire dans le bleu (3–5 %) |
   | Mélange K=2 (table actuelle) | 0,25 / **1,93** | exact (ADR-0006) |
   | Milieu unique, albédo imposé (Masse & Walster 2023) | 0,47 / 2,31 | exact ; profil comprimé × 0,36–0,67 sur peau mate ou foncée |

   Le filtre est proche du bruit sur peau claire, surtout dans le rouge, là où la diffusion se voit.
   Il est moins bon que K=2 sur peau foncée : il ne reproduit pas le pic de diffusion simple de
   l'épiderme pigmenté. Aucune représentation ne passe V2, et la RMSE du log-profil n'est pas une
   mesure visuelle (ADR-0006 le notait déjà).

## Décision

1. **La référence ne change pas** : Monte Carlo bicouche, modèle `SkinModel` de l'ADR-0007. La
   phase 2 (inversion, LUT, cartes) n'est pas touchée.
2. **Deux représentations candidates sont implémentées dans la référence Python** :
   - **A. Mélange K=2** (ADR-0004 tel quel).
   - **B. Filtre + derme.** À chaque traversée de la surface (entrée, sortie, et chaque
     réflexion interne), le poids est multiplié par exp(−τ_eff(λ)/|cos θ|). Le derme est un milieu
     homogène aux coefficients physiques (μa derme, μs, g(λ)) : le random walk est celui du moteur,
     sans table. τ_eff est choisi pour que l'albédo soit exactement celui de la table
     multicouche : une table 3D A_filtre(τ, μa/μs′, g) plus une inversion 1D au hit, comme la
     correction d'albédo de l'ADR-0006. Si aucun τ_eff ≥ 0 ne suffit (épiderme plus clair que le
     derme), on prend τ = 0 et on corrige l'albédo simple du derme.
3. **Nouveau critère bloquant de la phase 3, « V4-bord »**, qui remplace V2-forme :
   - image spectrale (400–700 nm, 31 λ au moins, D65) d'une dalle éclairée par un demi-plan
     (bord franc) et par un spot de 1 mm, calculée à partir des profils radiaux ;
   - ΔE00 par pixel contre la référence bicouche, sur un panel de peaux tiré de l'a priori de
     population (peaux foncées incluses) et deux épaisseurs d'épiderme ;
   - seuil : ΔE00 p95 < 1 sur la bande de ±5 mm autour du bord.
4. **Règle de choix** : on retient la représentation la plus simple qui passe V4-bord, soit B, puis
   A. Si les deux passent, B. Si aucune ne passe, A avec la table « farm ». Un hybride selon τ_e
   (B pour les peaux claires, A pour les peaux foncées) n'est envisagé que si la transition peut
   être rendue continue.
5. **Masse & Walster** (milieu unique, diffusion constante) sert de base de comparaison dans
   V4-bord et de repli pour les texels hors peau.

## Alternatives considérées

| Alternative | Raison du rejet |
|---|---|
| Table « farm » pour K=2, sans autre candidat | Coûteux, et ne règle ni la géométrie mince ni le MIS de sélection ; reste la solution de secours |
| Épiderme non diffusant d'épaisseur réelle | Mesuré : trou dans le profil en champ proche (écart de forme 1,1 en médiane) |
| Filtre avec τ physique, sans ajustement | Albédo faux de ±15 à 25 % (−80 % sur peau foncée dans le bleu) |
| Milieu unique (Masse & Walster) | Profil comprimé sur peau pigmentée : la peau foncée devient trop opaque |
| Vrai random walk bicouche dans le moteur | Exige la distance à la surface le long de la normale (coque décalée), coûteux et fragile sur la géométrie mince ; à garder comme référence hors ligne si le moteur le permet |

## Conséquences

- **Positives**
  - B utilise le random walk natif : géométrie mince correcte, pas de sélection de milieu, pas
    de MIS spectral du mélange. Le coût au hit est une exponentielle par traversée de la surface.
  - Stockage par vertex inchangé (ADR-0004). L'épaisseur de l'épiderme n'intervient plus que
    par τ_eff.
  - V4-bord juge ce qui se voit, et départage aussi le chantier « farm » : on ne le lance que si
    A est retenu.
- **Négatives / risques**
  - B a une limite physique (pic de diffusion simple de l'épiderme pigmenté) qu'aucune table ne
    corrigera. Si les peaux foncées échouent à V4-bord, on revient à A ou à l'hybride.
  - τ_eff dépend de la distribution angulaire à la surface : ajusté sur une dalle, il peut dériver
    sur les zones courbes ou minces. À mesurer dans V4.
  - Deux implémentations à maintenir pendant la phase 3a.

## Questions ouvertes

- Le random walk du moteur permet-il de pondérer chaque traversée de la surface par une
  transmittance spectrale (hors Fresnel) ?
- Faut-il un τ_eff dépendant de l'angle, ou un τ_eff moyen suffit-il sur V4-bord ?
- V4-bord sur dalle suffit-il, ou faut-il ajouter une géométrie mince (oreille) dès la phase 3a ?
- La diaphonie mélanine / sang (clôture de la phase 2) reste un sujet séparé : l'empaquetage
  vasculaire ne l'améliore pas ([PACKAGING_STUDY](../phase1/PACKAGING_STUDY.md)).
