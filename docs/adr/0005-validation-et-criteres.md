# ADR-0005 : Stratégie de validation et critères d'acceptation

- Statut : Accepté
- Date : 2026-09-25
- Décideurs : Y. Granier

## Contexte

Aliaga & Jarabo valident uniquement *in silico* : leur réseau est comparé à leur propre
Monte Carlo, en incidence normale. Nous n'avons pas de banc in vivo. Il faut pourtant des
critères chiffrés et opposables pour décider du passage d'une phase à l'autre (voir la
[feuille de route](../ROADMAP.md)).

## Décision

La validation se fait en quatre niveaux. Chacun a une référence et un seuil (les valeurs
sont initiales et révisables par un ADR ultérieur).

| Niveau | Quoi | Référence | Critère d'acceptation |
|---|---|---|---|
| V1 : table directe | Interpolation de A et des moments | Monte Carlo direct sur 200 couples (peau, λ) aléatoires | Erreur relative sur A : médiane < 1 %, p95 < 3 % ; moments : p95 < 5 % |
| V2 : mélange K=2 | Profil R(r) du mélange contre la peau bicouche | Monte Carlo bicouche | RMSE du log-profil < 0,05 ; ΔE00 d'albédo < 0,5 ; pire cas (peau foncée, vert) nettement meilleur que K=1 |
| V2b : banc commun | Les 10 tons représentatifs d'Aliaga & Jarabo (Fig. 2) | Leurs RGB publiés | ΔE00 < 2 entre notre RGB simulé et le leur (écarts dus aux modèles documentés) |
| V3 : inversion | LUT 3D contre optimisation par texel ; retrouver des paramètres connus | Optimisation de référence ; tons synthétiques ; NIST (100 spectres) ; base faciale hyperspectrale (29 visages) | Écart LUT/optimisation < 1 % sur les paramètres identifiables ; ΔE00 de reconstruction < 1 ; les spectres reconstruits restent dans l'enveloppe NIST |
| V4 : moteur | Rendu du mélange K=2 contre un random walk réellement multicouche **dans le moteur** | Même moteur, géométrie bicouche explicite | ΔE00 < 1 sous éclairage frontal ; écart documenté en incidence rasante et en rétroéclairage (oreilles) ; surcoût de variance mesuré |
| V5 : mesures (quand disponibles) | Spectre et profil radial ponctuels sur un acteur test | Spectrophotomètre + spot laser/LED en polarisation croisée | Profils mesurés et prédits cohérents après calibration de μs′ ; ΔE00 < 2 sur les spectres |

Règles :
- Chaque résultat est **scripté et reproductible** (seed, nombre de photons, version de la
  table), et les figures sont générées dans `docs/figures/`.
- Un seuil non atteint **bloque** la phase suivante, sauf décision explicite consignée
  (nouvel ADR ou note dans la feuille de route).

## Alternatives considérées

- *Validation uniquement in silico* (comme Aliaga & Jarabo) : insuffisante pour la
  production. Les niveaux V4 (dans notre moteur) et V5 (mesures) sont ajoutés.
- *Validation perceptuelle seule* (revues artistiques) : utile, mais pas opposable. Elle
  est conservée en complément de V4.

## Conséquences

- Il faut un Monte Carlo de référence **indépendant** de la table, et un rendu
  multicouche explicite dans le moteur.
- V5 demande un peu de matériel et l'accès à un acteur. Le niveau V5 est optionnel pour
  la phase 1, mais requis avant une mise en production.

## Questions ouvertes

- Seuils ΔE : quels sont ceux déjà utilisés en production pour la validation couleur des
  assets ?
- Existe-t-il en interne des mesures spectrales ou de diffusion réutilisables ?
