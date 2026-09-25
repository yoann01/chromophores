# Architecture Decision Records

Format : Nygard / MADR allégé. On écrit un ADR par décision structurante, et on ne le
modifie plus une fois accepté : on le remplace par un nouvel ADR (statut « Remplacé
par »).

| # | Titre | Statut | Date |
|---|---|---|---|
| [0001](0001-separer-direct-et-inverse.md) | Séparer le modèle direct (bien posé, tabulé) de l'inversion (mal posée, probabiliste) | Accepté | 2026-09-25 |
| [0002](0002-modele-direct-table-sans-dimension.md) | Modèle direct : peau bicouche → table sans dimension → mélange K=2 par cellule | Accepté | 2026-09-25 |
| [0003](0003-inversion-rgb-lut3d-a-priori.md) | Inversion RGB → chromophores par LUT 3D avec a priori explicite et incertitude | Accepté | 2026-09-25 |
| [0004](0004-representation-moteur-shade-before-hit.md) | Représentation moteur : chromophores par vertex, évaluation au hit, mélange K=2 en hero wavelength | Accepté | 2026-09-25 |
| [0005](0005-validation-et-criteres.md) | Stratégie de validation et critères d'acceptation | Accepté | 2026-09-25 |
| [0006](0006-ajustements-phase1.md) | Ajustements phase 1 : mélange à g partagé, interpolation préservant l'albédo, critère V2 relatif au bruit (amende 0002 et 0005) | Accepté | 2026-09-25 |

La mise en œuvre par phases est décrite dans [`../ROADMAP.md`](../ROADMAP.md).

## Modèle

```
# ADR-XXXX : <titre>
- Statut : Proposé | Accepté | Rejeté | Remplacé par ADR-YYYY
- Date :
- Décideurs :
## Contexte
## Décision
## Alternatives considérées
## Conséquences (positives / négatives / risques)
## Questions ouvertes
```
