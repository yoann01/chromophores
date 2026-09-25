# V2 : mélange K=2 vs peau bicouche

## (a) Qualité du fit, sur toutes les cellules de la table (4648 cellules ajustées, 56 trop sombres pour un profil)

| | médiane | p95 | max |
|---|---|---|---|
| RMSE log-profil K=1 (non pondérée) | 1.435 | 5.843 | 12.298 |
| RMSE log-profil K=2 (non pondérée) | 0.431 | 1.286 | 7.602 |
| RMSE log-profil K=1 (pondérée énergie) | 0.303 | 2.005 | 21.049 |
| RMSE log-profil K=2 (pondérée énergie) | 0.124 | 0.560 | 14.954 |
| erreur relative A, K=1 | 0.011 | 0.113 | 0.989 |
| erreur relative A, K=2 | 0.001 | 0.006 | 0.314 |

Les valeurs incluent le bruit Monte Carlo de la table (3·10⁴ photons par cellule).

## (b) Contrôle indépendant : 60 couples (peau, λ ∈ 380–780 nm), mélange **interpolé** vs Monte Carlo direct

| | médiane | p95 | max |
|---|---|---|---|
| RMSE log-profil K=1 (non pondérée) | 3.119 | 5.009 | 6.469 |
| RMSE log-profil K=2 (non pondérée) | 1.089 | 3.070 | 4.514 |
| plancher de bruit MC (non pondéré) | 0.183 | 0.264 | 0.359 |
| RMSE log-profil K=1 (pondérée énergie) | 0.387 | 1.386 | 1.762 |
| RMSE log-profil K=2 (pondérée énergie) | 0.150 | 0.622 | 3.644 |
| plancher de bruit MC (pondéré énergie) | 0.036 | 0.084 | 0.134 |
| erreur relative A, K=1 | 0.005 | 0.024 | 0.038 |
| erreur relative A, K=2 | 0.005 | 0.024 | 0.038 |

Référence : Monte Carlo direct à 10⁵ photons. Le plancher de bruit est l'écart entre deux tirages indépendants de la référence.

## (c) Couleur de l'albédo (12 peaux, spectre 380–780 nm, D65)

| | médiane | p95 | max |
|---|---|---|---|
| ΔE00 albédo du mélange K=2 (avec correction d'albédo) vs MC direct | 0.237 | 0.403 | 0.436 |
| ΔE00 albédo du mélange K=2 (sans correction) vs MC direct | 0.825 | 2.104 | 2.355 |
| ΔE00 albédo de la table vs MC direct | 0.237 | 0.403 | 0.436 |

Critères ADR-0005 : RMSE log-profil < 0,05 (pondérée énergie, p95, contrôle indépendant) → **NON ATTEINT** ;
écart au plancher de bruit < 0,05 (p95) → **NON ATTEINT** ;
ΔE00 albédo < 0,5 → **OK**.
