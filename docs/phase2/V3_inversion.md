# V3 : inversion RGB → chromophores, sur spectres mesurés ISSA

- 300 spectres ISSA **non utilisés** pour construire l'a priori ; 300 dans le gamut de la LUT.
- Entrée : albédo en polarisation croisée (SCI − R_sp) rendu en sRGB linéaire D65. LUT 33³, a priori « population ».
- Bruit supposé sur le RGB : 0,004 + 2 % (calibration).

## (a) LUT vs optimisation MAP directe

| | médiane | p95 |
|---|---|---|
| erreur relative mélanine | 2.80 % | 13.43 % |
| erreur relative sang | 4.32 % | 13.71 % |
| écart / σ a posteriori (tous paramètres) | 0.048 | 0.223 |

Critère ADR-0005 (écart LUT / optimisation < 1 % sur les paramètres identifiables) → **NON ATTEINT**.

## (b) Ce que le RGB permet de retrouver (vs ajustement du spectre complet 400–700 nm)

| Paramètre | corrélation (log) RGB vs spectre | erreur relative médiane | |z| médian | % |z| < 2 | σ LUT médian / σ spectre médian |
|---|---|---|---|---|---|
| melanin_fraction | 0.84 | 41 % | 0.67 | 98 % | 0.59 / 0.23 |
| eumelanin_ratio | 0.49 | 81 % | 0.68 | 100 % | 1.85 / 0.96 |
| blood_fraction | 0.78 | 46 % | 0.81 | 92 % | 0.40 / 0.08 |
| oxygen_saturation | 0.60 | 12 % | 0.49 | 100 % | 1.57 / 0.63 |
| epidermis_thickness_cm | 0.53 | 35 % | 0.75 | 99 % | 1.04 / 0.33 |
| scattering_scale | 0.08 | 50 % | 1.43 | 70 % | 0.40 / 0.27 |

z = (x_RGB − x_spectre) / √(σ²_RGB + σ²_spectre) dans l'espace non borné : |z| < 2 pour ~95 % des cas si les incertitudes sont bien calibrées.

## (c) Reconstruction de la couleur

ΔE00 (D65) entre l'albédo d'entrée et l'albédo reconstruit : médiane 0.56, p95 1.64, max 4.33.
Critère ADR-0005 (ΔE00 < 1) → **NON ATTEINT (p95)**.

## (d) Spectre reconstruit depuis le RGB vs spectre mesuré

| Méthode RGB → spectre | RMS spectral médian (400–700) | p95 | ΔE00 sous D65 méd. / p95 | ΔE00 sous A méd. / p95 | ΔE00 sous FL11 méd. / p95 | ΔE00 sous LED-B3 méd. / p95 |
|---|---|---|---|---|---|---|
| biophysique (notre LUT) | 0.0097 | 0.0217 | 0.56 / 1.64 | 0.60 / 1.63 | 0.76 / 2.05 | 0.56 / 1.55 |
| Jakob & Hanika 2019 | 0.0345 | 0.0702 | 0.03 / 0.05 | 1.42 / 2.43 | 2.68 / 4.16 | 1.36 / 2.13 |
| Mallett & Yuksel 2019 | 0.0650 | 0.0840 | 0.10 / 0.13 | 1.38 / 2.39 | 1.06 / 2.07 | 0.74 / 1.25 |
