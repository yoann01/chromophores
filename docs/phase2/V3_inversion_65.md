# V3 : inversion RGB → chromophores, sur spectres mesurés ISSA

- 300 spectres ISSA **non utilisés** pour construire l'a priori ; 300 dans le gamut de la LUT.
- Entrée : albédo en polarisation croisée (SCI − R_sp) rendu en sRGB linéaire D65. LUT 33³, a priori « population ».
- Bruit supposé sur le RGB : 0,004 + 2 % (calibration).

## (a) LUT vs optimisation MAP directe

| | médiane | p95 |
|---|---|---|
| erreur relative mélanine | 0.97 % | 8.12 % |
| erreur relative sang | 1.36 % | 5.66 % |
| écart / σ a posteriori (tous paramètres) | 0.016 | 0.128 |

Critère ADR-0005 (écart LUT / optimisation < 1 % sur les paramètres identifiables) → **NON ATTEINT**.

## (b) Ce que le RGB permet de retrouver (vs ajustement du spectre complet 400–700 nm)

| Paramètre | corrélation (log) RGB vs spectre | erreur relative médiane | |z| médian | % |z| < 2 | σ LUT médian / σ spectre médian |
|---|---|---|---|---|---|
| melanin_fraction | 0.84 | 38 % | 0.67 | 97 % | 0.58 / 0.23 |
| eumelanin_ratio | 0.49 | 80 % | 0.68 | 100 % | 1.82 / 0.96 |
| blood_fraction | 0.79 | 46 % | 0.84 | 92 % | 0.38 / 0.08 |
| oxygen_saturation | 0.60 | 11 % | 0.48 | 100 % | 1.57 / 0.63 |
| epidermis_thickness_cm | 0.51 | 35 % | 0.75 | 99 % | 1.03 / 0.33 |
| scattering_scale | 0.09 | 49 % | 1.43 | 70 % | 0.40 / 0.27 |

z = (x_RGB − x_spectre) / √(σ²_RGB + σ²_spectre) dans l'espace non borné : |z| < 2 pour ~95 % des cas si les incertitudes sont bien calibrées.

## (c) Reconstruction de la couleur

ΔE00 (D65) entre l'albédo d'entrée et l'albédo reconstruit : médiane 0.39, p95 1.44, max 4.35.
Critère ADR-0005 (ΔE00 < 1) → **NON ATTEINT (p95)**.

## (d) Spectre reconstruit depuis le RGB vs spectre mesuré

| Méthode RGB → spectre | RMS spectral médian (400–700) | p95 | ΔE00 sous D65 méd. / p95 | ΔE00 sous A méd. / p95 | ΔE00 sous FL11 méd. / p95 | ΔE00 sous LED-B3 méd. / p95 |
|---|---|---|---|---|---|---|
| biophysique (notre LUT) | 0.0094 | 0.0209 | 0.39 / 1.44 | 0.51 / 1.45 | 0.71 / 1.98 | 0.43 / 1.40 |
| biophysique + correction colorimétrique | 0.0094 | 0.0212 | 0.00 / 0.00 | 0.36 / 1.14 | 0.66 / 1.84 | 0.28 / 0.79 |
| Jakob & Hanika 2019 | 0.0345 | 0.0702 | 0.03 / 0.05 | 1.42 / 2.43 | 2.68 / 4.16 | 1.36 / 2.13 |
| Mallett & Yuksel 2019 | 0.0650 | 0.0840 | 0.10 / 0.13 | 1.38 / 2.39 | 1.06 / 2.07 | 0.74 / 1.25 |
