# Comparaison avec BioSkin (Aliaga et al. 2023)

Décodeur pré-entraîné `BioSkin.pt` (p → R(λ), 380–1000 nm), exécuté en numpy
(`scripts/external/bioskin_numpy.py`), avec les paramètres publiés des 10 tons d'Aliaga & Jarabo 2026.
« Mel. Blend » est lu comme la fraction d'eumélanine.

## 1. Reproduction des RGB publiés (ΔE00 vs RGB de la Fig. 2 d'Aliaga & Jarabo 2026)

| Ton | ΔE00 décodeur 2023 + conversion BioSkin (RGB) | idem (BGR) | décodeur 2023 + D65/sRGB | notre modèle + D65/sRGB | notre modèle vs décodeur 2023 (D65) |
|---|---|---|---|---|---|
| #01 | 3.9 | 27.0 | 2.1 | 10.4 | 11.4 |
| #02 | 5.5 | 12.0 | 9.6 | 6.6 | 6.9 |
| #03 | 0.9 | 34.0 | 2.1 | 2.8 | 2.5 |
| #04 | 2.6 | 21.2 | 2.9 | 8.8 | 8.9 |
| #05 | 0.5 | 37.2 | 4.3 | 10.6 | 8.7 |
| #06 | 1.4 | 34.0 | 5.3 | 8.2 | 5.1 |
| #07 | 3.2 | 30.4 | 7.4 | 11.7 | 12.4 |
| #08 | 2.7 | 39.5 | 3.6 | 9.9 | 7.6 |
| #09 | 5.2 | 29.9 | 6.7 | 11.5 | 9.3 |
| #10 | 4.6 | 36.8 | 4.2 | 11.3 | 9.0 |
| **médiane** | 2.9 | 32.2 | 4.2 | 10.2 | 8.8 |

## 2. Peau quasi sans chromophores (mélanine 0,1 %, sang 0,1 %, épiderme 100 µm)

| λ (nm) | 450 | 550 | 650 | 750 | 900 |
|---|---|---|---|---|---|
| BioSkin 2023 | 0.615 | 0.642 | 0.757 | 0.767 | 0.765 |
| notre modèle | 0.383 | 0.440 | 0.518 | 0.553 | 0.581 |
| notre modèle sans absorption de fond | 0.631 | 0.661 | 0.820 | 0.824 | 0.756 |
