# Absorption de fond : arbitrage par les spectres mesurés ISSA

- Données : ISSA (Leeds Skin Database), spectrophotomètre SCI, 400–700 nm / 10 nm.
- Échantillon stratifié : 960 spectres (60 visage + 60 corps max. par groupe ethnique).
- Modèle de mesure : SCI = R_sp (Fresnel, n = 1,4) + (1 − R_sp)·A_table. Carotène 0,4 µM, bilirubine 0.
- 6 paramètres ajustés par spectre (mélanine, β, sang, SO₂, épaisseur, échelle de μs′ ∈ [0,3 ; 3]), multi-départs.
- Plage « plausible » de l'échelle de μs′ (Jacques 2013, dispersion inter-individuelle) : [0,6 ; 1,6].

| Échelle du fond | RMS résidu médian | p95 | % RMS < 0,01 | μs′ fitté médian | % μs′ hors [0,6 ; 1,6] | % μs′ en butée |
|---|---|---|---|---|---|---|
| 1 | 0.0047 | 0.0093 | 96 % | 1.49 | 56 % | 17 % |
| 0.5 | 0.0046 | 0.0090 | 97 % | 0.95 | 53 % | 14 % |
| 0.25 | 0.0049 | 0.0115 | 91 % | 0.70 | 58 % | 17 % |
| 0 | 0.0077 | 0.0185 | 65 % | 2.63 | 76 % | 43 % |

## Par groupe

| Groupe | n | RMS méd. (fond 1) | RMS méd. (fond 0.5) | RMS méd. (fond 0.25) | RMS méd. (fond 0) | μs′ méd. (fond 1) | μs′ méd. (fond 0.5) | μs′ méd. (fond 0.25) | μs′ méd. (fond 0) |
|---|---|---|---|---|---|---|---|---|---|
| AB | 120 | 0.0051 | 0.0044 | 0.0045 | 0.0051 | 2.41 | 1.81 | 1.22 | 2.49 |
| AF | 120 | 0.0028 | 0.0022 | 0.0020 | 0.0017 | 1.10 | 0.88 | 0.82 | 0.91 |
| CA | 120 | 0.0059 | 0.0060 | 0.0074 | 0.0129 | 1.31 | 0.73 | 0.57 | 2.89 |
| CN | 120 | 0.0047 | 0.0048 | 0.0053 | 0.0099 | 1.35 | 0.77 | 0.59 | 2.96 |
| IQ | 120 | 0.0056 | 0.0062 | 0.0079 | 0.0120 | 0.87 | 0.52 | 0.51 | 2.79 |
| JP | 120 | 0.0055 | 0.0054 | 0.0054 | 0.0071 | 2.49 | 1.41 | 0.93 | 2.61 |
| SA | 120 | 0.0034 | 0.0034 | 0.0039 | 0.0068 | 1.28 | 0.87 | 0.58 | 2.87 |
| TH | 120 | 0.0041 | 0.0041 | 0.0045 | 0.0072 | 1.31 | 1.01 | 0.65 | 2.98 |

Codes : CA caucasien, CN chinois, SA sud-asiatique, AF africain, IQ irakien, TH thaï, JP japonais, AB arabe.
