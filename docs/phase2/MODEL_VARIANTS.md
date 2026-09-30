# Variantes du modèle direct contre les mesures (phéomélanine corrigée)

Script : `scripts/phase2/model_variants_check.py`. Ajustement libre (a priori très large) de 298 spectres
ISSA (400–700 nm) et des 100 spectres NIST (400–1000 nm) ; médianes des paramètres ajustés.

| Modèle | Données | n | RMS médian | RMS p95 | Mélanine | Sang | Épiderme (mm) | Échelle μs′ |
|---|---|---|---|---|---|---|---|---|
| fond ×0,5, eau (actuel) | ISSA | 298 | 0.0045 | 0.0100 | 0.021 | 0.0128 | 0.192 | 0.98 |
| fond ×0,5, eau (actuel) | NIST | 100 | 0.0132 | 0.0186 | 0.005 | 0.0165 | 0.197 | 0.96 |
| fond ×0,5, sans eau | ISSA | 298 | 0.0046 | 0.0102 | 0.021 | 0.0127 | 0.193 | 0.96 |
| fond ×0,5, sans eau | NIST | 100 | 0.0275 | 0.0335 | 0.010 | 0.1089 | 0.192 | 2.17 |
| fond ×1, eau | ISSA | 298 | 0.0044 | 0.0106 | 0.030 | 0.0225 | 0.115 | 1.61 |
| fond ×1, eau | NIST | 100 | 0.0179 | 0.0242 | 0.007 | 0.0301 | 0.112 | 1.62 |
| fond ×1, sans eau (Opus) | ISSA | 298 | 0.0045 | 0.0105 | 0.029 | 0.0225 | 0.116 | 1.58 |
| fond ×1, sans eau (Opus) | NIST | 100 | 0.0284 | 0.0353 | 0.013 | 0.1323 | 0.120 | 2.99 |
