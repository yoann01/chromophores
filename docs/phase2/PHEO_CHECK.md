# Vérification de la constante de phéomélanine

Donner & Jensen 2006, éq. 3 et 4 : eumélanine 6,6·10¹⁰·λ⁻³·³³ mm⁻¹, phéomélanine 2,9·10¹⁴·λ⁻⁴·⁷⁵ mm⁻¹
(λ en nm), soit 6,6·10¹¹ et **2,9·10¹⁵** cm⁻¹. `spectra.py` avait 2,9·10¹⁴ cm⁻¹ pour la phéomélanine
(10 fois trop bas). Script : `scripts/phase2/pheo_check.py` (médianes des paramètres ajustés).

| Phéomélanine | Fond | Données | n | RMS médian | RMS p95 | Mélanine Cm | β (eu) | Sang | SO₂ | Épiderme (mm) | Échelle μs′ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ×1 | ×0.5 | ISSA | 298 | 0.0045 | 0.0100 | 0.106 | 0.09 | 0.0140 | 0.76 | 0.184 | 1.01 |
| ×1 | ×0.5 | NIST | 100 | 0.0131 | 0.0186 | 0.050 | 0.00 | 0.0164 | 0.90 | 0.197 | 0.96 |
| ×10 | ×0.5 | ISSA | 298 | 0.0045 | 0.0100 | 0.021 | 0.30 | 0.0128 | 0.74 | 0.192 | 0.98 |
| ×10 | ×0.5 | NIST | 100 | 0.0132 | 0.0186 | 0.005 | 0.00 | 0.0165 | 0.90 | 0.197 | 0.96 |
| ×1 | ×1 | ISSA | 298 | 0.0045 | 0.0108 | 0.209 | 0.05 | 0.0236 | 0.77 | 0.106 | 1.75 |
| ×1 | ×1 | NIST | 100 | 0.0180 | 0.0251 | 0.071 | 0.00 | 0.0293 | 0.98 | 0.112 | 1.60 |
| ×10 | ×1 | ISSA | 298 | 0.0044 | 0.0106 | 0.030 | 0.12 | 0.0225 | 0.74 | 0.115 | 1.61 |
| ×10 | ×1 | NIST | 100 | 0.0179 | 0.0242 | 0.007 | 0.00 | 0.0301 | 0.99 | 0.112 | 1.62 |

## Lecture

1. **Le papier tranche** : Opus est fidèle à Donner & Jensen 2006. L'erreur était dans
   `spectra.py` (phéomélanine recopiée en cm⁻¹ sans conversion), corrigée le 2026-09-30.
2. **Les mesures ne départagent pas les deux constantes par la qualité d'ajustement** (RMS
   identique) : l'ajustement compense par Cm et β. C'est pour cela que l'erreur est passée
   inaperçue.
3. **Mais les paramètres ajustés deviennent plausibles** avec la bonne constante : la mélanine
   médiane passe de 10,6 % à 2,1 % (ISSA) et de 5 % à 0,5 % (NIST), dans la plage de Donner &
   Jensen et proche du défaut d'Opus (Vm 0,005). β remonte de 0,09 à 0,30 sur ISSA, mais reste à 0
   sur NIST : il n'est toujours pas interprétable comme un ratio biologique.
4. **La calibration du fond de l'ADR-0007 tient** : avec la phéomélanine corrigée, ×0,5 et ×1
   ajustent ISSA aussi bien, mais ×1 exige μs′ × 1,6 et ajuste moins bien NIST (RMS 0,018 contre
   0,013). Sang, SO₂, épaisseur et diffusion ne bougent presque pas.
5. **Conséquence** : les cartes de mélanine et de ratio eu/phéo, les a priori et les LUT dépendent
   de la constante. Ils sont recalculés (ajustements de population, a priori, LUT 33³ et 65³,
   covariance, clôture). Les rapports antérieurs (V3, VFace, Emily, benchmark BioSkin) datent
   d'avant la correction ; leurs conclusions sur le sang et la diffusion restent valables.
