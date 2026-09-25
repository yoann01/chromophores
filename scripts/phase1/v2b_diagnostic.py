"""Diagnostic only (the phase-1 model is NOT changed): which global model
differences could explain the colour gap with Aliaga & Jarabo's tones?

Grid search over two global factors applied to all tones: the scale of the
bloodless-tissue (baseline) absorption and the scale of the reduced
scattering, for both readings of the melanin blend.
"""

import colour
import numpy as np
from common1 import write_report
from v2b_aliaga_tones import M_RGB, TONES, lab

from chromophores import table
from chromophores.forward import SPECULAR
from chromophores.optics import SkinModel, SkinParams
from chromophores.spectra import WAVELENGTHS

T = table.LayeredTable(method="linear")  # speed; interpolation error << the gap studied here


def tone_dE(model, scatter, beta_mode):
    out = []
    for rgb, mel, blood, thick, so2, blend in TONES:
        beta = blend if beta_mode == 0 else 1 - blend
        p = SkinParams(melanin_fraction=mel, eumelanin_ratio=beta, blood_fraction=blood, oxygen_saturation=so2,
                       epidermis_thickness_cm=thick, bilirubin_umol=0.0, scattering_scale=scatter)
        tau, ad, de, g, _ = table.dimensionless(WAVELENGTHS, p, model)
        ours = M_RGB @ ((1 - SPECULAR) * np.exp(T._A(T.coords(tau, ad, de, g))))
        out.append(float(colour.delta_E(lab(ours), lab(rgb), method="CIE 2000")))
    return np.array(out)


rows = []
for beta_mode in (0, 1):
    for b in [1.0, 0.5, 0.25, 0.1, 0.0]:
        for s in [1.0, 1.5, 2.0, 3.0]:
            dE = tone_dE(SkinModel(baseline_scale=b), s, beta_mode)
            rows.append((beta_mode, b, s, np.median(dE), dE.max()))
rows.sort(key=lambda r: r[3])
lines = ["| Lecture du blend | Échelle fond | Échelle μs′ | ΔE00 médian | ΔE00 max |", "|---|---|---|---|---|"]
for bm, b, s, med, mx in rows[:8]:
    lines.append(f"| {'β = blend' if bm == 0 else 'β = 1 − blend'} | {b:g} | {s:g} | {med:.1f} | {mx:.1f} |")
ref = [r for r in rows if r[1] == 1.0 and r[2] == 1.0]
text = f"""# V2b — diagnostic (le modèle phase 1 n'est pas modifié)

Recherche sur grille de deux facteurs globaux (tous tons confondus) : échelle de l'absorption de fond
(Jacques, 7,84·10⁸ λ⁻³·²⁵⁵) et échelle de μs′. Modèle actuel : fond 1, μs′ 1 → ΔE00 médian
{min(r[3] for r in ref):.1f}.

Meilleures combinaisons :

{chr(10).join(lines)}
"""
write_report("V2b_diagnostic.md", text)
