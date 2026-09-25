"""V2b (ADR-0005): the 10 representative tones of Aliaga & Jarabo 2026 (Fig. 2).

Their published biophysical parameters are run through our phase-1 model and
the resulting D65 linear-sRGB albedo is compared with their published RGB
triplets (assumed linear sRGB, D65).

Unknowns in the mapping, reported rather than tuned:
- "Mel. Blend" is read either as the eumelanin fraction (beta = blend) or
  as the pheomelanin fraction (beta = 1 - blend);
- epidermal thickness is assumed in cm (0.0023-0.0341 cm = 23-341 um);
- their basal bilirubin / carotene / water levels and exact scattering are
  unknown; we use ours (scattering scale 1).
"""

import colour
import numpy as np
from common1 import SERIES, plt, save, write_report

from chromophores import color
from chromophores.forward import SkinForward
from chromophores.optics import SkinParams
from chromophores.spectra import WAVELENGTHS

# (rgb, melanin VF, blood VF, epi thickness, SO2, mel blend) from Fig. 2.
TONES = [
    ((0.069, 0.024, 0.013), 0.3963, 0.0241, 0.0171, 0.8982, 0.0970),
    ((0.034, 0.020, 0.014), 0.3257, 0.0355, 0.0334, 0.6706, 0.7814),
    ((0.170, 0.045, 0.021), 0.2290, 0.0552, 0.0115, 0.7770, 0.2421),
    ((0.111, 0.059, 0.042), 0.1049, 0.0590, 0.0341, 0.5594, 0.7935),
    ((0.409, 0.145, 0.074), 0.1212, 0.0334, 0.0041, 0.6389, 0.8214),
    ((0.373, 0.166, 0.089), 0.0512, 0.0079, 0.0245, 0.6023, 0.3877),
    ((0.632, 0.286, 0.219), 0.0199, 0.0341, 0.0023, 0.5469, 0.1447),
    ((0.472, 0.212, 0.084), 0.0997, 0.0057, 0.0065, 0.6632, 0.6112),
    ((0.783, 0.460, 0.344), 0.0058, 0.0077, 0.0144, 0.7498, 0.4623),
    ((0.626, 0.346, 0.183), 0.0287, 0.0022, 0.0137, 0.9463, 0.6038),
]

F = SkinForward()
M_RGB = color.reflectance_to_linear_srgb_matrix(WAVELENGTHS)
SRGB = colour.RGB_COLOURSPACES["sRGB"]


def lab(rgb_lin):
    xyz = colour.RGB_to_XYZ(np.asarray(rgb_lin), SRGB, apply_cctf_decoding=False)
    return colour.XYZ_to_Lab(xyz, SRGB.whitepoint)


rows = []
for i, (rgb, mel, blood, thick, so2, blend) in enumerate(TONES):
    res = {}
    for name, beta in (("β = blend", blend), ("β = 1 − blend", 1 - blend)):
        p = SkinParams(melanin_fraction=mel, eumelanin_ratio=beta, blood_fraction=blood, oxygen_saturation=so2,
                       epidermis_thickness_cm=thick, bilirubin_umol=0.0)
        ours = M_RGB @ F.albedo(WAVELENGTHS, p)
        res[name] = (ours, float(colour.delta_E(lab(ours), lab(rgb), method="CIE 2000")))
    rows.append((i + 1, rgb, res))

keys = ["β = blend", "β = 1 − blend"]
lines = ["| Ton | RGB publié | " + " | ".join(f"RGB nous ({k}) | ΔE00" for k in keys) + " |", "|---|---|" + "---|---|" * len(keys)]
for i, rgb, res in rows:
    cells = " | ".join(f"({res[k][0][0]:.3f}, {res[k][0][1]:.3f}, {res[k][0][2]:.3f}) | {res[k][1]:.1f}" for k in keys)
    lines.append(f"| #{i:02d} | ({rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f}) | {cells} |")
dE = {k: np.array([r[2][k][1] for r in rows]) for k in keys}
best = min(keys, key=lambda k: np.median(dE[k]))
summary = " ; ".join(f"{k} : médiane {np.median(dE[k]):.1f}, max {dE[k].max():.1f}" for k in keys)
text = f"""# V2b : les 10 tons représentatifs d'Aliaga & Jarabo (Fig. 2)

Paramètres publiés → notre modèle phase 1 (table multicouche) → sRGB linéaire D65, comparés aux RGB publiés.

{chr(10).join(lines)}

ΔE00 : {summary}.
Critère ADR-0005 : ΔE00 < 2 → **{"OK" if dE[best].max() < 2 else "NON ATTEINT"}** (meilleure lecture : {best}).
"""
write_report("V2b_aliaga_tones.md", text)

fig, ax = plt.subplots(figsize=(9, 2.6))
for i, rgb, res in rows:
    ax.add_patch(plt.Rectangle((i - 1, 1), 0.95, 0.95, color=color.encode_srgb(np.array(rgb))))
    ax.add_patch(plt.Rectangle((i - 1, 0), 0.95, 0.95, color=color.encode_srgb(res[best][0])))
    ax.text(i - 0.52, -0.25, f"#{i:02d}\nΔE {res[best][1]:.1f}", ha="center", va="top", fontsize=7)
ax.text(-0.1, 1.47, "publié", ha="right", va="center", fontsize=8)
ax.text(-0.1, 0.47, "notre modèle", ha="right", va="center", fontsize=8)
ax.set_xlim(-1.4, 10)
ax.set_ylim(-0.9, 2)
ax.axis("off")
ax.set_title(f"Tons d'Aliaga & Jarabo : publié vs notre modèle ({best})", loc="left", fontsize=10)
save(fig, "phase1_v2b_tones.png")
