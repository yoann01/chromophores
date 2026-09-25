"""Compare our phase-1 forward model with the BioSkin (Aliaga et al. 2023) decoder.

1. Can the 2023 decoder reproduce the RGB published for the 10 tones of
   Aliaga & Jarabo 2026 (Fig. 2)? Tested with their own RGB conversion
   (CIE CMFs, equal-energy, custom XYZ->RGB matrix, BGR order) and with ours
   (D65, linear sRGB).
2. Spectra: BioSkin decoder vs our table, tone by tone.
3. A nearly chromophore-free skin, to isolate baseline absorption and scattering.
"""

import os
import sys

import colour
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "phase1"))
from bioskin_numpy import BioSkinDecoder  # noqa: E402
from common import SERIES, plt, save  # noqa: E402
from common1 import write_report  # noqa: E402
from v2b_aliaga_tones import TONES, lab  # noqa: E402

from chromophores import color  # noqa: E402
from chromophores.forward import SkinForward  # noqa: E402
from chromophores.optics import SkinModel, SkinParams  # noqa: E402

BIOSKIN = os.environ.get("BIOSKIN_DIR", "/home/user/facebookresearch/bioskin")
D = BioSkinDecoder(os.path.join(BIOSKIN, "pretrained_models", "BioSkin"))
LAM = D.lam.astype(float)
VIS = LAM <= 780
F = SkinForward()

# Their RGB conversion (bioskin/spectrum/color_spectrum.py): CMFs on 380-780 (resampled), sum / n,
# custom "D65" matrix, output in BGR order.
M_THEIRS = np.array([[12.420749767575597, -3.7146545623403702, 0.21326517187700963],
                     [-5.8918942422391520, 7.1897002408297634, -0.78197265064351096],
                     [-1.9108846122515208, 0.15926100824153880, 4.0520420060904261]])
CMFS = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]
XYZBAR = np.stack([np.interp(LAM[VIS], CMFS.wavelengths, CMFS.values[:, k]) for k in range(3)])


def rgb_theirs(R):
    xyz = XYZBAR @ R[VIS] / VIS.sum()
    r, g, b = (M_THEIRS.T @ xyz)  # transposed matrix as in create_matrix_xyz_to_lrgb
    return np.array([b, g, r]), np.array([r, g, b])  # (as stored = BGR, re-ordered RGB)


M_OURS = color.reflectance_to_linear_srgb_matrix(LAM[VIS])


def rgb_ours(R):
    return M_OURS @ R[VIS]


def de(a, b):
    return float(colour.delta_E(lab(a), lab(b), method="CIE 2000"))


rows, spectra = [], []
for i, (rgb, mel, blood, thick, so2, blend) in enumerate(TONES):
    Rb = D.decode(D.unwarp(np.array([mel]), np.array([blood]), np.array([thick]), np.array([so2]), np.array([blend])))[0]
    p = SkinParams(melanin_fraction=mel, eumelanin_ratio=blend, blood_fraction=blood, oxygen_saturation=so2,
                   epidermis_thickness_cm=thick, bilirubin_umol=0.0)
    Ro = F.albedo(LAM, p)
    bgr_t, rgb_t = rgb_theirs(Rb)
    rows.append((i + 1, rgb, de(rgb_t, rgb), de(bgr_t, rgb), de(rgb_ours(Rb), rgb), de(rgb_ours(Ro), rgb), de(rgb_ours(Ro), rgb_ours(Rb)), rgb_t))
    spectra.append((Rb, Ro))

lines = ["| Ton | ΔE00 décodeur 2023 + conversion BioSkin (RGB) | idem (BGR) | décodeur 2023 + D65/sRGB | notre modèle + D65/sRGB | notre modèle vs décodeur 2023 (D65) |",
         "|---|---|---|---|---|---|"]
for i, rgb, a, b, c, d_, e_, _ in rows:
    lines.append(f"| #{i:02d} | {a:.1f} | {b:.1f} | {c:.1f} | {d_:.1f} | {e_:.1f} |")
med = np.median(np.array([r[2:7] for r in rows]), 0)
lines.append(f"| **médiane** | {med[0]:.1f} | {med[1]:.1f} | {med[2]:.1f} | {med[3]:.1f} | {med[4]:.1f} |")

# Near chromophore-free skin: baseline + scattering only.
Rb0 = D.decode(D.unwarp(np.array([0.001]), np.array([0.001]), np.array([0.01]), np.array([0.75]), np.array([0.5])))[0]
p0 = SkinParams(melanin_fraction=0.001, blood_fraction=0.001, epidermis_thickness_cm=0.01, carotene_umol=0.0, bilirubin_umol=0.0)
Ro0 = F.albedo(LAM, p0)
Ro0_nb = F_nb = SkinForward(model=SkinModel(baseline_scale=0.0)).albedo(LAM, p0)
probe = [450, 550, 650, 750, 900]
idx = [int(np.argmin(np.abs(LAM - l))) for l in probe]
lines0 = ["| λ (nm) | " + " | ".join(str(p) for p in probe) + " |", "|---|" + "---|" * len(probe),
          "| BioSkin 2023 | " + " | ".join(f"{Rb0[k]:.3f}" for k in idx) + " |",
          "| notre modèle | " + " | ".join(f"{Ro0[k]:.3f}" for k in idx) + " |",
          "| notre modèle sans absorption de fond | " + " | ".join(f"{Ro0_nb[k]:.3f}" for k in idx) + " |"]

text = f"""# Comparaison avec BioSkin (Aliaga et al. 2023)

Décodeur pré-entraîné `BioSkin.pt` (p → R(λ), 380–1000 nm), exécuté en numpy
(`scripts/external/bioskin_numpy.py`), avec les paramètres publiés des 10 tons d'Aliaga & Jarabo 2026.
« Mel. Blend » est lu comme la fraction d'eumélanine.

## 1. Reproduction des RGB publiés (ΔE00 vs RGB de la Fig. 2 d'Aliaga & Jarabo 2026)

{chr(10).join(lines)}

## 2. Peau quasi sans chromophores (mélanine 0,1 %, sang 0,1 %, épiderme 100 µm)

{chr(10).join(lines0)}
"""
write_report("BIOSKIN_COMPARISON.md", text)

fig, axes = plt.subplots(2, 5, figsize=(15, 5.6), sharex=True)
for ax, (i, *_), (Rb, Ro) in zip(axes.ravel(), rows, spectra):
    ax.plot(LAM, Rb, color=SERIES[0], label="BioSkin 2023")
    ax.plot(LAM, Ro, color=SERIES[1], label="notre modèle")
    ax.set_title(f"Ton #{i:02d}", loc="left", fontsize=9)
    ax.set_ylim(0, 0.8)
axes[0, 0].legend(fontsize=7)
for ax in axes[1]:
    ax.set_xlabel("λ (nm)")
for ax in axes[:, 0]:
    ax.set_ylabel("Réflectance diffuse")
save(fig, "bioskin_vs_ours_spectra.png")
