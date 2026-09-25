"""Production-oriented colour checks for RGB -> chromophore inversion.

1. Camera metamerism: a calibrated albedo (e.g. VFace, cross-polarised
   photogrammetry) is produced by a real camera, colour-corrected with a
   3x3 matrix fitted on a ColorChecker. How far is it from the CIE albedo
   the inversion assumes, on skin spectra?
2. Illuminant metamerism: skins that share the reference RGB albedo under
   D65 (i.e. equally valid answers of an RGB inversion) are rendered under
   set lights (tungsten, fluorescent, phosphor LED, RGB LED). The spread is
   the colour error the *prior* of the inversion is responsible for.
"""

import colour
import numpy as np
from common import SERIES, plt, save
from identifiability import RANGES, from_unit, to_unit

from chromophores import color, reflectance
from chromophores.optics import SkinParams
from chromophores.spectra import WAVELENGTHS as LAM

ILLUMINANTS = ["D65", "A", "FL11", "LED-B3", "LED-RGB1"]
ILL_LABELS = {"D65": "D65 (capture)", "A": "A (tungstène)", "FL11": "FL11 (fluo tribande)", "LED-B3": "LED-B3 (LED blanche)", "LED-RGB1": "LED-RGB1 (LED RGB)"}
CMFS = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]


def resample(sd_values, wl):
    return np.interp(LAM, wl, sd_values)


def xyz_matrix(ill):
    sd = colour.SDS_ILLUMINANTS[ill]
    e = resample(sd.values, sd.wavelengths)
    cm = np.stack([resample(CMFS.values[:, k], CMFS.wavelengths) for k in range(3)])
    m = cm * e
    return m / m[1].sum()


def lab(R, ill):
    M = xyz_matrix(ill)
    white = M.sum(1)
    return colour.XYZ_to_Lab(np.atleast_2d(R) @ M.T, colour.XYZ_to_xy(white))


def random_skins(n, seed=0):
    rng = np.random.default_rng(seed)
    U = rng.random((n, len(RANGES)))
    return U, np.array([reflectance.diffuse_reflectance(LAM, from_unit(u)) for u in U])


def camera_metamerism():
    cam = colour.MSDS_CAMERA_SENSITIVITIES["Nikon 5100 (NPL)"]
    S = np.stack([resample(cam.values[:, k], cam.wavelengths) for k in range(3)])
    d65 = resample(colour.SDS_ILLUMINANTS["D65"].values, colour.SDS_ILLUMINANTS["D65"].wavelengths)
    Mcam = S * d65
    Mcam = Mcam / Mcam.sum(1, keepdims=True)  # white balance
    Mxyz = xyz_matrix("D65")
    # ColorChecker calibration (least squares camera RGB -> XYZ).
    cc = colour.SDS_COLOURCHECKERS["ColorChecker N Ohta"]
    Rcc = np.array([resample(sd.values, sd.wavelengths) for sd in cc.values()])
    C = np.linalg.lstsq(Rcc @ Mcam.T, Rcc @ Mxyz.T, rcond=None)[0]
    _, Rsk = random_skins(4000, seed=3)
    xyz_true = Rsk @ Mxyz.T
    xyz_cam = Rsk @ Mcam.T @ C
    w = colour.XYZ_to_xy(Mxyz.sum(1))
    de = colour.delta_E(colour.XYZ_to_Lab(xyz_cam, w), colour.XYZ_to_Lab(xyz_true, w), method="CIE 2000")
    # Same calibration but fitted on skin spectra (a "skin-tuned" chart).
    C2 = np.linalg.lstsq(Rsk[::2] @ Mcam.T, Rsk[::2] @ Mxyz.T, rcond=None)[0]
    de2 = colour.delta_E(colour.XYZ_to_Lab(Rsk[1::2] @ Mcam.T @ C2, w), colour.XYZ_to_Lab(xyz_true[1::2], w), method="CIE 2000")
    print("\n| Calibration caméra (Nikon D5100) | ΔE00 médian | ΔE00 p95 | ΔE00 max |\n|---|---|---|---|")
    for name, d in [("ColorChecker 24", de), ("Spectres de peau", de2)]:
        print(f"| {name} | {np.median(d):.2f} | {np.percentile(d, 95):.2f} | {d.max():.2f} |")


def illuminant_metamerism(n=40000):
    U, R = random_skins(n)
    labs_d65 = lab(R, "D65")
    targets = {
        "claire": SkinParams(melanin_fraction=0.02, blood_fraction=0.015, bilirubin_umol=1.0),
        "moyenne": SkinParams(melanin_fraction=0.08, bilirubin_umol=1.0),
        "foncée": SkinParams(melanin_fraction=0.30, eumelanin_ratio=0.9, bilirubin_umol=1.0),
    }
    print("\n| Peau | # métamères (ΔE00<1 sous D65) | " + " | ".join(f"{ILL_LABELS[i]} méd / max" for i in ILLUMINANTS) + " |")
    print("|---|---|" + "---|" * len(ILLUMINANTS))
    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = np.arange(len(ILLUMINANTS))
    for j, ((name, p), c) in enumerate(zip(targets.items(), SERIES)):
        R0 = reflectance.diffuse_reflectance(LAM, p)
        de0 = colour.delta_E(labs_d65, lab(R0, "D65"), method="CIE 2000")
        idx = np.where(de0 < 1.0)[0]
        stats = []
        for ill in ILLUMINANTS:
            d = colour.delta_E(lab(R[idx], ill), lab(R0, ill), method="CIE 2000")
            stats.append((np.median(d), np.percentile(d, 95), d.max()))
        print(f"| {name} | {len(idx)} | " + " | ".join(f"{m:.2f} / {mx:.2f}" for m, _, mx in stats) + " |")
        s = np.array(stats)
        ax.bar(x + (j - 1) * 0.26, s[:, 1], width=0.24, color=c, label=f"Peau {name}")
    ax.axhline(1.0, color="#52514e", lw=1, ls="--")
    ax.text(-0.45, 1.02, "seuil de visibilité ≈ 1", ha="left", va="bottom", fontsize=8, color="#52514e")
    ax.set_xticks(x, [ILL_LABELS[i] for i in ILLUMINANTS], fontsize=8, rotation=20, ha="right")
    ax.set_ylim(0, 1.45)
    ax.set_ylabel("ΔE00 (95e centile) vs référence")
    ax.set_title("Métamères RGB (D65) rendus sous éclairages de plateau", loc="left")
    ax.legend(fontsize=8, loc="upper right", ncol=3)
    save(fig, "illuminant_metamerism.png")


if __name__ == "__main__":
    camera_metamerism()
    illuminant_metamerism()
