"""Which chromophore parameters can be recovered from which measurement?

Local (Fisher / Cramer-Rao) identifiability analysis of the biophysical
parameters under four acquisition setups:

  1. RGB diffuse albedo (the input of RGB->spectral inversion methods)
  2. RGB albedo + per-channel diffusion length (what a renderer's
     "scattering radius" encodes; measurable with spatially resolved
     reflectance / structured light)
  3. 8-band narrow-band multispectral (LED multiplexing, cross-polarised)
  4. Full spectrum 380-780 nm, 10 nm (spectrophotometer / hyperspectral)

Parameters are normalised to their plausible range; the prior is a Gaussian
with std = 0.5 range. A posterior std close to 0.5 means "the data tells
us nothing, only the prior does". Also runs a metamer search: skins with
the same RGB albedo (Delta E00 < 1) but very different chromophores.
"""

import colour
import numpy as np
from common import SERIES, plt, save

from chromophores import color, optics, reflectance
from chromophores.optics import SkinParams
from chromophores.spectra import WAVELENGTHS as LAM

# name, low, high, log-scale
RANGES = [
    ("melanin_fraction", 0.01, 0.43, True),
    ("eumelanin_ratio", 0.0, 1.0, False),
    ("blood_fraction", 0.002, 0.07, True),
    ("oxygen_saturation", 0.5, 1.0, False),
    ("carotene_umol", 0.0, 3.0, False),
    ("bilirubin_umol", 0.0, 10.0, False),
    ("epidermis_thickness_cm", 0.005, 0.02, True),
    ("scattering_scale", 0.6, 1.5, True),
]
LABELS = ["Mélanine Cm", "Ratio eu/phéo", "Sang Ch", "SO2", "Carotène", "Bilirubine", "Épaisseur épiderme", "Diffusion µs'"]
PRIOR_STD = 0.5
NOISE = 0.005  # absolute reflectance noise per channel / band
NOISE_LOG_LENGTH = 0.05  # 5 % on diffusion length

M_RGB = color.reflectance_to_linear_srgb_matrix(LAM)
M_XYZ = color.reflectance_to_xyz_matrix(LAM)
M_MS8 = color.gaussian_bands_matrix([450, 480, 520, 545, 575, 600, 640, 700], fwhm=20.0, lam=LAM)
W_CH = np.abs(M_RGB) / np.abs(M_RGB).sum(1, keepdims=True)


def to_unit(p):
    u = []
    for name, lo, hi, lg in RANGES:
        v = getattr(p, name)
        u.append((np.log(v) - np.log(lo)) / (np.log(hi) - np.log(lo)) if lg else (v - lo) / (hi - lo))
    return np.array(u)


def from_unit(u):
    kw = {}
    for (name, lo, hi, lg), x in zip(RANGES, u):
        kw[name] = float(np.exp(np.log(lo) + x * (np.log(hi) - np.log(lo))) if lg else lo + x * (hi - lo))
    return SkinParams(**kw)


def measurements(u):
    p = from_unit(u)
    R = reflectance.diffuse_reflectance(LAM, p)
    _, der = optics.layers(LAM, p)
    sigma_tr = np.sqrt(3 * der.mua * (der.mua + der.musp))
    log_len = np.log(W_CH @ (1.0 / sigma_tr))
    rgb = M_RGB @ R
    return {
        "RGB": (rgb, NOISE),
        "RGB + rayon de diffusion": (np.concatenate([rgb / NOISE, log_len / NOISE_LOG_LENGTH]), 1.0),
        "Multispectral 8 bandes": (M_MS8 @ R, NOISE),
        "Spectre 380–780 nm": (R, NOISE),
    }


def posterior_std(u, eps=1e-4):
    base = measurements(u)
    jac = {k: [] for k in base}
    for i in range(len(u)):
        du = np.zeros_like(u)
        du[i] = eps
        hi, lo = measurements(u + du), measurements(u - du)
        for k in base:
            jac[k].append((hi[k][0] - lo[k][0]) / (2 * eps))
    out = {}
    for k, (_, sigma) in base.items():
        J = np.array(jac[k]).T / sigma
        F = J.T @ J + np.eye(len(u)) / PRIOR_STD**2
        out[k] = np.sqrt(np.diag(np.linalg.inv(F)))
    return out


def identifiability_figure():
    skins = {
        "claire": SkinParams(melanin_fraction=0.02, blood_fraction=0.015, bilirubin_umol=1.0),
        "moyenne": SkinParams(melanin_fraction=0.08, bilirubin_umol=1.0),
        "foncée": SkinParams(melanin_fraction=0.30, eumelanin_ratio=0.9, bilirubin_umol=1.0),
    }
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    table = []
    for ax, (skin, p) in zip(axes, skins.items()):
        res = posterior_std(to_unit(p))
        x = np.arange(len(LABELS))
        w = 0.2
        for j, ((k, s), c) in enumerate(zip(res.items(), SERIES)):
            ax.bar(x + (j - 1.5) * w, s, width=w - 0.02, color=c, label=k)
            table.append((skin, k, s))
        ax.axhline(PRIOR_STD, color="#52514e", lw=1, ls="--")
        ax.text(len(LABELS) - 0.5, PRIOR_STD + 0.01, "a priori seul", ha="right", va="bottom", fontsize=8, color="#52514e")
        ax.set_xticks(x, LABELS, rotation=45, ha="right", fontsize=8)
        ax.set_title(f"Peau {skin}", loc="left")
    axes[0].set_ylabel("Écart-type a posteriori (fraction de la plage)")
    axes[0].legend(fontsize=8, loc="upper left")
    save(fig, "identifiability.png")

    print("\n| Peau | Mesure | " + " | ".join(LABELS) + " |")
    print("|---|---|" + "---|" * len(LABELS))
    for skin, k, s in table:
        print(f"| {skin} | {k} | " + " | ".join(f"{v:.2f}" for v in s) + " |")


def metamer_figure(n=30000, seed=0):
    rng = np.random.default_rng(seed)
    target = to_unit(SkinParams(melanin_fraction=0.08))
    U = rng.random((n, len(RANGES)))
    R = np.array([reflectance.diffuse_reflectance(LAM, from_unit(u)) for u in U])
    lab = colour.XYZ_to_Lab((M_XYZ @ R.T).T)
    R0 = reflectance.diffuse_reflectance(LAM, from_unit(target))
    lab0 = colour.XYZ_to_Lab(M_XYZ @ R0)
    de = colour.delta_E(lab, lab0, method="CIE 2000")
    idx = np.where(de < 1.0)[0]
    print(f"\n{len(idx)} / {n} random skins are RGB-metamers (ΔE00 < 1) of the medium skin")
    if len(idx) == 0:
        return
    spread = U[idx].max(0) - U[idx].min(0)
    print("| Paramètre | Plage couverte par les métamères (fraction de la plage totale) |\n|---|---|")
    for lab_, s in zip(LABELS, spread):
        print(f"| {lab_} | {s:.2f} |")
    # Show the most spectrally different metamers.
    diff = np.abs(R[idx] - R0).max(1)
    pick = idx[np.argsort(diff)[-3:]]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(LAM, R0, color=SERIES[0], label="Référence (peau moyenne)")
    for j, k in enumerate(pick):
        p = from_unit(U[k])
        axes[0].plot(LAM, R[k], color=SERIES[j + 1], label=f"Métamère {j + 1} (ΔE00={de[k]:.2f})")
        _, der0 = optics.layers(LAM, from_unit(target))
        _, der = optics.layers(LAM, p)
        l0 = 1 / np.sqrt(3 * der0.mua * (der0.mua + der0.musp))
        axes[1].plot(LAM, 10 * np.sqrt(3 * der.mua * (der.mua + der.musp)) ** -1, color=SERIES[j + 1], label=f"Métamère {j + 1}")
        print(f"metamer {j + 1}: " + ", ".join(f"{a}={getattr(p, a):.3g}" for a, *_ in RANGES))
    axes[1].plot(LAM, 10 * l0, color=SERIES[0], label="Référence")
    axes[0].set_title("Même albédo RGB, spectres différents", loc="left")
    axes[0].set_xlabel("Longueur d'onde (nm)")
    axes[0].set_ylabel("Réflectance diffuse")
    axes[0].legend(fontsize=8)
    axes[1].set_title("…et longueurs de diffusion différentes", loc="left")
    axes[1].set_xlabel("Longueur d'onde (nm)")
    axes[1].set_ylabel("1/σtr dans le derme (mm)")
    axes[1].legend(fontsize=8)
    save(fig, "rgb_metamers.png")


if __name__ == "__main__":
    identifiability_figure()
    metamer_figure()
