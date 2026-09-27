"""Digital Emily (WikiHuman / USC ICT Light Stage) : inversion and specular separation test.

Data (data/Emily/, research licence of the WikiHuman project):
- 00_diffuse_unlit_raw.exr : cross-polarised diffuse albedo, UV space, linear
- 00_specular_unlit_raw.exr : polarisation-difference specular intensity, UV space (arbitrary scale)
- cam1_mixed_w.exr : unpolarised photograph, other expression, image space (reference only)

1. Colour-space check against the measured skin gamut.
2. Chromophore maps from the true cross-polarised diffuse (reference maps).
3. Specular separation: mixed = diffuse + k * specular (white), for several k. The diffuse is
   recovered by removing the smallest amount of white that brings the colour back onto the
   measured skin gamut (dichromatic model). Compared with the true diffuse and with no separation.
"""

import os
import sys

import colour
import numpy as np
import OpenEXR
from PIL import Image
from scipy.spatial import cKDTree

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from common import SERIES, plt, save  # noqa: E402

from chromophores import inversion as inv, lut3d  # noqa: E402

ROOT = os.path.join(HERE, "..", "..")
DATA = os.path.join(ROOT, "data", "Emily")
LUT = os.path.join(ROOT, "chromophores", "data", "rgb_lut_face_65_v1.npz")
REPORT = os.path.join(ROOT, "docs", "phase2", "EMILY.md")
ID = {n: i for i, n in enumerate(inv.NAMES)}
N = 512  # working resolution for the separation study
K_LEVELS = [0.05, 0.1, 0.2]  # specular scale relative to the (normalised) specular map
TOL = 0.02  # gamut tolerance (encoded RGB distance)


def load(name):
    with OpenEXR.File(os.path.join(DATA, name + ".exr")) as e:
        ch = e.channels()
        return ch[list(ch)[0]].pixels.astype(np.float64)[..., :3]


def resize(a, n):
    return np.stack([np.asarray(Image.fromarray(a[..., c].astype(np.float32)).resize((n, n), Image.BILINEAR)) for c in range(3)], -1)


def de00(a, b):
    la = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(a, 0, None), apply_cctf_decoding=False))
    lb = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(b, 0, None), apply_cctf_decoding=False))
    return colour.delta_E(la, lb, method="CIE 2000")


def separate(mixed, tree, s_max=0.6, steps=121):
    """Remove the smallest white offset s >= 0 that brings the colour onto the skin gamut."""
    flat = mixed.reshape(-1, 3)
    grid = np.linspace(0.0, s_max, steps)
    dist = np.stack([tree.query(lut3d.encode(np.clip(flat - s, 0, None)))[0] for s in grid], 1)
    ok = dist <= TOL
    first = np.where(ok.any(1), ok.argmax(1), dist.argmin(1))
    s = grid[first]
    return np.clip(flat - s[:, None], 0, None).reshape(mixed.shape), s.reshape(mixed.shape[:-1])


def main():
    diff_raw = load("00_diffuse_unlit_raw")
    spec = load("00_specular_unlit_raw")[..., 0]
    ref_tree = cKDTree(lut3d.encode(lut3d.reference_skin_rgb()))
    lut = lut3d.RGBToChromophores(LUT)
    gains, in_before, in_after = lut3d.autocalibrate(diff_raw)
    diff = diff_raw * gains

    # 1. colour-space check
    rng = np.random.default_rng(0)
    pix = diff_raw.reshape(-1, 3)[rng.choice(diff.shape[0] * diff.shape[1], 200_000, replace=False)]
    gam_lin = np.mean(ref_tree.query(lut3d.encode(pix))[0] <= 0.06) * 100
    gam_enc = np.mean(ref_tree.query(lut3d.encode(lut3d.decode(pix)))[0] <= 0.06) * 100

    # 2. reference maps from the true diffuse
    d = resize(diff, N)
    sp = resize(np.repeat(spec[..., None], 3, -1), N)[..., 0]
    out_true = lut(d)
    skin = out_true["in_gamut"]

    # 3. separation study
    rows = []
    examples = {}
    for k in K_LEVELS:
        mixed = d + k * sp[..., None]
        rec, s_est = separate(mixed, ref_tree)
        out_mix, out_rec = lut(mixed), lut(rec)
        dE_mix = de00(d, mixed)[skin]
        dE_rec = de00(d, rec)[skin]
        err = {}
        for n in ("melanin_fraction", "blood_fraction"):
            t = out_true["map"][..., ID[n]][skin]
            err[n] = (np.median(np.abs(out_mix["map"][..., ID[n]][skin] / t - 1)), np.median(np.abs(out_rec["map"][..., ID[n]][skin] / t - 1)))
        corr = np.corrcoef(s_est[skin], (k * sp)[skin])[0, 1]
        rows.append((k, np.median((k * sp)[skin]), np.median(dE_mix), np.median(dE_rec), np.percentile(dE_rec, 95), err, corr,
                     100 * np.mean(out_mix["in_gamut"][skin])))
        examples[k] = (mixed, rec, s_est, out_mix, out_rec)

    # figures
    shown = [("melanin_fraction", "Mélanine"), ("blood_fraction", "Sang"), ("oxygen_saturation", "SO₂")]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.4))
    axes[0].imshow(lut3d.encode(d))
    axes[0].set_title("Diffuse cross-pol (Emily)", loc="left", fontsize=9)
    for ax, (n, lab) in zip(axes[1:], shown):
        v = np.where(skin, out_true["map"][..., ID[n]], np.nan)
        lo, hi = np.nanpercentile(v, [2, 98])
        im = ax.imshow(v, cmap="magma" if n != "oxygen_saturation" else "viridis", vmin=lo, vmax=hi)
        ax.set_title(lab, loc="left", fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.04)
    for ax in axes:
        ax.axis("off")
    save(fig, "phase2_emily_maps.png")

    k = 0.1
    mixed, rec, s_est, out_mix, out_rec = examples[k]
    fig, axes = plt.subplots(2, 4, figsize=(16, 8.4))
    for ax, img, title in zip(axes[0], [d, mixed, rec], ["Diffuse vraie (cross-pol)", f"Mixte synthétique (k = {k})", "Diffuse récupérée (séparation)"]):
        ax.imshow(lut3d.encode(img))
        ax.set_title(title, loc="left", fontsize=9)
    im = axes[0, 3].imshow(np.where(skin, s_est, np.nan), cmap="Greys")
    axes[0, 3].set_title("Spéculaire estimé (blanc retiré)", loc="left", fontsize=9)
    fig.colorbar(im, ax=axes[0, 3], fraction=0.04)
    n = "melanin_fraction"
    t = np.where(skin, out_true["map"][..., ID[n]], np.nan)
    lo, hi = np.nanpercentile(t, [2, 98])
    for ax, v, title in zip(axes[1], [t, np.where(skin, out_mix["map"][..., ID[n]], np.nan), np.where(skin, out_rec["map"][..., ID[n]], np.nan)],
                            ["Mélanine (vraie diffuse)", "Mélanine sans séparation", "Mélanine après séparation"]):
        ax.imshow(v, cmap="magma", vmin=lo, vmax=hi)
        ax.set_title(title, loc="left", fontsize=9)
    im = axes[1, 3].imshow(np.where(skin, de00(d, rec), np.nan), cmap="Greys", vmin=0, vmax=3)
    axes[1, 3].set_title("ΔE00 diffuse récupérée vs vraie", loc="left", fontsize=9)
    fig.colorbar(im, ax=axes[1, 3], fraction=0.04)
    for ax in axes.ravel():
        ax.axis("off")
    save(fig, "phase2_emily_separation.png")

    lines = ["| k (échelle spéculaire) | spéculaire médian ajouté | ΔE00 mixte vs vraie | ΔE00 récupérée vs vraie (méd. / p95) | "
             "erreur mélanine sans / avec séparation | erreur sang sans / avec séparation | corrélation spéculaire estimé / vrai | % texels mixtes dans le gamut |",
             "|---|---|---|---|---|---|---|---|"]
    for k, sm, dm, dr, dr95, err, corr, gm in rows:
        lines.append(f"| {k} | {sm:.3f} | {dm:.2f} | {dr:.2f} / {dr95:.2f} | {100 * err['melanin_fraction'][0]:.0f} % / {100 * err['melanin_fraction'][1]:.0f} % | "
                     f"{100 * err['blood_fraction'][0]:.0f} % / {100 * err['blood_fraction'][1]:.0f} % | {corr:.2f} | {gm:.0f} % |")
    med = {n: np.median(out_true["map"][..., ID[n]][skin]) for n in inv.NAMES}
    text = f"""# Digital Emily : inversion et séparation du spéculaire

Données : `data/Emily/` (projet WikiHuman, USC ICT ; licence recherche). Textures UV 1280², étude à {N}².

## 1. Espace colorimétrique

Texels (tous, cheveux et vêtements compris) dans le gamut peau : **{gam_lin:.0f} %** lus comme linéaires,
{gam_enc:.0f} % lus comme sRGB encodé. **La diffuse d'Emily est bien linéaire** (contrairement au fichier VFace).

Auto-calibration (exposition + balance R/B, `lut3d.autocalibrate`) : gains (R, G, B) = ({gains[0]:.3f}, {gains[1]:.3f}, {gains[2]:.3f}) ;
texels à moins de 0,03 du gamut peau : {100 * in_before:.0f} % → {100 * in_after:.0f} %. Le gain d'exposition (≈ ×{gains[1]:.2f}) est appliqué
dans la suite. Attention : une exposition globale et un teint globalement plus clair sont en partie confondus ;
ce gain doit être validé par un artiste ou une référence.

## 2. Cartes depuis la vraie diffuse (polarisation croisée)

Texels dans le gamut peau : {100 * skin.mean():.0f} % (cheveux, yeux, bouche et vêtements exclus automatiquement en grande partie).
Médianes : mélanine {med['melanin_fraction']:.3f}, sang {med['blood_fraction']:.4f}, SO₂ {med['oxygen_saturation']:.2f}.

![cartes](../figures/phase2_emily_maps.png)

## 3. Séparation du spéculaire (mixte synthétique = diffuse + k × spéculaire mesuré)

Méthode : retirer le plus petit blanc s ≥ 0 qui ramène la couleur sur le gamut de la peau mesurée
(tolérance {TOL} en RGB encodé). Erreurs relatives médianes sur les cartes, par rapport aux cartes de la vraie diffuse.

{chr(10).join(lines)}

![séparation](../figures/phase2_emily_separation.png)

Note : l'échelle absolue de la texture spéculaire est inconnue ; k couvre des niveaux de reflet faibles à forts.
"""
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    open(REPORT, "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
