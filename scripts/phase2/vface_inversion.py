"""Phase 2 on a real production albedo (TexturingXYZ VFace, UDIM 1001).

1. Colour-space check against the measured skin gamut (auto-calibration test).
2. Chromophore + uncertainty maps (face prior LUT, 65^3).
3. Full-resolution crop to check that the albedo detail is kept.
4. Closure: albedo re-rendered from the MAP chromophores (model only) vs input,
   i.e. the size of the residual needed for an exact reconstruction.
"""

import os
import sys

import colour
import numpy as np
from PIL import Image
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from common import plt, save  # noqa: E402

from chromophores import inversion as inv, lut3d  # noqa: E402

Image.MAX_IMAGE_PIXELS = None
ROOT = os.path.join(HERE, "..", "..")
SRC = os.path.join(ROOT, "data", "VFace", "XYZ_albedo_lin_srgb.1001.png")
LUT = os.path.join(ROOT, "chromophores", "data", "rgb_lut_face_65_v1.npz")
REPORT = os.path.join(ROOT, "docs", "phase2", "VFACE_1001.md")
ID = {n: i for i, n in enumerate(inv.NAMES)}
CROP = (1150, 1700, 640, 640)  # x, y, w, h at full resolution (left cheek)


def gamut_check(enc, ref_tree):
    rng = np.random.default_rng(0)
    pix = enc.reshape(-1, 3)[rng.choice(enc.shape[0] * enc.shape[1], 200_000, replace=False)]
    rows = []
    for label, lin in (("lu comme linéaire (nom du fichier)", pix), ("lu comme sRGB encodé", lut3d.decode(pix))):
        d, _ = ref_tree.query(lut3d.encode(lin))
        rows.append((label, 100 * np.mean(d <= 0.06), np.median(d), np.median(lin, 0)))
    return rows


def node_reconstruction(lut_path):
    """Albedo re-rendered from each LUT node's MAP (model only, no residual)."""
    t = np.load(lut_path)
    fwd = inv.Forward(specular=0.0)
    x = t["x"].reshape(-1, inv.N)
    rgb = np.full((len(x), 3), np.nan)
    for i in np.where(np.isfinite(x[:, 0]))[0]:
        rgb[i] = lut3d.M_SRGB @ fwd.spectrum(x[i])
    n = int(t["n"])
    grid = rgb.reshape(n, n, n, 3)
    valid = np.isfinite(grid[..., 0])
    idx_v = np.argwhere(valid)
    _, nn = cKDTree(idx_v).query(np.argwhere(~valid))
    grid[tuple(np.argwhere(~valid).T)] = grid[tuple(idx_v[nn].T)]
    return RegularGridInterpolator((t["axis"],) * 3, grid)


def de00(a, b):
    la = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(a, 0, None), apply_cctf_decoding=False))
    lb = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(b, 0, None), apply_cctf_decoding=False))
    return colour.delta_E(la, lb, method="CIE 2000")


def main():
    img = Image.open(SRC).convert("RGB")
    full = np.asarray(img, float) / 255.0  # sRGB-encoded values (see gamut check)
    ref_tree = cKDTree(lut3d.encode(lut3d.reference_skin_rgb()))
    gam = gamut_check(full, ref_tree)

    lut = lut3d.RGBToChromophores(LUT)
    recon = node_reconstruction(LUT)

    small = np.asarray(img.resize((1024, 1024), Image.LANCZOS), float) / 255.0
    lin = lut3d.decode(small)
    out = lut(lin)
    rec = recon(lut3d.encode(lin).reshape(-1, 3)).reshape(lin.shape)
    dE = de00(lin, rec)
    ing = out["in_gamut"]

    x, y, w, h = CROP
    crop_enc = full[y : y + h, x : x + w]
    crop = lut(lut3d.decode(crop_enc))

    shown = [("melanin_fraction", "Mélanine", "magma"), ("blood_fraction", "Sang", "magma"),
             ("oxygen_saturation", "SO₂", "viridis"), ("epidermis_thickness_cm", "Épaisseur épiderme (cm)", "viridis")]
    fig, axes = plt.subplots(2, 5, figsize=(17, 7.6))
    axes[0, 0].imshow(small)
    axes[0, 0].add_patch(plt.Rectangle((x / 4, y / 4), w / 4, h / 4, fill=False, ec="#0b0b0b", lw=1))
    axes[0, 0].set_title("Albédo VFace (sRGB)", loc="left", fontsize=9)
    im = axes[1, 0].imshow(np.where(ing, dE, np.nan), cmap="Greys", vmin=0, vmax=2)
    axes[1, 0].set_title("ΔE00 modèle seul vs albédo (résidu)", loc="left", fontsize=9)
    fig.colorbar(im, ax=axes[1, 0], fraction=0.035)
    for j, (n, lab, cmap) in enumerate(shown):
        k = ID[n]
        v = np.where(ing, out["map"][..., k], np.nan)
        lo, hi = np.nanpercentile(v, [2, 98])
        im = axes[0, j + 1].imshow(v, cmap=cmap, vmin=lo, vmax=hi)
        axes[0, j + 1].set_title(lab, loc="left", fontsize=9)
        fig.colorbar(im, ax=axes[0, j + 1], fraction=0.035)
        spread = np.where(ing, np.log(out["p84"][..., k] / out["p16"][..., k]) / 2, np.nan)
        im = axes[1, j + 1].imshow(spread, cmap="Greys", vmin=0, vmax=np.nanpercentile(spread, 98))
        axes[1, j + 1].set_title(f"{lab} — incertitude (σ rel.)", loc="left", fontsize=9)
        fig.colorbar(im, ax=axes[1, j + 1], fraction=0.035)
    for ax in axes.ravel():
        ax.axis("off")
    save(fig, "phase2_vface_maps.png")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    axes[0].imshow(crop_enc)
    axes[0].set_title("Albédo, pleine résolution (joue)", loc="left", fontsize=10)
    for ax, (n, lab) in zip(axes[1:], [("melanin_fraction", "Mélanine"), ("blood_fraction", "Sang")]):
        v = crop["map"][..., ID[n]]
        lo, hi = np.percentile(v, [1, 99])
        im = ax.imshow(v, cmap="magma", vmin=lo, vmax=hi)
        ax.set_title(f"{lab}, pleine résolution", loc="left", fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.04)
    for ax in axes:
        ax.axis("off")
    save(fig, "phase2_vface_crop.png")

    P = out["map"][ing]
    lines_g = "\n".join(f"| {lab} | {pct:.1f} % | {d:.3f} | ({m[0]:.3f}, {m[1]:.3f}, {m[2]:.3f}) |" for lab, pct, d, m in gam)
    med = {n: np.median(P[:, ID[n]]) for n in inv.NAMES}
    text = f"""# VFace UDIM 1001 : inversion sur un albédo de production

Fichier : `data/VFace/XYZ_albedo_lin_srgb.1001.png` (4096², RGBA 8 bits, alpha constant).

## 1. Espace colorimétrique (auto-calibration)

| Interprétation | Texels dans le gamut peau mesuré (ISSA + NIST) | Distance médiane (encodé) | RGB linéaire médian |
|---|---|---|---|
{lines_g}

Référence ISSA + NIST : RGB linéaire médian (0,374 ; 0,209 ; 0,132), p95 (0,487 ; 0,313 ; 0,223).
**Malgré son nom, le fichier est encodé en sRGB (gamma)** : lu comme linéaire, il est hors du gamut
de la peau réelle pour 93 % des texels.

## 2. Cartes (LUT 65³, a priori « visage »), sur 1024²

- Texels dans le gamut peau : {100 * ing.mean():.1f} %.
- Médianes : mélanine {med['melanin_fraction']:.4f}, sang {med['blood_fraction']:.4f}, SO₂ {med['oxygen_saturation']:.2f},
  épaisseur épiderme {med['epidermis_thickness_cm']:.4f} cm, échelle μs′ {med['scattering_scale']:.2f} (a priori).

![cartes](../figures/phase2_vface_maps.png)

## 3. Détail à pleine résolution

![joue](../figures/phase2_vface_crop.png)

## 4. Fermeture : albédo re-rendu par le modèle seul (sans résidu)

ΔE00 (D65) entre l'albédo d'entrée et l'albédo recalculé à partir des chromophores :
médiane {np.median(dE[ing]):.2f}, p95 {np.percentile(dE[ing], 95):.2f}, max {dE[ing].max():.2f}.
C'est l'amplitude du résidu (correction colorimétrique) nécessaire pour une reconstruction exacte.
"""
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    open(REPORT, "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
