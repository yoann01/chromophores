"""Phase 2: chromophore and uncertainty maps from albedo textures.

Test data: the two example face albedos of the BioSkin repository (MIT,
facebookresearch/BioSkin, test_images/), sRGB-encoded JPEG, with skin masks.
They are assumed calibrated, cross-polarised, linear-sRGB after decoding.
"""

import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from common import plt, save  # noqa: E402

from chromophores import inversion as inv, lut3d  # noqa: E402

BIOSKIN = os.environ.get("BIOSKIN_DIR", "/home/user/facebookresearch/bioskin")
DATA = os.path.join(HERE, "..", "..", "chromophores", "data")
HEIGHT = 900


def load(name, mask_name):
    im = Image.open(os.path.join(BIOSKIN, "test_images", name)).convert("RGB")
    w, h = im.size
    size = (int(w * HEIGHT / h), HEIGHT)
    rgb = np.asarray(im.resize(size, Image.LANCZOS), float) / 255.0
    mask = np.asarray(Image.open(os.path.join(BIOSKIN, "test_images", mask_name)).convert("L").resize(size, Image.BILINEAR), float) / 255.0 > 0.5
    return lut3d.decode(rgb), mask


def main():
    path = os.path.join(DATA, "rgb_lut_face_65_v1.npz")
    if not os.path.exists(path):
        path = os.path.join(DATA, "rgb_lut_face_v1.npz")
    lut = lut3d.RGBToChromophores(path)
    shown = [("melanin_fraction", "Mélanine", "magma"), ("blood_fraction", "Sang", "magma"), ("oxygen_saturation", "SO₂", "viridis"),
             ("epidermis_thickness_cm", "Épaisseur épiderme (cm)", "viridis")]
    ID = {n: i for i, n in enumerate(inv.NAMES)}
    stats = []
    for name, mname in (("1.jpg", "1_mask.jpg"), ("2.jpg", "2_mask.jpg")):
        lin, mask = load(name, mname)
        out = lut(lin)
        fig, axes = plt.subplots(2, 5, figsize=(16, 7.4))
        axes[0, 0].imshow(lut3d.encode(lin))
        axes[0, 0].set_title("Albédo d'entrée", loc="left", fontsize=9)
        gam = out["in_gamut"] & mask
        axes[1, 0].imshow(np.where(mask, np.where(out["in_gamut"], 1.0, 0.0), np.nan), cmap="gray", vmin=0, vmax=1)
        axes[1, 0].set_title(f"Dans le gamut peau : {100 * gam.sum() / max(mask.sum(), 1):.1f} % des texels", loc="left", fontsize=9)
        for j, (n, lab, cmap) in enumerate(shown):
            k = ID[n]
            v = np.where(gam, out["map"][..., k], np.nan)
            lo, hi = np.nanpercentile(v, [2, 98])
            im = axes[0, j + 1].imshow(v, cmap=cmap, vmin=lo, vmax=hi)
            axes[0, j + 1].set_title(lab, loc="left", fontsize=9)
            fig.colorbar(im, ax=axes[0, j + 1], fraction=0.035)
            spread = np.where(gam, np.log(out["p84"][..., k] / out["p16"][..., k]) / 2, np.nan)  # ~ relative 1-sigma
            im = axes[1, j + 1].imshow(spread, cmap="Greys", vmin=0, vmax=np.nanpercentile(spread, 98))
            axes[1, j + 1].set_title(f"{lab} — incertitude (σ relatif)", loc="left", fontsize=9)
            fig.colorbar(im, ax=axes[1, j + 1], fraction=0.035)
            stats.append((name, lab, np.nanmedian(v), np.nanmedian(spread)))
        for ax in axes.ravel():
            ax.axis("off")
        save(fig, f"phase2_maps_{name.split('.')[0]}.png")
    for s in stats:
        print(f"{s[0]} {s[1]:28s} médiane {s[2]:.4f}  σ relatif médian {s[3]:.2f}")


if __name__ == "__main__":
    main()
