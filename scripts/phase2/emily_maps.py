"""Full-resolution chromophore maps from the Digital Emily cross-polarised diffuse."""

import os
import sys

import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from common import plt, save  # noqa: E402
from emily import LUT, load  # noqa: E402

from chromophores import inversion as inv, lut3d  # noqa: E402

ID = {n: i for i, n in enumerate(inv.NAMES)}
CROP = (330, 470, 360, 300)  # x, y, w, h: under-eye / cheek area (freckles, pigmentation)
MAPS = [("melanin_fraction", "Mélanine (fraction épiderme)", "magma"),
        ("blood_fraction", "Hémoglobine (fraction sanguine derme)", "magma"),
        ("oxygen_saturation", "Oxygénation SO₂", "viridis"),
        ("epidermis_thickness_cm", "Épaisseur épiderme (cm)", "viridis"),
        ("eumelanin_ratio", "Ratio eu/phéo (β, peu fiable)", "cividis"),
        ("scattering_scale", "Échelle de diffusion μs′ (a priori)", "cividis")]


def show(ax, v, cmap, title, mask, pct=(2, 98)):
    v = np.where(mask, v, np.nan)
    lo, hi = np.nanpercentile(v, pct)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("#cfcfcf")
    im = ax.imshow(v, cmap=cm, vmin=lo, vmax=hi)
    ax.set_title(title, loc="left", fontsize=9)
    ax.axis("off")
    plt.colorbar(im, ax=ax, fraction=0.04)


def main():
    raw = load("00_diffuse_unlit_raw")
    gains, _, _ = lut3d.autocalibrate(raw)
    d = raw * gains
    lut = lut3d.RGBToChromophores(LUT)
    out = lut(d)
    skin = out["in_gamut"]

    fig, axes = plt.subplots(3, 4, figsize=(18, 14))
    axes[0, 0].imshow(lut3d.encode(d))
    x, y, w, h = CROP
    axes[0, 0].add_patch(plt.Rectangle((x, y), w, h, fill=False, ec="#0b0b0b", lw=1))
    axes[0, 0].set_title(f"Diffuse cross-pol (calibrée ×{gains[1]:.2f})", loc="left", fontsize=9)
    axes[0, 0].axis("off")
    for ax, (n, lab, cmap) in zip([axes[0, 1], axes[0, 2], axes[0, 3], axes[1, 0], axes[1, 1], axes[1, 2]], MAPS):
        show(ax, out["map"][..., ID[n]], cmap, lab, skin)
    rel = lambda n: np.log(out["p84"][..., ID[n]] / out["p16"][..., ID[n]]) / 2  # noqa: E731
    show(axes[1, 3], rel("melanin_fraction"), "Greys", "Incertitude mélanine (σ relatif)", skin)
    show(axes[2, 0], rel("blood_fraction"), "Greys", "Incertitude hémoglobine (σ relatif)", skin)
    show(axes[2, 1], rel("oxygen_saturation"), "Greys", "Incertitude SO₂ (σ relatif)", skin)
    axes[2, 2].imshow(np.where(skin, 1.0, 0.0), cmap="Greys_r")
    axes[2, 2].set_title(f"Gamut peau ({100 * skin.mean():.0f} % des texels ; gris = hors peau)", loc="left", fontsize=9)
    axes[2, 2].axis("off")
    axes[2, 3].axis("off")
    fig.suptitle("Digital Emily — cartes de chromophores depuis la diffuse (LUT 65³, a priori visage)", x=0.01, ha="left")
    save(fig, "phase2_emily_fullmaps.png")

    sl = (slice(y, y + h), slice(x, x + w))
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))
    axes[0].imshow(lut3d.encode(d[sl]))
    axes[0].set_title("Diffuse (zoom, pleine résolution)", loc="left", fontsize=9)
    axes[0].axis("off")
    for ax, (n, lab, cmap) in zip(axes[1:], MAPS[:3]):
        show(ax, out["map"][sl + (ID[n],)], cmap, lab, skin[sl], pct=(1, 99))
    save(fig, "phase2_emily_zoom.png")
    for n, *_ in MAPS:
        v = out["map"][..., ID[n]][skin]
        print(f"{n:24s} p5 {np.percentile(v, 5):.4f}  médiane {np.median(v):.4f}  p95 {np.percentile(v, 95):.4f}")


if __name__ == "__main__":
    main()
