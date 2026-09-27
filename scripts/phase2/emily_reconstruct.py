"""Rebuild the Emily diffuse albedo from the extracted chromophore maps only.

Each texel's physical parameters (melanin, blood, SO2, epidermal thickness, beta, scattering)
go through the forward model (layered table -> diffuse reflectance spectrum 380-780 nm ->
linear sRGB). No residual/colour correction is applied: this measures what the maps alone
reproduce. The spectra also give the albedo under another illuminant (tungsten A).
"""

import os
import sys

import colour
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from common import plt, save  # noqa: E402
from emily import LUT, load  # noqa: E402
from emily_maps import CROP  # noqa: E402

from chromophores import inversion as inv, lut3d, table  # noqa: E402
from chromophores.forward import SPECULAR  # noqa: E402
from chromophores.optics import SkinParams  # noqa: E402

LAM = inv.LAMBDA
T = table.LayeredTable(method="linear")


def spectra_from_maps(P, chunk=40_000):
    """Vectorised forward model: physical parameters (N, 6) -> diffuse albedo spectra (N, n_lambda)."""
    out = np.empty((len(P), len(LAM)))
    for i in range(0, len(P), chunk):
        c = P[i : i + chunk]
        p = SkinParams(bilirubin_umol=0.0, **{n: c[:, k : k + 1] for k, n in enumerate(inv.NAMES)})
        tau, ad, de, g, _ = table.dimensionless(LAM, p)
        out[i : i + chunk] = (1 - SPECULAR) * np.exp(T._A(T.coords(tau, ad, de, np.broadcast_to(g, tau.shape))))
    return out


def rgb_under(S, ill):
    shape = colour.SpectralShape(380, 780, 10)
    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"].copy().align(shape)
    E = colour.SDS_ILLUMINANTS[ill].copy().align(shape).values
    M = cmfs.values.T * E
    M = M / M[1].sum()
    xyz = S @ M.T
    # von Kries (CAT02) to D65 for display
    xyz_d65 = colour.chromatic_adaptation(xyz, colour.xy_to_XYZ(colour.XYZ_to_xy(M.sum(1))), colour.xy_to_XYZ(colour.CCS_ILLUMINANTS["CIE 1931 2 Degree Standard Observer"]["D65"]))
    return colour.XYZ_to_sRGB(xyz_d65, apply_cctf_encoding=False)


def de00(a, b):
    la = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(a, 0, None), apply_cctf_decoding=False))
    lb = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(b, 0, None), apply_cctf_decoding=False))
    return colour.delta_E(la, lb, method="CIE 2000")


def main():
    raw = load("00_diffuse_unlit_raw")
    gains, _, _ = lut3d.autocalibrate(raw)
    d = raw * gains
    out = lut3d.RGBToChromophores(LUT)(d)
    skin = out["in_gamut"]
    H, W = skin.shape
    P = out["map"].reshape(-1, inv.N)
    S = spectra_from_maps(P)
    rec = (S @ lut3d.M_SRGB.T).reshape(H, W, 3)
    rec_A = rgb_under(S, "A").reshape(H, W, 3)
    dE = de00(d, rec)
    shown = np.where(skin[..., None], rec, d * 0 + 0.8)

    fig, axes = plt.subplots(2, 3, figsize=(17, 11.5))
    axes[0, 0].imshow(lut3d.encode(d))
    axes[0, 0].set_title("Diffuse originale (calibrée)", loc="left", fontsize=10)
    axes[0, 1].imshow(lut3d.encode(shown))
    axes[0, 1].set_title("Reconstruite depuis les cartes seules (gris = hors peau)", loc="left", fontsize=10)
    im = axes[0, 2].imshow(np.where(skin, dE, np.nan), cmap="Greys", vmin=0, vmax=3)
    axes[0, 2].set_title("ΔE00 originale vs reconstruite", loc="left", fontsize=10)
    plt.colorbar(im, ax=axes[0, 2], fraction=0.04)
    x, y, w, h = CROP
    sl = (slice(y, y + h), slice(x, x + w))
    axes[1, 0].imshow(lut3d.encode(d[sl]))
    axes[1, 0].set_title("Zoom — originale", loc="left", fontsize=10)
    axes[1, 1].imshow(lut3d.encode(shown[sl]))
    axes[1, 1].set_title("Zoom — reconstruite", loc="left", fontsize=10)
    axes[1, 2].imshow(lut3d.encode(np.where(skin[..., None], rec_A, 0.8)))
    axes[1, 2].set_title("Bonus : mêmes cartes sous tungstène (A), adapté D65", loc="left", fontsize=10)
    for ax in axes.ravel():
        ax.axis("off")
    save(fig, "phase2_emily_reconstruction.png")
    v = dE[skin]
    print(f"texels peau : {100 * skin.mean():.0f} %  ΔE00 médiane {np.median(v):.2f}  p90 {np.percentile(v, 90):.2f}  p95 {np.percentile(v, 95):.2f}  p99 {np.percentile(v, 99):.2f}")
    print(f"part des texels peau avec ΔE00 < 1 : {100 * np.mean(v < 1):.1f} %, < 2 : {100 * np.mean(v < 2):.1f} %")


if __name__ == "__main__":
    main()
