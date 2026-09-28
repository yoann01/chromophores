"""Head-to-head on reflectance: our pipeline vs BioSkin (Aliaga et al. 2023; basis of the
reflectance part of Aliaga & Jarabo 2026), on measured data.

1. ISSA held-out spectra (the 300 of V3): RGB -> each pipeline -> reconstructed spectrum, compared
   with the MEASURED spectrum (400-700 nm) and in colour under D65 / A / FL11 / LED-B3.
   BioSkin is run in two conditions:
     - native: input RGB computed with BioSkin's own conversion (their best case);
     - production: standard linear sRGB D65 albedo (what a delivered texture is), fed in BGR as
       their loader does.
2. Digital Emily diffuse: maps and albedo reconstruction error for each pipeline.
"""

import os
import sys

import colour
import numpy as np

HERE = os.path.dirname(__file__)
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))
sys.path.insert(0, ROOT)
from bioskin_numpy import BioSkinDecoder, theirs_rgb, warp  # noqa: E402
from common import SERIES, plt, save  # noqa: E402

from chromophores import inversion as inv, lut3d  # noqa: E402
from chromophores.datasets import load_issa  # noqa: E402
from chromophores.forward import SPECULAR  # noqa: E402

BIOSKIN = os.environ.get("BIOSKIN_DIR", "/home/user/facebookresearch/bioskin")
D = BioSkinDecoder(os.path.join(BIOSKIN, "pretrained_models", "BioSkin"))
LAM_B = D.lam.astype(float)
LAM = inv.LAMBDA
SHAPE = colour.SpectralShape(380, 780, 10)
ILLUMS = ["D65", "A", "FL11", "LED-B3"]
REPORT = os.path.join(ROOT, "docs", "phase2", "BENCHMARK_BIOSKIN.md")
N_TEST = 300


def lab_under(R, ill):
    E = colour.SDS_ILLUMINANTS[ill].copy().align(SHAPE).values
    cm = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"].copy().align(SHAPE).values.T * E
    cm = cm / cm[1].sum()
    return colour.XYZ_to_Lab(np.atleast_2d(R) @ cm.T, colour.XYZ_to_xy(cm.sum(1)))


def bioskin_spectra(bgr):
    """BGR input (their convention) -> their spectrum resampled on 380-780 / 10 nm."""
    R = D.decode(D.encode(bgr))
    return np.array([np.interp(LAM, LAM_B, r) for r in R]), R


def ours(rgb_lin, lut, fwd):
    out = lut(rgb_lin)
    x = inv.to_unbounded(out["map"])
    S = np.array([fwd.spectrum(xi) for xi in x])
    S_exact = np.array([inv.colour_exact(s, y, lut3d.M_SRGB)[0] for s, y in zip(S, rgb_lin)])
    return S, S_exact, out


def held_out_issa():
    meta, lam_i, R_i = load_issa()
    used = set(np.load(os.path.join(ROOT, "data", "derived", "issa_fits.npz"))["index"].tolist())
    full = np.where(np.all(np.isfinite(R_i[:, (lam_i >= 400) & (lam_i <= 700)]), 1))[0]
    pool = np.array([i for i in full if i not in used])
    test = np.random.default_rng(7).choice(pool, N_TEST, replace=False)  # same draw as V3
    diff = np.clip(R_i[test] - SPECULAR, 0, None)
    return np.array([np.interp(LAM, lam_i[np.isfinite(r)], r[np.isfinite(r)]) for r in diff])


def metrics(S, ref):
    band = (LAM >= 400) & (LAM <= 700)
    rms = np.sqrt(np.mean((S[:, band] - ref[:, band]) ** 2, 1))
    dE = {ill: colour.delta_E(lab_under(ref, ill), lab_under(S, ill), method="CIE 2000") for ill in ILLUMS}
    return rms, dE


def main():
    lut = lut3d.RGBToChromophores(os.path.join(ROOT, "chromophores", "data", "rgb_lut_65_v1.npz"))
    fwd = inv.Forward(specular=0.0)

    # ---------------- 1. ISSA held-out ----------------
    spec = held_out_issa()
    rgb_std = spec @ lut3d.M_SRGB.T  # standard linear sRGB D65
    spec_b = np.array([np.interp(LAM_B, LAM, s) for s in spec])
    rgb_native = theirs_rgb(spec_b, LAM_B)  # their conversion (R, G, B)

    S_ours, S_ours_exact, out = ours(rgb_std, lut, fwd)
    S_b_native, _ = bioskin_spectra(rgb_native[:, ::-1])
    S_b_prod, _ = bioskin_spectra(rgb_std[:, ::-1])
    # Best-effort adaptation: 3x3 matrix sRGB D65 -> their RGB convention, fitted on ISSA spectra
    # disjoint from the test set (the prior-fit subset).
    meta_all, lam_all, R_all = load_issa()
    fit_idx = np.load(os.path.join(ROOT, "data", "derived", "issa_fits.npz"))["index"]
    fit_spec = np.array([np.interp(LAM, lam_all[np.isfinite(r)], r[np.isfinite(r)]) for r in np.clip(R_all[fit_idx] - SPECULAR, 0, None)])
    A_std = fit_spec @ lut3d.M_SRGB.T
    A_their = theirs_rgb(np.array([np.interp(LAM_B, LAM, s_) for s_ in fit_spec]), LAM_B)
    C = np.linalg.lstsq(A_std, A_their, rcond=None)[0]
    S_b_adapt, _ = bioskin_spectra((rgb_std @ C)[:, ::-1])
    methods = {
        "Notre chaîne (biophysique seule)": S_ours,
        "Notre chaîne + correction colorimétrique (production)": S_ours_exact,
        "BioSkin 2023, entrée native (leur conversion)": S_b_native,
        "BioSkin 2023, entrée sRGB D65 (texture livrée)": S_b_prod,
        "BioSkin 2023, sRGB D65 converti dans leur convention (3×3)": S_b_adapt,
    }
    lines = ["| Chaîne | RMS spectral méd. / p95 (400–700) | " + " | ".join(f"ΔE00 {i} méd. / p95" for i in ILLUMS) + " |",
             "|---|---|" + "---|" * len(ILLUMS)]
    res = {}
    for name, S in methods.items():
        rms, dE = metrics(S, spec)
        res[name] = (rms, dE)
        lines.append(f"| {name} | {np.median(rms):.4f} / {np.percentile(rms, 95):.4f} | "
                     + " | ".join(f"{np.median(dE[i]):.2f} / {np.percentile(dE[i], 95):.2f}" for i in ILLUMS) + " |")
    # closure in each pipeline's own colour convention
    rec_native = theirs_rgb(np.array([np.interp(LAM_B, LAM, s) for s in S_b_native]), LAM_B)
    rel = lambda a, b: np.median(np.abs(a / b - 1)) * 100  # noqa: E731

    # ---------------- 2. Emily ----------------
    from emily import load  # noqa: E402

    raw = load("00_diffuse_unlit_raw")
    gains, _, _ = lut3d.autocalibrate(raw)
    d = (raw * gains)[::2, ::2]  # 640^2
    H, W, _ = d.shape
    flat = d.reshape(-1, 3)
    lut_face = lut3d.RGBToChromophores(os.path.join(ROOT, "chromophores", "data", "rgb_lut_face_65_v1.npz"))
    out_e = lut_face(flat)
    skin = out_e["in_gamut"]
    sys.path.insert(0, os.path.join(ROOT, "scripts", "phase2"))
    from emily_reconstruct import spectra_from_maps  # noqa: E402

    rec_ours = spectra_from_maps(out_e["map"]) @ lut3d.M_SRGB.T
    xb = D.encode(flat[:, ::-1])
    Sb = D.decode(xb)
    rec_b_raw = np.array([np.interp(LAM, LAM_B, s) for s in Sb]) @ lut3d.M_SRGB.T
    xb = D.encode((flat @ C)[:, ::-1])
    Sb = D.decode(xb)
    rec_b_std = np.array([np.interp(LAM, LAM_B, s) for s in Sb]) @ lut3d.M_SRGB.T
    Pb = warp(xb)

    def de00(a, b):
        la = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(a, 0, None), apply_cctf_decoding=False))
        lb = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(b, 0, None), apply_cctf_decoding=False))
        return colour.delta_E(la, lb, method="CIE 2000")

    dE_o = de00(flat, rec_ours)[skin]
    dE_b = de00(flat, rec_b_std)[skin]
    dE_b_raw = de00(flat, rec_b_raw)[skin]

    fig, axes = plt.subplots(2, 4, figsize=(17, 8.6))
    panels = [(lut3d.encode(d), "Diffuse Emily (calibrée)"),
              (lut3d.encode(np.where(skin[:, None], rec_ours, 0.8).reshape(H, W, 3)), f"Reconstruite — notre chaîne (ΔE méd. {np.median(dE_o):.2f})"),
              (lut3d.encode(np.where(skin[:, None], rec_b_std, 0.8).reshape(H, W, 3)), f"Reconstruite — BioSkin, entrée adaptée (ΔE méd. {np.median(dE_b):.2f})")]
    for ax, (img, t) in zip(axes[0], panels):
        ax.imshow(img)
        ax.set_title(t, loc="left", fontsize=9)
    im = axes[0, 3].imshow(np.where(skin, de00(flat, rec_b_std), np.nan).reshape(H, W) - np.where(skin, de00(flat, rec_ours), np.nan).reshape(H, W),
                           cmap="RdBu_r", vmin=-3, vmax=3)
    axes[0, 3].set_title("ΔE BioSkin − ΔE nous (rouge : nous meilleurs)", loc="left", fontsize=9)
    plt.colorbar(im, ax=axes[0, 3], fraction=0.04)
    for ax, v, t in zip(axes[1], [out_e["map"][:, 0], Pb[:, 0], out_e["map"][:, 2], Pb[:, 1]],
                        ["Mélanine — nous", "Mélanine — BioSkin", "Sang — nous", "Hémoglobine — BioSkin"]):
        v = np.where(skin, v, np.nan).reshape(H, W)
        lo, hi = np.nanpercentile(v, [2, 98])
        cm = plt.get_cmap("magma").copy()
        cm.set_bad("#cfcfcf")
        im = ax.imshow(v, cmap=cm, vmin=lo, vmax=hi)
        ax.set_title(t, loc="left", fontsize=9)
        plt.colorbar(im, ax=ax, fraction=0.04)
    for ax in axes.ravel():
        ax.axis("off")
    save(fig, "benchmark_bioskin_emily.png")

    i_show = np.argsort(np.sum(spec, 1))[[10, N_TEST // 2, N_TEST - 10]]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, i in zip(axes, i_show):
        ax.plot(LAM, spec[i], color="#0b0b0b", lw=2.2, label="mesuré (ISSA)")
        for (name, S), c in zip(methods.items(), SERIES):
            ax.plot(LAM, S[i], color=c, lw=1.5, label=name)
        ax.set_xlabel("λ (nm)")
    axes[0].set_ylabel("Albédo diffus")
    axes[0].legend(fontsize=7)
    save(fig, "benchmark_bioskin_spectra.png")

    text = f"""# Benchmark réflectance : notre chaîne vs BioSkin (Aliaga et al. 2023)

BioSkin est la base de la partie réflectance d'Aliaga & Jarabo 2026 (leur partie translucidité
n'est pas publiée). Décodeur et encodeur pré-entraînés (`BioSkin.pt`, MIT), exécutés en numpy.

## 1. Spectres ISSA mesurés, non utilisés pour nos a priori ({N_TEST})

Référence : spectre **mesuré** (albédo diffus = SCI − R_sp). Chaque chaîne ne voit que le RGB.

{chr(10).join(lines)}

Fermeture dans leur propre convention (BioSkin natif, RGB reconstruit vs entrée) : erreur relative
médiane {rel(rec_native, rgb_native):.1f} %.

![spectres](../figures/benchmark_bioskin_spectra.png)

## 2. Digital Emily (diffuse cross-pol calibrée, 640²)

Albédo reconstruit à partir des paramètres estimés, comparé à l'entrée (sRGB D65), sur les texels peau :

| Chaîne | ΔE00 médiane | p90 | p95 | % < 1 | % < 2 |
|---|---|---|---|---|---|
| Notre chaîne (cartes seules, sans correction) | {np.median(dE_o):.2f} | {np.percentile(dE_o, 90):.2f} | {np.percentile(dE_o, 95):.2f} | {100 * np.mean(dE_o < 1):.0f} % | {100 * np.mean(dE_o < 2):.0f} % |
| BioSkin 2023, entrée sRGB D65 brute | {np.median(dE_b_raw):.2f} | {np.percentile(dE_b_raw, 90):.2f} | {np.percentile(dE_b_raw, 95):.2f} | {100 * np.mean(dE_b_raw < 1):.0f} % | {100 * np.mean(dE_b_raw < 2):.0f} % |
| BioSkin 2023, entrée convertie dans leur convention (3×3) | {np.median(dE_b):.2f} | {np.percentile(dE_b, 90):.2f} | {np.percentile(dE_b, 95):.2f} | {100 * np.mean(dE_b < 1):.0f} % | {100 * np.mean(dE_b < 2):.0f} % |

![emily](../figures/benchmark_bioskin_emily.png)

Remarque : les grandeurs « mélanine » et « hémoglobine » des deux modèles n'ont pas la même
définition ni la même échelle ; on compare leur structure spatiale, pas leurs valeurs.
"""
    open(REPORT, "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
