"""Phase 2 closure: run the production pipeline on Emily and VFace, measure the effect of the
non-skin mask and of the uncertainty-guided regularisation, export EXR maps.

Outputs: exports/<asset>/*.exr, docs/phase2/CLOTURE.md, docs/figures/phase2_close_*.png
"""

import os
import sys
import time

import colour
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter

HERE = os.path.dirname(__file__)
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)
from common import plt, save  # noqa: E402
from emily import load  # noqa: E402

from chromophores import inversion as inv, lut3d, pipeline as pl  # noqa: E402

Image.MAX_IMAGE_PIXELS = None
ID = {n: i for i, n in enumerate(inv.NAMES)}
VFACE = os.path.join(ROOT, "data", "VFace", "XYZ_albedo_lin_srgb.1001.png")
VFACE_RES = int(os.environ.get("VFACE_RES", "1024"))  # full 4096 locally; 1024 for the repository
VFACE_CROP = (1150, 1700, 640, 640)  # at 4096


def de00(a, b):
    la = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(a, 0, None), apply_cctf_decoding=False))
    lb = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(b, 0, None), apply_cctf_decoding=False))
    return colour.delta_E(la, lb, method="CIE 2000")


def bandpass(v, s):
    return gaussian_filter(v, s) - gaussian_filter(v, 4 * s)


def crosstalk(r, sl):
    """Correlation of log melanin and log blood per spatial band (fine: sigma 2 px, patch: sigma 16 px),
    before / after regularisation; plus the HF std reduction of every parameter (unbounded space)."""
    out = {}
    for label, x in (("avant", r.x_raw), ("après", r.x)):
        P = inv.to_physical(x[sl])
        m, b = np.log(P[..., ID["melanin_fraction"]]), np.log(P[..., ID["blood_fraction"]])
        out[label] = {s: float(np.corrcoef(bandpass(m, s).ravel(), bandpass(b, s).ravel())[0, 1]) for s in (2, 16)}
    hf = lambda x: np.std(x[sl] - gaussian_filter(x[sl], (2, 2, 0)), axis=(0, 1))  # noqa: E731
    out["hf_ratio"] = hf(r.x) / hf(r.x_raw)
    return out


def raw_closure(r):
    """Albedo from the unregularised maps (no residual)."""
    S = pl.spectra_from_params(inv.to_physical(r.x_raw.reshape(-1, inv.N)))
    return (S @ lut3d.M_SRGB.T).reshape(r.albedo.shape)


def figure(r, name, crop_sl, title):
    maps = r.maps
    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    axes[0, 0].imshow(lut3d.encode(r.albedo))
    axes[0, 0].set_title(f"{title} — albédo calibré", loc="left", fontsize=9)
    cm_mask = plt.get_cmap("Greys_r")
    axes[1, 0].imshow(np.where(r.skin, 1.0, np.where(r.in_gamut, 0.45, 0.0)), cmap=cm_mask, vmin=0, vmax=1)
    axes[1, 0].set_title("Masque : blanc = peau, gris = rejeté par l'a priori, noir = hors gamut", loc="left", fontsize=8)
    for j, (n, lab, cmap) in enumerate([("melanin_fraction", "Mélanine", "magma"), ("blood_fraction", "Hémoglobine", "magma"),
                                        ("oxygen_saturation", "SO₂", "viridis")]):
        cm = plt.get_cmap(cmap).copy()
        cm.set_bad("#cfcfcf")
        v = np.where(r.skin, maps[..., ID[n]], np.nan)
        lo, hi = np.nanpercentile(v, [2, 98])
        im = axes[0, j + 1].imshow(v, cmap=cm, vmin=lo, vmax=hi)
        axes[0, j + 1].set_title(lab, loc="left", fontsize=9)
        plt.colorbar(im, ax=axes[0, j + 1], fraction=0.04)
    # crop before / after regularisation
    ax_list = [axes[1, 1], axes[1, 2], axes[1, 3]]
    raw = inv.to_physical(r.x_raw[crop_sl])
    reg = inv.to_physical(r.x[crop_sl])
    m = r.skin[crop_sl]
    for ax, v, t in zip(ax_list, [raw[..., ID["melanin_fraction"]], reg[..., ID["melanin_fraction"]], reg[..., ID["blood_fraction"]]],
                        ["Zoom mélanine — avant régularisation", "Zoom mélanine — après", "Zoom hémoglobine — après"]):
        cm = plt.get_cmap("magma").copy()
        cm.set_bad("#cfcfcf")
        v = np.where(m, v, np.nan)
        lo, hi = np.nanpercentile(v, [1, 99])
        im = ax.imshow(v, cmap=cm, vmin=lo, vmax=hi)
        ax.set_title(t, loc="left", fontsize=9)
        plt.colorbar(im, ax=ax, fraction=0.04)
    for ax in axes.ravel():
        ax.axis("off")
    save(fig, f"phase2_close_{name}.png")


def run(name, albedo, crop_sl, title, export_dir):
    t0 = time.time()
    r = pl.process_albedo(albedo)
    dt = time.time() - t0
    written = pl.export(r, export_dir)
    corr = crosstalk(r, crop_sl)
    dE_raw = de00(r.albedo, raw_closure(r))[r.skin]
    dE_reg = de00(r.albedo, r.rgb_model)[r.skin]
    B = np.abs(lut3d.M_SRGB) / np.abs(lut3d.M_SRGB).max(1, keepdims=True)
    S = pl.spectra_from_params(inv.to_physical(r.x.reshape(-1, inv.N)))
    closed = (S * (1 + r.residual.reshape(-1, 3) @ B)) @ lut3d.M_SRGB.T
    closure_err = float(np.abs(closed - r.albedo.reshape(-1, 3)).max())
    figure(r, name, crop_sl, title)
    return {
        "name": title, "shape": r.albedo.shape[:2], "time": dt, "gains": r.gains, "in_gamut": 100 * r.in_gamut.mean(),
        "skin": 100 * r.skin.mean(), "corr": corr, "dE_raw": dE_raw, "dE_reg": dE_reg, "closure": closure_err,
        "files": [os.path.relpath(p, ROOT) for p in written],
    }


def main():
    rows = []
    em = load("00_diffuse_unlit_raw")
    rows.append(run("emily", em, (slice(470, 770), slice(330, 690)), "Digital Emily", os.path.join(ROOT, "exports", "emily")))
    img = Image.open(VFACE).convert("RGB").resize((VFACE_RES, VFACE_RES), Image.LANCZOS)
    vf = lut3d.decode(np.asarray(img, float) / 255.0)  # the file is sRGB-encoded (see VFACE_1001.md)
    s = VFACE_RES / 4096
    x, y, w, h = (int(v * s) for v in VFACE_CROP)
    rows.append(run("vface", vf, (slice(y, y + h), slice(x, x + w)), "VFace 1001", os.path.join(ROOT, "exports", f"vface_{VFACE_RES}")))

    lines = ["| Asset | Résolution | Temps | Gains (R, G, B) | Gamut | Peau (après masque a priori) | "
             "Corr. mél./sang fine (σ 2 px) avant → après | Corr. mél./sang plaques (σ 16 px) avant → après | ΔE00 cartes seules avant → après régul. (méd. / p95) | Fermeture avec résidu (err. max) |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['name']} | {r['shape'][1]}×{r['shape'][0]} | {r['time']:.0f} s | ({r['gains'][0]:.2f}, {r['gains'][1]:.2f}, {r['gains'][2]:.2f}) | "
                     f"{r['in_gamut']:.0f} % | {r['skin']:.1f} % | {r['corr']['avant'][2]:+.2f} → {r['corr']['après'][2]:+.2f} | "
                     f"{r['corr']['avant'][16]:+.2f} → {r['corr']['après'][16]:+.2f} | "
                     f"{np.median(r['dE_raw']):.2f} / {np.percentile(r['dE_raw'], 95):.2f} → {np.median(r['dE_reg']):.2f} / {np.percentile(r['dE_reg'], 95):.2f} | "
                     f"{r['closure']:.1e} |")
    hf = ["| Asset | " + " | ".join(inv.NAMES) + " |", "|---|" + "---|" * inv.N]
    for r in rows:
        hf.append(f"| {r['name']} | " + " | ".join(f"{v:.2f}" for v in r["corr"]["hf_ratio"]) + " |")
    files = "\n".join(f"- `{f}`" for r in rows for f in r["files"])
    text = f"""# Clôture de la phase 2 : pipeline albédo → chromophores

Module : `chromophores/pipeline.py` (`process_albedo`, `export`). Script : `scripts/phase2/close_phase2.py`
(`VFACE_RES=4096` pour la pleine résolution en local).

Étapes : auto-calibration → LUT 65³ (a priori visage) avec covariance complète par nœud →
masque non-peau (hors gamut, ou chromophores improbables sous l'a priori : distance de
Mahalanobis > χ²₆(99,9 %)) → régularisation guidée par l'incertitude des paramètres secondaires
(σ = 6 px, λ = 1 ; mélanine et sang gardés au MAP) →
modèle direct par texel → résidu couleur (3 coefficients) → export EXR.

{chr(10).join(lines)}

- **Corrélations mélanine / sang** : entre log(mélanine) et log(sang), par bande spatiale
  (différence de gaussiennes σ / 4σ), sur l'extrait de joue (rougeurs de VFace).
- **ΔE00 cartes seules** : albédo recalculé depuis les cartes, sans résidu. Le résidu garantit
  dans tous les cas la fermeture exacte (dernière colonne, précision machine).

### Effet de la régularisation : écart-type haute fréquence après / avant (espace non borné)

{chr(10).join(hf)}

## Constats

1. **La diaphonie mélanine / sang est un phénomène basse fréquence.** Dans les rougeurs de
   VFace, mélanine et sang sont fortement anticorrélés à l'échelle des plaques (σ 16 px), alors
   qu'ils sont positivement corrélés à l'échelle fine (grain de peau, taches : les deux
   assombrissent). Or l'a posteriori par texel corrèle positivement mélanine et sang (≈ +0,4) :
   l'échange observé ne suit **pas** la direction mal déterminée de l'inversion. C'est donc un
   **biais systématique du modèle direct** (la forme spectrale du sang des rougeurs n'est pas
   exactement celle du modèle, et l'ajustement compense en retirant de la mélanine), pas du bruit.
2. **Aucune régularisation locale ne peut le corriger.** Une première version (régularisation
   de tous les paramètres, pondérée par la covariance) le confirmait : diaphonie inchangée, mais
   25 % du détail fin de la mélanine perdu, ce qui est inacceptable pour les artistes. La version
   retenue **garde mélanine et sang au MAP par texel** (détail complet) et ne régularise que les
   paramètres mal contraints (ratio eu/phéo, SO₂, épaisseur, diffusion), conditionnellement aux
   premiers : cartes secondaires moins bruitées, sans toucher aux cartes principales. Seule
   exception : l'échelle de diffusion devient plus bruitée (elle absorbe ce que les autres cèdent).
   Elle n'est de toute façon pas identifiable depuis le RGB (corrélation 0,08, §2 du bilan) :
   **en production, utiliser une valeur constante (a priori) plutôt que la carte.**
3. **Le masque par l'a priori ne rejette presque rien** au-delà du gamut : la LUT tire déjà le MAP
   vers l'a priori, donc la distance de Mahalanobis reste faible. Le masque utile est le gamut
   (cheveux, yeux, intérieur de la bouche sur Emily).
4. **Pistes pour la diaphonie** (phase 3 / modèle direct) : effet de confinement du sang dans les
   vaisseaux (`vessel_radius_cm`, qui aplatit l'absorption de l'hémoglobine), sang superficiel
   dans le derme papillaire, ou a priori spatial explicite (mélanine « haute fréquence », sang
   « basse fréquence ») validé contre une capture multispectrale.

![Emily](../figures/phase2_close_emily.png)
![VFace](../figures/phase2_close_vface.png)

## Fichiers exportés (EXR half, ZIP)

{files}

Cartes : mélanine, ratio eu/phéo, hémoglobine, oxygénation, épaisseur d'épiderme (cm), échelle de
diffusion, chacune avec son incertitude relative (`*_sigma`, ≈ σ relatif à 1 écart-type) ;
`skin_mask` ; `albedo_input_calibrated` ; `albedo_from_maps` (modèle seul) ;
`colour_residual_coeffs` (3 coefficients : spectre final = spectre des cartes × (1 + Bᵀc), B = lobes
des fonctions colorimétriques, voir `inversion.colour_exact`).
"""
    path = os.path.join(ROOT, "docs", "phase2", "CLOTURE.md")
    open(path, "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
