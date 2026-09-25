"""V3 (ADR-0005): RGB -> chromophores inversion, validated on held-out measured spectra.

(a) LUT vs direct MAP optimisation (same prior).
(b) Measured ISSA spectra (not used to build the prior) -> cross-polarised RGB -> LUT,
    compared with the fit of the full spectrum (best available "truth" without in-vivo ground truth).
(c) Colour reconstruction error.
(d) Spectral uplift: spectrum predicted from RGB (biophysical) vs measured spectrum,
    compared with generic RGB -> spectrum methods (Jakob & Hanika 2019, Mallett & Yuksel 2019),
    including the colour error under other illuminants (what a spectral renderer sees).
"""

import os
import sys
from multiprocessing import Pool

import colour
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from common import SERIES, plt, save  # noqa: E402

from chromophores import inversion as inv, lut3d  # noqa: E402
from chromophores.datasets import load_issa  # noqa: E402
from chromophores.forward import SPECULAR  # noqa: E402

REPORT = os.path.join(HERE, "..", "..", "docs", "phase2")
os.makedirs(REPORT, exist_ok=True)
N_TEST = 300
LAM = inv.LAMBDA
ID = {n: i for i, n in enumerate(inv.NAMES)}
SHAPE = colour.SpectralShape(380, 780, 10)
ILLUMS = ["D65", "A", "FL11", "LED-B3"]
_W = {}


def worker_init():
    _W["fwd_rgb"] = inv.Forward(specular=0.0)
    _W["prior"] = inv.load_priors()["population"]
    _W["broad"] = inv.Prior(np.zeros(inv.N), 16.0 * np.eye(inv.N))


def direct_rgb(y):
    p = _W["prior"]
    starts = [p.mean, p.mean + np.array([2.0, 0, 0, 0, 0, 0]), p.mean + np.array([-2.0, 0, 1.0, 0, 0, 0])]
    return inv.map_estimate(y, lut3d.M_SRGB, _W["fwd_rgb"], p, lut3d.noise_model(y), starts).x


def spectral_fit(args):
    r, lam = args
    key = ("spec", tuple(lam))
    if key not in _W:
        _W[key] = inv.Forward(lam=lam, specular=0.0)
    starts = [np.array([m, 0.5, 0.0, 0.5, 0.0, 0.0]) for m in (-2.0, 0.0, 2.0)]
    est = inv.map_estimate(r, np.eye(len(lam)), _W[key], _W["broad"], 0.005, starts)
    return est.x, np.sqrt(np.diag(est.cov))


def lab_under(R, ill):
    sd = colour.SDS_ILLUMINANTS[ill].copy().align(SHAPE)
    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"].copy().align(SHAPE)
    E = sd.values
    M = cmfs.values.T * E
    M = M / M[1].sum()
    white = colour.XYZ_to_xy(M.sum(1))
    return colour.XYZ_to_Lab(np.atleast_2d(R) @ M.T, white)


def main():
    lut = lut3d.RGBToChromophores(os.environ.get("LUT_PATH"))
    meta, lam_i, R_i = load_issa()
    used = set(np.load(os.path.join(HERE, "..", "..", "data", "derived", "issa_fits.npz"))["index"].tolist())
    full = np.where(np.all(np.isfinite(R_i[:, (lam_i >= 400) & (lam_i <= 700)]), 1))[0]
    pool_idx = np.array([i for i in full if i not in used])
    rng = np.random.default_rng(7)
    test = rng.choice(pool_idx, N_TEST, replace=False)
    sel = (lam_i >= 400) & (lam_i <= 700)
    lam_fit = lam_i[sel]
    # Measured cross-polarised albedo spectrum on the 380-780 grid (edges extended) and on the fit range.
    diff = np.clip(R_i[test] - SPECULAR, 0, None)
    spec380 = np.array([np.interp(LAM, lam_i[np.isfinite(r)], r[np.isfinite(r)]) for r in diff])
    rgb = spec380 @ lut3d.M_SRGB.T

    out = lut(rgb)
    x_lut = inv.to_unbounded(out["map"])
    s_lut = inv.to_unbounded(out["p84"]) - x_lut
    with Pool(os.cpu_count(), initializer=worker_init) as pool:
        x_dir = np.array(pool.map(direct_rgb, list(rgb), chunksize=8))
        fits = pool.map(spectral_fit, [(d[sel], lam_fit) for d in diff], chunksize=8)
    x_spec = np.array([f[0] for f in fits])
    s_spec = np.array([f[1] for f in fits])
    ing = out["in_gamut"]

    # (a) LUT vs direct MAP, in units of posterior std
    za = np.abs(x_lut - x_dir) / np.maximum(s_lut, 1e-6)
    P_lut, P_dir, P_spec = inv.to_physical(x_lut), inv.to_physical(x_dir), inv.to_physical(x_spec)
    rel_a = {n: np.abs(P_lut[:, ID[n]] / P_dir[:, ID[n]] - 1) for n in ("melanin_fraction", "blood_fraction")}

    # (b) RGB inversion vs full-spectrum fit
    zb = (x_lut - x_spec) / np.sqrt(s_lut**2 + s_spec**2)
    lines_b = ["| Paramètre | corrélation (log) RGB vs spectre | erreur relative médiane | |z| médian | % |z| < 2 | σ LUT médian / σ spectre médian |",
               "|---|---|---|---|---|---|"]
    for n in inv.NAMES:
        k = ID[n]
        a, b = P_lut[ing, k], P_spec[ing, k]
        corr = np.corrcoef(np.log(a), np.log(b))[0, 1]
        rel = np.median(np.abs(a / b - 1))
        lines_b.append(f"| {n} | {corr:.2f} | {100 * rel:.0f} % | {np.median(np.abs(zb[ing, k])):.2f} | {100 * np.mean(np.abs(zb[ing, k]) < 2):.0f} % | "
                       f"{np.median(s_lut[ing, k]):.2f} / {np.median(s_spec[ing, k]):.2f} |")

    # (c) colour reconstruction: spectrum at LUT MAP -> RGB vs input
    fwd = inv.Forward(specular=0.0)
    rec = np.array([fwd.spectrum(x) for x in x_lut])
    dE_c = colour.delta_E(lab_under(spec380, "D65"), lab_under(rec, "D65"), method="CIE 2000")

    # (d) spectral uplift quality vs generic methods
    xyz_in = colour.sRGB_to_XYZ(np.clip(rgb, 1e-6, None), apply_cctf_decoding=False)
    jak = np.array([colour.XYZ_to_sd(x, method="Jakob 2019").copy().align(SHAPE).values for x in xyz_in])
    mal = np.array([colour.recovery.RGB_to_sd_Mallett2019(np.clip(c, 0, 1)).copy().align(SHAPE).values for c in rgb])
    rec_exact = np.array([inv.colour_exact(r_, y_, lut3d.M_SRGB)[0] for r_, y_ in zip(rec, rgb)])
    methods = {"biophysique (notre LUT)": rec, "biophysique + correction colorimétrique": rec_exact,
               "Jakob & Hanika 2019": jak, "Mallett & Yuksel 2019": mal}
    fit_band = (LAM >= 400) & (LAM <= 700)
    lines_d = ["| Méthode RGB → spectre | RMS spectral médian (400–700) | p95 | " + " | ".join(f"ΔE00 sous {i} méd. / p95" for i in ILLUMS) + " |",
               "|---|---|---|" + "---|" * len(ILLUMS)]
    for name, S in methods.items():
        rms = np.sqrt(np.mean((S[:, fit_band] - spec380[:, fit_band]) ** 2, 1))[ing]
        cells = []
        for ill in ILLUMS:
            d = colour.delta_E(lab_under(spec380[ing], ill), lab_under(S[ing], ill), method="CIE 2000")
            cells.append(f"{np.median(d):.2f} / {np.percentile(d, 95):.2f}")
        lines_d.append(f"| {name} | {np.median(rms):.4f} | {np.percentile(rms, 95):.4f} | " + " | ".join(cells) + " |")

    crit_a = all(np.median(v[ing]) < 0.01 for v in rel_a.values())
    text = f"""# V3 : inversion RGB → chromophores, sur spectres mesurés ISSA

- {N_TEST} spectres ISSA **non utilisés** pour construire l'a priori ; {int(ing.sum())} dans le gamut de la LUT.
- Entrée : albédo en polarisation croisée (SCI − R_sp) rendu en sRGB linéaire D65. LUT 33³, a priori « population ».
- Bruit supposé sur le RGB : 0,004 + 2 % (calibration).

## (a) LUT vs optimisation MAP directe

| | médiane | p95 |
|---|---|---|
| erreur relative mélanine | {100 * np.median(rel_a['melanin_fraction'][ing]):.2f} % | {100 * np.percentile(rel_a['melanin_fraction'][ing], 95):.2f} % |
| erreur relative sang | {100 * np.median(rel_a['blood_fraction'][ing]):.2f} % | {100 * np.percentile(rel_a['blood_fraction'][ing], 95):.2f} % |
| écart / σ a posteriori (tous paramètres) | {np.median(za[ing]):.3f} | {np.percentile(za[ing], 95):.3f} |

Critère ADR-0005 (écart LUT / optimisation < 1 % sur les paramètres identifiables) → **{"OK" if crit_a else "NON ATTEINT"}**.

## (b) Ce que le RGB permet de retrouver (vs ajustement du spectre complet 400–700 nm)

{chr(10).join(lines_b)}

z = (x_RGB − x_spectre) / √(σ²_RGB + σ²_spectre) dans l'espace non borné : |z| < 2 pour ~95 % des cas si les incertitudes sont bien calibrées.

## (c) Reconstruction de la couleur

ΔE00 (D65) entre l'albédo d'entrée et l'albédo reconstruit : médiane {np.median(dE_c[ing]):.2f}, p95 {np.percentile(dE_c[ing], 95):.2f}, max {dE_c[ing].max():.2f}.
Critère ADR-0005 (ΔE00 < 1) → **{"OK" if np.percentile(dE_c[ing], 95) < 1 else "NON ATTEINT (p95)"}**.

## (d) Spectre reconstruit depuis le RGB vs spectre mesuré

{chr(10).join(lines_d)}
"""
    path = os.path.join(REPORT, os.environ.get("V3_REPORT", "V3_inversion.md"))
    open(path, "w").write(text)
    print(text)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, n, lab in zip(axes, ["melanin_fraction", "blood_fraction", "scattering_scale"], ["Mélanine", "Sang", "Échelle de μs′"]):
        k = ID[n]
        ax.loglog(P_spec[ing, k], P_lut[ing, k], "o", ms=3.5, color=SERIES[0], alpha=0.6, mec="none")
        lo, hi = inv.BOUNDS[n][:2]
        ax.plot([lo, hi], [lo, hi], color="#52514e", lw=1)
        ax.set_xlabel(f"{lab} — ajustement du spectre complet")
        ax.set_ylabel(f"{lab} — depuis le RGB (LUT)")
        ax.set_title(lab, loc="left")
    save(fig, "phase2_v3_rgb_vs_spectrum.png")
    i_show = np.argsort(np.sum(spec380, 1))[[5, len(test) // 2, len(test) - 5]]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=False)
    for ax, i in zip(axes, i_show):
        ax.plot(LAM, spec380[i], color="#0b0b0b", lw=2.2, label="mesuré (ISSA)")
        for (name, S), c in zip(methods.items(), SERIES):
            ax.plot(LAM, S[i], color=c, lw=1.6, label=name)
        ax.set_xlabel("λ (nm)")
    axes[0].set_ylabel("Albédo diffus")
    axes[0].legend(fontsize=7)
    fig.suptitle("Spectre reconstruit depuis le seul RGB : exemples (peau foncée, moyenne, claire)", x=0.01, ha="left")
    save(fig, "phase2_v3_uplift_examples.png")


if __name__ == "__main__":
    main()
