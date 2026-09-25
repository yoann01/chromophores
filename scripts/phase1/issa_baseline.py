"""Arbitrate the baseline (bloodless tissue) absorption with measured skin spectra.

Data: International Skin Spectra Archive (ISSA / Leeds Skin Database), 15k SCI
spectrophotometer spectra, 400-700 nm. A stratified subsample (per ethnic
group, face and body sites) is fitted with the phase-1 forward model for
several baseline scales; scattering is left free over a wide range (0.3-3)
so that a wrong baseline can only be compensated by an implausible mu_s'.

Measurement model: SCI = specular (Fresnel at 8 deg, n = 1.4) + diffuse albedo.
"""

import os
import sys
from multiprocessing import Pool

import numpy as np
from common1 import SERIES, plt, save, write_report
from scipy.optimize import least_squares

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from chromophores import table  # noqa: E402
from chromophores.datasets import ISSA_FACE, load_issa  # noqa: E402
from chromophores.forward import SPECULAR  # noqa: E402
from chromophores.optics import SkinModel, SkinParams  # noqa: E402

BASELINES = [1.0, 0.5, 0.25, 0.0]
PER_GROUP = 120
BOUNDS = {  # name: (low, high, log)
    "melanin_fraction": (0.001, 0.6, True),
    "eumelanin_ratio": (0.0, 1.0, False),
    "blood_fraction": (0.0005, 0.2, True),
    "oxygen_saturation": (0.3, 1.0, False),
    "epidermis_thickness_cm": (0.002, 0.035, True),
    "scattering_scale": (0.3, 3.0, True),
}
NAMES = list(BOUNDS)
_T = {}


def to_params(x):
    kw = {}
    for v, (name, (lo, hi, lg)) in zip(x, BOUNDS.items()):
        s = 1 / (1 + np.exp(-v))
        kw[name] = float(np.exp(np.log(lo) + s * (np.log(hi) - np.log(lo)))) if lg else lo + s * (hi - lo)
    return SkinParams(bilirubin_umol=0.0, **kw)


def model_sci(x, lam, model):
    p = to_params(x)
    tau, ad, de, g, _ = table.dimensionless(lam, p, model)
    T = _T["table"]
    return SPECULAR + (1 - SPECULAR) * np.exp(T._A(T.coords(tau, ad, de, g)))


def fit_one(args):
    R, lam, bs = args[:3]
    extra = args[3] if len(args) > 3 else {}
    if "table" not in _T:
        _T["table"] = table.LayeredTable(method="linear")  # speed; relative comparison between baselines
    model = SkinModel(baseline_scale=bs, **extra)
    best = None
    for mel0 in (-2.0, 0.0, 2.0):
        for sc0 in (-1.0, 0.5):
            x0 = np.array([mel0, 0.5, 0.0, 0.5, 0.0, sc0])
            sol = least_squares(lambda x: model_sci(x, lam, model) - R, x0, max_nfev=300)
            if best is None or sol.cost < best.cost:
                best = sol
    p = to_params(best.x)
    return np.sqrt(np.mean(best.fun**2)), best.fun, [getattr(p, n) for n in NAMES]


def main():
    meta, lam_all, R_all = load_issa()
    sel = (lam_all >= 400) & (lam_all <= 700)
    lam = lam_all[sel]
    R = R_all[:, sel]
    ok = np.all(np.isfinite(R), 1)
    rng = np.random.default_rng(0)
    idx = []
    for eth, grp in meta[ok].groupby("ethnicity"):
        face = grp[grp.location.isin(ISSA_FACE)].index.to_numpy()
        body = grp[~grp.location.isin(ISSA_FACE)].index.to_numpy()
        idx += list(rng.choice(face, min(len(face), PER_GROUP // 2), replace=False))
        idx += list(rng.choice(body, min(len(body), PER_GROUP // 2), replace=False))
    idx = np.array(idx)
    sub_meta = meta.loc[idx].reset_index(drop=True)
    results = {}
    with Pool(os.cpu_count()) as pool:
        for bs in BASELINES:
            out = pool.map(fit_one, [(R[i], lam, bs) for i in idx], chunksize=8)
            results[bs] = (np.array([o[0] for o in out]), np.array([o[1] for o in out]), np.array([o[2] for o in out]))
            print(f"baseline {bs}: done", flush=True)

    si = NAMES.index("scattering_scale")
    lines = ["| Échelle du fond | RMS résidu médian | p95 | % RMS < 0,01 | μs′ fitté médian | % μs′ hors [0,6 ; 1,6] | % μs′ en butée |",
             "|---|---|---|---|---|---|---|"]
    for bs, (rms, res, P) in results.items():
        sc = P[:, si]
        out_plaus = np.mean((sc < 0.6) | (sc > 1.6)) * 100
        at_bound = np.mean((sc < 0.31) | (sc > 2.9)) * 100
        lines.append(f"| {bs:g} | {np.median(rms):.4f} | {np.percentile(rms, 95):.4f} | {100 * np.mean(rms < 0.01):.0f} % | "
                     f"{np.median(sc):.2f} | {out_plaus:.0f} % | {at_bound:.0f} % |")
    eth_lines = ["| Groupe | n | " + " | ".join(f"RMS méd. (fond {bs:g})" for bs in BASELINES) + " | " +
                 " | ".join(f"μs′ méd. (fond {bs:g})" for bs in BASELINES) + " |",
                 "|---|---|" + "---|" * (2 * len(BASELINES))]
    for eth in sorted(sub_meta.ethnicity.unique()):
        m = (sub_meta.ethnicity == eth).to_numpy()
        eth_lines.append(f"| {eth} | {m.sum()} | " + " | ".join(f"{np.median(results[bs][0][m]):.4f}" for bs in BASELINES) + " | " +
                         " | ".join(f"{np.median(results[bs][2][m, si]):.2f}" for bs in BASELINES) + " |")
    text = f"""# Absorption de fond : arbitrage par les spectres mesurés ISSA

- Données : ISSA (Leeds Skin Database), spectrophotomètre SCI, 400–700 nm / 10 nm.
- Échantillon stratifié : {len(idx)} spectres ({PER_GROUP // 2} visage + {PER_GROUP // 2} corps max. par groupe ethnique).
- Modèle de mesure : SCI = R_sp (Fresnel, n = 1,4) + (1 − R_sp)·A_table. Carotène 0,4 µM, bilirubine 0.
- 6 paramètres ajustés par spectre (mélanine, β, sang, SO₂, épaisseur, échelle de μs′ ∈ [0,3 ; 3]), multi-départs.
- Plage « plausible » de l'échelle de μs′ (Jacques 2013, dispersion inter-individuelle) : [0,6 ; 1,6].

{chr(10).join(lines)}

## Par groupe

{chr(10).join(eth_lines)}

Codes : CA caucasien, CN chinois, SA sud-asiatique, AF africain, IQ irakien, TH thaï, JP japonais, AB arabe.
"""
    write_report("ISSA_BASELINE.md", text)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for (bs, (rms, res, P)), c in zip(results.items(), SERIES):
        axes[0].plot(lam, np.median(res, 0), color=c, label=f"fond × {bs:g}")
        axes[0].fill_between(lam, np.percentile(res, 25, 0), np.percentile(res, 75, 0), color=c, alpha=0.12, lw=0)
        axes[1].hist(P[:, si], bins=np.geomspace(0.3, 3, 30), histtype="step", color=c, lw=2, label=f"fond × {bs:g}")
    axes[0].axhline(0, color="#52514e", lw=1)
    axes[0].set_xlabel("λ (nm)")
    axes[0].set_ylabel("Résidu modèle − mesure (médiane, IQR)")
    axes[0].set_title("Résidu spectral systématique", loc="left")
    axes[0].legend(fontsize=8)
    axes[1].set_xscale("log")
    axes[1].axvspan(0.6, 1.6, color="#e4e3df", alpha=0.6, lw=0)
    axes[1].set_xlabel("Échelle de μs′ ajustée")
    axes[1].set_ylabel("Nombre de spectres")
    axes[1].set_title("Diffusion nécessaire (zone grise : plausible)", loc="left")
    save(fig, "issa_baseline.png")


if __name__ == "__main__":
    main()
