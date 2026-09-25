"""Vessel packaging of blood (van Veen et al. 2002): does it remove the systematic
haemoglobin-band residual seen on ISSA and NIST? Baseline scale fixed at 0.5."""

import os
from multiprocessing import Pool

import numpy as np
from common1 import SERIES, plt, save, write_report
from issa_baseline import NAMES, fit_one

from chromophores.datasets import ISSA_FACE, load_issa, load_nist

RADII_UM = [0, 10, 25, 50, 100]


def run(R, lam, label):
    res = {}
    with Pool(os.cpu_count()) as pool:
        for r_um in RADII_UM:
            out = pool.map(fit_one, [(r, lam, 0.5, {"vessel_radius_cm": r_um * 1e-4}) for r in R], chunksize=4)
            res[r_um] = (np.array([o[0] for o in out]), np.array([o[1] for o in out]), np.array([o[2] for o in out]))
            print(f"{label} radius {r_um} um done", flush=True)
    return res


def main():
    meta, lam_i, R_i = load_issa()
    sel = (lam_i >= 400) & (lam_i <= 700)
    ok = np.all(np.isfinite(R_i[:, sel]), 1)
    rng = np.random.default_rng(1)
    idx = np.concatenate([rng.choice(g.index.to_numpy(), min(len(g), 60), replace=False) for _, g in meta[ok].groupby("ethnicity")])
    lam_n, R_n = load_nist()
    lam_nist = np.arange(400.0, 1000.0, 10.0)
    runs = {
        "ISSA 400–700 nm": (lam_i[sel], R_i[idx][:, sel]),
        "NIST 400–1000 nm": (lam_nist, np.array([np.interp(lam_nist, lam_n, r) for r in R_n])),
    }
    si = NAMES.index("scattering_scale")
    text = "# Empaquetage vasculaire du sang (fond × 0,5)\n\n"
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, (label, (lam, R)) in zip(axes, runs.items()):
        res = run(R, lam, label)
        text += f"## {label} ({len(R)} spectres)\n\n| Rayon des vaisseaux | RMS médian | p95 | μs′ médian | sang médian |\n|---|---|---|---|---|\n"
        for (r_um, (rms, r, P)), c in zip(res.items(), SERIES):
            text += f"| {r_um} µm | {np.median(rms):.4f} | {np.percentile(rms, 95):.4f} | {np.median(P[:, si]):.2f} | {np.median(P[:, NAMES.index('blood_fraction')]):.4f} |\n"
            ax.plot(lam, np.median(r, 0), color=c, label=f"{r_um} µm")
        text += "\n"
        ax.axhline(0, color="#52514e", lw=1)
        ax.set_title(f"{label} : résidu médian", loc="left")
        ax.set_xlabel("λ (nm)")
    axes[0].set_ylabel("Résidu modèle − mesure")
    axes[0].legend(fontsize=8, title="rayon vaisseaux")
    save(fig, "packaging_study.png")
    write_report("PACKAGING_STUDY.md", text)


if __name__ == "__main__":
    main()
