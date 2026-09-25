"""Same arbitration as issa_baseline.py on the NIST reference spectra (100 forearms).

NIST extends to the NIR, which also constrains scattering and water: fits are
run on 400-700 nm (comparable with ISSA) and 400-1000 nm.
"""

import os
from multiprocessing import Pool

import numpy as np
from common1 import SERIES, plt, save, write_report
from issa_baseline import BASELINES, NAMES, fit_one

from chromophores.datasets import load_nist


def main():
    lam_n, R_n = load_nist()
    runs = {}
    for name, lam in (("400–700 nm", np.arange(400.0, 701.0, 10.0)), ("400–1000 nm", np.arange(400.0, 1000.0, 10.0))):
        R = np.array([np.interp(lam, lam_n, r) for r in R_n])
        res = {}
        with Pool(os.cpu_count()) as pool:
            for bs in BASELINES:
                out = pool.map(fit_one, [(r, lam, bs) for r in R], chunksize=4)
                res[bs] = (np.array([o[0] for o in out]), np.array([o[1] for o in out]), np.array([o[2] for o in out]))
                print(f"{name} baseline {bs}: done", flush=True)
        runs[name] = (lam, res)

    si = NAMES.index("scattering_scale")
    text = "# Absorption de fond : arbitrage par les 100 spectres NIST (avant-bras)\n\n"
    text += "Même protocole que ISSA (6 paramètres libres, μs′ ∈ [0,3 ; 3], SCI = R_sp + (1 − R_sp)·A).\n\n"
    for name, (lam, res) in runs.items():
        text += f"## {name}\n\n| Échelle du fond | RMS résidu médian | p95 | % RMS < 0,01 | μs′ fitté médian | % μs′ hors [0,6 ; 1,6] | % μs′ en butée |\n|---|---|---|---|---|---|---|\n"
        for bs, (rms, _, P) in res.items():
            sc = P[:, si]
            text += (f"| {bs:g} | {np.median(rms):.4f} | {np.percentile(rms, 95):.4f} | {100 * np.mean(rms < 0.01):.0f} % | {np.median(sc):.2f} | "
                     f"{100 * np.mean((sc < 0.6) | (sc > 1.6)):.0f} % | {100 * np.mean((sc < 0.31) | (sc > 2.9)):.0f} % |\n")
        text += "\n"
    write_report("NIST_BASELINE.md", text)

    lam, res = runs["400–1000 nm"]
    fig, ax = plt.subplots(figsize=(7, 4))
    for (bs, (rms, r, P)), c in zip(res.items(), SERIES):
        ax.plot(lam, np.median(r, 0), color=c, label=f"fond × {bs:g}")
        ax.fill_between(lam, np.percentile(r, 25, 0), np.percentile(r, 75, 0), color=c, alpha=0.12, lw=0)
    ax.axhline(0, color="#52514e", lw=1)
    ax.set_xlabel("λ (nm)")
    ax.set_ylabel("Résidu modèle − mesure (médiane, IQR)")
    ax.set_title("NIST 400–1000 nm : résidu spectral systématique", loc="left")
    ax.legend(fontsize=8)
    save(fig, "nist_baseline.png")


if __name__ == "__main__":
    main()
