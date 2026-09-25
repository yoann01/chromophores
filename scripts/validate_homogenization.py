"""Validate the dimensionless homogenisation table.

1. Table vs direct Monte Carlo on random (skin, wavelength) pairs:
   multiple-scattering albedo A and mean exit radius.
2. Radial profile of the layered skin vs the homogeneous random-walk proxy
   fitted by moment matching (what the renderer would actually trace).
"""

import numpy as np
from common import SERIES, plt, save
from identifiability import RANGES, from_unit

from chromophores import montecarlo, optics
from chromophores.homogenize import SkinHomogenizer
from chromophores.optics import SkinParams

H = SkinHomogenizer()
DERMIS_CM = 5.0  # effectively semi-infinite, consistent with the table


def layered_mc(p, lam, n_photons=100_000, r_max=1.0, n_bins=200, seed=5):
    epi, der = optics.layers(np.array([lam]), p, DERMIS_CM)
    return montecarlo.run_layers(
        [epi.mua[0], der.mua[0]], [epi.mus[0], der.mus[0]], [epi.g, der.g], [epi.thickness_cm, der.thickness_cm], n_photons, r_max, n_bins, seed
    )


def table_accuracy(n=40, seed=11):
    rng = np.random.default_rng(seed)
    errA, errR = [], []
    for k in range(n):
        p = from_unit(rng.random(len(RANGES)))
        lam = float(rng.uniform(400, 780))
        A, rbar, _ = layered_mc(p, lam, seed=100 + k)
        m = H.medium(np.array([lam]), p)
        errA.append(m.albedo_ms[0] / A - 1)
        errR.append(m.mean_radius_cm[0] / rbar - 1)
    errA, errR = np.abs(errA), np.abs(errR)
    print("\n| Grandeur | Erreur relative médiane | p90 | max |\n|---|---|---|---|")
    for name, e in [("Albédo multiple A", errA), ("Rayon moyen de sortie", errR)]:
        print(f"| {name} | {100 * np.median(e):.1f} % | {100 * np.percentile(e, 90):.1f} % | {100 * e.max():.1f} % |")


def profiles():
    cases = [("moyenne", SkinParams(melanin_fraction=0.08)), ("foncée", SkinParams(melanin_fraction=0.30, eumelanin_ratio=0.9))]
    lams = [450.0, 550.0, 650.0, 750.0]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    print("\n| Peau | λ (nm) | A | rayon moyen (mm) | α homogène | 1/σt (mm) | écart L1 de r·R(r) |\n|---|---|---|---|---|---|---|")
    n_bins, r_max = 200, 1.0
    edges = np.linspace(0, r_max, n_bins + 1)
    rc = 0.5 * (edges[1:] + edges[:-1])
    areas = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    for ax, (name, p) in zip(axes, cases):
        for lam, c in zip(lams, SERIES):
            A, rbar, h = layered_mc(p, lam, n_bins=n_bins, r_max=r_max)
            m = H.medium(np.array([lam]), p)
            st, ss = m.sigma_t[0], m.sigma_s[0]
            Ah, rh, hh = montecarlo.run_layers([st - ss], [ss], [m.g], [1e3 / st], 100_000, r_max, n_bins, seed=77)
            f_l, f_h = rc * h / areas, rc * hh / areas
            l1 = np.abs(f_l - f_h).sum() / f_l.sum()
            print(f"| {name} | {lam:.0f} | {A:.3f} | {10 * rbar:.2f} | {ss / st:.5f} | {10 / st:.3f} | {100 * l1:.0f} % |")
            ax.semilogy(10 * rc, np.maximum(h / areas, 1e-5), color=c, label=f"{lam:.0f} nm multicouche")
            ax.semilogy(10 * rc, np.maximum(hh / areas, 1e-5), color=c, ls="--", lw=1.5, label=f"{lam:.0f} nm homogène")
        ax.set_xlim(0, 5)
        ax.set_xlabel("r (mm)")
        ax.set_title(f"Peau {name} : multicouche vs milieu homogène équivalent", loc="left", fontsize=10)
    axes[0].set_ylabel("R(r) (cm⁻²)")
    axes[0].legend(fontsize=7, ncol=2)
    save(fig, "homogenization_profiles.png")


if __name__ == "__main__":
    table_accuracy()
    profiles()
