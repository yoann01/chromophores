"""Figure 2: analytic (diffusion) model vs Monte Carlo random walk.

Shows (a) the bias of cheap analytic models commonly used for chromophore
estimation and (b) the per-wavelength radial profile R(r) that defines the
"scattering radius" a renderer needs, which an RGB albedo does not contain.
"""

import numpy as np
from common import SERIES, plt, save

from chromophores import montecarlo, reflectance
from chromophores.optics import SkinParams

SKINS = {
    "claire (Cm=2 %)": SkinParams(melanin_fraction=0.02, blood_fraction=0.015),
    "moyenne (Cm=8 %)": SkinParams(melanin_fraction=0.08),
    "foncée (Cm=30 %)": SkinParams(melanin_fraction=0.30, eumelanin_ratio=0.9),
}
lam_mc = np.arange(400.0, 781.0, 20.0)
lam = np.arange(380.0, 781.0, 5.0)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
rows = []
for (name, p), c in zip(SKINS.items(), SERIES):
    ra = reflectance.diffuse_reflectance(lam, p)
    rm, r, prof = montecarlo.simulate(lam_mc, p, n_photons=100_000)
    ax.plot(lam, ra, color=c, label=f"{name} — analytique")
    ax.plot(lam_mc, rm, "o", color=c, ms=5, mec="#fcfcfb", mew=1, label=f"{name} — Monte Carlo")
    rel = np.abs(np.interp(lam_mc, lam, ra) - rm) / rm
    rows.append((name, rel.mean(), rel.max()))
    if name.startswith("moyenne"):
        prof_mid, r_mid = prof, r
ax.set_xlabel("Longueur d'onde (nm)")
ax.set_ylabel("Réflectance diffuse")
ax.set_title("Modèle analytique vs marche aléatoire", loc="left")
ax.legend(fontsize=7)

ax = axes[1]
for k, (wl, c) in enumerate(zip([440, 560, 660, 760], SERIES)):
    i = int(np.argmin(np.abs(lam_mc - wl)))
    ax.semilogy(r_mid * 10, np.maximum(prof_mid[i], 1e-6), color=c, label=f"{lam_mc[i]:.0f} nm")
ax.set_xlabel("Distance radiale r (mm)")
ax.set_ylabel("R(r) (cm⁻²)")
ax.set_xlim(0, 6)
ax.set_title("Profil radial (peau moyenne) — le « rayon de diffusion »", loc="left")
ax.legend(fontsize=8)
save(fig, "forward_model_validation.png")

print("\n| Peau | Écart relatif moyen | Écart relatif max |\n|---|---|---|")
for name, m, x in rows:
    print(f"| {name} | {100 * m:.1f} % | {100 * x:.1f} % |")
