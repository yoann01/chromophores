"""V2 (ADR-0005): K=2 mixture vs two-layer skin.

(a) Fit quality over all table cells (K=1 vs K=2).
(b) Independent check: random (skin, wavelength) pairs, direct Monte Carlo
    profile vs the profile of the *interpolated* mixture (what the engine gets).
(c) Albedo colour error: full 380-780 nm spectra of random skins, direct Monte
    Carlo albedo vs mixture albedo (Delta E00, D65).
(d) Figure: profiles of representative skins, MC vs K=1 vs K=2.
"""

import colour
import numpy as np
from common1 import SERIES, plt, random_skin, save, write_report

from chromophores import color, table
from chromophores.forward import SPECULAR, SkinForward
from chromophores.mixture import MixtureModel
from chromophores.optics import SkinParams
from chromophores.spectra import WAVELENGTHS

F = SkinForward()
MM = MixtureModel(F.lut)
EDGES = table.RHO_EDGES
CENTERS = table.rho_centers(EDGES)
AREAS = table.bin_areas(EDGES)


def mc_profile(lam, p, n_photons=100_000, seed=1):
    tau_e, a_d, d_e, g, musp = (float(v[0]) for v in table.dimensionless(np.array([lam]), p))
    A, r1, r2, hist = table.simulate_cell(tau_e, a_d, d_e, g, n_photons, seed)
    return A, hist, musp


def mixture_log_profile(lam, p, K):
    m = F.media(np.array([lam]), p, K)
    musp = table.dimensionless(np.array([lam]), p)[4][0]
    comps = [(np.log(1 - m.alpha[k, 0]), m.g[k, 0], m.sigma_t[k, 0] / musp) for k in range(K)]
    return MM.evaluate(comps, m.weight[:, 0], CENTERS)


def log_rmse(A, hist, logP, energy=False):
    valid = (hist / A > 1e-5) & (hist > 0) & np.isfinite(logP)
    ref = np.log(hist[valid] / AREAS[valid] / A)
    d2 = (logP[valid] - ref) ** 2
    if energy:
        w = hist[valid] / hist[valid].sum()
        return float(np.sqrt(np.sum(w * d2)))
    return float(np.sqrt(np.mean(d2)))


def noise_floor(A, hist, A2, hist2, energy=False):
    """Per-run log-RMSE between two independent MC runs of the same cell (/sqrt 2)."""
    both = (hist > 0) & (hist2 > 0)
    logP2 = np.full(len(hist), np.nan)
    logP2[both] = np.log(hist2[both] / AREAS[both] / A2)
    return log_rmse(A, np.where(both, hist, 0), logP2, energy) / np.sqrt(2)


# (a) fit statistics over cells
S = F.mixture.meta["stats"].reshape(-1, 6)
fin = np.isfinite(S[:, 0])
rm1, re1, ae1, rm2, re2, ae2 = S[fin].T

# (b) random pairs
rng = np.random.default_rng(7)
pairs = []
for k in range(60):
    p = random_skin(rng)
    lam = float(rng.uniform(380, 780))
    A, hist, _ = mc_profile(lam, p, seed=5000 + k)
    if A < 1e-3:
        continue
    Ab, histb, _ = mc_profile(lam, p, seed=9000 + k)
    A1, l1 = mixture_log_profile(lam, p, 1)
    A2, l2 = mixture_log_profile(lam, p, 2)
    pairs.append((log_rmse(A, hist, l1), log_rmse(A, hist, l2), log_rmse(A, hist, l1, True), log_rmse(A, hist, l2, True),
                  noise_floor(A, hist, Ab, histb), noise_floor(A, hist, Ab, histb, True), A1 / A - 1, A2 / A - 1))
P = np.array(pairs)

# (c) albedo colour error on full spectra
M_XYZ = color.reflectance_to_xyz_matrix(WAVELENGTHS)
white = colour.XYZ_to_xy(M_XYZ.sum(1))
dE_mix, dE_raw, dE_tab = [], [], []
for k in range(12):
    p = random_skin(rng)
    ref = np.array([mc_profile(l, p, 30_000, seed=777 + 50 * k + i)[0] for i, l in enumerate(WAVELENGTHS)]) * (1 - SPECULAR)
    lab_ref = colour.XYZ_to_Lab(M_XYZ @ ref, white)
    dE_mix.append(colour.delta_E(colour.XYZ_to_Lab(M_XYZ @ F.mixture_albedo(WAVELENGTHS, p), white), lab_ref, method="CIE 2000"))
    dE_raw.append(colour.delta_E(colour.XYZ_to_Lab(M_XYZ @ F.mixture_albedo(WAVELENGTHS, p, preserve_albedo=False), white), lab_ref, method="CIE 2000"))
    dE_tab.append(colour.delta_E(colour.XYZ_to_Lab(M_XYZ @ F.albedo(WAVELENGTHS, p), white), lab_ref, method="CIE 2000"))
dE_mix, dE_raw, dE_tab = np.array(dE_mix), np.array(dE_raw), np.array(dE_tab)


def q(v):
    v = np.abs(v)
    return f"{np.median(v):.3f} | {np.percentile(v, 95):.3f} | {v.max():.3f}"


crit_rmse = np.percentile(P[:, 3], 95) < 0.05 if len(P) else False
crit_floor = np.percentile(P[:, 3] - P[:, 5], 95) < 0.05 if len(P) else False
crit_dE = dE_mix.max() < 0.5
text = f"""# V2 : mélange K=2 vs peau bicouche

## (a) Qualité du fit, sur toutes les cellules de la table ({int(fin.sum())} cellules ajustées, {int((~fin).sum())} trop sombres pour un profil)

| | médiane | p95 | max |
|---|---|---|---|
| RMSE log-profil K=1 (non pondérée) | {q(rm1)} |
| RMSE log-profil K=2 (non pondérée) | {q(rm2)} |
| RMSE log-profil K=1 (pondérée énergie) | {q(re1)} |
| RMSE log-profil K=2 (pondérée énergie) | {q(re2)} |
| erreur relative A, K=1 | {q(ae1)} |
| erreur relative A, K=2 | {q(ae2)} |

Les valeurs incluent le bruit Monte Carlo de la table (3·10⁴ photons par cellule).

## (b) Contrôle indépendant : {len(P)} couples (peau, λ ∈ 380–780 nm), mélange **interpolé** vs Monte Carlo direct

| | médiane | p95 | max |
|---|---|---|---|
| RMSE log-profil K=1 (non pondérée) | {q(P[:, 0])} |
| RMSE log-profil K=2 (non pondérée) | {q(P[:, 1])} |
| plancher de bruit MC (non pondéré) | {q(P[:, 4])} |
| RMSE log-profil K=1 (pondérée énergie) | {q(P[:, 2])} |
| RMSE log-profil K=2 (pondérée énergie) | {q(P[:, 3])} |
| plancher de bruit MC (pondéré énergie) | {q(P[:, 5])} |
| erreur relative A, K=1 | {q(P[:, 6])} |
| erreur relative A, K=2 | {q(P[:, 7])} |

Référence : Monte Carlo direct à 10⁵ photons. Le plancher de bruit est l'écart entre deux tirages indépendants de la référence.

## (c) Couleur de l'albédo (12 peaux, spectre 380–780 nm, D65)

| | médiane | p95 | max |
|---|---|---|---|
| ΔE00 albédo du mélange K=2 (avec correction d'albédo) vs MC direct | {q(dE_mix)} |
| ΔE00 albédo du mélange K=2 (sans correction) vs MC direct | {q(dE_raw)} |
| ΔE00 albédo de la table vs MC direct | {q(dE_tab)} |

Critères ADR-0005 : RMSE log-profil < 0,05 (pondérée énergie, p95, contrôle indépendant) → **{"OK" if crit_rmse else "NON ATTEINT"}** ;
écart au plancher de bruit < 0,05 (p95) → **{"OK" if crit_floor else "NON ATTEINT"}** ;
ΔE00 albédo < 0,5 → **{"OK" if crit_dE else "NON ATTEINT"}**.
"""
write_report("V2_mixture.md", text)

# (d) figure
cases = [("moyenne", SkinParams(melanin_fraction=0.08, blood_fraction=0.02, epidermis_thickness_cm=0.01)),
         ("foncée", SkinParams(melanin_fraction=0.35, eumelanin_ratio=0.9, epidermis_thickness_cm=0.01))]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, (name, p) in zip(axes, cases):
    for lam, c in zip([450.0, 550.0, 650.0], SERIES):
        A, hist, musp = mc_profile(lam, p, 200_000, seed=31)
        r_mm = 10 * CENTERS / musp
        ax.loglog(r_mm, np.maximum(hist / AREAS, 1e-9) * musp**2, color=c, lw=2, label=f"{lam:.0f} nm Monte Carlo")
        for K, ls in ((1, ":"), (2, "--")):
            AK, lK = mixture_log_profile(lam, p, K)
            ax.loglog(r_mm, AK * np.exp(lK) * musp**2, color=c, lw=1.3, ls=ls, label=f"{lam:.0f} nm K={K}")
    ax.set_xlim(1e-3, 20)
    ax.set_ylim(1e-4, 1e5)
    ax.set_xlabel("r (mm)")
    ax.set_title(f"Peau {name}", loc="left")
axes[0].set_ylabel("R(r) (cm⁻²)")
axes[0].legend(fontsize=7, ncol=3, loc="lower left")
save(fig, "phase1_v2_profiles.png")
