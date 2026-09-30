"""ADR-0008 feasibility: can the epidermis be rendered as an absorbing interface filter over a
homogeneous dermis (plain random walk), instead of a K=2 mixture emulating the two-layer skin?

Reference: two-layer Monte Carlo (same scattering in both layers, ADR-0002). Candidates, all
compared on the radial profile shape (energy-weighted log RMSE, ADR-0006) and on the albedo:
  filter   epidermis = non-scattering absorbing layer (tau_e), dermis simulated exactly;
  homog.   one homogeneous medium whose albedo is matched to the reference, dermal scattering kept
           ("constant scattering", Masse & Walster 2023): melanin spread through the volume;
  K=2      current engine representation (table + mixture, ADR-0004 / ADR-0006).

Output: docs/phase3/FILTER_FEASIBILITY.md, docs/figures/phase3_filter_profiles.png
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(__file__)
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)
from common import plt, save  # noqa: E402

from chromophores import montecarlo, table  # noqa: E402
from chromophores.forward import SkinForward  # noqa: E402
from chromophores.mixture import MixtureModel  # noqa: E402
from chromophores.optics import SkinParams  # noqa: E402

N = 200_000
EDGES = table.RHO_EDGES
CENTERS = table.rho_centers(EDGES)
AREAS = table.bin_areas(EDGES)
F = SkinForward()
MM = MixtureModel(F.lut)

SKINS = {  # melanin volume fraction (epidermis 100 um, blood 2 %, SO2 75 %)
    "très claire": 0.01,
    "claire": 0.04,
    "mate": 0.12,
    "foncée": 0.30,
}
LAMS = [450.0, 550.0, 650.0]


def log_rmse(A, hist, logP, energy=True):
    valid = (hist / A > 1e-5) & (hist > 0) & np.isfinite(logP)
    ref = np.log(hist[valid] / AREAS[valid] / A)
    d2 = (logP[valid] - ref) ** 2
    if energy:
        w = hist[valid] / hist[valid].sum()
        return float(np.sqrt(np.sum(w * d2)))
    return float(np.sqrt(np.mean(d2)))


def logP_of(A, hist):
    with np.errstate(divide="ignore"):
        return np.where(hist > 0, np.log(np.maximum(hist, 1e-300) / AREAS / A), np.nan)


def mean_r(hist):
    return float(np.sum(hist * CENTERS) / np.sum(hist))


def case(lam, p, seed):
    tau_e, a_d, d_e, g, _ = (float(v[0]) for v in table.dimensionless(np.array([lam]), p))
    mus = 1.0 / (1.0 - g)
    depth = table.DERMIS_DEPTH
    out = {"tau_e": tau_e, "a_d": a_d, "d_e": d_e, "g": g}
    out["ref"] = montecarlo.run_binned([tau_e / d_e, a_d], [mus, mus], [g, g], [d_e, depth], EDGES, N, seed)
    out["ref2"] = montecarlo.run_binned([tau_e / d_e, a_d], [mus, mus], [g, g], [d_e, depth], EDGES, N, seed + 7)
    out["filter"] = montecarlo.run_binned([tau_e / d_e, a_d], [1e-9, mus], [g, g], [d_e, depth], EDGES, N, seed + 13)
    eps = 1e-4 * d_e  # zero-thickness interface filter: exp(-tau_e / |cos|) at each crossing
    out["filter0"] = montecarlo.run_binned([tau_e / eps, a_d], [1e-9, mus], [g, g], [eps, depth], EDGES, N, seed + 17)
    A_ref = out["ref"][0]

    def A_filter(t, n=60_000):
        return montecarlo.run_binned([t / eps, a_d], [1e-9, mus], [g, g], [eps, depth], EDGES, n, seed + 23)[0]

    # tau_eff such that the filter reproduces the reference albedo (secant on log A, A decreasing in tau)
    t0, A0 = tau_e, out["filter0"][0]
    t1 = max(t0 + np.log(A0 / A_ref) / 1.5, 1e-4)
    A1 = A_filter(t1)
    for _ in range(5):
        if abs(np.log(A1 / A_ref)) < 3e-3 or A1 == A0:
            break
        t2 = max(t1 - np.log(A1 / A_ref) * (t1 - t0) / (np.log(A1) - np.log(A0)), 1e-5)
        t0, A0, t1, A1 = t1, A1, t2, A_filter(t2)
    out["tau_eff"] = t1
    out["filterC"] = montecarlo.run_binned([t1 / eps, a_d], [1e-9, mus], [g, g], [eps, depth], EDGES, N, seed + 29)
    u = float(F.lut.invert_A(A_ref, g))
    alpha = 1.0 - np.exp(u)
    out["homog"] = montecarlo.run_binned([mus * (1 - alpha) / alpha], [mus], [g], [depth], EDGES, N, seed + 19)
    m = F.media(np.array([lam]), p, 2)
    musp = table.dimensionless(np.array([lam]), p)[4][0]
    comps = [(np.log(1 - m.alpha[k, 0]), m.g[k, 0], m.sigma_t[k, 0] / musp) for k in range(2)]
    out["k2_A"], out["k2_logP"] = MM.evaluate(comps, m.weight[:, 0], CENTERS)
    return out


def main():
    rows, examples = [], {}
    seed = 9001
    for name, mel in SKINS.items():
        p = SkinParams(melanin_fraction=mel)
        for lam in LAMS:
            c = case(lam, p, seed)
            seed += 101
            A, _, _, h = c["ref"]
            A2, _, _, h2 = c["ref2"]
            both = (h > 0) & (h2 > 0)
            floor = log_rmse(A, np.where(both, h, 0), np.where(both, logP_of(A2, h2), np.nan)) / np.sqrt(2)
            r = {"skin": name, "mel": mel, "lam": lam, "tau_e": c["tau_e"], "A": A, "floor": floor, "r_ref": mean_r(h)}
            for k in ("filter", "filter0", "filterC", "homog"):
                Ak, _, _, hk = c[k]
                r[f"{k}_rmse"] = log_rmse(A, h, logP_of(Ak, hk))
                r[f"{k}_dA"] = Ak / A - 1
                r[f"{k}_r"] = mean_r(hk) / r["r_ref"]
            r["tau_ratio"] = c["tau_eff"] / c["tau_e"]
            r["k2_rmse"] = log_rmse(A, h, c["k2_logP"])
            r["k2_dA"] = c["k2_A"] / A - 1
            rows.append(r)
            examples[(name, lam)] = c
            print(f"{name:12s} {lam:.0f} nm  tau_e={c['tau_e']:.3f}  A={A:.3f}  floor={floor:.3f}  "
                  f"filter {r['filter_rmse']:.3f} ({r['filter_dA']:+.1%})  filter0 {r['filter0_rmse']:.3f} ({r['filter0_dA']:+.1%})  filterC {r['filterC_rmse']:.3f} (tau_eff/tau {c['tau_eff'] / c['tau_e']:.2f}, {r['filterC_dA']:+.1%})  homog {r['homog_rmse']:.3f} (<r> x{r['homog_r']:.2f})  "
                  f"K2 {r['k2_rmse']:.3f}", flush=True)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for ax, key in zip(axes, [("claire", 550.0), ("foncée", 450.0), ("foncée", 650.0)]):
        c = examples[key]
        A, _, _, h = c["ref"]
        ok = h > 0
        ax.loglog(CENTERS[ok], (h / AREAS / A)[ok], "k-", lw=2, label="Référence bicouche (MC)")
        for k, lab, st in (("filter", "Filtre épais (100 µm) + derme", "C0--"), ("filter0", "Filtre d'interface + derme", "C1--"), ("filterC", "Filtre d'interface, τ_eff (albédo exact)", "C4-"), ("homog", "Milieu unique (Masse & Walster)", "C3:")):
            Ak, _, _, hk = c[k]
            okk = hk > 0
            ax.loglog(CENTERS[okk], (hk / AREAS / Ak)[okk], st, lw=1.8, label=lab)
        ax.loglog(CENTERS, np.exp(c["k2_logP"]), "C2-.", lw=1.4, label="Mélange K=2 (actuel)")
        ax.set_title(f"Peau {key[0]}, {key[1]:.0f} nm (τ_e = {c['tau_e']:.2f})", loc="left", fontsize=10)
        ax.set_xlabel("rayon (unités 1/μs′)")
        ax.set_ylim(1e-6, None)
    axes[0].set_ylabel("R(r) / A")
    axes[0].legend(fontsize=8)
    save(fig, "phase3_filter_profiles.png")

    def q(k):
        v = np.array([r[k] for r in rows])
        return np.median(v), np.percentile(v, 95) if v.ndim else v

    lines = ["| Peau | λ | τ_e | A réf. | Plancher de bruit | Filtre épais : RMSE forme | Filtre épais : ΔA | Filtre d'interface : RMSE forme | Filtre d'interface : ΔA | Filtre τ_eff : RMSE forme | τ_eff / τ_e | Milieu unique : RMSE forme | Milieu unique : ⟨r⟩ / réf. | K=2 : RMSE forme | K=2 : ΔA |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['skin']} | {r['lam']:.0f} | {r['tau_e']:.3f} | {r['A']:.3f} | {r['floor']:.3f} | {r['filter_rmse']:.3f} | "
                     f"{r['filter_dA']:+.1%} | {r['filter0_rmse']:.3f} | {r['filter0_dA']:+.1%} | {r['filterC_rmse']:.3f} | {r['tau_ratio']:.2f} | {r['homog_rmse']:.3f} | {r['homog_r']:.2f} | {r['k2_rmse']:.3f} | {r['k2_dA']:+.1%} |")
    summ = ["| Méthode | RMSE forme − plancher (médiane / max) | \\|ΔA\\| (médiane / max) |", "|---|---|---|"]
    fl = np.array([r["floor"] for r in rows])
    for k, lab in (("filter", "Filtre épais (100 µm, non diffusant) + derme"), ("filter0", "Filtre d'interface (épaisseur nulle) + derme"), ("filterC", "Filtre d'interface, τ_eff ajusté sur l'albédo"), ("homog", "Milieu unique (albédo imposé)"), ("k2", "Mélange K=2")):
        e = np.array([r[f"{k}_rmse"] for r in rows]) - fl
        dA = np.abs([r.get(f"{k}_dA", 0.0) for r in rows])
        dA_txt = "0 (imposé)" if k == "homog" else f"{np.median(dA):.1%} / {dA.max():.1%}"
        summ.append(f"| {lab} | {np.median(e):.3f} / {e.max():.3f} | {dA_txt} |")
    text = f"""# Faisabilité ADR-0008 : épiderme en filtre d'interface

Script : `scripts/phase3/filter_feasibility.py` ({N:,} photons par simulation). Référence : Monte Carlo
bicouche (épiderme 100 µm, sang 2 %, SO₂ 75 %, même diffusion dans les deux couches). Métrique de
forme : RMSE du log-profil pondérée par l'énergie (ADR-0006), à comparer au plancher de bruit.

{chr(10).join(summ)}

{chr(10).join(lines)}

![profils](../figures/phase3_filter_profiles.png)

## Constats

1. **Un épiderme non diffusant de 100 µm est exclu** : il crée un trou dans le profil en champ
   proche (plus aucune diffusion dans les 100 premiers µm).
2. **Le filtre d'interface (épaisseur nulle) a une bonne forme mais un albédo faux** (±15 à 25 %,
   −80 % sur peau foncée dans le bleu). Deux effets de signes opposés : les trajets obliques
   traversent le filtre sur une longueur d/|cos θ| (le filtre absorbe trop), et dans la vraie peau la
   mélanine agit aussi sur les photons qui diffusent plusieurs fois dans l'épiderme (le filtre
   absorbe trop peu dans le rouge).
3. **Avec une épaisseur optique effective τ_eff ajustée sur l'albédo**, le filtre bat le mélange
   K=2 en médiane (0,14 contre 0,25 au-dessus du bruit). Il est proche du plancher de bruit sur les
   peaux claires, en particulier dans le rouge, où la diffusion se voit (0,03–0,05 contre
   0,14–0,26). Il est moins bon sur peau foncée (bleu et rouge) : l'épiderme pigmenté produit un pic
   de diffusion simple en champ proche qu'aucun filtre ne peut reproduire. τ_eff / τ_e va de 0,65
   à 1,2. Sur peau très claire dans le bleu, aucun τ_eff ≥ 0 ne suffit : l'épiderme presque sans
   mélanine absorbe *moins* que le derme sanguin, et le filtre ne peut pas représenter une couche
   plus claire que ce qu'elle recouvre (écart d'albédo résiduel de 3 à 5 %).
4. **Le milieu unique (Masse & Walster)** est excellent sur peau très claire, mais comprime le
   profil sur peau mate ou foncée (⟨r⟩ × 0,36 à 0,67) : la peau foncée devient trop opaque.
5. **Aucune méthode ne passe le critère V2** (écart au plancher < 0,05). La RMSE du log-profil
   n'est pas une mesure visuelle : l'arbitrage doit se faire sur un test d'image (ADR-0008).
"""
    os.makedirs(os.path.join(ROOT, "docs", "phase3"), exist_ok=True)
    open(os.path.join(ROOT, "docs", "phase3", "FILTER_FEASIBILITY.md"), "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
