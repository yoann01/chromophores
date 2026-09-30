"""Refit ISSA (400-700 nm) and NIST (400-1000 nm) under SkinModel variants (corrected pheomelanin).

Question from the Opus concordance: do background x0.5 and dermal water survive a refit?
Output: docs/phase2/MODEL_VARIANTS.md
"""

import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(__file__)
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, ROOT)
from chromophores import inversion as inv  # noqa: E402
from chromophores.datasets import load_issa, load_nist  # noqa: E402
from chromophores.optics import SkinModel  # noqa: E402

BROAD = inv.Prior(np.zeros(inv.N), 16.0 * np.eye(inv.N), "broad")
STARTS = [np.array([m, 0.5, 0.0, 0.5, 0.0, 0.0]) for m in (-2.0, 0.0, 2.0)]
VARIANTS = {
    "fond ×0,5, eau (actuel)": dict(baseline_scale=0.5, dermis_water_fraction=0.65),
    "fond ×0,5, sans eau": dict(baseline_scale=0.5, dermis_water_fraction=0.0),
    "fond ×1, eau": dict(baseline_scale=1.0, dermis_water_fraction=0.65),
    "fond ×1, sans eau (Opus)": dict(baseline_scale=1.0, dermis_water_fraction=0.0),
}
_F = {}


def fit(args):
    y, lam, key, name = args
    k = (key, name)
    if k not in _F:
        _F[k] = inv.Forward(lam=lam, model=SkinModel(**VARIANTS[name]), specular=1.0)
    est = inv.map_estimate(y, np.eye(len(lam)), _F[k], BROAD, 0.005, STARTS)
    return est.x, est.residual_norm / np.sqrt(len(lam)) * 0.005


def main(per_cell=4):
    meta, lam_i, R_i = load_issa()
    sel = (lam_i >= 400) & (lam_i <= 700)
    ok = np.all(np.isfinite(R_i[:, sel]), 1)
    rng = np.random.default_rng(5)
    idx = []
    for _, g in meta[ok].groupby(["ethnicity", "location"]):
        idx += list(rng.choice(g.index.to_numpy(), min(len(g), per_cell), replace=False))
    lam_n, R_n = load_nist()
    lam = np.arange(400.0, 1000.0, 10.0)
    data = {"ISSA": [(R_i[i, sel], lam_i[sel], "issa") for i in idx],
            "NIST": [(np.interp(lam, lam_n, r), lam, "nist") for r in R_n]}
    lines = ["| Modèle | Données | n | RMS médian | RMS p95 | Mélanine | Sang | Épiderme (mm) | Échelle μs′ |", "|---|---|---|---|---|---|---|---|---|"]
    for name in VARIANTS:
        for dname, items in data.items():
            with Pool(os.cpu_count()) as pool:
                out = pool.map(fit, [it + (name,) for it in items], chunksize=4)
            X = inv.to_physical(np.array([o[0] for o in out]))
            rms = np.array([o[1] for o in out])
            m = dict(zip(inv.NAMES, np.median(X, 0)))
            row = (f"| {name} | {dname} | {len(items)} | {np.median(rms):.4f} | {np.percentile(rms, 95):.4f} | {m['melanin_fraction']:.3f} | "
                   f"{m['blood_fraction']:.4f} | {10 * m['epidermis_thickness_cm']:.3f} | {m['scattering_scale']:.2f} |")
            print(row, flush=True)
            lines.append(row)
    text = f"""# Variantes du modèle direct contre les mesures (phéomélanine corrigée)

Script : `scripts/phase2/model_variants_check.py`. Ajustement libre (a priori très large) de 298 spectres
ISSA (400–700 nm) et des 100 spectres NIST (400–1000 nm) ; médianes des paramètres ajustés.

{chr(10).join(lines)}
"""
    open(os.path.join(ROOT, "docs", "phase2", "MODEL_VARIANTS.md"), "w").write(text)


if __name__ == "__main__":
    main()
