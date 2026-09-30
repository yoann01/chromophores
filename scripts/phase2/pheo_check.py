"""Pheomelanin constant check (Donner & Jensen 2006, eq. 4: 2.9e14 lambda^-4.75 mm^-1 = 2.9e15 cm^-1).

spectra.py used 2.9e14 in cm^-1 (10x too low, now fixed). Refit a stratified ISSA subset (400-700 nm) and the
NIST spectra (400-1000 nm) with pheomelanin x1 (former) and x10 (paper), for background absorption
x0.5 (ADR-0007) and x1, and compare fit quality and fitted parameters.

Output: docs/phase2/PHEO_CHECK.md
"""

import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(__file__)
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, ROOT)
from chromophores import inversion as inv, spectra  # noqa: E402
from chromophores.datasets import load_issa, load_nist  # noqa: E402
from chromophores.optics import SkinModel  # noqa: E402

BROAD = inv.Prior(np.zeros(inv.N), 16.0 * np.eye(inv.N), "broad")
STARTS = [np.array([m, 0.5, 0.0, 0.5, 0.0, 0.0]) for m in (-2.0, 0.0, 2.0)]
VARIANTS = [(1.0, 0.5), (10.0, 0.5), (1.0, 1.0), (10.0, 1.0)]  # (pheomelanin factor, background scale)


def _ORIG(lam):  # the former (10x too low) constant, so that factors are relative to it
    return 2.9e14 * np.asarray(lam, float) ** -4.75


_F = {}


def fit(args):
    y, lam, key, pheo, base = args
    spectra.mua_pheomelanin = lambda l, f=pheo: f * _ORIG(l)  # noqa: E731
    k = (key, pheo, base)
    if k not in _F:
        _F[k] = inv.Forward(lam=lam, model=SkinModel(baseline_scale=base), specular=1.0)
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
    rows = []
    for pheo, base in VARIANTS:
        for name, items in data.items():
            with Pool(os.cpu_count()) as pool:
                out = pool.map(fit, [it + (pheo, base) for it in items], chunksize=4)
            X = inv.to_physical(np.array([o[0] for o in out]))
            rms = np.array([o[1] for o in out])
            med = dict(zip(inv.NAMES, np.median(X, 0)))
            rows.append((pheo, base, name, len(items), np.median(rms), np.percentile(rms, 95), med))
            print(f"pheo x{pheo:g} base x{base:g} {name}: RMS {np.median(rms):.4f} / {np.percentile(rms, 95):.4f}  "
                  f"Cm {med['melanin_fraction']:.3f} beta {med['eumelanin_ratio']:.2f} blood {med['blood_fraction']:.4f} "
                  f"d_e {10 * med['epidermis_thickness_cm']:.3f} mm  mus' x{med['scattering_scale']:.2f}", flush=True)
    lines = ["| Phéomélanine | Fond | Données | n | RMS médian | RMS p95 | Mélanine Cm | β (eu) | Sang | SO₂ | Épiderme (mm) | Échelle μs′ |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for pheo, base, name, n, m, p, med in rows:
        lines.append(f"| ×{pheo:g} | ×{base:g} | {name} | {n} | {m:.4f} | {p:.4f} | {med['melanin_fraction']:.3f} | {med['eumelanin_ratio']:.2f} | "
                     f"{med['blood_fraction']:.4f} | {med['oxygen_saturation']:.2f} | {10 * med['epidermis_thickness_cm']:.3f} | {med['scattering_scale']:.2f} |")
    text = f"""# Vérification de la constante de phéomélanine

Donner & Jensen 2006, éq. 3 et 4 : eumélanine 6,6·10¹⁰·λ⁻³·³³ mm⁻¹, phéomélanine 2,9·10¹⁴·λ⁻⁴·⁷⁵ mm⁻¹
(λ en nm), soit 6,6·10¹¹ et **2,9·10¹⁵** cm⁻¹. `spectra.py` avait 2,9·10¹⁴ cm⁻¹ pour la phéomélanine
(10 fois trop bas). Script : `scripts/phase2/pheo_check.py` (médianes des paramètres ajustés).

{chr(10).join(lines)}
"""  # interpretation appended by hand in PHEO_CHECK.md
    open(os.path.join(ROOT, "docs", "phase2", "PHEO_CHECK.md"), "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
