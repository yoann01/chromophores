"""V1 (ADR-0005): layered table interpolation vs direct Monte Carlo.

200 random (skin, wavelength) pairs over 250-1000 nm. For each pair the exact
dimensionless coordinates are simulated directly (1e5 photons) and compared
with the table interpolation (A, <rho>, <rho^2>). Adding-doubling gives an
independent check of A. Pairs with A < 1e-3 (visually black, dominated by MC
noise) and pairs outside the table domain are reported separately.
"""

import numpy as np
from common1 import random_skin, write_report

from chromophores import adding_doubling, table

N = 200
rng = np.random.default_rng(2026)
T = table.LayeredTable()

rows = []
for k in range(N):
    p = random_skin(rng)
    lam = float(rng.uniform(250, 1000))
    tau_e, a_d, d_e, g, _ = (float(v[0]) for v in table.dimensionless(np.array([lam]), p))
    raw = np.array([np.log(tau_e), np.log(a_d), np.log(d_e), g])
    inside = bool(np.all(np.abs(T.coords(tau_e, a_d, d_e, g) - raw) < 1e-9))
    A, r1, r2, _ = table.simulate_cell(tau_e, a_d, d_e, g, 100_000, seed=90001 + 17 * k)
    t = T(tau_e, a_d, d_e, g)
    ad = adding_doubling.two_layer_reflectance(tau_e, a_d, d_e, g) if inside else np.nan
    rows.append((lam, A, t.A.item(), r1, t.r1.item(), r2, t.r2.item(), ad, inside))

R = np.array(rows, dtype=float)
lam, A, tA, r1, tr1, r2, tr2, ad, inside = R.T
inside = inside.astype(bool)
ok = inside & (A >= 1e-3)


def stats(e):
    e = np.abs(e)
    return f"{100 * np.median(e):.2f} % | {100 * np.percentile(e, 95):.2f} % | {100 * e.max():.2f} %"


eA, e1, e2, eAD = tA / A - 1, tr1 / r1 - 1, tr2 / r2 - 1, ad / A - 1
crit_A = np.median(np.abs(eA[ok])) < 0.01 and np.percentile(np.abs(eA[ok]), 95) < 0.03
crit_m = np.percentile(np.abs(e1[ok]), 95) < 0.05 and np.percentile(np.abs(e2[ok]), 95) < 0.05
text = f"""# V1 : table multicouche vs Monte Carlo direct

- {N} couples (peau, λ) aléatoires, λ ∈ [250, 1000] nm ; référence : Monte Carlo direct, 10⁵ photons.
- Hors domaine de la table (coordonnées écrêtées) : {int((~inside).sum())} couples, exclus des statistiques.
- A < 10⁻³ (noir visuel, bruit MC) : {int((inside & (A < 1e-3)).sum())} couples, exclus.
- Couples évalués : {int(ok.sum())}.

| Grandeur | médiane | p95 | max |
|---|---|---|---|
| A (table vs MC) | {stats(eA[ok])} |
| ⟨ρ⟩ (table vs MC) | {stats(e1[ok])} |
| ⟨ρ²⟩ (table vs MC) | {stats(e2[ok])} |
| A (adding-doubling vs MC) | {stats(eAD[ok])} |

Critères ADR-0005 : A médiane < 1 %, p95 < 3 % → **{"OK" if crit_A else "NON ATTEINT"}** ;
moments p95 < 5 % → **{"OK" if crit_m else "NON ATTEINT"}**.
"""
if (~inside).any():
    text += f"\nλ des couples hors domaine : {', '.join(f'{v:.0f}' for v in np.sort(lam[~inside]))} nm.\n"
write_report("V1_table.md", text)
