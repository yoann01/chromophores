"""Phase 2, step 2: Gaussian priors (unbounded space) from the population fits.

- ISSA fits give the chromophore distribution (population, per ethnic group,
  per face region).
- The scattering scale is weakly constrained by 400-700 nm: its marginal is
  taken from the NIST fits (400-1000 nm) and decorrelated from the rest.
Output: chromophores/data/priors_v1.npz
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from chromophores import inversion as inv  # noqa: E402

DERIVED = os.path.join(HERE, "..", "..", "data", "derived")
OUT = os.path.join(HERE, "..", "..", "chromophores", "data", "priors_v1.npz")
REGIONS = {"forehead": [6], "cheek": [2, 3], "chin": [4], "nose": [9], "face": [2, 3, 4, 6, 9], "body": [1, 5, 7, 8, 10, 11, 12]}
SI = inv.NAMES.index("scattering_scale")


def build(X, nist_x, label):
    p = inv.Prior.from_samples(X, label)
    mean, cov = p.mean.copy(), p.cov.copy()
    mean[SI] = nist_x[:, SI].mean()
    cov[SI, :] = cov[:, SI] = 0.0
    cov[SI, SI] = nist_x[:, SI].var() + 1e-3
    return inv.Prior(mean, cov, label, len(X))


def main():
    d = np.load(os.path.join(DERIVED, "issa_fits.npz"))
    nist = np.load(os.path.join(DERIVED, "nist_fits.npz"))["x"]
    X, eth, loc = d["x"], d["ethnicity"], d["location"]
    good = d["rms"] < 0.012  # drop failed fits
    priors = {"population": build(X[good], nist, "population")}
    for e in np.unique(eth):
        priors[f"ethnicity:{e}"] = build(X[good & (eth == e)], nist, f"ethnicity:{e}")
    for r, locs in REGIONS.items():
        priors[f"region:{r}"] = build(X[good & np.isin(loc, locs)], nist, f"region:{r}")
    np.savez_compressed(OUT, names=np.array(inv.NAMES), labels=np.array(list(priors)),
                        means=np.array([p.mean for p in priors.values()]), covs=np.array([p.cov for p in priors.values()]),
                        counts=np.array([p.n for p in priors.values()]))
    for k, p in priors.items():
        print(f"{k:22s} n={p.n:5d} median phys = {np.round(inv.to_physical(p.mean), 4)}")


if __name__ == "__main__":
    main()
