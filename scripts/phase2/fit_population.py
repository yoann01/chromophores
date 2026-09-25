"""Phase 2, step 1: spectral fits of measured skin -> population prior data.

ISSA (400-700 nm, SCI) stratified by ethnic group x body location, and NIST
(400-1000 nm) for the scattering scale (the visible range alone does not
constrain it). A very broad prior (sigma = 4 in unbounded space) only
regularises degenerate directions. Output: data/derived/{issa,nist}_fits.npz.
"""

import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from chromophores import inversion as inv  # noqa: E402
from chromophores.datasets import load_issa, load_nist  # noqa: E402

OUT = os.path.join(HERE, "..", "..", "data", "derived")
BROAD = inv.Prior(np.zeros(inv.N), 16.0 * np.eye(inv.N), "broad")
STARTS = [np.array([m, 0.5, 0.0, 0.5, 0.0, 0.0]) for m in (-2.0, 0.0, 2.0)]
_F = {}


def fit(args):
    y, lam, key = args
    if key not in _F:
        _F[key] = inv.Forward(lam=lam, specular=1.0)
    est = inv.map_estimate(y, np.eye(len(lam)), _F[key], BROAD, 0.005, STARTS)
    return est.x, np.sqrt(np.diag(est.cov)), est.residual_norm / np.sqrt(len(lam)) * 0.005


def main(per_cell=40):
    meta, lam_i, R_i = load_issa()
    sel = (lam_i >= 400) & (lam_i <= 700)
    ok = np.all(np.isfinite(R_i[:, sel]), 1)
    rng = np.random.default_rng(2)
    idx = []
    for _, g in meta[ok].groupby(["ethnicity", "location"]):
        idx += list(rng.choice(g.index.to_numpy(), min(len(g), per_cell), replace=False))
    idx = np.array(sorted(idx))
    with Pool(os.cpu_count()) as pool:
        out = pool.map(fit, [(R_i[i, sel], lam_i[sel], "issa") for i in idx], chunksize=8)
    np.savez_compressed(os.path.join(OUT, "issa_fits.npz"), index=idx, x=np.array([o[0] for o in out]),
                        std=np.array([o[1] for o in out]), rms=np.array([o[2] for o in out]),
                        ethnicity=meta.ethnicity.to_numpy()[idx].astype(str), location=meta.location.to_numpy()[idx].astype(int),
                        names=np.array(inv.NAMES))
    print(f"ISSA: {len(idx)} spectra fitted", flush=True)
    lam_n, R_n = load_nist()
    lam = np.arange(400.0, 1000.0, 10.0)
    with Pool(os.cpu_count()) as pool:
        out = pool.map(fit, [(np.interp(lam, lam_n, r), lam, "nist") for r in R_n], chunksize=4)
    np.savez_compressed(os.path.join(OUT, "nist_fits.npz"), x=np.array([o[0] for o in out]), std=np.array([o[1] for o in out]),
                        rms=np.array([o[2] for o in out]), names=np.array(inv.NAMES))
    print("NIST: 100 spectra fitted", flush=True)


if __name__ == "__main__":
    main()
