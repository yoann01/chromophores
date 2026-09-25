"""K-media mixture fitted once per cell of the dimensionless table (ADR-0002).

For each cell of the layered table we fit a mixture of K homogeneous
random-walk media (K=1 and K=2), in units where mu_s' = 1:

    R_mix(rho) = sum_k w_k A_h(alpha_k, g_k) s_k^2 P_h(s_k rho | alpha_k, g_k)

with s_k = sigma_t,k / mu_s'. Loss: log normalised profile over the bins that carry energy (as in Aliaga &
Jarabo 2026), each bin weighted by its Monte Carlo noise (relative error
~ 1/sqrt(counts), plus a 2 % floor), plus a penalty on the total reflectance
such that a 1 % error on A costs as much as a 0.05 mean log error.

Reported metrics per cell: unweighted log-RMSE (all bins equal, tail-noise
dominated), energy-weighted log-RMSE (bins weighted by their share of the
reflected energy, i.e. what is visible at render time) and the error on A.

Parameters are stored in an interpolation-friendly space:
u_k = log(1 - alpha_k), g_k, log s_k, logit w.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from importlib import resources
from multiprocessing import Pool

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import least_squares

from .homogeneous_lut import HomogeneousLUT
from .table import bin_areas, rho_centers

VERSION = 2
FIT_FILE = "mixture_fit_v2.npz"
A_WEIGHT = 5.0  # 1 % error on A ~ 0.05 log-RMSE
MIN_BIN_ENERGY = 1e-5
G_MAX = 0.92


def _sig(x):
    return 1.0 / (1.0 + np.exp(-x))


def _logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


NOISE_FLOOR = 0.02


class CellTarget:
    def __init__(self, A, r1, hist, edges, n_photons=30_000):
        self.A, self.r1 = A, r1
        frac = hist / max(A, 1e-300)
        self.valid = (frac > MIN_BIN_ENERGY) & (hist > 0)
        self.rho = rho_centers(edges)[self.valid]
        self.logP = np.log(hist[self.valid] / bin_areas(edges)[self.valid] / A)
        counts = hist[self.valid] * n_photons
        w = 1.0 / np.sqrt(1.0 / np.maximum(counts, 1.0) + NOISE_FLOOR**2)
        self.w_noise = w / np.sqrt(np.mean(w**2)) if len(w) else w
        e = frac[self.valid]
        self.w_energy = e / e.sum() if len(e) else e


class MixtureModel:
    """tie_g: all components use the tissue anisotropy of the cell (fewer
    degenerate solutions, smoother across cells) instead of a free g_k."""

    def __init__(self, lut: HomogeneousLUT, tie_g=False):
        self.lut = lut
        self.u0, self.u1 = lut.u[0], lut.u[-1]
        self.tie_g = tie_g
        self.g_cell = 0.8

    # ----- parametrisation -------------------------------------------------
    def decode(self, x, K):
        comps = []
        n = 2 if self.tie_g else 3
        for k in range(K):
            u = self.u0 + (self.u1 - self.u0) * _sig(x[n * k])
            g = self.g_cell if self.tie_g else G_MAX * _sig(x[n * k + 1])
            s = np.exp(x[n * k + n - 1])
            comps.append((u, g, s))
        w = np.array([1.0]) if K == 1 else np.array([_sig(x[n * K]), 1 - _sig(x[n * K])])
        return comps, w

    def encode(self, comps, w):
        x = []
        for u, g, s in comps:
            x += [_logit((u - self.u0) / (self.u1 - self.u0))]
            if not self.tie_g:
                x += [_logit(np.clip(g, 1e-3, G_MAX - 1e-3) / G_MAX)]
            x += [np.log(s)]
        if len(comps) == 2:
            x.append(_logit(w[0]))
        return np.array(x)

    # ----- forward ----------------------------------------------------------
    def evaluate(self, comps, w, rho):
        """(A_mix, log normalised profile at rho)."""
        Ak = np.array([self.lut.A(u, g) for u, g, _ in comps]).ravel()
        dens = sum(w[k] * Ak[k] * s**2 * np.exp(self.lut.log_density(u, g, s * rho)) for k, (u, g, s) in enumerate(comps))
        A = float(np.dot(w, Ak))
        return A, np.log(np.maximum(dens, 1e-300)) - np.log(A)

    def residuals(self, x, K, target: CellTarget):
        comps, w = self.decode(x, K)
        A, logP = self.evaluate(comps, w, target.rho)
        n = len(target.rho)
        return np.concatenate([target.w_noise * (logP - target.logP), [np.sqrt(n) * A_WEIGHT * (A / target.A - 1)]])

    # ----- fitting ----------------------------------------------------------
    def moment_init(self, target: CellTarget, g):
        u = float(self.lut.invert_A(target.A, g))
        s = float(np.ravel(self.lut.r1(u, g))[0]) / target.r1
        return u, g, s

    def fit(self, target: CellTarget, g_cell, neighbour=None):
        self.g_cell = g_cell
        n = len(target.rho)
        u, g, s = self.moment_init(target, g_cell)
        if n < 5 or not np.isfinite(s):
            # Almost no light comes back (very dark epidermis): moment matching only.
            s = s if np.isfinite(s) else 1.0
            x1 = self.encode([(u, g, s)], [1.0])
            x2 = neighbour if neighbour is not None else self.encode([(u, g, s), (u, g, s)], [0.5])
            return {1: self.report(x1, 1, target), 2: self.report(x2, 2, target)}
        # K = 1
        best1 = least_squares(self.residuals, self.encode([(u, g, s)], [1.0]), args=(1, target), max_nfev=200)
        # K = 2, multi-start
        u_hi = min(u + 2.0, self.u1 - 0.1)
        starts = [
            self.encode([(u, g, s), (u_hi, g, 3 * s)], [0.8]),
            self.encode([(u - 0.5, g, 0.8 * s), (min(u + 1.0, self.u1 - 0.1), 0.5, 2 * s)], [0.6]),
            self.encode([(u, 0.3, s), (u, 0.85, 0.5 * s)], [0.5]),
        ]
        if neighbour is not None:
            starts.append(neighbour)
        sols = [least_squares(self.residuals, x0, args=(2, target), max_nfev=300) for x0 in starts]
        costs = np.array([sol.cost for sol in sols])
        pick = int(np.argmin(costs))
        if neighbour is not None and costs[-1] <= 1.05 * costs[pick]:
            pick = len(sols) - 1  # prefer continuity with the neighbouring cell
        return {1: self.report(best1.x, 1, target), 2: self.report(self.canonical(sols[pick].x), 2, target)}

    def refit_regularised(self, target: CellTarget, g_cell, x0, prior, mu):
        """Refit K=2 with a penalty sqrt(mu) * (x - prior) (smoothing across cells)."""
        self.g_cell = g_cell
        if len(target.rho) < 5:
            return prior

        def res(x):
            return np.concatenate([self.residuals(x, 2, target), np.sqrt(mu * len(target.rho)) * (x - prior)])

        return self.canonical(least_squares(res, x0, max_nfev=200).x)

    def canonical(self, x2):
        comps, w = self.decode(x2, 2)
        if comps[0][0] > comps[1][0]:  # component 1 = higher albedo
            return self.encode(comps[::-1], w[::-1])
        return x2

    def report(self, x, K, target: CellTarget):
        cK, wK = self.decode(x, K)
        n = len(target.rho)
        if n < 5:
            return {"x": x, "rmse": np.nan, "rmse_energy": np.nan, "A_err": np.nan, "n_bins": n}
        A, logP = self.evaluate(cK, wK, target.rho)
        return {
            "x": x,
            "rmse": float(np.sqrt(np.mean((logP - target.logP) ** 2))),
            "rmse_energy": float(np.sqrt(np.sum(target.w_energy * (logP - target.logP) ** 2))),
            "A_err": A / target.A - 1,
            "n_bins": n,
        }


# ---------------------------------------------------------------------------
# Table-wide fitting (multiprocessing)
# ---------------------------------------------------------------------------
_W = {}


def _init_worker(table_path, lut_path, tie_g):
    t = np.load(table_path)
    _W["t"] = {k: t[k] for k in t.files}
    _W["model"] = MixtureModel(HomogeneousLUT(lut_path), tie_g=tie_g)


def _target(idx):
    t = _W["t"]
    return CellTarget(t["A"][idx], t["r1"][idx], t["hist"][idx], t["rho_edges"], int(t["n_photons"]))


def _fit_slice(i):
    t, model = _W["t"], _W["model"]
    shape = t["A"].shape[1:]
    n2 = 5 if model.tie_g else 7
    x1 = np.zeros(shape + (2 if model.tie_g else 3,))
    x2 = np.zeros(shape + (n2,))
    for m in range(shape[1]):
        for q in range(shape[2]):
            prev = None
            for j in range(shape[0]):
                res = model.fit(_target((i, j, m, q)), float(t["g"][q]), prev)
                prev = res[2]["x"]
                x1[j, m, q], x2[j, m, q] = res[1]["x"], res[2]["x"]
    return i, x1, x2


def _refit_cells(args):
    cells, X2, mu = args
    t, model = _W["t"], _W["model"]
    out = []
    for idx in cells:
        nb = []
        for d in range(4):
            for step in (-1, 1):
                j = list(idx)
                j[d] += step
                if 0 <= j[d] < X2.shape[d]:
                    nb.append(X2[tuple(j)])
        prior = np.mean(nb, 0)
        out.append((idx, model.refit_regularised(_target(idx), float(t["g"][idx[3]]), X2[idx], prior, mu)))
    return out


def _stats(idx):
    t, model = _W["t"], _W["model"]
    model.g_cell = float(t["g"][idx[3]])
    tgt = _target(idx)
    r1 = model.report(_W["X1"][idx], 1, tgt)
    r2 = model.report(_W["X2"][idx], 2, tgt)
    return [r1["rmse"], r1["rmse_energy"], r1["A_err"], r2["rmse"], r2["rmse_energy"], r2["A_err"]]


def _stats_chunk(args):
    cells, X1, X2 = args
    _W["X1"], _W["X2"] = X1, X2
    return [(idx, _stats(idx)) for idx in cells]


def fit_table(table_path, lut_path, out_path, processes=4, tie_g=True, smoothing_passes=3, mu=0.02, verbose=True):
    """Fit every cell, then smooth the K=2 parameters with regularised red-black passes."""
    t = np.load(table_path)
    shape = t["A"].shape
    g_axis = t["g"]
    n1, n2 = (2, 5) if tie_g else (3, 7)
    X1, X2 = np.zeros(shape + (n1,)), np.zeros(shape + (n2,))
    cells = list(np.ndindex(shape))

    def chunks(seq, n):
        k = max(1, len(seq) // (4 * n))
        return [seq[i : i + k] for i in range(0, len(seq), k)]

    with Pool(processes, initializer=_init_worker, initargs=(table_path, lut_path, tie_g)) as pool:
        for i, x1, x2 in pool.imap_unordered(_fit_slice, range(shape[0])):
            X1[i], X2[i] = x1, x2
            if verbose:
                print(f"mixture fit: slice {i} done", flush=True)
        for it in range(smoothing_passes):
            for parity in (0, 1):
                subset = [c for c in cells if sum(c) % 2 == parity]
                for res in pool.imap_unordered(_refit_cells, [(ch, X2, mu) for ch in chunks(subset, processes)]):
                    for idx, x in res:
                        X2[idx] = x
            if verbose:
                print(f"smoothing pass {it + 1}/{smoothing_passes} done", flush=True)
        S = np.zeros(shape + (6,))
        for res in pool.imap_unordered(_stats_chunk, [(ch, X1, X2) for ch in chunks(cells, processes)]):
            for idx, st in res:
                S[idx] = st
    model = MixtureModel(HomogeneousLUT(lut_path), tie_g=tie_g)
    # Store decoded, interpolation-friendly parameters: u, g, log s (, logit w).
    P1 = np.zeros(shape + (3,))
    P2 = np.zeros(shape + (7,))
    for idx in cells:
        model.g_cell = float(g_axis[idx[3]])
        (c,), _ = model.decode(X1[idx], 1)
        P1[idx] = [c[0], c[1], np.log(c[2])]
        comps, w = model.decode(X2[idx], 2)
        P2[idx] = [comps[0][0], comps[0][1], np.log(comps[0][2]), comps[1][0], comps[1][1], np.log(comps[1][2]), _logit(w[0])]
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    np.savez_compressed(
        out_path,
        version=VERSION,
        tie_g=tie_g,
        smoothing_passes=smoothing_passes,
        mu=mu,
        tau_e=t["tau_e"],
        a_d=t["a_d"],
        d_e=t["d_e"],
        g=g_axis,
        k1=P1,
        k2=P2,
        stats=S,
    )


def default_path():
    return resources.files("chromophores") / "data" / FIT_FILE


@dataclass
class Mixture:
    """Random-walk media for a set of wavelengths. Arrays of shape (K, n_lambda)."""

    alpha: np.ndarray
    g: np.ndarray
    sigma_t: np.ndarray  # cm^-1
    weight: np.ndarray

    @property
    def sigma_s(self):
        return self.alpha * self.sigma_t

    @property
    def sigma_a(self):
        return (1 - self.alpha) * self.sigma_t


class MixtureTable:
    """Interpolates the fitted mixture parameters over the 4D grid."""

    def __init__(self, path=None):
        with np.load(path or default_path()) as f:
            t = {k: f[k] for k in f.files}
        self.meta = t
        axes = (np.log(t["tau_e"]), np.log(t["a_d"]), np.log(t["d_e"]), t["g"])
        self.axes = axes
        self._k1 = RegularGridInterpolator(axes, t["k1"])
        self._k2 = RegularGridInterpolator(axes, t["k2"])

    def __call__(self, coords, musp, K=2):
        """coords: (n, 4) log-space coordinates from LayeredTable.coords; musp in cm^-1."""
        if K == 1:
            v = self._k1(coords)
            return Mixture(1 - np.exp(v[None, :, 0]), v[None, :, 1], np.exp(v[None, :, 2]) * musp, np.ones((1, len(musp))))
        v = self._k2(coords)
        w1 = _sig(v[:, 6])
        return Mixture(
            np.stack([1 - np.exp(v[:, 0]), 1 - np.exp(v[:, 3])]),
            np.stack([v[:, 1], v[:, 4]]),
            np.stack([np.exp(v[:, 2]), np.exp(v[:, 5])]) * musp,
            np.stack([w1, 1 - w1]),
        )
