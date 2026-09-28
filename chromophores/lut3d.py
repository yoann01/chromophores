"""Phase 2 (ADR-0003): 3D LUT  calibrated linear RGB -> chromophores (MAP) + uncertainty.

The grid is defined on sRGB-*encoded* values (as usual for 3D LUTs, better
resolution in the darks), 33^3 nodes by default. Only nodes close to the
gamut of real skin (measured ISSA / NIST spectra rendered as cross-polarised
albedo) are solved; the others are flagged invalid. Each valid node stores
the MAP in unbounded parameter space, the Laplace standard deviations, and
the reconstruction error.

Input convention: linear sRGB (D65) diffuse albedo, cross-polarised (no
specular). Other working spaces (ACEScg, camera RGB) go through a 3x3
conversion before the lookup.
"""

from __future__ import annotations

import os
from importlib import resources
from multiprocessing import Pool

import colour
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

from . import color
from . import inversion as inv

LUT_FILE = "rgb_lut_v1.npz"
M_SRGB = color.reflectance_to_linear_srgb_matrix(inv.LAMBDA)


def encode(lin):
    lin = np.clip(lin, 0.0, 1.0)
    return np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * np.power(lin, 1 / 2.4) - 0.055)


def decode(enc):
    enc = np.clip(enc, 0.0, 1.0)
    return np.where(enc <= 0.04045, enc / 12.92, np.power((enc + 0.055) / 1.055, 2.4))


def noise_model(rgb_lin, absolute=0.004, relative=0.02):
    """Per-channel measurement std: sensor/calibration floor + relative calibration error."""
    return absolute + relative * np.abs(rgb_lin)


def reference_skin_rgb():
    """Linear sRGB of measured skin (ISSA + NIST) as cross-polarised albedo (SCI minus specular)."""
    from .datasets import load_issa, load_nist
    from .forward import SPECULAR

    out = []
    meta, lam, R = load_issa()
    for r in R:
        ok = np.isfinite(r)
        if ok.sum() < 20:
            continue
        spec = np.interp(inv.LAMBDA, lam[ok], r[ok])  # edge values extended (tiny CMF weight)
        out.append(M_SRGB @ np.clip(spec - SPECULAR, 0, None))
    lam_n, R_n = load_nist()
    for r in R_n:
        out.append(M_SRGB @ np.clip(np.interp(inv.LAMBDA, lam_n, r) - SPECULAR, 0, None))
    return np.array(out)


def gamut_mask(grid_enc, reference_rgb, radius=0.06):
    """Nodes (encoded coordinates) within `radius` of a measured skin colour (encoded)."""
    tree = cKDTree(encode(reference_rgb))
    d, _ = tree.query(grid_enc)
    return d <= radius


_W = {}


def _solve(args):
    enc, prior_mean, prior_cov = args
    if "fwd" not in _W:
        _W["fwd"] = inv.Forward(specular=0.0)
    prior = inv.Prior(prior_mean, prior_cov)
    y = decode(enc)
    starts = [prior.mean, prior.mean + np.array([2.0, 0, 0, 0, 0, 0]), prior.mean + np.array([-2.0, 0, 1.0, 0, 0, 0])]
    est = inv.map_estimate(y, M_SRGB, _W["fwd"], prior, noise_model(y), starts)
    lab_in = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(y, apply_cctf_decoding=False))
    lab_out = colour.XYZ_to_Lab(colour.sRGB_to_XYZ(np.clip(est.y_model, 0, None), apply_cctf_decoding=False))
    return est.x, np.sqrt(np.diag(est.cov)), float(colour.delta_E(lab_in, lab_out, method="CIE 2000"))


def build(prior: inv.Prior, path, n=33, radius=0.06, processes=None, verbose=True):
    axis = np.linspace(0.0, 1.0, n)
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), -1).reshape(-1, 3)
    valid = gamut_mask(grid, reference_skin_rgb(), radius)
    nodes = np.where(valid)[0]
    if verbose:
        print(f"LUT {n}^3: {len(nodes)} nodes in the skin gamut", flush=True)
    X = np.full((len(grid), inv.N), np.nan)
    S = np.full((len(grid), inv.N), np.nan)
    dE = np.full(len(grid), np.nan)
    with Pool(processes or os.cpu_count()) as pool:
        for i, (x, s, de) in zip(nodes, pool.imap(_solve, [(grid[i], prior.mean, prior.cov) for i in nodes], chunksize=16)):
            X[i], S[i], dE[i] = x, s, de
    np.savez_compressed(path, n=n, axis=axis, valid=valid.reshape(n, n, n), x=X.reshape(n, n, n, -1), std=S.reshape(n, n, n, -1),
                        dE=dE.reshape(n, n, n), prior_mean=prior.mean, prior_cov=prior.cov, prior_label=prior.label,
                        names=np.array(inv.NAMES), radius=radius)


def default_path():
    return resources.files("chromophores") / "data" / LUT_FILE


class RGBToChromophores:
    """Applies the LUT to linear-sRGB pixels: returns physical MAP, 16-84 % bounds, flags."""

    def __init__(self, path=None):
        with np.load(path or default_path()) as f:
            t = {k: f[k] for k in f.files}
        self.meta = t
        axis = t["axis"]
        valid = t["valid"]
        x, s = t["x"].copy(), t["std"].copy()
        # Fill invalid nodes with the nearest valid one so that trilinear weights stay usable at the gamut border.
        idx_valid = np.argwhere(valid)
        tree = cKDTree(idx_valid)
        _, nn = tree.query(np.argwhere(~valid))
        src = tuple(idx_valid[nn].T)
        dst = tuple(np.argwhere(~valid).T)
        x[dst], s[dst] = x[src], s[src]
        dE = t["dE"].copy()
        dE[dst] = np.nan
        self._x = RegularGridInterpolator((axis,) * 3, x)
        self._s = RegularGridInterpolator((axis,) * 3, s)
        self._valid = RegularGridInterpolator((axis,) * 3, valid.astype(float))
        self._dE = RegularGridInterpolator((axis,) * 3, np.nan_to_num(dE, nan=99.0))
        self._cov = None
        if "cov_tri" in t:
            c = t["cov_tri"].astype(np.float64)
            c[dst] = c[src]
            self._cov = RegularGridInterpolator((axis,) * 3, c)

    def __call__(self, rgb_lin):
        shape = np.shape(rgb_lin)[:-1]
        e = encode(np.reshape(rgb_lin, (-1, 3)))
        x, s = self._x(e), self._s(e)
        lo, hi = inv.to_physical(x - s), inv.to_physical(x + s)
        out = {
            "map": inv.to_physical(x).reshape(shape + (inv.N,)),
            "p16": lo.reshape(shape + (inv.N,)),
            "p84": hi.reshape(shape + (inv.N,)),
            "in_gamut": (self._valid(e) > 0.999).reshape(shape),
            "dE_node": self._dE(e).reshape(shape),
            "x": x.reshape(shape + (inv.N,)),
        }
        if self._cov is not None:
            out["cov"] = tri_to_full(self._cov(e)).reshape(shape + (inv.N, inv.N))
        return out


TRI = np.triu_indices(inv.N)


def tri_to_full(tri):
    full = np.zeros(tri.shape[:-1] + (inv.N, inv.N))
    full[..., TRI[0], TRI[1]] = tri
    full[..., TRI[1], TRI[0]] = tri
    return full


def _node_cov(args):
    x, y, mean, cov = args
    if "fwd" not in _W:
        _W["fwd"] = inv.Forward(specular=0.0)
    prior = inv.Prior(mean, cov)
    noise = noise_model(y)

    def res(v):
        return np.concatenate([(M_SRGB @ _W["fwd"].spectrum(v) - y) / noise, prior.residual(v)])

    r0 = res(x)
    J = np.stack([(res(x + h) - r0) / 1e-4 for h in np.eye(inv.N) * 1e-4], 1)
    C = np.linalg.pinv(J.T @ J)
    return C[TRI]


def add_covariance(path, processes=None):
    """Add the full Laplace covariance (upper triangle) of every valid node to an existing LUT file."""
    with np.load(path) as f:
        t = {k: f[k] for k in f.files}
    valid = t["valid"]
    n = int(t["n"])
    axis = t["axis"]
    idx = np.argwhere(valid)
    args = [(t["x"][tuple(i)], decode(axis[i]), t["prior_mean"], t["prior_cov"]) for i in idx]
    with Pool(processes or os.cpu_count()) as pool:
        tri = pool.map(_node_cov, args, chunksize=32)
    cov_tri = np.full((n, n, n, len(TRI[0])), np.nan, np.float32)
    cov_tri[tuple(idx.T)] = np.array(tri, np.float32)
    t["cov_tri"] = cov_tri
    np.savez_compressed(path, **t)


def autocalibrate(rgb_lin, reference_rgb=None, tol=0.03, n_samples=100_000, seed=0):
    """Global exposure + white balance that best bring a texture onto the measured skin gamut.

    Low-dimensional on purpose (3 gains: exposure and R/B balance relative to G) so that it
    corrects a delivery/calibration error without absorbing the actual skin tone. The
    objective is robust (truncated distance), so non-skin texels (hair, eyes, clothes) weigh
    little. Returns (gains[3], fraction in gamut before, after).
    """
    from scipy.optimize import minimize

    ref = reference_skin_rgb() if reference_rgb is None else reference_rgb
    tree = cKDTree(encode(ref))
    flat = np.reshape(rgb_lin, (-1, 3))
    rng = np.random.default_rng(seed)
    pix = flat[rng.choice(len(flat), min(n_samples, len(flat)), replace=False)]

    def gains(v):
        e, r, b = np.exp(v)
        return np.array([e * r, e, e * b])

    def cost(v):
        d, _ = tree.query(encode(pix * gains(v)))
        return np.mean(np.minimum(d, 0.1))

    best = min(
        (minimize(cost, np.array([np.log(e0), 0.0, 0.0]), method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-5}) for e0 in (0.5, 0.75, 1.0, 1.3)),
        key=lambda r: r.fun,
    )
    g = gains(best.x)
    frac = lambda k: float(np.mean(tree.query(encode(pix * k))[0] <= tol))  # noqa: E731
    return g, frac(np.ones(3)), frac(g)
