"""Phase 2 (ADR-0003): RGB (or spectrum) -> chromophores, MAP with explicit prior + Laplace uncertainty.

Parameter space: 6 biophysical parameters mapped to an unbounded vector x
(logit of the position in [low, high], on a log scale for positive
quantities). Priors are Gaussians in x (mean + full covariance), built from
measured spectra (ISSA / NIST, see scripts/phase2).

Measurement model (albedo map or spectrophotometer):
    y = M @ [R_sp * specular + (1 - R_sp) * A_table(lambda)]
with M a 3 x n_lambda camera/colour matrix (or identity for spectra) and
`specular` = 1 for SCI spectrophotometer data, 0 for a cross-polarised albedo.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from . import optics, table
from .forward import SPECULAR
from .optics import SkinModel, SkinParams

BOUNDS = {  # name: (low, high, log scale)
    "melanin_fraction": (0.001, 0.6, True),
    "eumelanin_ratio": (0.0, 1.0, False),
    "blood_fraction": (0.0005, 0.2, True),
    "oxygen_saturation": (0.3, 1.0, False),
    "epidermis_thickness_cm": (0.002, 0.035, True),
    "scattering_scale": (0.3, 3.0, True),
}
NAMES = list(BOUNDS)
N = len(NAMES)
LAMBDA = np.arange(380.0, 781.0, 10.0)


def _sig(v):
    return 1.0 / (1.0 + np.exp(-v))


def to_physical(x):
    """Unbounded vector(s) x[..., 6] -> physical values[..., 6]."""
    x = np.asarray(x, float)
    out = np.empty_like(x)
    for k, (lo, hi, lg) in enumerate(BOUNDS.values()):
        s = _sig(x[..., k])
        out[..., k] = np.exp(np.log(lo) + s * (np.log(hi) - np.log(lo))) if lg else lo + s * (hi - lo)
    return out


def to_unbounded(v):
    v = np.asarray(v, float)
    out = np.empty_like(v)
    for k, (lo, hi, lg) in enumerate(BOUNDS.values()):
        s = (np.log(v[..., k]) - np.log(lo)) / (np.log(hi) - np.log(lo)) if lg else (v[..., k] - lo) / (hi - lo)
        s = np.clip(s, 1e-6, 1 - 1e-6)
        out[..., k] = np.log(s / (1 - s))
    return out


def to_params(x, **extra):
    v = to_physical(x)
    return SkinParams(bilirubin_umol=0.0, **{n: float(val) for n, val in zip(NAMES, v)}, **extra)


@dataclass
class Prior:
    mean: np.ndarray  # (6,) in unbounded space
    cov: np.ndarray  # (6, 6)
    label: str = ""
    n: int = 0

    def __post_init__(self):
        self.L = np.linalg.cholesky(np.linalg.inv(self.cov))  # whitening: r = L^T (x - mean)

    def residual(self, x):
        return self.L.T @ (x - self.mean)

    @classmethod
    def from_samples(cls, X, label="", shrink=1e-3, inflate=1.0):
        """Gaussian fit in unbounded space; `inflate` widens it (e.g. to cover inter-dataset shift)."""
        mean = X.mean(0)
        cov = np.cov(X.T) * inflate**2 + shrink * np.eye(X.shape[1])
        return cls(mean, cov, label, len(X))

    def to_dict(self):
        return {"mean": self.mean, "cov": self.cov, "label": self.label, "n": self.n}


class Forward:
    """Albedo/spectrum forward model on a fixed wavelength grid, using the layered table."""

    def __init__(self, lam=LAMBDA, model: SkinModel = optics.DEFAULT_MODEL, specular=0.0, table_method="linear"):
        self.lam = np.asarray(lam, float)
        self.model = model
        self.specular = specular
        self.T = table.LayeredTable(method=table_method)

    def spectrum(self, x):
        p = to_params(x)
        tau, ad, de, g, _ = table.dimensionless(self.lam, p, self.model)
        A = np.exp(self.T._A(self.T.coords(tau, ad, de, g)))
        return self.specular * SPECULAR + (1 - SPECULAR) * A


@dataclass
class Estimate:
    x: np.ndarray  # MAP in unbounded space
    cov: np.ndarray  # Laplace covariance in unbounded space
    residual_norm: float  # data misfit (in noise units)
    y_model: np.ndarray

    @property
    def physical(self):
        return to_physical(self.x)

    def interval(self, k=1.0):
        """Per-parameter physical interval from x +/- k sigma (asymmetric, exact through the monotonic map)."""
        s = np.sqrt(np.diag(self.cov))
        return to_physical(self.x - k * s), to_physical(self.x + k * s)


def map_estimate(y, M, forward: Forward, prior: Prior, noise, x0s=None, max_nfev=200):
    """MAP of x given measurement y = M @ spectrum(x) (+ noise), Gaussian prior; Laplace covariance."""
    y = np.asarray(y, float)
    noise = np.broadcast_to(np.asarray(noise, float), y.shape)

    def res(x):
        return np.concatenate([(M @ forward.spectrum(x) - y) / noise, prior.residual(x)])

    starts = [prior.mean] if x0s is None else list(x0s)
    best = None
    for x0 in starts:
        sol = least_squares(res, x0, max_nfev=max_nfev, x_scale=1.0)
        if best is None or sol.cost < best.cost:
            best = sol
    J = best.jac
    cov = np.linalg.pinv(J.T @ J)
    r = best.fun[: len(y)]
    return Estimate(best.x, cov, float(np.sqrt(np.sum(r**2))), M @ forward.spectrum(best.x))


def load_priors(path=None):
    """{label: Prior} from chromophores/data/priors_v1.npz (scripts/phase2/build_priors.py)."""
    from importlib import resources

    with np.load(path or resources.files("chromophores") / "data" / "priors_v1.npz") as f:
        return {str(lab): Prior(m, c, str(lab), int(n)) for lab, m, c, n in zip(f["labels"], f["means"], f["covs"], f["counts"])}


def colour_exact(spectrum, y, M):
    """Smooth multiplicative correction so that M @ corrected = y exactly.

    corrected(lambda) = spectrum(lambda) * (1 + sum_k c_k B_k(lambda)), with B the
    (normalised, positive) rows of M: 3 broad smooth lobes. The biophysical spectral
    shape (haemoglobin bands, melanin slope) is kept; only a low-frequency tint is
    added. Returns (corrected spectrum, c).
    """
    B = np.abs(M) / np.abs(M).max(1, keepdims=True)
    A = M @ (spectrum[:, None] * B.T)
    c = np.linalg.solve(A, y - M @ spectrum)
    return spectrum * (1 + B.T @ c), c
