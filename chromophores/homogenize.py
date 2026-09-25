"""Per-wavelength homogenisation of two-layer skin for a hero-wavelength random walk.

PROTOTYPE (3D table, constant g, single medium). Superseded by the phase-1
modules `table` and `mixture` (ADR-0002); kept to reproduce docs/PRODUCTION.md.

Key observation: at a fixed wavelength, with index-matched layers sharing the
same scattering (mu_s, g) and a (quasi) semi-infinite dermis, the diffuse
response of the two-layer model depends only on three dimensionless numbers:

    tau_e = mu_a,epi * d_epi          epidermal absorption optical depth
    a_d   = mu_a,derm / mu_s'         dermal absorption / reduced scattering
    D     = d_epi * mu_s'             epidermal thickness in transport mfp

and lengths scale as 1 / mu_s'. A small 3D table, precomputed once by
Monte Carlo, therefore gives the multiple-scattering albedo A and mean exit
radius of *any* skin at *any* wavelength -- independently of the chromophore
spectra. At shade/hit time (hero wavelengths) we evaluate the chromophore
absorption analytically, look up (A, r_mean), and invert them into the
single-scattering albedo and extinction of a homogeneous random-walk medium
(moment matching, a generalisation of Chiang et al. 2016 albedo inversion).
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from . import montecarlo, optics
from .optics import SkinModel, SkinParams

LEGACY = SkinModel.legacy()

TABLE_FILE = "skin_homogenization_table.npz"

TAU_E = np.geomspace(1e-3, 20.0, 12)
A_D = np.geomspace(1e-3, 30.0, 14)
D_E = np.geomspace(0.05, 3.0, 6)
ALPHA = 1.0 - np.geomspace(3e-5, 0.999, 40)[::-1]  # homogeneous single-scattering albedo grid
DERMIS_DEPTH = 200.0  # in transport mfp: effectively semi-infinite


def _layered_point(tau_e, a_d, d_e, n_photons, seed):
    g = optics.ANISOTROPY
    mus = 1.0 / (1.0 - g)  # mu_s' = 1
    A, rbar, _ = montecarlo.run_layers(
        [tau_e / d_e, a_d], [mus, mus], [g, g], [d_e, DERMIS_DEPTH], n_photons, r_max=100.0, n_bins=10, seed=seed
    )
    return A, rbar


def _homogeneous_point(alpha, g_h, n_photons, seed):
    # sigma_t = 1: lengths in units of the mean free path.
    A, rbar, _ = montecarlo.run_layers(
        [1.0 - alpha], [alpha], [g_h], [1e4], n_photons, r_max=1e3, n_bins=10, seed=seed
    )
    return A, rbar


def build_table(n_photons=20_000, g_h=optics.ANISOTROPY, verbose=True):
    layered = np.zeros((len(TAU_E), len(A_D), len(D_E), 2))
    k = 0
    for i, t in enumerate(TAU_E):
        for j, a in enumerate(A_D):
            for m, d in enumerate(D_E):
                layered[i, j, m] = _layered_point(t, a, d, n_photons, seed=7 + k)
                k += 1
        if verbose:
            print(f"layered table: {i + 1}/{len(TAU_E)}")
    homog = np.array([_homogeneous_point(al, g_h, n_photons, seed=99991 + q) for q, al in enumerate(ALPHA)])
    return {"tau_e": TAU_E, "a_d": A_D, "d_e": D_E, "layered": layered, "alpha": ALPHA, "homogeneous": homog, "g_h": g_h}


def save_table(table, path):
    np.savez_compressed(path, **table)


def load_table(path=None):
    if path is None:
        path = resources.files("chromophores") / "data" / TABLE_FILE
    with np.load(path) as f:
        return {k: f[k] for k in f.files}


@dataclass
class HomogeneousMedium:
    """Random-walk medium for a set of wavelengths (cm^-1)."""

    sigma_s: np.ndarray
    sigma_a: np.ndarray
    g: float
    albedo_ms: np.ndarray  # target multiple-scattering albedo (diffuse reflectance)
    mean_radius_cm: np.ndarray

    @property
    def sigma_t(self):
        return self.sigma_s + self.sigma_a


class SkinHomogenizer:
    def __init__(self, table=None):
        t = table if table is not None else load_table()
        axes = (np.log(t["tau_e"]), np.log(t["a_d"]), np.log(t["d_e"]))
        self._logA = RegularGridInterpolator(axes, np.log(t["layered"][..., 0]))
        self._logr = RegularGridInterpolator(axes, np.log(t["layered"][..., 1]))
        self._axes = axes
        h = t["homogeneous"]
        self._alpha, self._hA, self._hr = t["alpha"], h[:, 0], h[:, 1]
        self.g_h = float(t["g_h"])

    def dimensionless(self, lam, p: SkinParams):
        epi, der = optics.layers(lam, p, model=LEGACY)
        musp = der.musp
        x = np.stack([epi.mua * epi.thickness_cm, der.mua / musp, epi.thickness_cm * musp], -1)
        return x, musp

    def lookup(self, x):
        """(A, mean radius in 1/mu_s' units) for dimensionless coordinates x[..., 3]."""
        lx = np.log(np.asarray(x))
        lx = np.stack([np.clip(lx[..., k], ax[0], ax[-1]) for k, ax in enumerate(self._axes)], -1)
        return np.exp(self._logA(lx)), np.exp(self._logr(lx))

    def invert(self, A, rbar_cm):
        """Moment matching: homogeneous (alpha, sigma_t) with same A and mean exit radius."""
        alpha = np.interp(A, self._hA, self._alpha)
        rbar_mfp = np.interp(alpha, self._alpha, self._hr)
        sigma_t = rbar_mfp / rbar_cm
        return alpha, sigma_t

    def medium(self, lam, p: SkinParams) -> HomogeneousMedium:
        """Shade-time evaluation: chromophores + hero wavelengths -> random-walk coefficients."""
        x, musp = self.dimensionless(np.atleast_1d(lam), p)
        A, rbar = self.lookup(x)
        rbar_cm = rbar / musp
        alpha, sigma_t = self.invert(A, rbar_cm)
        return HomogeneousMedium(alpha * sigma_t, (1 - alpha) * sigma_t, self.g_h, A, rbar_cm)
