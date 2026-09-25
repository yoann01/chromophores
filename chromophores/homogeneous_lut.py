"""Scale-invariant LUT of a homogeneous semi-infinite random-walk medium.

For a homogeneous medium the profile shape depends only on (alpha, g);
sigma_t only rescales the radius (Chiang et al. 2016; Aliaga & Jarabo 2026):

    R(r | alpha, g, sigma_t) = sigma_t^2 * R(sigma_t r | alpha, g, 1)

This LUT stores, with sigma_t = 1, the total diffuse reflectance, the exit
radius moments and the radial profile on log bins. It is the building block
of the K-media mixture fit (chromophores.mixture). Same boundary as the
layered table: smooth Fresnel interface, n = 1.4, normal incidence,
specular excluded.
"""

from __future__ import annotations

import os
from importlib import resources

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from . import montecarlo
from .table import PROFILE_FLOOR, bin_areas, rho_centers

VERSION = 2
LUT_FILE = "homogeneous_lut_v2.npz"

# u = log(1 - alpha): alpha from 0.001 to 1 - 2e-5 (v2: extended to very absorbing media).
U = np.linspace(np.log(2e-5), np.log(0.999), 48)
G = np.array([0.0, 0.25, 0.5, 0.65, 0.75, 0.82, 0.88, 0.92])
RHO_EDGES = np.concatenate([[0.0], np.geomspace(1e-3, 1e3, 64)])  # units of 1/sigma_t
DEPTH = 1e5


def alpha_of_u(u):
    return 1.0 - np.exp(u)


def build(path, n_photons=50_000, verbose=True):
    shape = (len(U), len(G))
    A, r1, r2 = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    hist = np.zeros(shape + (len(RHO_EDGES) - 1,))
    for i, u in enumerate(U):
        a = alpha_of_u(u)
        for j, g in enumerate(G):
            A[i, j], r1[i, j], r2[i, j], hist[i, j] = montecarlo.run_binned(
                [1.0 - a], [a], [g], [DEPTH], RHO_EDGES, n_photons, seed=104729 * (1 + i * len(G) + j)
            )
        if verbose:
            print(f"homogeneous LUT: {i + 1}/{len(U)}", flush=True)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    np.savez_compressed(path, version=VERSION, n_photons=n_photons, u=U, g=G, rho_edges=RHO_EDGES, A=A, r1=r1, r2=r2, hist=hist)


def default_path():
    return resources.files("chromophores") / "data" / LUT_FILE


class HomogeneousLUT:
    def __init__(self, path=None):
        with np.load(path or default_path()) as f:
            t = {k: f[k] for k in f.files}
        self.meta = t
        self.u, self.g = t["u"], t["g"]
        self.edges = t["rho_edges"]
        self.log_centers = np.log(rho_centers(self.edges))
        prof = np.maximum(t["hist"] / bin_areas(self.edges), PROFILE_FLOOR)
        # Normalised profile density P = R(rho) / A (integrates to 1 over the plane).
        logP = np.log(prof) - np.log(t["A"])[..., None]
        self._logA = RegularGridInterpolator((self.u, self.g), np.log(t["A"]))
        self._logr1 = RegularGridInterpolator((self.u, self.g), np.log(t["r1"]))
        self._logP = RegularGridInterpolator((self.u, self.g, self.log_centers), logP, bounds_error=False, fill_value=None)

    def clip(self, u, g):
        return np.clip(u, self.u[0], self.u[-1]), np.clip(g, self.g[0], self.g[-1])

    def A(self, u, g):
        u, g = self.clip(u, g)
        return np.exp(self._logA(np.stack(np.broadcast_arrays(u, g), -1)))

    def r1(self, u, g):
        u, g = self.clip(u, g)
        return np.exp(self._logr1(np.stack(np.broadcast_arrays(u, g), -1)))

    def log_density(self, u, g, rho):
        """log P(rho | alpha, g) with sigma_t = 1 (per unit area, normalised)."""
        u, g = self.clip(u, g)
        lr = np.clip(np.log(np.maximum(rho, 1e-12)), self.log_centers[0], self.log_centers[-1])
        u, g, lr = np.broadcast_arrays(u, g, lr)
        return self._logP(np.stack([u, g, lr], -1))

    def invert_A(self, A, g):
        """u such that A(u, g) = A (A is monotonic in u)."""
        grid = self.A(self.u, np.full_like(self.u, g))  # decreasing with u
        return np.interp(-np.log(A), -np.log(grid), self.u)
