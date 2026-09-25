"""Dimensionless 4D table of the two-layer skin response (ADR-0002).

At a given wavelength, with index-matched layers sharing the same scattering
(mu_s, g) and a semi-infinite dermis, the diffuse response depends only on

    tau_e = mu_a,epi * d_epi      epidermal absorption optical depth
    a_d   = mu_a,derm / mu_s'     dermal absorption / reduced scattering
    D     = d_epi * mu_s'         epidermal thickness in transport mfp
    g     = anisotropy

with all lengths scaling as 1/mu_s'. Convention: reflectances are per photon
transmitted into the tissue (entry Fresnel excluded); the diffuse albedo seen
by a cross-polarised camera is (1 - R_sp) * A. Each cell stores, in units where
mu_s' = 1: the total diffuse reflectance A, the moments <rho>, <rho^2> of
the exit radius, and the radial profile on logarithmic bins.

The build is split in slices along tau_e so that it can be distributed
(render farm) and resumed; `merge_slices` assembles the final table.
"""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass
from importlib import resources

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from . import montecarlo

VERSION = 1
TABLE_FILE = "layered_table_v1.npz"

TAU_E = np.geomspace(1e-3, 30.0, 12)
A_D = np.geomspace(1e-3, 30.0, 14)
D_E = np.geomspace(0.02, 8.0, 7)
G = np.array([0.68, 0.76, 0.84, 0.92])
# Radial bins in units of 1/mu_s': [0, 1e-3] then 47 log bins up to 100.
RHO_EDGES = np.concatenate([[0.0], np.geomspace(1e-3, 100.0, 48)])
DERMIS_DEPTH = 200.0  # transport mfp: effectively semi-infinite
PROFILE_FLOOR = 1e-12
A_FLOOR = 1e-9


def rho_centers(edges=RHO_EDGES):
    c = np.sqrt(np.maximum(edges[:-1], edges[1] * 1e-1) * edges[1:])
    c[0] = 0.5 * edges[1]
    return c


def bin_areas(edges=RHO_EDGES):
    return np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)


def simulate_cell(tau_e, a_d, d_e, g, n_photons, seed):
    """Monte Carlo for one cell, lengths in units of 1/mu_s'."""
    mus = 1.0 / (1.0 - g)
    return montecarlo.run_binned(
        [tau_e / d_e, a_d], [mus, mus], [g, g], [d_e, DERMIS_DEPTH], RHO_EDGES, n_photons, seed
    )


def build_slice(i_tau, n_photons):
    shape = (len(A_D), len(D_E), len(G))
    A, r1, r2 = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    hist = np.zeros(shape + (len(RHO_EDGES) - 1,))
    for j, a in enumerate(A_D):
        for m, d in enumerate(D_E):
            for q, g in enumerate(G):
                seed = 1 + ((i_tau * len(A_D) + j) * len(D_E) + m) * len(G) + q
                A[j, m, q], r1[j, m, q], r2[j, m, q], hist[j, m, q] = simulate_cell(TAU_E[i_tau], a, d, g, n_photons, seed * 7919)
    return {"A": A, "r1": r1, "r2": r2, "hist": hist}


def build(out_dir, n_photons=30_000, slices=None, verbose=True):
    """Build (or resume) the per-slice files in out_dir."""
    os.makedirs(out_dir, exist_ok=True)
    for i in range(len(TAU_E)) if slices is None else slices:
        path = os.path.join(out_dir, f"slice_{i:03d}.npz")
        if os.path.exists(path):
            continue
        np.savez_compressed(path, n_photons=n_photons, **build_slice(i, n_photons))
        if verbose:
            print(f"slice {i + 1}/{len(TAU_E)} done", flush=True)


def merge_slices(out_dir, path):
    files = sorted(glob.glob(os.path.join(out_dir, "slice_*.npz")))
    if len(files) != len(TAU_E):
        raise RuntimeError(f"expected {len(TAU_E)} slices, found {len(files)}")
    parts = [np.load(f) for f in files]
    data = {k: np.stack([p[k] for p in parts]) for k in ("A", "r1", "r2", "hist")}
    np.savez_compressed(
        path,
        version=VERSION,
        n_photons=int(parts[0]["n_photons"]),
        tau_e=TAU_E,
        a_d=A_D,
        d_e=D_E,
        g=G,
        rho_edges=RHO_EDGES,
        dermis_depth=DERMIS_DEPTH,
        **data,
    )


def default_path():
    return resources.files("chromophores") / "data" / TABLE_FILE


@dataclass
class LayeredResponse:
    A: np.ndarray  # total diffuse reflectance
    r1: np.ndarray  # <rho> in 1/mu_s' units
    r2: np.ndarray  # <rho^2>
    profile: np.ndarray | None  # R(rho) per unit area, shape (..., n_bins)


class LayeredTable:
    """Multilinear interpolation (log axes, log values) of the 4D table."""

    def __init__(self, path=None, method="pchip"):
        with np.load(path or default_path()) as f:
            t = {k: f[k] for k in f.files}
        self.meta = t
        self.axes = (np.log(t["tau_e"]), np.log(t["a_d"]), np.log(t["d_e"]), t["g"])
        self.edges = t["rho_edges"]
        self.centers = rho_centers(self.edges)
        areas = bin_areas(self.edges)
        prof = np.maximum(t["hist"] / areas, PROFILE_FLOOR)
        # Scalars use monotone PCHIP (tensor product): 2-3x lower interpolation error than
        # multilinear on this grid (see docs/phase1). The profile stays multilinear.
        interp = lambda v, method=method: RegularGridInterpolator(self.axes, v, method=method)  # noqa: E731
        # Floors guard cells where (almost) no photon came back (very dark epidermis).
        self._A = interp(np.log(np.maximum(t["A"], A_FLOOR)))
        self._r1 = interp(np.log(np.maximum(t["r1"], RHO_EDGES[1])))
        self._r2 = interp(np.log(np.maximum(t["r2"], RHO_EDGES[1] ** 2)))
        self._prof = interp(np.log(prof), "linear")

    def coords(self, tau_e, a_d, d_e, g):
        x = [np.log(np.asarray(tau_e, float)), np.log(np.asarray(a_d, float)), np.log(np.asarray(d_e, float)), np.asarray(g, float)]
        x = np.broadcast_arrays(*x)
        return np.stack([np.clip(v, ax[0], ax[-1]) for v, ax in zip(x, self.axes)], -1)

    def __call__(self, tau_e, a_d, d_e, g, profile=False):
        x = self.coords(tau_e, a_d, d_e, g)
        return LayeredResponse(
            np.exp(self._A(x)),
            np.exp(self._r1(x)),
            np.exp(self._r2(x)),
            np.exp(self._prof(x)) if profile else None,
        )


def dimensionless(lam, p, model=None):
    """Skin parameters + wavelengths -> (tau_e, a_d, d_e, g, mu_s' [cm^-1])."""
    from . import optics

    model = model or optics.DEFAULT_MODEL
    epi, der = optics.layers(np.atleast_1d(lam), p, model=model)
    musp = der.musp
    return epi.mua * epi.thickness_cm, der.mua / musp, epi.thickness_cm * musp, der.g, musp
