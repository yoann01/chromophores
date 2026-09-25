"""Reference Monte Carlo random walk in a two-layer slab (MCML-like).

Same kind of random walk as a path tracer uses for subsurface scattering:
exponential free flights, Henyey-Greenstein phase function, Fresnel
boundary at the air/skin interface, index-matched internal interface.
Returns the diffuse reflectance and its radial profile R(r), which is
what the "scattering radius" of rendering models tries to capture.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from . import optics
from .optics import SkinModel, SkinParams


@njit(cache=True)
def _fresnel(n1, n2, cos_i):
    sin_t2 = (n1 / n2) ** 2 * (1.0 - cos_i * cos_i)
    if sin_t2 >= 1.0:
        return 1.0
    cos_t = np.sqrt(1.0 - sin_t2)
    rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    rp = (n1 * cos_t - n2 * cos_i) / (n1 * cos_t + n2 * cos_i)
    return 0.5 * (rs * rs + rp * rp)


@njit(cache=True)
def _sample_hg(g, ux, uy, uz):
    if abs(g) < 1e-6:
        cost = 2.0 * np.random.random() - 1.0
    else:
        t = (1.0 - g * g) / (1.0 - g + 2.0 * g * np.random.random())
        cost = (1.0 + g * g - t * t) / (2.0 * g)
    cost = min(1.0, max(-1.0, cost))
    sint = np.sqrt(1.0 - cost * cost)
    phi = 2.0 * np.pi * np.random.random()
    cosp, sinp = np.cos(phi), np.sin(phi)
    if abs(uz) > 0.99999:
        return sint * cosp, sint * sinp, cost * np.sign(uz)
    tmp = np.sqrt(1.0 - uz * uz)
    nx = sint * (ux * uz * cosp - uy * sinp) / tmp + ux * cost
    ny = sint * (uy * uz * cosp + ux * sinp) / tmp + uy * cost
    nz = -sint * cosp * tmp + uz * cost
    return nx, ny, nz


@njit(parallel=True, cache=True)
def _run(n_photons, mua, mus, g, z_bounds, n_tissue, edges, seed):
    n_layers = mua.shape[0]
    n_bins = edges.shape[0] - 1
    hist = np.zeros(n_bins)
    n_chunks = 64
    per = n_photons // n_chunks
    partial = np.zeros((n_chunks, n_bins + 3))
    for c in prange(n_chunks):
        np.random.seed(seed + c)
        for _ in range(per):
            x = y = z = 0.0
            ux, uy, uz = 0.0, 0.0, 1.0
            w = 1.0  # specular reflection is excluded (diffuse albedo)
            layer = 0
            while True:
                tau = -np.log(np.random.random() + 1e-300)
                # Free flight, possibly crossing index-matched layer interfaces.
                alive = True
                while True:
                    mut = mua[layer] + mus[layer]
                    s = tau / mut
                    if uz > 0:
                        d_b = (z_bounds[layer + 1] - z) / uz
                    elif uz < 0:
                        d_b = (z_bounds[layer] - z) / uz
                    else:
                        d_b = 1e30
                    if s < d_b:
                        x += s * ux
                        y += s * uy
                        z += s * uz
                        break
                    x += d_b * ux
                    y += d_b * uy
                    z += d_b * uz
                    tau -= d_b * mut
                    if uz < 0 and layer == 0:
                        # Air/skin interface.
                        if np.random.random() < _fresnel(n_tissue, 1.0, -uz):
                            uz = -uz
                            z = 0.0
                        else:
                            r = np.sqrt(x * x + y * y)
                            partial[c, n_bins] += w
                            partial[c, n_bins + 1] += w * r
                            partial[c, n_bins + 2] += w * r * r
                            b = np.searchsorted(edges, r, side="right") - 1
                            if 0 <= b < n_bins:
                                partial[c, b] += w
                            alive = False
                            break
                    elif uz > 0 and layer == n_layers - 1:
                        alive = False  # transmitted out of the bottom
                        break
                    else:
                        layer += 1 if uz > 0 else -1
                if not alive:
                    break
                w *= mus[layer] / (mua[layer] + mus[layer])
                ux, uy, uz = _sample_hg(g[layer], ux, uy, uz)
                if w < 1e-4:
                    if np.random.random() < 0.1:
                        w *= 10.0
                    else:
                        break
    total = 0.0
    total_r = 0.0
    total_r2 = 0.0
    for c in range(n_chunks):
        total += partial[c, n_bins]
        total_r += partial[c, n_bins + 1]
        total_r2 += partial[c, n_bins + 2]
        for b in range(n_bins):
            hist[b] += partial[c, b]
    n = per * n_chunks
    norm = max(total, 1e-300)
    return total / n, total_r / norm, total_r2 / norm, hist / n


def run_binned(mua, mus, g, thickness, edges, n_photons=100_000, seed=1, n_tissue=optics.REFRACTIVE_INDEX):
    """Random walk in index-matched layers with arbitrary radial bin edges.

    Returns (A, <r>, <r^2>, hist) where hist[i] is the reflected weight per
    launched photon exiting with edges[i] <= r < edges[i+1].
    """
    z_bounds = np.concatenate([[0.0], np.cumsum(np.asarray(thickness, dtype=float))])
    return _run(n_photons, np.asarray(mua, float), np.asarray(mus, float), np.asarray(g, float), z_bounds, n_tissue, np.asarray(edges, float), seed)


def run_layers(mua, mus, g, thickness, n_photons=100_000, r_max=1.0, n_bins=100, seed=1, n_tissue=optics.REFRACTIVE_INDEX):
    """Random walk in a stack of index-matched layers (units are arbitrary but consistent).

    Returns (total diffuse reflectance, mean exit radius, radial histogram per photon).
    """
    A, r1, _, hist = run_binned(mua, mus, g, thickness, np.linspace(0.0, r_max, n_bins + 1), n_photons, seed, n_tissue)
    return A, r1, hist


def simulate(lam, p: SkinParams, n_photons=200_000, dermis_thickness_cm=0.2, r_max_cm=1.0, n_bins=100, seed=1, model: SkinModel = SkinModel.legacy()):
    """Diffuse reflectance per wavelength + radial profile R(r) [cm^-2].

    Returns (R, r_centers, profile) with profile of shape (len(lam), n_bins).
    """
    lam = np.atleast_1d(np.asarray(lam, dtype=float))
    epi, der = optics.layers(lam, p, dermis_thickness_cm, model)
    edges = np.linspace(0.0, r_max_cm, n_bins + 1)
    areas = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    R = np.empty(len(lam))
    prof = np.empty((len(lam), n_bins))
    for i in range(len(lam)):
        mua = np.array([epi.mua[i], der.mua[i]])
        mus = np.array([epi.mus[i], der.mus[i]])
        g = np.array([np.broadcast_to(epi.g, lam.shape)[i], np.broadcast_to(der.g, lam.shape)[i]])
        R[i], _, h = run_layers(mua, mus, g, [epi.thickness_cm, der.thickness_cm], n_photons, r_max_cm, n_bins, seed + 1000 * i)
        prof[i] = h / areas
    return R, 0.5 * (edges[1:] + edges[:-1]), prof
