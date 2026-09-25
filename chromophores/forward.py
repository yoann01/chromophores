"""Phase-1 forward model: chromophores + wavelengths -> albedo and random-walk media.

This is the reference implementation of what the engine evaluates at hit
time (ADR-0004): analytic absorption, dimensionless coordinates, lookup in
the layered table (albedo) and in the fitted mixture table (media).
"""

from __future__ import annotations

import numpy as np

from . import optics, table
from .homogeneous_lut import HomogeneousLUT
from .mixture import Mixture, MixtureTable
from .optics import SkinModel, SkinParams

SPECULAR = ((optics.REFRACTIVE_INDEX - 1) / (optics.REFRACTIVE_INDEX + 1)) ** 2


class SkinForward:
    def __init__(self, model: SkinModel = optics.DEFAULT_MODEL, table_path=None, mixture_path=None, lut_path=None):
        self.model = model
        self.table = table.LayeredTable(table_path)
        self.mixture = MixtureTable(mixture_path)
        self.lut = HomogeneousLUT(lut_path)

    def coords(self, lam, p: SkinParams):
        tau_e, a_d, d_e, g, musp = table.dimensionless(lam, p, self.model)
        raw = np.stack([np.log(tau_e), np.log(a_d), np.log(d_e), g], -1)
        clipped = self.table.coords(tau_e, a_d, d_e, g)
        return clipped, musp, np.any(np.abs(raw - clipped) > 1e-9, -1)

    def albedo(self, lam, p: SkinParams, include_entry_fresnel=True):
        """Diffuse albedo from the layered table (reference)."""
        x, _, _ = self.coords(lam, p)
        A = np.exp(self.table._A(x))
        return (1 - SPECULAR) * A if include_entry_fresnel else A

    def media(self, lam, p: SkinParams, K=2, preserve_albedo=True) -> Mixture:
        """Random-walk media at the given (hero) wavelengths.

        With preserve_albedo, the interpolated mixture is corrected so that its
        albedo equals the layered-table albedo exactly: all u_k = log(1 - alpha_k)
        are shifted by a common delta (A is monotonic in u), sigma_t and the
        weights are unchanged. This decouples colour accuracy from the
        (nonlinear) interpolation of the mixture parameters.
        """
        x, musp, _ = self.coords(lam, p)
        m = self.mixture(x, musp, K)
        if preserve_albedo:
            m = self._preserve_albedo(m, np.exp(self.table._A(x)))
        return m

    def _preserve_albedo(self, m: Mixture, A_target, iterations=40):
        u = np.log(np.maximum(1 - m.alpha, 1e-12))
        lo = np.full(A_target.shape, -4.0)
        hi = np.full(A_target.shape, 4.0)
        for _ in range(iterations):  # bisection on the common shift (A decreases with u)
            mid = 0.5 * (lo + hi)
            A = np.sum(m.weight * self.lut.A(u + mid, m.g), 0)
            too_bright = A > A_target
            lo = np.where(too_bright, mid, lo)
            hi = np.where(too_bright, hi, mid)
        u_new = np.clip(u + 0.5 * (lo + hi), self.lut.u[0], self.lut.u[-1])
        return Mixture(1 - np.exp(u_new), m.g, m.sigma_t, m.weight)

    def mixture_albedo(self, lam, p: SkinParams, K=2, include_entry_fresnel=True, preserve_albedo=True):
        """Albedo actually produced by the fitted mixture in a random walk."""
        m = self.media(lam, p, K, preserve_albedo)
        u = np.log(np.maximum(1 - m.alpha, 1e-12))
        A = np.sum(m.weight * self.lut.A(u, m.g), 0)
        return (1 - SPECULAR) * A if include_entry_fresnel else A
