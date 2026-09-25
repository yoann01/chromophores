"""Spectrum -> CIE XYZ -> linear sRGB, and simple multispectral camera models."""

from __future__ import annotations

import colour
import numpy as np

from .spectra import WAVELENGTHS

XYZ_TO_LINEAR_SRGB = np.array(
    [
        [3.2404542, -1.5371385, -0.4985314],
        [-0.9692660, 1.8760108, 0.0415560],
        [0.0556434, -0.2040259, 1.0572252],
    ]
)


def _cmfs_and_d65(lam):
    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"].copy().align(colour.SpectralShape(360, 830, 1))
    d65 = colour.SDS_ILLUMINANTS["D65"].copy().align(colour.SpectralShape(300, 830, 1))
    xyz = np.stack([np.interp(lam, cmfs.wavelengths, cmfs.values[:, k]) for k in range(3)])
    e = np.interp(lam, d65.wavelengths, d65.values)
    return xyz, e


def reflectance_to_xyz_matrix(lam=WAVELENGTHS):
    """3 x len(lam) matrix M such that XYZ = M @ R (Y of the perfect white = 1)."""
    xyz, e = _cmfs_and_d65(lam)
    m = xyz * e
    return m / m[1].sum()


def reflectance_to_linear_srgb_matrix(lam=WAVELENGTHS):
    """3 x len(lam) matrix such that linear sRGB albedo = M @ R (white -> ~(1,1,1))."""
    return XYZ_TO_LINEAR_SRGB @ reflectance_to_xyz_matrix(lam)


def gaussian_bands_matrix(centers, fwhm, lam=WAVELENGTHS):
    """Idealised narrow-band multispectral camera (e.g. LED multiplexing)."""
    sigma = fwhm / 2.3548
    m = np.exp(-0.5 * ((lam[None, :] - np.asarray(centers)[:, None]) / sigma) ** 2)
    return m / m.sum(1, keepdims=True)


def encode_srgb(linear):
    linear = np.clip(linear, 0.0, 1.0)
    return np.where(linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055)
