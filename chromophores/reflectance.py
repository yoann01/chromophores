"""Fast analytic diffuse reflectance of two-layer skin.

The epidermis is treated as a thin absorbing filter (Beer-Lambert, one
collimated pass in and one diffuse pass out) on top of a semi-infinite
dermis whose total diffuse reflectance comes from the dipole diffusion
model with a Fresnel-mismatched boundary (Farrell et al. 1992, integrated
over the surface). This is the kind of cheap model used in many
chromophore-estimation papers; montecarlo.py provides the reference it
should be checked against.
"""

from __future__ import annotations

import numpy as np

from . import optics
from .optics import SkinModel, SkinParams


def internal_reflection_A(n=optics.REFRACTIVE_INDEX):
    r_d = -1.440 / n**2 + 0.710 / n + 0.668 + 0.0636 * n
    return (1 + r_d) / (1 - r_d)


def semi_infinite_diffuse_reflectance(mua, musp, n=optics.REFRACTIVE_INDEX):
    """Total diffuse reflectance of a semi-infinite medium (dipole model)."""
    a = musp / (mua + musp)
    s = np.sqrt(3.0 * (1.0 - a))
    A = internal_reflection_A(n)
    return 0.5 * a * (1.0 + np.exp(-4.0 / 3.0 * A * s)) * np.exp(-s)


def diffuse_reflectance(lam, p: SkinParams, path_in=1.0, path_out=2.0, model: SkinModel = SkinModel.legacy()):
    """Diffuse (specular excluded) reflectance spectrum of skin (prototype model by default)."""
    epi, der = optics.layers(lam, p, model=model)
    t_epi = np.exp(-epi.mua * epi.thickness_cm * (path_in + path_out))
    return t_epi * semi_infinite_diffuse_reflectance(der.mua, der.musp)
