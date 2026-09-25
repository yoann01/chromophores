"""Analytic cross-check: adding-doubling for the two-layer slab (iadpython, S. Prahl).

Plane-parallel, so it only gives the total diffuse reflectance A (no radial
profile). Same assumptions as the Monte Carlo tables: normal collimated
incidence, smooth Fresnel air/tissue boundary, index-matched layers,
Henyey-Greenstein phase function.

Convention (same as montecarlo / table): A is the diffuse reflectance per
photon *transmitted* into the tissue, i.e. the specular reflection is removed
and the result is divided by the entry Fresnel transmission (1 - R_sp). The
diffuse albedo seen by a cross-polarised camera is (1 - R_sp) * A.
"""

from __future__ import annotations

import copy

import iadpython as iad

from . import optics


def two_layer_reflectance(tau_e, a_d, d_e, g, n=optics.REFRACTIVE_INDEX, dermis_depth=200.0, quad_pts=16):
    """Total diffuse reflectance in the dimensionless units of chromophores.table."""
    mus = 1.0 / (1.0 - g)
    mua_e = tau_e / d_e
    epi = iad.Sample(a=mus / (mus + mua_e), b=(mus + mua_e) * d_e, g=g, d=d_e, n=n, n_above=1.0, n_below=n, quad_pts=quad_pts)
    epi.update_quadrature()
    der = copy.deepcopy(epi)
    der.a, der.b, der.d = mus / (mus + a_d), (mus + a_d) * dermis_depth, dermis_depth
    der.update_quadrature()
    R_e, T_e = iad.simple_layer_matrices(epi)
    R_d, T_d = iad.simple_layer_matrices(der)
    R02, R20, T02, T20 = iad.add_layers(epi, R_e, R_e, T_e, T_e, R_d, R_d, T_d, T_d)
    R01, R10, T01, T10 = iad.boundary_layer(epi, top=True)
    R03, _, T03, _ = iad.add_slide_above(epi, R01, R10, T01, T10, R02, R20, T02, T20)
    ur1, _, _, _ = epi.UX1_and_UXU(R03, T03)
    specular = ((n - 1.0) / (n + 1.0)) ** 2
    return (ur1 - specular) / (1.0 - specular)
