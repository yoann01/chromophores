"""Absorption spectra of the main skin chromophores.

All absorption coefficients are returned in cm^-1, wavelengths are in nm.

Sources
-------
- Hemoglobin (HbO2, Hb): S. Prahl, OMLC tabulated molar extinction
  coefficients (see data/hemoglobin_prahl.csv).
- Eumelanin / pheomelanin: power laws from Jacques (1998) and
  Donner & Jensen (2006), "A Spectral BSSRDF for Shading Human Skin".
- Water: Segelstein (1981) / Hale & Querry (1973), see data/water_segelstein.csv.
- Baseline (bloodless, melanin-free) tissue: Jacques (1998),
  mu_a = 7.84e8 * lambda^-3.255.
- beta-carotene and bilirubin: *approximate* Gaussian-mixture fits of the
  shape of published solution spectra, normalised to literature peak molar
  extinction. They must be replaced by tabulated data (e.g. PhotochemCAD)
  before being used for quantitative work. See docs/ETAT_DE_L_ART.md.
"""

from __future__ import annotations

from importlib import resources

import numpy as np

LN10 = np.log(10.0)
HB_MOLAR_MASS = 64500.0  # g/mol
HB_BLOOD_CONCENTRATION = 150.0  # g/L of whole blood

# Visible range used for colour computations.
WAVELENGTHS = np.arange(380.0, 781.0, 10.0)
# Full simulated range (aligned with Aliaga & Jarabo 2026).
WAVELENGTHS_FULL = np.arange(250.0, 1001.0, 2.0)


def _load_hemoglobin():
    path = resources.files("chromophores") / "data" / "hemoglobin_prahl.csv"
    with path.open() as f:
        table = np.loadtxt(f, delimiter=",", comments="#", skiprows=3)
    return table[:, 0], table[:, 1], table[:, 2]


_HB_LAMBDA, _EPS_HBO2, _EPS_HB = _load_hemoglobin()


def _load_water():
    path = resources.files("chromophores") / "data" / "water_segelstein.csv"
    with path.open() as f:
        table = np.loadtxt(f, delimiter=",", comments="#", skiprows=3)
    return table[:, 0], table[:, 1]


_WATER_LAMBDA, _WATER_MUA = _load_water()


def _interp(lam, x, y):
    lam = np.asarray(lam, dtype=float)
    if lam.min() < x.min() or lam.max() > x.max():
        raise ValueError(f"wavelength outside tabulated range [{x.min()}, {x.max()}] nm")
    return np.interp(lam, x, y)


def mua_oxyhemoglobin(lam):
    """Absorption of whole blood, fully oxygenated (cm^-1)."""
    eps = _interp(lam, _HB_LAMBDA, _EPS_HBO2)
    return LN10 * eps * HB_BLOOD_CONCENTRATION / HB_MOLAR_MASS


def mua_deoxyhemoglobin(lam):
    """Absorption of whole blood, fully deoxygenated (cm^-1)."""
    eps = _interp(lam, _HB_LAMBDA, _EPS_HB)
    return LN10 * eps * HB_BLOOD_CONCENTRATION / HB_MOLAR_MASS


def mua_blood(lam, so2):
    """Absorption of whole blood with oxygen saturation `so2` in [0, 1]."""
    return so2 * mua_oxyhemoglobin(lam) + (1.0 - so2) * mua_deoxyhemoglobin(lam)


def mua_water(lam):
    """Absorption of pure liquid water (cm^-1), Segelstein (1981)."""
    return _interp(lam, _WATER_LAMBDA, _WATER_MUA)


def mua_eumelanin(lam):
    """Absorption inside a eumelanin-filled melanosome (cm^-1)."""
    lam = np.asarray(lam, dtype=float)
    return 6.6e11 * lam ** -3.33


def mua_pheomelanin(lam):
    """Absorption inside a pheomelanin-filled melanosome (cm^-1)."""
    lam = np.asarray(lam, dtype=float)
    return 2.9e14 * lam ** -4.75


def mua_baseline(lam):
    """Background absorption of bloodless, melanin-free tissue (cm^-1)."""
    lam = np.asarray(lam, dtype=float)
    return 7.84e8 * lam ** -3.255


def _gaussian_mixture(lam, peaks, heights, widths):
    lam = np.asarray(lam, dtype=float)[..., None]
    g = np.asarray(heights) * np.exp(-0.5 * ((lam - np.asarray(peaks)) / np.asarray(widths)) ** 2)
    return g.sum(-1)


def eps_beta_carotene(lam):
    """APPROXIMATE molar extinction of beta-carotene (cm^-1 / M).

    Three-band vibronic structure (~425, 450, 478 nm) scaled to a peak of
    ~1.4e5 at 450 nm. Placeholder until tabulated data is imported.
    """
    bands = ([425.0, 451.0, 478.0, 400.0], [0.62, 1.0, 0.86, 0.25], [14.0, 14.0, 14.0, 30.0])
    return 1.39e5 * _gaussian_mixture(lam, *bands) / _gaussian_mixture(451.0, *bands)


def eps_bilirubin(lam):
    """APPROXIMATE molar extinction of bilirubin (cm^-1 / M), broad band at ~460 nm."""
    return 5.5e4 * _gaussian_mixture(lam, [458.0], [1.0], [32.0])


def mua_from_molar(eps, concentration_mol_per_l):
    return LN10 * eps * concentration_mol_per_l
