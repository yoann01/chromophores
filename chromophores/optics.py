"""Biophysical parameters -> per-layer optical coefficients (two-layer skin).

Layer 1 (epidermis): melanin (eu/pheo mixture), carotene, baseline.
Layer 2 (dermis):    blood (HbO2/Hb), carotene, bilirubin, baseline.

Parametrisation follows Donner & Jensen (2006) and Jacques (1998, 2013).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields, replace

import numpy as np

from . import spectra

# Reduced scattering of average skin, Jacques (2013): mu_s' = a (lambda/500)^-b.
SCATTER_A = 46.0  # cm^-1
SCATTER_B = 1.421
ANISOTROPY = 0.8
REFRACTIVE_INDEX = 1.4


@dataclass(frozen=True)
class SkinParams:
    melanin_fraction: float = 0.05  # volume fraction of melanosomes in epidermis
    eumelanin_ratio: float = 0.7  # eumelanin / (eu + pheo)
    blood_fraction: float = 0.02  # volume fraction of blood in dermis
    oxygen_saturation: float = 0.75  # SO2 of dermal blood
    carotene_umol: float = 0.4  # beta-carotene, umol/L of tissue
    bilirubin_umol: float = 0.0  # bilirubin, umol/L of tissue
    epidermis_thickness_cm: float = 0.01
    scattering_scale: float = 1.0  # multiplier on mu_s'

    def to_dict(self):
        return asdict(self)

    def with_(self, **kw):
        return replace(self, **kw)

    @classmethod
    def names(cls):
        return [f.name for f in fields(cls)]


@dataclass
class LayerOptics:
    mua: np.ndarray  # cm^-1
    mus: np.ndarray  # cm^-1
    g: float
    thickness_cm: float

    @property
    def musp(self):
        return self.mus * (1.0 - self.g)


def reduced_scattering(lam, scale=1.0):
    lam = np.asarray(lam, dtype=float)
    return scale * SCATTER_A * (lam / 500.0) ** -SCATTER_B


def epidermis_mua(lam, p: SkinParams):
    mel = p.eumelanin_ratio * spectra.mua_eumelanin(lam) + (1 - p.eumelanin_ratio) * spectra.mua_pheomelanin(lam)
    car = spectra.mua_from_molar(spectra.eps_beta_carotene(lam), p.carotene_umol * 1e-6)
    return p.melanin_fraction * mel + (1 - p.melanin_fraction) * spectra.mua_baseline(lam) + car


def dermis_mua(lam, p: SkinParams):
    blood = spectra.mua_blood(lam, p.oxygen_saturation)
    car = spectra.mua_from_molar(spectra.eps_beta_carotene(lam), p.carotene_umol * 1e-6)
    bil = spectra.mua_from_molar(spectra.eps_bilirubin(lam), p.bilirubin_umol * 1e-6)
    return p.blood_fraction * blood + (1 - p.blood_fraction) * spectra.mua_baseline(lam) + car + bil


def layers(lam, p: SkinParams, dermis_thickness_cm=0.2):
    musp = reduced_scattering(lam, p.scattering_scale)
    mus = musp / (1.0 - ANISOTROPY)
    epi = LayerOptics(epidermis_mua(lam, p), mus, ANISOTROPY, p.epidermis_thickness_cm)
    der = LayerOptics(dermis_mua(lam, p), mus, ANISOTROPY, dermis_thickness_cm)
    return epi, der
