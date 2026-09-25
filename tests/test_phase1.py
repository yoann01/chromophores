import numpy as np
import pytest

from chromophores import SkinParams, adding_doubling, table
from chromophores.optics import SkinModel, layers


def test_phase1_model_has_aliaga_anisotropy_and_water():
    lam = np.array([500.0, 970.0])
    epi, der = layers(lam, SkinParams())
    np.testing.assert_allclose(der.g, 0.62 + 0.00029 * lam)
    _, der_dry = layers(lam, SkinParams(), model=SkinModel(dermis_water_fraction=0.0))
    assert der.mua[1] - der_dry.mua[1] > 0.2  # water band at 970 nm
    assert abs(der.mua[0] - der_dry.mua[0]) < 1e-3


def test_adding_doubling_agrees_with_monte_carlo():
    cell = (0.5, 0.05, 0.3, 0.84)
    A_mc = table.simulate_cell(*cell, 40_000, seed=3)[0]
    assert abs(adding_doubling.two_layer_reflectance(*cell) / A_mc - 1) < 0.03


@pytest.fixture(scope="module")
def forward():
    from chromophores.forward import SkinForward

    try:
        return SkinForward()
    except FileNotFoundError:
        pytest.skip("phase-1 tables not built (scripts/phase1/build_tables.py, fit_mixture.py)")


def test_hero_wavelength_media_are_valid(forward):
    lam = np.array([450.0, 550.0, 650.0, 750.0])
    m = forward.media(lam, SkinParams())
    assert m.alpha.shape == (2, 4)
    np.testing.assert_allclose(m.weight.sum(0), 1.0)
    assert np.all((m.alpha > 0) & (m.alpha < 1) & (m.sigma_t > 0) & (m.g >= 0) & (m.g < 1))


def test_mixture_albedo_matches_table_albedo(forward):
    lam = np.arange(400.0, 781.0, 20.0)
    p = SkinParams()
    np.testing.assert_allclose(forward.mixture_albedo(lam, p), forward.albedo(lam, p), rtol=0.05)
