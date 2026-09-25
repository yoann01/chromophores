import numpy as np

from chromophores import SkinParams, color, montecarlo, reflectance, spectra


def test_hemoglobin_isosbestic_near_800nm():
    lam = np.arange(780.0, 820.0, 2.0)
    diff = np.abs(spectra.mua_oxyhemoglobin(lam) - spectra.mua_deoxyhemoglobin(lam))
    assert 790 <= lam[np.argmin(diff)] <= 812


def test_white_maps_to_white():
    m = color.reflectance_to_linear_srgb_matrix()
    np.testing.assert_allclose(m.sum(1), 1.0, atol=2e-3)


def test_more_melanin_is_darker():
    lam = spectra.WAVELENGTHS
    light = reflectance.diffuse_reflectance(lam, SkinParams(melanin_fraction=0.02))
    dark = reflectance.diffuse_reflectance(lam, SkinParams(melanin_fraction=0.3))
    assert np.all(dark < light)


def test_monte_carlo_matches_diffusion_for_thick_bright_medium():
    # Without epidermal absorption and with a thick dermis, the dipole model
    # should be within ~10 % of the random walk at red wavelengths.
    p = SkinParams(melanin_fraction=0.0, carotene_umol=0.0, epidermis_thickness_cm=1e-6)
    lam = np.array([650.0])
    r_mc, _, _ = montecarlo.simulate(lam, p, n_photons=64_000, dermis_thickness_cm=2.0)
    r_an = reflectance.diffuse_reflectance(lam, p)
    assert abs(r_mc[0] - r_an[0]) / r_mc[0] < 0.1
