"""Figure 1: absorption spectra of the skin chromophores at typical concentrations."""

import numpy as np
from common import SERIES, plt, save

from chromophores import spectra
from chromophores.optics import SkinParams

p = SkinParams()
lam = np.arange(380.0, 1000.0, 2.0)
curves = [
    ("Eumélanine (épiderme, Cm=5 %)", p.melanin_fraction * spectra.mua_eumelanin(lam)),
    ("Phéomélanine (épiderme, Cm=5 %)", p.melanin_fraction * spectra.mua_pheomelanin(lam)),
    ("HbO2 (derme, 2 % sang)", p.blood_fraction * spectra.mua_oxyhemoglobin(lam)),
    ("Hb (derme, 2 % sang)", p.blood_fraction * spectra.mua_deoxyhemoglobin(lam)),
    ("β-carotène 0.4 µM (approx.)", spectra.mua_from_molar(spectra.eps_beta_carotene(lam), 0.4e-6)),
    ("Bilirubine 5 µM (approx.)", spectra.mua_from_molar(spectra.eps_bilirubin(lam), 5e-6)),
    ("Fond tissulaire", spectra.mua_baseline(lam)),
]
fig, ax = plt.subplots(figsize=(8, 4.5))
for (name, y), c in zip(curves, SERIES):
    ax.semilogy(lam, y, color=c, label=name)
ax.set_xlabel("Longueur d'onde (nm)")
ax.set_ylabel("µa (cm⁻¹)")
ax.set_ylim(1e-3, 1e3)
ax.set_title("Absorption des chromophores de la peau (concentrations typiques)", loc="left")
ax.legend(fontsize=8, ncol=2)
save(fig, "chromophores_absorption.png")
