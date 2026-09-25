# chromophores

Recherche sur la simulation biophysique de la peau (chromophores → coefficients optiques
spectraux → subsurface scattering path-tracé) et sur l'inversion RGB/multispectral →
chromophores.

- [`docs/ETAT_DE_L_ART.md`](docs/ETAT_DE_L_ART.md) : synthèse bibliographique et jeux de données.
- [`docs/AXE_DE_RECHERCHE.md`](docs/AXE_DE_RECHERCHE.md) : axe proposé (identifiabilité et
  acquisition des chromophores), premiers résultats, programme.

## Code

| Module | Rôle |
|---|---|
| `chromophores/spectra.py` | Spectres d'absorption : HbO₂/Hb (Prahl), eu/phéomélanine, fond tissulaire, carotène et bilirubine (approchés) |
| `chromophores/optics.py` | Paramètres biophysiques → μa, μs, g par couche (épiderme / derme) |
| `chromophores/reflectance.py` | Modèle analytique rapide (filtre épidermique + dipôle) |
| `chromophores/montecarlo.py` | Marche aléatoire de référence (numba) : réflectance + profil radial R(r) |
| `chromophores/color.py` | Spectre → XYZ → sRGB linéaire (D65), caméras multispectrales idéalisées |

```bash
pip install -r requirements.txt
python scripts/plot_chromophores.py        # figure des spectres
python scripts/validate_forward_model.py   # analytique vs Monte Carlo
python scripts/identifiability.py          # Cramér-Rao + recherche de métamères RGB
python -m pytest -q tests
```

Les figures sont écrites dans `docs/figures/`.
