# chromophores

Recherche sur la simulation biophysique de la peau (chromophores → coefficients optiques
spectraux → subsurface scattering path-tracé) et sur l'inversion RGB/multispectral →
chromophores.

- [`docs/ETAT_DE_L_ART.md`](docs/ETAT_DE_L_ART.md) : synthèse bibliographique et jeux de données.
- [`docs/TOPO_ALIAGA_JARABO_2026.md`](docs/TOPO_ALIAGA_JARABO_2026.md) : analyse du papier d'Aliaga &
  Jarabo (2026) et mise en parallèle avec nos résultats.
- [`docs/AXE_DE_RECHERCHE.md`](docs/AXE_DE_RECHERCHE.md) : axe proposé (identifiabilité et
  acquisition des chromophores), premiers résultats, programme.
- [`docs/PRODUCTION.md`](docs/PRODUCTION.md) : pipeline pour digital doubles (albédo VFace →
  chromophores par vertex → homogénéisation par λ héroïque en shade-before-hit).

- [`docs/adr/`](docs/adr/) : décisions d'architecture (ADR 0001–0005).
- [`docs/ROADMAP.md`](docs/ROADMAP.md) : phases, livrables et critères de passage.

## Code

| Module | Rôle |
|---|---|
| `chromophores/spectra.py` | Spectres d'absorption : HbO₂/Hb (Prahl), eu/phéomélanine, fond tissulaire, carotène et bilirubine (approchés) |
| `chromophores/optics.py` | Paramètres biophysiques → μa, μs, g par couche (épiderme / derme) |
| `chromophores/reflectance.py` | Modèle analytique rapide (filtre épidermique + dipôle) |
| `chromophores/montecarlo.py` | Marche aléatoire de référence (numba) : réflectance + profil radial R(r) |
| `chromophores/homogenize.py` | Table 3D sans dimension (Monte Carlo) → milieu homogène de random walk par longueur d'onde |
| `chromophores/color.py` | Spectre → XYZ → sRGB linéaire (D65), caméras multispectrales idéalisées |

```bash
pip install -r requirements.txt
python scripts/plot_chromophores.py        # figure des spectres
python scripts/validate_forward_model.py   # analytique vs Monte Carlo
python scripts/identifiability.py          # Cramér-Rao + recherche de métamères RGB
python scripts/production_metamerism.py    # calibration caméra + métamérisme d'illuminant
python scripts/build_homogenization_table.py  # ~5 min, régénère chromophores/data/*.npz
python scripts/validate_homogenization.py  # table vs MC, profils multicouche vs homogène
python -m pytest -q tests
```

Les figures sont écrites dans `docs/figures/`.
