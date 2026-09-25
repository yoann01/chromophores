# État de l'art : chromophores de la peau et subsurface scattering spectral

> Synthèse de départ. Le papier d'Aliaga & Jarabo (2026) est maintenant analysé en
> détail dans [`TOPO_ALIAGA_JARABO_2026.md`](TOPO_ALIAGA_JARABO_2026.md). Les références marquées ⚠ sont citées de mémoire et doivent
> être vérifiées (auteurs, année, venue).

## 1. Le problème

La peau est un milieu multicouche : stratum corneum, épiderme, derme papillaire, derme
réticulaire et hypoderme. Son apparence résulte de :

- **l'absorption par les chromophores** :
  - mélanines (eumélanine brun-noir, phéomélanine rouge-jaune) dans l'épiderme ;
  - hémoglobine oxygénée (HbO₂) et désoxygénée (Hb) dans le derme ;
  - chromophores mineurs : β-carotène, bilirubine, eau et lipides (surtout dans le proche IR) ;
  - absorption de fond du tissu ;
- **la diffusion** : de type Mie (fibres de collagène, organites) et de type Rayleigh
  (petites structures). Le μs′ décroît avec λ, et l'anisotropie est g ≈ 0,8–0,9.

En production, le SSS path-tracé (random walk) est paramétré par un **albédo diffus RGB**,
un **rayon de diffusion** RGB et une **anisotropie**. L'albédo de diffusion simple est
obtenu par *albedo inversion* (Christensen & Burley 2015 ; Chiang, Kutz & Burley 2016).
Le milieu est homogénéisé et le rayon est réglé à la main.

## 2. Lignées de travaux

### 2.1 Modèles biophysiques directs (graphics)
| Travail | Idée clé |
|---|---|
| Krishnaswamy & Baranoski 2004, *BioSpec* ⚠ | Modèle spectral stochastique multicouche, référence « physiologique » |
| Donner & Jensen 2006, *A Spectral BSSRDF for Shading Human Skin* | Deux couches, paramètres Cm, βm, Ch ; multipôle spectral. **Paramétrisation reprise ici** |
| Donner et al. 2008, *Layered heterogeneous reflectance model…* (SIGGRAPH Asia) | Cartes de chromophores spatialement variables, estimées à partir de photos multispectrales |
| Jimenez et al. 2010, *A Practical Appearance Model for Dynamic Facial Color* | Cartes mélanine/hémoglobine, dynamique du sang (émotions, pression) |
| Chen, Baranoski et al. 2015, *Hyperspectral Modeling of Skin Appearance* (TOG) | Modèle hyperspectral (HyLIoS), étend au proche IR/UV |
| Iglesias-Guitian et al. 2015, *Biophysically-based model of skin aging* (CGF) | Évolution des chromophores et de la structure avec l'âge |

### 2.2 Inversion : de l'image vers les chromophores
| Travail | Idée clé |
|---|---|
| Tsumura et al. 2003 (SIGGRAPH) | ICA sur densités optiques → cartes mélanine / hémoglobine depuis RGB |
| Alotaibi & Smith 2017 / BioFaceNet ⚠ | Modèle morphable biophysique, inversion par réseau |
| Gitlina et al. 2020, *Practical measurement and reconstruction of spectral skin reflectance* (EGSR) | Acquisition pratique (illumination multispectrale + polarisation) |
| Gevaux et al. 2021 ⚠ | Estimation temps réel des chromophores depuis de l'hyperspectral, par réseau |
| Aliaga et al. 2022, *Estimation of Spectral Biophysical Skin Properties from Captured RGB Albedo* ([arXiv 2201.10695](https://arxiv.org/pdf/2201.10695)) | Encodeur/décodeur RGB ↔ biophysique, entraîné sur des simulations |
| Aliaga et al. 2023, *A Hyperspectral Space of Skin Tones for Inverse Rendering of Biophysical Skin Properties* ([CGF](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.14887)) | Espace latent des tons de peau appris sur des spectres simulés |
| Biophysical skin model for heterogeneous volume rendering, CVM 2025 ([lien](https://www.sciopen.com/article/10.26599/CVM.2025.9450360)) | Deux couches volumétriques hétérogènes, coefficients appris depuis une texture d'albédo |
| **Aliaga & Jarabo 2026, *Spectral SSS from RGB via Biophysical Skin Inversion*** ([arXiv 2606.27604](https://arxiv.org/html/2606.27604v1)) | Mélange de 3 milieux non corrélés, décodeur neuronal chaîné RGB → (g, rayon, albédo) spectraux, random walk avec sélection aléatoire du milieu. Voir le [topo](TOPO_ALIAGA_JARABO_2026.md) |

### 2.3 Optique biomédicale (mesure quantitative)
- Jacques 2013, *Optical properties of biological tissues: a review* (Phys. Med. Biol.) : lois μs′(λ) et spectres de la mélanine et du fond.
- Wang, Jacques & Zheng 1995, *MCML* : Monte Carlo multicouche de référence (notre `montecarlo.py` en est une version minimale).
- Zonios et al. 2001 : mélanine, hémoglobine et diffusion *in vivo* par spectroscopie de réflectance diffuse (DRS).
- Cuccia et al. 2009, **SFDI** (*spatial frequency domain imaging*) : cartes de μa **et** μs′ à grand champ, obtenues en projetant des motifs sinusoïdaux. C'est la mesure la plus proche du « rayon de diffusion » d'un moteur de rendu.
- Tomographie des chromophores par imagerie hyperspectrale ([PubMed 41365253](https://pubmed.ncbi.nlm.nih.gov/41365253/)).
- Limites optiques mélanine/érythème ([bioRxiv 2025](https://www.biorxiv.org/content/10.64898/2025.12.22.696093v1.full)) : le signal de l'hémoglobine est atténué quand la mélanine augmente, ce qui rend la détection de l'érythème difficile sur peau foncée. Nos résultats d'identifiabilité (§4 d'`AXE_DE_RECHERCHE.md`) retrouvent ce comportement.

### 2.4 Données disponibles
| Jeu | Contenu |
|---|---|
| NIST, *Reference Data Set of Human Skin Reflectance* ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7339745/)) | 100 spectres traçables, 250–2500 nm, avec incertitudes |
| Hyper-Skin ([arXiv 2310.17911](https://arxiv.org/pdf/2310.17911)) | Visages hyperspectraux + RGB (VIS + NIR) : reconstruction spectrale depuis RGB |
| *Hyperspectral Imaging Database of Human Facial Skin* ([PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11823275/)) | 29 visages, 400–720 nm par pas de 10 nm, tons de peau variés |
| Tables OMLC (Prahl) | Hémoglobine (incluse dans le dépôt), eau, lipides |
| PhotochemCAD ⚠ | Spectres du β-carotène et de la bilirubine (à importer ; approximés pour l'instant) |

## 3. Limites identifiées dans l'approche RGB → spectral

1. **Le problème est mal posé.** 3 mesures pour 6 à 8 paramètres biophysiques. Le
   réseau apprend donc un *a priori* implicite, fixé par la distribution des données
   synthétiques d'entraînement.
2. **L'albédo seul ne contient pas le rayon de diffusion.** Deux peaux de même albédo
   RGB peuvent avoir des μs′ et des longueurs de diffusion très différentes (voir la
   figure `rgb_metamers.png`). Le rayon est alors « halluciné » par l'a priori.
3. **Chromophores mineurs quasi invisibles** en RGB : SO₂, carotène, bilirubine,
   ratio eu/phéomélanine.
4. **Peaux foncées** : la mélanine masque le derme, et l'inversion y dépend presque
   entièrement de l'a priori. C'est un enjeu d'équité de représentation.
5. **Données d'entraînement 100 % simulées** : le biais du modèle direct (diffusion
   analytique vs Monte Carlo, spectres approchés) se propage à l'inversion.
