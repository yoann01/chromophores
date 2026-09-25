# Feuille de route

Les phases suivent l'ordre « le moins cher d'abord ». Chaque phase se termine par un
critère de passage tiré de l'[ADR-0005](adr/0005-validation-et-criteres.md).

| Phase | Contenu | ADR | Livrables | Critère de passage |
|---|---|---|---|---|
| **1. Modèle direct** | Modèle de peau aligné sur Aliaga & Jarabo (250–1000 nm, g(λ), eau, carotène et bilirubine tabulés) ; vérification analytique (adding-doubling, diffusion fréquentielle) ; table 4D sans dimension (Monte Carlo CPU/farm, bins log) ; fit K=2 par cellule | 0001, 0002 | Module direct + table + LUT homogène (α, g) + scripts de validation | V1, V2, V2b |
| **2. Inversion** | Modèle caméra / calibration peau ; LUT 3D RGB → chromophores (MAP + covariance) ; a priori population (bayésien empirique sur VFace) et par région ; gestion de l'AO | 0003 | LUT 3D + fichiers d'a priori + cartes de chromophores et d'incertitude sur 2 ou 3 visages VFace | V3 |
| **3. Intégration moteur** | Attributs par vertex ; évaluation au hit ; mélange K=2 avec MIS spectral ; re-remplissage de la table par le moteur ; filtrage en espace d'absorption | 0004 | Implémentation moteur + rendus de comparaison | V4 |
| **4. Mesures** | Protocole de mesure ponctuelle (spectre + profil radial) sur un acteur test ; calibration de l'a priori μs′ par acteur | 0003, 0005 | Protocole + données + a priori acteur | V5 |

## Prérequis à lever avant ou pendant la phase 1

- [x] Valider les ADR 0001 à 0005 (acceptés le 2026-09-25 ; questions ouvertes conservées comme risques).
- [ ] Récupérer des spectres tabulés du carotène et de la bilirubine (PhotochemCAD).
      Eau : fait (Segelstein).
- [ ] Choisir l'espace RGB de travail (sRGB linéaire ou ACEScg) et documenter la
      livraison VFace.
- [ ] Vérifier l'accès à la farm (ou à un poste multicœur) pour remplir la table
      (`scripts/phase1/build_tables.py layered --slices …` est distribuable).

## État

| Phase | État |
|---|---|
| 1 | **Prototype livré** (2026-09-25). V1-A ✅, V2-albédo ✅ ; V1-moments, V2-forme et V2b ❌. Voir [bilan](phase1/README.md) et [ADR-0006](adr/0006-ajustements-phase1.md) (proposé). |
| 2 | En attente de décision : démarrage possible (dépend seulement de l'albédo, validé). |
| 3 | Bloquée par V2-forme (table « farm » requise). |
| 4 | — |
