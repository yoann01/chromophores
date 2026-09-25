# V1 : table multicouche vs Monte Carlo direct

- 200 couples (peau, λ) aléatoires, λ ∈ [250, 1000] nm ; référence : Monte Carlo direct, 10⁵ photons.
- Hors domaine de la table (coordonnées écrêtées) : 4 couples, exclus des statistiques.
- A < 10⁻³ (noir visuel, bruit MC) : 0 couples, exclus.
- Couples évalués : 196.

| Grandeur | médiane | p95 | max |
|---|---|---|---|
| A (table vs MC) | 0.46 % | 2.43 % | 5.80 % |
| ⟨ρ⟩ (table vs MC) | 0.36 % | 5.33 % | 16.72 % |
| ⟨ρ²⟩ (table vs MC) | 0.64 % | 6.01 % | 26.01 % |
| A (adding-doubling vs MC) | 0.22 % | 0.98 % | 4.68 % |

Critères ADR-0005 : A médiane < 1 %, p95 < 3 % → **OK** ;
moments p95 < 5 % → **NON ATTEINT**.

λ des couples hors domaine : 707, 870, 949, 1000 nm.
