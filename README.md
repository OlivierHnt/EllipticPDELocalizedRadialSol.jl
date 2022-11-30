# EllipticPDELocalizedRadialSol.jl

EllipticPDELocalizedRadialSol.jl is a Julia package to prove the existence of localized radial solutions of semilinear elliptic equations of the form

$$
\Delta U + \mathbf{N}(U) = 0, \qquad U = U(x) \in \mathbb{R}^q, \quad x \in \mathbb{R}^d,
$$

as described in the paper [Constructive proofs for localized radial solutions of semilinear elliptic systems on $\mathbb{R}^d$]().

The package exports one function `prove` which applies Proposition 2.2, Theorem 3.2 and Corollary 3.4 for the input data (cf. folder [data](https://github.com/OlivierHnt/EllipticPDELocalizedRadialSol.jl/tree/main/data)).

Each file in the folder [data](https://github.com/OlivierHnt/EllipticPDELocalizedRadialSol.jl/tree/main/data) contains the following variables (cf. Section 2 and 3 for more details):
- `equation`
- `x`
- `n_T_pad`
- `n_C_pad`
- `c`
- `Λ` (`\Lambda<tab>`)
- `Γ` (`\Gamma<tab>`)
- `l`
- `r_star`
- `L`
- `ℒ_x` (`\scrL<tab>_x`)
- `ℒ_y` (`\scrL<tab>_y`)
- `ν` (`\nu<tab>`)
- `ϱ` (`\varrho<tab>`)

## Executing the proofs

[Download this repository](https://github.com/OlivierHnt/EllipticPDELocalizedRadialSol.jl/archive/refs/heads/main.zip) and execute the following command in the downloaded directory:

```
EllipticPDELocalizedRadialSol_CompiledApp/bin/EllipticPDELocalizedRadialSol [FILENAME ...]
```

where each command-line argument given in `[FILENAME ...]` must be the path of a file contained in the folder [data](https://github.com/OlivierHnt/EllipticPDELocalizedRadialSol.jl/tree/main/data).
