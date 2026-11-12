# StandardHaloModelsPP
*A `julia` package encoding the Standard Halo Model ([SHM](https://arxiv.org/abs/2105.00599)) and anisotropic extension, the [SHM++](https://arxiv.org/abs/1810.11468).*

[![Build Status](https://github.com/tgorordo/StandardHaloModelsPP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tgorordo/StandardHaloModelsPP.jl/actions/workflows/CI.yml?query=branch%3Amain)

Possible Dark Matter (DM) direct detection signals depend on the local DM velocity distribution.
As a tool for experimental & phenomenological modeling, this julia package implements a
parameterized family of models of that velocity distribution based on the superposition of a spherical primary component
and an additional "sausage" anisotropic component. The package also provides a simple model of the Earth's motion
through the galaxy, and so can provide the DM velocity distribution in the lab as well as galactic frames of reference.

Utilities for commonly useful integrals of the distribution are also available, along with an
[Integrals.jl](https://github.com/SciML/Integrals.jl)-based problem interface for integrating
arbitrary kernels against it (with sweeps over parameters as threaded
SciML ensembles). The integrals are primarily implemented numerically, but in certain
circumstances (a zeroed anisotropic component) analytic forms are used.

Provides preset halo models:
- `SHM` (default): the SHM following the recommended conventions of Baxter et al. (2021).
- `SHMpp`: the SHM++ of Evans, O'Hare & McCabe (2019).
- `SHMlegacy`, `SHMnoesc` and `SHMerr`: respectively provide historical v₀ = 220 km/s SHM, an untruncated
  SHM, and the SHM *with* measurement uncertainties (against the Baxter recommendation) on its speeds for 
  sensitivity analyses.

See the [documentation](https://tgorordo.github.io/StandardHaloModelsPP.jl/) for more information.

## Installation
The package requires Julia 1.12 or later. Install it from this repository with the package manager:
```julia
using Pkg
Pkg.add(url = "https://github.com/tgorordo/StandardHaloModelsPP.jl")
```
or, in the REPL's package mode (press `]`), `add https://github.com/tgorordo/StandardHaloModelsPP.jl`.
Then `using StandardHaloModelsPP`. Physical units come from `NaturallyUnitful` (re-exporting
`Unitful`), which you will usually want to load too.

## Citations and License
If you use this package, you should likely cite at least (see [`CITATION.bib`](CITATION.bib)):

[1]: [N. Wyn Evans, Ciaran A.J. O'Hare, Christopher McCabe, *Phys. Rev. D.* **99**, 023012 (2019)](https://link.aps.org/doi/10.1103/PhysRevD.99.023012) ([arXiv:1810.11468](https://arxiv.org/abs/1810.11468)), for the SHM++ (`SHMpp`) parameterization;

[2]: [D. Baxter et al., *Eur. Phys. J. C* **81**, 907 (2021)](https://doi.org/10.1140/epjc/s10052-021-09655-y) ([arXiv:2105.00599](https://arxiv.org/abs/2105.00599)), for the standard conventions (`SHM` and derivatives) and the local standard of rest of the solar ephemeris;

[3]: [C. McCabe, *JCAP* **02**, 027 (2014)](https://doi.org/10.1088/1475-7516/2014/02/027) ([arXiv:1312.1355](https://arxiv.org/abs/1312.1355)), for the Earth ephemeris.

See also the [References in the documentation](https://tgorordo.github.io/StandardHaloModelsPP.jl/ref/) for some further attribution.
Otherwise, the source of this package is licensed under the very permissive [MIT License](LICENSE) - do with it what you will! 

## Development
A [Nix](https://nixos.org) flake provides a development environment: Julia, the [`just`](https://github.com/casey/just) taskrunner,
and a TeX Live/Asymptote toolchain used to generate documentation. Enter it with `nix develop`, or
automatically with [direnv](https://direnv.net) (`.envrc`). Then:

- `just test`, `just aqua`: run the test suite, or only the Aqua.jl checks.
- `just bench`: run the benchmarks (`test/runbench.jl`).
- `just format`: format `src` and `test` with JuliaFormatter.
- `just docs build`, `just docs serve`: build the documentation (PDF and HTML), and serve it locally.
