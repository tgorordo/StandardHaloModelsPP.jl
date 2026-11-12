# Introduction - Galactic-Frame Velocity Distributions

```@meta
CurrentModule = StandardHaloModelsPP
```

[`StandardHaloModelsPP.jl`](https://github.com/tgorordo/StandardHaloModelsPP.jl)
implements Dark Matter velocity distributions, and integrals of them, for the Standard Halo
Model (**SHM**) and an extension of the distribution, the **SHM++**. 
The SHM is an isotropic, truncated Maxwellian [Lewin:1995rx](@cite); 
by default its parameters follow the recommended conventions for reporting direct-detection results
 [Baxter:2021pqo](@cite). 
 The SHM++ adds an anisotropic component to better model local measurements, 
 accounting for the "Gaia Sausage" [10.1093/mnras/sty982, Necib_2019](@cite), as defined by:
> [PhysRevD.99.023012](@cite): [N. Wyn Evans, Ciaran A.J. O'Hare, Christopher McCabe, *Phys. Rev. D.* **99**, 023012 (2019)](@cite PhysRevD.99.023012)
> ([arXiv:1810.11468](https://arxiv.org/abs/1810.11468))

A [`HaloModel`](@ref) struct encodes the SHM++ parameterization, and [Presets](@ref) are provided 
for common choices. This page covers the galactic-frame distributions; their
[lab-frame counterparts](@ref "Galactic & Solar Ephemerides - the Earth-Frame Distributions")
and [integrals](@ref "Velocity Distribution Integrals") follow.
```@docs; canonical=false
HaloModel
HaloModel()
```

Every function takes a model as an optional first argument, defaulting to the [`SHM`](@ref)
(e.g. `fχ(v) == fχ(SHM, v)`), and every lab-frame function an optional time after it,
defaulting to [`T0`](@ref) (June 1st.).

!!! info "Physical Units"
    Quantities carry units using [Unitful.jl](https://juliaphysics.github.io/Unitful.jl/stable/),
    as re-exported by [NaturallyUnitful.jl](https://github.com/MasonProtter/NaturallyUnitful.jl):
    densities in e.g. `u"GeV/cm^3"`, speeds in e.g. `u"km/s"`. Any compatible units are
    accepted, so e.g. you may use [UnitfulAstro.jl](https://juliaastro.org/UnitfulAstro.jl/stable/).

## Presets

| Preset | ``\rho_0`` (GeV cm``^{-3}``) | ``v_0`` (km s``^{-1}``) | ``v_\text{esc}`` (km s``^{-1}``) | ``\beta`` | ``\eta`` | Use |
|:--|:--|:--|:--|:--|:--|:--|
| [`SHM`](@ref) | 0.3 | 238 | 544 | 0 | 0 | The SHM of [Baxter:2021pqo](@cite), and default model |
| [`SHMlegacy`](@ref) | 0.3 | 220 | 544 | 0 | 0 | A historical SHM, for comparison with older results |
| [`SHMnoesc`](@ref) | 0.3 | 238 | ``\infty`` | 0 | 0 | The SHM without escape-speed truncation |
| [`SHMerr`](@ref) | 0.3 | 238 ± 1.5 | 544 ± 33 | 0 | 0 | The SHM with measurement uncertainties, for sensitivity studies |
| [`SHMpp`](@ref) | 0.55 ± 0.17 | 233 ± 3 | 528 ± 25 | 0.9 ± 0.05 | 0.2 ± 0.1 | The SHM++ of [PhysRevD.99.023012](@cite) |

Lab-frame boosts (see more on [Ephemerides](@ref)) use the LSR speed ``v_\text{LSR} = v_0`` by default. 
Other models can be built with the [`HaloModel`](@ref) constructors:
e.g. `HaloModel(SHM; vₑ = 500u"km/s")`, or with a different LSR speed,
`HaloModel(SHMpp; VLSR = 238u"km/s")`.

## The Standard Halo Model

The SHM presets follow the conventions of [Baxter:2021pqo; Table 1](@cite), with the
historical ``v_0 = 220`` km/s variant of [PhysRevD.99.023012; Table I](@cite), an
untruncated variant, and a variant with measurement uncertainties:
```@docs; canonical=false
SHM
SHMlegacy
SHMnoesc
SHMerr
```

Its velocity distribution is a Maxwellian truncated at the escape speed
[PhysRevD.99.023012; Eqs. (1, 2)](@cite):
```@docs; canonical=false
fR(m::HaloModel, v::Unitful.Velocity)
```

## Anisotropic Refinement - the SHM++

The SHM++ adds two parameters for the so-called "Gaia Sausage", its anisotropy ``\beta`` and its
fraction ``\eta`` of the local DM [PhysRevD.99.023012; Table I](@cite):
```@docs; canonical=false
SHMpp
```

Its distribution mixes the round SHM component ``f_R`` with the anisotropic sausage component
``f_S`` [PhysRevD.99.023012; Eq. (3)](@cite):
```@docs; canonical=false
fχ(m::HaloModel, v::AbstractVector{<:Unitful.Velocity})
```
with the sausage [PhysRevD.99.023012; Eqs. (5-7)](@cite):
```@docs; canonical=false
fS(m::HaloModel, v::AbstractVector{<:Unitful.Velocity})
```

!!! tip "The Galactic-Frame Distribution"
    [`fχ`](@ref) is the full galactic-frame distribution, and usually what you want to use. 
    To evaluate it at many velocities, see [Performance Tips](@ref).

## Measurement Uncertainties

Parameters with measurement uncertainties use
[Measurements.jl](https://juliaphysics.github.io/Measurements.jl/stable/), as do the ephemeris
constants (see [Ephemerides](@ref)).

- Pointwise functions ([`fχ`](@ref), [`fχlab`](@ref), ...) follow the model: a model
  carrying `Measurement`s gives `Measurement` results, and an exact one plain numbers (with
  the central-value ephemeris).
- Integrals ([`fχspeed`](@ref), [`g`](@ref), [`velint`](@ref), [`HaloIntegralProblem`](@ref),
  ...) evaluate central values by default, since propagating through a quadrature is
  typically 10-130x slower; pass `uncertainty = true` to propagate.
- [`SHM`](@ref) is exact, following the reporting conventions of [Baxter:2021pqo](@cite);
  [`SHMerr`](@ref) attaches measurement uncertainties to its speeds for sensitivity studies
  (not for reported limits), and [`SHMpp`](@ref) carries the uncertainties of the SHM++.

The central-value model of any `HaloModel` is `value(m)` (from Measurements.jl), which
evaluates as fast as an exact model:
```@docs; canonical=false
Measurements.value(m::HaloModel)
```
