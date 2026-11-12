# Velocity Distribution Integrals

```@meta
CurrentModule = StandardHaloModelsPP
```
[`StandardHaloModelsPP.jl`](https://github.com/tgorordo/StandardHaloModelsPP.jl) provides
a variety of utilities for integrals of the lab-frame distribution; its marginals, 
the mean speed and mean reciprocal speed above a threshold, and integrals of general kernels against it. 
They are computed in closed form where one exists, and otherwise with
[Integrals.jl](https://github.com/SciML/Integrals.jl). See also the [Performance Tips](@ref).

!!! note "Uncertainties"
    By default these integrals evaluate the central values of the model's parameters, even
    for a model carrying `Measurement`s such as [`SHMpp`](@ref). Pass `uncertainty = true`
    to propagate them instead. That is as much as 10-130x slower, since Measurements.jl
    tracks derivatives through every evaluation of the integrand.

!!! todo "TODO: Faster uncertainty propagation options"
    Dual numbers (e.g. [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)) would give the result and its gradient in one
    integration, with the escape-boundary terms still added separately. Monte Carlo over an
    `EnsembleProblem` of sampled models (e.g. [MonteCarloMeasurements.jl](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/)) 
    would give the full distribution, where first order is poor (e.g. ``v_\text{min}`` near [`vmaxlab`](@ref)).

## Marginalized Distributions

The lab-frame distribution marginalized over ``v_z \equiv v_\phi`` (the planar distribution)
and over solid angle (the speed distribution). For ``\eta = 0`` models the speed distribution
has a closed form [DelNobile:2021wmp; Eq. (I.55)](@cite), used by default:
```@docs; canonical=false
fχplanar(m::HaloModel, t::TimeLike, vr::Unitful.Velocity, vθ::Unitful.Velocity)
fχspeed(m::HaloModel, t::TimeLike, v::Unitful.Velocity)
```

These recreate plots akin to [PhysRevD.99.023012; Fig. 1](@cite), at [`T0`](@ref). The
paper's SHM is [`SHMlegacy`](@ref) here, and its SHM++ band spans ``\eta``, while the ribbon
here propagates all of `SHMpp`'s uncertainties:
```@example margin_dist_plots
using Plots, NaturallyUnitful, StandardHaloModelsPP
using Measurements: value, uncertainty

# Marginalized planar velocity distribution
vrrange = LinRange(-600, 600, 72)u"km/s"
vθrange = LinRange(-800, 200, 60)u"km/s"
# (the colorbar needs plain numbers, so scale and unit-strip by-hand)
z = [1e6 * ustrip(u"s^2/km^2", fχplanar(SHMpp, vr, vθ)) for vθ in vθrange, vr in vrrange]
p1 = contourf(vrrange, vθrange, z; color = cgrad(:RdYlBu, rev = true), lw = 0.2,
              xlabel = "v_r", ylabel = "v_θ", colorbar_title = "10⁶ f (s² km⁻²)")

# Speed distribution, with its uncertainty as a ribbon (Unitful values give axis units)
vrange = LinRange(0, 850, 60)u"km/s"
y = fχspeed.(SHMpp, vrange; uncertainty = true)
p2 = plot(vrange, value.(y); ribbon = uncertainty.(y), label = "SHM++",
          xlabel = "v", ylabel = "f(v)")
plot!(p2, vrange, fχspeed.(SHM, vrange); label = "SHM (analytic)", ls = :dash)
plot!(p2, vrange, fχspeed.(SHM, vrange; analytic = false); label = "SHM (numeric)", ls = :dot)
plot!(p2, vrange, fχspeed.(SHMlegacy, vrange); label = "SHMlegacy", ls = :dashdot)

plot(p1, p2; layout = 2, size = (950, 360), margin = 6Plots.mm)
```

## Mean Speed and Reciprocal

The mean reciprocal speed ``g`` (elsewhere often ``\eta(v_\text{min})``) and the mean speed
``h`` above a threshold ``v_\text{min}`` are the halo integrals of velocity-independent and
velocity-dependent scattering rates [DelNobile:2021wmp](@cite). For ``\eta = 0`` models they
have closed forms [DelNobile:2021wmp; Eqs. (I.60, I.63)](@cite), used by default:
```@docs; canonical=false
g(m::HaloModel, t::TimeLike, vmin::Unitful.Velocity)
h(m::HaloModel, t::TimeLike, vmin::Unitful.Velocity)
```

For the presets:
```@example mean_vel_plots
using Plots, NaturallyUnitful, StandardHaloModelsPP
using Measurements: value, uncertainty

vrange = LinRange(1, 850, 25)u"km/s"
y = g.(SHMpp, vrange; uncertainty = true)
p1 = plot(vrange, value.(y); ribbon = uncertainty.(y), label = "SHM++",
          title = "g(v_min)", xlabel = "v_min")
plot!(p1, vrange, g.(SHM, vrange); label = "SHM (analytic)", ls = :dash)
plot!(p1, vrange, g.(SHM, vrange; analytic = false); label = "SHM (numeric)", ls = :dot)
plot!(p1, vrange, g.(SHMlegacy, vrange); label = "SHMlegacy", ls = :dashdot)

y = h.(SHMpp, vrange; uncertainty = true)
p2 = plot(vrange, value.(y); ribbon = uncertainty.(y), label = "SHM++",
          title = "h(v_min)", xlabel = "v_min")
plot!(p2, vrange, h.(SHM, vrange); label = "SHM (analytic)", ls = :dash)
plot!(p2, vrange, h.(SHM, vrange; analytic = false); label = "SHM (numeric)", ls = :dot)
plot!(p2, vrange, h.(SHMlegacy, vrange); label = "SHMlegacy", ls = :dashdot)

plot(p1, p2; layout = 2, size = (950, 360), margin = 6Plots.mm)
```

## Arbitrary Convolutions

[`velint`](@ref) integrates a kernel against the full 3-D lab-frame distribution, and
[`speedint`](@ref) a scalar kernel against the speed distribution. A kernel may be any
callable taking `(v, p)` or just `(v)`. Both return the integral itself; for control over the
algorithm and error estimate, or for sweeps, use the [Integration Problems](@ref) they are
built on. For which to use, see [Performance Tips](@ref).
```@docs; canonical=false
velint(m::HaloModel, t::TimeLike, kk)
speedint(m::HaloModel, t::TimeLike, kk)
```

For example, ``g`` is `velint` with the kernel ``1/|\mathbf{v}|``. The two agree, here for
the SHM++, where `g` is itself numerical:
```@example velint_meanvel
using Plots, NaturallyUnitful, LinearAlgebra, StandardHaloModelsPP
using Measurements: value, uncertainty

gg(m, vmin) = velint(m, v -> 1 / norm(v); vmin)

vrange = LinRange(1, 850, 25)u"km/s"
y, yy = g.(SHMpp, vrange; uncertainty = true), gg.(SHMpp, vrange)
p = plot(vrange, value.(y); ribbon = uncertainty.(y), label = "g", title = "SHM++ g(v_min)",
         xlabel = "v_min")
scatter!(p, vrange, yy; label = "via velint", markersize = 3)

p # hide
```

## Integration Problems

`velint` and `speedint` wrap a [`HaloIntegralProblem`](@ref): an Integrals.jl
[`IntegralProblem`](@extref Integrals :std:doc:`basics/IntegralProblem`) for the
integral, with units handled for you. Use the problem directly to choose an
[algorithm](@extref Integrals :std:doc:`solvers/IntegralSolvers`), inspect the error
estimate or return code, or `remake` it with a parameter changed. `solve`, `remake` and the
SciML ensemble types (`EnsembleProblem`, `EnsembleSerial`, `EnsembleThreads`,
`EnsembleDistributed`) are re-exported.

```@docs; canonical=false
HaloIntegralProblem
HaloIntegralSolution
```

### Sweeps

A sweep over time, ``v_\text{min}``, kernel parameters or models is a SciML
`EnsembleProblem` whose `prob_func` remakes the problem for each point. Rebuilding costs
about 3% of a solve, and the result is an `EnsembleSolution` of
[`HaloIntegralSolution`](@ref)s. With `EnsembleThreads()` the points are solved in parallel
(start Julia with e.g. `-t auto`):
```julia
using Dates
ts = [DateTime(2020, 1, 1) + Day(d) for d in 0:364]            # annual modulation
ens = EnsembleProblem(prob; prob_func = (p, ctx) -> remake(p; t = ts[ctx.sim_id]))
sols = solve(ens, nothing, EnsembleThreads(); trajectories = length(ts))
[s.u for s in sols]
```
Passing `nothing` as the algorithm uses the problem's default.

### Building Your Own Integrals

To integrate against the halo in your own `IntegralProblem` (e.g. jointly with a recoil
energy), use the distributions frozen at a time, which fix the lab frame and normalizations
once and evaluate cheaply: `fχspeed(m, t)` for the speed distribution (analytic where
possible, and quiet inside an integrand) and `fχlab(m, t)` for the full velocity
distribution. [`vmaxlab`](@ref) gives the largest lab-frame speed with support, and
`vE(m, t)` / `VE(m, t)` the model's Earth velocity and speed.
```@docs; canonical=false
vmaxlab
```

!!! todo "TODO: Vector space methods for anisotropic detectors"
    Provide basis projections (spherical harmonics × radial functions) of the distributions
    for the vector space methods of
    [lillard2025vectorspaceintegrationdark, lillard2025waveletharmonicintegrationmethods](@cite),
    as implemented in VSDM [lillard2025vectorspacesdarkmatter](@cite): rates become matrix
    products, and detector orientations rotations of the coefficients.
