module StandardHaloModelsPP

using NaturallyUnitful, Measurements
using Measurements: measurement, value, uncertainty
import ConstructionBase


_isvelocity(x) = dimension(x) == dimension(u"km/s")
_isdensity(x) = dimension(x) in (dimension(u"GeV/cm^3"), dimension(u"kg/m^3"))

export HaloModel
"""
    HaloModel(ρ₀::Union{Unitful.Pressure,Unitful.Density}, v₀::Unitful.Velocity, vₑ::Unitful.Velocity, β::Real, η::Real[, VLSR::Unitful.Velocity])

An SHM++-parameterized halo model, per the definition in
[N. Wyn Evans, Ciaran A.J. O'Hare, Christopher McCabe, *Phys. Rev. D.* **99**, 023012 (2019)](@cite PhysRevD.99.023012)
([arXiv](@cite evans2018shm)).

- `localdensity`: `ρ₀`, local DM mass density: with mass in energy units (e.g.
    `u"GeV/cm^3"`, hence a `Unitful` pressure), or as a mass density 
    (e.g. [`UnitfulAstro`](https://juliaastro.org/UnitfulAstro.jl/stable/): `u"Msun/pc^3"`); 
    units are stored as given.
- `circspeed`: `v₀`, circular speed of the flat asymptote of the galactic rotation curve,
    which sets the DM velocity dispersion ``\\sigma = v_0 / \\sqrt{2}``.
- `escspeed`: `vₑ`, galactic escape speed.
- `anisotropy`: `β`, Gaia sausage anisotropy, `0 ≤ β < 1`.
- `sausagefrac`: `η`, Gaia sausage prominence, `0 ≤ η ≤ 1`.

- `lsrspeed`: `VLSR`, speed of the Local Standard of Rest, which sets the boost to the
    lab frame. By default (omitted, or `nothing`) it is `v₀`: the SHM's isothermal halo has
    a flat rotation curve with circular speed ``v_0 = \\sqrt{2}\\sigma``, so the two coincide
    [PhysRevD.99.023012](@cite). Give it to boost with a different speed than the
    dispersion implies. `m.VLSR` is always the speed in use, even when tied to `v₀`.
"""
struct HaloModel{Tρ,Tv,Tβ,Tl<:Union{Nothing,Tv}}
    localdensity::Tρ
    circspeed::Tv
    escspeed::Tv
    anisotropy::Tβ
    sausagefrac::Tβ
    lsrspeed::Tl # `nothing`: the same as `circspeed`

    function HaloModel(ρ₀, v₀, vₑ, β, η, VLSR = nothing)
        _isdensity(ρ₀) || throw(
            ArgumentError(
                "Local density ρ₀ must be an energy or mass density (e.g. in GeV/cm^3 " *
                "or Msun/pc^3), got $ρ₀",
            ),
        )
        for (name, v) in (("circular speed v₀", v₀), ("escape speed vₑ", vₑ))
            _isvelocity(v) ||
                throw(ArgumentError("Galactic $name must be a velocity, got $v"))
        end
        VLSR === nothing ||
            _isvelocity(VLSR) ||
            throw(ArgumentError("LSR speed VLSR must be a velocity, got $VLSR"))
        β isa Real ||
            throw(ArgumentError("Sausage anisotropy β must be a plain number, got $β"))
        η isa Real ||
            throw(ArgumentError("Sausage fraction η must be a plain number, got $η"))

        ρ₀ > zero(ρ₀) ||
            throw(ArgumentError("Local density ρ₀ must be positive, got $ρ₀"))
        0 ≤ β < 1 ||
            throw(ArgumentError("Sausage anisotropy β must be 0 ≤ β < 1, got $β"))
        0 ≤ η ≤ 1 ||
            throw(ArgumentError("Sausage fraction η must be 0 ≤ η ≤ 1, got $η"))
        v₀ > zero(v₀) ||
            throw(ArgumentError("Galactic circular speed v₀ must be positive, got $v₀"))
        vₑ > zero(vₑ) ||
            throw(ArgumentError("Galactic escape speed vₑ must be positive, got $vₑ"))
        VLSR === nothing ||
            VLSR > zero(VLSR) ||
            throw(ArgumentError("LSR speed VLSR must be positive, got $VLSR"))

        if VLSR === nothing
            v₀, vₑ = promote(float(v₀), float(uconvert(unit(v₀), vₑ)))
        else
            v₀, vₑ, VLSR = promote(
                float(v₀),
                float(uconvert(unit(v₀), vₑ)),
                float(uconvert(unit(v₀), VLSR)),
            )
        end
        β, η = promote(float(β), float(η))
        return new{typeof(ρ₀),typeof(v₀),typeof(β),typeof(VLSR)}(ρ₀, v₀, vₑ, β, η, VLSR)
    end
end
Base.broadcastable(m::HaloModel) = Ref(m)

_fieldalias(s::Symbol) =
    s === :ρ₀ ? :localdensity :
    s === :v₀ ? :circspeed :
    s === :vₑ ? :escspeed :
    s === :β ? :anisotropy : s === :η ? :sausagefrac : s === :VLSR ? :lsrspeed : s

_lsrspeed(m::HaloModel) = something(getfield(m, :lsrspeed), getfield(m, :circspeed))

@inline function Base.getproperty(m::HaloModel, s::Symbol)
    f = _fieldalias(s)
    return f === :lsrspeed ? _lsrspeed(m) : getfield(m, f)
end
Base.propertynames(m::HaloModel, private::Bool = false) =
    (fieldnames(HaloModel)..., :ρ₀, :v₀, :vₑ, :β, :η, :VLSR)

_fields(m::HaloModel) =
    NamedTuple{fieldnames(HaloModel)}(ntuple(i -> getfield(m, i), fieldcount(typeof(m))))
function _canonical(kwargs)
    nt = NamedTuple()
    for (k, v) in pairs(kwargs)
        f = _fieldalias(k)
        f in fieldnames(HaloModel) || throw(
            ArgumentError(
                "unknown HaloModel field `$k`; expected ρ₀/localdensity, v₀/circspeed, " *
                "vₑ/escspeed, β/anisotropy, η/sausagefrac or VLSR/lsrspeed",
            ),
        )
        haskey(nt, f) && throw(ArgumentError("HaloModel field `$f` given twice"))
        nt = merge(nt, NamedTuple{(f,)}((v,)))
    end
    return nt
end

"""
    HaloModel(; ρ₀, v₀, vₑ, β = 0, η = 0, VLSR = nothing)
    HaloModel(m::HaloModel; kwargs...)

Keyword forms of the [`HaloModel`](@ref) constructor. Each parameter may be named by its
field name or a symbol shorthand (`v₀` or `circspeed`, `vₑ` or `escspeed`, ...).

Also a copying form which copies `m` with the given parameters changed, e.g.
`HaloModel(SHMpp; v₀ = 220u"km/s")`. Unchanged `Measurement` parameters are carried over
as the same objects to maintain correlations. 

This also allows e.g. [Accessors.jl](https://juliaobjects.github.io/Accessors.jl/stable/) 
to work automatically on symbol shorthands: `@set SHMpp.v₀ = 220u"km/s"`.
"""
function HaloModel(; kwargs...)
    p = merge((; anisotropy = 0, sausagefrac = 0, lsrspeed = nothing), _canonical(kwargs))
    for f in (:localdensity, :circspeed, :escspeed)
        haskey(p, f) ||
            throw(ArgumentError("HaloModel requires `$f` (or its symbol alias)"))
    end
    return HaloModel(
        p.localdensity,
        p.circspeed,
        p.escspeed,
        p.anisotropy,
        p.sausagefrac,
        p.lsrspeed,
    )
end
HaloModel(m::HaloModel; kwargs...) = HaloModel(merge(_fields(m), _canonical(kwargs))...)

ConstructionBase.getproperties(m::HaloModel) = _fields(m)
ConstructionBase.setproperties(m::HaloModel, patch::NamedTuple) = HaloModel(m; patch...)

# Sun and Earth ephemerides
include("Ephemeris.jl")

export SHM
"""
The default `HaloModel` for the SHM, following the recommended conventions for reporting
direct-detection results in
[D. Baxter et al., *Eur. Phys. J. C* **81**, 907 (2021)](@cite Baxter:2021pqo)
([arXiv](https://arxiv.org/abs/2105.00599))

| SHM | | |
|-------------------------|---------------|---------------------------------|
| Local DM density        | ``\\rho_0``      | ``0.3\\text{ GeV}\\text{ cm}^{-3}`` |
| Circular rotation speed | ``v_0``         | ``238\\text{ km}\\text{ s}^{-1}``   |
| Escape speed            | ``v_\\text{esc}`` | ``544\\text{ km}\\text{ s}^{-1}``   |
| Velocity distribution   | ``f_R(\\mathbf{v})`` | Gaussian |
(additional SHM++ parameters are zeroed.) [Baxter:2021pqo; Table 1](@cite)

`v₀`, and so the LSR speed of the lab-frame boost, is the reference ephemeris'
[`VLSR`](@ref) (its central value). For the pre-2021 `v₀ = 220` km/s convention see
[`SHMlegacy`](@ref).
"""
const SHM = HaloModel(0.3u"GeV/cm^3", 238u"km/s", 544u"km/s", 0, 0)

export SHMlegacy
"""
A historical SHM, with `v₀ = 220` km/s with lab-frame boosts commensurate 
(otherwise as [`SHM`](@ref)), for comparison with results that predate 
[Baxter:2021pqo](@cite); see e.g. [Lewin:1995rx](@cite) and
[evans2018shm; Table 1](@cite).
"""
const SHMlegacy = HaloModel(0.3u"GeV/cm^3", 220u"km/s", 544u"km/s", 0, 0)

export SHMnoesc
"""
An untruncated variant of [`SHM`](@ref): `vₑ = Inf`, so the Maxwellian has no escape-speed
cutoff and normalizes to unity over all velocities.
"""
const SHMnoesc = HaloModel(0.3u"GeV/cm^3", 238u"km/s", (Inf)u"km/s", 0, 0)

export SHMerr
"""
[`SHM`](@ref) with measurement uncertainties on its two characteristic speeds *included*, 
for sensitivity studies of Halo Model conventions and uncertainties.

| SHMerr | | |
|-------------------------|---------------|---------------------------------|
| Local DM density        | ``\\rho_0``      | ``0.3\\text{ GeV}\\text{ cm}^{-3}`` (fixed) |
| Circular rotation speed | ``v_0``         | ``238.0 \\pm 1.5\\text{ km}\\text{ s}^{-1}`` |
| Escape speed            | ``v_\\text{esc}`` | ``544 \\pm 33\\text{ km}\\text{ s}^{-1}`` |

- `v₀`: using uncertainty quoted but not used in [Baxter:2021pqo; Sec. 3.1.5](@cite). 
- `vₑ`: using the RAVE 90% confidence interval ``498 < v_\\text{esc} < 608`` km/s of
  [Smith:2006ym](@cite) (median 544 km/s, agreeing with the `SHM` convention), symmetrized and converted
  to ``1\\sigma``: ``(608 - 498) / (2 \\times 1.645) \\approx 33`` km/s.
- `ρ₀` is a fixed value rather than a measurement with uncertainty, as it is not a directly measured quantity; 
    set another value (usually in the range 0.2-0.6 GeV/cm³) by creating your own 
    `HaloModel` from this one via [`HaloModel`](@ref).
- ``v_\\text{LSR} = v_0`` as `Measurements.jl` objects so the dispersion
    and the boost are correctly correlated.

!!! note
    [Baxter:2021pqo; Sec. 3.1](@cite) recommends against attaching uncertainties to the SHM
    parameters when reporting results: they "do not accurately reflect the uncertainties in
    the dark matter halo, nor do they represent likelihood distributions that can be
    meaningfully profiled over". [`SHM`](@ref) therefore stays exact (and is the default
    model); use `SHMerr` for sensitivity studies, but not for reporting limits.
    Propagating Measurements uncertainties also has a significant performance impact on integrals.
"""
const SHMerr = HaloModel(0.3u"GeV/cm^3", VLSR, measurement(544.0u"km/s", 33.0u"km/s"), 0, 0)

export SHMpp
"""
The default `HaloModel` for the SHM++.

| SHM++ | | |
|-------------------------|---------------|---------------------------------|
| Local DM density        | ``\\rho_0``      | ``0.55 \\pm 0.17\\text{ GeV}\\text{ cm}^{-3}`` |
| Circular rotation speed | ``v_0``         | ``233 \\pm 3\\text{ km}\\text{ s}^{-1}`` |
| Escape speed            | ``v_\\text{esc}`` | ``528 \\pm 25\\text{ km}\\text{ s}^{-1}`` |
| Sausage anisotropy      | ``\\beta``       | ``0.9\\pm 0.05``                   |
| Sausage fraction        | ``\\eta``        | ``0.2\\pm 0.1``                    |
| Velocity distribution   | ``f_\\chi(\\mathbf{v})`` | Gaussian+Sausage |
[evans2018shm; Table 1](@cite)

As in the defining SHM++ [paper](@cite evans2018shm), the lab-frame boost uses
``v_\\text{LSR} = v_0 = 233 \\pm 3`` km/s, with the boost and the dispersion correlated.
Create your own modified copy if you want a SHM++ halo model with Baxter [Baxter:2021pqo](@cite) parameters
instead.

You may want to strip Measurement uncertainties via `Measurements.value(m)` if you are not
interested in sensitivity questions wrt. the Gaia sausage and just want to use this as a halo model.
"""
const SHMpp = HaloModel(
    measurement(0.55u"GeV/cm^3", 0.17u"GeV/cm^3"),
    measurement(233u"km/s", 3u"km/s"),
    measurement(528u"km/s", 25u"km/s"),
    measurement(0.9, 0.05),
    measurement(0.2, 0.1),
)

# Base.show helpers
const _PRESETS = (; SHM, SHMlegacy, SHMnoesc, SHMerr, SHMpp)
_presetname(m::HaloModel) = findfirst(p -> p === m, _PRESETS) # Symbol, or nothing

const _SHOWFIELDS = (
    (:localdensity, "ρ₀"),
    (:circspeed, "v₀"),
    (:escspeed, "vₑ"),
    (:anisotropy, "β"),
    (:sausagefrac, "η"),
    (:lsrspeed, "VLSR"),
)
_istied(m::HaloModel) = getfield(m, :lsrspeed) === nothing

function Base.show(io::IO, m::HaloModel)
    name = _presetname(m)
    name === nothing || return print(io, name)
    print(io, "HaloModel(")
    shown = (fs for fs in _SHOWFIELDS if !(first(fs) === :lsrspeed && _istied(m)))
    join(io, ("$s = $(getfield(m, f))" for (f, s) in shown), ", ")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", m::HaloModel)
    name = _presetname(m)
    print(io, "HaloModel", name === nothing ? "" : " $name", ":")
    w = maximum(length ∘ string ∘ first, _SHOWFIELDS)
    ws = maximum(textwidth ∘ last, _SHOWFIELDS)
    for (f, s) in _SHOWFIELDS
        val = f === :lsrspeed && _istied(m) ? "v₀" : getfield(m, f)
        print(io, "\n  ", rpad(string(f), w), "  ", rpad(s, ws), " = ", val)
    end
end

"""
    Measurements.value(m::HaloModel)

Copy `m` with every `Measurement` parameter replaced by its central value.
"""
function Measurements.value(m::HaloModel)
    l = getfield(m, :lsrspeed)
    return HaloModel(
        value(m.ρ₀),
        value(m.v₀),
        value(m.vₑ),
        value(m.β),
        value(m.η),
        l === nothing ? nothing : value(l),
    )
end

# Velocity distribution definitions
include("VelocityDistributions.jl")

# Analytic integral forms
include("AnalyticSHMIntegrals.jl")

# Helper/utils and Integrated marginal distributions
include("MarginalDistributions.jl")

# General velocity integration interface
include("VelocityProblems.jl")

# Convenience and quick user functions
include("VelocityIntegration.jl")

# Precompile potentially hot entry points.
using PrecompileTools
@compile_workload begin
    v = [300.0, 0.0, 0.0]u"km/s"
    tol = (reltol = 1e-2, abstol = 1e-4)

    # Pointwise distributions, exact and Measurement, galactic and lab frame.
    fχ(SHM, v)
    fχ(SHMpp, v)
    fR(SHMpp, v)
    fS(SHMpp, v)
    fχlab(SHM, T0, v)
    fχlab(SHMpp, T0, v)

    # Analytic fast paths (η = 0).
    fχspeed(SHM, 300u"km/s")
    g(SHM, 200u"km/s")
    h(SHM, 200u"km/s")

    # Numerical paths on central values (the integrators' default): the angular cubature
    # (_labframe/_fχspeed), the planar quadrature and the 3D cubature (_velocityproblem).
    fχspeed(SHM, 300u"km/s"; analytic = false, tol...)
    fχspeed(SHMpp, 300u"km/s"; tol...)
    fχplanar(SHM, 100u"km/s", -50u"km/s"; tol...)
    g(SHM, 200u"km/s"; analytic = false, tol...)
    g(SHMpp, 200u"km/s"; tol...)

    # Specializations for propagating Measurement uncertainties.
    fχspeed(SHMpp, 300u"km/s"; uncertainty = true, tol...)
    fχplanar(SHMpp, 100u"km/s", -50u"km/s"; uncertainty = true, tol...)
    g(SHMpp, 200u"km/s"; uncertainty = true, tol...)

    # Frozen speed distributions, as used in various integrands.
    fχspeed(SHM, T0)(300u"km/s")
    fχspeed(SHMpp, T0; tol...)(300u"km/s")

    # velint/speedint/HaloIntegralProblem cannot be precompiled: 
    # their kernel is a user closure, of a distinct type per call site,
    # so cannot be generically precompiled.
end

end # module StandardHaloModelsPP
