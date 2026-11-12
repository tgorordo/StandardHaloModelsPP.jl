using NaturallyUnitful, Dates
using Measurements
using Measurements: value, uncertainty
using LinearAlgebra, Integrals, StaticArrays

# Integration problems ============================================================

"""
    _radialsupport(Q, V, n, vmin) -> (lo, Δ)

Find the radial speed range with nonzero support from the halo.

`V` Earth galactic velocity,
`n` lab-frame DM velocity unit vector, 
`Q = vₑ² - |V|²`, 
`vmin` the kinematic lower cutoff.

Returns (`lo`, `Δ`), lower bound and width of the range of speeds with support along `n`.
(when `Δ==0` `lo` is moved `oneunit` away from zero to avoid coordinate singularity).
"""
@inline function _radialsupport(Q, V, n, vmin)
    c = n ⋅ V
    D = Q + c^2

    r = value(D) > zero(value(D)) ? sqrt(D) : zero(D) / oneunit(c)
    lo = max(promote(vmin, -c - r)...)
    Δ = max(-c + r - lo, zero(lo))

    return (iszero(value(Δ)) ? lo + oneunit(lo) : lo), Δ
end

"""
    _capangle(Q, VEn, vmin) -> γ0

Polar angle from `V̂` beyond which a lab-frame direction `n̂` has support at some speed
`v ≥ vmin`, i.e. `|v n̂ + V| < vₑ`; all of the integrand's support lies in `γ ∈ [γ0, π]`.

`Q = vₑ² - |V|²`, `VEn = |V|`, `vmin` the kinematic lower cutoff. Along `n̂` at angle `γ`
the condition is `cos γ < (Q - v²)/(2v|V|)`, largest at `v = max(vmin, √(-Q))`.
"""
function _capangle(Q, VEn, vmin)
    q, V = value(Q), value(VEn)
    vs = max(promote(value(vmin), sqrt(max(-q, zero(q))))...)
    iszero(vs) && return 0.0 # vmin = 0 inside the escape sphere: every direction
    return acos(clamp(uconvert(NoUnits, (q - vs^2) / (2vs * V)), -1.0, 1.0))
end

# argument wrapping for kernels
struct _Kernel{F} <: Function
    f::F
end
@inline (k::_Kernel)(v, ::Tuple{}) = applicable(k.f, v) ? k.f(v) : k.f(v, ())
@inline (k::_Kernel)(v, p) = applicable(k.f, v, p) ? k.f(v, p) : k.f(v)

_takes(f, n) = any(m -> m.nargs - 1 == n || (m.isva && m.nargs - 2 ≤ n), methods(f))
function _kernel2(kk)
    _takes(kk, 1) ||
        _takes(kk, 2) ||
        throw(ArgumentError("a kernel must be callable like f(v) or f(v, p); got $kk"))
    return _Kernel(kk)
end

"""
    _velocityproblem(m::HaloModel, t, kk; vmin = 0u"km/s", p = ())

Construct the 3D velocity integrand over (s, w, α) ∈ [0,1] × [0,1] × [0,2π].

Kernel `kk` receives `(vxyz, vv, p)` with `vv = norm(vxyz)`.

Returns the unit-free problem and the unit of its integral. As in `_fχspeed`, directions
n̂ are polar about V̂, at angle γ, and restricted to the cap γ ∈ [γ0, π] holding all the
support above `vmin` (see `_capangle`), rescaled onto w ∈ [0,1]; α is the full azimuth,
since the kernel need not be symmetric. The radial coordinate is rescaled onto s ∈ [0,1]
spanning the support [lo, lo + Δ] along n̂. So the domain is the support, and adaptive
cubature can't miss even a small cap near `vmaxlab`.

The support along n̂ shrinks to nothing at γ0, where the integrand vanishes, so the
central-value bound γ0 needs no first-order boundary term for `Measurement` models.
"""
function _velocityproblem(
    m::HaloModel,
    t::TimeLike,
    kk::K;
    vmin = 0u"km/s",
    p = (),
) where {K}
    _checkanisotropy(m)
    V = _vE(m, t)              # hoist
    fp = prepare(m)            # hoist
    Q = _vesc_eff(m)^2 - V ⋅ V # direction-independent

    VEn = norm(V)
    ê3 = value.(V ./ VEn)
    ê1, ê2 = _orthobsym(ê3)
    γ0 = _capangle(Q, VEn, vmin)
    Δγ = π - γ0

    k(swα, q) =
        let s = swα[1], w = swα[2], α = swα[3], fp = q[1], V = q[2], pp = q[3]
            sγ, cγ = sincos(γ0 + w * Δγ)
            sα, cα = sincos(α)
            n = (sγ * cα) .* ê1 .+ (sγ * sα) .* ê2 .+ cγ .* ê3
            lo, Δ = _radialsupport(Q, V, n, vmin) # Δ = 0, lo ≠ 0 where the support is empty
            vv = lo + s * Δ
            vxyz = vv .* n

            # Jacobians: vv²sin(γ) from the spherical measure, Δγ and Δ from the rescalings.
            return fp(vxyz + V) * kk(vxyz, vv, pp) * vv^2 * sγ * Δγ * Δ
        end
    uk = unit(k(SA[0.5, 0.5, π], (fp, V, p)))

    domain = (SA[0.0, 0.0, 0.0], SA[1.0, 1.0, 2π])
    return IntegralProblem{false}((x, q) -> ustrip(uk, k(x, q)), domain, (fp, V, p)), uk
end

"""
    _speedproblem(m::HaloModel, t, kk; vmin = 0u"km/s", vmax = nothing, p = (), kwargs...)

Construct the 1D speed integrand over v ∈ [vmin, vmax], in the units of vmin.

Returns the unit-free problem and the unit of its integral.
"""
function _speedproblem(
    m::HaloModel,
    t::TimeLike,
    kk::K;
    vmin = 0u"km/s",
    vmax = nothing,
    p = (),
    analytic = true,
    inneralg = _CUBALG,
    innerreltol = 1e-7,
    innerabstol = 0.0,
) where {K}
    vmax = something(vmax, vmaxlab(m, t))

    inner = _frozenspeed(
        m,
        t;
        analytic,
        alg = inneralg,
        reltol = innerreltol,
        abstol = innerabstol,
    )
    return _speedradial(inner, kk, vmin, vmax, p)
end

function _speedradial(inner::F, kk::K, vmin, vmax, p) where {F,K}
    ui = unit(vmin)
    k(v, q) =
        let pp = q[1]
            return inner(v) * kk(v, pp) * ui
        end

    # probe at a non-degenerate speed: vmin may be 0, where kernels like 1/v are singular.
    uk = unit(k(300u"km/s", (p,)))

    domain = (value(ustrip(ui, vmin)), value(ustrip(ui, vmax)))
    return IntegralProblem{false}((v, q) -> ustrip(uk, k(v * ui, q)), domain, (p,)), uk
end

export HaloIntegralSolution
"""
    HaloIntegralSolution

The result of `solve`-ing a [`HaloIntegralProblem`](@ref).

- `u`: the integral, with units.
- `resid`: the solver's error estimate, in the same units, or `nothing` if the algorithm
  does not provide one.
- `retcode`: the SciML return code; checked by `SciMLBase.successful_retcode(sol)`.
- `sol`: the underlying (unit-free) Integrals.jl solution.
"""
struct HaloIntegralSolution{T,R,C,S}
    u::T
    resid::R
    retcode::C
    sol::S
end
SciMLBase.successful_retcode(s::HaloIntegralSolution) =
    SciMLBase.successful_retcode(s.retcode)
Base.show(io::IO, s::HaloIntegralSolution) =
    print(io, "HaloIntegralSolution(u = ", s.u, ", retcode = ", s.retcode, ")")

_withunit(::Nothing, u) = nothing
_withunit(x, u) = x * u

# Solve a unit-free problem whose integral has unit `uk`; shared by `solve` and g/h.
function _solve(prob, uk, alg; reltol, abstol, kwargs...)
    at = _abstol(abstol, uk)
    sol = solve(prob, alg; reltol, abstol = at, kwargs...)
    rc = _retcode(sol, reltol, at)
    return HaloIntegralSolution(sol.u * uk, _withunit(sol.resid, uk), rc, sol)
end

export HaloIntegralProblem, solve, remake
"""
    HaloIntegralProblem([m::HaloModel = SHM,] [t::TimeLike = T0,] kk; measure = :velocity, uncertainty = false, kwargs...)

An integral problem with a kernel `kk` against the lab-frame velocity distribution of `m` at time
`t`, set up for Integrals.jl. The kernel `kk` may be any callable (a function, or a
callable object) taking `(v, p)` or just `(v)`: with the default `p = ()` it is called as
`kk(v)` if it has a method for that (so `norm` or `inv ∘ norm` work), otherwise as
`kk(v, ())`; with a `p`, as `kk(v, p)` if it has a method for that, otherwise as `kk(v)`.
`solve` it with any Integrals.jl algorithm:

```julia
using LinearAlgebra: norm
using Integrals: HCubatureJL
using Dates: DateTime, Day

prob = HaloIntegralProblem(SHM, v -> 1 / norm(v); vmin = 300u"km/s")
sol = solve(prob)
sol = solve(prob, HCubatureJL(); reltol = 1e-6, maxiters = 10^6)
sol.u                                     # integral solution with units
remake(prob; t = T0 + Day(90))            # integral at another time
```

## Measures

- `measure = :velocity`: ``\\int_{v>v_\\text{min}} f_{\\chi,\\text{lab}}(\\mathbf{v})\\,
  k(\\mathbf{v}, p)\\,d^3v``, with `kk(v::AbstractVector, p)`. Keywords: `vmin`, `p`. The
  problem is 3D over rescaled spherical coordinates ``(s, w, \\alpha) \\in [0,1] \\times
  [0,1] \\times [0,2\\pi]``, polar about the Earth's velocity, spanning exactly the
  distribution's support above ``v_\\text{min}``: the escape-speed boundary is the domain
  boundary, and near ``v_\\text{max}`` the domain shrinks onto the small cap of directions
  that still reach it; default algorithm `HCubatureJL()`, `reltol = 1e-4`, `abstol = 0`. See [`velint`](@ref).

- `measure = :speed`: ``\\int_{v_\\text{min}}^{v_\\text{max}} k(v, p)\\,
  f_{\\chi,\\text{lab}}^{\\text{speed}}(v)\\,dv``, with `kk(v::Unitful.Velocity, p)`.
  Keywords: `vmin`, `vmax` (default ``v_\\text{esc} + V_E``), `p`, `analytic`, and the
  inner angular quadrature's `inneralg`, `innerreltol` (`1e-7`), `innerabstol` (`0`).
  The problem is 1D over ``v``; the inner integral, [`fχspeed`](@ref), is fixed at
  construction, and is analytic when ``\\eta = 0`` and ``V_E < v_\\text{esc}``. Default
  algorithm `QuadGKJL()`, `reltol = 1e-6`, `abstol = 0`. See [`speedint`](@ref).

## Solving

`solve(prob[, alg]; reltol, abstol, kwargs...)` returns a [`HaloIntegralSolution`](@ref);
further `kwargs` go to Integrals.jl's `solve`. `abstol` may be a `Quantity` in the units of
the result; a plain number is taken in the integrand's internal (unit-stripped) units.
`remake(prob; kwargs...)` rebuilds the problem with any constructor argument changed
(`m`, `t`, `kernel`, `measure`, `uncertainty`, or the measure's keywords), at a few percent
of the cost of a solve.

## Sweeps

A sweep is a SciML `EnsembleProblem` whose `prob_func` remakes the problem per point; it
returns an `EnsembleSolution` of [`HaloIntegralSolution`](@ref)s, and runs on threads with
`EnsembleThreads()` (start Julia with `-t auto`):

```julia
ts = [DateTime(2020, 1, 1) + Day(d) for d in 0:364]
ens = EnsembleProblem(prob; prob_func = (p, ctx) -> remake(p; t = ts[ctx.sim_id]))
sols = solve(ens, nothing, EnsembleThreads(); trajectories = length(ts))
```

Pass `nothing` as the algorithm for the problem's default. Kernel parameters sweep the
same way, e.g. `remake(p; p = Es[ctx.sim_id])` for a kernel `kk(v, p)`.

## Uncertainties

By default only the central values of `m`'s parameters (and of the ephemeris) are integrated.
`uncertainty = true` additionally propagates a `Measurement` model's uncertainties. That is
typically 10-130x slower, and only works with pure-Julia algorithms (e.g. `QuadGKJL`,
`HCubatureJL`): C-backed ones (Cuba, Cubature) need plain floating-point integrands. The
defaults measure error on central values (`norm = value ∘ norm`); a user-supplied
adaptive algorithm without that norm still works, but possibly more slowly.
"""
struct HaloIntegralProblem{P<:IntegralProblem,U,A,S}
    prob::P         # unit-free Integrals.jl problem
    unit::U         # unit of the integral: result = solution.u * unit
    alg::A          # default algorithm for the measure
    reltol::Float64 # default tolerances for the measure
    abstol::Float64
    spec::S         # constructor arguments, for remake and show
end

function HaloIntegralProblem(
    m::HaloModel,
    t::TimeLike,
    kk;
    measure::Symbol = :velocity,
    uncertainty::Bool = false,
    kwargs...,
)
    spec = (; m, t, kernel = kk, measure, uncertainty, kwargs = NamedTuple(kwargs))
    return _haloproblem(
        Val(measure),
        _precision(m, uncertainty),
        t,
        _kernel2(kk),
        spec;
        kwargs...,
    )
end
HaloIntegralProblem(t::TimeLike, kk; kwargs...) = HaloIntegralProblem(SHM, t, kk; kwargs...)
HaloIntegralProblem(m::HaloModel, kk; kwargs...) = HaloIntegralProblem(m, T0, kk; kwargs...)
HaloIntegralProblem(kk; kwargs...) = HaloIntegralProblem(SHM, T0, kk; kwargs...)
# A model and/or time but no kernel: a missing argument, not a kernel (as for velint).
HaloIntegralProblem(m::HaloModel; kwargs...) = throw(MethodError(HaloIntegralProblem, (m,)))
HaloIntegralProblem(t::TimeLike; kwargs...) = throw(MethodError(HaloIntegralProblem, (t,)))
HaloIntegralProblem(m::HaloModel, t::TimeLike; kwargs...) =
    throw(MethodError(HaloIntegralProblem, (m, t)))

_haloproblem(::Val{:velocity}, m, t, kk, spec; kwargs...) = HaloIntegralProblem(
    _velocityproblem(m, t, (vxyz, vv, pp) -> kk(vxyz, pp); kwargs...)...,
    _CUBALG,
    1e-4,
    0.0,
    spec,
)
_haloproblem(::Val{:speed}, m, t, kk, spec; kwargs...) =
    HaloIntegralProblem(_speedproblem(m, t, kk; kwargs...)..., _QUADALG, 1e-6, 0.0, spec)
_haloproblem(::Val{M}, args...; kwargs...) where {M} =
    throw(ArgumentError("measure must be :velocity or :speed, got :$M"))

function _solvechecked(hp, alg, name; kwargs...)
    sol = solve(hp, alg; kwargs...)
    _checkretcode(sol.retcode, name)
    return sol.u
end

SciMLBase.solve(
    hp::HaloIntegralProblem,
    alg = nothing;            # use problem default
    reltol = hp.reltol,
    abstol = hp.abstol,
    kwargs...,
) = _solve(hp.prob, hp.unit, something(alg, hp.alg); reltol, abstol, kwargs...)

using Integrals.SciMLBase:
    EnsembleProblem, EnsembleSerial, EnsembleThreads, EnsembleDistributed
export EnsembleProblem, EnsembleSerial, EnsembleThreads, EnsembleDistributed

function SciMLBase.remake(
    hp::HaloIntegralProblem;
    m = hp.spec.m,
    t = hp.spec.t,
    kernel = hp.spec.kernel,
    measure = hp.spec.measure,
    uncertainty = hp.spec.uncertainty,
    kwargs...,
)
    return HaloIntegralProblem(
        m,
        t,
        kernel;
        measure,
        uncertainty,
        hp.spec.kwargs...,
        kwargs...,
    )
end

function Base.show(io::IO, hp::HaloIntegralProblem)
    (; m, t, measure, uncertainty, kwargs) = hp.spec
    print(io, "HaloIntegralProblem(", m, ", ", t, ", kernel; measure = :", measure)
    uncertainty && print(io, ", uncertainty = true")
    for (k, v) in pairs(kwargs)
        print(io, ", ", k, " = ", v)
    end
    print(io, ")")
end
