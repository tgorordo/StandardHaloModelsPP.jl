using NaturallyUnitful, Dates
using Measurements
using Measurements: value, uncertainty
using LinearAlgebra, Integrals, StaticArrays

export velint
"""
    velint([m::HaloModel = SHM,] [t::TimeLike = T0,] kk; uncertainty = false, vmin = 0u"km/s", p = (), alg = HCubatureJL(), reltol = 1e-4, abstol = 0.0, kwargs...)

Arbitrary convolution ``\\int_{v>v_\\text{min}} f_{\\chi,\\text{lab}}(\\mathbf{v})\\,
k(\\mathbf{v},p)\\,d^3v`` of a kernel `kk(v, p)` (or `kk(v)`) with the lab-frame velocity distribution:
`solve(HaloIntegralProblem(m, t, kk; ...), alg; ...).u`. See
[`HaloIntegralProblem`](@ref) for the problem itself, e.g. to choose an algorithm, inspect
the error estimate or sweep `t`.

$(_ABSTOL_DOC)

$(_UNCERTAINTY_DOC)

See also [`speedint`](@ref), which convolves against the marginalized speed distribution.
"""
velint(
    m::HaloModel,
    t::TimeLike,
    kk;
    uncertainty::Bool = false,
    vmin = 0u"km/s",
    p = (),
    alg = _CUBALG,
    reltol = 1e-4,
    abstol = 0.0,
    kwargs...,
) = _solvechecked(
    HaloIntegralProblem(m, t, kk; uncertainty, vmin, p),
    alg,
    :velint;
    reltol,
    abstol,
    kwargs...,
)
velint(t::TimeLike, kk; kwargs...) = velint(SHM, t, kk; kwargs...)
velint(m::HaloModel, kk; kwargs...) = velint(m, T0, kk; kwargs...)
velint(kk; kwargs...) = velint(SHM, T0, kk; kwargs...)

export g
"""
    g([m::HaloModel = SHM,] [t::TimeLike = T0,] vmin::Unitful.Velocity; uncertainty = false, analytic = true, alg = HCubatureJL(), reltol = 1e-4, abstol = 0, kwargs...)

The Mean Reciprocal Speed ``g(v_\\text{min}, t) = \\int_{v>v_\\text{min}}
\\frac{f(\\mathbf{v} + \\mathbf{v}_E(t))}{v}\\,d^3v``.

$(_ANALYTIC_DOC) $(_ABSTOL_DOC)

$(_UNCERTAINTY_DOC)
"""
function g(
    m::HaloModel,
    t::TimeLike,
    vmin::Unitful.Velocity;
    uncertainty::Bool = false,
    kwargs...,
)
    return _g_impl(_precision(m, uncertainty), t, vmin; kwargs...)
end
function _g_impl(
    m::HaloModel,
    t::TimeLike,
    vmin::Unitful.Velocity;
    analytic = true,
    alg = _CUBALG,
    reltol = 1e-4,
    abstol = 0.0,
    kwargs...,
)
    analytic && _hasanalytic(m, t) && return _gSHM(m, t, vmin)
    # radial kernel: the velocity integrand already holds the speed, so norm(v) is not
    # recomputed
    prob, uk = _velocityproblem(m, t, (vxyz, vv, pp) -> 1 / vv; vmin)
    sol = _solve(prob, uk, alg; reltol, abstol, kwargs...)
    _checkretcode(sol.retcode, :g)
    return sol.u
end
g(t::TimeLike, vmin::Unitful.Velocity; kwargs...) = g(SHM, t, vmin; kwargs...)
g(m::HaloModel, vmin::Unitful.Velocity; kwargs...) = g(m, T0, vmin; kwargs...)
g(vmin::Unitful.Velocity; kwargs...) = g(SHM, T0, vmin; kwargs...)

export h
"""
    h([m::HaloModel = SHM,] [t::TimeLike = T0,] vmin::Unitful.Velocity; uncertainty = false, analytic = true, alg = HCubatureJL(), reltol = 1e-4, abstol = 0, kwargs...)

The Mean Speed ``h(v_\\text{min}, t) = \\int_{v>v_\\text{min}} v\\,
f(\\mathbf{v} + \\mathbf{v}_E(t))\\,d^3v``.

$(_ANALYTIC_DOC) $(_ABSTOL_DOC)

$(_UNCERTAINTY_DOC)
"""
function h(
    m::HaloModel,
    t::TimeLike,
    vmin::Unitful.Velocity;
    uncertainty::Bool = false,
    kwargs...,
)
    return _h_impl(_precision(m, uncertainty), t, vmin; kwargs...)
end
function _h_impl(
    m::HaloModel,
    t::TimeLike,
    vmin::Unitful.Velocity;
    analytic = true,
    alg = _CUBALG,
    reltol = 1e-4,
    abstol = 0.0,
    kwargs...,
)
    analytic && _hasanalytic(m, t) && return _hSHM(m, t, vmin)
    prob, uk = _velocityproblem(m, t, (vxyz, vv, pp) -> vv; vmin)
    sol = _solve(prob, uk, alg; reltol, abstol, kwargs...)
    _checkretcode(sol.retcode, :h)
    return sol.u
end
h(t::TimeLike, vmin::Unitful.Velocity; kwargs...) = h(SHM, t, vmin; kwargs...)
h(m::HaloModel, vmin::Unitful.Velocity; kwargs...) = h(m, T0, vmin; kwargs...)
h(vmin::Unitful.Velocity; kwargs...) = h(SHM, T0, vmin; kwargs...)

export speedint
"""
    speedint([m::HaloModel = SHM,] [t::TimeLike = T0,] kk; uncertainty = false, vmin = 0u"km/s", vmax = nothing, p = (), alg = QuadGKJL(), reltol = 1e-6, abstol = 0.0, analytic = true, inneralg = HCubatureJL(), innerreltol = reltol / 10, innerabstol = abstol / 10, kwargs...)

Convolution ``\\int k(v,p)\\,f_{\\chi,\\text{lab}}^{\\text{speed}}(v)\\,dv`` of a scalar kernel
`kk(v, p)` (or `kk(v)`) with the marginalized speed distribution:
`solve(HaloIntegralProblem(m, t, kk; measure = :speed, ...), alg; ...).u`. See
[`HaloIntegralProblem`](@ref) for the problem itself.

`abstol` defaults to `0`, a purely relative tolerance, as for [`velint`](@ref). It may be a
`Quantity` in the units of the result; a plain number is taken in the integrand's internal
(unit-stripped) units. The default `innerabstol` is a 1/10 of plain `abstol`, and
`0` (the inner quadrature then controlled by `innerreltol` alone) for the default or a
`Quantity` `abstol`, since the inner integral has different units; pass `innerabstol`
explicitly, as a number or a `Quantity` in the units of [`fχspeed`](@ref), to override.
Further `kwargs` are passed to the outer Integrals.jl `solve`.

$(_UNCERTAINTY_DOC)

See also [`velint`](@ref), which convolves against the full 3-D velocity distribution.
"""
function speedint(
    m::HaloModel,
    t::TimeLike,
    kk;
    uncertainty::Bool = false,
    vmin = 0u"km/s",
    vmax = nothing,
    p = (),
    alg = _QUADALG,
    reltol = 1e-6,
    abstol = 0.0,
    analytic = true,
    inneralg = _CUBALG,
    innerreltol = reltol / 10,
    innerabstol = abstol isa Quantity ? 0.0 : abstol / 10,
    kwargs...,
)
    prob = HaloIntegralProblem(
        m,
        t,
        kk;
        measure = :speed,
        uncertainty,
        vmin,
        vmax,
        p,
        analytic,
        inneralg,
        innerreltol,
        innerabstol,
    )
    return _solvechecked(prob, alg, :speedint; reltol, abstol, kwargs...)
end
speedint(t::TimeLike, kk; kwargs...) = speedint(SHM, t, kk; kwargs...)
speedint(m::HaloModel, kk; kwargs...) = speedint(m, T0, kk; kwargs...)
speedint(kk; kwargs...) = speedint(SHM, T0, kk; kwargs...)

for f in (:velint, :speedint)
    @eval begin
        $f(m::HaloModel; kwargs...) = throw(MethodError($f, (m,)))
        $f(t::TimeLike; kwargs...) = throw(MethodError($f, (t,)))
        $f(m::HaloModel, t::TimeLike; kwargs...) = throw(MethodError($f, (m, t)))
    end
end
