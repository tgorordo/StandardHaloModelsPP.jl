using NaturallyUnitful, Dates
using Measurements
using Measurements: value, uncertainty
using LinearAlgebra, Integrals, StaticArrays

## Integration helpers needed by: 
# MarginalDistributions.jl (inc. first),
# VelocityProblems.jl, and 
# VelocityIntegration.jl.

_abstol(abstol::Quantity, uk; scale = 1) = value(ustrip(uk, abstol)) / scale
_abstol(abstol, uk; scale = 1) = abstol

_precision(m::HaloModel, uncertainty::Bool) = uncertainty ? m : value(m)

_hasanalytic(m::HaloModel, t) = iszero(m.η) && value(_VE(m, t)) < value(m.vₑ)

# Escape speed for integrators to cut off: vₑ, or for an untruncated model (vₑ = ∞) an
# effective n·v₀, beyond which the Maxwellian exp(-v²/v₀²) is below e^{-n²} (e^{-36}).
_vesc_eff(m::HaloModel; n = 6) = isfinite(m.vₑ) ? m.vₑ : n * m.v₀

function _retcode(sol, reltol, abstol)
    SciMLBase.successful_retcode(sol.retcode) || return sol.retcode
    sol.resid === nothing && return sol.retcode
    tol = max(abstol, reltol * abs(value(sol.u)))
    return value(sol.resid) ≤ tol * (1 + 1e-6) ? sol.retcode : SciMLBase.ReturnCode.MaxIters
end

function _checkretcode(retcode, name)
    SciMLBase.successful_retcode(retcode) && return nothing
    @warn "$name: the integral did not converge (retcode = $retcode), so the result may " *
          "be inaccurate. Loosen `reltol`/`abstol` or raise `maxiters`; " *
          "`solve(HaloIntegralProblem(...))` returns the full solution."
    return nothing
end

_enorm(x::Number) = abs(value(x))
_enorm(x) = value(norm(x))

# Default algorithms. 
const _CUBALG = HCubatureJL(norm = _enorm)
const _QUADALG = QuadGKJL(norm = _enorm)

# Docstring paragraphs common to the integrators; interpolated to be consistent.
const _ABSTOL_DOC = """
`abstol` defaults to `0`, a purely relative tolerance. It may be a `Quantity` in the units of
the result; a plain number is taken in the integrand's internal (unit-stripped) units.
Further `kwargs` are passed to Integrals.jl's `solve`."""
const _UNCERTAINTY_DOC = """
By default the central values of `m`'s parameters (and of the ephemeris) are integrated.
Pass `uncertainty = true` to propagate a `Measurement` model's uncertainties
(typically 30-90x slower due to derivative tracking/propagation)."""
const _ANALYTIC_DOC = """
When ``\\eta = 0`` and the Earth is inside the escape sphere (``V_E < v_\\text{esc}``, true of
all the presets) an exact analytic form is used; pass `analytic=false` to force numerical
integration."""

## Marginal Distributions ====================

export fχplanar
"""
    fχplanar([m::HaloModel = SHM,] [t::TimeLike = T0,] vr::Unitful.Velocity, vθ::Unitful.Velocity; uncertainty = false, alg = QuadGKJL(), reltol = 1e-6, abstol = 0, kwargs...)

The galactic-planar velocity distribution, marginalized over ``v_z \\sim v_\\phi``:
```math
f_{\\chi,\\text{lab}}^{\\text{planar}}(v_r, v_\\theta) = \\int f_{\\chi,\\text{lab}}(\\mathbf{v})\\, dv_z
```

$(_ABSTOL_DOC)

$(_UNCERTAINTY_DOC)

See also [`fχspeed`](@ref) for the marginalization over solid angle.
"""
function fχplanar(
    m::HaloModel,
    t::TimeLike,
    vr::Unitful.Velocity,
    vθ::Unitful.Velocity;
    uncertainty::Bool = false,
    kwargs...,
)
    return _fχplanar_impl(_precision(m, uncertainty), t, vr, vθ; kwargs...)
end
function _fχplanar_impl(
    m::HaloModel,
    t::TimeLike,
    vr::Unitful.Velocity,
    vθ::Unitful.Velocity;
    alg = _QUADALG,
    reltol = 1e-6,
    abstol = 0.0,
    kwargs...,
)
    ui = unit(vθ)
    vr, vθ = promote(uconvert(ui, vr), vθ)

    fp = prepare(m)
    V = _vE(m, t)
    k(vϕ, q) =
        let fp = q[1], V = q[2], vr = q[3], vθ = q[4]
            fp(SA[vr, vθ, vϕ] + V) * ui
        end
    kprobe = k(100u"km/s", (fp, V, vr, vθ))
    uk = unit(kprobe)

    D = _vesc_eff(m)^2 - (vr + V[1])^2 - (vθ + V[2])^2
    value(ustrip(ui^2, D)) ≤ 0 && return zero(kprobe)

    r = sqrt(D)
    lom, him = -V[3] - r, -V[3] + r # the limits, with their uncertainties
    lo = value(ustrip(ui, lom))
    hi = value(ustrip(ui, him))

    domain = (lo, hi)
    prob = IntegralProblem((vϕ, q) -> ustrip(uk, k(vϕ * ui, q)), domain, (fp, V, vr, vθ))
    at = _abstol(abstol, uk)
    sol = solve(prob, alg; reltol, abstol = at, kwargs...)
    _checkretcode(_retcode(sol, reltol, at), :fχplanar)
    I = sol.u * uk
    _ismeas(typeof(D)) || return I

    # The quadrature runs between the central-value limits; to first order their
    # uncertainty adds f(hi)·δhi - f(lo)·δlo, where δ = limit - its central value (zero
    # value, so only the propagated uncertainty changes). f is evaluated just inside each
    # limit, since exactly on the escape sphere rounding can land outside it, where fχ = 0.
    ε = 1e-9 * (hi - lo)
    fb(z) = fp(SA[vr, vθ, z*ui] + V)
    return I + fb(hi - ε) * (him - hi * ui) - fb(lo + ε) * (lom - lo * ui)
end
fχplanar(t::TimeLike, vr::Unitful.Velocity, vθ::Unitful.Velocity; kwargs...) =
    fχplanar(SHM, t, vr, vθ; kwargs...)
fχplanar(m::HaloModel, vr::Unitful.Velocity, vθ::Unitful.Velocity; kwargs...) =
    fχplanar(m, T0, vr, vθ; kwargs...)
fχplanar(vr::Unitful.Velocity, vθ::Unitful.Velocity; kwargs...) =
    fχplanar(T0, vr, vθ; kwargs...)

export fχspeed
function _fχspeed(
    L,
    v;
    alg = _CUBALG,
    reltol = 1e-4,
    abstol = 0.0,
    checkretcode = false,
    kwargs...,
)
    (; fp, Vp, VEn, ê1, ê2, ê3, vesc2) = L
    # Opening angle of the escape-speed boundary: directions with γ > γ0 lie inside the
    # escape sphere. The quadrature runs at the central value; uncertainty enters
    # through the boundary term added at the end.
    cm = uconvert(NoUnits, (vesc2 - v^2 - VEn^2) / (2v * VEn))
    c = clamp(value(cm), -1.0, 1.0)
    γ0 = acos(c)
    Δγ = π - γ0

    k(sα, q) =
        let s = sα[1], α = sα[2], fp = q[1], Vp = q[2], vv = q[3]
            γ = γ0 + s * Δγ
            sγ, cγ = sincos(γ)
            sα_, cα_ = sincos(α)
            n = (sγ * cα_) .* ê1 .+ (sγ * sα_) .* ê2 .+ cγ .* ê3
            # Jacobians: sin(γ) from dΩ, Δγ from the polar rescaling.
            return vv^2 * fp(vv .* n + Vp) * sγ * Δγ
        end
    uk = unit(k(SA[0.5, 0.0], (fp, Vp, v)))

    domain = (SA[0.0, 0.0], SA[1.0, π])
    prob = IntegralProblem((x, q) -> ustrip(uk, k(x, q)), domain, (fp, Vp, v))
    at = _abstol(abstol, uk; scale = 2)
    sol = solve(prob, alg; reltol, abstol = at, kwargs...)
    checkretcode && _checkretcode(_retcode(sol, reltol, at), :fχspeed)
    I = 2 * sol.u * uk

    _ismeas(typeof(cm)) || return I
    -1 < value(cm) < 1 || return I
    return I + _speedboundaryterm(L, v, cm, γ0, Δγ, reltol)
end

# First-order effect of the uncertain escape boundary on the angular integral
#     I(γ0) = ∫_{γ0}^{π} F(γ) dγ,   F(γ) = ∫ v² f(v n̂(γ, α) + V) sin γ dα,
# namely δI = -F(γ̄0) δγ0 with δγ0 = acos(c) - γ̄0 w/ no central value. 
# F is a 1D α quadrature at fixed γ, cheap next to the 2D integral. 
function _speedboundaryterm(L, v, cm, γ0, Δγ, reltol)
    (; fp, Vp, ê1, ê2, ê3) = L
    sγ, cγ = sincos(γ0 + 1e-9 * Δγ)
    F(α, q) =
        let sα_ = sin(α), cα_ = cos(α)
            n = (sγ * cα_) .* ê1 .+ (sγ * sα_) .* ê2 .+ cγ .* ê3
            v^2 * fp(v .* n + Vp) * sγ
        end
    uF = unit(F(0.0, nothing))
    sol = solve(
        IntegralProblem((α, q) -> value(ustrip(uF, F(α, q))), (0.0, Float64(π))),
        QuadGKJL();
        reltol,
    )
    return -(2 * sol.u * uF) * (acos(cm) - γ0)
end

@inline function _labframe(m::HaloModel, t)
    V = _vE(m, t)
    VEn = norm(V)
    ê3 = value.(V ./ VEn)
    ê1, ê2 = _orthobsym(ê3)
    return (;
        fp = prepare(m),
        Vp = V .- (ê2 ⋅ V) .* ê2,
        VEn,
        ê1,
        ê2,
        ê3,
        vesc2 = _vesc_eff(m)^2,
    )
end

"""
    fχspeed(m::HaloModel, t::TimeLike; uncertainty = false, analytic = true, alg = HCubatureJL(), reltol = 1e-4, abstol = 0, kwargs...) -> v -> fχspeed(m, t, v)
    fχspeed([m::HaloModel = SHM,] [t::TimeLike = T0,] v::Unitful.Velocity; uncertainty = false, analytic = true, alg = HCubatureJL(), reltol = 1e-4, abstol = 0, kwargs...)

Marginalized speed distribution, integrated over solid angle:
```math
f_{\\chi,\\text{lab}}^{\\text{speed}}(v) = \\int f_{\\chi,\\text{lab}}(\\mathbf{v})\\, d\\Omega
```

$(_ANALYTIC_DOC) $(_ABSTOL_DOC)

$(_UNCERTAINTY_DOC)

The two-argument form `fχspeed(m, t)` returns the distribution at `t` as a callable
`v -> fχspeed(m, t, v)`, with the lab frame, normalisations and the analytic-or-numerical
choice fixed. Use it to evaluate repeatedly, e.g. inside your own integrand. It does
not warn on non-convergence per evaluation (check caller's outer integral instead).

See also [`fχplanar`](@ref) for the marginalization over ``v_z``.
"""
function fχspeed(
    m::HaloModel,
    t::TimeLike,
    v::Unitful.Velocity;
    uncertainty::Bool = false,
    kwargs...,
)
    return _fχspeed_impl(_precision(m, uncertainty), t, v; kwargs...)
end
function _fχspeed_impl(
    m::HaloModel,
    t::TimeLike,
    v::Unitful.Velocity;
    alg = _CUBALG,
    reltol = 1e-4,
    abstol = 0.0,
    analytic = true,
    kwargs...,
)
    analytic && _hasanalytic(m, t) && return _fχspeedSHM(m, t, v)
    return _fχspeed(_labframe(m, t), v; alg, reltol, abstol, checkretcode = true, kwargs...)
end
fχspeed(t::TimeLike, v::Unitful.Velocity; kwargs...) = fχspeed(SHM, t, v; kwargs...)
fχspeed(m::HaloModel, v::Unitful.Velocity; kwargs...) = fχspeed(m, T0, v; kwargs...)
fχspeed(v::Unitful.Velocity; kwargs...) = fχspeed(T0, v; kwargs...)

struct _FrozenSpeed{M<:HaloModel,T,L,A,K} <: Function
    m::M
    t::T
    L::L    # lab frame from _labframe, or `nothing` on analytic codepath
    alg::A
    kw::K   # (; reltol, abstol, kwargs...)
end
(d::_FrozenSpeed{<:HaloModel,<:Any,Nothing})(v::Unitful.Velocity) = _fχspeedSHM(d.m, d.t, v)
(d::_FrozenSpeed)(v::Unitful.Velocity) = _fχspeed(d.L, v; alg = d.alg, d.kw...)
Base.show(io::IO, d::_FrozenSpeed) = print(io, "fχspeed(", d.m, ", ", d.t, ")")

function _frozenspeed(
    m::HaloModel,
    t::TimeLike;
    analytic = true,
    alg = _CUBALG,
    reltol = 1e-4,
    abstol = 0.0,
    kwargs...,
)
    L = analytic && _hasanalytic(m, t) ? nothing : _labframe(m, t)
    return _FrozenSpeed(m, t, L, alg, (; reltol, abstol, kwargs...))
end
fχspeed(m::HaloModel, t::TimeLike; uncertainty::Bool = false, kwargs...) =
    _frozenspeed(_precision(m, uncertainty), t; kwargs...)

export vmaxlab
"""
    vmaxlab([m::HaloModel = SHM,] [t::TimeLike = T0])

The largest lab-frame speed at which `m`'s distribution is nonzero at time `t`:
``v_\\text{esc} + V_E(t)``. For ``v_\\text{esc} = \\infty`` (e.g. [`SHMnoesc`](@ref)) the
escape speed is the effective cutoff ``6 v_0`` (beyond which the Maxwellian is below
``e^{-36}`` effective cutoff).

Like the pointwise distributions a model `m` with uncertainties returns a
`Measurement` (with the ephemeris' uncertainties included, see [`VE`](@ref)); use
`value(vmaxlab(m, t))` for a plain integration limit.
"""
vmaxlab(m::HaloModel, t::TimeLike) = _vesc_eff(m) + _VE(m, t)
vmaxlab(m::HaloModel) = vmaxlab(m, T0)
vmaxlab(t::TimeLike) = vmaxlab(SHM, t)
vmaxlab() = vmaxlab(SHM, T0)
