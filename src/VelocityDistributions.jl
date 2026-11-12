using NaturallyUnitful
using Measurements: value #, uncertainty
using LinearAlgebra, SpecialFunctions

# patch a removable singularity
function _erfid(x, β; delta = 1e-4)
    b = β / (1 - β)
    u = b * x^2
    value(u) < delta && return (2/sqrt(π)) * x * (1 + u/3 + u^2/10 + u^3/42)
    return erfi(x * sqrt(b)) / sqrt(b)
end

function _dawsond(x, β; delta = 1e-4)
    b = β / (1 - β)
    u = b * x^2
    value(u) < delta && return x * (1 - 2u/3 + 4u^2/15 - 8u^3/105)
    return dawson(x * sqrt(b)) / sqrt(b)
end

"""
    PreparedHalo

Internal: callable cache of the `HaloModel` normalizations for hoisting out of integrands.
`p(v)` evaluates ``f_\\chi``. Users can get the same effect through the frozen forms
`fχ(m)`, `fR(m)`, ... (see [`fχ`](@ref)) without depending on this type.
"""
struct PreparedHalo{TAr,TAs,TWr,TWs,TSv,TSr,TSθ,TV2,Tη}
    Ar::TAr         # fR prefactor
    As::TAs         # fS prefactor
    wR::TWr         # (1-η)·Ar, the fR weight in fχ
    wS::TWs         # η·As,     the fS weight in fχ
    invσv2::TSv     # 1/(2σv²)
    invσr2::TSr     # 1/(2σr²)
    invσθ2::TSθ     # 1/(2σθ²) = 1/(2σϕ²)
    vₑ2::TV2        # escape speed squared
    η::Tη
end

"""
    prepare(m::HaloModel) -> PreparedHalo

Internal: hoist the normalizations out of hot loops. Reuse when evaluating repeatedly.
"""
function prepare(m::HaloModel)
    v₀, vₑ, β, η = m.v₀, m.vₑ, m.β, m.η

    # Maxwellian component.
    σv2 = v₀^2 / 2
    zr = vₑ / v₀
    NResc = isfinite(vₑ) ? erf(zr) - 2zr * exp(-zr^2) / sqrt(π) : one(zr)
    Ar = 1 / ((2π)^(3 / 2) * σv2 * sqrt(σv2) * NResc)

    # Sausage component.
    σr2 = 3v₀^2 / (2(3 - 2β))
    σθ2 = 3v₀^2 * (1 - β) / (2(3 - 2β))
    zs = vₑ / sqrt(2σr2)

    NSesc = if iszero(β)
        NResc
    elseif isfinite(vₑ)
        #erf(zs) - exp(-vₑ^2 / (2σθ2)) * _erfid(zs, β)
        erf(zs) - (2 / sqrt(π)) * exp(-zs^2) * _dawsond(zs, β)
    else
        one(zs)
    end
    As = 1 / ((2π)^(3 / 2) * sqrt(σr2) * σθ2 * NSesc)
    return PreparedHalo(
        Ar,
        As,
        (1 - η) * Ar,
        η * As,
        inv(2σv2),
        inv(2σr2),
        inv(2σθ2),
        vₑ^2,
        η,
    )
end

@inline _zeroexp(A, x) = zero(A) * zero(Base.promote_op(exp, typeof(x)))

export fR
@inline function _fR(p::PreparedHalo, v2)
    x = -v2 * p.invσv2
    v2 > p.vₑ2 && return _zeroexp(p.Ar, x)
    return p.Ar * exp(x)
end
"""
    fR([m::HaloModel = SHM,] v::AbstractVector{<:Unitful.Velocity})
    fR([m::HaloModel = SHM,] v::Unitful.Velocity)
    fR(m::HaloModel) -> v -> fR(m, v)

Galactic-frame Maxwellian component of ``f_\\chi``.

``f_R(\\mathbf{v}) = \\frac{1}{(2\\pi\\sigma_v^2)^{3/2}N_{R,\\text{esc}}}
    \\exp\\left(-\\frac{|\\mathbf{v}|^2}{2\\sigma_v^2}\\right)\\Theta(v_\\text{esc} - |\\mathbf{v}|)``

``N_{R,\\text{esc}} = \\text{erf}\\left(\\frac{v_\\text{esc}}{\\sqrt{2}\\sigma_v}\\right) 
    - \\sqrt{\\frac{2}{\\pi}}\\frac{v_\\text{esc}}{\\sigma_v}\\exp\\left(-\\frac{v_\\text{esc}^2}{2\\sigma_v^2}\\right)``

with ``\\sigma_v = v_0/\\sqrt{2}``. Isotropic, so `v` may be given as a scalar speed.

For `fR(m)` see "Repeated evaluation" under [`fχ`](@ref).
"""
fR(m::HaloModel, v::Unitful.Velocity) = _fR(prepare(m), v^2)
fR(m::HaloModel, v::AbstractVector{<:Unitful.Velocity}) = _fR(prepare(m), v ⋅ v)
fR(v::Unitful.Velocity) = fR(SHM, v)
fR(v::AbstractVector{<:Unitful.Velocity}) = fR(SHM, v)
export fRlab
"""
    fRlab([m::HaloModel = SHM,] [t::TimeLike = T0,] v::AbstractVector{<:Unitful.Velocity})
    fRlab(m::HaloModel) -> ([t,] v) -> fRlab(m, t, v)
    fRlab(m::HaloModel, t::TimeLike) -> v -> fRlab(m, t, v)

Lab-frame [`fR`](@ref): `fR(m, v + vE(m, t))`. For `fRlab(m)` see "Repeated evaluation"
under [`fχ`](@ref); `fRlab(m, t)` also fixes the time, and so the Earth's velocity.
"""
fRlab(m::HaloModel, t::TimeLike, v::AbstractVector{<:Unitful.Velocity}) =
    fR(m, v + _vE(m, t))
fRlab(t::TimeLike, v::AbstractVector{<:Unitful.Velocity}) = fRlab(SHM, t, v)
fRlab(m::HaloModel, v::AbstractVector{<:Unitful.Velocity}) = fRlab(m, T0, v)
fRlab(v::AbstractVector{<:Unitful.Velocity}) = fRlab(SHM, v)

export fS
@inline function _fS(p::PreparedHalo, v::AbstractVector{<:Unitful.Velocity})
    # Spherical polar at Earth -> cartesian LSR (arXiv:1810.11468 footnote 1):
    # v[1] ≡ vr, v[2] ≡ vθ, v[3] ≡ vϕ.
    v2 = v ⋅ v
    x = -v[1]^2 * p.invσr2 - (v[2]^2 + v[3]^2) * p.invσθ2
    v2 > p.vₑ2 && return _zeroexp(p.As, x)
    return p.As * exp(x)
end
"""
    fS([m::HaloModel = SHM,] v::AbstractVector{<:Unitful.Velocity})
    fS(m::HaloModel) -> v -> fS(m, v)

Galactic-frame anisotropic (Gaia sausage) component of ``f_\\chi``.

``f_S(\\mathbf{v}) = \\frac{1}{(2\\pi)^{3/2}\\sigma_r\\sigma_\\theta^2 N_{S,\\text{esc}}}
    \\exp\\left(-\\frac{v_r^2}{2\\sigma_r^2} - \\frac{v_\\theta^2}{2\\sigma_\\theta^2} - \\frac{v_\\phi^2}{2\\sigma_\\phi^2}\\right)
    \\Theta(v_\\text{esc} - |\\mathbf{v}|)``

``N_{S,\\text{esc}} = \\text{erf}\\left(\\frac{v_\\text{esc}}{\\sqrt{2}\\sigma_r}\\right) - \\sqrt{\\frac{1-\\beta}{\\beta}}
\\exp\\left(-\\frac{v_\\text{esc}^2}{2\\sigma_\\theta^2}\\right)\\text{erfi}\\left(\\frac{v_\\text{esc}}{\\sqrt{2}\\sigma_r}\\sqrt{\\frac{\\beta}{1-\\beta}}\\right)``

Anisotropic, so `v` must be a vector; `v[1]` is the radial (sausage) axis ``v_r``.

For `fS(m)` see "Repeated evaluation" under [`fχ`](@ref).
"""
fS(m::HaloModel, v::AbstractVector{<:Unitful.Velocity}) = _fS(prepare(m), v)
fS(v::AbstractVector{<:Unitful.Velocity}) = fS(SHM, v)
export fχ
@inline function (p::PreparedHalo)(v::AbstractVector{<:Unitful.Velocity})
    v2 = v ⋅ v
    x = -v2 * p.invσv2

    v2 > p.vₑ2 && return _zeroexp(p.wR, x)
    a = p.wR * exp(x)

    return iszero(p.η) ? a :
           a + p.wS * exp(-v[1]^2 * p.invσr2 - (v[2]^2 + v[3]^2) * p.invσθ2)
end
"""
    fχ([m::HaloModel = SHM,] v::AbstractVector{<:Unitful.Velocity})
    fχ(m::HaloModel) -> v -> fχ(m, v)

Galactic-frame SHM++ velocity distribution: a mixture of the round [`fR`](@ref) and
anisotropic [`fS`](@ref) components.

``f_\\chi(\\mathbf{v}) = (1 - \\eta) f_R(\\mathbf{v}) + \\eta f_S(\\mathbf{v})``

## Repeated evaluation

Each `fχ(m, v)` call recomputes the model's normalizations (`erf`, `erfi` and `exp` of the
parameters), which for a `Measurement` model costs several times the evaluation itself.

The one-argument form `f = fχ(m)` computes them once and returns a callable `f(v)`, to use
in hot loops, closures or integrands. Broadcasting over velocities with a single model, as
in `fχ.(m, vs)` or `fχ.(vs)`, does this automatically.

The same applies to [`fR`](@ref), [`fS`](@ref), [`fRlab`](@ref) and [`fχlab`](@ref).
"""
fχ(m::HaloModel, v::AbstractVector{<:Unitful.Velocity}) = prepare(m)(v)
fχ(v::AbstractVector{<:Unitful.Velocity}) = fχ(SHM, v)
export fχlab
"""
    fχlab([m::HaloModel = SHM,] [t::TimeLike = T0,] v::AbstractVector{<:Unitful.Velocity})
    fχlab(m::HaloModel) -> ([t,] v) -> fχlab(m, t, v)
    fχlab(m::HaloModel, t::TimeLike) -> v -> fχlab(m, t, v)

Lab-frame [`fχ`](@ref): `fχ(m, v + vE(m, t))`. For `fχlab(m)` see "Repeated evaluation"
under [`fχ`](@ref); `fχlab(m, t)` also fixes the time, and so the Earth's velocity.
"""
fχlab(m::HaloModel, t::TimeLike, v::AbstractVector{<:Unitful.Velocity}) =
    fχ(m, v + _vE(m, t))
fχlab(t::TimeLike, v::AbstractVector{<:Unitful.Velocity}) = fχlab(SHM, t, v)
fχlab(m::HaloModel, v::AbstractVector{<:Unitful.Velocity}) = fχlab(m, T0, v)
fχlab(v::AbstractVector{<:Unitful.Velocity}) = fχlab(SHM, v)

# Frozen distributions ==============================================================
struct _Frozen{F,M<:HaloModel,P<:PreparedHalo} <: Function
    f::F
    m::M
    p::P
end
_Frozen(f, m::HaloModel) = _Frozen(f, m, prepare(m))
Base.show(io::IO, d::_Frozen) = print(io, d.f, "(", d.m, ")")

const _VelVec = AbstractVector{<:Unitful.Velocity}

fR(m::HaloModel) = _Frozen(fR, m)
(d::_Frozen{typeof(fR)})(v::Unitful.Velocity) = _fR(d.p, v^2)
(d::_Frozen{typeof(fR)})(v::_VelVec) = _fR(d.p, v ⋅ v)

fS(m::HaloModel) = _Frozen(fS, m)
(d::_Frozen{typeof(fS)})(v::_VelVec) = _fS(d.p, v)

fχ(m::HaloModel) = _Frozen(fχ, m)
(d::_Frozen{typeof(fχ)})(v::_VelVec) = d.p(v)

fRlab(m::HaloModel) = _Frozen(fRlab, m)
(d::_Frozen{typeof(fRlab)})(t::TimeLike, v::_VelVec) =
    (w = v + _vE(d.m, t); _fR(d.p, w ⋅ w))
(d::_Frozen{typeof(fRlab)})(v::_VelVec) = d(T0, v)

fχlab(m::HaloModel) = _Frozen(fχlab, m)
(d::_Frozen{typeof(fχlab)})(t::TimeLike, v::_VelVec) = d.p(v + _vE(d.m, t))
(d::_Frozen{typeof(fχlab)})(v::_VelVec) = d(T0, v)

struct _FrozenLab{F,M<:HaloModel,P<:PreparedHalo,T,V} <: Function
    f::F
    m::M
    p::P
    t::T
    V::V    # _vE(m, t)
end
_FrozenLab(f, m::HaloModel, t) = _FrozenLab(f, m, prepare(m), t, _vE(m, t))
Base.show(io::IO, d::_FrozenLab) = print(io, d.f, "(", d.m, ", ", d.t, ")")

fRlab(m::HaloModel, t::TimeLike) = _FrozenLab(fRlab, m, t)
(d::_FrozenLab{typeof(fRlab)})(v::_VelVec) = (w = v + d.V; _fR(d.p, w ⋅ w))

fχlab(m::HaloModel, t::TimeLike) = _FrozenLab(fχlab, m, t)
(d::_FrozenLab{typeof(fχlab)})(v::_VelVec) = d.p(v + d.V)

for f in (:fR, :fS, :fχ, :fRlab, :fχlab)
    @eval Base.Broadcast.broadcasted(::typeof($f), m::HaloModel, v, args...) =
        Base.Broadcast.broadcasted($f(m), v, args...)
end
for f in (:fRlab, :fχlab)
    @eval Base.Broadcast.broadcasted(::typeof($f), m::HaloModel, t::TimeLike, v) =
        Base.Broadcast.broadcasted($f(m, t), v)
    @eval Base.Broadcast.broadcasted(::typeof($f), m::HaloModel, t::TimeLike) = $f(m, t)
end

const _VelArray = AbstractArray{<:Union{Unitful.Velocity,_VelVec}}
const _VelVecArray = AbstractArray{<:_VelVec}
Base.Broadcast.broadcasted(::typeof(fR), v::_VelArray) =
    Base.Broadcast.broadcasted(fR(SHM), v)
for f in (:fS, :fχ, :fRlab, :fχlab)
    @eval Base.Broadcast.broadcasted(::typeof($f), v::_VelVecArray) =
        Base.Broadcast.broadcasted($f(SHM), v)
end
for f in (:fRlab, :fχlab)
    @eval Base.Broadcast.broadcasted(
        ::typeof($f),
        t::Union{TimeLike,AbstractArray{<:TimeLike}},
        v,
    ) = Base.Broadcast.broadcasted($f(SHM), t, v)
end
