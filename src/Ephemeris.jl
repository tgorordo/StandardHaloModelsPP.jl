using LinearAlgebra, NaturallyUnitful, UnitfulAstro, Dates, StaticArrays

# Loose naming conventions in this file:
#   lowercase (vLSR, vSUN, vES, vE) - velocity *vectors*
#   uppercase (VLSR, VSUN, VES, VE) - the corresponding *speeds* (norms) or scalar/time constants
#   leading underscore (_vE, _VE)   - model-aware: dispatches on `isexact(m)` to use the
#                                     central-value ephemeris for models without uncertainties



## Solar motion in the Galactic and LSR frames =====================================

export d0
"""
Distance of the sun from the Galactic Center (Sgr. A*). 

Value: ``d_0 =`` $(d0).

[R. Abuter and others (GRAVITY collaboration). Astron. Astrophys. 647, A59 (2021)](@cite Abuter:2021yys):
``8275 \\pm 9_\\text{stat} \\pm 33_\\text{sys}`` pc, combined in quadrature. This is the
distance underlying [`VLSR`](@ref) in [Baxter:2021pqo](@cite); it supersedes the
``8.122 \\pm 0.031`` kpc of [GRAVITY:2018ofz](@cite) used in the defining SHM++ paper.
"""
const d0 = measurement(8.275, 0.034)u"kpc"

export θ0
"""
Azimuthal angle of the solar system in Galactocentric Spherical Coordinates. 

Value: ``\\theta_0 =`` $(θ0)

(Zero by definition.)
"""
const θ0 = 0.0

export z0
"""
Normal distance of the sun above the Galactic Plane. 

Value: ``z_0 =`` $(z0)

[J. Bland-Hawthorn and O. Gerhard, Ann. Rev. Astron. Astrophys. 54, 529 (2016)](@cite Bland-Hawthorn:2016lwg),
adopting the SDSS estimate of [M. Jurić and others (SDSS collaboration), Astrophys. J. 673, 864 (2008)](@cite SDSS:2005kst).
"""
const z0 = measurement(25.0, 5.0)u"pc"

export η0
"""
Altitude (angle) of the sun above the Galactic Plane. 

Value: ``\\eta_0 = \\sin^{-1}{z_0 / d_0} =`` $(η0) (radians)
"""
const η0 = asin(z0 / d0)

export ϕ0
"""
Polar angle of the sun from the galactic north zenith.

Value: ``\\phi_0 = \\frac{\\pi}{2} - \\eta_0 =`` $(ϕ0) (radians)
"""
const ϕ0 = π/2 - η0

export R0
"""
Radial (in-galactic-plane/cylindrical) distance of the sun from the Galactic Center (Sgr. A*). 

Value: ``R_0 = d_0 \\cos{\\eta_0} =`` $(R0) 
"""
const R0 = d0 * cos(η0)

export VLSR
"""
Circular speed of the Local Standard of Rest (LSR).

Value: $(VLSR)

[D. Baxter and others, Eur. Phys. J. C 81, 907 (2021); Sec. 3.1.5](@cite Baxter:2021pqo).

!!! danger
    `VLSR` is part of the ephemeris, independent of the circular speed `v₀` of any
    [`HaloModel`](@ref). It fixes the Galactic → lab coordinate transformation, while `v₀`
    sets the dispersion of the DM velocity distribution. They coincide for [`SHM`](@ref)
    but [`SHMpp`](@ref) keeps `v₀ = 233` km/s, and [`SHMlegacy`](@ref) keeps the
    historical `v₀ = 220` km/s; both are still boosted with this `VLSR` for lab-frame
    distributions.
"""
const VLSR = measurement(238.0, 1.5)u"km/s"

export vLSR
"""
Circular velocity (vector) of the Local Standard of Rest (LSR) in galactic rectangular coordinates.
 
Value: ``v_\\text{LSR} = (v_x, v_y, v_z) =`` ($(vLSR[1]), $(vLSR[2]), $(vLSR[3]))
"""
const vLSR = SA[0.0u"km/s", VLSR, 0.0u"km/s"]

export vPEC
"""
Peculiar velocity of the sun in the Local Standard of Rest (LSR) in galactic rectangular coordinates.

Value: ``v_\\text{PEC} = (U, V, W) =`` ($(vPEC[1]), $(vPEC[2]), $(vPEC[3]))

[R. Schönrich, J. Binney, W. Dehnen, Mon. Not. Roy. Astron. Soc. 403, 1829 (2010)](@cite Sch_nrich_2010),
with the (statistical and systematic) uncertainties as adopted by [the SHM++](@cite PhysRevD.99.023012).
"""
const vPEC = SA[
    measurement(11.1, 1.5)u"km/s", # UO, 
    measurement(12.2, 2.0)u"km/s", # VO, 
    measurement(7.3, 1)u"km/s",    # WO
]

export vSUN
"""
Velocity of the sun in the Galactic frame.

Value: ``v_\\odot = (v_x, v_y, v_z) = (v_r, v_\\theta, v_\\phi) =`` ($(vSUN[1]), $(vSUN[2]), $(vSUN[3]))
"""
const vSUN = vLSR + vPEC

export VSUN
"""
Speed of the sun in the galactic frame. A good ballpark typical DM velocity ignoring Earth ephemeris.
"""
const VSUN = norm(vSUN)

## Ephemeris of the Earth about the Sun ============================================

"Angular frequency of the Earth's Solar orbit (rad/yr). Value: $(ωES)"
const ωES = 2π / 1u"yr"

export J2000
"The J2000.0 epoch, 2000-01-01 12:00 (TT, taken as UTC at this precision). Value: $(J2000)"
const J2000 = DateTime(2000, 1, 1, 12)

"March 20th, 2020 Equinox `DateTime`. Value: $(Teq)"
const Teq = DateTime(2020, 3, 20, 3, 50)

"Julian centuries from [`J2000`](@ref) to `t`: the epoch of date ``T`` of [McCabe_2014; Eq. (15)](@cite)."
_jcenturies(t::DateTime) = Dates.value(t - J2000) / (86_400_000 * 36525)


"""
Direction of the Earth's velocity at the March equinox [`Teq`](@ref) (equivalently the
Sun → Earth direction at the June solstice), in galactic coordinates.

``\\hat\\epsilon_1 = \\hat{y}(T)`` of [McCabe_2014; Eq. (14)](@cite), at the epoch of `Teq`;
cf. ``\\hat\\epsilon_1`` of [Lee:2013xxa; Eq. (B.7)](@cite).

Value: $(ϵ1)
"""
const ϵ1 =
    SA[0.993821, 0.110992, 0.000352] + SA[0.001316, -0.011851, 0.021267] * _jcenturies(Teq)

"""
Earth → Sun direction at the March equinox [`Teq`](@ref) (equivalently the direction of
the Earth's velocity at the June solstice), in galactic coordinates.

``\\hat\\epsilon_2 = -\\hat{x}(T)`` of [McCabe_2014; Eq. (13)](@cite), at the epoch of `Teq`;
cf. ``\\hat\\epsilon_2`` of [Lee:2013xxa; Eq. (B.7)](@cite).

Value: $(ϵ2)
"""
const ϵ2 = -(
    SA[0.054876, -0.494109, 0.867666] +
    SA[-0.024232, -0.002689, 1.546e-6] * _jcenturies(Teq)
)

_hat(v) = v / norm(v) # unit vector

export T0
"June 1st reference vS-ve alignment. Value: $(T0)"
const T0 = Teq + round(Day, atan(_hat(vSUN) ⋅ _hat(ϵ2), _hat(vSUN) ⋅ _hat(ϵ1)) / ωES)

"Radius of Earth's orbit. Value: $(RES)"
const RES = 149.598e6u"km"
# neglects eccentricity

"""
Tangential speed of the Earth in the Solar frame, ``a\\omega``, for a circular orbit.
cf. ``\\langle u_E \\rangle = 29.79`` km/s of [McCabe_2014; Eq. (12)](@cite).

Value: $(VES)
"""
const VES = uconvert(u"km/s", RES * ωES)

export vES
"""
Velocity of the Earth in the Solar frame as a function of time, for a circular orbit:

``\\mathbf{v}_{ES}(t) = V_{ES} \\left[\\hat\\epsilon_1 \\cos\\omega(t - t_\\text{eq}) + \\hat\\epsilon_2 \\sin\\omega(t - t_\\text{eq})\\right]``

the zeroth order in the eccentricity of [Lee:2013xxa; Eq. (B.7)](@cite) and
[McCabe_2014; Eq. (21)](@cite).
"""
vES(t::DateTime) = VES * (cos(ωES * (t - Teq)) .* ϵ1 + sin(ωES * (t - Teq)) .* ϵ2)
vES(t::Unitful.Time) = vES(T0 + round(Millisecond, t))

export vE
"Velocity of the Earth in the galactic frame as a function of time."
vE(t) = vES(t) + vSUN

export VE
"Speed of the Earth in the galactic frame as a function of time."
VE(t::DateTime) = norm(vE(t))
VE(t::Unitful.Time) = VE(T0 + round(Millisecond, t))

export DayOfYear
"""
    DayOfYear(d::Real[, year::Integer = Dates.year(T0)])

Day `d` of `year` as a time argument: `d = 1` is January 1st, 00:00 UTC, and fractional
days are allowed, so `DayOfYear(152.5)` is noon on June 1st 2020.

A plain number is deliberately *not* accepted as a time; instead use a `Unitful` time such as `100u"d"`
to mean the time elapsed since [`T0`](@ref).
"""
struct DayOfYear
    t::DateTime
end
DayOfYear(d::Real, year::Integer = Dates.year(T0)) =
    DayOfYear(DateTime(year) + Millisecond(round(Int64, (d - 1) * 86_400_000)))
Dates.DateTime(d::DayOfYear) = d.t
Base.broadcastable(d::DayOfYear) = Ref(d)
function Base.show(io::IO, d::DayOfYear)
    y = Dates.year(d.t)
    print(io, "DayOfYear(", Dates.value(d.t - DateTime(y)) / 86_400_000 + 1, ", ", y, ")")
end

vES(t::Union{Date,DayOfYear}) = vES(DateTime(t))
VE(t::Union{Date,DayOfYear}) = VE(DateTime(t))

export TimeLike
"""
    TimeLike = Union{DateTime, Date, DayOfYear, Unitful.Time}
"""
const TimeLike = Union{DateTime,Date,DayOfYear,Unitful.Time}

"""
Annual mean of [`VE`](@ref).

Value: $(VEavg)
"""
const VEavg = let N = 360
    # trapezoid rule
    sum(VE((i / N) * 1u"yr") for i = 0:(N-1)) / N
end

"""
Time after [`T0`](@ref) at which [`VE`](@ref) rises through its annual mean
[`VEavg`](@ref), so `VE(Tavg) ≈ VEavg`. About March 2nd. 

`VE` also falls through its mean roughly half a year later, around September 1st.

Value: $(Tavg)
"""
const Tavg =
    let f = d -> value(ustrip(u"km/s", VE(d * u"d"))) - value(ustrip(u"km/s", VEavg))
        # Bisection on central values between `VE` min (~`T0 + yr/2`) and max (~`T0 + yr`).

        a, b = 365.25 / 2, 365.25
        for _ = 1:40
            c = (a + b) / 2
            f(c) < 0 ? (a = c) : (b = c)
        end
        ((a + b) / 2)u"d"
    end

# Models with no Measurement parameters get the central-value solar velocity; 
# i.e. exact models do not inherit vSUN uncertainties.
_ismeas(::Type{<:Measurements.Measurement}) = true
_ismeas(::Type{<:Quantity{T}}) where {T} = _ismeas(T)
_ismeas(::Type) = false

isexact(::Type{<:HaloModel{Tρ,Tv,Tβ}}) where {Tρ,Tv,Tβ} =
    !(_ismeas(Tρ) || _ismeas(Tv) || _ismeas(Tβ))
isexact(m::HaloModel) = isexact(typeof(m))

"Central-value solar velocity, for exact-model evaluation. See `isexact`."
const vSUNc = value.(vSUN)

"""
Ephemeris with the uncertainty-propagating `vE` for a model that carries uncertainties, and central-values otherwise.
"""
_vE(m::HaloModel, t) = _vE(Val(isexact(m)), t)
_vE(::Val{true}, t) = vES(t) + vSUNc
_vE(::Val{false}, t) = vES(t) + vSUN

_VE(m::HaloModel, t) = norm(_vE(m, t))

"""
    vE(m::HaloModel, t::TimeLike)
    VE(m::HaloModel, t::TimeLike)

The Earth's galactic velocity and speed at `t`.

If `m` does not contain uncertainties then ephemeris uncertainties are
not added to the result (for performance reasons), 
otherwise uncertainties will be included.

Unqualified `vE(t)`/`VE(t)` always carry uncertainties.
"""
vE(m::HaloModel, t::TimeLike) = _vE(m, t)
VE(m::HaloModel, t::TimeLike) = _VE(m, t)

# Sausage anisotropy axis: v[1] ≡ vr (arXiv:1810.11468 footnote 1).
const _AXIS = SA[1.0, 0.0, 0.0]

# Orthonormal frame with ê3 as given and ê1 in the (ê3, x̂) plane, so ê2 ⋅ x̂ = 0.
# Makes the angular integrand even in α; see _labframe.
function _orthobsym(ê3)
    x⊥ = _AXIS - ê3[1] .* ê3
    n = norm(x⊥)

    n < 1e-8 && return SA[0.0, 1.0, 0.0], ê3 × SA[0.0, 1.0, 0.0]
    ê1 = x⊥ ./ n
    return ê1, ê3 × ê1
end
