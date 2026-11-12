using SpecialFunctions

# Closed form for an SHM-like isotropic (i.e. η = 0 -- no sausage) model
# based on exposition in DelNobile:2021wmp (see Documenter references)

# Dimensionless speeds throughout, in units of v₀: z = v/v₀, zE = V_E/v₀, zesc = vₑ/v₀.

# Normalization of the truncated Maxwellian (NResc in prepare).
function _NescSHM(m::HaloModel)
    z = NoUnits(m.vₑ / m.v₀)
    isfinite(z) || return one(z)
    return erf(z) - (2 / sqrt(π)) * z * exp(-z^2)
end

# frozen constants
struct _SHMConsts{V,Z,E,N}
    v0::V   # v₀
    zesc::Z # vₑ/v₀
    zE::E   # V_E(t)/v₀
    Nesc::N # _NescSHM(m)
end
_SHMConsts(m::HaloModel, t) =
    _SHMConsts(m.v₀, NoUnits(m.vₑ / m.v₀), NoUnits(_VE(m, t) / m.v₀), _NescSHM(m))

function _analyticconsts(m::HaloModel, t, analytic::Bool)
    (analytic && iszero(m.η)) || return nothing
    c = _SHMConsts(m, t)
    return value(c.zE) < value(c.zesc) ? c : nothing
end

# The lab-frame speed distribution (fχspeed).
_fχspeedSHM(m::HaloModel, t::TimeLike, v::Unitful.Velocity) =
    _fχspeedSHM(_SHMConsts(m, t), v)
function _fχspeedSHM(c::_SHMConsts, v::Unitful.Velocity)
    (; v0, zesc, zE, Nesc) = c
    z = NoUnits(v / v0)
    xp = min(z + zE, zesc)
    xm = min(z - zE, zesc)
    return z * (exp(-xm^2) - exp(-xp^2)) / (sqrt(π) * v0 * zE * Nesc)
end

# The mean reciprocal speed g(vmin).
_gSHM(m::HaloModel, t::TimeLike, vmin::Unitful.Velocity) = _gSHM(_SHMConsts(m, t), vmin)
function _gSHM(c::_SHMConsts, vmin::Unitful.Velocity)
    (; v0, zesc, zE, Nesc) = c
    zmin = NoUnits(vmin / v0)
    zps = min(zmin + zE, zesc)
    zms = min(zmin - zE, zesc)
    esc = isfinite(zesc) ? (2 / sqrt(π)) * exp(-zesc^2) * (zps - zms) : zero(zps)
    return (erf(zps) - erf(zms) - esc) / (2 * Nesc * v0 * zE)
end

# The mean speed h(vmin).
_hSHM(m::HaloModel, t::TimeLike, vmin::Unitful.Velocity) = _hSHM(_SHMConsts(m, t), vmin)
function _hSHM(c::_SHMConsts, vmin::Unitful.Velocity)
    (; v0, zesc, zE, Nesc) = c
    zmin = NoUnits(vmin / v0)

    zp = min(zmin, zesc + zE)
    zm = min(zmin, zesc - zE)

    zps = min(zmin + zE, zesc)
    zms = min(zmin - zE, zesc)

    esc =
        isfinite(zesc) ?
        exp(-zesc^2) * ((zm^3 - zp^3) / 3 + 2zE * (1 + zE^2 / 3 + zesc^2)) : zero(zps)
    return (
        (sqrt(π) / 4) * (1 + 2zE^2) * (erf(zps) - erf(zms)) -
        ((zps - 2zE) * exp(-zps^2) - (zms + 2zE) * exp(-zms^2)) / 2 - esc
    ) * v0 / (sqrt(π) * Nesc * zE)
end
