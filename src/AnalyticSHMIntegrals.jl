using SpecialFunctions

# Closed form for an SHM-like isotropic (i.e. η = 0 -- no sausage) model 
# based on exposition in DelNobile:2021wmp (see Documenter references)

# Dimensionless speeds throughout, in units of v₀: z = v/v₀, zE = V_E/v₀, zesc = vₑ/v₀.

# Normalisation of the truncated Maxwellian (NResc in prepare).
function _NescSHM(m::HaloModel)
    z = NoUnits(m.vₑ / m.v₀)
    isfinite(z) || return one(z)
    return erf(z) - (2 / sqrt(π)) * z * exp(-z^2)
end

# The lab-frame speed distribution (fχspeed).
function _fχspeedSHM(m::HaloModel, t::TimeLike, v::Unitful.Velocity)
    v0 = m.v₀
    zesc = NoUnits(m.vₑ / v0)
    z = NoUnits(v / v0)
    zE = NoUnits(_VE(m, t) / v0)
    xp = min(z + zE, zesc)
    xm = min(z - zE, zesc)
    return z * (exp(-xm^2) - exp(-xp^2)) / (sqrt(π) * v0 * zE * _NescSHM(m))
end

# The mean reciprocal speed g(vmin).
function _gSHM(m::HaloModel, t::TimeLike, vmin::Unitful.Velocity)
    v0 = m.v₀
    zmin = NoUnits(vmin / v0)
    zesc = NoUnits(m.vₑ / v0)
    zE = NoUnits(_VE(m, t) / v0)
    zps = min(zmin + zE, zesc)
    zms = min(zmin - zE, zesc)
    esc = isfinite(zesc) ? (2 / sqrt(π)) * exp(-zesc^2) * (zps - zms) : zero(zps)
    return (erf(zps) - erf(zms) - esc) / (2 * _NescSHM(m) * v0 * zE)
end

# The mean speed h(vmin).
function _hSHM(m::HaloModel, t::TimeLike, vmin::Unitful.Velocity)
    v0 = m.v₀
    zmin = NoUnits(vmin / v0)
    zesc = NoUnits(m.vₑ / v0)
    zE = NoUnits(_VE(m, t) / v0)

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
    ) * v0 / (sqrt(π) * _NescSHM(m) * zE)
end
