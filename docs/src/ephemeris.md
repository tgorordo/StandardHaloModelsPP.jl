# Galactic & Solar Ephemerides - the Earth-Frame Distributions

```@meta
CurrentModule = StandardHaloModelsPP
```

The Earth-frame distributions follow from the
[galactic-frame ones](@ref "Introduction - Galactic-Frame Velocity Distributions") by a boost
with the Earth's galactic velocity as a composition of the Earth-Sun and Sun-(Sag. A⋆) orbits.

[`StandardHaloModelsPP.jl`](https://github.com/tgorordo/StandardHaloModelsPP.jl) provides
a simple model of that motion in which the Sun's ephemeris follows the
conventions of [Baxter:2021pqo](@cite), with its peculiar velocity from
[Sch_nrich_2010](@cite) and its position from [Abuter:2021yys](@cite), and the Earth's orbit
follows the conventions of [McCabe_2014, Lee:2013xxa](@cite).

## Ephemerides

The Sun's galactic velocity ``\vec{v}_\odot = \vec{v}_\text{LSR} + \vec{v}_\text{PEC}`` is that
of the Local Standard of Rest (LSR), the circular orbit about the Galactic Center at the Sun's
position, plus the Sun's peculiar velocity relative to it. Schematically (not to scale):
```@raw html
<iframe src="../figs/gal.html" width="800" height="600"></iframe>
```
```@raw latex
\begin{figure}[htbp]
\centering
\includegraphics[width=0.8\textwidth]{figs/gal.pdf}
\caption{The Sun's galactocentric position, the Galactic-Rectangular axes (along which the halo velocity components are taken), and the LSR and peculiar velocities of the Sun.}
\end{figure}
```
([Fig. `gal.asy` source](figs/gal.asy))

Galactic-Rectangular coordinates (green, centered in the plane below the Sun) are defined
with ``\hat{x}`` pointing from the Sun to Sgr A*, ``\hat{y}`` along the disk's rotation
(clockwise seen from the North Galactic Pole), and ``\hat{z}`` to the North Galactic Pole
[McCabe_2014](@cite) [BinneyMerrifield:1998; Sec. 10.3.1](@cite).

The DM halo velocity distribution is usually expressed in Galactocentric-Spherical 
coordinates [PhysRevD.99.023012](@cite). Its components at the sun are taken
along the Galactic-Rectangular axes [^1], and at the location of the solar system 
``(v_r, v_\theta, v_\phi) \equiv (v_x, v_y, v_z)`` (labeled in the figure). 

The lab-frame distribution is related by a Galilean boost in that frame
[McCabe_2014](@cite): a DM velocity ``\vec{v}`` in the lab is ``\vec{v} + \vec{v}_E`` in the
galaxy; ``f_\text{lab}(\vec{v}, t) = f(\vec{v} + \vec{v}_E(t))``.

[^1]: The angle of the sun above the galactic plane is small and can be neglected -- [`η0`](@ref) is *extremely* exaggerated in the figure.

Position of the Sun:
```@docs; canonical=false
d0
z0
θ0
η0
ϕ0
R0
```

Velocity of the Sun:
```@docs; canonical=false
VLSR
vLSR
vPEC
vSUN
VSUN
```

In most halo models the DM (isothermal) distribution is assumed to be stationary with respect to
the galaxy, so it is the motion of the solar system through the galaxy
that sets ``v_0``, the DM velocity dispersion [PhysRevD.99.023012](@cite) -- 
i.e. the LSR speed is essentially that parameter of a halo model. Each [`HaloModel`](@ref) therefore
carries its own LSR speed `m.VLSR`, by default equivalent to its `v₀`, and ephemeris functions are
model-aware; i.e. [`vE(m, t)`](@ref vE), the lab-frame distributions and the integrals all 
use it for their boost to the lab. 
The constants here are a reference ephemeris, with [`VLSR`](@ref) from [Baxter:2021pqo](@cite)
(that of [`SHM`](@ref)), used whenever no model is given, e.g. by `VE(t)`.

The Earth's orbit is taken as circular, in the ecliptic, i.e. at zeroth order in its
eccentricity [Lee:2013xxa, McCabe_2014](@cite):
```math
\vec{v}_{ES}(t) = V_{ES}\left[\hat\epsilon_1 \cos\omega(t - t_\text{eq}) + \hat\epsilon_2 \sin\omega(t - t_\text{eq})\right],
```
- ``\omega = 2\pi/\text{yr}`` is the orbital angular frequency.
- ``V_{ES} = a\omega \approx 29.79`` km/s, with ``a`` = 1 au.
- ``t_\text{eq}`` is the March equinox.
- ``\hat\epsilon_1`` is the direction of the Earth's velocity, and ``\hat\epsilon_2`` the
  Earth → Sun direction, at ``t_\text{eq}``.

Both basis vectors are expressed in Galactic-Rectangular coordinates, to first order in the
slow precession of the equinoxes (measured from [`J2000`](@ref)), and evaluated at the epoch
of ``t_\text{eq}`` [McCabe_2014; Eqs. (13-15)](@cite).

!!! todo "TODO: Include an elliptical model"
    [McCabe_2014](@cite) gives ``\vec{v}_{ES}`` to first order in the eccentricity, 
    but Analyses relying on e.g. an annual modulation might need a more-accurate elliptic 
    model of the Earth's orbit.

The Earth's galactic velocity is then
```math
\vec{v}_E(t) = \vec{v}_\odot + \vec{v}_{ES}(t),
```
whose speed peaks at ``T_0`` (about June 1st), when ``\vec{v}_{ES}`` is most closely aligned
with ``\vec{v}_\odot``. Times may be given as a `DateTime`, a `Date`, a
[`DayOfYear`](@ref), or a `Unitful` time, which is read as the time elapsed since ``T_0``:
```@docs; canonical=false
TimeLike
DayOfYear
```
```@raw html
<iframe src="../figs/eph.html" width="800" height="600"></iframe>
```
```@raw latex
\begin{figure}[htbp]
\centering
\includegraphics[width=0.8\textwidth]{figs/eph.pdf}
\caption{The Earth's orbit (the ecliptic, at 1 au) in Galactic-Rectangular coordinates, with the Earth at the March equinox $t_\mathrm{eq}$ and at $T_0$. The velocities $\vec{v}_\odot$ and $\vec{v}_{ES}$ are drawn to a common scale.}
\end{figure}
```
([Fig. `eph.asy` source](figs/eph.asy))

```@docs; canonical=false
T0
J2000
vES
vE
VE
vE(m::HaloModel, t::TimeLike)
```

The annual mean of the Earth's speed, and the time at which it is reached, are also
available, though not exported (as `StandardHaloModelsPP.VEavg` and
`StandardHaloModelsPP.Tavg`):
```@docs; canonical=false
VEavg
Tavg
```

## Lab-Frame Distributions

The lab-frame distributions, transforming the velocity argument with the model's ephemeris
([`vE(m, t)`](@ref vE)):
```@docs; canonical=false
fχlab(m::HaloModel, t::TimeLike, v::AbstractVector{<:Unitful.Velocity})
fRlab(m::HaloModel, t::TimeLike, v::AbstractVector{<:Unitful.Velocity})
```
To evaluate either at many velocities, see [Performance Tips](@ref).

!!! todo "TODO: Sidereal Motion of the Lab"
    Analyses involving daily (sidereal) modulation and directional signals need the lab's motion and
    orientation with the Earth's rotation, as a function of location on its surface and the time of day.
