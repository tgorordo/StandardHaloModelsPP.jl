# Performance Tips

```@meta
CurrentModule = StandardHaloModelsPP
```

## Evaluating a Distribution at Many Points

Each call like `fχ(m, v)` first computes the model's normalizations, which for a model
carrying `Measurement`s costs several times the evaluation itself. To evaluate repeatedly,
*freeze* the model once, and call the result:
```julia
f = fχ(SHMpp)           # normalizations computed once
f.(vs)                  # or f(v) in a loop, a closure, your own integrand, ...

flab = fχlab(SHMpp, t)  # also fixes the time, and so the Earth's velocity
flab(v)
```
The same works for [`fR`](@ref), [`fS`](@ref), [`fRlab`](@ref) and [`fχlab`](@ref)
(`fχlab(m)` takes `([t,] v)`; `fχlab(m, t)` takes `v`).

Broadcasting over velocities with a single model does this for you: `fχ.(m, vs)`,
`fχlab.(m, t, vs)` and the default-model `fχ.(vs)` all freeze the model once. For 1000
points of [`SHMpp`](@ref) this is about 7x faster than calling `fχ(m, v)` per point.
Broadcasting over an array of models still evaluates each model separately, as usual.

## Uncertainties Cost Time

A model carrying `Measurement`s (such as [`SHMpp`](@ref)) makes every evaluation propagate
uncertainties: about 10x slower pointwise. When you don't need them, use its central-value
model `value(m)` instead, which evaluates as fast as an exact one.

Integrals evaluate central values by default; `uncertainty = true` propagates the model's
uncertainties, at roughly 10-130x the cost (most for [`g`](@ref) and [`h`](@ref)).

## Analytic Fast Paths

For a model without the anisotropic component (``\eta = 0``, i.e. every preset but
`SHMpp`), with the Earth inside its escape sphere, [`fχspeed`](@ref), [`g`](@ref) and
[`h`](@ref) use exact closed forms [DelNobile:2021wmp; Sec. I.4](@cite), which take under
100 ns. Pass `analytic = false` to force the numerical integrals instead (e.g. to check them).

## Choosing Between `velint` and `speedint`

For a kernel that depends only on the speed ``|\mathbf{v}|``, both [`velint`](@ref) (the
3-D velocity integral) and [`speedint`](@ref) (a 1-D integral over the speed distribution)
apply:
- For ``\eta = 0`` models, `speedint` uses the analytic speed distribution and is the
  faster: around 20 µs, against 0.04-0.2 ms for `velint`.
- For anisotropic models (e.g. `SHMpp`), every point of the speed distribution is itself an
  angular integral, and `velint` is faster (about 0.04-2 ms against 0.4-5 ms, both at a
  relative tolerance of ``10^{-4}``). Both get cheaper as ``v_\text{min}`` rises.

`speedint` defaults to a tighter `reltol` (``10^{-6}``) than `velint` (``10^{-4}``). For an
anisotropic model comparing defaults is a misleading benchmark (a `speedint` takes around 50 ms).
Pass `reltol` explicitly when comparing the two, or when the default precision isn't needed.

A kernel depending on the direction of ``\mathbf{v}`` needs `velint`.

## Tolerances

Every integrator defaults to a purely relative tolerance: a `reltol` suited to its method, and
`abstol = 0`. To set an absolute tolerance, give `abstol` as a `Quantity` in the units of the
result. A plain number is taken in the integrand's internal units, so it only makes sense for
integrals of order one: a small integral (e.g. a rate for a small cross-section, or ``g`` far
into its tail) with a plain-number `abstol` can stop well short of its `reltol`.

## Many Integrals: Sweeps and Threads

To compute an integral at many times, ``v_\text{min}`` values or kernel parameters, solve an
`EnsembleProblem` of [`HaloIntegralProblem`](@ref)s on threads, as described under
[Sweeps](@ref). Start Julia with threads (e.g. `julia -t auto`) for `EnsembleThreads()` to
run in parallel.

To integrate against the distribution inside your own integrals, use the forms frozen at a
time, `fχspeed(m, t)` and `fχlab(m, t)`, rather than calling `fχspeed(m, t, v)` or
`fχlab(m, t, v)` at every point (see [Building Your Own Integrals](@ref)). Likewise, to
evaluate [`g`](@ref) or [`h`](@ref) at many ``v_\text{min}`` at one time (e.g. across a
recoil-energy spectrum), use `g(m, t)` and `h(m, t)`: on the analytic path they take around
25-40 ns per call, against 100-130 ns for `g(m, t, vmin)`.
