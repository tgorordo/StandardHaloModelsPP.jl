using BenchmarkTools, StandardHaloModelsPP, NaturallyUnitful
using LinearAlgebra: norm
using Measurements: value
using Dates: Day

const v3 = [100.0, -50.0, 25.0]u"km/s"
const vs = [v3 .* (i / 200) for i = 1:1000]  # 1000 points, in and out of the escape sphere
const ts = [T0 + Day(d) for d = 0:30:330]    # 12 times, for sweeps
const Mc = value(SHMpp)                      # β ≠ 0, η ≠ 0, but no Measurement

const SUITE = BenchmarkGroup()

## Cases:
#   "SHM"      exact, η = 0: the integrals take the closed forms
#   "central"  value(SHMpp): the numerical paths, on plain numbers
#   "SHMpp"    Measurement + anisotropy, the worst case; integrals pass `uncertainty = true`

# Pointwise distributions, galactic and lab frame.
for (name, m) in (("SHM", SHM), ("central", Mc), ("SHMpp", SHMpp))
    SUITE["fχ"][name] = @benchmarkable fχ($m, $v3)
    SUITE["fχlab"][name] = @benchmarkable fχlab($m, $v3)
end
# A call per point recomputes the normalizations; a broadcast freezes them once.
for (name, m) in (("SHM", SHM), ("SHMpp", SHMpp))
    SUITE["fχ × 1000"][name*" per-call"] = @benchmarkable [fχ($m, v) for v in $vs]
    SUITE["fχ × 1000"][name*" broadcast"] = @benchmarkable fχ.($m, $vs)
end

# Integrals.
const CASES =
    (("SHM", SHM, (;)), ("central", Mc, (;)), ("SHMpp", SHMpp, (; uncertainty = true)))
for (name, m, kw) in CASES
    SUITE["g"][name] = @benchmarkable g($m, $(200u"km/s"); $kw...)
    SUITE["fχspeed"][name] = @benchmarkable fχspeed($m, $(300u"km/s"); $kw...)
    SUITE["fχplanar"][name] =
        @benchmarkable fχplanar($m, $(100u"km/s"), $(-50u"km/s"); $kw...)
end
# g far into its tail: a tiny integral, over a small cap of directions near vmaxlab
SUITE["g"]["central tail"] = @benchmarkable g($Mc, $(650u"km/s"))
# 3D cubature vs nested quadrature, whose inner fχspeed is analytic (SHM) or not
for (name, m) in (("SHM", SHM), ("central", Mc))
    SUITE["velint"][name] = @benchmarkable velint($m, (v, p) -> 1.0)
    SUITE["speedint"][name] = @benchmarkable speedint($m, (v, p) -> 1.0)
end

# Problem API: remaking a problem vs solving it, and a sweep over times as an
# EnsembleProblem, serially and on threads.
const kg = v -> 1 / norm(v)
const hp = HaloIntegralProblem(SHM, kg; vmin = 200u"km/s")
const ens = EnsembleProblem(hp; prob_func = (p, ctx) -> remake(p; t = ts[ctx.sim_id]))
SUITE["problem"]["remake t"] = @benchmarkable remake($hp; t = $(T0 + Day(90)))
SUITE["problem"]["solve"] = @benchmarkable solve($hp)
SUITE["sweep × 12"]["serial"] =
    @benchmarkable solve($ens, nothing, EnsembleSerial(); trajectories = $(length(ts)))
SUITE["sweep × 12"]["threads"] =
    @benchmarkable solve($ens, nothing, EnsembleThreads(); trajectories = $(length(ts)))

"""
    runbench(; save = nothing, compare = nothing, tune = true, kwargs...)

Tune and run `SUITE`, and print allocation count, memory and minimum time per case,
laid out like `Test`'s summary.

- `save`: a file to save this run's minimum estimates to, to compare a later run against.
- `compare`: a file saved by an earlier run. Adds each case's time relative to it, marked
  `+` (regression) or `-` (improvement) when beyond BenchmarkTools' default 5% tolerance.
- `tune`: tune the evaluations per sample first, which sub-microsecond cases need to be
  timed accurately. Further `kwargs` (e.g. `seconds = 1` for a quick pass) go to `run`.

Start Julia with threads (`-t auto`, as `just bench` does) for threaded sweeps to
correctly run in parallel.
"""
function runbench(; save = nothing, compare = nothing, tune = true, kwargs...)
    tune && tune!(SUITE)
    res = minimum(run(SUITE; verbose = false, kwargs...))
    save === nothing || BenchmarkTools.save(save, res)
    ref = compare === nothing ? nothing : only(BenchmarkTools.load(compare))

    head = ["Benchmark Summary:", "Allocs", "Memory", "Min Time"]
    ref === nothing || push!(head, "vs Ref")
    rows = Vector{String}[]
    for grp in sort(collect(keys(res)))
        push!(rows, [grp])
        for case in sort(collect(keys(res[grp])))
            t = res[grp][case]
            row = [
                "  " * case,
                string(t.allocs),
                BenchmarkTools.prettymemory(t.memory),
                BenchmarkTools.prettytime(t.time),
            ]
            ref === nothing || push!(row, _vsref(t, ref, grp, case))
            push!(rows, row)
        end
    end

    # Column widths from content, as `Test` does, rather than fixed padding.
    w = [
        maximum(r -> length(r) ≥ i ? textwidth(r[i]) : 0, [[head]; rows]) for
        i in eachindex(head)
    ]
    line(r) =
        length(r) == 1 ? rpad(r[1], w[1]) * " |" :
        rpad(r[1], w[1]) * " | " * join((lpad(r[i], w[i]) for i = 2:length(r)), "  ")

    println(line(head))
    foreach(r -> println(line(r)), rows)
    return res
end

# A case's time relative to the reference run: "1.03x", marked + / - beyond the tolerance.
function _vsref(t, ref, grp, case)
    (haskey(ref, grp) && haskey(ref[grp], case)) || return "new"
    r = ref[grp][case]
    j = judge(t, r)
    mark = isregression(j) ? " +" : isimprovement(j) ? " -" : "  "
    return string(round(t.time / r.time, digits = 2), "x", mark)
end
