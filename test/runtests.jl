using StandardHaloModelsPP
using NaturallyUnitful
using LinearAlgebra: norm, ⋅
using Measurements: Measurement, uncertainty, value
using Dates: Date, DateTime, Day, Millisecond, month
using StaticArrays: SA
using Test

const S = StandardHaloModelsPP

# callable-object kernels: one arity, the other, and both
struct InvSpeed end
(::InvSpeed)(v) = 1 / norm(v)
struct ScaledInvSpeed
    a::Float64
end
(k::ScaledInvSpeed)(v, p) = k.a * p / norm(v)
struct InvSpeedOr end
(k::InvSpeedOr)(v) = k(v, 1.0)
(::InvSpeedOr)(v, a) = a / norm(v)

@testset "StandardHaloModelsPP.jl" begin
    include("runaqua.jl")

    @testset "Distributions" begin
        @testset "constructor validation" begin
            ρ₀, v₀, vₑ = 0.3u"GeV/cm^3", 220u"km/s", 544u"km/s"
            for args in (
                (0.3, v₀, vₑ, 0, 0),             # ρ₀ not a density
                (ρ₀, 220, vₑ, 0, 0),             # v₀ not a velocity
                (ρ₀, v₀, 544u"km", 0, 0),        # vₑ not a velocity
                (ρ₀, v₀, vₑ, 0, 0, 230),         # VLSR not a velocity
                (ρ₀, v₀, vₑ, 0.5u"km/s", 0),     # β not a plain number
                (ρ₀, v₀, vₑ, 0, 0.5u"km/s"),     # η not a plain number
                (-ρ₀, v₀, vₑ, 0, 0),             # ρ₀ non-positive
                (ρ₀, v₀, vₑ, 1, 0),              # β out of range
                (ρ₀, v₀, vₑ, 0, 1.1),            # η out of range
                (ρ₀, -v₀, vₑ, 0, 0),             # v₀ non-positive
                (ρ₀, v₀, -vₑ, 0, 0),             # vₑ non-positive
                (ρ₀, v₀, vₑ, 0, 0, -v₀),         # VLSR non-positive
            )
                @test_throws ArgumentError HaloModel(args...)
            end
            # speeds are promoted to v₀'s units; a mass density is kept as given
            @test HaloModel(ρ₀, v₀, 544_000u"m/s", 0, 0).vₑ === 544.0u"km/s"
            @test HaloModel(5e-22u"kg/m^3", v₀, vₑ, 0, 0).ρ₀ === 5e-22u"kg/m^3"
        end

        @testset "keyword and copy constructors" begin
            @test HaloModel(; ρ₀ = 0.3u"GeV/cm^3", v₀ = 238u"km/s", vₑ = 544u"km/s") === SHM
            @test HaloModel(;
                localdensity = 0.3u"GeV/cm^3",
                circspeed = 238u"km/s",
                escspeed = 544u"km/s",
            ) === SHM
            @test HaloModel(SHM; v₀ = 220u"km/s") === SHMlegacy
            # a missing, unknown or repeated parameter
            @test_throws ArgumentError HaloModel(; v₀ = 238u"km/s")
            @test_throws ArgumentError HaloModel(SHM; vesc = 500u"km/s")
            @test_throws ArgumentError HaloModel(SHM; v₀ = 1u"km/s", circspeed = 2u"km/s")
            # unchanged Measurement parameters are the same objects (correlation kept)
            m = HaloModel(SHMpp; v₀ = 220u"km/s")
            @test m.β === SHMpp.β && m.vₑ === SHMpp.vₑ && value(m.v₀) == 220u"km/s"
            # ConstructionBase (Accessors.jl's @set) sees only the fields, but takes aliases
            CB = S.ConstructionBase
            @test keys(CB.getproperties(SHM)) == fieldnames(HaloModel)
            @test CB.setproperties(SHM, (; v₀ = 220u"km/s")) === SHMlegacy
        end

        @testset "LSR speed" begin
            # tied to v₀ unless given
            @test SHM.VLSR === SHM.v₀ && SHMpp.VLSR === SHMpp.v₀
            @test HaloModel(SHM; v₀ = 250u"km/s").VLSR == 250u"km/s"
            m = HaloModel(SHM; VLSR = 230u"km/s")
            @test m.VLSR == m.lsrspeed == 230u"km/s" && m.v₀ == SHM.v₀
            @test HaloModel(m; v₀ = 250u"km/s").VLSR == 230u"km/s"
            @test HaloModel(m; VLSR = nothing) === SHM
            # promoted to v₀'s units and type, and stripped by value
            @test HaloModel(SHM; VLSR = 230_000u"m/s") === m
            @test HaloModel(SHMpp; VLSR = 238u"km/s").VLSR isa Quantity{<:Measurement}
            @test value(HaloModel(SHMpp; VLSR = VLSR)).VLSR === 238.0u"km/s"
        end

        @testset "show" begin
            @test sprint(show, SHMpp) == "SHMpp"      # presets print by name
            @test sprint(show, fχ(SHM)) == "fχ(SHM)"
            # a modified copy is not a preset, so it prints its fields
            c = HaloModel(SHMpp; ρ₀ = 0.3u"GeV/cm^3")
            @test startswith(sprint(show, c), "HaloModel(ρ₀ = ")
            full = sprint(show, MIME"text/plain"(), SHMpp)
            @test startswith(full, "HaloModel SHMpp:") && occursin("sausagefrac", full)
            # a tied LSR speed shows as v₀, an untied one by value
            @test occursin("VLSR = v₀", full) && !occursin("VLSR", sprint(show, c))
            u = HaloModel(SHM; VLSR = 230u"km/s")
            @test occursin("VLSR = 230", sprint(show, u))
            @test occursin("VLSR = 230", sprint(show, MIME"text/plain"(), u))
        end

        @testset "escape-speed truncation" begin
            for m in (SHM, SHMpp)
                @test value(ustrip(fR(m, 500.0u"km/s"))) > 0
                @test iszero(fR(m, 544.1u"km/s"))
            end
            @test ustrip(fR(SHMnoesc, 1000.0u"km/s")) > 0
        end

        @testset "fχ composition identity" begin
            v = [100.0, -50.0, 25.0]u"km/s"
            for m in (SHM, SHMpp)
                @test fχ(m, v) ≈ (1 - m.η) * fR(m, v) + m.η * fS(m, v)
            end
        end

        # β = 0 makes the sausage isotropic, so it must coincide with the round component;
        # β ≠ 0 must not (or the β = 0 checks would pass vacuously)
        @testset "fS vs fR" begin
            v = [100.0, -50.0, 25.0]u"km/s"
            @test fS(SHM, v) ≈ fR(SHM, v)
            @test fS(SHMnoesc, v) ≈ fR(SHMnoesc, v)
            @test !(
                fS(SHMpp, [300.0, 0.0, 0.0]u"km/s") ≈ fS(SHMpp, [0.0, 300.0, 0.0]u"km/s")
            )
        end

        # the sausage normalization once overflowed (erfi) as β → 1
        @testset "sausage normalization at β → 1" begin
            v = [100.0, 0.0, 0.0]u"km/s"
            for (β, vₑ) in ((0.9, 3400), (0.9976, 544), (0.999, 544))
                m = HaloModel(0.3u"GeV/cm^3", 238u"km/s", vₑ * u"km/s", β, 1.0)
                @test 0 < ustrip(fS(m, v)) < Inf
            end
        end

        # The frozen forms, and the broadcast specializations built on them, must agree with
        # the plain calls: they only reorder when the normalizations are computed.
        @testset "frozen forms and broadcasting" begin
            vs = [[100.0, -50.0, 25.0]u"km/s", [600.0, 0.0, 0.0]u"km/s"] # in, out of sphere
            ss = [150.0u"km/s", 700.0u"km/s"]
            ts = [T0, T0 + Day(100)]
            for m in (SHM, SHMpp)
                for f in (fR, fS, fχ, fRlab, fχlab)
                    @test f(m).(vs) == f.(m, vs) == [f(m, v) for v in vs]
                end
                @test fR(m).(ss) == fR.(m, ss) == [fR(m, s) for s in ss]
                for f in (fRlab, fχlab)
                    @test f(m)(ts[2], vs[1]) == f(m, ts[2], vs[1])
                    @test f.(m, ts, Ref(vs[1])) == [f(m, t, vs[1]) for t in ts]
                    # frozen at a time too: ≈, since the same arithmetic compiled in a
                    # different context can differ in the last bit
                    @test f(m, ts[2]).(vs) ≈ [f(m, ts[2], v) for v in vs] rtol = 1e-14
                    @test f.(m, ts[2], vs) ≈ [f(m, ts[2], v) for v in vs] rtol = 1e-14
                end
            end
            # the default model
            @test fχ.(vs) == [fχ(v) for v in vs]
            @test fR.(ss) == [fR(s) for s in ss]
            @test fχlab.(ts, Ref(vs[1])) == [fχlab(t, vs[1]) for t in ts]
            # with no velocities, a broadcast is the frozen form
            for f in (fR, fS, fχ, fRlab, fχlab)
                @test f.(SHM)(vs[1]) == f(SHM, vs[1])
            end
            for f in (fRlab, fχlab)
                @test f.(SHM, T0)(vs[1]) == f(SHM, T0, vs[1])
            end
        end

        # plain Vector and SVector velocities, and other units, give the same answer
        @testset "Vector/SVector and unit invariance" begin
            vv, sv = [100.0, -50.0, 25.0]u"km/s", SA[100.0, -50.0, 25.0]u"km/s"
            for m in (SHM, SHMpp)
                @test fχlab(m, vv) ≈ fχlab(m, sv)
                @test fχ(m, vv) ≈ fχ(m, uconvert.(u"m/s", vv))
            end
        end

        # fχlab/fRlab are the galactic-frame forms boosted by the model's v_E; exact models
        # use the central-value ephemeris
        @testset "lab frame is the boosted galactic frame" begin
            v, t = [100.0, -50.0, 25.0]u"km/s", T0
            for m in (SHM, SHMpp)
                V = vE(m, t)
                @test fχlab(m, t, v) ≈ fχ(m, v + V)
                @test fRlab(m, t, v) ≈ fR(m, v + V)
                @test uncertainty(ustrip(fχlab(m, t, v))) ≈
                      uncertainty(ustrip(fχ(m, v + V)))
            end
            @test iszero(uncertainty(ustrip(fχlab(SHM, t, v))))
        end

        # fχplanar marginalizes fχlab over v_z; checked against an independent trapezoid sum
        @testset "fχplanar" begin
            function trapz_vz(m, vr, vθ; N = 1601, L = 800.0)
                dz = 2L / (N - 1)
                s = zero(fχlab(m, [vr, vθ, 0.0u"km/s"]))
                for i = 1:N
                    w = (i == 1 || i == N) ? 0.5 : 1.0
                    s += w * fχlab(m, [vr, vθ, (-L + (i - 1) * dz) * u"km/s"])
                end
                return s * dz * u"km/s"
            end
            for m in (SHM, SHMpp)
                for (vr, vθ) in ((0.0, 0.0), (100.0, -50.0), (200.0, 150.0))
                    @test fχplanar(m, vr * u"km/s", vθ * u"km/s") ≈
                          trapz_vz(m, vr * u"km/s", vθ * u"km/s") rtol = 1e-3
                end
                # outside the escape sphere the discriminant is negative and returns early
                @test iszero(fχplanar(m, 900.0u"km/s", 900.0u"km/s"))
                # the same point in other units gives the same result, to solver precision
                a = fχplanar(m, 100u"km/s", -50u"km/s")
                b = fχplanar(m, 100_000u"m/s", -50_000u"m/s")
                @test a ≈ uconvert(unit(a), b) rtol = 1e-5
            end
            # vr and vθ need not share a numeric type
            @test fχplanar(SHM, 100u"km/s", -50.0u"km/s") ≈
                  fχplanar(SHM, 100u"km/s", -50u"km/s")
        end
    end

    @testset "Ephemeris" begin
        @testset "isexact and value" begin
            @test S.isexact(SHM) && !S.isexact(SHMpp) && !S.isexact(SHMerr)
            @test S.isexact(value(SHMpp))
            @test value(SHMerr) == SHM
            @test value(SHMpp) ==
                  HaloModel(0.55u"GeV/cm^3", 233u"km/s", 528u"km/s", 0.9, 0.2)
        end

        # vE(m, t) follows the model's precision: an exact model gets the central-value
        # ephemeris, so its results stay plain numbers
        @testset "model-aware ephemeris" begin
            t = T0 + Day(45)
            @test !(eltype(vE(SHM, t)) <: Quantity{<:Measurement})
            @test eltype(vE(SHMpp, t)) <: Quantity{<:Measurement}
            @test vE(SHM, t) ≈ value.(vE(t))
        end

        # boost uses the model's LSR speed; the reference ephemeris (no model) is SHM's
        @testset "model LSR speed" begin
            t = T0 + Day(45)
            ŷ = SA[0.0, 1.0, 0.0]
            @test vE(SHMlegacy, t) ≈ vE(SHM, t) - 18u"km/s" * ŷ
            @test vE(HaloModel(SHM; VLSR = 230u"km/s"), t) ≈ vE(SHM, t) - 8u"km/s" * ŷ
            # SHMpp boosts with its own v₀, so the boost is correlated with the dispersion
            Vy, v₀ = ustrip(u"km/s", vE(SHMpp, t)[2]), ustrip(u"km/s", SHMpp.v₀)
            @test S.Measurements.derivative(Vy, v₀) ≈ 1
            # SHMerr's v₀ is the reference VLSR itself, so its boost is the reference one,
            # uncertainties included
            @test SHMerr.v₀ === VLSR
            a, b = ustrip.(u"km/s", vE(SHMerr, t)), ustrip.(u"km/s", vE(t))
            @test value.(a) ≈ value.(b) && uncertainty.(a) ≈ uncertainty.(b)
            # the analytic path boosts the same way as the numerical one
            @test fχspeed(SHMlegacy, t, 300u"km/s") ≈
                  fχspeed(SHMlegacy, t, 300u"km/s"; analytic = false) rtol = 1e-4
        end

        @testset "time arguments" begin
            @test VE(Date(2020, 6, 1)) == VE(DateTime(2020, 6, 1))
            @test VE(DayOfYear(153)) == VE(DateTime(2020, 6, 1))  # 2020 is a leap year
            @test VE(45u"d") == VE(T0 + Day(45))  # Unitful times count from T0
            @test DayOfYear(1.5, 2021) == DayOfYear(DateTime(2021, 1, 1, 12))
            @test sprint(show, DayOfYear(153)) == "DayOfYear(153.0, 2020)"
            # all accepted wherever a time is
            v = [100.0, -50.0, 25.0]u"km/s"
            @test fχlab(SHM, Date(2020, 6, 1), v) ==
                  fχlab(SHM, DayOfYear(153), v) ==
                  fχlab(SHM, DateTime(2020, 6, 1), v)
        end

        # T0 is *defined* as the phase at which the Earth's speed through the halo peaks
        # (Ephemeris.jl). The amplitude is external ground truth: the Earth's 29.8 km/s
        # orbital speed projected onto v̂SUN gives ~29 km/s peak-to-peak.
        @testset "annual modulation of VE" begin
            VEd(d) = value(ustrip(u"km/s", VE(d * u"d")))
            lo, hi = extrema(VEd(d) for d = 0:364)
            @test VEd(0) ≈ hi                                    # T0 is the maximum
            @test 25 < hi - lo < 32                              # cf. 2·VES·cos(~61°)
            # VEavg is the annual mean, not merely some value in range
            @test value(ustrip(u"km/s", S.VEavg)) ≈ sum(VEd, 0:364) / 365 rtol = 1e-4
            # Tavg is the rising crossing of the mean, about March 2nd (Lewin & Smith)
            @test VEd(ustrip(u"d", S.Tavg)) ≈ value(ustrip(u"km/s", S.VEavg)) rtol = 1e-6
            @test month(T0 + round(Millisecond, S.Tavg)) in (2, 3)
        end
    end

    @testset "Integrals" begin
        vmin = 200u"km/s"

        @testset "analytic vs numerical" begin
            @test g(SHM, vmin; analytic = true) ≈ g(SHM, vmin; analytic = false) rtol = 1e-3
            @test h(SHM, vmin; analytic = true) ≈ h(SHM, vmin; analytic = false) rtol = 1e-3
            @test fχspeed(SHM, 300u"km/s"; analytic = true) ≈
                  fχspeed(SHM, 300u"km/s"; analytic = false) rtol = 1e-3
            # vₑ = Inf: the escape correction must be skipped, not evaluated as 0 * Inf
            for f in (g, h)
                @test f(SHMnoesc, vmin; analytic = true) ≈
                      f(SHMnoesc, vmin; analytic = false) rtol = 1e-3
            end
            # η ≠ 0 has no closed form, so there the flag is a silent no-op
            @test g(SHMpp, vmin; analytic = true) ≈ g(SHMpp, vmin; analytic = false)
            # Near vmaxlab only a small cap of directions has support. The 3D cubature spans
            # just that cap; over the whole sphere it used to miss it and return exactly 0.
            vmax = vmaxlab(SHM, T0)
            for vm in (700u"km/s", vmax - 20u"km/s", vmax - 0.5u"km/s"), f in (g, h)
                @test f(SHM, vm; analytic = false) ≈ f(SHM, vm) rtol = 1e-4
            end
            # ... and converges to tight tolerances there
            a = g(SHM, 770u"km/s")
            @test g(SHM, 770u"km/s"; analytic = false, reltol = 1e-8) ≈ a rtol = 1e-8
        end

        # HaloIntegralProblem is the single code path under velint/speedint, so these check
        # the problem API itself: construction, solve, remake and the solution.
        @testset "HaloIntegralProblem" begin
            kk(v, p) = 1 / norm(v)
            pr = HaloIntegralProblem(SHM, kk; vmin)
            sol = solve(pr)
            @test S.SciMLBase.successful_retcode(sol)
            @test unit(sol.u) == u"s/km" && unit(sol.resid) == u"s/km"
            @test sol.u == velint(SHM, kk; vmin)
            # a solver stopped short of the tolerance is not reported as Success (the raw
            # Integrals.jl QuadGK/HCubature solutions say Success regardless)
            @test solve(pr; maxiters = 1).retcode == S.SciMLBase.ReturnCode.MaxIters
            # remake: equal to a fresh build, and actually changes what it is given
            t2 = T0 + Day(180)
            r = solve(remake(pr; t = t2)).u
            @test r == solve(HaloIntegralProblem(SHM, t2, kk; vmin)).u
            @test r != sol.u
            @test solve(remake(pr; vmin = 500u"km/s")).u ≈ g(SHM, 500u"km/s") rtol = 1e-3
            # the speed measure
            ps = HaloIntegralProblem(SHM, (v, p) -> 1 / v; measure = :speed, vmin)
            @test solve(ps).u == speedint(SHM, (v, p) -> 1 / v; vmin)
            @test_throws ArgumentError HaloIntegralProblem(SHM, kk; measure = :planar)
            # uncertainty is fixed at construction, and survives remake
            pu = HaloIntegralProblem(SHMpp, kk; vmin, uncertainty = true)
            @test uncertainty(solve(pu).u) > 0u"s/km"
            @test uncertainty(solve(remake(pu; t = t2)).u) > 0u"s/km"
            @test startswith(sprint(show, pu), "HaloIntegralProblem(SHMpp, ")
            @test occursin("uncertainty = true", sprint(show, pu))
        end

        # Sweeps are SciML ensembles of remade HaloIntegralProblems; they must agree with
        # solving each point directly, serially or on threads.
        @testset "sweeps (EnsembleProblem)" begin
            prob = HaloIntegralProblem(SHM, v -> 1 / norm(v); vmin)
            ts = [T0 + Day(d) for d = 0:60:300]
            ens =
                EnsembleProblem(prob; prob_func = (p, ctx) -> remake(p; t = ts[ctx.sim_id]))
            direct = [solve(remake(prob; t)).u for t in ts]
            for mode in (EnsembleSerial(), EnsembleThreads())
                sols = solve(ens, nothing, mode; trajectories = length(ts))
                @test [s.u for s in sols.u] == direct
            end
            # kernel parameters sweep the same way
            ps = [1.0, 2.0, 3.0]
            pk = HaloIntegralProblem(SHM, (v, p) -> p / norm(v); vmin, p = 1.0)
            ensp =
                EnsembleProblem(pk; prob_func = (q, ctx) -> remake(q; p = ps[ctx.sim_id]))
            sp = solve(ensp, nothing, EnsembleSerial(); trajectories = 3)
            @test [s.u for s in sp.u] ≈ ps .* solve(pk).u
        end

        @testset "frozen speed distribution and moments" begin
            t, v = T0 + Day(45), 300u"km/s"
            # the analytic (SHM) and numerical (value(SHMpp)) paths
            for m in (SHM, value(SHMpp))
                @test fχspeed(m, t)(v) == fχspeed(m, t, v)
                for f in (g, h)
                    @test f(m, t)(vmin) == f(m, t, vmin)
                end
            end
            @test fχspeed(SHM, t; analytic = false)(v) ==
                  fχspeed(SHM, t, v; analytic = false)
            # Measurement models: SHMerr is analytic, so its hoisted constants must keep the
            # model/ephemeris correlation
            for m in (SHMerr, SHMpp)
                a = g(m, t; uncertainty = true)(v)
                @test uncertainty(a) > 0u"s/km"
                @test a == g(m, t, v; uncertainty = true)
            end
            @test fχspeed(SHMpp, t; uncertainty = true)(v) ==
                  fχspeed(SHMpp, t, v; uncertainty = true)
            # quiet: a frozen form inside someone else's integrand must not warn per call
            @test_logs fχspeed(value(SHMpp), t; maxiters = 1)(v)
            @test sprint(show, fχspeed(SHM, t)) == "fχspeed(SHM, $t)"
            @test sprint(show, g(SHM, t)) == "g(SHM, $t)"
        end

        # the largest lab-frame speed with support, with the 6v₀ cutoff when vₑ = ∞
        @testset "vmaxlab" begin
            t = T0 + Day(45)
            @test iszero(fχspeed(SHM, t, vmaxlab(SHM, t) + 1u"km/s"))
            @test !iszero(fχspeed(SHM, t, vmaxlab(SHM, t) - 1u"km/s"))
            @test vmaxlab(SHMnoesc, t) ≈ 6SHMnoesc.v₀ + VE(SHMnoesc, t)
        end

        @testset "convergence warnings" begin
            # a solver that gives up is reported by the wrappers; no silent return
            @test_logs (:warn, r"did not converge") g(SHMpp, 200u"km/s"; maxiters = 1)
            @test_logs (:warn, r"did not converge") velint(SHM, (v, p) -> 1.0; maxiters = 1)
            @test_logs (:warn, r"did not converge") fχspeed(SHMpp, 300u"km/s"; maxiters = 1)
            vr, vθ = 100u"km/s", -50u"km/s"
            @test_logs (:warn, r"did not converge") fχplanar(SHM, vr, vθ; maxiters = 1)
            @test_logs g(SHMpp, 200u"km/s")  # and quiet when it converges

            # a needle-thin sausage warns from the numerical integrals, but not from the
            # pointwise forms, nor when there is no sausage to resolve (η = 0)
            mβ = HaloModel(0.3u"GeV/cm^3", 238u"km/s", 544u"km/s", 0.995, 0.2)
            @test_logs (:warn, r"β > 0.99") g(mβ, 300u"km/s")
            @test_logs (:warn, r"β > 0.99") fχspeed(mβ, 300u"km/s")
            @test_logs fχ(mβ, [100.0, 0.0, 0.0]u"km/s")
            @test_logs g(HaloModel(mβ; η = 0), 300u"km/s"; analytic = false)
        end

        @testset "kernels" begin
            # kernels may omit p
            @test velint(SHM, v -> 1 / norm(v); vmin) ==
                  velint(SHM, (v, p) -> 1 / norm(v); vmin)
            @test speedint(SHM, v -> 1 / v; vmin) == speedint(SHM, (v, p) -> 1 / v; vmin)
            @test_throws ArgumentError HaloIntegralProblem(SHM, () -> 1.0)
            # a model and/or time with no kernel is a missing argument
            for f in (velint, speedint, HaloIntegralProblem),
                args in ((SHM,), (T0,), (SHM, T0))

                @test_throws MethodError f(args...)
            end
            # callable objects, functions with unrelated 2-argument methods (norm(v, p)) and
            # compositions (any arity) all work
            @test velint(SHM, InvSpeed(); vmin) == velint(SHM, v -> 1 / norm(v); vmin)
            @test velint(SHM, norm) == velint(SHM, v -> norm(v))
            @test velint(SHM, inv ∘ norm; vmin) == velint(SHM, v -> inv(norm(v)); vmin)
            # p selects the arity: omitted → kk(v), given → kk(v, p)
            @test velint(SHM, ScaledInvSpeed(2.0); vmin, p = 3.0) ≈
                  6 * velint(SHM, InvSpeed(); vmin) rtol = 1e-12
            @test velint(SHM, InvSpeedOr(); vmin) == velint(SHM, InvSpeed(); vmin)
            @test velint(SHM, InvSpeedOr(); vmin, p = 3.0) ≈
                  3 * velint(SHM, InvSpeed(); vmin) rtol = 1e-12
        end

        # a Quantity abstol is converted to the result's units, so it means the same thing
        # whatever units the inputs are in; a plain number means that number in those units
        @testset "abstol units" begin
            a = g(SHMpp, 200u"km/s"; abstol = 1e-9u"s/km")
            b = g(SHMpp, 200_000u"m/s"; abstol = 1e-9u"s/km")
            @test unit(a) == u"s/km"
            @test a ≈ uconvert(unit(a), b) rtol = 1e-5
            # including for fχspeed, whose cubature covers only half the angular domain
            @test fχspeed(value(SHMpp), 300u"km/s"; abstol = 1e-6) ==
                  fχspeed(value(SHMpp), 300u"km/s"; abstol = 1e-6u"s/km")
        end

        # A legal model with vₑ below the Earth's speed: some rays miss the escape sphere
        # and others cross it between two positive roots. This used to throw a DomainError.
        @testset "escape speed below Earth's speed" begin
            m = HaloModel(0.3u"GeV/cm^3", 220u"km/s", 150u"km/s", 0, 0)
            a = g(m, 0u"km/s"; analytic = false)
            @test 0 < ustrip(a) < Inf
            @test velint(m, (v, p) -> 1.0; reltol = 1e-6) ≈ 1 rtol = 1e-4
            # the closed forms assume V_E < vₑ (and are ~35% off), so the default
            # `analytic = true` must fall back to the numerical path, which agrees with an
            # independent integral
            @test g(m, 0u"km/s") ≈ a rtol = 1e-3
            @test speedint(m, (v, p) -> 1 / v; analytic = false) ≈ a rtol = 1e-3
            @test fχspeed(m, 100u"km/s") ≈ fχspeed(m, 100u"km/s"; analytic = false)
        end

        @testset "normalization" begin
            for m in (SHM, SHMpp)
                @test velint(m, (v, p) -> 1.0) ≈ 1 rtol = 1e-3
                @test speedint(m, (v, p) -> 1.0) ≈ 1 rtol = 1e-3
            end
            @test speedint(SHM, (v, p) -> 1.0; analytic = false) ≈ 1 rtol = 1e-3
        end

        # the speed measure's inner quadrature, analytic or not, against the closed form
        @testset "speedint vs g" begin
            @test speedint(SHM, (v, p) -> 1 / v; vmin) ≈
                  speedint(SHM, (v, p) -> 1 / v; vmin, analytic = false) rtol = 1e-3
            @test speedint(SHM, (v, p) -> 1 / v; vmin) ≈ g(SHM, vmin) rtol = 1e-3
        end

        # g/h use an internal radial kernel that reuses the speed instead of norm(v); at
        # matched tolerances the two must agree to solver precision.
        @testset "radial kernel matches vector kernel" begin
            tol = (; reltol = 1e-6, abstol = 1e-10)
            @test g(SHM, vmin; analytic = false, tol...) ≈
                  velint(SHM, (v, p) -> 1 / norm(v); vmin, tol...) rtol = 1e-6
            @test h(SHM, vmin; analytic = false, tol...) ≈
                  velint(SHM, (v, p) -> norm(v); vmin, tol...) rtol = 1e-6
        end

        # Integrals.jl infers an IntegralProblem's in-place-ness at runtime unless told, so
        # an untyped construction anywhere made every integrator return Any.
        @testset "inference" begin
            t, v = T0 + Day(45), 300.0u"km/s"
            for m in (SHM, value(SHMpp))
                @inferred fχspeed(m, t, v)
                for f in (g, h)
                    @inferred f(m, t, vmin)
                end
            end
            @inferred fχplanar(SHM, t, 100u"km/s", -50u"km/s")
            @inferred velint(SHM, t, v -> 1 / norm(v); vmin)
            @inferred speedint(SHM, t, v -> 1 / v; vmin)
        end

        # the default integrates central values: same value, no uncertainty, for every
        # integrator; `uncertainty = true` propagates them
        @testset "uncertainty propagation" begin
            for (f, args) in (
                (g, (vmin,)),
                (h, (vmin,)),
                (fχspeed, (300u"km/s",)),
                (fχplanar, (100u"km/s", -50u"km/s")),
            )
                a, b = f(SHMpp, args...), f(SHMpp, args...; uncertainty = true)
                @test iszero(uncertainty(a)) && !iszero(uncertainty(b))
                @test value(a) ≈ value(b) rtol = 1e-3
            end
            @test iszero(uncertainty(velint(SHMpp, (v, p) -> 1.0)))
            @test iszero(uncertainty(speedint(SHMpp, (v, p) -> 1.0)))
        end

        # The escape boundary moves with vₑ, v₀ and the ephemeris, and fχ is nonzero on it
        # (the Maxwellian is truncated, not tapered), so an integral's uncertainty must
        # include the boundary's motion. velint gets this from its Measurement-valued radial
        # Jacobian; _fχspeed and fχplanar integrate between central-value limits and add a
        # first-order boundary term.
        @testset "escape-boundary uncertainty" begin
            # normalization carries no uncertainty: ∂(∫f)/∂θ = ∂1/∂θ = 0
            @test uncertainty(speedint(SHMpp, (v, p) -> 1.0; uncertainty = true)) < 1e-4
            @test uncertainty(velint(SHMpp, (v, p) -> 1.0; uncertainty = true)) < 1e-4
            # two independent constructions of g agree in uncertainty, not just in value
            for vm in (0u"km/s", 500u"km/s")
                a = g(SHMpp, vm; uncertainty = true)
                b = speedint(SHMpp, (v, p) -> 1 / v; vmin = vm, uncertainty = true)
                @test uncertainty(a) > 0u"s/km"                  # not vacuously 0 ≈ 0
                @test value(a) ≈ value(b) rtol = 1e-3
                @test uncertainty(a) ≈ uncertainty(b) rtol = 2e-2
            end
            # fχplanar against a direct Measurement-limit rule: v_z = lo + sΔ, s ∈ [0,1],
            # lo and Δ carrying their uncertainties, Simpson on the smooth interior
            vr, vθ = -300.0u"km/s", -400.0u"km/s"
            V = vE(SHMpp, T0)
            D = S._vesc_eff(SHMpp)^2 - (vr + V[1])^2 - (vθ + V[2])^2
            lo, Δ = -V[3] - sqrt(D), 2sqrt(D)
            N = 2000
            fz(s) = fχ(SHMpp, SA[vr, vθ, lo+s*Δ] + V)
            w(i) = i == 0 || i == N ? 1 : (isodd(i) ? 4 : 2)
            ref = Δ * sum(w(i) * fz(clamp(i / N, 1e-9, 1 - 1e-9)) for i = 0:N) / (3N)
            got = fχplanar(SHMpp, vr, vθ; uncertainty = true)
            @test value(got) ≈ value(ref) rtol = 1e-5
            @test uncertainty(got) ≈ uncertainty(ref) rtol = 1e-3
        end

        # velint builds its own fixed frame and halves nothing, so it independently checks
        # the α-halving in _fχspeed where β ≠ 0 and η ≠ 0 make the angular
        # integrand α-dependent. value(SHMpp) keeps both nonzero.
        @testset "symmetric frame vs velint (anisotropic)" begin
            M = value(SHMpp)
            @test speedint(M, inv; vmin) ≈ velint(M, inv ∘ norm; vmin) rtol = 1e-4
            @test speedint(M, identity; vmin) ≈ velint(M, norm; vmin) rtol = 1e-4
        end

        # The α-halving in _fχspeed is only valid if the angular integrand is even in α.
        # That premise is invisible in central values: the odd term enters through ê2 ⋅ V,
        # which is 0 ± 1.0 km/s for SHMpp. So this must run on SHMpp and compare
        # uncertainties, or it proves nothing.
        @testset "α-evenness premise (SHMpp)" begin
            L = S._labframe(SHMpp, 0.0u"d")
            V = vE(SHMpp, 0.0u"d")
            # the frame premise, and the projection killing the odd component in uncertainty
            @test L.ê2 ⋅ SA[1.0, 0.0, 0.0] ≈ 0 atol = 1e-15
            @test uncertainty(L.ê2 ⋅ V) > 0.1u"km/s"       # the term being removed is real
            @test uncertainty(L.ê2 ⋅ L.Vp) < 1e-12u"km/s"
            @test norm(L.Vp) ≈ norm(V)                      # V' preserves the norm ...
            # ... and the anisotropy axis (componentwise: the arithmetic rebuilds the
            # Measurement's derivative tree, so not ===)
            @test value(L.Vp[1]) == value(V[1]) && uncertainty(L.Vp[1]) == uncertainty(V[1])
            # the consequence: the integrand is even in α, in value AND uncertainty ...
            kern(α, W) =
                let (sγ, cγ) = sincos(0.7), (sα, cα) = sincos(α)
                    L.fp(
                        200.0u"km/s" .*
                        ((sγ * cα) .* L.ê1 .+ (sγ * sα) .* L.ê2 .+ cγ .* L.ê3) + W,
                    )
                end
            a, b = kern(1.1, L.Vp), kern(2π - 1.1, L.Vp)
            @test ustrip(value(a)) > 0 && ustrip(uncertainty(a)) > 0 # inside the sphere
            @test value(a) ≈ value(b) rtol = 1e-12
            @test uncertainty(a) ≈ uncertainty(b) rtol = 1e-12
            # ... and without the projection, evenness fails in uncertainty while still
            # holding in value: exactly the regression this guards
            a, b = kern(1.1, V), kern(2π - 1.1, V)
            @test value(a) ≈ value(b) rtol = 1e-12
            @test !isapprox(uncertainty(a), uncertainty(b); rtol = 1e-6)
        end
    end
end
