# Comprehensive unit tests for the MAPPLE spectral model.
#
#   Core model:        src/Mapple.jl    (`mapple`, `mapple!`, `MAPPLE`, `mapple_sort!`,
#                                         `fit`, `predict`)
#   Optim refinement:  ext/OptimExt.jl  (`fit_mapple` 3-arg, `fit!`)
#
# Characterization tests pin down the current, correct behaviour of the model. The later
# items also cover the robustness/reporting overhaul: the conservative peak-amplitude init
# and the `log_A` do-block-leak fix, the Hz reporting accessors, the `rsquared` flat-spectrum
# guard, the peak-width guard, data-driven (`:auto`) peak counting, BIC component
# selection, and the outer-edge windowing fix (final items in this file). These tests are deterministic; every item that fits a model also saves a log–log
# figure of the fit to `test/mapple_figs/` via the shared `MapplePlots` setup module.
#
# Construction helpers are repeated inside each `@testitem` because every item runs in
# its own isolated module.

# Shared plotting helper. Defined once as a `@testmodule` so CairoMakie loads a single
# time and is reused by every item that requests `setup = [MapplePlots]`. Only the
# `ComponentArray` *type* is imported (not all of ComponentArrays) so its `Axis` export does not
# shadow Makie's `Axis`.
@testmodule MapplePlots begin
    using CairoMakie
    using ForwardDiff
    import Fathom: fathom, OnePanel
    import TimeseriesTools: MAPPLE, mapple
    import ComponentArrays: ComponentArray
    set_theme!(fathom())

    const FIGDIR = joinpath(@__DIR__, "mapple_figs")

    _curve(f, m::MAPPLE) = mapple(f, m.params)
    _curve(f, p::ComponentArray) = mapple(f, p)

    """
        save_fit(name, f, data, fit; init = nothing)

    Save a log-log figure of a MAPPLE fit to `test/mapple_figs/<name>.png`: the measured `data`
    (linear spectral density at linear frequencies `f`) as points, the rough `init` (initial
    estimate) curve in grey when supplied, and the refined `fit` (full fit) curve in red on top.
    `fit`/`init` may each be a `MAPPLE` or a raw parameter `ComponentArray`.
    """
    function save_fit(name, f, data, fit; init = nothing, subdir = "")
        outdir = joinpath(FIGDIR, subdir)
        isdir(outdir) || mkpath(outdir)
        fv = collect(Float64, f)
        fig = OnePanel()
        ax = Axis(
            fig[1, 1]; xscale = log10, yscale = log10,
            xlabel = "frequency", ylabel = "power spectral density"
        )
        scatter!(ax, fv, collect(Float64, data); label = "data", markersize = 4, color = (:black, 0.4))
        init === nothing ||
            lines!(ax, fv, _curve(fv, init); label = "initial estimate", color = :gray, linewidth = 2, linestyle = :dash)
        lines!(ax, fv, _curve(fv, fit); label = "full fit", color = :crimson, linewidth = 2, linestyle = :dash)
        # Horizontal legend in a row above the axis: spans the width, never overlaps the data.
        Legend(
            fig[0, 1], ax; orientation = :horizontal, framevisible = false,
            tellheight = true, tellwidth = false
        )
        save(joinpath(outdir, "$(name).png"), fig)
        return nothing
    end
end

@testitem "mapple model: power law, peaks, positivity" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    f = exp10.(range(0, 3, length = 200))

    # A single component with no peaks is exactly the power law A·f^β.
    p1 = ComponentArray(;
        log_A = 0.5, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    s1 = mapple(f, p1)
    @test s1 ≈ exp10(0.5) .* f .^ (-2.0)
    @test all(>(0), s1)                # strictly positive => log10 in the loss is safe
    @test all(isfinite, s1)

    # The two `mapple` call forms agree.
    @test mapple(f, p1) ≈ mapple(f, p1[[:log_A, :transition_width, :components]], p1[[:peaks]])

    # A Gaussian peak adds power near its centre frequency (10^log_f).
    p2 = ComponentArray(;
        log_A = 0.5,
        peaks = [mkpeak(; log_f = 1.0, log_σ = 0.3, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    s2 = mapple(f, p2)
    icentre = argmin(abs.(f .- 10.0))
    @test s2[icentre] > s1[icentre]
    @test all(>(0), s2)

    # Output is invariant to the order components are supplied (sorted internally).
    pfwd = ComponentArray(;
        log_A = 1.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 5.0, β = -3.0)],
        transition_width = 0.05
    )
    prev = ComponentArray(;
        log_A = 1.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = -3.0), mkcomp(; log_f_stop = 1.5, β = -1.0)],
        transition_width = 0.05
    )
    @test mapple(f, pfwd) ≈ mapple(f, prev)

    # The hot path is type stable.
    @test @inferred(mapple(f, p1)) isa Vector{Float64}
end

@testitem "mapple model: multi-component finite & positive" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    f = exp10.(range(-1, 4, length = 400))
    p = ComponentArray(;
        log_A = 1.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 5.0, β = -3.0)],
        transition_width = 0.05
    )
    s = mapple(f, p)
    @test all(isfinite, s)
    @test all(>(0), s)
    @test @inferred(mapple(f, p)) isa Vector{Float64}

    # A peak added on top only increases power.
    pp = ComponentArray(;
        log_A = 1.0,
        peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 2.0)],
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 5.0, β = -3.0)],
        transition_width = 0.05
    )
    @test all(mapple(f, pp) .≥ s .- 1.0e-8)
end

@testitem "MAPPLE: construction & predict" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

    # The keyword constructor mirrors the documented field layout.
    m = MAPPLE(;
        log_A = 1.0,
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5)],
        components = [mkcomp(; log_f_stop = 1.5, β = 2.0)],
        transition_width = 0.2
    )
    @test m isa MAPPLE
    @test length(m.params.peaks) == 1
    @test length(m.params.components) == 1
    @test m.params.log_A == 1.0
    @test m.params.transition_width == 0.2

    # `predict` is exactly the model evaluation at the requested frequencies.
    f = exp10.(range(0, 3, length = 100))
    @test predict(m, f) == mapple(f, m.params)
end

@testitem "MAPPLE: sort reorders blocks without corruption" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    mkparams() = ComponentArray(;
        log_A = 1.0,
        peaks = [
            mkpeak(; log_f = 2.5, log_σ = 0.11, log_A = 0.51),
            mkpeak(; log_f = 0.7, log_σ = 0.12, log_A = 1.52),
        ],
        components = [
            mkcomp(; log_f_stop = 5.0, β = 4.0),
            mkcomp(; log_f_stop = 1.5, β = 2.0),
        ],
        transition_width = 0.2
    )

    m = MAPPLE(mkparams())
    @test !issorted(m.params.peaks.log_f)
    @test !issorted(m.params.components.log_f_stop)

    f = exp10.(range(0, 3, length = 100))
    s = mapple(f, m.params)

    # Non-mutating `sort` must not touch the original.
    m2 = sort(m)
    @test m2 isa MAPPLE
    @test !issorted(m.params.components.log_f_stop)

    # Reordering is a pure relabeling: it preserves the predicted spectrum and the
    # *set* of peaks/components, with each sub-component's fields kept together.
    @test predict(m2, f) ≈ s
    @test Set(collect(Float64, m2.params.peaks.log_A)) == Set([0.51, 1.52])
    @test Set(collect(Float64, m2.params.components.β)) == Set([2.0, 4.0])
    # the 0.7-peak must carry its own σ and A with it to the front
    i07 = findfirst(≈(0.7), collect(Float64, m2.params.peaks.log_f))
    @test collect(Float64, m2.params.peaks.log_σ)[i07] ≈ 0.12
    @test collect(Float64, m2.params.peaks.log_A)[i07] ≈ 1.52

    # When the model is already sorted, sorting is a no-op and is safe.
    msorted = MAPPLE(
        ComponentArray(;
            log_A = 1.0,
            peaks = [
                mkpeak(; log_f = 0.7, log_σ = 0.12, log_A = 1.52),
                mkpeak(; log_f = 2.5, log_σ = 0.11, log_A = 0.51),
            ],
            components = [
                mkcomp(; log_f_stop = 1.5, β = 2.0),
                mkcomp(; log_f_stop = 5.0, β = 4.0),
            ],
            transition_width = 0.2
        )
    )
    ssorted = mapple(f, msorted.params)
    sort!(msorted)
    @test mapple(f, msorted.params) ≈ ssorted
end

@testitem "mapple: fit recovers a known spectrum (peaks=0)" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

    truth = ComponentArray(;
        components = [mkcomp(; β = 2.0, log_f_stop = 1.5), mkcomp(; β = 4.0, log_f_stop = 5.0)],
        peaks = [
            mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5),
            mkpeak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5),
        ],
        transition_width = 0.15, log_A = 1.0
    )

    log_f = range(0, 3, length = 500)
    f = exp10.(log_f)
    Random.seed!(42)
    s_clean = mapple(f, truth)
    log_s = log10.(s_clean) .+ 0.1 .* randn(length(s_clean))

    # Rough fit (peak-finding + regression) then Optim refinement; components only.
    init = fit_mapple(log_f, log_s; components = 2, peaks = 0, w = 50)
    refined = fit_mapple(log_f, log_s, init; autodiff = Optim.ADTypes.AutoForwardDiff())

    loss(p) = sum((log_s .- log10.(mapple(f, p))) .^ 2)
    @test loss(refined) ≤ loss(init)            # refinement improves the fit
    @test loss(refined) ≤ 1.5 * loss(truth)     # lands near the optimum
    @test cor(mapple(f, refined), s_clean) > 0.99
    @test sort(collect(Float64, refined.components.β)) ≈ [2.0, 4.0] atol = 0.2
    @test refined.log_A ≈ 1.0 atol = 0.1

    MapplePlots.save_fit("fit_recovers_peaks0", f, exp10.(log_s), refined; init = init, subdir = "diagnostics")
end

@testitem "mapple: type-form entry points & empty-block construction" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    f = exp10.(range(0, 3, length = 200))
    p = ComponentArray(;
        log_A = 0.5, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    s = mapple(f, p)

    # `fit(MAPPLE, freqs, spectrum)` wraps the frequency vector in a frequency dimension
    # and fits.
    m = fit(MAPPLE, f, s; peaks = 0, components = 1)
    @test m isa MAPPLE
    full = fit!(deepcopy(m), ToolsArray(s, 𝑓(f)))   # refine the rough fit for the figure
    MapplePlots.save_fit("type_form_fit", f, s, full; init = m, subdir = "diagnostics")

    # `fit!(::Type{MAPPLE}, ...)` was a never-working misnomer for `fit`; it is removed,
    # so calling it raises a MethodError.
    @test_throws MethodError fit!(MAPPLE, f, s)

    # A zero-peak model builds via the `MAPPLE` constructor even from a plain empty
    # `Vector{<:ComponentArray}` (the constructor normalises empty blocks).
    PeakT = typeof(mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0))
    m0 = MAPPLE(;
        log_A = 0.5, peaks = Vector{PeakT}(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    @test m0 isa MAPPLE
    @test length(m0.params.peaks) == 0
    @test all(>(0), mapple(f, m0.params))
end

@testitem "mapple: peak finder misses peaks on a rising background" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    import TimeseriesTools: findpeaks, Log10𝑓

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

    log_f = range(0, 3, length = 500)
    f = exp10.(log_f)
    truth = ComponentArray(;
        log_A = 1.0,
        peaks = [
            mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5),
            mkpeak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5),
        ],
        components = [mkcomp(; log_f_stop = 1.5, β = 2.0), mkcomp(; log_f_stop = 5.0, β = 4.0)],
        transition_width = 0.15
    )
    log_s = log10.(mapple(f, truth))

    logspec = ToolsArray(log_s, Log10𝑓(collect(log_f)))
    minprom = (maximum(log_s) - minimum(log_s)) / 50
    _, proms, _ = findpeaks(logspec, 50; minprom)
    # Both true peaks sit on a steeply *rising* slope, so neither is a local maximum of
    # the raw spectrum: stage-1 detection finds fewer than the 2 true peaks. This is the
    # motivation for detrending the spectrum before peak finding.
    @test length(proms) < 2
end

@testitem "MAPPLE: typed accessors" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    m = MAPPLE(;
        log_A = 1.0,
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5)],
        components = [mkcomp(; log_f_stop = 1.5, β = 2.0), mkcomp(; log_f_stop = 5.0, β = 4.0)],
        transition_width = 0.2
    )
    @test betas(m) isa Vector{Float64}
    @test betas(m) == [2.0, 4.0]
    @test breakpoints(m) == [1.5, 5.0]
    @test peakfreqs(m) == [0.7]
    @test peaksigmas(m) == [0.1]
    @test peakamplitudes(m) == [1.5]
    # empty peaks must not error and stay typed
    PeakT = typeof(mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0))
    m0 = MAPPLE(;
        log_A = 1.0, peaks = Vector{PeakT}(),
        components = [mkcomp(; log_f_stop = 5.0, β = 2.0)], transition_width = 0.2
    )
    @test peakfreqs(m0) == Float64[]
    @test peakfreqs(m0) isa Vector{Float64}
end

@testitem "MAPPLE: diagnostics" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    truth = MAPPLE(;
        log_A = 1.0,
        peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    f = exp10.(range(0, 3, length = 300))
    spectrum = ToolsArray(mapple(f, truth.params), 𝑓(f))
    # the exact model on its own spectrum: zero residual, perfect fit
    @test maximum(abs, mapple_residuals(truth, spectrum)) < 1.0e-8
    @test mapple_loss(truth, spectrum) < 1.0e-12
    @test rsquared(truth, spectrum) ≈ 1.0 atol = 1.0e-8
    # a worse model has higher loss and lower R²
    wrong = MAPPLE(;
        log_A = 1.0,
        peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 5.0, β = -1.0)], transition_width = 0.1
    )
    @test mapple_loss(wrong, spectrum) > mapple_loss(truth, spectrum)
    @test rsquared(wrong, spectrum) < rsquared(truth, spectrum)
end

@testitem "MAPPLE: precomputed log_f fast path" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    p = ComponentArray(;
        log_A = 0.5, peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 2.0, β = -1.0), mkcomp(; log_f_stop = 5.0, β = -3.0)],
        transition_width = 0.1
    )
    f = exp10.(range(0, 3, length = 400))
    log_f = log10.(f)
    cps = p[[:log_A, :transition_width, :components]]
    pks = p[[:peaks]]
    s_ref = mapple(f, cps, pks)
    s_pre = similar(s_ref)
    TimeseriesTools.mapple!(s_pre, f, cps, pks; log_f)
    @test s_pre ≈ s_ref
    # passing log_f avoids allocating a fresh log10.(f) each call
    TimeseriesTools.mapple!(s_pre, f, cps, pks; log_f)   # warm up
    TimeseriesTools.mapple!(s_pre, f, cps, pks)
    a = @allocated TimeseriesTools.mapple!(s_pre, f, cps, pks; log_f)
    b = @allocated TimeseriesTools.mapple!(s_pre, f, cps, pks)
    @test a < b
end

@testitem "mapple: robust fit on a steep rising spectrum (regression)" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    # Steep rising background (β = 2, 4) with two peaks: the case that used to collapse to a
    # degenerate β / log_A optimum. The zero-peak detrend + single joint refine fits it cleanly.
    truth = ComponentArray(;
        components = [mkcomp(; β = 2.0, log_f_stop = 1.5), mkcomp(; β = 4.0, log_f_stop = 5.0)],
        peaks = [
            mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5),
            mkpeak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5),
        ],
        transition_width = 0.15, log_A = 1.0
    )
    log_f = range(0, 3, length = 500)
    f = exp10.(log_f)
    Random.seed!(0)
    s_clean = mapple(f, truth)
    log_s = log10.(s_clean) .+ 0.05 .* randn(length(s_clean))
    loss(p) = sum(abs2, log_s .- log10.(mapple(f, p)))
    init = fit_mapple(log_f, log_s; components = 2, peaks = 2, w = 50)
    refined = fit_mapple(log_f, log_s, init; autodiff = Optim.ADTypes.AutoForwardDiff())
    @test length(refined.peaks) == 2                 # peaks detected, not the fallback
    @test loss(refined) < 1.5 * loss(truth)          # no degenerate collapse
    @test cor(mapple(f, refined), s_clean) > 0.99
    @test sort(collect(Float64, refined.components.β)) ≈ [2.0, 4.0] atol = 0.5

    MapplePlots.save_fit("robust_steep_rising", f, exp10.(log_s), refined; init = init, subdir = "diagnostics")
end

@testitem "mapple: rough-fit base amplitude is not clobbered by the peak loop" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Random, Optim, ForwardDiff
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    # The peak-construction `do` block must use fresh names: assigning `log_A` inside it would
    # leak to the enclosing `log_A = first(log_s)` and corrupt the base amplitude to the last
    # peak's value (a catastrophic init that biases the subsequent refine).
    truth = ComponentArray(;
        components = [mkcomp(; β = 2.0, log_f_stop = 1.5), mkcomp(; β = 4.0, log_f_stop = 5.0)],
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5), mkpeak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5)],
        transition_width = 0.15, log_A = 1.0
    )
    log_f = range(0, 3, length = 500); f = exp10.(log_f)
    Random.seed!(0); log_s = log10.(mapple(f, truth)) .+ 0.05 .* randn(length(f))
    init = fit_mapple(log_f, log_s; components = 2, peaks = 2, w = 50)
    @test init.log_A ≈ first(log_s)              # base amplitude is the DC estimate, not a peak's
    @test all(isfinite, mapple(f, init))         # init is sane, not a 10^(peak log_A) blow-up

    refined = fit_mapple(log_f, log_s, init)     # full fit, for the figure
    MapplePlots.save_fit("rough_init_base_amplitude", f, exp10.(log_s), refined; init = init, subdir = "diagnostics")
end

@testitem "MAPPLE: Hz reporting accessors" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    m = MAPPLE(;
        log_A = 1.0, peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 0.5)],
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 2.5, β = -3.0)],
        transition_width = 0.1
    )
    # Physical-unit accessors convert OUT of log-10 space at the reporting boundary.
    @test peakcentres(m) ≈ exp10.(peakfreqs(m)) ≈ [10.0]
    @test breakfrequencies(m) ≈ exp10.(breakpoints(m)) ≈ [10.0^1.5, 10.0^2.5]
    @test peakheights(m) ≈ exp10.(peakamplitudes(m)) ≈ [10.0^0.5]
    # FWHM in Hz: σ = f·tanh(log_σ), FWHM = 2√(2ln2)·σ.
    @test peakbandwidths(m) ≈ [2 * sqrt(2 * log(2)) * 10.0 * tanh(0.2)]
    @test all(>(0), peakbandwidths(m))
    # Empty peaks stay typed and empty.
    PeakT = typeof(mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0))
    m0 = MAPPLE(;
        log_A = 1.0, peaks = Vector{PeakT}(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    @test peakcentres(m0) == Float64[] && peakbandwidths(m0) == Float64[] && peakheights(m0) == Float64[]
end

@testitem "MAPPLE: rsquared guards a flat spectrum" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    m = MAPPLE(;
        log_A = 1.0, peaks = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0),
        components = [mkcomp(; log_f_stop = 5.0, β = 0.0)], transition_width = 0.1
    )
    f = exp10.(range(0, 3, length = 100))
    flat = ToolsArray(fill(10.0, 100), 𝑓(f))   # constant spectrum => ss_tot = 0
    @test isnan(rsquared(m, flat))             # undefined, not ±Inf or a DivideError
end

@testitem "mapple: peak-width guard rejects out-of-range detections" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Random, Optim, ForwardDiff
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)
    # Falling 1-component background with one clear, narrow peak.
    truth = ComponentArray(;
        log_A = 1.0, peaks = [mkpeak(; log_f = 1.5, log_σ = 0.1, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    log_f = range(0, 3, length = 400); f = exp10.(log_f)
    Random.seed!(7); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))
    # Default limits: the narrow peak is found.
    init = fit_mapple(log_f, log_s; components = 1, w = 50)
    @test length(init.peaks) ≥ 1
    # Demanding very wide peaks rejects every detection before counting.
    @test length(fit_mapple(log_f, log_s; components = 1, w = 50, peak_width_limits = (2.0, 3.0)).peaks) == 0

    # Refine for a clean full fit: the recovered peak lands at the true centre and height.
    refined = fit_mapple(log_f, log_s, init)
    m = MAPPLE(refined)
    @test only(peakcentres(m)) ≈ 10.0^1.5 rtol = 0.05
    @test only(peakamplitudes(m)) ≈ 1.0 atol = 0.15

    MapplePlots.save_fit("peak_width_guard", f, exp10.(log_s), refined; init = init, subdir = "low")
end

@testitem "mapple: :auto detects a prominent peak and rejects noise" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Random, Optim, ForwardDiff
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    truth = ComponentArray(;
        log_A = 1.0, peaks = [mkpeak(; log_f = 1.5, log_σ = 0.12, log_A = 1.2)],
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    log_f = range(0, 3, length = 400); f = exp10.(log_f)
    Random.seed!(11); log_s = log10.(mapple(f, truth)) .+ 0.03 .* randn(length(f))
    # Rough init: `:auto` finds the one real peak and seeds its amplitude from the LOCAL background
    # under it, so the start is already near the true height (anchoring to the global spectral floor
    # instead would start it orders of magnitude too low for the refine to grow — see the peak-init
    # note in `fit_mapple`).
    init = fit_mapple(log_f, log_s; components = 1, w = 50)   # peaks = :auto
    @test length(init.peaks) == 1                              # the one real peak, not noise bumps
    @test only(collect(Float64, init.peaks.log_f)) ≈ 1.5 atol = 0.2
    @test only(collect(Float64, init.peaks.log_A)) ≈ 1.2 atol = 0.6   # seeded near the true height

    # The Optim refinement then polishes the centre and height to the data.
    refined = fit_mapple(log_f, log_s, init)
    m = MAPPLE(refined)
    @test only(peakcentres(m)) ≈ 10.0^1.5 rtol = 0.05
    @test only(peakamplitudes(m)) ≈ 1.2 atol = 0.15
    @test only(peakheights(m)) ≈ 10.0^1.2 rtol = 0.25

    MapplePlots.save_fit("auto_detect_peak", f, exp10.(log_s), refined; init = init, subdir = "low")
end

@testitem "mapple: BIC selects the simplest adequate component count" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Random, Optim, ForwardDiff
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)
    # A genuine single power law: the free-knot BIC penalty must not over-select components.
    truth = ComponentArray(;
        log_A = 1.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    f = exp10.(range(0, 3, length = 150))
    Random.seed!(5); s = exp10.(log10.(mapple(f, truth)) .+ 0.03 .* randn(length(f)))
    spec = ToolsArray(s, 𝑓(f))
    m = fit(MAPPLE, spec)                       # :auto BIC selection (refined when Optim loaded)
    @test length(m.params.components) == 1
    @test rsquared(m, spec) > 0.99

    init = fit_mapple(log10.(f), log10.(s); components = length(m.params.components))
    MapplePlots.save_fit("bic_component_selection", f, s, m; init = init, subdir = "diagnostics")
end

@testitem "mapple: fit on a flat / low-dynamic-range spectrum does not crash (regression)" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)
    f = exp10.(range(0, 3, length = 200))
    # A perfectly flat spectrum at a high baseline. The log_A upper bound used to be
    # range-relative (100*(max-min log_s)), collapsing to 0 while the init log_A = first(log_s) = 2
    # sat above it, so Fminbox rejected the start. The bound is now absolute, so fit and fit! work.
    flat = ToolsArray(fill(100.0, 200), 𝑓(f))
    m = fit(MAPPLE, flat)
    @test m isa MAPPLE
    @test all(>(0), predict(m, f))
    m2 = MAPPLE(;
        log_A = 2.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = 0.0)], transition_width = 0.1
    )
    @test fit!(m2, flat) isa MAPPLE
    @test m2.params.log_A ≈ 2.0 atol = 0.5      # base amplitude near the flat level (10^2)

    init = fit_mapple(log10.(f), log10.(parent(flat)); components = 1)
    MapplePlots.save_fit("flat_spectrum_fit", f, parent(flat), m; init = init, subdir = "diagnostics")
end

@testitem "mapple: fit! is robust to degenerate spectra (finite + bounded)" tags = [:mapple] begin
    using TimeseriesTools, Optim, ForwardDiff

    f = exp10.(range(0, 3; length = 200))
    fitβ(x; kw...) = (m = fit(MAPPLE, x; components = 2, peaks = 0); fit!(m, x; kw...); first(m.params.components.β))
    @test isfinite(fitβ(ToolsArray(f .^ -1.5, 𝑓(f))))   # sanity: a clean power law fits

    # (A) A zero/silent bin → `log10(0) = -Inf` in the objective must not poison the whole fit into
    # NaN/Inf. Regression: before the safe-log10 guard, one zero bin drove the fitted exponent to NaN.
    zb = collect(f .^ -1.5); zb[80:90] .= 0.0
    @test isfinite(fitβ(ToolsArray(zb, 𝑓(f))))

    # (B) A real low-dynamic-range MAD curve (a WRCircuit working-regime `itot` input) on which the
    # UNCAPPED Fminbox refine wandered ~30-60 s even though the optimum is reached early. The bounded
    # default budget must reach the SAME exponent in a small fraction of that time. Regression: with no
    # default iteration cap, curves like this turned a per-neuron parameter sweep into an effective hang.
    lags = unique(round.(Int, exp10.(range(1, 4, 100)))) ./ 10
    s = [1.45763, 1.57037, 1.6742, 1.7685, 1.85299, 1.92805, 1.9947, 2.05432, 2.15678, 2.20155,
        2.2818, 2.31816, 2.3851, 2.44564, 2.47383, 2.55167, 2.59857, 2.64202, 2.70158, 2.73814,
        2.77541, 2.69008, 2.41115, 2.12487, 1.80093, 1.78285, 2.01123, 2.25617, 2.42017, 2.51865,
        2.58963, 2.63558, 2.64861, 2.40036, 1.98147, 2.02651, 2.31385, 2.54985, 2.69065, 2.59615,
        2.25737, 2.3472, 2.53125, 2.60808, 2.34354, 2.47621, 2.70026, 2.59311, 2.59517, 2.72218,
        2.58298, 2.69988, 2.61941, 2.70915, 2.71232, 2.78858, 2.80357, 2.80875, 2.87425, 2.8744,
        2.91203, 2.93889, 2.94648, 3.0124, 3.03716, 3.08576, 3.12241, 3.15742, 3.1982, 3.25253,
        3.26843, 3.28528, 3.31252, 3.31327, 3.37337, 3.37107, 3.3838, 3.40228, 3.40549, 3.39346,
        3.37416, 3.32477, 3.29771, 3.28583, 3.29731, 3.32491, 3.36026, 3.37732, 3.39327, 3.42337,
        3.4288, 3.44047, 3.47026, 3.42306, 3.4001, 3.34346, 3.29111, 3.29352, 3.32632]
    hard = ToolsArray(s, 𝑓(Float64.(lags)))
    fitβ(hard); fitβ(hard; iterations = 5000, outer_iterations = 50)        # warmup both budgets
    β_cap = fitβ(hard)                                   # default (bounded budget)
    t_cap = @elapsed fitβ(hard)
    β_big = fitβ(hard; iterations = 5000, outer_iterations = 50)            # an effectively unbounded budget
    @test β_cap ≈ β_big atol = 0.05                      # the cap reaches the same exponent...
    @test t_cap < 10                                     # ...far faster (uncapped was ~tens of seconds)
end

@testitem "mapple: refine options route through the `refine` channel" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)
    f = exp10.(range(0, 3, length = 200))
    truth = ComponentArray(;
        log_A = 1.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    Random.seed!(3); s = exp10.(log10.(mapple(f, truth)) .+ 0.05 .* randn(length(f)))
    spec = ToolsArray(s, 𝑓(f))
    # Optim refine options (iterations, multistart) pass through `refine` without reaching the
    # peak finder. Previously these leaked into findpeaks and raised a MethodError.
    @test fit(MAPPLE, spec; refine = (; iterations = 50)) isa MAPPLE
    m = fit(MAPPLE, spec; refine = (; multistart = 2))
    @test m isa MAPPLE

    init = fit_mapple(log10.(f), log10.(s); components = 1)
    MapplePlots.save_fit("refine_channel", f, s, m; init = init, subdir = "diagnostics")
end

@testitem "mapple: degenerate fits raise a clear ArgumentError" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)
    f = exp10.(range(0, 3, length = 100)); s = exp10.((-2.0) .* log10.(f))
    spec = ToolsArray(s, 𝑓(f))
    # A non-positive component count is rejected up front, not via a cryptic range error.
    @test_throws ArgumentError fit(MAPPLE, spec; components = 0)
    @test_throws ArgumentError fit_mapple(log10.(f), log10.(s); components = 0)
    # Fewer than two frequencies (or a zero-span grid) cannot be fit.
    @test_throws ArgumentError fit_mapple([1.0], [0.5]; components = 1)
    @test_throws ArgumentError fit(MAPPLE, [10.0], [3.0])
    # `peak_width_limits` must be ordered (wmin ≤ wmax). A swapped tuple is caught up front rather
    # than silently rejecting every detection.
    @test_throws ArgumentError fit_mapple(
        log10.(f), log10.(s); components = 1, peak_width_limits = (2.0, 0.5)
    )
    # The MAPPLE constructor requires at least one component.
    @test_throws ArgumentError MAPPLE(;
        log_A = 1.0, peaks = nopeaks(),
        components = ComponentArray[], transition_width = 0.1
    )
end

@testitem "MAPPLE: show reports Hz consistently with the accessors" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    m = sort(
        MAPPLE(;
            log_A = 1.0, peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 0.5)],
            components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 2.5, β = -3.0)],
            transition_width = 0.1
        )
    )
    str = sprint(show, m)
    @test occursin("Hz", str) && occursin("FWHM", str)
    # The figures printed by `show` are exactly the accessor values (single source of truth).
    @test occursin("$(round(peakcentres(m)[1], sigdigits = 4))", str)
    @test occursin("$(round(breakfrequencies(m)[1], sigdigits = 4))", str)
    @test occursin("$(round(peakbandwidths(m)[1], sigdigits = 4))", str)
end

@testitem "mapple: multistart refinement is robust and never worsens the fit" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    truth = ComponentArray(;
        components = [mkcomp(; β = 2.0, log_f_stop = 1.5), mkcomp(; β = 4.0, log_f_stop = 5.0)],
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5), mkpeak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5)],
        transition_width = 0.15, log_A = 1.0
    )
    log_f = range(0, 3, length = 500); f = exp10.(log_f)
    Random.seed!(0); log_s = log10.(mapple(f, truth)) .+ 0.05 .* randn(length(f))
    loss(p) = sum(abs2, log_s .- log10.(mapple(f, p)))
    init = fit_mapple(log_f, log_s; components = 2, peaks = 2, w = 50)
    single = fit_mapple(log_f, log_s, init)
    Random.seed!(1); multi = fit_mapple(log_f, log_s, init; multistart = 4)
    @test multi isa ComponentArray
    @test loss(multi) ≤ loss(single) + 1.0e-8     # restarts never worsen the unperturbed result
    # Multistart survives a degenerate over-parameterised noiseless target without crashing
    # (a perturbed restart where Fminbox cannot build a barrier is skipped, not fatal).
    clean = log10.(mapple(f, truth))
    cinit = fit_mapple(log_f, clean; components = 3, peaks = 2, w = 50)
    Random.seed!(2); @test fit_mapple(log_f, clean, cinit; multistart = 3) isa ComponentArray

    MapplePlots.save_fit("multistart_refine", f, exp10.(log_s), multi; init = init, subdir = "diagnostics")
end

# --- Complex-spectrum accuracy ----------------------------------------------------------------
# These items exercise the refinement on harder spectra (multiple peaks, >2 power-law segments).
# They judge accuracy in LOG-10 space — the space the loss is defined in — rather than by linear
# `cor`, which is dominated by the high-power low-frequency samples and stays ≈1 even when a
# mid-band peak is missed entirely. The conservative rough init only seeds peaks and slopes; the
# Optim refine is what grows peaks to their true height, so every accuracy claim is on `refined`.

@testitem "mapple: refined fit recovers peak heights on a multi-peak, multi-component spectrum" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

    # Three broken-power-law segments with two well-separated peaks (one in the first segment,
    # one in the last). This is the case the misleading `auto_detect_peak` figure failed to show:
    # the rough init's peaks start near the spectral floor; only after refinement do they reach
    # the data height.
    truth = ComponentArray(;
        log_A = 2.0,
        peaks = [
            mkpeak(; log_f = 0.5, log_σ = 0.1, log_A = 1.5),
            mkpeak(; log_f = 2.5, log_σ = 0.12, log_A = 1.0),
        ],
        components = [
            mkcomp(; log_f_stop = 1.2, β = -1.0), mkcomp(; log_f_stop = 2.0, β = -2.5),
            mkcomp(; log_f_stop = 5.0, β = -3.5),
        ], transition_width = 0.08
    )
    log_f = range(0, 3, length = 600); f = exp10.(log_f)
    Random.seed!(2); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))

    init = fit_mapple(log_f, log_s; components = 3, peaks = 2, w = 40)
    refined = fit_mapple(log_f, log_s, init)
    m = MAPPLE(refined)

    # Both peaks are kept (neither driven to log_A → −∞ and dropped).
    @test length(refined.peaks) == 2
    # Centres and heights recovered. Heights are the crux: the rough init starts them an order of
    # magnitude low, so this fails unless the refine genuinely grows the peaks.
    @test sort(peakfreqs(m)) ≈ [0.5, 2.5] atol = 0.15
    @test sort(peakamplitudes(m)) ≈ [1.0, 1.5] atol = 0.2
    @test sort(peakheights(m)) ≈ exp10.([1.0, 1.5]) rtol = 0.3

    # The whole spectrum is reproduced in log space (not just correlated in linear space).
    logfit = log10.(mapple(f, refined)); logtruth = log10.(mapple(f, truth))
    logr2 = 1 - sum(abs2, logtruth .- logfit) / sum(abs2, logtruth .- mean(logtruth))
    @test logr2 > 0.99
    @test maximum(abs, logtruth .- logfit) < 0.15

    MapplePlots.save_fit("complex_multipeak", f, exp10.(log_s), refined; init = init, subdir = "med")
end

@testitem "mapple: >2 components — three power-law segments recover their slopes" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    # A genuine three-segment background (steepening: −1, −2.5, −4). With no peaks to trade against,
    # the slopes are identifiable, so we can pin them down directly.
    truth = ComponentArray(;
        log_A = 1.0, peaks = nopeaks(),
        components = [
            mkcomp(; log_f_stop = 1.0, β = -1.0), mkcomp(; log_f_stop = 2.0, β = -2.5),
            mkcomp(; log_f_stop = 5.0, β = -4.0),
        ], transition_width = 0.1
    )
    log_f = range(0, 3, length = 500); f = exp10.(log_f)
    Random.seed!(1); log_s = log10.(mapple(f, truth)) .+ 0.03 .* randn(length(f))

    init = fit_mapple(log_f, log_s; components = 3, peaks = 0, w = 50)
    refined = fit_mapple(log_f, log_s, init)
    m = MAPPLE(refined)

    @test length(refined.components) == 3
    @test sort(betas(m)) ≈ [-4.0, -2.5, -1.0] atol = 0.2
    # `n` components carry `n - 1` breakpoints: the last `log_f_stop` is inert (the final component
    # is never windowed), so only the two interior knots are asserted.
    @test sort(breakpoints(m))[1:2] ≈ [1.0, 2.0] atol = 0.2

    logfit = log10.(mapple(f, refined)); logtruth = log10.(mapple(f, truth))
    logr2 = 1 - sum(abs2, logtruth .- logfit) / sum(abs2, logtruth .- mean(logtruth))
    @test logr2 > 0.99

    MapplePlots.save_fit("three_component_slopes", f, exp10.(log_s), refined; init = init, subdir = "med")
end

@testitem "mapple: four-component spectrum is reproduced accurately" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    # Four segments over four decades. The individual (β, breakpoint) values are NOT uniquely
    # identifiable here — several decompositions yield near-identical curves — so we assert the
    # spectrum is reproduced, not the parameters. This is the honest claim for deep stacks of
    # power laws: the model is expressive enough to capture the shape.
    truth = ComponentArray(;
        log_A = 2.0, peaks = nopeaks(),
        components = [
            mkcomp(; log_f_stop = 0.5, β = -0.5), mkcomp(; log_f_stop = 1.5, β = -1.5),
            mkcomp(; log_f_stop = 2.5, β = -2.5), mkcomp(; log_f_stop = 6.0, β = -3.5),
        ],
        transition_width = 0.07
    )
    log_f = range(-0.5, 3.5, length = 700); f = exp10.(log_f)
    Random.seed!(11); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))

    init = fit_mapple(log_f, log_s; components = 4, peaks = 0, w = 40)
    refined = fit_mapple(log_f, log_s, init)

    @test length(refined.components) == 4
    logfit = log10.(mapple(f, refined)); logtruth = log10.(mapple(f, truth))
    logr2 = 1 - sum(abs2, logtruth .- logfit) / sum(abs2, logtruth .- mean(logtruth))
    @test logr2 > 0.99
    @test maximum(abs, logtruth .- logfit) < 0.1

    MapplePlots.save_fit("four_component_spectrum", f, exp10.(log_s), refined; init = init, subdir = "high")
end

# --- Graduated complexity sweep ---------------------------------------------------------------
# A staircase of increasingly complex spectra (more breaks, more peaks) so the figures show where
# the fit starts to break. Figures are bucketed into `mapple_figs/{low,med,high,very_high}/` for
# navigation; the recovery tests above also file into these tiers (auto_detect_peak/peak_width_guard
# → low, complex_multipeak/three_component_slopes → med, four_component_spectrum → high).
#
# With peaks seeded from the LOCAL background (see the peak-init note in `fit_mapple`), peak recovery
# is now robust well past a handful of peaks: the low/med/high tiers recover every peak. The frontier
# has moved to the BACKGROUND — at ~6 components the broken-power-law has too many near-degenerate
# (β, breakpoint) combinations, so the optimizer wanders (occasionally to a NaN gradient, which Optim
# absorbs) and the recovered slopes/knots drift even though the peaks are still found. The `very_high`
# tier exercises that regime; we assert only the robustness invariants there (finite, positive
# spectrum; a loose log-R² floor) and record the worst-case residual via `@info`, with the figure
# carrying the rest of the story. Lower tiers additionally assert that all peaks are recovered.

@testitem "mapple: complexity sweep (low) — simple spectra fit cleanly" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    log_f = range(0, 3.5, length = 700); f = exp10.(log_f)
    logr2(refined, truth) = (
        lp = log10.(mapple(f, refined)); lt = log10.(mapple(f, truth));
        1 - sum(abs2, lt .- lp) / sum(abs2, lt .- mean(lt))
    )

    # Two segments, one peak: firmly inside the regime the refine handles well.
    truth = ComponentArray(;
        log_A = 2.0, peaks = [mkpeak(; log_f = 0.6, log_σ = 0.1, log_A = 1.3)],
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 10.0, β = -2.8)],
        transition_width = 0.06
    )
    Random.seed!(1); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))
    init = fit_mapple(log_f, log_s; components = 2, peaks = 1, w = 40)
    refined = fit_mapple(log_f, log_s, init)
    sfit = mapple(f, refined)
    @test all(isfinite, sfit) && all(>(0), sfit)
    @test logr2(refined, truth) > 0.99

    MapplePlots.save_fit("two_components_one_peak", f, exp10.(log_s), refined; init = init, subdir = "low")
end

@testitem "mapple: complexity sweep (medium) — three peaks on three segments" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    log_f = range(0, 3.5, length = 700); f = exp10.(log_f)
    logr2(refined, truth) = (
        lp = log10.(mapple(f, refined)); lt = log10.(mapple(f, truth));
        1 - sum(abs2, lt .- lp) / sum(abs2, lt .- mean(lt))
    )

    # Three segments, three peaks (one per segment). All three are recovered — both centres and
    # heights — thanks to the local-background amplitude seeding.
    truth = ComponentArray(;
        log_A = 2.0,
        peaks = [
            mkpeak(; log_f = 0.45, log_σ = 0.09, log_A = 1.2),
            mkpeak(; log_f = 1.6, log_σ = 0.09, log_A = 1.0),
            mkpeak(; log_f = 2.75, log_σ = 0.09, log_A = 1.0),
        ],
        components = [
            mkcomp(; log_f_stop = 1.0, β = -1.0), mkcomp(; log_f_stop = 2.2, β = -2.0),
            mkcomp(; log_f_stop = 10.0, β = -3.2),
        ], transition_width = 0.06
    )
    Random.seed!(1); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))
    init = fit_mapple(log_f, log_s; components = 3, peaks = 3, w = 40)
    refined = fit_mapple(log_f, log_s, init)
    sfit = mapple(f, refined)
    m = MAPPLE(refined)
    @test all(isfinite, sfit) && all(>(0), sfit)
    @test logr2(refined, truth) > 0.99

    # All three peaks are kept (none collapsed to log_A → −∞) and land at their true centres/heights.
    @test count(>(-3.0), collect(Float64, refined.peaks.log_A)) == 3
    @test sort(peakfreqs(m)) ≈ [0.45, 1.6, 2.75] atol = 0.15
    @test sort(peakamplitudes(m)) ≈ [1.0, 1.0, 1.2] atol = 0.2

    @info "medium-complexity fit" peaks_recovered = "3/3" maxres = round(maximum(abs, log10.(mapple(f, truth)) .- log10.(sfit)), digits = 3)

    MapplePlots.save_fit("three_components_three_peaks", f, exp10.(log_s), refined; init = init, subdir = "med")
end

@testitem "mapple: complexity sweep (high) — dense peaks on four segments" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    log_f = range(0, 3.5, length = 700); f = exp10.(log_f)
    logr2(refined, truth) = (
        lp = log10.(mapple(f, refined)); lt = log10.(mapple(f, truth));
        1 - sum(abs2, lt .- lp) / sum(abs2, lt .- mean(lt))
    )

    # Four segments with three then four peaks. Even on this dense field every peak is recovered at
    # its true centre — the boxed refine (one of the two bound sets `fit_mapple` tries) stops the
    # peaks from sliding together or one ballooning into a background-like blob, which is what used
    # to swallow the smaller peaks next to large ones.
    components = [
        mkcomp(; log_f_stop = 0.8, β = -0.8), mkcomp(; log_f_stop = 1.7, β = -1.6),
        mkcomp(; log_f_stop = 2.6, β = -2.6), mkcomp(; log_f_stop = 10.0, β = -3.6),
    ]
    configs = [
        (
            "four_components_three_peaks", 3,
            [
                mkpeak(; log_f = 0.35, log_σ = 0.08, log_A = 1.2),
                mkpeak(; log_f = 1.25, log_σ = 0.08, log_A = 1.1),
                mkpeak(; log_f = 2.2, log_σ = 0.08, log_A = 1.0),
            ],
        ),
        (
            "four_components_four_peaks", 4,
            [
                mkpeak(; log_f = 0.35, log_σ = 0.08, log_A = 1.2),
                mkpeak(; log_f = 1.2, log_σ = 0.08, log_A = 1.1),
                mkpeak(; log_f = 2.05, log_σ = 0.08, log_A = 1.0),
                mkpeak(; log_f = 2.95, log_σ = 0.08, log_A = 1.0),
            ],
        ),
    ]
    for (name, np, pks) in configs
        truth = ComponentArray(; log_A = 2.5, peaks = pks, components = components, transition_width = 0.05)
        Random.seed!(1); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))
        init = fit_mapple(log_f, log_s; components = 4, peaks = np, w = 40)
        refined = fit_mapple(log_f, log_s, init)
        sfit = mapple(f, refined)
        m = MAPPLE(refined)
        @test all(isfinite, sfit) && all(>(0), sfit)   # never returns a non-physical spectrum
        @test logr2(refined, truth) > 0.99             # background and peaks both captured
        @test count(>(-3.0), collect(Float64, refined.peaks.log_A)) == np   # every peak kept
        @test sort(peakfreqs(m)) ≈ sort([p.log_f for p in pks]) atol = 0.15   # every peak at its centre

        @info "high-complexity fit" name peaks_recovered = "$np/$np" maxres = round(maximum(abs, log10.(mapple(f, truth)) .- log10.(sfit)), digits = 3)

        MapplePlots.save_fit(name, f, exp10.(log_s), refined; init = init, subdir = "high")
    end
end

@testitem "mapple: complexity sweep (very high) — background identifiability breaks down" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    lo, hi = -0.5, 4.5
    log_f = range(lo, hi, length = 1000); f = exp10.(log_f)
    logr2(refined, truth) = (
        lp = log10.(mapple(f, refined)); lt = log10.(mapple(f, truth));
        1 - sum(abs2, lt .- lp) / sum(abs2, lt .- mean(lt))
    )

    # Six segments, six peaks over five decades — the regime where the broken-power-law background
    # becomes non-identifiable: many (β, breakpoint) sets give near-identical curves, the optimizer
    # wanders (Optim absorbs the occasional NaN gradient), and the recovered slopes/knots drift even
    # though the peaks are still found. We assert only that the result stays physical and broadly
    # tracks the data; the figure shows the local background mismatches that define this frontier.
    edges = collect(range(lo, hi, length = 7)); stops = edges[2:end]; stops[end] = hi + 1.0
    betas = collect(range(-0.8, -0.8 - 0.7 * 5, length = 6))
    components = [mkcomp(; log_f_stop = stops[i], β = betas[i]) for i in 1:6]
    centres = collect(range(lo + 0.6, hi - 0.6, length = 6)) .+ 0.17
    peaks = [mkpeak(; log_f = centres[i], log_σ = 0.08, log_A = 1.1) for i in 1:6]
    truth = ComponentArray(; log_A = 2.5, peaks = peaks, components = components, transition_width = 0.05)
    Random.seed!(1); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))

    init = fit_mapple(log_f, log_s; components = 6, peaks = 6, w = 40)
    refined = fit_mapple(log_f, log_s, init)
    sfit = mapple(f, refined)
    @test all(isfinite, sfit) && all(>(0), sfit)   # still physical despite the degeneracy
    @test logr2(refined, truth) > 0.9              # broadly tracks the data, but degraded vs lower tiers

    alive = count(>(-3.0), collect(Float64, refined.peaks.log_A))
    @info "very-high-complexity fit" peaks_recovered = "$alive/6" maxres = round(maximum(abs, log10.(mapple(f, truth)) .- log10.(sfit)), digits = 3)

    MapplePlots.save_fit("six_components_six_peaks", f, exp10.(log_s), refined; init = init, subdir = "very_high")
end

# --- Noise-robustness sweeps ------------------------------------------------------------------
# How well does the fit hold up as measurement noise grows? Two noise models are exercised on a
# fixed, moderately-complex truth (two power-law segments + one peak):
#
#   1. UNIFORM (homoscedastic): constant-variance Gaussian noise in log-10 space, the same at
#      every frequency. The fit minimises an *unweighted* sum of squared log-residuals
#      (`OptimExt.mapple_loss`), so this is the noise model the loss is matched to. Because the
#      fitted curve is smooth and least squares is unbiased, the noise averages out over the grid
#      and recovery of the clean truth degrades only gently with σ.
#
#   2. FREQUENCY-INCREASING (heteroscedastic): the noise standard deviation ramps up with
#      frequency, so the high-frequency tail is far noisier than the low-frequency body — the
#      realistic shape for an un-averaged periodogram approaching its noise floor. The unweighted
#      loss does NOT down-weight the noisy tail, so the test checks that the (clean) low-frequency
#      structure is still recovered while the residual carries the injected frequency structure.
#
# Both items save one log–log figure per noise level, bucketed into `mapple_figs/noise/{uniform,
# increasing}/` with σ-sorted filenames, matching the per-tier figure convention used above.

@testitem "mapple: robustness to uniform spectral noise" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    logr2(fit, truth) = (
        lp = log10.(mapple(f, fit)); lt = log10.(mapple(f, truth));
        1 - sum(abs2, lt .- lp) / sum(abs2, lt .- mean(lt))
    )

    # Two falling segments with a single mid-band peak: comfortably inside the regime the refine
    # handles, so any loss of accuracy is attributable to the added noise, not the model complexity.
    truth = ComponentArray(;
        log_A = 2.0, peaks = [mkpeak(; log_f = 1.2, log_σ = 0.1, log_A = 1.2)],
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 10.0, β = -2.5)],
        transition_width = 0.08
    )
    log_f = range(0, 3, length = 600); f = exp10.(log_f)
    s_clean = mapple(f, truth)

    # A staircase of log-space noise amplitudes, from barely-perturbed to severe (σ = 0.4 ⇒ a
    # typical point sits a factor ~10^0.4 ≈ 2.5× off the true spectrum).
    sigmas = [0.02, 0.05, 0.1, 0.2, 0.4]
    scores = Float64[]
    for (i, σ) in enumerate(sigmas)
        Random.seed!(100 + i)
        log_s = log10.(s_clean) .+ σ .* randn(length(f))
        init = fit_mapple(log_f, log_s; components = 2, peaks = 1, w = 50)
        refined = fit_mapple(log_f, log_s, init)
        sfit = mapple(f, refined)
        @test all(isfinite, sfit) && all(>(0), sfit)        # a physical spectrum at every level
        push!(scores, logr2(refined, truth))
        MapplePlots.save_fit(
            "uniform_s$(lpad(round(Int, σ * 100), 2, '0'))", f, exp10.(log_s), refined;
            init = init, subdir = "noise/uniform"
        )
    end

    @test scores[1] > 0.99                                   # low noise: the clean truth is recovered
    @test all(>(0.8), scores)                               # never collapses, even at σ = 0.4
    # Graceful degradation: the heaviest two levels fit the clean truth no better, on average, than
    # the lightest level. (A trend, not strict per-step monotonicity, since each level is one draw.)
    @test mean(scores[(end - 1):end]) < scores[1]

    @info "uniform-noise sweep" sigmas = sigmas logr2 = round.(scores, digits = 4)
end

@testitem "mapple: robustness to frequency-increasing noise variance" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

    truth = ComponentArray(;
        log_A = 2.0, peaks = [mkpeak(; log_f = 1.2, log_σ = 0.1, log_A = 1.2)],
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 10.0, β = -2.5)],
        transition_width = 0.08
    )
    log_f = range(0, 3, length = 600); f = exp10.(log_f)
    s_clean = mapple(f, truth)
    n = length(f); half = n ÷ 2

    # Coefficient of determination of the whole fit against the CLEAN truth, in log-10 space.
    logr2(fit) = (
        lp = log10.(mapple(f, fit)); lt = log10.(s_clean);
        1 - sum(abs2, lt .- lp) / sum(abs2, lt .- mean(lt))
    )

    # σ ramps linearly with frequency from a small floor to a growing ceiling: the noise variance
    # is concentrated in the high-frequency tail. Three ceilings, from mild to severe.
    tail_res = Float64[]
    for (σ_lo, σ_hi) in [(0.02, 0.1), (0.02, 0.2), (0.02, 0.35)]
        Random.seed!(round(Int, σ_hi * 100))
        ramp = collect(range(σ_lo, σ_hi, length = n))
        log_s = log10.(s_clean) .+ ramp .* randn(n)
        init = fit_mapple(log_f, log_s; components = 2, peaks = 1, w = 50)
        refined = fit_mapple(log_f, log_s, init)
        sfit = mapple(f, refined)
        resid = log_s .- log10.(sfit)
        lo, hi = mean(abs, resid[1:half]), mean(abs, resid[(half + 1):end])

        @test all(isfinite, sfit) && all(>(0), sfit)            # physical despite the noisy tail
        # The injected variance — and so the fit residual — is concentrated at high frequency: the
        # heteroscedastic structure survives the fit rather than being smeared uniformly.
        @test hi > lo
        # Yet the unweighted fit still tracks the overall clean spectral shape: it is not dragged
        # into a wild curve by the noisy tail (this is the robustness claim the figures illustrate;
        # `logsample`ing the spectrum first sharpens it further — see the `fit(MAPPLE, …)` docstring).
        @test logr2(refined) > 0.9
        push!(tail_res, hi)

        MapplePlots.save_fit(
            "increasing_hi$(lpad(round(Int, σ_hi * 100), 2, '0'))", f, exp10.(log_s), refined;
            init = init, subdir = "noise/increasing"
        )
    end
    # Heavier tail noise leaves a larger high-frequency residual: degradation is graceful and
    # localised to the tail, not a global blow-up.
    @test tail_res[end] > tail_res[1]

    @info "frequency-increasing-noise sweep" tail_residual = round.(tail_res, digits = 4)
end

@testitem "mapple: degenerate peak parameters stay finite" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

    f = exp10.(range(0, 3, length = 200))

    # A hand-built model can carry a degenerate peak width log_σ ≤ 0 (the fit bounds keep log_σ > 0,
    # but nothing stops a user constructing one). Evaluation floors σ at a tiny positive fraction of
    # the centre, so the Gaussian never divides by zero / produces NaN at the peak centre.
    for bad_log_σ in (0.0, -0.5)
        p = ComponentArray(;
            log_A = 0.5, peaks = [mkpeak(; log_f = 1.0, log_σ = bad_log_σ, log_A = 1.0)],
            components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
        )
        s = mapple(f, p)
        @test all(isfinite, s)
        @test all(>(0), s)
    end

    # A low-prominence peak seeds its amplitude through the conservative init
    # `log_A = bg + log10(max(10^prom - 1, eps()))`; the `eps()` floor keeps log_A finite as the
    # prominence → 0, so the rough fit never returns NaN/Inf parameters on a faint bump.
    log_f = log10.(f)
    truth = ComponentArray(;
        log_A = 0.5, peaks = [mkpeak(; log_f = 1.5, log_σ = 0.08, log_A = -1.7)],
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1
    )
    log_s = log10.(mapple(f, truth))
    init = fit_mapple(log_f, log_s; components = 1, w = 50, peak_threshold = 1.0)
    @test all(isfinite, init)
end

# --- Outer-edge windowing ----------------------------------------------------------------------
# Regression cover for a windowing bug in `mapple!`: every component was multiplied by BOTH a
# `tanh` shoulder opening at the previous knot and one closing at its own, including the outermost
# components, which have nothing beyond them to crossfade into. The closing shoulder on the LAST
# component attenuated the model across the top of its own domain --- to exactly half at that knot
# --- so `last(β)` stopped being a slope and became a nuisance parameter absorbing the fade. Fits
# whose top segment ended inside the data compensated with a runaway β (+2.5, +8.5 observed) or
# collapsed the final knot onto its neighbour to switch the segment off entirely.
#
# The suite missed it because every multi-component case placed the truth model's final
# `log_f_stop` well beyond the sampled grid, where the fade falls outside the data.

@testitem "mapple: outer components are unwindowed at the edges of the model" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    f = exp10.(range(0, 3, length = 400))

    # Equal slopes across every component must reproduce a SINGLE power law exactly, wherever the
    # knots sit: the amplitudes chain to a common value, and unwindowed outer edges make the
    # shoulders a partition of unity (tanh is odd, so the two halves of a crossfade sum to 1).
    # A closing shoulder on the last component breaks this by tapering the top toward zero.
    for stop in (1.5, 2.0, 3.0)      # inside the grid, at its top edge, and at the last sample
        p = ComponentArray(;
            log_A = 0.7, peaks = nopeaks(),
            components = [mkcomp(; log_f_stop = 1.0, β = -1.3), mkcomp(; log_f_stop = stop, β = -1.3)],
            transition_width = 0.1
        )
        @test mapple(f, p) ≈ exp10(0.7) .* f .^ -1.3 rtol = 1.0e-10
    end

    # Same invariant with three components and well-separated knots. Interior crossfades only
    # partition to unity when they do not overlap, so this holds to a tolerance rather than exactly.
    p3 = ComponentArray(;
        log_A = 0.7, peaks = nopeaks(),
        components = [
            mkcomp(; log_f_stop = 1.0, β = -1.3), mkcomp(; log_f_stop = 2.0, β = -1.3),
            mkcomp(; log_f_stop = 2.5, β = -1.3),
        ], transition_width = 0.05
    )
    @test mapple(f, p3) ≈ exp10(0.7) .* f .^ -1.3 rtol = 1.0e-6

    # The specific symptom: at the last component's own knot the model was halved.
    p = ComponentArray(;
        log_A = 0.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 1.0, β = 0.0), mkcomp(; log_f_stop = 2.0, β = 0.5)],
        transition_width = 0.1
    )
    s = mapple(f, p)
    i = argmin(abs.(f .- 100.0))                    # f = 10^2, the final knot
    @test s[i] ≈ mapple([100.0], p)[1]
    @test s[i] > 0.9 * (s[i - 1] + s[i + 1]) / 2    # no half-weight notch at the knot
    # ... and the top decade follows the second segment's slope, not a taper toward zero.
    hi = findall(>(exp10(2.2)), f)
    slope = (log10(s[hi[end]]) - log10(s[hi[1]])) / (log10(f[hi[end]]) - log10(f[hi[1]]))
    @test slope ≈ 0.5 atol = 1.0e-6

    # The first component is likewise unwindowed below: the bottom decade follows ITS slope. The
    # tolerance is looser than the top decade's because the INTERIOR crossfade still leaks here ---
    # 0.4 decades below the knot component 1's weight is (1 + tanh 4)/2 = 0.99966, not 1, drifting
    # ~3e-4 across the decade. A closing shoulder on an outer edge would instead taper by O(0.1).
    lo = findall(<(exp10(0.6)), f)
    slope_lo = (log10(s[lo[end]]) - log10(s[lo[1]])) / (log10(f[lo[end]]) - log10(f[lo[1]]))
    @test slope_lo ≈ 0.0 atol = 1.0e-3
end

@testitem "mapple: the final component's log_f_stop is inert" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    f = exp10.(range(0, 3, length = 300))
    # `n` components carry `n - 1` breakpoints; the last `log_f_stop` bounds nothing and neither
    # windows the segment nor enters the amplitude chaining, so moving it must not touch the model.
    build(stop) = ComponentArray(;
        log_A = 0.4, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 1.2, β = -0.8), mkcomp(; log_f_stop = stop, β = -2.4)],
        transition_width = 0.08
    )
    reference = mapple(f, build(1.3))
    for stop in (1.25, 2.0, 3.0, 12.0, 500.0)
        @test mapple(f, build(stop)) == reference
    end

    # Sorting is by `log_f_stop`, so a final knot BELOW its neighbour reorders the components; the
    # segments then swap roles and the curve legitimately differs. Guard the inert claim against
    # being read as "any value", and confirm sorting still yields a finite, positive spectrum.
    swapped = mapple(f, build(0.5))
    @test all(isfinite, swapped) && all(>(0), swapped)

    # A single component has neither edge, so its knot is inert too.
    one_comp(stop) = ComponentArray(;
        log_A = 0.4, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = stop, β = -1.1)], transition_width = 0.08
    )
    @test mapple(f, one_comp(0.5)) == mapple(f, one_comp(9.0))
end

@testitem "mapple: the inert final knot is held at the top of the fitted band" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    # The final `log_f_stop` has exactly zero loss-gradient (it windows nothing), so left free it is
    # still boxed and Fminbox's barrier drifts it through its box --- 1.25 decades on this fit ---
    # towards its neighbour, where `sortperm` reads a crossing as a reordering of the components.
    # It is pinned to the top of the fitted band instead, so the accessor keeps a uniform shape and
    # a meaningful value.
    log_f = collect(range(0, 3, length = 400)); f = exp10.(log_f)
    truth = ComponentArray(;
        log_A = 0.5, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 9.0, β = -2.5)],
        transition_width = 0.1
    )
    Random.seed!(3); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))

    for n in 2:4
        init = fit_mapple(log_f, log_s; components = n, peaks = 0)
        refined = fit_mapple(log_f, log_s, init)
        m = MAPPLE(refined)
        @test length(breakpoints(m)) == n                       # accessor shape is unchanged
        @test breakpoints(m)[end] ≈ maximum(log_f) atol = 1.0e-5 # ... and the last entry is the band top
        @test breakpoints(m)[end] == maximum(breakpoints(m))    # the pinned knot stays sorted last
    end

    # An explicit `fix` is applied after the automatic pin, so a caller can still override it.
    init = fit_mapple(log_f, log_s; components = 2, peaks = 0)
    over = fit_mapple(log_f, log_s, init; fix = ["components[2].log_f_stop" => 2.0])
    @test MAPPLE(over).params.components[2].log_f_stop ≈ 2.0 atol = 1.0e-5
end

@testitem "mapple: a top segment ending inside the data recovers its slope" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    # Truth built ANALYTICALLY (not through `mapple`), so the test cannot be satisfied by
    # reproducing the model's own windowing. A Fano-factor shape: a flat shot-noise floor, a rise,
    # then saturation --- the case that exposed the bug, where the top segment's knot is at 10^2
    # with a full decade of data above it.
    log_f = range(0, 3, length = 400); f = exp10.(log_f)
    truth_log(x) = x < 1.0 ? 0.0 : (x < 2.0 ? 0.28 * (x - 1.0) : 0.28)
    Random.seed!(7); log_s = truth_log.(log_f) .+ 0.004 .* randn(length(f))

    init = fit_mapple(log_f, log_s; components = 3, peaks = 0)
    refined = fit_mapple(log_f, log_s, init; fix = ["components[1].β" => 0.0])
    m = MAPPLE(refined)
    MapplePlots.save_fit("top_segment_inside_data", f, exp10.(log_s), refined; init = init, subdir = "med")

    βs, knots = betas(m), breakpoints(m)
    @test βs[1] ≈ 0.0 atol = 1.0e-6            # pinned floor
    @test βs[2] ≈ 0.28 atol = 0.05             # the rise
    @test βs[3] ≈ 0.0 atol = 0.05              # the plateau: a slope, not a fade compensator
    @test knots[1] ≈ 1.0 atol = 0.15
    @test knots[2] ≈ 2.0 atol = 0.15

    # The fitted curve must track the data ACROSS the top decade, above the final knot. This is the
    # assertion the bug failed: the closing shoulder attenuated the model there by up to log10(2).
    top = findall(>(exp10(knots[2])), f)
    @test length(top) > 50
    @test maximum(abs.(log_s[top] .- log10.(mapple(f, refined))[top])) < 0.02
end

@testitem "mapple: BIC charges only the parameters actually searched" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random
    using TimeseriesTools: _held, _dof, _bic
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

    log_f = collect(range(0, 3, length = 400)); f = exp10.(log_f)
    truth = ComponentArray(;
        log_A = 0.5, peaks = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0),
        components = [
            ComponentArray(; log_f_stop = 1.0, β = -1.0), ComponentArray(; log_f_stop = 2.0, β = -2.0),
            ComponentArray(; log_f_stop = 9.0, β = -3.0),
        ], transition_width = 0.1
    )
    Random.seed!(3); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))
    p = fit_mapple(log_f, log_s; components = 3, peaks = 0)

    # 3 dof per component (β + a 2-dof free knot) + 2 global, minus what the fit never searched:
    # a held knot returns 2, any other held parameter 1.
    @test _dof(p, Dict{Int, Float64}()) == 2 + 3 * 3
    @test _dof(p, _held(log_f, p, nothing)) == 11 - 2            # the inert knot is never searched
    @test _dof(p, _held(log_f, p, ["components[1].β" => 0.0])) == 11 - 2 - 1
    @test _dof(p, _held(log_f, p, ["components[1].β" => 0.0, "transition_width" => 0.1])) == 11 - 2 - 2
    @test _dof(p, _held(log_f, p, ["components[1].log_f_stop" => 1.0])) == 11 - 2 - 2  # a knot is 2

    # The automatic pin refunds the same 2 dof at every candidate count, so it cannot shift the
    # component selection --- it only makes the absolute BIC right.
    for nc in 1:4
        q = fit_mapple(log_f, log_s; components = nc, peaks = 0)
        @test (2 + 3nc) - _dof(q, _held(log_f, q, nothing)) == 2
    end

    # BIC must be scored under the same weights the candidates were refined with (`_select_mapple`
    # resolves `w` out of `refine`); an unweighted score on a weighted fit is a different criterion.
    w = TimeseriesTools.logweights(log_f)
    @test _bic(p, f, log_s, w) != _bic(p, f, log_s)
    @test _bic(p, f, log_s, nothing, _held(log_f, p, nothing)) < _bic(p, f, log_s)  # fewer dof, lower BIC

    # The property that must not regress: a genuine three-segment spectrum still selects three.
    m = fit(MAPPLE, f, exp10.(log_s))
    @test length(betas(m)) == 3
    @test sort(betas(m)) ≈ [-3.0, -2.0, -1.0] atol = 0.15
    MapplePlots.save_fit("bic_dof_accounting", f, exp10.(log_s), m; subdir = "med")
end

@testitem "mapple: vcov/stderror are calibrated and skip held parameters" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics, LinearAlgebra
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    log_f = collect(range(0, 3, length = 300)); f = exp10.(log_f)
    truth = ComponentArray(;
        log_A = 0.5, peaks = nopeaks(),
        components = [
            mkcomp(; log_f_stop = 1.0, β = -1.0), mkcomp(; log_f_stop = 2.0, β = -2.0),
            mkcomp(; log_f_stop = 9.0, β = -3.0),
        ], transition_width = 0.1
    )
    function fitone(seed; fix = nothing)
        Random.seed!(seed)
        ls = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))
        m = MAPPLE(fit_mapple(log_f, ls, fit_mapple(log_f, ls; components = 3, peaks = 0); fix))
        return m, Timeseries(exp10.(ls), f)
    end

    m, spec = fitone(1)
    se = stderror(m, spec)
    Σ = vcov(m, spec)

    @test se isa ComponentArray && length(se) == length(m.params)  # shaped like the parameters
    @test all(≥(0), se)
    @test se.components[3].log_f_stop == 0                          # the inert knot is not searched
    @test size(Σ) == (length(freelabels(m, spec)), length(freelabels(m, spec)))
    @test Σ ≈ Σ' rtol = 1.0e-8
    @test all(>(0), diag(Σ))
    @test length(freelabels(m, spec)) == length(m.params) - 1       # everything but the inert knot
    @test "components[3].log_f_stop" ∉ freelabels(m, spec)

    # Calibration: on independent noise the Laplace SE must predict the estimator's real scatter.
    βs = [betas(first(fitone(100 + s)))[2] for s in 1:15]
    @test se.components[2].β ≈ std(βs) rtol = 0.4

    # A held parameter is a constant: zero variance, and dropped from the covariance.
    fix = ["components[1].β" => -1.0]
    mh, spech = fitone(1; fix)
    seh = stderror(mh, spech; fix)
    @test seh.components[1].β == 0
    @test "components[1].β" ∉ freelabels(mh, spech; fix)
    @test length(freelabels(mh, spech; fix)) == length(m.params) - 2
end

@testitem "mapple: :auto skips fix labels a candidate does not have" setup = [MapplePlots] tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    log_f = collect(range(0, 3, length = 300)); f = exp10.(log_f)
    truth = ComponentArray(;
        log_A = 0.5, peaks = nopeaks(),
        components = [
            mkcomp(; log_f_stop = 1.0, β = -1.0), mkcomp(; log_f_stop = 2.0, β = -2.0),
            mkcomp(; log_f_stop = 9.0, β = -3.0),
        ], transition_width = 0.1
    )
    Random.seed!(3); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))
    spec = Timeseries(exp10.(log_s), f)

    # `"components[3].β"` does not exist for the 1- and 2-component candidates. Before, the sweep
    # threw there; now those candidates are fitted without it and only the 3-component one holds it.
    m = @test_nowarn fit(
        MAPPLE, spec; peaks = 0, components = :auto, max_components = 3,
        refine = (; fix = ["components[1].β" => -1.0, "components[3].β" => -3.0])
    )
    @test length(betas(m)) == 3
    @test betas(m)[1] == -1.0        # count-independent: held exactly
    @test betas(m)[3] == -3.0        # count-dependent: held exactly where it applies

    # A count-independent pin alone still works across the whole sweep.
    m2 = fit(
        MAPPLE, spec; peaks = 0, components = :auto, max_components = 3,
        refine = (; fix = ["transition_width" => 0.1])
    )
    @test m2.params.transition_width == 0.1

    # A label no candidate has is a typo, and must still be rejected rather than silently ignored.
    @test_throws ArgumentError fit(
        MAPPLE, spec; peaks = 0, components = :auto, max_components = 3,
        refine = (; fix = ["components[1].beta" => 0.0])
    )
    # ... including one that is merely out of range for EVERY candidate.
    @test_throws ArgumentError fit(
        MAPPLE, spec; peaks = 0, components = :auto, max_components = 3,
        refine = (; fix = ["components[9].β" => 0.0])
    )
end
