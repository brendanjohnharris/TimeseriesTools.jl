# Comprehensive unit tests for the MAPPLE spectral model.
#
#   Core model:        src/Mapple.jl    (`mapple`, `mapple!`, `MAPPLE`, `mapple_sort!`,
#                                         `fit`, `predict`)
#   Optim refinement:  ext/OptimExt.jl  (`fit_mapple` 3-arg, `fit!`)
#
# Characterization tests pin down the current, correct behaviour of the model. The later
# items also cover the robustness/reporting overhaul: the conservative peak-amplitude init
# and the `log_A` do-block-leak fix, the Hz reporting accessors, the `rsquared` flat-spectrum
# guard, the peak-width guard, data-driven (`:auto`) peak counting, and BIC component
# selection. These tests are deliberately plot-free and deterministic.
#
# Construction helpers are repeated inside each `@testitem` because every item runs in
# its own isolated module.

@testitem "mapple model: power law, peaks, positivity" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    f = exp10.(range(0, 3, length = 200))

    # A single component with no peaks is exactly the power law A·f^β.
    p1 = ComponentArray(; log_A = 0.5, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1)
    s1 = mapple(f, p1)
    @test s1 ≈ exp10(0.5) .* f .^ (-2.0)
    @test all(>(0), s1)                # strictly positive => log10 in the loss is safe
    @test all(isfinite, s1)

    # The two `mapple` call forms agree.
    @test mapple(f, p1) ≈ mapple(f, p1[[:log_A, :transition_width, :components]], p1[[:peaks]])

    # A Gaussian peak adds power near its centre frequency (10^log_f).
    p2 = ComponentArray(; log_A = 0.5,
        peaks = [mkpeak(; log_f = 1.0, log_σ = 0.3, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1)
    s2 = mapple(f, p2)
    icentre = argmin(abs.(f .- 10.0))
    @test s2[icentre] > s1[icentre]
    @test all(>(0), s2)

    # Output is invariant to the order components are supplied (sorted internally).
    pfwd = ComponentArray(; log_A = 1.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 5.0, β = -3.0)],
        transition_width = 0.05)
    prev = ComponentArray(; log_A = 1.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = -3.0), mkcomp(; log_f_stop = 1.5, β = -1.0)],
        transition_width = 0.05)
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
    p = ComponentArray(; log_A = 1.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 5.0, β = -3.0)],
        transition_width = 0.05)
    s = mapple(f, p)
    @test all(isfinite, s)
    @test all(>(0), s)
    @test @inferred(mapple(f, p)) isa Vector{Float64}

    # A peak added on top only increases power.
    pp = ComponentArray(; log_A = 1.0,
        peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 2.0)],
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 5.0, β = -3.0)],
        transition_width = 0.05)
    @test all(mapple(f, pp) .≥ s .- 1e-8)
end

@testitem "MAPPLE: construction & predict" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

    # The keyword constructor mirrors the documented field layout.
    m = MAPPLE(; log_A = 1.0,
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5)],
        components = [mkcomp(; log_f_stop = 1.5, β = 2.0)],
        transition_width = 0.2)
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
    mkparams() = ComponentArray(; log_A = 1.0,
        peaks = [mkpeak(; log_f = 2.5, log_σ = 0.11, log_A = 0.51),
            mkpeak(; log_f = 0.7, log_σ = 0.12, log_A = 1.52)],
        components = [mkcomp(; log_f_stop = 5.0, β = 4.0),
            mkcomp(; log_f_stop = 1.5, β = 2.0)],
        transition_width = 0.2)

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
    msorted = MAPPLE(ComponentArray(; log_A = 1.0,
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.12, log_A = 1.52),
            mkpeak(; log_f = 2.5, log_σ = 0.11, log_A = 0.51)],
        components = [mkcomp(; log_f_stop = 1.5, β = 2.0),
            mkcomp(; log_f_stop = 5.0, β = 4.0)],
        transition_width = 0.2))
    ssorted = mapple(f, msorted.params)
    sort!(msorted)
    @test mapple(f, msorted.params) ≈ ssorted
end

@testitem "mapple: fit recovers a known spectrum (peaks=0)" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

    truth = ComponentArray(;
        components = [mkcomp(; β = 2.0, log_f_stop = 1.5), mkcomp(; β = 4.0, log_f_stop = 5.0)],
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5),
            mkpeak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5)],
        transition_width = 0.15, log_A = 1.0)

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
end

@testitem "mapple: type-form entry points & empty-block construction" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim

    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)

    f = exp10.(range(0, 3, length = 200))
    p = ComponentArray(; log_A = 0.5, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1)
    s = mapple(f, p)

    # `fit(MAPPLE, freqs, spectrum)` wraps the frequency vector in a frequency dimension
    # and fits.
    @test fit(MAPPLE, f, s; peaks = 0, components = 1) isa MAPPLE

    # `fit!(::Type{MAPPLE}, ...)` was a never-working misnomer for `fit`; it is removed,
    # so calling it raises a MethodError.
    @test_throws MethodError fit!(MAPPLE, f, s)

    # A zero-peak model builds via the `MAPPLE` constructor even from a plain empty
    # `Vector{<:ComponentArray}` (the constructor normalises empty blocks).
    PeakT = typeof(mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0))
    m0 = MAPPLE(; log_A = 0.5, peaks = Vector{PeakT}(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1)
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
    truth = ComponentArray(; log_A = 1.0,
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5),
            mkpeak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5)],
        components = [mkcomp(; log_f_stop = 1.5, β = 2.0), mkcomp(; log_f_stop = 5.0, β = 4.0)],
        transition_width = 0.15)
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
    m = MAPPLE(; log_A = 1.0,
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5)],
        components = [mkcomp(; log_f_stop = 1.5, β = 2.0), mkcomp(; log_f_stop = 5.0, β = 4.0)],
        transition_width = 0.2)
    @test betas(m) isa Vector{Float64}
    @test betas(m) == [2.0, 4.0]
    @test breakpoints(m) == [1.5, 5.0]
    @test peakfreqs(m) == [0.7]
    @test peaksigmas(m) == [0.1]
    @test peakamplitudes(m) == [1.5]
    # empty peaks must not error and stay typed
    PeakT = typeof(mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0))
    m0 = MAPPLE(; log_A = 1.0, peaks = Vector{PeakT}(),
        components = [mkcomp(; log_f_stop = 5.0, β = 2.0)], transition_width = 0.2)
    @test peakfreqs(m0) == Float64[]
    @test peakfreqs(m0) isa Vector{Float64}
end

@testitem "MAPPLE: diagnostics" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    truth = MAPPLE(; log_A = 1.0,
        peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1)
    f = exp10.(range(0, 3, length = 300))
    spectrum = ToolsArray(mapple(f, truth.params), 𝑓(f))
    # the exact model on its own spectrum: zero residual, perfect fit
    @test maximum(abs, mapple_residuals(truth, spectrum)) < 1.0e-8
    @test mapple_loss(truth, spectrum) < 1.0e-12
    @test rsquared(truth, spectrum) ≈ 1.0 atol = 1.0e-8
    # a worse model has higher loss and lower R²
    wrong = MAPPLE(; log_A = 1.0,
        peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 5.0, β = -1.0)], transition_width = 0.1)
    @test mapple_loss(wrong, spectrum) > mapple_loss(truth, spectrum)
    @test rsquared(wrong, spectrum) < rsquared(truth, spectrum)
end

@testitem "MAPPLE: precomputed log_f fast path" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    p = ComponentArray(; log_A = 0.5, peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 2.0, β = -1.0), mkcomp(; log_f_stop = 5.0, β = -3.0)],
        transition_width = 0.1)
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

@testitem "mapple: robust fit on a steep rising spectrum (regression)" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Optim, ForwardDiff, Random, Statistics
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    # Steep rising background (β = 2, 4) with two peaks: the case that used to collapse to a
    # degenerate β / log_A optimum. The zero-peak detrend + single joint refine fits it cleanly.
    truth = ComponentArray(;
        components = [mkcomp(; β = 2.0, log_f_stop = 1.5), mkcomp(; β = 4.0, log_f_stop = 5.0)],
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5),
            mkpeak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5)],
        transition_width = 0.15, log_A = 1.0)
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
end

@testitem "mapple: rough-fit base amplitude is not clobbered by the peak loop" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Random
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    # The peak-construction `do` block must use fresh names: assigning `log_A` inside it would
    # leak to the enclosing `log_A = first(log_s)` and corrupt the base amplitude to the last
    # peak's value (a catastrophic init that biases the subsequent refine).
    truth = ComponentArray(;
        components = [mkcomp(; β = 2.0, log_f_stop = 1.5), mkcomp(; β = 4.0, log_f_stop = 5.0)],
        peaks = [mkpeak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5), mkpeak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5)],
        transition_width = 0.15, log_A = 1.0)
    log_f = range(0, 3, length = 500); f = exp10.(log_f)
    Random.seed!(0); log_s = log10.(mapple(f, truth)) .+ 0.05 .* randn(length(f))
    init = fit_mapple(log_f, log_s; components = 2, peaks = 2, w = 50)
    @test init.log_A ≈ first(log_s)              # base amplitude is the DC estimate, not a peak's
    @test all(isfinite, mapple(f, init))         # init is sane, not a 10^(peak log_A) blow-up
end

@testitem "MAPPLE: Hz reporting accessors" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    m = MAPPLE(; log_A = 1.0, peaks = [mkpeak(; log_f = 1.0, log_σ = 0.2, log_A = 0.5)],
        components = [mkcomp(; log_f_stop = 1.5, β = -1.0), mkcomp(; log_f_stop = 2.5, β = -3.0)],
        transition_width = 0.1)
    # Physical-unit accessors convert OUT of log-10 space at the reporting boundary.
    @test peakcentres(m) ≈ exp10.(peakfreqs(m)) ≈ [10.0]
    @test breakfrequencies(m) ≈ exp10.(breakpoints(m)) ≈ [10.0^1.5, 10.0^2.5]
    @test peakheights(m) ≈ exp10.(peakamplitudes(m)) ≈ [10.0^0.5]
    # FWHM in Hz: σ = f·tanh(log_σ), FWHM = 2√(2ln2)·σ.
    @test peakbandwidths(m) ≈ [2 * sqrt(2 * log(2)) * 10.0 * tanh(0.2)]
    @test all(>(0), peakbandwidths(m))
    # Empty peaks stay typed and empty.
    PeakT = typeof(mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0))
    m0 = MAPPLE(; log_A = 1.0, peaks = Vector{PeakT}(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1)
    @test peakcentres(m0) == Float64[] && peakbandwidths(m0) == Float64[] && peakheights(m0) == Float64[]
end

@testitem "MAPPLE: rsquared guards a flat spectrum" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    m = MAPPLE(; log_A = 1.0, peaks = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0),
        components = [mkcomp(; log_f_stop = 5.0, β = 0.0)], transition_width = 0.1)
    f = exp10.(range(0, 3, length = 100))
    flat = ToolsArray(fill(10.0, 100), 𝑓(f))   # constant spectrum => ss_tot = 0
    @test isnan(rsquared(m, flat))             # undefined, not ±Inf or a DivideError
end

@testitem "mapple: peak-width guard rejects out-of-range detections" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Random
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)
    # Falling 1-component background with one clear, narrow peak.
    truth = ComponentArray(; log_A = 1.0, peaks = [mkpeak(; log_f = 1.5, log_σ = 0.1, log_A = 1.0)],
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1)
    log_f = range(0, 3, length = 400); f = exp10.(log_f)
    Random.seed!(7); log_s = log10.(mapple(f, truth)) .+ 0.02 .* randn(length(f))
    # Default limits: the narrow peak is found.
    @test length(fit_mapple(log_f, log_s; components = 1, w = 50).peaks) ≥ 1
    # Demanding very wide peaks rejects every detection before counting.
    @test length(fit_mapple(log_f, log_s; components = 1, w = 50, peak_width_limits = (2.0, 3.0)).peaks) == 0
end

@testitem "mapple: :auto detects a prominent peak and rejects noise" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Random
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    truth = ComponentArray(; log_A = 1.0, peaks = [mkpeak(; log_f = 1.5, log_σ = 0.12, log_A = 1.2)],
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1)
    log_f = range(0, 3, length = 400); f = exp10.(log_f)
    Random.seed!(11); log_s = log10.(mapple(f, truth)) .+ 0.03 .* randn(length(f))
    p = fit_mapple(log_f, log_s; components = 1, w = 50)   # peaks = :auto
    @test length(p.peaks) == 1                               # the one real peak, not noise bumps
    @test only(collect(Float64, p.peaks.log_f)) ≈ 1.5 atol = 0.2
end

@testitem "mapple: BIC selects the simplest adequate component count" tags = [:mapple] begin
    using TimeseriesTools, ComponentArrays, Random, Optim, ForwardDiff
    mkcomp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
    mkpeak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))
    nopeaks() = map(_ -> mkpeak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)
    # A genuine single power law: the free-knot BIC penalty must not over-select components.
    truth = ComponentArray(; log_A = 1.0, peaks = nopeaks(),
        components = [mkcomp(; log_f_stop = 5.0, β = -2.0)], transition_width = 0.1)
    f = exp10.(range(0, 3, length = 150))
    Random.seed!(5); s = exp10.(log10.(mapple(f, truth)) .+ 0.03 .* randn(length(f)))
    spec = ToolsArray(s, 𝑓(f))
    m = fit(MAPPLE, spec)                       # :auto BIC selection (refined when Optim loaded)
    @test length(m.params.components) == 1
    @test rsquared(m, spec) > 0.99
end
