@testitem "AutocorrelationsExt" begin # Optimize this some more?
    using Autocorrelations, StatsBase, BenchmarkTools, MeanSquaredDisplacement, CairoMakie,
          Unitful, TimeseriesMakie
    import TimeseriesTools: Timeseries
    x = colorednoise(1:10)
    @test Autocorrelations.default_lags(x) == 0:1:9

    lags = 1:5
    out1 = Vector{eltype(x)}(undef, length(lags))
    out2 = Vector{eltype(x)}(undef, length(lags))
    acf!(out1, x, lags)
    acf!(out2, parent(x), lags)
    @test out1 == out2
    @test parent(acf(x)) == acf(parent(x))
    @test times(acf(x)) == range(0, length(x) - 1) * samplingperiod(x)
    @test fftacf(x) == fftacf(parent(x))
    @test acf(x) == dotacf(x)
    @test acf(x) != fftacf(x)
    @test acf(x) ≈ fftacf(x)

    x = colorednoise(0.01:0.01:100)
    lags = 1:1000
    out1 = Vector{eltype(x)}(undef, length(lags))
    out2 = Vector{eltype(x)}(undef, length(lags))
    acf!(out1, x, lags)
    acf!(out2, parent(x), lags)
    @test out1 == out2
    @test parent(acf(x)) == acf(parent(x))
    @test times(acf(x)) == range(0, length(x) - 1) * samplingperiod(x)
    @test fftacf(x) == fftacf(parent(x))
    @test acf(x) == fftacf(x)
    @test acf(x) != dotacf(x)
    @test acf(x) ≈ dotacf(x)
    @test !(acf(x .+ 10; demean = true) ≈ acf(x .+ 10; demean = false))

    # ? StatsBase is biased
    y = set(x, sin.(times(x)))
    r1 = autocor(y, 800:1000; demean = true) # StatsBase method is biased
    r2 = dotacf(y, 800:1000; demean = false, normalize = true)
    @test minimum(r2)≈-1 rtol=1e-3
    @test minimum(r1) > -0.95 # This is no bueno

    m = @test_nowarn msdist(x)
    @test m == imsd(parent(x))
    @test m == msdist(parent(x))
    a = @benchmark msdist($x)
    b = @benchmark msdist($(parent(x)))
    c = @benchmark imsd($(parent(x))) # Compare to MeanSquaredDisplacement
    @test median(a.times) < median(c.times) .* 2
    @test b.allocs ≤ c.allocs

    x = Timeseries(cumsum(randn(10000, 100), dims = 1), 0.1:0.1:1000, 1:100)
    m = msdist(x, 1:1000)
    traces(m; linecolor = (:gray, 0.1), axis = (; xscale = log10, yscale = log10))
    m = dropdims(mean(m, dims = 2), dims = 2)
    plot!(decompose(m)...)
    current_figure()
end

@testitem "DSPExt" begin
    using DSP
    using CairoMakie
    using TimeseriesMakie
    using TimeseriesTools
    import TimeseriesTools.Timeseries # or TS
    using StatsBase

    N = 100000
    dt = 0.005
    t = dt:dt:10
    x = [0.00 .* colorednoise(t) .+ sin.(2 * t .+ 2 * randn()) for _ in 1:200]
    y = ToolsArray(x, Var(1:200)) |> stack
    x̂ = Timeseries(vcat(collect.(x)...), dt:dt:(sum(length.(x)) * dt))
    x = phasestitch(x)

    p = @test_nowarn heatmap(y) # Should default to the DimensionalData recipe
    @test p.plot isa CairoMakie.Heatmap

    pargs = Makie.convert_arguments(Traces, y)
    @test pargs[1] == lookup(y, 1) |> collect
    @test pargs[2] == lookup(y, 2) |> collect
    @test pargs[3] isa Matrix

    f = Figure()
    ax = Axis(f[1, 1])
    p = @test_nowarn traces!(ax, y[Var(1:5)])
    f

    @test_nowarn plot(x[𝑡(1:10000)])
    plot(x̂[𝑡(1500:(length(t) * 5))])

    # And a power spectrum of a 'perfect' signal
    _t = dt:dt:(dt * N)
    p = Timeseries(sin.(2 * _t), _t)
    S′ = powerspectrum(p, dt * 4)
    @test_nowarn spectrumplot(S′)

    # Power spectrum of the concatenated time series
    Ŝ = powerspectrum(x̂[1:N], dt * 4)
    @test_nowarn spectrumplot(Ŝ)

    # Power spectrum of the phasestitched time series
    S = powerspectrum(x[1:N], dt * 4)
    fax = @test_nowarn spectrumplot(S)

    pac = autocor(p, [10])[1]
    @test ≈(pac, autocor(x[𝑡(1:10000)] |> collect, [10])[1]; rtol = 1e-2)
    # @test pac - autocor(x̂[ 𝑡(1:10000)] |> collect, [10])[1] >   pac -
    # autocor(x[ 𝑡(1:10000)] |> collect, [10])[1]
end

@testitem "ContinuousWaveletsExt" begin
    import TimeseriesTools: Timeseries
    using ContinuousWavelets, BenchmarkTools
    # Define a test time series
    fs = 200
    t = range(0, stop = 5, length = 100 * fs + 1)
    x = (0.8 .* sin.(2 * π * 40 * t) + 1.1 .* sin.(2 * π * 100 * t)) .^ 2
    ts = x = TimeseriesTools.Timeseries(x, t)
    f_min = fs / 100
    S = waveletspectrogram(x)
    @test S isa TimeseriesTools.TimeseriesBase.Spectra.RegularSpectrogram

    # Multivariate
    x = cat(ts, ts .* randn(length(ts)), dims = Var(1:2))
    S = @test_nowarn waveletspectrogram(x)
    @test all(isa.(dims(S), (𝑡, 𝑓, Var)))

    # GPU test
    if false
        using CUDA
        using BenchmarkTools
        BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60
        @benchmark waveletspectrogram(x)
        x = CuArray(x)
        @benchmark CUDA.@sync waveletspectrogram(x)
    end

    if false
        using CUDA
        x = cat(ts, ts .* randn(length(ts)), dims = Var(1:2))
        S = @test_nowarn waveletspectrogram(x)
        @test all(isa.(dims(S), (𝑡, 𝑓, Var)))

        y = set(x, CuArray(x.data))
        S = @test_nowarn waveletspectrogram(y)
        @test all(isa.(dims(S), (𝑡, 𝑓, Var)))

        @test all(x .== y)
        @test dims(x) == dims(y)
    end
end

@testitem "TimeseriesSurrogatesExt" begin
    using StatsBase, TimeseriesSurrogates
    θ = 3 # A Fano-factor of 3
    μ = 1
    α = μ / θ # A mean of 1
    N = 500000
    x = gammarenewal(N, α, θ)
    dt = diff(times(x))
    F = var(dt) / mean(dt)
    @test x isa SpikeTrain
    @test F≈θ rtol=5e-2
    @test mean(dt)≈α * F rtol=5e-2

    # Jitter surrogate
    y = set(x, 𝑡 => surrogate(times(x), RandomJitter(0.1, 0.1)))
    @test y isa SpikeTrain
    @test issorted(times(y))
    @test minimum(times(y))≈minimum(times(x)) atol=0.5
    @test maximum(times(y))≈maximum(times(x)) atol=0.5
    @test x != y
    sur = @test_nowarn surrogenerator(times(x), RandomJitter(0.1, 0.1))
    @test all(copy(sur()) .!= sur())

    # Gamma renewal surrogate
    y = set(x, 𝑡 => surrogate(times(x), GammaRenewal()))
    dt̂ = diff(times(y))
    F̂ = var(dt̂) / mean(dt̂)
    @test y isa SpikeTrain
    @test issorted(times(y))
    @test F̂≈θ rtol=5e-2
    @test mean(dt̂)≈α * F̂ rtol=5e-2
    @test minimum(times(y))≈minimum(times(x)) atol=10 * μ
    @test maximum(times(y))≈maximum(times(x)) atol=0.01 * N
end

# @testitem "DiffEqBaseExt" begin using DifferentialEquations f(u, p, t) = 1.01 * u u0 = 1 /
#     2 tspan = (0.0, 1.0) prob = ODEProblem(f, u0, tspan, saveat=0.1) sol = solve(prob)

#     x = Timeseries(sol)
# end

@testitem "GeneralizedPhaseExt" begin
    using GeneralizedPhase, Unitful
    x = bandpass(colorednoise(0.01:0.01:10), (10, 15))
    X = [bandpass(colorednoise(0.1:0.1:100), (0.1, 0.5)) for _ in 1:10]
    X = ToolsArray(X, Var(1:10)) |> stack
    _ϕ = @test_nowarn _generalized_phase(x)
    ϕ = @test_nowarn _generalized_phase(X)

    x = set(x, 𝑡 => lookup(x, 𝑡).data * u"s")
    X = set(X, 𝑡 => lookup(X, 𝑡).data * u"s")

    ϕ = @test_nowarn _generalized_phase(x)
    ϕ = @test_nowarn _generalized_phase(X)
end

@testitem "ComplexityMeasuresExt" begin
    using Distributions, LinearAlgebra, ComplexityMeasures
    μ = [1.0, -4.0]
    σ = [2.0, 2.0]
    𝒩 = MvNormal(μ, LinearAlgebra.Diagonal(map(abs2, σ)))
    N = 500
    D = Timeseries(hcat(sort([rand(𝒩) for i in 1:N])...)', 1:N, 1:2)
    p = probabilities(NaiveKernel(1.5), StateSpaceSet(D))

    ComplexityMeasures.entropy(Shannon(), ValueBinning(RectangularBinning(100)),
                               StateSpaceSet(D))
end

@testitem "Upsampling" begin
    using DataInterpolations
    using Unitful

    # * Vector
    ts = (rand(100) .- 0.5) * 2π
    x = Timeseries(sinc.(ts), ts)
    @test_throws "`interpolate` only supports forward" TimeseriesTools.interpolate(x)

    sort!(ts)
    x = Timeseries(sinc.(ts), ts)
    itp = TimeseriesTools.interpolate(x)
    y = itp(dims(x) |> only)
    @test x ≈ y

    ts = (-π):0.1:π
    x = Timeseries(sinc.(ts), ts)
    itp = TimeseriesTools.interpolate(x)
    z = @test_nowarn upsample(x, 2)
    @test all(z[𝑡(At(ts))] .≈ x)
    @test z≈sinc.(lookup(z) |> only) atol=1e-2

    # * Matrix
    x = Timeseries(randn(100, 100), 0.1:0.1:10, Var(1:100))
    y = @test_nowarn upsample(x, 2)
    z = @test_nowarn upsample(x, 2, dims = 1)
    @test y == z
    z = @test_nowarn upsample(x, 2, dims = (1, 2))
    @test all(y .≈ z[𝑡(At(lookup(y, 𝑡))), Var(At(lookup(y, Var)))])
    @test all(x .≈ z[𝑡(At(lookup(x, 𝑡))), Var(At(lookup(x, Var)))])
    @test DimensionalData.name.(dims(x)) == DimensionalData.name.(dims(z))
    @test length(dims(z, 1)) == length(dims(z, 2)) == 199

    # * 3D array
    x = Timeseries(randn(100, 100, 100), 0.1:0.1:10, Var(1:100), X(1:100))
    y = @test_nowarn upsample(x, 2)
    z = @test_nowarn upsample(x, 2, dims = 1)
    @test y == z
    z = @test_nowarn upsample(x, 2, dims = (3, 2, 1))
    @test all(y .≈ z[𝑡(At(lookup(y, 𝑡))), Var(At(lookup(y, Var))), X(At(lookup(y, X)))])
    @test all(x .≈ z[𝑡(At(lookup(x, 𝑡))), Var(At(lookup(x, Var))), X(At(lookup(x, X)))])
    @test DimensionalData.name.(dims(x)) == DimensionalData.name.(dims(z))
    @test length(dims(z, 1)) == length(dims(z, 2)) == 199

    # * Unitful data
    ts = ((-π):0.1:π) * u"s"
    x = Timeseries(sinc.(ustrip(ts)), ts) * u"V"
    itp = TimeseriesTools.interpolate(x)
    @test unit(eltype(itp(dims(x) |> only))) == NoUnits
    z = @test_nowarn upsample(x, 2)
    @test unit(eltype(z)) == u"V"
    @test unit(eltype(lookup(z) |> only)) == u"s"
    @test all(z[𝑡(At(ts))] .≈ x)
end

@testitem "Unit Power" begin
    using Unitful
    using Statistics
    N = UnitPower
    _X = Timeseries(rand(100), 0.01:0.01:1)
    y = copy(_X)
    T = fit(N, y)
    x = normalize(y, T)
    @test mean(x .^ 2) ≈ 1
    @test !isnothing(T.p)
    @test denormalize(x, T) ≈ y
    @test_nowarn normalize!(y, T)
    @test y == x
    @test_nowarn denormalize!(x, T)
    @test all(x .≈ _X)

    _X = set(y, 𝑡 => times(y) * u"s") * u"V"
    # _X = X * u"V"
    y = copy(_X)
    T = fit(N, y)
    x = normalize(y, T)
    @test ustripall(mean(x .^ 2)) ≈ 1
    @test !isnothing(T.p)
    @test denormalize(x, T) == y
end

@testitem "SpectralFittingExt" begin
    using SpectralFitting
    using SpectralFitting, CairoMakie, TimeseriesMakie
    using Statistics

    x = collect(range(-5, 5, 200))

    A_true = 3.0
    pos_true = 1.3
    sigma_true = 0.8
    err_true = 0.2

    y = @. A_true * exp(-(x - pos_true)^2 / (2 * sigma_true^2))

    y_noisy = y .+ (0.2 * randn(length(y)))

    z = ToolsArray(y_noisy, 𝑓(x))

    data = @inferred InjectiveData(x, y_noisy) # Time series gets passed in here
    zdata = @inferred InjectiveData(z)

    model = GaussianLine(μ = FitParam(0.0))
    prob = FittingProblem(model => data)
    result = SpectralFitting.fit(prob, LevenbergMarquadt())

    zprob = FittingProblem(model => zdata)
    zresult = SpectralFitting.fit(zprob, LevenbergMarquadt())
    @test result.u == zresult.u

    # * Multivariate
    Y = [y_noisy .+ 0.1 * randn(length(y)) for _ in 1:5] |> stack
    Z = ToolsArray(Y, (𝑓(x), Var(1:size(Y, 2))))
    zdata = @inferred InjectiveData(Z)
    zprob = @inferred FittingProblem(model => zdata)
    zresult = SpectralFitting.fit(zprob, LevenbergMarquadt())

    # * Try 1/f fitting
    α = 2.0
    x = colorednoise(0.001:0.001:10; α)
    s = spectrum(x)[𝑓 = OpenInterval(1, Inf)]
    model = PowerLaw()
    data = InjectiveData(s)
    prob = FittingProblem(model => data)
    result = SpectralFitting.fit(prob, LevenbergMarquadt()) # * Need to have a balanced optimizer that accounts for loglog spacing/errors
    update_model!(model, result)

    f = make_model_domain(ContiguouslyBinned(), data)
    y = invokemodel(f, model)

    plotspectrum(s)
    lines!(data.domain, y)
    display(current_figure())
end

@testitem "OptimExt" begin
    # For fitting power spectra
    include("./optim.jl")
end

@testitem "MAPPLE" begin
    using ComponentArrays
    components = ComponentVector([ComponentVector(; β = 4.0, log_f_stop = 5.0),
                                     ComponentVector(; β = 2.0, log_f_stop = 1.5)])
    peaks = ComponentVector([ComponentArray(; log_f = 2.5, log_σ = 0.1, log_A = 0.5),
                                ComponentVector(; log_f = 0.7, log_σ = 0.1, log_A = 1.5)])
    params = ComponentVector(; components, peaks, transition_width = 0.2, log_A = 1.0)

    log_f = range(0, 3, length = 500)
    f = map(exp10, log_f)
    s = mapple(f, params)
    log_s = log10.(s) .+ 0.1 * randn(length(s))
    s = exp10.(log_s)

    m = MAPPLE(params)
    @test !issorted(m.params.peaks.log_f)
    @test !issorted(m.params.components.log_f_stop)

    m2 = sort(m)
    @test issorted(m2.params.peaks.log_f)
    @test issorted(m2.params.components.log_f_stop)

    sort!(m)
    @test issorted(m.params.peaks.log_f)
    @test issorted(m.params.components.log_f_stop)
end

@testitem "MAPPLE MAD" begin
    using CairoMakie
    using Optim

    x = zeros(100000)
    x[1] = 0.0
    noise_strength = 0.1
    for i in 2:length(x)
        x[i] = x[i - 1] - 0.01 * x[i - 1] + 0.0001 * randn()
    end
    x = ToolsArray(x, 𝑡(1:length(x)))

    mad = madev(x)

    m = fit(MAPPLE, logsample(mad[10:end]); peaks = 0, components = 2)
    mad_guess = predict(m, mad)
    _mad = logsample(mad[1:end])
    fit!(m, _mad; show_trace = true)
    mad_refined = predict(m, mad)

    f = Figure()
    ax = CairoMakie.Axis(f[1, 1]; xscale = log10, yscale = log10)
    lines!(ax, mad; label = "MAD", color = :green)
    scatter!(ax, _mad; color = :green, markersize = 10)
    lines!(ax, mad_guess; label = "MAPPLE fit", color = :red)
    lines!(ax, mad_refined; label = "MAPPLE refined", color = :black)
    display(f)

    @test m.params.components[1].β≈-0.5 atol=0.05
    @test m.params.components[2].β≈0.0 atol=0.05
end
