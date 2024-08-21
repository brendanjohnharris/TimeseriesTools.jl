@testset "AutocorrelationsExt" begin # Optimize this some more?
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
    @test median(a.times) < median(c.times) .* 1.1
    @test b.allocs ≤ c.allocs

    x = TimeSeries(0.1:0.1:1000, 1:100, cumsum(randn(10000, 100), dims = 1))
    m = msdist(x, 1:1000)
    traces(m; color = (:gray, 0.4))
    m = dropdims(mean(m, dims = 2), dims = 2)
    plot!(decompose(m)...)
    current_figure()
end

@testset "DSPExt" begin
    using DSP
    using TimeseriesTools
    import TimeseriesTools.TimeSeries # or TS

    N = 100000
    dt = 0.005
    t = dt:dt:10
    x = [0.00 .* colorednoise(t) .+ sin.(2 * t .+ 2 * randn()) for _ in 1:200]
    y = hcat(Var(1:200), x...)
    x̂ = TimeSeries(dt:dt:(sum(length.(x)) * dt), vcat(collect.(x)...))
    x = phasestitch(x)

    @test_nowarn stackedtraces(y[Var(1:10)], spacing = :even, linewidth = 5, offset = 1.3;
                               axis = (; xlabel = "Time"))
    @test_nowarn plot(x[𝑡(1:10000)])
    plot(x̂[𝑡(1500:(length(t) * 5))])

    # And a power spectrumof a 'perfect' signal
    _t = dt:dt:(dt * N)
    p = TimeSeries(_t, sin.(2 * _t))
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

@testset "ContinuousWaveletsExt" begin
    # Define a test time series
    fs = 200
    t = range(0, stop = 5, length = 100 * fs + 1)
    x = (0.8 .* sin.(2 * π * 40 * t) + 1.1 .* sin.(2 * π * 100 * t)) .^ 2
    ts = x = TimeseriesTools.TimeSeries(t, x)
    f_min = fs / 100
    S = waveletspectrogram(x)
    @test S isa RegularSpectrogram

    # Multivariate
    x = cat(Var(1:2), ts, ts .* randn(length(ts)))
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
        x = cat(Var(1:2), ts, ts .* randn(length(ts)))
        S = @test_nowarn waveletspectrogram(x)
        @test all(isa.(dims(S), (𝑡, 𝑓, Var)))

        y = set(x, CuArray(x.data))
        S = @test_nowarn waveletspectrogram(y)
        @test all(isa.(dims(S), (𝑡, 𝑓, Var)))

        @test all(x .== y)
        @test dims(x) == dims(y)
    end
end

@testset "TimeseriesSurrogatesExt" begin
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
    @test minimum(times(y))≈minimum(times(x)) atol=6 * μ
    @test maximum(times(y))≈maximum(times(x)) atol=0.01 * N
end

# @testset "DiffEqBaseExt" begin using DifferentialEquations f(u, p, t) = 1.01 * u u0 = 1 /
#     2 tspan = (0.0, 1.0) prob = ODEProblem(f, u0, tspan, saveat=0.1) sol = solve(prob)

#     x = TimeSeries(sol)
# end

@testset "GeneralizedPhaseExt" begin
    x = bandpass(colorednoise(0.01:0.01:10), (10, 15))
    X = cat(Var(1:10), [bandpass(colorednoise(0.1:0.1:100), (0.1, 0.5)) for _ in 1:10]...)
    _ϕ = @test_nowarn _generalized_phase(x)
    ϕ = @test_nowarn _generalized_phase(X)

    x = set(x, 𝑡 => lookup(x, 𝑡).data * u"s")
    X = set(X, 𝑡 => lookup(X, 𝑡).data * u"s")

    ϕ = @test_nowarn _generalized_phase(x)
    ϕ = @test_nowarn _generalized_phase(X)
end

@testset "ComplexityMeasuresExt" begin
    μ = [1.0, -4.0]
    σ = [2.0, 2.0]
    𝒩 = MvNormal(μ, LinearAlgebra.Diagonal(map(abs2, σ)))
    N = 500
    D = Timeseries(1:N, 1:2, hcat(sort([rand(𝒩) for i in 1:N])...)')
    p = probabilities(NaiveKernel(1.5), StateSpaceSet(D))

    ComplexityMeasures.entropy(Shannon(), ValueBinning(RectangularBinning(100)),
                               StateSpaceSet(D))
end

@testset "Upsampling" begin
    x = TimeSeries(0.1:0.1:10, Var(1:100), randn(100, 100))
    itp = TimeseriesTools.interpolate(x)
    y = itp(dims(x)...)
    @test x ≈ y
    z = @test_nowarn upsample(x, 2)
    @test length(dims(z, 1)) == length(dims(z, 2)) == 199
end
