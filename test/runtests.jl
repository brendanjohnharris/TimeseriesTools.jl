using Test
using Unitful
import Unitful.unit
using FFTW
using CairoMakie
using DSP
using ContinuousWavelets
using StatsBase
using TimeseriesSurrogates
using IntervalSets

using TimeseriesTools
import TimeseriesTools: TimeSeries, name, rectifytime

using Documenter
using ImageMagick
using BenchmarkTools
using Foresight

@testset "TimeseriesTools.jl" begin
    ts = 1:100
    x = @test_nowarn TimeSeries(ts, randn(100))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa UnivariateTimeSeries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1 / step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test Interval(x) == first(extrema(ts)) .. last(extrema(ts))
    @test x[Ti(1 .. 10)] == x[1:10]
    @test all(x[Ti(At(1:10))] .== x[1:10])
    # @test x[Ti(At(1:10))] != x[1:10]
end

@testset "Multivariate time series" begin
    ts = 1:100
    x = @test_nowarn TimeSeries(ts, 1:5, randn(100, 5))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa MultivariateTimeSeries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1 / step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test Interval(x) == first(extrema(ts)) .. last(extrema(ts))
    @test x[Ti(1 .. 10)] == x[1:10, :]
end

@testset "Multidimensional time series" begin
    x = @test_nowarn TimeSeries(Ti(1:100), X(1:10), randn(100, 10))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa MultidimensionalTimeSeries

    x = @test_nowarn TimeSeries(Ti(1:100), X(1:10), Y(1:10), randn(100, 10, 10))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa MultidimensionalTimeSeries

    x = @test_nowarn TimeSeries(Ti(1:100), X(randn(10) |> sort), Y(1:10),
                                randn(100, 10, 10))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test !(x isa MultidimensionalTimeSeries)

    x = @test_nowarn TimeSeries(Ti(sort(randn(100))), randn(100))
    @test x isa AbstractTimeSeries
    @test !(x isa RegularTimeSeries)
    @test !(x isa MultidimensionalTimeSeries)
end

@testset "Makie" begin
    x = TimeSeries(0.01:0.01:10, randn(1000))
    p = @test_nowarn plot(x)
    @test p.plot isa Lines
    @test 10 ≤ p.axis.finallimits.val.widths[1] < 12
    x = TimeSeries(0.01:0.01:10, 1:2, randn(1000, 2))
    p = @test_nowarn plot(x)
    @test p.plot isa Heatmap
    @test 10 ≤ p.axis.finallimits.val.widths[1] < 12
    @test 2 ≤ p.axis.finallimits.val.widths[2] < 3
end

@testset "Spectra" begin
    # Define a test time series
    fs = 1000
    t = range(0, stop = 1, length = fs + 1)
    x = 0.8 .* sin.(2 * π * 50 * t) + 1.1 .* sin.(2 * π * 100 * t)
    ts = x = TimeseriesTools.TimeSeries(t, x)
    f_min = fs / 100
    Pxx = powerspectrum(ts, f_min)
    @test Pxx isa RegularSpectrum

    # Plotting
    p = @test_nowarn lines(Pxx)

    freqs = dims(Pxx, Freq)
    peaks = findall(x -> x > maximum(Pxx) / 2, Pxx)
    @test collect(freqs[peaks])≈[50.0, 100.0] rtol=1e-2

    X = hcat(ts, ts)
    mts = DimArray(X, (Ti(t), Var(:)))
    Pxx_mts = powerspectrum(mts, f_min)
    @test Pxx_mts isa MultivariateSpectrum
    @test Pxx_mts[:, 1] == Pxx_mts[:, 2] == Pxx

    for i in 1:size(Pxx_mts, 2)
        Pxx = Pxx_mts[:, i]
        freqs = dims(Pxx, Freq)
        peaks = findall(x -> x > maximum(Pxx) / 2, Pxx)
        @test collect(freqs[peaks])≈[50.0, 100.0] rtol=1e-2
    end

    #  !!!Test padding
    fs = 1000
    t = range(0, stop = 1, length = fs + 1)
    x = 0.8 .* sin.(2 * π * 50 * t) + 1.1 .* sin.(2 * π * 100 * t)
    ts = x = TimeseriesTools.TimeSeries(t, x)
    f_min = fs / 100
    Pa = powerspectrum(ts, f_min; padding = 0)
    Pb = powerspectrum(ts, f_min / 10; padding = 100)
    @test Pb isa RegularSpectrum

    freqs = dims(Pb, Freq)
    peaks = findall(x -> x > maximum(Pb) / 2, Pb)
    @test collect(freqs[peaks])≈[50.0, 100.0] rtol=1e-2

    # @test 2 * sum(energyspectrum(x) .^ 2) .= sum(x .^ 2)
    @test sum(x .^ 2) .* samplingperiod(x)≈sum(Pa) .* step(TimeseriesTools.freqs(Pa)) * 2 rtol=1e-3
    @test sum(x .^ 2) .* samplingperiod(x)≈sum(Pb) .* step(TimeseriesTools.freqs(Pb)) * 2 rtol=1e-5
    # # Plotting
    # f = Figure()
    # ax = Axis(f[1, 1])
    # @test_nowarn lines!(ax, TimeseriesTools.freqs(Pa), Pa)
    # @test_nowarn lines!(ax, TimeseriesTools.freqs(Pb), Pb)
    # save("tmp.pdf", f)
end

@testset "Unitful" begin
    ts = (1:1000)u"s"
    x = @test_nowarn TimeSeries(ts, randn(1000))
    @test TimeSeries(ustrip.(ts), collect(x), u"s") == x
    @test x isa AbstractTimeSeries
    @test x isa UnitfulTimeSeries
    @test x isa RegularTimeSeries
    @test x isa UnivariateTimeSeries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1 / step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test x[Ti(1u"s" .. 10u"s")] == x[1:10]
    @test x[Ti = 1:10] == x[1:10]
    @test_nowarn spectrum(x, 0.1)
end

@testset "Twice unitful" begin
    ts = ((-100 + 0.01):0.0005:100) * u"s"
    f = rfftfreq(length(ts), 1 / step(ts))
    x = 4.2u"V" .* sin.(2 * π * 50u"Hz" * ts) .+ 3.1u"V" .* sin.(2 * π * 100u"Hz" * ts)
    x = TimeSeries(ts, x)
    S = energyspectrum(x, 0.0)
    P = powerspectrum(x, 0.0)

    @test sum(x .^ 2) .* samplingperiod(x) ≈ sum(S) .* step(dims(S, Freq)) * 2
    @test sum(x .^ 2) .* samplingperiod(x) ./ duration(x) ≈
          sum(S) .* step(dims(S, Freq)) ./ duration(x) * 2
    @test unit(eltype(S)) == u"V^2*s^2" # Units of energy spectrum

    peaks = findall(x -> x > maximum(P) / 3, P)
    peakfs = f[peaks]
    peakamps = P[peaks]
    @test all(round.(ustrip.(peakfs)) .∈ ([50, 100],))
    @test first(peakamps) / last(peakamps)≈4.2 / 3.1 rtol=1e-1

    x = exp.(-ustrip.(ts) .^ 2)
    x = TimeSeries(ts, x * u"V")
    ℱ = sqrt(π) .* exp.(-π^2 .* ustrip.(f) .^ 2)
    _S = abs.(ℱ) .^ 2 * u"V^2*s^2"
    S = energyspectrum(x, 0.0)
    @test sum(_S) .* step(f)≈sum(S) .* step(dims(S, Freq)) rtol=0.05

    lines(ustrip.(f), ustrip.(_S), axis = (; limits = ((0, 1), (0, 4))))
    plot!(collect(ustrip.(dims(S, Freq))), collect(ustrip.(S)))
    current_figure()
end

# @testset "Doctests" begin
#     using TimeseriesTools
#     using Unitful
#     using Documenter

#     DocMeta.setdocmeta!(TimeseriesTools, :DocTestSetup, :(using Unitful, TimeseriesTools); recursive=true)

#     doctest(TimeseriesTools)
# end

@testset "Readme" begin
    using TimeseriesTools, CairoMakie, Unitful
    import TimeseriesTools.TimeSeries # or TS

    t = 0.005:0.005:1e5
    x = colorednoise(t, u"s") * u"V"

    # Plot the time series
    f = Figure(; size = (720, 480))
    ax = Axis(f[1, 1])
    @test_nowarn plot!(ax, x[1:10000])
    save("./timeseries.png", f; px_per_unit = 3)

    # Calculate the power spectrum
    S = _powerspectrum(x, 0.0001)
    f = Figure(; size = (720, 480))
    ax = Axis(f[1, 1])
    @test_nowarn plot!(ax, S, linewidth = 1)
    @test_nowarn save("./powerspectrum.png", f; px_per_unit = 3)

    # Shadows
    x = loadtimeseries("./test_timeseries.csv")

    f = Figure(; size = (500, 480))
    ax = Axis3(f[1, 1])
    trajectory!(ax, collect.(eachcol(x))...; colormap = :turbo, linewidth = 0.1)
    ax.xlabelvisible = ax.ylabelvisible = ax.zlabelvisible = ax.xticksvisible = ax.yticksvisible = ax.zticksvisible = ax.xticklabelsvisible = ax.yticklabelsvisible = ax.zticklabelsvisible = false
    ax.azimuth = ax.azimuth[] + 0.25
    ax.elevation = ax.elevation[] + 0.25
    shadows!(ax, collect.(eachcol(x))...; color = (:slategray, 0.5), linewidth = 0.05)
    save("./shadows.png", f; px_per_unit = 3)
end

@testset "Readme_dark" begin
    using CairoMakie, TimeseriesTools, Unitful
    import TimeseriesTools.TimeSeries # or TS
    set_theme!(foresight(:dark, :transparent))

    t = 0.005:0.005:1e5
    x = colorednoise(t, u"s") * u"V"

    # Plot the time series
    f = Figure(; size = (720, 480))
    ax = Axis(f[1, 1])
    @test_nowarn plot!(ax, x[1:10000])
    save("./timeseries_dark.png", f; px_per_unit = 3)

    # Calculate the power spectrum
    S = _powerspectrum(x, 0.0001)
    f = Figure(; size = (720, 480))
    ax = Axis(f[1, 1])
    @test_nowarn plot!(ax, S, linewidth = 1)
    @test_nowarn save("./powerspectrum_dark.png", f; px_per_unit = 3)

    # Shadows
    x = loadtimeseries("./test_timeseries.csv")

    f = Figure(; size = (500, 480))
    ax = Axis3(f[1, 1])
    trajectory!(ax, collect.(eachcol(x))...; colormap = :turbo, linewidth = 0.1)
    ax.xlabelvisible = ax.ylabelvisible = ax.zlabelvisible = ax.xticksvisible = ax.yticksvisible = ax.zticksvisible = ax.xticklabelsvisible = ax.yticklabelsvisible = ax.zticklabelsvisible = false
    ax.azimuth = ax.azimuth[] + 0.25
    ax.elevation = ax.elevation[] + 0.25
    shadows!(ax, collect.(eachcol(x))...; color = (:white, 0.5), linewidth = 0.05)
    save("./shadows_dark.png", f; px_per_unit = 3)
end

@testset "Unit Power" begin
    N = UnitPower
    _X = TimeSeries(0.01:0.01:1, rand(100))
    X = copy(_X)
    T = fit(N, X)
    Y = normalize(X, T)
    @test sum(Y .^ 2) / duration(Y) ≈ 1
    @test !isnothing(T.p)
    @test denormalize(Y, T) ≈ X
    @test_nowarn normalize!(X, T)
    @test X == Y
    @test_nowarn denormalize!(Y, T)
    @test all(Y .≈ _X)

    _X = TimeseriesTools.unitfultimeseries(X, u"s") * u"V"
    X = copy(_X)
    T = fit(N, X)
    Y = normalize(X, T)
    @test ustrip(sum(Y .^ 2) / duration(Y)) ≈ 1
    @test !isnothing(T.p)
    @test_throws "Denormalization of unitful arrays currently not supported" denormalize(Y,
                                                                                         T)
    X = @test_nowarn normalize(X, T)
    @test X == Y
    Y = @test_throws "Denormalization of unitful arrays currently not supported" denormalize(Y,
                                                                                             T)
    # @test all(Y .≈ _X)
end

@testset "IO" begin
    x = TimeSeries(0.001:0.001:1, 1:3, rand(1000, 3); metadata = Dict(:a => :test),
                   name = "name")

    f = tempname() * ".jld2"
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x == _x

    f = tempname() * ".csv"
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x ≈ _x

    x = x[:, 1]
    savetimeseries(f, x)
    _x = @test_logs (:warn, "Cannot load refdims yet") loadtimeseries(f)
    @test refdims(_x) == ()
    @test all(x .≈ _x)

    x = TimeSeries(0.001:0.001:1, 1:3, rand(1000, 3); metadata = Dict(:a => :test))
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x ≈ _x

    # Currently not the greatest way of handling non-serializable metadata
    x = TimeSeries(0.001:0.001:1, 1:3, rand(1000, 3);
                   metadata = Dict(:a => DimensionalData.NoName())) # Something that can't be serialized
    @test_logs (:warn, ErrorException("Cannot serialize type DimensionalData.NoName")) savetimeseries(f,
                                                                                                      x)
    _x = loadtimeseries(f)
    @test metadata(_x) == DimensionalData.Dimensions.LookupArrays.NoMetadata()
    @test x ≈ _x

    x = TimeSeries(0.001:0.001:1, 1:3, rand(1000, 3); name = TimeSeries) # Something that can't be serialized
    @test_logs (:warn, ErrorException("Cannot serialize type typeof(TimeSeries)")) savetimeseries(f,
                                                                                                  x)
    _x = loadtimeseries(f)
    @test name(_x) == DimensionalData.NoName()
    @test x ≈ _x

    x = TimeSeries(0.001:0.001:1, [TimeSeries, TimeSeries, TimeSeries], rand(1000, 3))
    @test_logs (:warn, ErrorException("Cannot serialize type typeof(TimeSeries)")) savetimeseries(f,
                                                                                                  x)
    _x = loadtimeseries(f)
    @test x ≈ _x

    x = TimeSeries(0.001:0.001:1, rand(1000))
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x ≈ _x

    x = TimeSeries((0.001:0.001:1) * u"s", 1:3, rand(1000, 3); metadata = Dict(:a => :test),
                   name = "name") * u"V"

    f = tempname() * ".jld2"
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x == _x
end

@testset "Traces" begin
    using CairoMakie, TimeseriesTools, Unitful
    import TimeseriesTools.TimeSeries # or TS

    t = 0.005:0.005:1e4
    x = colorednoise(t, u"s") * u"V"
    X = cat(Var(1:2), x, x .+ 1.0 * u"V", dims = 2)

    # Calculate the power spectrum
    S = _powerspectrum(x, 0.0005)[2:end, :]
    f = Figure(; size = (720, 480))
    ax = Axis(f[1, 1], xscale = log10, yscale = log10)
    # x, y, z = collect.(ustrip.(decompose(S)))
    @test_nowarn traces!(ax, S; colormap = :turbo)
end

@testset "Spectrum plot" begin
    using DSP
    using CairoMakie, TimeseriesTools, Unitful
    import TimeseriesTools.TimeSeries # or TS

    t = 0.005:0.005:1e4
    x = colorednoise(t, u"s") * u"V"
    X = cat(Var(1:2), x, x .+ 1.0 * u"V", dims = 2)

    # Calculate the power spectrum
    S = _powerspectrum(x, 0.0005)[2:end, :]
    f = Figure(; size = (720, 480))
    ax = Axis(f[1, 1], xscale = log10, yscale = log10)
    @test_nowarn spectrumplot!(ax, S, linewidth = 2)

    # * Test peaks
    x = bandpass(x, (0.1u"Hz", 0.2u"Hz"))
    S = powerspectrum(x, 0.0005)
    spectrumplot(S; peaks = true)
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
    @test_nowarn plot(x[Ti(1:10000)])
    plot(x̂[Ti(1500:(length(t) * 5))])

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
    @test ≈(pac, autocor(x[Ti(1:10000)] |> collect, [10])[1]; rtol = 1e-2)
    # @test pac - autocor(x̂[Ti(1:10000)] |> collect, [10])[1] >
    #   pac - autocor(x[Ti(1:10000)] |> collect, [10])[1]
end

@testset "Interlace" begin
    x = TimeSeries(0:0.1:1, randn(11))
    y = TimeSeries(0.05:0.1:1, randn(10))
    z = @test_nowarn interlace(x, y)
    @test all(collect(times(z)) .== 0.0:0.05:1.0)
end

@testset "Buffer" begin
    N = 10
    x = TimeSeries(0.1:0.1:10, randn(100))
    y = @test_nowarn buffer(x, 10)
    @test length(y) == N
    @test y[1] == x[1:(length(x) ÷ N)]
    @test cat(y..., dims = Ti) == x[1:((length(x) ÷ N) * N)]

    y = @test_nowarn buffer(x, 10, 0; discard = false)
    @test cat(y..., dims = Ti) == x

    y = @test_nowarn buffer(x, 10, N ÷ 2)
    @test length(y) == 2 * N - 1

    x = TimeSeries(0:0.1:10, 1:10, randn(101, 10))
    y = buffer(x, 10)
    @test length(y) == 10

    x = TimeSeries(0.1:0.1:10, randn(100))
    y = @test_nowarn window(x, 2, 1)
    @test all(length.(y) .== 2)
    y = @test_nowarn delayembed(x, 2, 1, 1)
    y = @test_nowarn delayembed(x, 2, 1, 2)
    @test length(y) == length(x)
    @test samplingperiod(y) == 2 * samplingperiod(x)
    y = @test_nowarn delayembed(x, 2, 2, 1)
end

@testset "Spike FFT" begin
    ts = 0:0.01:100
    t = [abs(_t - round(_t)) < 0.05 ? 1 : 0 for _t in ts][1:(end - 1)]
    t = findall(t .> 0) ./ 100 # Should have a period of 1 second
    t = TimeSeries(t, trues(length(t)))
    @test t isa SpikeTrain

    p = @test_nowarn spikefft(0:0.1:10, t)
    fs = (0.01, 50)
    e = @test_nowarn energyspectrum(t, fs; method = :schild)
    @test (2 * sum(e[2:end]) + e[1]) * fs[1] ≈ sum(t)
    p = @test_nowarn powerspectrum(t, fs; method = :schild)

    # Test this returns an identical result for spikes measured at regular intervals
    x = TimeSeries(ts, zeros(length(ts)))
    x[Ti(At(times(t)))] .= 1.0 / sqrt(samplingperiod(x))
    et = energyspectrum(x, 0.01)
    @test (2 * sum(et[2:end]) + et[1]) * fs[1] ≈ sum(t)

    # if false # Yep works, better, even
    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, et, color = :crimson)
    lines!(ax, e)
    f
    # end

    # Multivariate
    T = hcat(Var(1:4), t, t, t, t)
    P = @test_nowarn powerspectrum(T, fs; method = :schild)
    @test P[:, 1] == p

    # * Autocovariance spectrum
    p = energyspectrum(t, fs; method = stoic(; σ = 0.01))
    @test (2 * sum(p[2:end]) + p[1]) * fs[1] ≈ sum(t)

    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, decompose(et)..., color = :crimson)
    lines!(ax, decompose(p)...)
    f
end

@testset "Spike-time tiling coefficient" begin
    using IntervalSets
    using LinearAlgebra
    using Distributions
    ts = 0:0.01:100
    t = [abs(_t - round(_t)) < 0.05 ? 1 : 0 for _t in ts][1:(end - 1)]
    t = findall(t .> 0) ./ 100 # Should have a period of 1 second
    t = TimeSeries(t, trues(length(t)))
    Δt = 0.025
    c = @test_nowarn sttc(t, t .+ 0.02; Δt)

    t1 = rand(0 .. 10, 200)
    η = 8000.0
    t2 = t1 .+ η .* randn(length(t1)) * Δt
    sort!.([t1, t2])
    @test_nowarn sttc(t1, t2; Δt)

    # Positive semi-definite? Not quite, but almost.
    η = 1
    t1 = findall(rand(Poisson(0.01), 100000) .> 0) ./ 1000
    t2 = findall(rand(Poisson(0.05), 100000) .> 0) ./ 1000
    t3 = vcat(t1 .* η .* randn(length(t1)) * Δt, t2 .* η .* randn(length(t2)) * Δt)
    t4 = vcat(t2 .* η .* randn(length(t2)) * Δt, t1 .* η .* randn(length(t1)) * Δt)
    t5 = vcat(t1 .* η .* randn(length(t1)) * Δt, t1 .* η .* randn(length(t1)) * Δt)
    t6 = vcat(t2 .* η .* randn(length(t2)) * Δt, t2 .* η .* randn(length(t2)) * Δt)
    sort!.([t3, t4, t5, t6])

    ts = [t1, t2, t3, t4, t5, t6]

    is = Iterators.product(ts, ts)
    Λ = [sttc(i...) for i in is]
    λ = eigvals(Λ)
end

@testset "Spike-time overlap-integral coefficient (stoic)" begin
    x = randn(1000) |> sort
    y = randn(1000) |> sort
    σ = 0.25
    Δt = σ * 10

    # * Verify the integral of the product of two gaussians
    G = TimeseriesTools.normal
    t1 = 0.025
    t2 = 0.0
    f(x) = G(t1, σ)(x) * G(t2, σ)(x)
    I1 = sum(f.(-1:0.001:1)) * 0.001
    I2 = G(t1, sqrt(2) * σ)(t2)
    @test I1≈I2 rtol=1e-6
    # Aw yeah

    D = @test_nowarn closeneighbours(x, y; Δt)
    @test stoic(x, y; Δt, σ)≈1.0 rtol=5e-2

    x = y
    @test stoic(x, y; Δt, σ) == 1.0

    x = rand(0 .. 1000, 1000) |> sort
    y = rand(0 .. 1000, 1000) |> sort
    σ = 100
    Δt = σ * 10
    @test stoic(x, y; Δt, σ)≈1.0 rtol=0.02
    σ = 0.001
    Δt = σ * 10
    @test stoic(x, y; σ)≈0.0 atol=1e-2

    @test stoic([0.0], [0.0]; σ = 1, normalize = false) ≈ 0.5 / sqrt(π)
    @test stoic([0.0], [1.0]; σ = 1, normalize = false) ≈ 1 / (2 * exp(1 / 4) * sqrt(π))
    @test stoic([0.0, 10.0], [1.0, 10.0]; σ = 1, normalize = false)≈0.50179 rtol=1e-4

    # * Is it positive semi-definite?
    x = [rand(0 .. 100, 100) |> sort for _ in 1:100]
    [_x .= x[22] for _x in x[round.(Int, rand(1 .. length(x), 20))]]
    [_x .= sort(x[22] .+ 0.01 .* randn(100))
     for _x in x[round.(Int, rand(1 .. length(x), 20))]]
    [_x .= sort(x[22] .+ 0.05 .* randn(100))
     for _x in x[round.(Int, rand(1 .. length(x), 20))]]
    ρ = @test_nowarn pairwise(stoic(; σ = 0.01), x)
    e = eigvals(ρ)
    @test minimum(real.(e)) + 1e-10 > 0.0
    @test all(isapprox.(imag.(e), 0.0; atol = 1e-10))
end

@testset "Stoic spike-train length" begin
    # * Set up independent gamma renewal processes and verify stoic scaling with length vs. kernel width
    Ns = range(start = 100, step = 100, length = 100)
    xs = [gammarenewal(N, 1, 1) for N in Ns]
    ρ = @test_nowarn pairwise(stoic(; σ = 0.01), xs; symmetric = true)
    @test mean(ρ[ρ .!= 1])≈0 atol=0.05
    ρ[ρ .== 1] .= NaN

    f = Figure()
    ax = Axis(f[1, 1]; aspect = 1, xlabel = "N₁", ylabel = "N₂")
    p = heatmap!(ax, Ns, Ns, ρ)
    Colorbar(f[1, 2], p, label = "stoic")
end
@testset "Stoic spike-train fano" begin
    # * Set up independent gamma renewal processes and verify stoic scaling with length vs. kernel width
    θs = range(start = 0.1, step = 0.01, length = 150)
    xs = [gammarenewal(10000, 1, θ) for θ in θs]
    ρ = pairwise(stoic(; σ = 0.01), xs; symmetric = true)
    ρ[ρ .== 1] .= NaN

    f = Figure(size = (720, 360))
    ax = Axis(f[1, 1]; aspect = 1, xlabel = "θ₁", ylabel = "θ₂")
    p = heatmap!(ax, θs, θs, ρ)
    Colorbar(f[1, 2], p, label = "stoic")
    f

    ρ2 = pairwise(sttc(; Δt = 0.03), xs; symmetric = true)
    ρ2[ρ2 .== 1] .= NaN

    ax = Axis(f[1, 3]; aspect = 1, xlabel = "θ₁", ylabel = "θ₂")
    p = heatmap!(ax, θs, θs, abs.(ρ2))
    Colorbar(f[1, 4], p, label = "|sttc|")

    rowsize!(f.layout, 1, Relative(0.6))
    f
end

# @testset "Stoic firing rate and fano factor" begin
#     N = 5000

#     fr = range(start = 0.1, stop = 500, length = 100) # 1 the mean ISI
#     θs = range(start = 0.1, stop = 2, length = 100) # 1 the mean ISI
#     αs = 1 ./ (fr .* θs)

#     ps = Iterators.product(αs, θs)

#     ρ = zeros(length(fr), length(θs))
#     n = 50
#     Threads.@threads for i in 1:n
#         _ρ = map(ps) do (α, θ)
#             x = gammarenewal(N, α, θ)
#             y = gammarenewal(N, α, θ)
#             # sttc(x, y; Δt = 0.025) |> abs
#             stoic(x, y; σ = 0.01) |> abs
#         end
#         ρ .+= _ρ
#     end
#     ρ ./= n # Mean

#     f = Figure()
#     ax = Axis(f[1, 1]; aspect = 1, xlabel = "Mean firing rate", ylabel = "Fano factor")
#     p = heatmap!(ax, fr, θs, ρ)
#     Colorbar(f[1, 2], p, label = "stoic")
#     f
# end

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
    @test all(isa.(dims(S), (Ti, Freq, Var)))

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
        @test all(isa.(dims(S), (Ti, Freq, Var)))

        y = set(x, CuArray(x.data))
        S = @test_nowarn waveletspectrogram(y)
        @test all(isa.(dims(S), (Ti, Freq, Var)))

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
    y = set(x, Ti => surrogate(times(x), RandomJitter(0.1, 0.1)))
    @test y isa SpikeTrain
    @test issorted(times(y))
    @test minimum(times(y))≈minimum(times(x)) atol=0.5
    @test maximum(times(y))≈maximum(times(x)) atol=0.5
    @test x != y
    sur = @test_nowarn surrogenerator(times(x), RandomJitter(0.1, 0.1))
    @test all(copy(sur()) .!= sur())

    # Gamma renewal surrogate
    y = set(x, Ti => surrogate(times(x), GammaRenewal()))
    dt̂ = diff(times(y))
    F̂ = var(dt̂) / mean(dt̂)
    @test y isa SpikeTrain
    @test issorted(times(y))
    @test F̂≈θ rtol=5e-2
    @test mean(dt̂)≈α * F̂ rtol=5e-2
    @test minimum(times(y))≈minimum(times(x)) atol=6 * μ
    @test maximum(times(y))≈maximum(times(x)) atol=0.01 * N
end

@testset "Central differences" begin
    x = colorednoise(0.01:0.01:10)
    X = cat(Var(1:10), [colorednoise(0.1:0.1:100) for _ in 1:10]...)

    dx = @test_nowarn centraldiff(x)
    @test all(dx[2:(end - 1)] .== (x[3:end] - x[1:(end - 2)]) / 2)
    @test times(dx) == times(x)

    dX = @test_nowarn centraldiff(X)
    @test all(dX[2:(end - 1), :] .== (X[3:end, :] - X[1:(end - 2), :]) / 2)
    @test times(dX) == times(X)
    @test dims(dX, Var) == dims(X, Var)

    dX = @test_nowarn centralderiv(X)
    @test all(dX[2:(end - 1), :] .==
              ((X[3:end, :] - X[1:(end - 2), :]) / 2) ./ samplingperiod(X))

    x = @test_nowarn Timeseries(0.1:0.1:1000, sin)
    𝑓 = instantaneousfreq(x)
    @test std(𝑓[2500:(end - 2500)]) < 0.001
    @test mean(𝑓[2500:(end - 2500)])≈1 / 2π rtol=1e-5

    ϕ = analyticphase(x)[1000:(end - 1000)]
    dϕ = @test_nowarn centraldiff(ϕ)
end

@testset "Rectification" begin
    ts = 0.1:0.1:1000
    x = TimeSeries(ts .+ randn(length(ts)) .* 1e-8, sin)
    @test issorted(times(x))
    _x = @test_nowarn rectifytime(x)
    @test all(x .== _x)
    @test ts == times(_x)

    y = TimeSeries(ts .+ randn(length(ts)) .* 1e-8, cos)
    @test issorted(times(y))
    _x, _y = (rectifytime([x, y])...,)

    @test all(x .== _x)
    @test ts == times(_x)
    @test all(y .== _y)
    @test ts == times(_y)

    x = @test_nowarn TimeSeries(Ti(1:100), X((1:10) .+ 1e-9 .* randn(10)), randn(100, 10))
    y = @test_nowarn rectify(x, dims = X)
    @test dims(y, X) == X(1:10)

    x = @test_nowarn TimeSeries(Ti(1:100), X((1:10) .+ 1e-9 .* randn(10)),
                                Y((1:5) .+ 1e-9 .* randn(5)), randn(100, 10, 5))
    y1 = @test_nowarn rectify(x, dims = X)
    y2 = @test_nowarn rectify(x, dims = Y)
    y3 = @test_nowarn rectify(x, dims = [X, Y])
    @test dims(y1, X) == dims(y3, X) == X(1:10)
    @test dims(y2, Y) == dims(y3, Y) == Y(1:5)
    @test dims(y1, Y) == dims(x, Y)
    @test dims(y2, X) == dims(x, X)
end

# @testset "DiffEqBaseExt" begin
#     using DifferentialEquations
#     f(u, p, t) = 1.01 * u
#     u0 = 1 / 2
#     tspan = (0.0, 1.0)
#     prob = ODEProblem(f, u0, tspan, saveat=0.1)
#     sol = solve(prob)

#     x = TimeSeries(sol)
# end
