using Unitful
import Unitful.unit
using FFTW
using CairoMakie
using TimeseriesTools
import TimeseriesTools.TimeSeries
using TimeseriesSurrogates
using Test
using Documenter
using ImageMagick
using Foresight

@testset "TimeseriesTools.jl" begin
    ts = 1:100
    x = @test_nowarn TimeSeries(ts, randn(100))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa UnivariateTimeSeries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1/step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test Interval(x) == first(extrema(ts))..last(extrema(ts))
    @test x[At(dims(x, Ti)[1:10])] == x[1:10]
end

@testset "Multivariate time series" begin
    ts = 1:100
    x = @test_nowarn TimeSeries(ts, 1:5, randn(100, 5))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa MultivariateTimeSeries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1/step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test Interval(x) == first(extrema(ts))..last(extrema(ts))
    @test x[At(dims(x, Ti)[1:10]), :] == x[1:10, :]
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
    t = range(0, stop=1, length=fs+1)
    x = 0.8.*sin.(2 * π * 50 * t) + 1.1.*sin.(2 * π * 100 * t)
    ts = x = TimeseriesTools.TimeSeries(t, x)
    f_min = fs/100
    Pxx = powerspectrum(ts, f_min)
    @test Pxx isa RegularSpectrum

    # Plotting
    p = @test_nowarn lines(Pxx)


    freqs = dims(Pxx, Freq)
    peaks = findall(x -> x > maximum(Pxx) / 2, Pxx)
    @test collect(freqs[peaks]) ≈ [50.0, 100.0] rtol=1e-2

    X = hcat(ts, ts)
    mts = DimArray(X, (Ti(t), Var(:)))
    Pxx_mts = powerspectrum(mts, f_min)
    @test Pxx_mts isa MultivariateSpectrum
    @test Pxx_mts[:, 1] == Pxx_mts[:, 2] == Pxx

    for i in 1:size(Pxx_mts, 2)
        Pxx = Pxx_mts[:, i]
        freqs = dims(Pxx, Freq)
        peaks = findall(x -> x > maximum(Pxx) / 2, Pxx)
        @test collect(freqs[peaks]) ≈ [50.0, 100.0] rtol=1e-2
    end
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
    @test samplingrate(x) == 1/step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test x[At(dims(x, Ti)[1:10])] == x[1:10]
    @test_nowarn spectrum(x)
end

@testset "Twice unitful" begin
    ts = (-100+0.01:0.0005:100)*u"s"
    f = rfftfreq(length(ts), 1/step(ts))
    x = 4.2u"V".*sin.(2 * π * 50u"Hz" * ts) .+ 3.1u"V".*sin.(2 * π * 100u"Hz" * ts)
    x = TimeSeries(ts, x)
    S = energyspectrum(x, 0.0)
    P = powerspectrum(x, 0.0)

    @test sum(x.^2).*samplingperiod(x) ≈ sum(S).*step(dims(S, Freq))*2
    @test sum(x.^2).*samplingperiod(x)./duration(x) ≈ sum(S).*step(dims(S, Freq))./duration(x)*2
    @test unit(eltype(S)) == u"V^2*s^2" # Units of energy spectrum

    peaks = findall(x -> x > maximum(P) / 3, P)
    peakfs = f[peaks]
    peakamps = P[peaks]
    @test all(round.(ustrip.(peakfs)) .∈ ([50, 100],))
    @test first(peakamps)/last(peakamps) ≈ 4.2/3.1 rtol=1e-1


    x = exp.(-ustrip.(ts).^2)
    x = TimeSeries(ts, x*u"V")
    ℱ = sqrt(π).*exp.(-π^2 .* ustrip.(f).^2)
    _S = abs.(ℱ).^2*u"V^2*s^2"
    S = energyspectrum(x, 0.0)
    @test sum(_S).*step(f) ≈ sum(S).*step(dims(S, Freq)) rtol=0.05

    lines(ustrip.(f), ustrip.(_S), axis=(; limits=((0, 1),(0, 4))))
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
    x = colorednoise(t, u"s")*u"V"

    # Plot the time series
    f = Figure(; resolution=(720, 480))
    ax = Axis(f[1, 1])
    @test_nowarn plot!(ax, x[1:10000])
    save("./timeseries.png", f)

    # Calculate the power spectrum
    S = _powerspectrum(x, 0.0001)
    f = Figure(; resolution=(720, 480))
    ax = Axis(f[1, 1])
    @test_nowarn plot!(ax, S)
    @test_nowarn save("./powerspectrum.png", f)
end

@testset "Readme_dark" begin
    using CairoMakie, TimeseriesTools, Unitful
    import TimeseriesTools.TimeSeries # or TS
    set_theme!(foresight(:dark, :transparent))

    t = 0.005:0.005:1e5
    x = colorednoise(t, u"s")*u"V"

    # Plot the time series
    f = Figure(; resolution=(720, 480))
    ax = Axis(f[1, 1])
    @test_nowarn plot!(ax, x[1:10000])
    save("./timeseries_dark.png", f)

    # Calculate the power spectrum
    S = _powerspectrum(x, 0.0001)
    f = Figure(; resolution=(720, 480))
    ax = Axis(f[1, 1])
    @test_nowarn plot!(ax, S)
    @test_nowarn save("./powerspectrum_dark.png", f)
end

@testset "Unit Power" begin
    N = UnitPower
    _X = TimeSeries(0.01:0.01:1, rand(100))
    X = copy(_X)
    T = fit(N, X)
    Y = normalize(X, T)
    @test sum(Y.^2)/duration(Y) ≈ 1
    @test !isnothing(T.p)
    @test denormalize(Y, T) ≈ X
    @test_nowarn normalize!(X, T)
    @test X == Y
    @test_nowarn denormalize!(Y, T)
    @test all(Y .≈ _X)

    _X = TimeseriesTools.unitfultimeseries(X, u"s")*u"V"
    X = copy(_X)
    T = fit(N, X)
    Y = normalize(X, T)
    @test ustrip(sum(Y.^2)/duration(Y)) ≈ 1
    @test !isnothing(T.p)
    @test_throws "Denormalization of unitful arrays currently not supported" denormalize(Y, T)
    X = @test_nowarn normalize(X, T)
    @test X == Y
    Y = @test_throws "Denormalization of unitful arrays currently not supported" denormalize(Y, T)
    # @test all(Y .≈ _X)
end


@testset "TimeseriesSurrogatesExt" begin
    x = TimeSeries(0.001:0.001:1, rand(1000))
    @test_nowarn s = surrogate(x, FT())
end

@testset "IO" begin
    x = TimeSeries(0.001:0.001:1, 1:3, rand(1000, 3); metadata=Dict(:a=>:test), name="name")

    f = tempname()*".jld2"
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x == _x

    f = tempname()*".csv"
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x ≈ _x

    x = x[:, 1]
    savetimeseries(f, x)
    _x = @test_logs (:warn, "Cannot load refdims yet") loadtimeseries(f)
    @test refdims(_x) == ()
    @test all(x .≈ _x)

    x = TimeSeries(0.001:0.001:1, 1:3, rand(1000, 3); metadata=Dict(:a=>:test))
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    display(_x)
    @test x ≈ _x

    # Currently not the greatest way of handling non-serializable metadata
    x = TimeSeries(0.001:0.001:1, 1:3, rand(1000, 3); metadata=Dict(:a=>DimensionalData.NoName())) # Something that can't be serialized
    @test_logs (:warn, ErrorException("Cannot serialize type DimensionalData.NoName")) savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test metadata(_x) == DimensionalData.Dimensions.LookupArrays.NoMetadata()
    @test x ≈ _x

    x = TimeSeries(0.001:0.001:1, 1:3, rand(1000, 3); name=TimeSeries) # Something that can't be serialized
    @test_logs (:warn, ErrorException("Cannot serialize type typeof(TimeSeries)")) savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test name(_x) == DimensionalData.NoName()
    @test x ≈ _x

    x = TimeSeries(0.001:0.001:1, [TimeSeries, TimeSeries, TimeSeries], rand(1000, 3))
    @test_logs (:warn, ErrorException("Cannot serialize type typeof(TimeSeries)")) savetimeseries(f, x)
    _x = loadtimeseries(f)
    display(_x)
    @test x ≈ _x

    x = TimeSeries(0.001:0.001:1, rand(1000))
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x ≈ _x


    x = TimeSeries((0.001:0.001:1)*u"s", 1:3, rand(1000, 3); metadata=Dict(:a=>:test), name="name")*u"V"

    f = tempname()*".jld2"
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x == _x
end


@testset "Traces" begin
    using CairoMakie, TimeseriesTools, Unitful
    import TimeseriesTools.TimeSeries # or TS

    t = 0.005:0.005:1e4
    x = colorednoise(t, u"s")*u"V"
    X = cat(Var(1:2), x, x.+1.0*u"V", dims=2)

    # Calculate the power spectrum
    S = _powerspectrum(x, 0.0005)
    f = Figure(; resolution=(720, 480))
    ax = Axis(f[1, 1])
    x, y, z = collect.(ustrip.(decompose(X)))
    @test_nowarn traces!(ax, x, y, z)
    @test_nowarn save("./powerspectrum.png", f)
end
