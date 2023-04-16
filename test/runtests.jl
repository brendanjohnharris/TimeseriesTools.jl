using Unitful
using FFTW
using TimeseriesTools
using Test

@testset "TimeseriesTools.Vl" begin
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

@testset "Spectra" begin
    # Define a test time series
    fs = 1000
    t = range(0, stop=1, length=fs+1)
    x = sin.(2 * π * 50 * t) + sin.(2 * π * 100 * t)
    ts = x = TimeSeries(t, x)
    f_min = fs/100
    Pxx = powerspectrum(ts, f_min)
    @test Pxx isa RegularSpectrum

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
        peaks = findall(x -> x > maximum(Pxx) / 1.5, Pxx)
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
    ts = (-50+0.0005:0.0005:50)*u"s"
    f = rfftfreq(length(ts), 1/step(ts))
    x = 4.2u"V".*sin.(2 * π * 50u"Hz" * ts) .+ 3.1u"V".*sin.(2 * π * 100u"Hz" * ts)
    x = TimeSeries(ts, x)
    S = energyspectrum(x, 0.0)
    P = powerspectrum(x, 0.0)

    @test sum(x.^2).*samplingperiod(x) ≈ sum(S).*step(dims(S, Freq))*2
    @test sum(x.^2).*samplingperiod(x)./duration(x) ≈ sum(S).*step(dims(S, Freq))./duration(x)*2
    @test unit(eltype(S)) == u"V^2*s^2" # Units of energy spectrum

    peaks = findall(x -> x > maximum(S) / 3, S)
    peakfs = f[peaks]
    peakamps = S[peaks]
    @test ustrip.(peakfs) ≈ [50, 100] rtol=1e-3
    @test first(peakamps)/last(peakamps) ≈ 4.2.^2/3.1.^2 rtol=1e-3


    x = exp.(-ustrip.(ts).^2)
    x = TimeSeries(ts, x*u"V")
    ℱ = sqrt(π).*exp.(-π^2 .* ustrip.(f).^2)
    _S = abs.(ℱ).^2*u"V^2*s^2"
    S = energyspectrum(x, 10.0)
    sum(x.^2).*samplingperiod(x)
    @test sum(_S).*step(f) ≈ sum(S).*step(dims(S, Freq)) rtol=0.05
end
