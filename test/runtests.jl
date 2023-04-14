using Unitful
using TimeseriesTools
using Test

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

@testset "Spectra" begin
    # Define a test time series
    fs = 1000
    t = range(0, stop=1, length=fs+1)
    x = sin.(2 * π * 50 * t) + sin.(2 * π * 100 * t)
    ts = x = TimeSeries(t, x)
    f_min = fs/100
    Pxx = spectrum(ts, f_min)
    @test Pxx isa RegularSpectrum
    @test length(Pxx) == expected_nfft ÷ 2 + 1
    freqs, Pxx_values = coordinates(Pxx), Pxx.data
    peaks = findall(x -> x > maximum(Pxx_values) / 2, Pxx_values)
    @test freqs[peaks] ≈ [50, 100] rtol=1e-2

    X = hcat(ts, ts)
    mts = DimArray(X, (Ti(t), Var(:)))
    Pxx_mts = spectrum(mts)
    @test Pxx_mts isa MultivariateDimArray{<:Number, 2, Tuple{FreqDim, VarDim}, DimensionalData.Categorical{DimArray}}

    for i in 1:size(Pxx_mts, Var)
        Pxx_values = Pxx_mts[:, i]
        freqs, Pxx_values = coordinates(Pxx_mts, dims=Freq), Pxx_values.data
        peaks = findall(x -> x > maximum(Pxx_values) / 2, Pxx_values)
        @test freqs[peaks] ≈ [50, 100] rtol=1e-2
    end
end

@testset "Unitful" begin
    ts = (1:1000)u"s"
    x = @test_nowarn TimeSeries(ts, randn(1000))
    @test TimeSeries(ustrip(ts), collect(x), u"s") == x
    @test x isa AbstractTimeSeries
    @test x isa UnitfulTimeSeries
    @test x isa RegularTimeSeries
    @test x isa UnivariateTimeSeries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1/step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test x[At(dims(x, Ti)[1:10])] == x[1:10]
    @test spectrum(x)
end

@testset "Twice unitful" begin
    ts = (-50+0.001:0.001:50)*u"s"
    f = rfftfreq(length(ts), 1/step(ts))
    x = 4.2u"J".*sin.(2 * π * 50u"Hz" * ts) .+ 3.1u"J".*sin.(2 * π * 100u"Hz" * ts)
    x = TimeSeries(ts, x)
    S = spectrum(x, 0)

    # Why do we need to strip the units from the signal?
    @test sum(x.^2).*ustrip(samplingperiod(x)) == sum(S)./step(dims(S, Freq))*2
    @test unit(eltype(S)) == u"J^2/s"

    peaks = findall(x -> x > maximum(S) / 3, S)
    peakfs = f[peaks]
    peakamps = S[peaks]
    @test ustrip.(peakfs) ≈ [50, 100] rtol=1e-3
    @test first(peakamps)/last(peakamps) ≈ 4.2.^2/3.1.^2 rtol=1e-3


    x = exp.(-ustrip(ts).^2)
    x = TimeSeries(ts, x*u"J")
    ℱ = sqrt(π).*exp.(-π^2 .* ustrip(f).^2)
    _S = abs.(ℱ).^2*u"J^2"
    S = spectrum(x, 0)
    sum(x.^2).*samplingperiod(x)
    @test sum(_S)./ustrip(step(f)) == sum(S)./step(dims(S, Freq))
end
