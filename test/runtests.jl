using TestItems
using TestItemRunner

@run_package_tests

@testitem "Dates" begin
    using Dates, Unitful
    x = 1:100
    t = DateTime(1901):Year(1):DateTime(2000)
    y = @test_nowarn Timeseries(x, t)
    @test y isa RegularTimeseries
    @test samplingperiod(y) == Year(1)
    @test times(y) == t
    @test duration(y) == last(t) - first(t)
    @test unit(y) == NoUnits
end

@testitem "Spectra" begin
    using Unitful, CairoMakie
    # Define a test time series
    fs = 1000
    t = range(0, stop = 1, length = fs + 1)
    x = 0.8 .* sin.(2 * π * 50 * t) + 1.1 .* sin.(2 * π * 100 * t)
    ts = x = TimeseriesTools.Timeseries(x, t)
    f_min = fs / 100
    Pxx = powerspectrum(ts, f_min)
    @test Pxx isa RegularSpectrum

    xu = set(x, 𝑡 => ustripall(t) * u"s")
    Pxu = @test_nowarn powerspectrum(xu, f_min)
    @test unit(eltype(Pxu)) == u"s"
    @test unit(eltype(lookup(Pxu, 1))) == u"s^-1"
    @test all(ustripall(Pxu) .≈ Pxx)

    @test_throws "DomainError" powerspectrum(x, 1e-6)

    # Plotting
    p = @test_nowarn lines(Pxx)

    freqs = dims(Pxx, 𝑓)
    peaks = findall(x -> x > maximum(Pxx) / 2, Pxx)
    @test collect(freqs[peaks])≈[50.0, 100.0] rtol=1e-2

    xx = hcat(ts, ts)
    mts = ToolsArray(xx, (𝑡(t), Var(:)))
    Pxx_mts = powerspectrum(mts, f_min)
    @test Pxx_mts isa MultivariateSpectrum
    @test Pxx_mts[:, 1] == Pxx_mts[:, 2] == Pxx

    for i in axes(Pxx_mts, 2)
        Pxx = Pxx_mts[:, i]
        freqs = dims(Pxx, 𝑓)
        peaks = findall(x -> x > maximum(Pxx) / 2, Pxx)
        @test collect(freqs[peaks])≈[50.0, 100.0] rtol=1e-2
    end

    #  !!!Test padding
    fs = 1000
    t = range(0, stop = 1, length = fs + 1)
    x = 0.8 .* sin.(2 * π * 50 * t) + 1.1 .* sin.(2 * π * 100 * t)
    ts = x = TimeseriesTools.Timeseries(x, t)
    f_min = fs / 100
    Pa = powerspectrum(ts, f_min; padding = 0)
    Pb = powerspectrum(ts, f_min / 10; padding = 100)
    @test Pb isa RegularSpectrum

    freqs = dims(Pb, 𝑓)
    peaks = findall(x -> x > maximum(Pb) / 2, Pb)
    @test collect(freqs[peaks])≈[50.0, 100.0] rtol=1e-2

    # @test 2 * sum(energyspectrum(x) .^ 2) .= sum(x .^ 2)
    @test sum(x .^ 2) .* samplingperiod(x)≈sum(Pa) .* step(TimeseriesTools.freqs(Pa)) * 2 rtol=1e-3
    @test sum(x .^ 2) .* samplingperiod(x)≈sum(Pb) .* step(TimeseriesTools.freqs(Pb)) * 2 rtol=1e-5
    # # Plotting
    # f = Figure() ax = Axis(f[1, 1]) @test_nowarn lines!(ax, TimeseriesTools.freqs(Pa), Pa)
    # @test_nowarn lines!(ax, TimeseriesTools.freqs(Pb), Pb) save("tmp.pdf", f)
end

include("Types.jl")
include("Utils.jl")
include("IO.jl")
include("Unitful.jl")
include("SpikeTrains.jl")
include("Operators.jl")
include("TimeseriesMakie.jl")
include("TimeseriesSurrogatesExt.jl")
include("Extensions.jl")
