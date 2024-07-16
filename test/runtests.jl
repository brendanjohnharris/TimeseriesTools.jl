using Test
using Term
using Unitful
import Unitful.unit
using FFTW
using CairoMakie
using DSP
using Dates
using ContinuousWavelets
using StatsBase
using TimeseriesSurrogates
using IntervalSets
using Dierckx
using GeneralizedPhase
using Autocorrelations
using MeanSquaredDisplacement

using TimeseriesTools
import TimeseriesTools: TimeSeries, name, rectifytime,
                        leftdiff, rightdiff,
                        NDAAFT, NDIAAFT, MVFT

using Documenter
using ImageMagick
using BenchmarkTools
using Foresight
using ComplexityMeasures
using Distributions
using LinearAlgebra

@testset "progressmap" begin
    S = [1:1000 for _ in 1:3]
    L = @test_nowarn progressmap(S; description = "outer") do s
        progressmap(x -> (sleep(0.0001); randn()), s; description = "inner")
    end
    L = @test_nowarn progressmap(S; description = "outer", parallel = true) do s
        progressmap(x -> (sleep(0.0001); randn()), s; description = "inner")
    end
    L = @test_nowarn progressmap(S; description = "outer", parallel = true) do s
        progressmap(x -> (sleep(0.0001); randn()), s; description = "inner",
                    parallel = true)
    end

    s = [x for x in 1:3]
    _L = map(exp10, collect(s))
    L = @test_nowarn progressmap(exp10, collect(s))
    @test _L == L
    @test typeof(L) == typeof(_L)

    S = [stack(X(1:3), [colorednoise(0:0.1:100) for _ in 1:1000]) for _ in 1:3]

    _L = map(exp10, collect(S[1]))
    L = @test_nowarn progressmap(exp10, collect(S[1]))
    @test _L == L
    @test typeof(L) == typeof(_L)

    _L = map(S) do s
        map(exp10, collect(s))
    end
    L = @test_nowarn progressmap(S) do s
        progressmap(exp10, collect(s))
    end
    @test _L == L
    @test L isa Vector{Any}

    _L = map(S) do s
        map(sum, collect(eachslice(s, dims = (2,))))
    end
    L = @test_nowarn progressmap(S) do s
        progressmap(sum, collect(eachslice(s, dims = (2,))))
    end
    @test _L == L
    @test L isa Vector{Any}

    _L = map(S) do s
        map(sum, eachslice(s, dims = (2,)))
    end
    L = @test_nowarn progressmap(S; parallel = false) do s
        progressmap(sum, eachslice(s, dims = (2,)); parallel = false)
    end
    @test _L == L
    @test typeof(L[1]) != typeof(_L[1]) # This is a limitation of the current implementation; does not handle rebuilding slices of a dim array
    L = @test_nowarn progressmap(S; parallel = false) do s
        map(sum, eachslice(s, dims = (2,)))
    end
    @test _L == L
    @test typeof(L[1]) == typeof(_L[1])
end

@testset "Dates" begin
    x = 1:100
    t = DateTime(1901):Year(1):DateTime(2000)
    y = @test_nowarn TimeSeries(t, x)

    @test samplingperiod(y) == Year(1)
    @test times(y) == t
    @test duration(y) == last(t) - first(t)
    @test unit(y) == NoUnits
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

    xu = set(x, Ti => ustripall(t) * u"s")
    Pxu = @test_nowarn powerspectrum(xu, f_min)
    @test unit(eltype(Pxu)) == u"s"
    @test unit(eltype(lookup(Pxu, 1))) == u"s^-1"
    @test all(ustripall(Pxu) .== Pxx)

    @test_throws "DomainError" powerspectrum(x, 1e-6)

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

    for i in axes(Pxx_mts, 2)
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
    # f = Figure() ax = Axis(f[1, 1]) @test_nowarn lines!(ax, TimeseriesTools.freqs(Pa), Pa)
    # @test_nowarn lines!(ax, TimeseriesTools.freqs(Pb), Pb) save("tmp.pdf", f)
end

include("Types.jl")
include("Utils.jl")
include("IO.jl")
include("Unitful.jl")
include("SpikeTrains.jl")
include("MakieExt.jl")
include("TimeseriesSurrogatesExt.jl")
include("Extensions.jl")
