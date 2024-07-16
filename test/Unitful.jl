
@testset "Unitful" begin
    ts = (1:1000)u"s"
    x = @test_nowarn TimeSeries(ts, randn(1000))
    @test TimeSeries(ustripall(ts), collect(x), u"s") == x
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

    @test_nowarn rectify(x; dims = Ti)
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
    @test all(round.(ustripall.(peakfs)) .∈ ([50, 100],))
    @test first(peakamps) / last(peakamps)≈4.2 / 3.1 rtol=1e-1

    x = exp.(-ustripall(ts) .^ 2)
    x = TimeSeries(ts, x * u"V")
    ℱ = sqrt(π) .* exp.(-π^2 .* ustripall(f) .^ 2)
    _S = abs.(ℱ) .^ 2 * u"V^2*s^2"
    S = energyspectrum(x, 0.0)
    @test sum(_S) .* step(f)≈sum(S) .* step(dims(S, Freq)) rtol=0.05

    lines(ustripall(f), ustripall(_S), axis = (; limits = ((0, 1), (0, 4))))
    plot!(collect(ustripall(dims(S, Freq))), collect(ustripall(S)))
    current_figure()
end
