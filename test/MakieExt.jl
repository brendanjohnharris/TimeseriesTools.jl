
@testitem "Makie" begin
    using CairoMakie
    import TimeseriesTools: TimeSeries
    Traces = Base.get_extension(TimeseriesTools, :MakieExt).Traces
    x = TimeSeries(0.01:0.01:10, randn(1000))

    p = @test_nowarn plot(x)
    @test p.plot isa Lines
    # @test 10 ≤ p.axis.finallimits.val.widths[1] < 12
    x = TimeSeries(0.01:0.01:10, 1:2, randn(1000, 2))
    p = @test_nowarn plot(x)
    @test p.plot isa Traces
    # @test 10 ≤ p.axis.finallimits.val.widths[1] < 12
    # @test 2 ≤ p.axis.finallimits.val.widths[2] < 3

    x = TimeSeries(0.01:0.01:10, [1, 3], randn(1000, 2))
    p = @test_nowarn plot(x)
    @test p.plot isa Heatmap
end

@testitem "Readme" begin
    using TimeseriesTools, CairoMakie, Unitful, Foresight
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
    x = loadtimeseries("./test_timeseries.tsv")

    f = Figure(; size = (500, 480))
    ax = Axis3(f[1, 1])
    trajectory!(ax, collect.(eachcol(x))...; colormap = :turbo, linewidth = 0.1)
    shadows!(ax, collect.(eachcol(x))...; color = (:slategray, 0.5), linewidth = 0.05)
    ax.xlabelvisible = ax.ylabelvisible = ax.zlabelvisible = ax.xticksvisible = ax.yticksvisible = ax.zticksvisible = ax.xticklabelsvisible = ax.yticklabelsvisible = ax.zticklabelsvisible = false
    ax.azimuth = ax.azimuth[] + 0.1
    ax.elevation = ax.elevation[] + 0.1
    save("./shadows.png", f; px_per_unit = 3)

    @test_nowarn trajectory!(ax, x)
    @test_nowarn shadows!(ax, x; color = (:slategray, 0.5), linewidth = 0.05)
end

@testitem "Readme_dark" begin
    using CairoMakie, TimeseriesTools, Unitful, Foresight
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
    x = loadtimeseries("./test_timeseries.tsv")

    f = Figure(; size = (500, 480))
    ax = Axis3(f[1, 1])
    trajectory!(ax, collect.(eachcol(x))...; colormap = :turbo, linewidth = 0.1)
    shadows!(ax, collect.(eachcol(x))...; color = (:white, 0.5), linewidth = 0.05)
    ax.xlabelvisible = ax.ylabelvisible = ax.zlabelvisible = ax.xticksvisible = ax.yticksvisible = ax.zticksvisible = ax.xticklabelsvisible = ax.yticklabelsvisible = ax.zticklabelsvisible = false
    ax.azimuth = ax.azimuth[] + 0.1
    ax.elevation = ax.elevation[] + 0.1
    save("./shadows_dark.png", f; px_per_unit = 3)
end

@testitem "Traces" begin
    using CairoMakie, TimeseriesTools, Unitful
    import TimeseriesTools.TimeSeries # or TS

    t = 0.005:0.005:1e4
    x = colorednoise(t, u"s") * u"V"
    X = cat(Var(1:2), x, x .+ 1.0 * u"V", dims = 2)

    # Calculate the power spectrum
    S = _powerspectrum(x, 0.0005)[2:end, :]
    f = Figure(; size = (720, 480))
    ax = Axis(f[1, 1], xscale = log10, yscale = log10)
    # x, y, z = collect.(ustripall(decompose(S)))
    @test_nowarn traces!(ax, S; colormap = :turbo)
end

@testitem "Spectrum plot" begin
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
    spectrumplot(S; peaks = true) # * Log-log plot

    # * Test recipes
    S = ustripall(S)
    f, ax, p = plot(S) # * Log-log plot
    tightlimits!(ax)
    @test ax isa Axis
    @test p isa SpectrumPlot
    @test ax.finallimits.val.widths[1]≈100 atol=1e-2 # Uses freq values
end
