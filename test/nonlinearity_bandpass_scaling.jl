using CairoMakie
using FFTW
using Statistics
using TimeseriesSurrogates
using ComplexityMeasures
using Term
using TimeseriesTools

file = "test/test_timeseries.tsv"
X = loadtimeseries(file)

begin
    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, x[5:end], x[1:(end - 4)], linewidth = 0.1)
    f
end

begin
    x = X[Ti = 500 .. 600, Var = At(:x)]
    x = rectify(x; dims = Ti)
    # method = WLS(IAAFT(), rescale = true)
    method = IAAFT()
    function F(x)
        y = @views x[100:(end - 100)] # Remove edge effects
        est = MissingDispersionPatterns(Dispersion(m = 3, c = 7))
        complexity_normalized(est, y)
        # entropy_sample(y)
    end
    T = SurrogateTest(F, parent(x), method; n = 100)
    p = pvalue(T, tail = :right)
end
begin
    x = X[Ti = 500 .. 10000, Var = At(:x)]
    x = rectify(x, dims = Ti)
    # fs = range(log10(0.1), log10(9), length = 10) .|> exp10
    fs = 1:9
    # fs = 1:2:20
    ps = progressmap(fs; parallel = true) do f
        y = bandpass(x, (10 - f / 2) .. (10 + f / 2))
        # y = x[1:f:end]
        T = SurrogateTest(F, parent(y), method; n = 500)
        pvalue(T, tail = :right)
        # return (F(parent(y)) - mean(T.vals)) ./ std(T.vals)
        return F(parent(y)), mean(T.vals), std(T.vals)
    end
    ts = first.(ps)
    ms = getindex.(ps, 2)

    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, fs, ts)
    lines!(ax, fs, ms)
    f
end

begin # * Test with varying surrogate threshold (partial randomization)
    ts = 0.01:0.01:10
    αs = 0:0.1:3

    method = IAAFT()
    function F(x)
        y = @views x[100:(end - 100)] # Remove edge effects
        est = MissingDispersionPatterns(Dispersion(m = 3, c = 7))
        complexity_normalized(est, y)
        # entropy_sample(y)
    end

    ps = progressmap(αs; parallel = true) do α
        y = colorednoise(ts; α)
        T = SurrogateTest(F, parent(y), method; n = 1000)
        pvalue(T, tail = :right)
        return F(parent(y)), mean(T.vals), std(T.vals)
    end

    ss = first.(ps)
    ms = getindex.(ps, 2)

    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, αs, ss)
    lines!(ax, αs, ms)
    f
end

begin # Scaling of complexity with 1/f exponent
    ts = 0.01:0.01:10
    αs = 0:0.1:3

    method = IAAFT()
    function F(x)
        y = @views x[100:(end - 100)] # Remove edge effects
        est = MissingDispersionPatterns(Dispersion(m = 3, c = 7))
        complexity_normalized(est, y)
        # entropy_sample(y)
    end

    ps = progressmap(αs; parallel = true) do α
        y = colorednoise(ts; α)
        T = SurrogateTest(F, parent(y), method; n = 1000)
        pvalue(T, tail = :right)
        return F(parent(y)), mean(T.vals), std(T.vals)
    end

    ss = first.(ps)
    ms = getindex.(ps, 2)

    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, αs, ss)
    lines!(ax, αs, ms)
    f
end
