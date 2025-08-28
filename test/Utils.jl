@testitem "Interlace" begin
    x = Timeseries(randn(11), 0:0.1:1)
    y = Timeseries(randn(10), 0.05:0.1:1)
    z = @test_nowarn interlace(x, y)
    @test all(collect(times(z)) .== 0.0:0.05:1.0)
end

@testitem "Stack" begin
    x = Timeseries(randn(100, 100), 0.1:0.1:10, Var(1:100))
    z = ToolsArray([x, x], 𝑓(1:2)) |> stack
    y = stack(DimArray([x, x], 𝑓(1:2)); dims = 1)
    @test dims(y, 1) == 𝑓(1:2)
end

@testitem "Buffer" begin
    N = 10
    x = Timeseries(randn(100), 0.1:0.1:10)
    y = @test_nowarn buffer(x, 10)
    @test length(y) == N
    @test y[1] == x[1:(length(x) ÷ N)]
    @test cat(y..., dims = 𝑡) == x[1:((length(x) ÷ N) * N)]

    y = @test_nowarn buffer(x, 10, 0; discard = false)
    @test cat(y..., dims = 𝑡) == x

    y = @test_nowarn buffer(x, 10, N ÷ 2)
    @test length(y) == 2 * N - 1

    x = Timeseries(randn(101, 10), 0:0.1:10, 1:10)
    y = buffer(x, 10)
    @test length(y) == 10

    x = Timeseries(randn(100), 0.1:0.1:10)
    y = @test_nowarn window(x, 2, 1)
    @test all(length.(y) .== 2)
    y = @test_nowarn delayembed(x, 2, 1, 1)
    y = @test_nowarn delayembed(x, 2, 1, 2)
    @test length(y) == length(x)
    @test samplingperiod(y) == 2 * samplingperiod(x)
    y = @test_nowarn delayembed(x, 2, 2, 1)
end

@testitem "Rectification" begin
    import TimeseriesTools: rectifytime
    ts = 0.1:0.1:1000
    x = Timeseries(sin, ts .+ randn(length(ts)) .* 1e-10)
    @test issorted(times(x))
    _x = @test_nowarn rectifytime(x)
    @test all(x .== _x)
    @test ts == times(_x)

    y = Timeseries(cos, ts .+ randn(length(ts)) .* 1e-10)
    @test issorted(times(y))
    _x, _y = (rectifytime(x, y)...,)

    @test all(x .== parent(_x))
    @test ts == times(_x)
    @test all(y .== parent(_y))
    @test ts == times(_y)

    x = @test_nowarn Timeseries(randn(100, 10), 𝑡(1:100), X((1:10) .+ 1e-10 .* randn(10)))
    y = @test_nowarn rectify(x, dims = X)
    @test dims(y, X) == X(1:10)

    x = @test_nowarn Timeseries(randn(100, 10, 5), 𝑡(1:100),
                                X((1:10) .+ 1e-10 .* randn(10)),
                                Y((1:5) .+ 1e-10 .* randn(5)))
    y1 = @test_nowarn rectify(x, dims = X)
    y2 = @test_nowarn rectify(x, dims = Y)
    y3 = @test_nowarn rectify(x, dims = [X, Y])
    @test dims(y1, X) == dims(y3, X) == X(1:10)
    @test dims(y2, Y) == dims(y3, Y) == Y(1:5)
    @test dims(y1, Y) == dims(x, Y)
    @test dims(y2, X) == dims(x, X)
end

@testitem "Central differences" begin
    using DSP, StatsBase
    x = colorednoise(0.01:0.01:10)
    X = ToolsArray([colorednoise(0.1:0.1:100) for _ in 1:10], Var(1:10)) |> stack

    dx = @test_nowarn centraldiff(x)
    @test all(dx[2:(end - 1)] .== (parent(x)[3:end] - parent(x)[1:(end - 2)]) / 2)
    @test times(dx) == times(x)

    dX = @test_nowarn centraldiff(X)
    @test all(dX[2:(end - 1), :] .== (parent(X)[3:end, :] - parent(X)[1:(end - 2), :]) / 2)
    @test times(dX) == times(X)
    @test dims(dX, Var) == dims(X, Var)

    dX = @test_nowarn centralderiv(X)
    @test all(dX[2:(end - 1), :] .==
              ((parent(X)[3:end, :] - parent(X)[1:(end - 2), :]) / 2) ./ samplingperiod(X))

    x = @test_nowarn Timeseries(sin, 0.1:0.1:1000)
    𝑓 = instantaneousfreq(x)
    @test std(𝑓[2500:(end - 2500)]) < 0.001
    @test mean(𝑓[2500:(end - 2500)])≈1 / 2π rtol=1e-5

    ϕ = analyticphase(x)[1000:(end - 1000)]
    dϕ = @test_nowarn centraldiff(ϕ; grad = phasegrad)
    @test sum(dϕ .!= centraldiff(ϕ)) > 4000
end

@testitem "Left and right derivatives" begin
    import TimeseriesTools: leftdiff, rightdiff
    x = colorednoise(0.01:0.01:10)
    X = ToolsArray([colorednoise(0.1:0.1:100) for _ in 1:10], Var(1:10)) |> stack

    dx = @test_nowarn leftdiff(x)
    @test all(parent(dx)[2:(end)] .== (parent(x)[2:end] - parent(x)[1:(end - 1)]))
    @test times(dx) == times(x)

    dX = @test_nowarn leftdiff(X)
    @test all(parent(dX)[2:(end), :] .== (parent(X)[2:end, :] - parent(X)[1:(end - 1), :]))
    @test times(dX) == times(X)
    @test dims(dX, Var) == dims(X, Var)

    dx = @test_nowarn rightdiff(x)
    @test all(parent(dx)[1:(end - 1)] .== (parent(x)[2:end] - parent(x)[1:(end - 1)]))
    @test times(dx) == times(x)

    dX = @test_nowarn rightdiff(X)
    @test all(parent(dX)[1:(end - 1), :] .==
              (parent(X)[2:end, :] - parent(X)[1:(end - 1), :]))
    @test times(dX) == times(X)
    @test dims(dX, Var) == dims(X, Var)
end

# @testitem "Irregular central derivative" begin
#     ts = 0.1:0.1:1000
#     x = Timeseries(sin, ts)
#     y = Timeseries(parent(x), ts .+ randn(length(ts)) .* 1e-10)
#     @test centralderiv(x) ≈ centralderiv(y)
# end

@testitem "Unitful derivative" begin
    using Unitful
    ts = 0.1:0.1:1000
    x = Timeseries(sin, ts)
    y = set(x, 𝑡 => ts .* u"s")
    @test ustripall(centralderiv(x)) == ustripall(centralderiv(y))
    @test unit(eltype(centralderiv(y))) == unit(u"1/s")
end

@testitem "coarsegrain" begin
    using StatsBase
    X = repeat(1:11, 1, 100)
    C = coarsegrain(X, dims = 1)
    M = mean(C, dims = 3)
    @test all(M[:, 1] .== 1.5:2:9.5)
    @test size(C, 1) == size(X, 1) ÷ 2

    C = coarsegrain(X)
    @test size(C) == (5, 50, 4)
    M = mean(C, dims = 3)
    @test all(M[:, 1] .== 1.5:2:9.5)

    C = coarsegrain(X; newdim = 2)
    M = mean(C, dims = 2)
    @test size(C) == (5, 200)
    @test all(M[:, 1] .== 1.5:2:9.5)

    X = cat(X, X; dims = 3)
    C = coarsegrain(X; dims = 1, newdim = 2)
    @test size(C) == (5, 200, 2)

    X = Timeseries(repeat(1:11, 1, 100), 1:11, 1:100)
    C = coarsegrain(X, dims = 1)
    M = dropdims(mean(C, dims = 3), dims = 3)
    @test all(M[:, 1] .== 1.5:2:9.5)
    @test size(C, 1) == size(X, 1) ÷ 2

    C = coarsegrain(X)
    @test size(C) == (5, 50, 4)
    M = dropdims(mean(C, dims = 3), dims = 3)
    @test all(M[:, 1] .== 1.5:2:9.5)

    C = coarsegrain(X; dims = 𝑡, newdim = Var)
    @test length(lookup(C, 1)) == size(C, 1)
    @test length(lookup(C, 2)) == size(C, 2)
    M = mean(C.data, dims = 2)
    @test size(C) == (5, 200)
    @test all(M[:, 1] .== 1.5:2:9.5)

    X = cat(X, X; dims = 3)
    C = coarsegrain(X; dims = 1, newdim = 2)
    @test size(C) == (5, 200, 2)
    @test_nowarn C[𝑡(Near(0.1))]
end

@testitem "matchdim" begin
    ts = 0:1:100
    X = [Timeseries(sin, ts .+ 1e-6 .* randn(101)) for _ in 1:10]
    X = Timeseries(X, 1:10)
    Y = matchdim(X)

    @test length(unique(dims.(Y))) == 1
    @test dims(Y[1], 𝑡) == 𝑡(ts)
end

@testitem "findpeaks" begin
    x = TimeseriesTools.Timeseries(x -> sin(x .* 2π / 4), 0.1:0.1:100)
    peaks = spiketrain(range(start = 1, stop = 100, step = 4))
    pks, proms = findpeaks(x)
    @test times(pks) == times(peaks)
    @test pks isa ToolsArray{<:Float64}
    @test proms isa ToolsArray{<:Float64}

    xx = cat(x, x .+ 1.0, dims = Var(1:2))

    # * test rebuild
    identity.(eachslice(xx; dims = 1)) isa ToolsArray

    pks, proms = findpeaks(xx)
    @test pks isa ToolsArray{<:ToolsArray}
    @test proms isa ToolsArray{<:ToolsArray}
    @test times(pks[1]) == times(peaks)
    @test times(pks[2]) == times(peaks)
    @test all(proms[1][2:(end - 1)] .== 2)

    m = @test_nowarn maskpeaks(x)
    M = @test_nowarn maskpeaks(xx)
    @test M[:, 2] == m

    # * With no-peak signal
    xx[:, 2] .= 1
    @test_nowarn findpeaks(xx)
    M = @test_nowarn maskpeaks(xx)
end

@testitem "ProgressLogging progressmap" begin
    using DimensionalData
    TimeseriesTools.PROGRESSMAP_BACKEND = :ProgressLogging
    S = 1:1000
    g = x -> (sleep(0.001); randn())
    out = progressmap(g, S)

    # * Check for matrix
    S = randn(10, 10)
    h(x) = (randn(1000, 1000)^10)^(-10)
    out = progressmap(h, S)
    @test out isa Matrix{Matrix{Float64}}

    # * Check for DimArray
    S = DimArray(randn(10, 10), (X(1:10), Y(1:10)))
    out = progressmap(h, S)
    @test out isa DimMatrix{Matrix{Float64}}
end
@testitem "progressmap schedulers" begin
    TimeseriesTools.PROGRESSMAP_BACKEND = :Threads

    # * Compare to Threads.@threads. Guess you should have some threads.
    function f(X)
        return sum(X * X)
    end
    Ns = 1:100:1000
    Xs = [rand(n, n) for n in Ns]
    @test progressmap(f, Xs, schedule = :dynamic) == progressmap(f, Xs, schedule = :static)

    Ns = range(0, 1, length = 100)

    if Threads.nthreads() > 2 && VERSION ≥ v"1.11"
        _, bs = @timed progressmap(sleep, Ns, schedule = :static)
        _, bd = @timed progressmap(sleep, Ns, schedule = :dynamic)
        # _, bg = @timed progressmap(sleep, Ns, schedule = :greedy)

        @test bs > bd # > bg
    end
end

# begin
#     D = DimensionalData.Dim{:x}(1:100)
#     T = TimeseriesTools.Dim{:x}(1:100)
#     @test D != T
#     @test all(D .== T)
#     @test DimensionalData.name(T) == DimensionalData.name(D)
#     𝑡(1:10)
#     𝑥(1:10)
#     𝑦(1:10)
#     𝑧(1:10)
# end
