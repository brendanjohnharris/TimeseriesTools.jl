@testset "Interlace" begin
    x = TimeSeries(0:0.1:1, randn(11))
    y = TimeSeries(0.05:0.1:1, randn(10))
    z = @test_nowarn interlace(x, y)
    @test all(collect(times(z)) .== 0.0:0.05:1.0)
end

@testset "Cat" begin
    x = TimeSeries(0.1:0.1:10, Var(1:100), randn(100, 100))
    y = cat(Freq(1:2), x, x)
    @test dims(y, 3) == Freq(1:2)
    z = stack(Freq(1:2), [x, x])
    @test y == z
    y = stack(Freq(1:2), [x, x]; dims = 1)
    @test dims(y, 1) == Freq(1:2)
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

@testset "Rectification" begin
    ts = 0.1:0.1:1000
    x = TimeSeries(ts .+ randn(length(ts)) .* 1e-9, sin)
    @test issorted(times(x))
    _x = @test_nowarn rectifytime(x)
    @test all(x .== _x)
    @test ts == times(_x)

    y = TimeSeries(ts .+ randn(length(ts)) .* 1e-9, cos)
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
    dϕ = @test_nowarn centraldiff(ϕ; grad = phasegrad)
    @test sum(dϕ .!= centraldiff(ϕ)) > 4000
end

@testset "Left and right derivatives" begin
    x = colorednoise(0.01:0.01:10)
    X = cat(Var(1:10), [colorednoise(0.1:0.1:100) for _ in 1:10]...)

    dx = @test_nowarn leftdiff(x)
    @test all(dx[2:(end)] .== (x[2:end] - x[1:(end - 1)]))
    @test times(dx) == times(x)

    dX = @test_nowarn leftdiff(X)
    @test all(dX[2:(end), :] .== (X[2:end, :] - X[1:(end - 1), :]))
    @test times(dX) == times(X)
    @test dims(dX, Var) == dims(X, Var)

    dx = @test_nowarn rightdiff(x)
    @test all(dx[1:(end - 1)] .== (x[2:end] - x[1:(end - 1)]))
    @test times(dx) == times(x)

    dX = @test_nowarn rightdiff(X)
    @test all(dX[1:(end - 1), :] .== (X[2:end, :] - X[1:(end - 1), :]))
    @test times(dX) == times(X)
    @test dims(dX, Var) == dims(X, Var)
end

# @testset "Irregular central derivative" begin
#     ts = 0.1:0.1:1000
#     x = TimeSeries(ts, sin)
#     y = TimeSeries(ts .+ randn(length(ts)) .* 1e-10, parent(x))
#     @test centralderiv(x) ≈ centralderiv(y)
# end

@testset "Unitful derivative" begin
    ts = 0.1:0.1:1000
    x = TimeSeries(ts, sin)
    y = set(x, Ti => ts .* u"s")
    @test ustripall(centralderiv(x)) == ustripall(centralderiv(y))
    @test unit(eltype(centralderiv(y))) == unit(u"1/s")
end

@testset "coarsegrain" begin
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

    X = TimeSeries(1:11, 1:100, repeat(1:11, 1, 100))
    C = coarsegrain(X, dims = 1)
    M = dropdims(mean(C, dims = 3), dims = 3)
    @test all(M[:, 1] .== 1.5:2:9.5)
    @test size(C, 1) == size(X, 1) ÷ 2

    C = coarsegrain(X)
    @test size(C) == (5, 50, 4)
    M = dropdims(mean(C, dims = 3), dims = 3)
    @test all(M[:, 1] .== 1.5:2:9.5)

    C = coarsegrain(X; dims = Ti, newdim = Var)
    @test length(lookup(C, 1)) == size(C, 1)
    @test length(lookup(C, 2)) == size(C, 2)
    M = mean(C.data, dims = 2)
    @test size(C) == (5, 200)
    @test all(M[:, 1] .== 1.5:2:9.5)

    X = cat(X, X; dims = 3)
    C = coarsegrain(X; dims = 1, newdim = 2)
    @test size(C) == (5, 200, 2)
    @test_nowarn C[Ti(Near(0.1))]
end

@testset "matchdim" begin
    ts = 0:1:100
    X = [Timeseries(ts .+ 1e-6 .* randn(101), sin) for _ in 1:10]
    X = TimeSeries(1:10, X)
    Y = matchdim(X)

    @test length(unique(dims.(Y))) == 1
    @test dims(Y[1], Ti) == Ti(ts)
end

@testset "findpeaks" begin
    x = TimeSeries(0.1:0.1:100, x -> sin(x .* 2π / 4))
    peaks = spiketrain(range(start = 1, stop = 100, step = 4))
    pks, proms = findpeaks(x)
    @test times(pks) == times(peaks)

    X = cat(Var(1:2), x, x .+ 1.0)
    pks, proms = findpeaks(X)
    @test pks isa DimArray{<:DimArray}
    @test proms isa DimArray{<:DimArray}
    @test times(pks[1]) == times(peaks)
    @test times(pks[2]) == times(peaks)
    @test all(proms[1][2:(end - 1)] .== 2)

    m = @test_nowarn maskpeaks(x)
    M = @test_nowarn maskpeaks(X)
    @test M[:, 2] == m

    # * With no-peak signal
    X[:, 2] .= 1
    @test_nowarn findpeaks(X)
    M = @test_nowarn maskpeaks(X)
end
