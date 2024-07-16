@testset "TimeseriesTools.jl" begin
    ts = 1:100
    x = @test_nowarn TimeSeries(ts, randn(100))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa UnivariateTimeSeries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1 / step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test Interval(x) == first(extrema(ts)) .. last(extrema(ts))
    @test x[Ti(1 .. 10)] == x[1:10]
    @test all(x[Ti(At(1:10))] .== x[1:10])
    # @test x[Ti(At(1:10))] != x[1:10]
end

@testset "Multivariate time series" begin
    ts = 1:100
    x = @test_nowarn TimeSeries(ts, 1:5, randn(100, 5))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa MultivariateTimeSeries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1 / step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test Interval(x) == first(extrema(ts)) .. last(extrema(ts))
    @test x[Ti(1 .. 10)] == x[1:10, :]
end

@testset "Multidimensional time series" begin
    x = @test_nowarn TimeSeries(Ti(1:100), X(1:10), randn(100, 10))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa MultidimensionalTimeSeries

    x = @test_nowarn TimeSeries(Ti(1:100), X(1:10), Y(1:10), randn(100, 10, 10))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test x isa MultidimensionalTimeSeries
    @test_nowarn x[Ti(Near(4:10))]

    x = @test_nowarn TimeSeries(Ti(1:100), X(randn(10) |> sort), Y(1:10),
                                randn(100, 10, 10))
    @test x isa AbstractTimeSeries
    @test x isa RegularTimeSeries
    @test !(x isa MultidimensionalTimeSeries)

    x = @test_nowarn TimeSeries(Ti(sort(randn(100))), randn(100))
    @test x isa AbstractTimeSeries
    @test !(x isa RegularTimeSeries)
    @test !(x isa MultidimensionalTimeSeries)
end
