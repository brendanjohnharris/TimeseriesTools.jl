@testitem "ToolsArrays" begin
    import DimensionalData: span
    x = ToolsArray(randn(10), (𝑡(1:10),))
    @test x isa ToolsArray
    @test !(x isa DimensionalData.DimArray)
    @test x isa DimensionalData.AbstractDimArray

    using DimensionalData
    import DimensionalData: ForwardOrdered, Regular, Points, Sampled, Metadata, order,
                            sampling, layerdims, index, locus, Intervals, intervalbounds
    a = [1 2; 3 4]
    a2 = [1 2 3 4
          3 4 5 6
          4 5 6 7]
    xmeta = Metadata(:meta => "X")
    ymeta = Metadata(:meta => "Y")
    tmeta = Metadata(:meta => "T")
    ameta = Metadata(:meta => "da")
    dimz = (X(Sampled(143.0:2.0:145.0; order = ForwardOrdered(), metadata = xmeta)),
            Y(Sampled(-38.0:2.0:-36.0; order = ForwardOrdered(), metadata = ymeta)))
    dimz2 = (Dim{:row}(10:10:30), Dim{:column}(-20:10:10))

    refdimz = (𝑡(1:1; metadata = tmeta),)
    da = @test_nowarn ToolsArray(a, dimz; refdims = refdimz, name = :test, metadata = ameta)
    val(dims(da, 1)) |> typeof
    da2 = ToolsArray(a2, dimz2; refdims = refdimz, name = :test2)
    lx = Sampled(143.0:2.0:145.0, ForwardOrdered(), Regular(2.0), Points(), xmeta)
    ly = Sampled(-38.0:2.0:-36.0, ForwardOrdered(), Regular(2.0), Points(), ymeta)
    db = DimArray(da)
    @test db isa DimArray
    @test dims(da) == dims(db)
    @test dims(db, X) == dims(da, X)
    @test refdims(db) == refdims(da)
    @test name(db) == name(da)
    @test metadata(db) == metadata(da)
    @test lookup(db) == lookup(da)
    @test order(db) == order(da)
    @test sampling(db) == sampling(da)
    @test span(db) == span(da)
    @test locus(db) == locus(da)
    @test bounds(db) == bounds(da)
    @test layerdims(db) == layerdims(da)
    @test index(db, Y) == index(da, Y)
    da_intervals = set(da, X => Intervals, Y => Intervals)
    db_intervals = set(db, X => Intervals, Y => Intervals)
    @test intervalbounds(da_intervals) == intervalbounds(db_intervals)
end

@testitem "Printing" begin
    # * Vector
    x = ToolsArray(randn(10), (𝑡(1:10),))
    @test_nowarn show(x)
    @test_nowarn display(x)

    # * Matrix
    x = ToolsArray(randn(10, 10), (𝑡(1:10), 𝑥(1:10)))
    @test_nowarn show(x)
    @test_nowarn display(x)

    # * Array
    x = ToolsArray(randn(2, 2, 2), (𝑡(1:2), 𝑥(1:2), 𝑦(1:2)))
    @test_nowarn show(x)
    @test_nowarn display(x)

    # * 2D Array of arrays.
    x = ToolsArray(fill(randn(3, 2), 10, 10), (𝑡(1:10), 𝑥(1:10)))
    @test_nowarn show(x)
    @test_nowarn display(x)

    # * 3D Array of arrays.
    x = ToolsArray(fill(randn(2, 2), 10, 10, 10), (𝑡(1:10), 𝑥(1:10), 𝑦(1:10)))
    @test_nowarn show(x)
    @test_nowarn display(x)

    # * 4D Array of arrays.
    x = ToolsArray(fill(randn(2, 2), 3, 3), (𝑡(1:3), 𝑥(1:3)))
    x = ToolsArray(fill(x, 2, 1, 4, 1),
                   (𝑡(1:2), 𝑥(1:1), 𝑦(1:4), 𝑧(1:1)))
    @test_nowarn show(x)
    @test_nowarn display(x)

    # * Spike train
    x = ToolsArray(rand(Bool, 5), (𝑡(1:5),))
    @test_nowarn show(x)
    @test_nowarn display(x)
end

@testitem "TimeseriesTools.jl" begin
    ts = 1:100
    x = @test_nowarn Timeseries(ts, randn(100))
    @test x isa AbstractTimeseries
    @test x isa RegularTimeseries
    @test x isa UnivariateTimeseries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1 / step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test Interval(x) == first(extrema(ts)) .. last(extrema(ts))
    @test x[𝑡(1 .. 10)] == x[1:10]
    @test all(x[𝑡(At(1:10))] .== x[1:10])
    # @test x[ 𝑡(At(1:10))] != x[1:10]
end

@testitem "Dim queries" begin
    ts = 1:100
    cs = 1:10
    x = Timeseries(ts, TDim{:channel}(cs), randn(100, 10))
    @test x isa ToolsArray
    @test all(isa.(dims(x), ToolsDimension))
    @test dims(x) == (𝑡(ts), TDim{:channel}(cs))
    @test dims(x, 1) == 𝑡(ts)
    @test dims(x, 𝑡) == 𝑡(ts)
    @test dims(x, TDim{:channel}) == TDim{:channel}(cs)
    @test !(dims(x, :channel) == TDim{:channel}(cs)) # You CAN'T use symbols for TDim{}s because DimensionalData.name2dim always returns a Dim{}

    DimensionalData.@dim U ToolsDim "U"
    @test U <: ToolsDimension

    x = ToolsArray(randn(10), (𝑥(1:10),))
    @test all(x[At(dims(x, 1))] .== x)
    @test lookup(x[At(dims(x, 1))]) != lookup(x) # One is Regular, one is Irregular
    @test all(lookup(x[At(dims(x, 1))]) .== lookup(x[At(1:10)])) # But same elements
end

@testitem "Multivariate time series" begin
    ts = 1:100
    x = @test_nowarn Timeseries(ts, 1:5, randn(100, 5))
    @test x isa AbstractTimeseries
    @test x isa RegularTimeseries
    @test x isa MultivariateTimeseries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1 / step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test Interval(x) == first(extrema(ts)) .. last(extrema(ts))
    @test x[𝑡(1 .. 10)] == x[1:10, :]
end

@testitem "Multidimensional time series" begin
    x = @test_nowarn Timeseries(𝑡(1:100), X(1:10), randn(100, 10))
    @test x isa AbstractTimeseries
    @test x isa RegularTimeseries
    @test x isa MultidimensionalTimeseries

    x = @test_nowarn Timeseries(𝑡(1:100), X(1:10), Y(1:10), randn(100, 10, 10))
    @test x isa AbstractTimeseries
    @test x isa RegularTimeseries
    @test x isa MultidimensionalTimeseries
    @test_nowarn x[𝑡(Near(4:10))]

    x = @test_nowarn Timeseries(𝑡(1:100), X(randn(10) |> sort), Y(1:10),
                                randn(100, 10, 10))
    @test x isa AbstractTimeseries
    @test x isa RegularTimeseries
    @test !(x isa MultidimensionalTimeseries)

    x = @test_nowarn Timeseries(𝑡(sort(randn(100))), randn(100))
    @test x isa AbstractTimeseries
    @test !(x isa RegularTimeseries)
    @test !(x isa MultidimensionalTimeseries)
end
