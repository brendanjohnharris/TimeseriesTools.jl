@testitem "IO" begin
    import TimeseriesTools: Timeseries
    using Unitful
    using JLD2
    x = Timeseries(rand(1000, 3), 0.001:0.001:1, 1:3; metadata = Dict(:a => :test),
                   name = "name")

    f = tempname() * ".jld2"
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x == _x

    f = tempname() * ".tsv"
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test all(x .≈ _x)

    x = x[:, 1]
    savetimeseries(f, x)
    _x = @test_logs (:warn, "Cannot load refdims yet") loadtimeseries(f)
    @test refdims(_x) == ()
    @test all(x .≈ _x)

    x = Timeseries(rand(1000, 3), 0.001:0.001:1, 1:3; metadata = Dict(:a => :test))
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x ≈ _x

    x = Timeseries(rand(1000, 3), 0.001:0.001:1, X(1:3); metadata = Dict(:a => :test))
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x ≈ _x
    @test [all(d .≈ _d) for (d, _d) in zip(dims(x), dims(_x))] |> all
    @test parent(lookup(_x, 1)) isa Vector{Float64}

    # Currently not the greatest way of handling non-serializable metadata
    x = Timeseries(rand(1000, 3), 0.001:0.001:1, 1:3;
                   metadata = Dict(:a => DimensionalData.NoName())) # Something that can't be serialized
    @test_logs (:warn, ErrorException("Cannot serialize type DimensionalData.NoName")) savetimeseries(f,
                                                                                                      x)
    _x = loadtimeseries(f)
    @test metadata(_x) == DimensionalData.Dimensions.LookupArrays.NoMetadata()
    @test all(x .≈ _x)
    @test [all(d .≈ _d) for (d, _d) in zip(dims(x), dims(_x))] |> all
    @test parent(lookup(_x, 1)) isa Vector{Float64}

    x = Timeseries(rand(1000, 3), 0.001:0.001:1, 1:3; name = Timeseries) # Something that can't be serialized
    @test_logs (:warn,
                ErrorException("Cannot serialize type typeof(TimeseriesTools.Timeseries)")) savetimeseries(f,
                                                                                                           x)
    _x = loadtimeseries(f)
    @test name(_x) == DimensionalData.NoName()
    @test all(x .≈ _x)
    @test [all(d .≈ _d) for (d, _d) in zip(dims(x), dims(_x))] |> all
    @test parent(lookup(_x, 1)) isa Vector{Float64}

    x = Timeseries(rand(1000, 3), 0.001:0.001:1, [Timeseries, Timeseries, Timeseries])
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test all(parent(x) .≈ parent(_x))
    @test parent(lookup(_x, 1)) isa Vector{Float64}
    @test parent(lookup(_x, 2)) isa Vector{Symbol}

    x = Timeseries(rand(1000), 0.001:0.001:1)
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test all(x .≈ _x)
    @test [all(d .≈ _d) for (d, _d) in zip(dims(x), dims(_x))] |> all
    @test parent(lookup(_x, 1)) isa Vector{Float64}

    x = Timeseries(rand(1000, 3), (0.001:0.001:1) * u"s", 1:3; metadata = Dict(:a => :test),
                   name = "name") * u"V"

    f = tempname() * ".jld2"
    savetimeseries(f, x)
    _x = loadtimeseries(f)
    @test x == _x
end
