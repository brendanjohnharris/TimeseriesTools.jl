@testitem "Dagger progressmap" begin
    using Distributed
    addprocs(2)
    @everywhere using DimensionalData
    @everywhere using Dagger
    @everywhere using TimeseriesTools
    TimeseriesTools.PROGRESSMAP_BACKEND = :Dagger

    S = 1:1000
    @everywhere f(x) = (sleep(0.001); DimArray(randn(10), (X(1:10),)))
    out = progressmap(f, S; numblocks = 100)

    # * Check for matrix
    S = randn(10, 10)
    @everywhere h(x) = DimArray((randn(1000, 1000)^5), (X(1:1000), Y(1:1000)))
    out = progressmap(h, S; numblocks = 100)
    @test out isa Matrix{<:DimMatrix}
    @test size(out) == size(S)

    # * Check for DimArray
    S = DimArray(randn(10, 10), (X(1:10), Y(1:10)))
    out = progressmap(h, S; numblocks = 100)
    @test out isa DimMatrix{<:DimMatrix}
    @test size(out) == size(S)

    S = ToolsArray(randn(10, 10), (TDim{:X}(1:10), TDim{:Y}(1:10))) # Can't do this yet. Don't use TDim.
    @test_throws "DimensionMismatch" progressmap(h, S; numblocks = 100)
    rmprocs()
end
