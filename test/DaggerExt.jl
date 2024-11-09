@testitem "Dagger progressmap" begin
    using Dagger
    using DimensionalData
    TimeseriesTools.PROGRESSMAP_BACKEND = :Dagger

    S = 1:1000
    g = x -> (sleep(0.001); randn())
    out = progressmap(g, S; numblocks = 100)

    # * Check for matrix
    S = randn(10, 10)
    h(x) = (randn(1000, 1000)^10)^(-10)
    out = progressmap(h, S; numblocks = 100)
    @test out isa Matrix{Matrix{Float64}}

    # * Check for DimArray
    S = DimArray(randn(10, 10), (X(1:10), Y(1:10)))
    out = progressmap(h, S; numblocks = 100)
    @test out isa DimMatrix{Matrix{Float64}}
end
