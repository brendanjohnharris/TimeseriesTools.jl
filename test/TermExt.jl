@testitem "Term progressmap" begin
    using Term

    S = [1:100 for _ in 1:3]
    TimeseriesTools.PROGRESSMAP_BACKEND = nothing
    @test_throws "ArgumentError" TimeseriesTools._progressmap(nothing, S) do s
        sleep(0.00001)
    end
    @test_throws "ArgumentError" progressmap(S) do s
        sleep(0.00001)
    end

    TimeseriesTools.PROGRESSMAP_BACKEND = :Term
    L = @test_nowarn progressmap(S; description = "outer") do s
        progressmap(x -> (sleep(0.00001); randn()), s; description = "inner")
    end
    L = @test_nowarn progressmap(S; description = "outer", parallel = true) do s
        progressmap(x -> (sleep(0.00001); randn()), s; description = "inner")
    end
    L = @test_nowarn progressmap(S; description = "outer", parallel = true) do s
        progressmap(x -> (sleep(0.00001); randn()), s; description = "inner",
                    parallel = true)
    end

    s = [x for x in 1:3]
    _L = map(exp10, collect(s))
    L = @test_nowarn progressmap(exp10, collect(s))
    @test _L == L
    @test typeof(L) == typeof(_L)

    S = [stack(X(1:3), [colorednoise(0:0.1:100) for _ in 1:1000]) for _ in 1:3]

    _L = map(exp10, collect(S[1]))
    L = @test_nowarn progressmap(exp10, collect(S[1]))
    @test _L == L
    @test typeof(L) == typeof(_L)

    _L = map(S) do s
        map(exp10, collect(s))
    end
    L = @test_nowarn progressmap(S) do s
        progressmap(exp10, collect(s))
    end
    @test _L == L
    @test L isa Vector{Any}

    _L = map(S) do s
        map(sum, collect(eachslice(s, dims = (2,))))
    end
    L = @test_nowarn progressmap(S) do s
        progressmap(sum, collect(eachslice(s, dims = (2,))))
    end
    @test _L == L
    @test L isa Vector{Any}

    _L = map(S) do s
        map(sum, eachslice(s, dims = (2,)))
    end
    L = @test_nowarn progressmap(S; parallel = false) do s
        progressmap(sum, eachslice(s, dims = (2,)); parallel = false)
    end
    @test _L == L
    @test typeof(L[1]) != typeof(_L[1]) # This is a limitation of the current implementation; does not handle rebuilding slices of a dim array
    L = @test_nowarn progressmap(S; parallel = false) do s
        map(sum, eachslice(s, dims = (2,)))
    end
    @test _L == L
    @test typeof(L[1]) == typeof(_L[1])
end
