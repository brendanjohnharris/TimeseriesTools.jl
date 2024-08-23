using TimeseriesTools.Operators
using Test

@testset "Operators" begin
    x = colorednoise(1:1000)
    _x = deepcopy(x)
    @test_nowarn ℬ!(x)
    @test all(x[2:end] .== _x[1:(end - 1)])
    ℬ!(x, 3)
    @test all(x[5:end] .== _x[1:(end - 4)])
    @test ℬ(_x, 4) == x

    x = deepcopy(_x)
    @test_nowarn ℒ!(x)
    @test all(_x[2:end] .== x[1:(end - 1)])
    ℒ!(x, 3)
    @test all(_x[5:end] .== x[1:(end - 4)])
    @test ℒ(_x, 4) == x

    x = colorednoise(0.0:0.01:1)
    T = 𝒯(-1)
    @test times(T(x)) == -1:step(x):0
end
