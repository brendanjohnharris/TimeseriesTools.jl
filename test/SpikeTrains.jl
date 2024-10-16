using Distances
using TimeseriesFeatures
using Test
using TimeseriesTools
using LinearAlgebra

@testset "STOIC feature" begin
    x = gammarenewal(100, 1.0, 1.0)
    y = gammarenewal(100, 1.0, 1.0)
    xy = ToolsArray([x, y], (Var(1:2),))
    xx = ToolsArray([x, x], (Var(1:2),))

    @test_nowarn TimeseriesTools.STOIC([x, y]) isa Matrix
    @test STOIC(xy) == STOIC([x, y])
    STOIC(xx) == ones(length(xx), length(xx))
    @test STOIC(xx) isa AbstractToolsArray
    @test dims(xx, 1) == dims(STOIC(xx), 1)
end

@testset "STOIC distance" begin
    x = [gammarenewal(100, 1.0, 1.0) for _ in 1:10]
    @test pairwise(StoicDist(; σ = 1.0, Δt = 1.0), x) isa Matrix

    x = ToolsArray(x, Var(1:length(x)))
    A = pairwise(StoicDist(), x)
    B = pairwise(StoicDist(), x, x)
    @test A == B
end

@testset "Spike FFT" begin
    ts = 0:0.01:100
    t = [abs(_t - round(_t)) < 0.05 ? 1 : 0 for _t in ts][1:(end - 1)]
    t = findall(t .> 0) ./ 100 # Should have a period of 1 second
    t = TimeSeries(t, trues(length(t)))
    @test t isa SpikeTrain

    p = @test_nowarn spikefft(0:0.1:10, t)
    @test_nowarn spikefft(0.1, times(t), Val(:schild))
    fs = (0.01, 50)
    e = @test_nowarn energyspectrum(t, fs; method = :schild)
    @test (2 * sum(e[2:end]) + e[1]) * fs[1] ≈ sum(t)
    p = @test_nowarn powerspectrum(t, fs; method = :schild)

    # Test this returns an identical result for spikes measured at regular intervals
    x = TimeSeries(ts, zeros(length(ts)))
    x[𝑡(Near(times(t)))] .= 1.0 / sqrt(samplingperiod(x))
    et = energyspectrum(x, 0.01)
    @test (2 * sum(et[2:end]) + et[1]) * fs[1] ≈ sum(t)

    # if false # Yep works, better, even
    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, et, color = :crimson)
    lines!(ax, e)
    f
    # end

    # Multivariate
    T = hcat(Var(1:4), t, t, t, t)
    P = @test_nowarn powerspectrum(T, fs; method = :schild)
    @test P[:, 1] == p

    # * Autocovariance spectrum
    p = energyspectrum(t, fs; method = stoic(; σ = 0.01))
    @test (2 * sum(p[2:end]) + p[1]) * fs[1] ≈ sum(t)

    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, decompose(et)..., color = :crimson)
    lines!(ax, decompose(p)...)
    f
end

@testset "Spike-time tiling coefficient" begin
    using IntervalSets
    using LinearAlgebra
    using Distributions
    ts = 0:0.01:100
    t = [abs(_t - round(_t)) < 0.05 ? 1 : 0 for _t in ts][1:(end - 1)]
    t = findall(t .> 0) ./ 100 # Should have a period of 1 second
    t = TimeSeries(t, trues(length(t)))
    Δt = 0.025
    c = @test_nowarn sttc(t, t .+ 0.02; Δt)

    t1 = rand(0 .. 10, 200)
    η = 8000.0
    t2 = t1 .+ η .* randn(length(t1)) * Δt
    sort!.([t1, t2])
    @test_nowarn sttc(t1, t2; Δt)

    # Positive semi-definite? Not quite, but almost.
    η = 1
    t1 = findall(rand(Poisson(0.01), 100000) .> 0) ./ 1000
    t2 = findall(rand(Poisson(0.05), 100000) .> 0) ./ 1000
    t3 = vcat(t1 .* η .* randn(length(t1)) * Δt, t2 .* η .* randn(length(t2)) * Δt)
    t4 = vcat(t2 .* η .* randn(length(t2)) * Δt, t1 .* η .* randn(length(t1)) * Δt)
    t5 = vcat(t1 .* η .* randn(length(t1)) * Δt, t1 .* η .* randn(length(t1)) * Δt)
    t6 = vcat(t2 .* η .* randn(length(t2)) * Δt, t2 .* η .* randn(length(t2)) * Δt)
    sort!.([t3, t4, t5, t6])

    ts = [t1, t2, t3, t4, t5, t6]

    is = Iterators.product(ts, ts)
    Λ = [sttc(i...) for i in is]
    λ = eigvals(Λ)
end

@testset "Spike-time overlap-integral coefficient (stoic)" begin
    x = randn(1000) |> sort
    y = randn(1000) |> sort
    σ = 0.25
    Δt = σ * 10

    # * Verify the integral of the product of two gaussians
    G = TimeseriesTools.normal
    t1 = 0.025
    t2 = 0.0
    ff(x) = G(t1, σ)(x) * G(t2, σ)(x)
    I1 = sum(ff.(-1:0.001:1)) * 0.001
    I2 = G(t1, sqrt(2) * σ)(t2)
    @test I1≈I2 rtol=1e-6
    # Aw yeah

    D = @test_nowarn closeneighbours(x, y; Δt)
    @test stoic(x, y; Δt, σ)≈1.0 rtol=5e-2

    x = y
    @test stoic(x, y; Δt, σ) == 1.0

    # * Verify the normalization
    @test sum(TimeseriesTools.normal(1.0).(-10:0.01:10)) .* 0.01 ≈ 1
    @test sum(TimeseriesTools.npi(0.2).(-10:0.01:10)) .* 0.01 ≈ 1
    # a = x
    # E1 = stoic(a, a; Δt = σ * 100, σ, normalize = false)
    # E2 = stoic(a, a; Δt = σ, σ = σ / 1000, normalize = false)

    x = rand(0 .. 1000, 1000) |> sort
    y = rand(0 .. 1000, 1000) |> sort
    σ = 100
    Δt = σ * 10
    @test stoic(x, y; Δt, σ)≈1.0 rtol=0.02
    σ = 0.001
    Δt = σ * 10
    @test stoic(x, y; σ)≈0.0 atol=1e-2

    @test stoic([0.0], [0.0]; σ = 1, normalize = false) ≈ 0.5 / sqrt(π)
    @test stoic([0.0], [1.0]; σ = 1, normalize = false) ≈ 1 / (2 * exp(1 / 4) * sqrt(π))
    @test stoic([0.0, 10.0], [1.0, 10.0]; σ = 1, normalize = false)≈0.50179 rtol=1e-4

    # * Is it positive semi-definite?
    x = [rand(0 .. 100, 100) |> sort for _ in 1:100]
    [_x .= x[22] for _x in x[round.(Int, rand(1 .. length(x), 20))]]
    [_x .= sort(x[22] .+ 0.01 .* randn(100))
     for _x in x[round.(Int, rand(1 .. length(x), 20))]]
    [_x .= sort(x[22] .+ 0.05 .* randn(100))
     for _x in x[round.(Int, rand(1 .. length(x), 20))]]
    ρ = @test_nowarn pairwise(stoic(; σ = 0.01), x)
    e = eigvals(ρ)
    @test minimum(real.(e)) + 1e-10 > 0.0
    @test all(isapprox.(imag.(e), 0.0; atol = 1e-10))
end

@testset "Stoic spike-train length" begin
    # * Set up independent gamma renewal processes and verify stoic scaling with length vs.
    #   kernel width
    Ns = range(start = 100, step = 100, length = 100)
    xs = [gammarenewal(N, 1, 1) for N in Ns]
    ρ = @test_nowarn pairwise(stoic(; σ = 0.01), xs; symmetric = true)
    @test mean(ρ[ρ .!= 1])≈0 atol=0.05
    ρ[ρ .== 1] .= NaN

    f = Figure()
    ax = Axis(f[1, 1]; aspect = 1, xlabel = "N₁", ylabel = "N₂")
    p = heatmap!(ax, Ns, Ns, ρ)
    Colorbar(f[1, 2], p, label = "stoic")
end
@testset "Stoic spike-train fano" begin
    # * Set up independent gamma renewal processes and verify stoic scaling with length vs.
    #   kernel width
    θs = range(start = 0.1, step = 0.01, length = 150)
    xs = [gammarenewal(10000, 1, θ) for θ in θs]
    ρ = pairwise(stoic(; σ = 0.01), xs; symmetric = true)
    ρ[ρ .== 1] .= NaN

    f = Figure(size = (720, 360))
    ax = Axis(f[1, 1]; aspect = 1, xlabel = "θ₁", ylabel = "θ₂")
    p = heatmap!(ax, θs, θs, ρ)
    Colorbar(f[1, 2], p, label = "stoic")
    f

    ρ2 = pairwise(sttc(; Δt = 0.03), xs; symmetric = true)
    ρ2[ρ2 .== 1] .= NaN

    ax = Axis(f[1, 3]; aspect = 1, xlabel = "θ₁", ylabel = "θ₂")
    p = heatmap!(ax, θs, θs, abs.(ρ2))
    Colorbar(f[1, 4], p, label = "|sttc|")

    rowsize!(f.layout, 1, Relative(0.6))
    f
end

# @testset "Stoic firing rate and fano factor" begin N = 5000

#     fr = range(start = 0.1, stop = 500, length = 100) # 1 the mean ISI θs = range(start =
#     0.1, stop = 2, length = 100) # 1 the mean ISI αs = 1 ./ (fr .* θs)

#     ps = Iterators.product(αs, θs)

#     ρ = zeros(length(fr), length(θs)) n = 50 Threads.@threads for i in 1:n _ρ = map(ps) do
#     (α, θ) x = gammarenewal(N, α, θ) y = gammarenewal(N, α, θ) # sttc(x, y; Δt = 0.025) |>
#     abs stoic(x, y; σ = 0.01) |> abs end ρ .+= _ρ end ρ ./= n # Mean

#     f = Figure()
#     ax = Axis(f[1, 1]; aspect = 1, xlabel = "Mean firing rate", ylabel = "Fano factor")
#     p = heatmap!(ax, fr, θs, ρ)
#     Colorbar(f[1, 2], p, label = "stoic")
#     f
# end
