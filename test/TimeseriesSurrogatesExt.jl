@testset "ND phase randomization" begin
    f = xy -> sin.(0.5 * 2π * sum(xy)) + cos.(0.1 * 2π * xy[2])

    # * Odd
    x = f.(Iterators.product(range(0, 1, length=5), range(0, 1, length=5)))
    ϕ = angle.(fft(x))
    phaserand!(ϕ)
    @test sum(fftshift((ϕ)) .+ reverse(fftshift((ϕ))) .== 0) == length(ϕ) - 1

    x = randn(11, 11, 11)
    ϕ = angle.(fft(x))
    phaserand!(ϕ)
    @test sum(fftshift((ϕ)) .+ reverse(fftshift((ϕ))) .== 0) == length(ϕ) - 1

    # * Even
    x = f.(Iterators.product(range(0, 1, length=6), range(0, 1, length=6)))
    ϕ = angle.(fft(x))
    phaserand!(ϕ)
    matchn = sum(fftshift((ϕ))[2:end, 2:end] .+ reverse(fftshift((ϕ))[2:end, 2:end]) .== 0)
    @test matchn == length(ϕ[2:end, 2:end]) - 1

    x = randn(10, 10, 10)
    ϕ = angle.(fft(x))
    phaserand!(ϕ)
    @test sum(fftshift((ϕ))[2:end, 2:end, 2:end] .+
              reverse(fftshift((ϕ))[2:end, 2:end, 2:end]) .== 0) ==
          length(ϕ[2:end, 2:end, 2:end]) - 1

    # * Mixed
    x = f.(Iterators.product(range(0, 1, length=6), range(0, 1, length=5)))
    ϕ = angle.(fft(x))
    phaserand!(ϕ)
    @test sum(fftshift((ϕ))[2:end, :] .+ reverse(fftshift((ϕ))[2:end, :]) .== 0) ==
          length(ϕ[2:end, :]) - 1

    x = randn(11, 10, 11)
    ϕ = angle.(fft(x))
    phaserand!(ϕ)
    @test sum(fftshift((ϕ))[:, 2:end, :] .+ reverse(fftshift((ϕ))[:, 2:end, :]) .== 0) ==
          length(ϕ[:, 2:end, :]) - 1
end

@testset "1D ND surrogates" begin
    x = loadtimeseries("./test_timeseries.tsv")[:, 1]
    x = bandpass(x, 1000, [1, 20])
    S = abs.(fft(x)) .^ 2
    s = spectrum(rectify(x, dims=𝑡))

    x̂ = deepcopy(x)
    x̂ .= surrogate(collect(x), FT())
    Ŝ = abs.(fft(x̂)) .^ 2
    ŝ = spectrum(rectify(x̂, dims=𝑡))
    @test length(Ŝ) == length(S)
    @test sum(abs.(S .- Ŝ)) ./ sum(S) ≈ 0 atol = 1e-9
    @test sum(abs.(s .- ŝ)) ./ sum(s) ≈ 0 atol = 1e-1

    x̂ = deepcopy(x)
    x̂ .= surrogate(collect(x), NDFT())
    Ŝ = abs.(fft(x̂)) .^ 2
    ŝ = spectrum(rectify(x̂, dims=𝑡))
    @test length(Ŝ) == length(S)
    @test sum(abs.(S .- Ŝ)) ./ sum(S) ≈ 0 atol = 1e-9
    @test sum(abs.(s .- ŝ)) ./ sum(s) ≈ 0 atol = 1e-1

    x̂ = deepcopy(x)
    x̂ .= surrogate(collect(x), NDAAFT())
    ŝ = spectrum(rectify(x̂, dims=𝑡))
    @test sum(abs.(s .- ŝ)) ./ sum(s) ≈ 0 atol = 0.15

    x̂ = deepcopy(x)
    x̂ .= surrogate(collect(x), NDIAAFT())

    ŝ = spectrum(rectify(x̂, dims=𝑡))
    @test sum(abs.(s .- ŝ)) ./ sum(s) ≈ 0 atol = 0.1
end

@testset "2D ND surrogates" begin
    f = xy -> sin.(2 * 2π * sum(xy)) + cos.(1 * 2π * xy[2])

    # * Odd
    x = f.(Iterators.product(range(0, 1, length=101), range(0, 1, length=101)))

    S = abs.(fft(x)) .^ 2

    ϕ = angle.(fft(x))
    phaserand!(ϕ)
    @test sum(fftshift((ϕ)) .+ reverse(fftshift((ϕ))) .== 0) == length(ϕ) - 1

    x̂ = deepcopy(x)
    x̂ .= surrogate(collect(x), NDFT())
    Ŝ = abs.(fft(x̂)) .^ 2
    @test length(Ŝ) == length(S)
    @test sum(abs.(S .- Ŝ)) ./ sum(S) ≈ 0 atol = 1e-9

    x̂ = deepcopy(x)
    x̂ .= surrogate(collect(x), NDAAFT())
    Ŝ = abs.(fft(x̂)) .^ 2
    @test length(Ŝ) == length(S)
    @test sum(abs.(S .- Ŝ)) ./ sum(S) ≈ 0 atol = 2e-3

    x̂ = deepcopy(x)
    x̂ .= surrogate(collect(x), NDIAAFT())
    Ŝ = abs.(fft(x̂)) .^ 2
    @test length(Ŝ) == length(S)
    @test sum(abs.(S .- Ŝ)) ./ sum(S) ≈ 0 atol = 2e-3

    x̂ = deepcopy(x)
    x̂ .= surrogate(collect(x), MVFT())
    Ŝ = abs.(fft(x̂)) .^ 2
    @test length(Ŝ) == length(S)
    @test sum(abs.(S .- Ŝ)) ./ sum(S) ≈ 0 atol = 1e-10

    # * Even
    x = f.(Iterators.product(range(0, 1, length=4), range(0, 1, length=4)))

    S = abs.(fft(x)) .^ 2

    ϕ = angle.(fft(x))
    phaserand!(ϕ)
    @test sum(fftshift((ϕ))[2:end, 2:end] .+ reverse(fftshift((ϕ))[2:end, 2:end]) .== 0) ==
          length(ϕ[2:end, 2:end]) - 1 # The zero frequency phase should be non-zero, although this doesn't matter

    x̂ = deepcopy(x)
    x̂ .= surrogate(collect(x), NDFT())
    Ŝ = abs.(fft(x̂)) .^ 2
    @test length(Ŝ) == length(S)
    @test sum(abs.(S .- Ŝ)) ./ sum(S) ≈ 0 atol = 1e-9
end

@testset "ND Fourier transform surrogates" begin
    xs = -0.6:0.01:0.6
    x = [stack(X(xs), [colorednoise(0:0.01:1) for _ in xs]) for _ in xs]
    x = stack(Y(xs), x)

    S = abs.(rfft(x)) .^ 2

    x̂ = deepcopy(x)
    x̂ .= surrogate(x, NDFT())

    Ŝ = abs.(rfft(x̂)) .^ 2

    @test S ≈ Ŝ rtol = 1e-10

    # * Larger array, smoothed
    xs = -0.6:0.01:0.6
    x = [stack(X(xs), [colorednoise(0:0.05:50) for _ in xs]) for _ in xs]
    x = stack(Y(xs), x)

    function G(x, μ, σ)
        d = length(μ)
        exponent = -0.5 * dot(x .- μ, (x .- μ) ./ σ)
        coeff = 1 / ((2π)^(d / 2) * σ)
        return coeff * exp(exponent)
    end
    G(μ, σ) = x -> G(x, μ, σ)

    G1(x) = G(x, [-0.25, -0.25], 0.05)
    G2(x) = G(x, [0.25, 0.25], 0.05)
    M1 = G1.(Iterators.product(lookup(x)[2:3]...))
    M2 = G2.(Iterators.product(lookup(x)[2:3]...))
    # x[X = 0 .. 0.5] .= reverse(x[X = 0 .. 0.5]) # A little boundary
    _x = 5.0 * mean(std(x, dims=1)) .* sin.(times(x))  # A slowly varying "true" signal
    _x = zeros(size(x, 1), 1, 1) .+ _x
    MM1 = permutedims(stack([M1 for _ in 1:size(x, 1)]), [3, 1, 2])
    MM1 = mapslices(x -> x .* _x, MM1, dims=1)

    _x = 10.0 * mean(std(x, dims=1)) .* cos.(times(x .* 1.5))  # A slowly varying "true" signal
    _x = zeros(size(x, 1), 1, 1) .+ _x
    MM2 = permutedims(stack([M2 for _ in 1:size(x, 1)]), [3, 1, 2])
    MM2 = mapslices(x -> x .* _x, MM2, dims=1)

    x = x .+ MM1 .+ MM2

    fg(x) = bandpass(x, 1 / step(xs), (1, 5))
    x = mapslices(fg, x, dims=2)
    x = mapslices(fg, x, dims=3)

    x = bandpass(x, 0.1 .. 0.5)
    y = angle.(hilbert(x))

    x = x[X=-0.5 .. 0.5, Y=-0.5 .. 0.5]
    y = y[X=-0.5 .. 0.5, Y=-0.5 .. 0.5]
    # if !haskey(ENV, "CI")
    #     f = Figure()
    #     ax = Axis(f[1, 1])
    #     xx = Observable(x[𝑡= 1])
    #     heatmap!(ax, xx; colorrange = extrema(x))
    #     record(f, "./MultidimModel_x.mp4", 1:2:900) do i
    #         xx[] = parent(x[𝑡= i])
    #     end
    # end
    # if !haskey(ENV, "CI")
    #     f = Figure()
    #     ax = Axis(f[1, 1])
    #     xx = Observable(y[𝑡= 1])
    #     heatmap!(ax, xx; colorrange = extrema(y), colormap = :twilight)
    #     record(f, "./MultidimModel_phi.mp4", 1:2:900) do i
    #         xx[] = y[𝑡= i]
    #     end
    # end

    S = abs.(fft(x)) .^ 2

    x̂ = deepcopy(x)
    x̂ .= surrogate(x, NDFT())

    Ŝ = abs.(fft(x̂)) .^ 2

    # Ss = periodogram(collect(x[50, :, :]), radialavg = true)
    # plot(Ss.freq, Ss.power)
    # Ss = periodogram(collect(x̂[50, :, :]), radialavg = true)
    # plot!(Ss.freq, Ss.power)
    # current_figure()

    @test S ≈ Ŝ rtol = 1e-10

    x = bandpass(x̂, 0.1 .. 0.5)
    y = angle.(hilbert(x))

    x = x[X=-0.5 .. 0.5, Y=-0.5 .. 0.5]
    y = y[X=-0.5 .. 0.5, Y=-0.5 .. 0.5]
    # if !haskey(ENV, "CI")
    #     f = Figure()
    #     ax = Axis(f[1, 1])
    #     xx = Observable(x[𝑡= 1])
    #     heatmap!(ax, xx; colorrange = extrema(x))
    #     record(f, "./MultidimModel_x_s.mp4", 1:2:900) do i
    #         xx[] = x[𝑡= i]
    #     end
    # end
    # if !haskey(ENV, "CI")
    #     f = Figure()
    #     ax = Axis(f[1, 1])
    #     xx = Observable(y[𝑡= 1])
    #     heatmap!(ax, xx; colorrange = extrema(y), colormap = :twilight)
    #     record(f, "./MultidimModel_phi_s.mp4", 1:2:900) do i
    #         xx[] = y[𝑡= i]
    #     end
    # end
end
