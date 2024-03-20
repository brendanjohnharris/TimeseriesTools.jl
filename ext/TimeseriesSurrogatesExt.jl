# module TimeseriesSurrogatesExt
import ..TimeseriesSurrogates
import ..TimeseriesSurrogates: Surrogate, SurrogateGenerator, surrogenerator
using TimeseriesTools
using Statistics
using Distributions
using DimensionalData

export RandomJitter, GammaRenewal, NDFT

# Spike train surrogates
struct RandomJitter <: Surrogate
    Δt::Real # The minimum jitter
    σ::Real # The Standard deviation of a half-normal describing the jitter distribution
end

RandomJitter(; Δt, σ) = RandomJitter(Δt, σ)

function surrogenerator(x::AbstractVector, rf::RandomJitter, rng = Random.default_rng())
    D = truncated(Normal(rf.Δt, rf.σ), rf.Δt, nothing)
    init = (; D)

    return SurrogateGenerator(rf, x, deepcopy(x), init, rng)
end

function (sg::SurrogateGenerator{<:RandomJitter})()
    _, s, rng = sg.x, sg.s, sg.rng
    s .+= rand(rng, [-1, 1], length(s)) .* rand(rng, sg.init.D, length(s))
    sort!(s)
    return s
end

struct GammaRenewal <: Surrogate end

function surrogenerator(x::AbstractVector, rf::GammaRenewal, rng = Random.default_rng())
    dt = diff(x)
    μ = mean(dt)
    θ = var(dt) / μ
    α = μ / θ
    D = Distributions.Gamma(α, θ)
    init = (; D)

    return SurrogateGenerator(rf, x, similar(x), init, rng)
end

function (sg::SurrogateGenerator{<:GammaRenewal})()
    x, s, rng = sg.x, sg.s, sg.rng
    pointprocess!(s, sg.init.D; rng)
    s .+= first(x) # To roughly align time itnervals
    return s
end

"""
Generate a surrogate for a multi-dimensional time series `X` using the multidimensional
phase-randomisation algorithm. If the array contains any boundary NaNs, the Fourier
transform subsamples the largest non-NaN rectangular array.
"""
struct NDFT <: Surrogate
end

function surrogenerator(x::AbstractArray{<:Real}, rf::NDFT, rng = Random.default_rng())
    # @assert !any(Bool.(mod.(size(x), 2))) # All even Assumes the largest valid
    # subrectangle is at least half the size of the array in each dim.
    NaNs = isnan.(x)
    any(NaNs) && (x = nansubarray(x))
    m = mean(x)
    σ = std(x) # Weird scaling issue
    oddsizes = size(x) .- Bool.(mod.(size(x), 2))
    x = getindex(x, [1:n for n in oddsizes]...)
    # ! We get periodic funniness here; better to extend by 1?

    ffreqs = fftfreq.(size(x))
    ds = [Dim{Symbol(i)}(ix) for (i, ix) in enumerate(ffreqs)]
    F = fft(x .- m)
    F = DimArray(F, Tuple(ds))

    s = zeros(2 .* size(x)) # Double the size of the surrogate array to have the correct frequencies but a buffer for resizing. Pretty slow, but what do you want?
    inverse = plan_ifft(s)
    sfreqs = fftfreq.(size(s))
    ds = [Dim{Symbol(i)}(ix) for (i, ix) in enumerate(sfreqs)]
    shuffledF = DimArray(zeros(Complex, length.(sfreqs)), Tuple(ds))

    n = size(F)
    r = abs.(F)
    ϕ = angle.(F)

    init = (; inverse, m, n, r, ϕ, shuffledF, NaNs, σ)
    return SurrogateGenerator(rf, x, s, init, rng)
end

function (sg::SurrogateGenerator{<:NDFT})()
    inverse, m, n, r, ϕ, shuffledF, NaNs, σ = getfield.(Ref(sg.init),
                                                        (:inverse, :m, :n, :r, :ϕ,
                                                         :shuffledF, :NaNs, :σ))
    s, rng = sg.s, sg.rng

    phaserand!(ϕ, rng)

    shuffledF[At.(dims(ϕ))...] .= r .* exp.(ϕ .* 1im)

    _s = inverse * shuffledF
    @assert all(isapprox.(real(_s), _s; atol = 1e-6))
    s .= σ .* real(_s) ./ std(real(_s)) .+ m
    s = s[axes(parent(NaNs))...] # Crop to the original size. Doesn't really matter where we crop, the surrogate is totally stationary
    s[NaNs] .= NaN

    return s
end

function nansubarray(X::AbstractMatrix{<:AbstractFloat})
    # Get a rectangular sub array that contains no nans
    Y = similar(X)
    idxs = isnan.(X)
    Y[idxs] .= -Inf
    Y[.!idxs] .= 1.0
    _, a, b = maxrect(Y)
    return X[a[1]:b[1], a[2]:b[2]]
end

function nansubarray(X::AbstractArray{<:AbstractFloat, 3})
    nansum = sum(isnan.(X), dims = 3)
    @assert all(nansum .∈ ([0, size(X, 3)],))
    Y = deepcopy(X[:, :, 1])
    idxs = isnan.(Y)
    Y[idxs] .= -Inf
    Y[.!idxs] .= 1.0
    _, a, b = maxrect(Y)
    return X[a[1]:b[1], a[2]:b[2], :]
end;
export nansubarray

function kadane!(start, fin, x::AbstractVector)
    x = deepcopy(x)
    S = 0
    maxS = -Inf
    tempStart = 1
    for (i, _x) in enumerate(x)
        S += _x
        if S < 0
            S = 0
            tempStart = i + 1
        elseif S > maxS
            maxS = S
            start .= tempStart
            fin .= i
        end
    end
    return maxS
end;

function maxrect(X::AbstractArray)
    maxS = -Inf
    tmp = ones(Float64, length(X[:, 1]))
    start = [1]
    fin = [1]
    endl = endr = endt = endb = 1
    for l in 1:lastindex(tmp)
        tmp .= 1
        for r in l:lastindex(X, 2)
            for i in 1:lastindex(tmp)
                tmp[i] += X[i, r]
            end
            S = kadane!(start, fin, tmp)
            if S > maxS
                maxS = S
                endl = l
                endr = r
                endt = start[1]
                endb = fin[1]
            end
        end
    end
    return maxS, (endt, endl), (endb, endr)
end

function phaserand!(ϕ, rng = Random.default_rng()) # ! Works!!! But only for gridded data
    length(ϕ) == 1 && return # This will correspond to (0, 0, ...)
    i2f(x) = [i - Int(size(ϕ, n) ÷ 2 + 1) for (n, i) in enumerate(x)]
    f2i(x) = [i + Int(size(ϕ, n) ÷ 2 + 1) for (n, i) in enumerate(x)]
    # First do the edges, recursively
    N = ndims(ϕ)
    edg = [fill(Colon(), N) for _ in 1:N] .|> Array{Union{Colon, Int}}
    [edg[n][n] = 1 for n in 1:N]
    for e in edg
        _ϕ = view(ϕ, e...)
        phaserand!(_ϕ, rng)
    end
    # Then the remainder
    for i in CartesianIndices(ϕ)
        i = Tuple(i)
        if !any(i .== 1)
            if all(i2f(i) .== 0)
                if mod(ϕ[i...], π) ≈ 0 || mod(ϕ[i...], π) ≈ π
                    # Then pass. These zero frequencies we don't want to touch.
                else
                    @warn "The zero-frequency phases are neither 0 nor π. Leaving untouched since I don't know the symmetry here."
                end
            else
                _i = f2i(.-i2f(i))
                ϕ[_i...] = rand(rng, Uniform(0, 2π))
                ϕ[i...] = -ϕ[_i...] # Phase symmetry
            end
        end
    end
end

# end # module
