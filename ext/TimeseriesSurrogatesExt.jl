# module TimeseriesSurrogatesExt
import ..TimeseriesSurrogates
import ..TimeseriesSurrogates: Surrogate, SurrogateGenerator, surrogenerator
using TimeseriesTools
using Statistics
using Distributions
using DimensionalData

export RandomJitter, GammaRenewal, NDFT, phaserand!

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

# NDFT surrogates

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

# * Only 100% accurate for ODD sized arrays
function phaserand!(ϕ, rng = Random.default_rng(), n = size(ϕ))
    if any(iseven.(n))
        ds = findall(iseven, n)
        Is = collect.(axes(ϕ))
        for d in ds
            idxs = collect(Any, axes(ϕ))
            idxs[d] = n[d] ÷ 2 + 1
            _ϕ = view(ϕ, idxs...)
            # phaserand!(_ϕ)
            popat!(Is[d], n[d] ÷ 2 + 1)
        end
        ϕ = view(ϕ, Is...)
    end
    fs = fftfreq.(size(ϕ))
    for _fs in Iterators.product(fs...)
        idx = map(_fs, fs) do _f, f
            findfirst(_f .== f)
        end
        any(isnothing.(idx)) && continue
        ϕ[idx...] = rand(rng, Uniform(-π, π))
        idx₋ = map(_fs, fs) do _f, f
            findfirst(-_f .== f)
        end
        any(isnothing.(idx₋)) && continue
        ϕ[idx₋...] = -ϕ[idx...]
    end
end

struct NDFT <: Surrogate
end

function surrogenerator(x, method::NDFT, rng = Random.default_rng())
    n = size(x)
    m = mean(x)
    forward = plan_fft(x)
    inverse = plan_ifft(forward * x)
    𝓕 = forward * (x .- m)

    init = (inverse = inverse,
            m = m,
            𝓕 = 𝓕,
            r = abs.(𝓕),
            ϕ = angle.(𝓕),
            shuffled𝓕 = similar(𝓕),
            coeffs = zeros(size(𝓕)),
            n = n)

    return SurrogateGenerator(method, x, similar(x), init, rng)
end

function (sg::SurrogateGenerator{<:NDFT})()
    s, rng = sg.s, sg.rng

    init_fields = (:inverse, :m, :r, :ϕ, :shuffled𝓕, :coeffs, :n)
    inverse, m, r, ϕ, shuffled𝓕, coeffs, n = getfield.(Ref(sg.init),
                                                       init_fields)
    coeffs .= ϕ
    phaserand!(coeffs, rng)
    shuffled𝓕 .= r .* exp.(coeffs .* 1im)
    _s = inverse * shuffled𝓕
    @assert all(isapprox.(imag.(_s), 0; atol = 1e-3))
    s .= real.(_s) .+ m
    return s
end

# end # module
