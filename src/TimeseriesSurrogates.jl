# module TimeseriesSurrogatesExt
import TimeseriesSurrogates
import TimeseriesSurrogates: Surrogate, SurrogateGenerator, surrogenerator
# using TimeseriesTools
using Statistics
import Distributions: Gamma, Normal, Uniform, truncated
# using TimeseriesTools.Statistics
# using TimeseriesTools.Distributions
# using TimeseriesTools.DimensionalData

export RandomJitter, GammaRenewal, phaserand!, NDFT

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
    s .= map(real, parent(_s)) .+ m
    return s
end

struct NDAAFT <: Surrogate
end

function surrogenerator(x, method::NDAAFT, rng = Random.default_rng())
    n = size(x)
    m = mean(x)
    forward = plan_fft(x)
    inverse = plan_ifft(forward * x)
    𝓕 = forward * (x .- m)

    x_sorted = deepcopy(x)
    sort!(view(x_sorted, :))

    init = (inverse = inverse,
            x_sorted = x_sorted,
            ix = zeros(Int, size(x)),
            m = m,
            𝓕 = 𝓕,
            r = abs.(𝓕),
            ϕ = angle.(𝓕),
            shuffled𝓕 = similar(𝓕),
            coeffs = zeros(size(𝓕)),
            n = n)

    return SurrogateGenerator(method, x, similar(x), init, rng)
end

function (sg::SurrogateGenerator{<:NDAAFT})()
    s, rng = sg.s, sg.rng

    init_fields = (:inverse, :x_sorted, :ix, :m, :r, :ϕ, :shuffled𝓕, :coeffs, :n)
    inverse, x_sorted, ix, m, r, ϕ, shuffled𝓕, coeffs, n = getfield.(Ref(sg.init),
                                                                     init_fields)
    coeffs .= ϕ
    phaserand!(coeffs, rng)
    shuffled𝓕 .= r .* exp.(coeffs .* 1im)
    _s = inverse * shuffled𝓕
    @assert all(isapprox.(imag.(_s), 0; atol = 1e-3))
    s .= real.(_s) .+ m

    sortperm!(view(ix, :), view(s, :))
    s[ix] .= x_sorted
    return s
end

struct NDIAAFT <: Surrogate
    M::Int
    tol::Real

    function NDIAAFT(; M::Int = 100, tol::Real = 1e-9)
        new(M, tol)
    end
end

Base.show(io::IO, x::NDIAAFT) = print(io, "NDIAAFT(M = $(x.M), tol = $(x.tol))")

function surrogenerator(x, method::NDIAAFT, rng = Random.default_rng())
    n = size(x)
    m = mean(x)
    forward = plan_fft(x)
    inverse = plan_ifft(forward * x)
    𝓕 = forward * (x .- m)

    x_sorted = deepcopy(x)
    sort!(view(x_sorted, :))

    xpower = abs.(similar(𝓕)) .^ 2
    spower = copy(xpower)

    init = (forward = forward,
            inverse = inverse,
            x_sorted = x_sorted,
            ix = zeros(Int, size(x)),
            xpower = xpower,
            spower = spower,
            m = m,
            𝓕 = 𝓕,
            r = abs.(𝓕),
            ϕ = angle.(𝓕),
            shuffled𝓕 = similar(𝓕),
            coeffs = zeros(size(𝓕)),
            n = n)

    return SurrogateGenerator(method, x, similar(x), init, rng)
end

function spectraloss(x, y)
    sum(abs.((x .- y) ./ (x .+ y))) / 2
end

function (sg::SurrogateGenerator{<:NDIAAFT})()
    x, s, rng = sg.x, sg.s, sg.rng

    init_fields = (:forward, :inverse, :x_sorted, :ix, :xpower, :spower, :m, :r, :ϕ, :𝓕,
                   :shuffled𝓕,
                   :coeffs, :n)
    forward, inverse, x_sorted, ix, xpower, spower, m, r, ϕ, 𝓕, shuffled𝓕, coeffs, n = getfield.(Ref(sg.init),
                                                                                                 init_fields)
    M = sg.method.M
    tol = sg.method.tol

    n = length(x)
    view(s, :) .= x[sample(rng, 1:n, n)]

    sum_old, sum_new = 1.0, 0.0
    iter = 1
    while iter <= M
        @debug "$iter: $sum_old → $sum_new ($(sum_new - sum_old))"
        𝓕 .= forward * s
        ϕ .= angle.(𝓕)
        𝓕 .= r .* exp.(ϕ .* 1im)

        s .= real.(inverse * 𝓕) # Assume complex part is negligible
        sortperm!(view(ix, :), view(s, :))
        s[ix] .= x_sorted

        spower = abs.(𝓕) .^ 2 # ! Need an average power periodogram here...

        sum_new = spectraloss(x, s)

        if abs(sum_old - sum_new) ./ sum_old < tol
            iter = M + 1
        else
            sum_old = sum_new
        end

        iter += 1
    end

    return s
end

struct MVFT <: Surrogate
end

function surrogenerator(X::AbstractMatrix, rf::MVFT, rng = Random.default_rng())
    x = X[:, 1]
    forward = plan_rfft(x)
    inverse = plan_irfft(forward * x, length(x))
    m = mean(X, dims = 1)
    𝓕 = map(eachcol(X), m) do x, _m
        forward * (x .- _m)
    end
    𝓕 = stack(𝓕, dims = 2)
    shuffled𝓕 = zero(𝓕)
    S = similar(X)
    n = size(𝓕, 1)
    r = abs.(𝓕)
    ϕ = angle.(𝓕)
    coeffs = zero(r[:, 1])

    init = (inverse = inverse, m = m, coeffs = coeffs, n = n, r = r,
            ϕ = ϕ, shuffled𝓕 = shuffled𝓕)
    return SurrogateGenerator(rf, X, S, init, rng)
end

function (sg::SurrogateGenerator{<:MVFT})()
    inverse, m, coeffs, n, r, ϕ, shuffled𝓕 = getfield.(Ref(sg.init),
                                                       (:inverse, :m, :coeffs, :n, :r, :ϕ,
                                                        :shuffled𝓕))
    S, rng = sg.s, sg.rng

    rand!(rng, Uniform(0, 2π), coeffs)

    map(eachcol(S), eachcol(r), eachcol(ϕ), eachcol(shuffled𝓕), m) do s, _r, _ϕ, 𝓕, _m
        𝓕 .= _r .* exp.((coeffs .+ _ϕ) .* 1im)
        s .= inverse * 𝓕
        s .+= _m
    end
    return S
end

# end # module
