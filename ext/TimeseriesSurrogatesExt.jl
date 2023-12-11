# module TimeseriesSurrogatesExt
import ..TimeseriesSurrogates
import ..TimeseriesSurrogates: Surrogate, SurrogateGenerator, surrogenerator
using TimeseriesTools
using Statistics
using Distributions

export RandomJitter, GammaRenewal

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
# end # module
