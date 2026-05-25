module DistributionsExt

using TimeseriesTools
using Random
import Distributions
import Distributions: Distribution, Gamma
import TimeseriesTools: pointprocess!, gammarenewal!, gammarenewal

function pointprocess!(spikes, D::Distribution; rng = Random.default_rng(), t0 = 0.0)
    N = length(spikes)
    t = deepcopy(t0)
    isis = rand(rng, D, N)
    for i in 1:N
        t += isis[i]
        spikes[i] = t
    end
    return
end

function gammarenewal!(
        spikes::AbstractVector, α, θ;
        t0 = randn(Random.default_rng()) * α * θ, kwargs...
    )
    D = Gamma(α, θ)
    return pointprocess!(spikes, D; t0, kwargs...)
end

function gammarenewal(N, args...; kwargs...)
    spikes = zeros(Float64, N)
    gammarenewal!(spikes, args...; kwargs...)
    return spiketrain(spikes)
end
function gammarenewal(spikes::SpikeTrain, args...; kwargs...)
    t = collect(times(spikes))
    gammarenewal!(t, args...; kwargs...)
    return SpikeTrain(t)
end

end # module
