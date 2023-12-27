using SparseArrays
using Random
using Distributions
export spikefft, sttc, convolve, closeneighbours, stoic, pointprocess!, gammarenewal!,
       gammarenewal

normal(σ) = x -> (1 / (σ * sqrt(2π))) .* exp.(-0.5 .* x .^ 2 ./ σ^2)
normal(μ, σ) = x -> (1 / (σ * sqrt(2π))) .* exp.(-0.5 .* (x .- μ) .^ 2 ./ σ^2)
npi(σ) = normal(sqrt(2) * σ) # The integral of the product of two gaussians with separation `x` and equal variance σ²
function convolve(t::SpikeTrain; kernel::Function, range = 0.0)
    @assert all(t .== true)
    fs = [(x -> kernel(x .- _t)) for _t in times(t)]
    isempty(fs) && return x -> 0.0
    if range > 0
        function f(x)
            inrange = [abs(x - b) < range for b in times(t)] # Only consider spikes that are within a reasonable time of one another; the others should be negligible if the kernel decays
            _x = [g(x) for (i, g) in enumerate(fs) if inrange[i]]
            isempty(_x) && return x -> 0.0
            sum(_x)
        end
    else
        function h(x)
            _x = [g(x) for g in fs]
            isempty(_x) && return 0.0
            sum(_x)
        end
    end
end

function convolve(t::SpikeTrain, p; kernel::Function = normal, kwargs...)
    convolve(t; kernel = kernel(p), kwargs...)
end
"""
    sttc(a, b; Δt = 0.025)

The spike-time tiling coefficient, a measure of correlation between spike trains [1].

# Arguments
- `a::Vector{<:Real}`: A sorted vector of spike times.
- `b::Vector{<:Real}`: A second sorted vector of spike times .
- `Δt::Real=0.025`: The time window for calculating the STTC.

# Returns
- `sttc::Real`: The STTC value.

# References
    [1] [Cutts & Eglen 2014](https://doi.org/10.1523%2FJNEUROSCI.2767-14.2014)
"""
function sttc(a, b; Δt = 0.025)
    if !issorted(a) || !issorted(b)
        error("Spike trains must be sorted")
    end

    Ta = 0
    ext = 0
    for _a in a
        Ta += min(_a + Δt - ext, 2 * Δt) # If the window overlaps the previous window, add the remainder. Otherwise, add the full window
        ext = _a + Δt
        # Assume the first and last spikes with their overhanging windows are negligible
    end
    Ta = Ta / (last(a) - first(a) + 2 * Δt)
    Tb = 0
    ext = 0
    for _b in b
        Tb += min(_b + Δt - ext, 2 * Δt)
        ext = _b + Δt
    end
    Tb = Tb / (last(b) - first(b) + 2 * Δt)

    i = 1 # Keep track of which spikes are behind us
    Na = 0
    for _a in a
        while _a > b[i] + Δt && i < length(b)
            i += 1
        end
        if b[i] - Δt < _a ≤ b[i] + Δt
            Na += 1
        end
    end
    i = 1
    Nb = 0
    for _b in b
        while _b > a[i] + Δt && i < length(a)
            i += 1
        end
        if a[i] - Δt < _b ≤ a[i] + Δt
            Nb += 1
        end
    end
    Pa = Na / length(a)
    Pb = Nb / length(b)
    return 0.5 * ((Pa - Tb) / (1 - Pa * Tb) + (Pb - Ta) / (1 - Pb * Ta))
end

function sttc(a::UnivariateTimeSeries, b::UnivariateTimeSeries; τ = 0.0, kwargs...)
    if τ != 0.0
        b = 𝒯(τ)(b)
    end
    sttc(times(a), times(b); kwargs...)
end
sttc(; kwargs...) = (x, y) -> sttc(x, y; kwargs...)

function mapneighbours!(x, y, f!; Δt)
    if !issorted(x) || !issorted(y)
        error("Spike trains must be sorted")
    end

    # Iterate through the train with the smallest number of spikes, looking for neighbours
    c = length(y) > length(x)
    a = c ? x : y
    b = c ? y : x
    la = length(a)
    lb = length(b)

    _j = 1 # Keep track of which spikes are behind us
    j = 1
    for i in eachindex(a)
        while _j < lb && b[_j] < a[i] - Δt # Catch up to this window
            _j += 1
        end
        j = _j # Catch up
        while j ≤ lb && a[i] - Δt ≤ b[j] ≤ a[i] + Δt
            f!(a[i], b[j], i, j)
            j += 1
        end
    end
end

"""
    closeneighbours(x, y; Δt)

Constructs a sparse matrix of distances between neighbouring spikes in two sorted spike trains.

# Arguments
- `x`: A sorted array representing the first spike train.
- `y`: A sorted array representing the second spike train.
- `Δt`: The maximum time difference allowed for two spikes to be considered neighbours.

# Returns
A sparse matrix `D` where `D[i, j]` represents the distance between the `i`-th spike in `x` and the `j`-th spike in `y`, for pairs of spikes within `Δt` of each other.
"""
function closeneighbours(x::AbstractVector{T}, y::AbstractVector{T};
                         kwargs...) where {T <: Real}
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{T}()
    function f!(a, b, i, j)
        push!(V, abs(a - b))
        push!(I, i)
        push!(J, j)
    end
    mapneighbours!(x, y, f!; kwargs...)
    lx = length(x)
    ly = length(y)
    D = ly > lx ? sparse(I, J, V, lx, ly) : sparse(J, I, V, lx, ly)
end

"""
    stoic(a, b; kpi = npi, σ = 0.025, Δt = σ * 10)

Compute the spike-train overlap-integral coefficient between two spike trains, after normalizing both convolutions to unit energy

# Arguments
- `a`: Spike train a.
- `b`: Spike train b.
- `kpi`: Kernel product integral, a function of the distance between two spikes. Default is `npi`, the integral of two gaussians with equal variance at a given distance from each other.
- `σ`: Width parameter of the kernel. For `npi`, this is the width of the unit-mass Gaussian kernels. Default is `0.025`.
- `Δt`: Time window for considering spikes as close. Default is `σ * 10`.
"""
function stoic(a, b; kpi = npi, σ = 0.025, Δt = σ * 10, normalize = true)
    if normalize
        𝐸a = stoic(a, a; kpi, σ, Δt, normalize = false)
        𝐸b = stoic(b, b; kpi, σ, Δt, normalize = false)
    else # Assume normalized
        𝐸a = 1.0
        𝐸b = 1.0
    end
    𝐶 = [0.0]
    function f!(a, b, i, j)
        𝐶[1] = 𝐶[1] + kpi(σ)(abs(a - b))
    end
    mapneighbours!(a, b, f!; Δt)
    𝐶[1] ./ sqrt(𝐸a * 𝐸b)
end

function stoic(a::UnivariateTimeSeries, b::UnivariateTimeSeries; τ = 0.0, kwargs...)
    if τ != 0.0
        b = 𝒯(τ)(b)
    end
    stoic(times(a), times(b); kwargs...)
end
stoic(; kwargs...) = (x, y) -> stoic(x, y; kwargs...)

"""
    pointprocess!(spikes, D::Distribution; rng = Random.default_rng(), t0 = 0.0)

Simulate a point process by sampling inter-spike intervals from a given distribution.

## Arguments
- `spikes`: An array to store the generated spike times.
- `D::Distribution`: The distribution from which to sample inter-spike intervals.
- `rng`: (optional) The random number generator to use. Defaults to `Random.default_rng()`.
- `t0`: (optional) The initial time. Defaults to `0.0`.
"""
function pointprocess!(spikes, D::Distribution; rng = Random.default_rng(), t0 = 0.0)
    N = length(spikes)
    t = deepcopy(t0)
    isis = rand(rng, D, N)
    for i in 1:N
        t += isis[i]
        spikes[i] = t
    end
end

"""
    gammarenewal!(spikes, α, θ; t0 = randn() * α / θ, kwargs...)

Generate a sequence of gamma-distributed renewal spikes.

Arguments:
- `spikes::AbstractVector`: The vector to store the generated spikes.
- `α`: The shape parameter of the gamma distribution.
- `θ`: The scale parameter of the gamma distribution.
- `t0`: The initial time of the spike train. Defaults to a random value drawn from a normal
  distribution with mean of 0 and standard deviation equal to the mean firing rate.
- `kwargs...`: Additional keyword arguments to be passed to [`pointprocess!`](@ref).
"""
function gammarenewal!(spikes::AbstractVector, α, θ;
                       t0 = randn(Random.default_rng()) * α * θ, kwargs...)
    N = length(spikes)
    D = Distributions.Gamma(α, θ)
    pointprocess!(spikes, D; t0, kwargs...)
end

"""
    gammarenewal(N, α, θ; t0)

Generate a spike train with inter-spike intervals drawn from a Gamma process.

# Arguments
- `N`: Number of spikes to generate.
- `α`: Shape parameter of the gamma distribution (equivalent to the mean ISI divided by the
 Fano factor).
- `θ`: Scale parameter of the gamma distribution (equivalent to the Fano factor).
- `t0`: The initial time of the spike train, prior to the first sampled spike. Defaults to a
 random value drawn from a normal distribution with mean of 0 and standard deviation equal
 to the mean firing rate.

# Returns
- A [`SpikeTrain`](@ref) containing the generated spike times.

See also [`gammarenewal!`](@ref).
"""
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
