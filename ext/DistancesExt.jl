# module DSPExt

using TimeseriesTools
using .Distances

struct StoicDist <: Metric
    σ::Real
    Δt::Real
    normalize::Bool
    kpi::Function
end
function StoicDist(; kpi = npi, σ = 0.025, Δt = σ * 10, normalize = true)
    StoicDist(σ, Δt, normalize, kpi)
end
function (D::StoicDist)(a::AbstractVector, b::AbstractVector)
    1 - stoic(a, b; σ = D.σ, Δt = D.Δt, kpi = D.kpi, normalize = D.normalize)
end # A distance, like cosine distance
export StoicDist

function Distances.pairwise(metric::PreMetric, a::AbstractDimVector)
    rebuild(a; data = pairwise(metric, parent(a)), dims = (dims(a, 1), dims(a, 1)))
end

function Distances.pairwise(metric::PreMetric, a::AbstractDimVector, b::AbstractDimVector)
    rebuild(a; data = pairwise(metric, parent(a), parent(b)),
            dims = (dims(a, 1), dims(b, 1)))
end

function Distances.result_type(metric::StoicDist, a::AbstractVector{<:SpikeTrain},
                               b::AbstractVector{<:SpikeTrain})
    eltype(lookup(first(a), 𝑡)) # ! Should check if all lookups have the same eltype
end
# end # module
