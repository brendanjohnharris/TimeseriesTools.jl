module SpectralFittingExt
# Use this by `] registry add https://github.com/astro-group-bristol/AstroRegistry`
# then `] add SpectralFitting`
using TimeseriesTools
using SpectralFitting
using Statistics
import SpectralFitting: InjectiveData

function InjectiveData(x::ToolsArray{T, 1}, args...; kwargs...) where {T}
    InjectiveData((Vector ∘ parent ∘ only ∘ lookup)(x), parent(x); kwargs...)
end

function InjectiveData(x::ToolsArray, args...; dims = 1, kwargs...)
    length(dims) > 1 && throw(ArgumentError("Can only fit spectra along one dimension"))
    dimidx = dimnum(x, dims)
    negdims = ntuple(i -> i < dimidx ? i : i + 1, ndims(x) - 1) # Type stable

    μ = mean(x, dims = negdims) |> parent
    μ = dropdims(μ, dims = negdims)

    σ = std(x, dims = negdims) |> parent
    σ = dropdims(σ, dims = negdims)

    f = parent(lookup(x, only(dims)))
    InjectiveData(f, μ, args...; codomain_variance = σ, kwargs...)
end

end
