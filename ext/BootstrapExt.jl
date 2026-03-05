module BootstrapExt
using Bootstrap
using Statistics
import TimeseriesTools: bootstrapaverage, bootstrapmedian, bootstrapmean, nansafe

function bootstrapaverage(average, x::AbstractVector; confint = 0.95,
                          N = 10000)#::Tuple{T, Tuple{T, T}} where {T}
    sum(!isnan, x) < 5 && return (NaN, (NaN, NaN))

    # * Estimate a sampling distribution of the average
    x = filter(!isnan, x)
    b = Bootstrap.bootstrap(nansafe(average), x, Bootstrap.BalancedSampling(N))
    μ, σ... = only(Bootstrap.confint(b, Bootstrap.BCaConfInt(confint)))
    return (; average = μ, confint = (; lower = σ[1], upper = σ[2]))
end

# function bootstrapaverage(average, X::AbstractArray; dims = 1, kwargs...)
#     if length(dims) > 1
#         error("Only one dimension can be specified")
#     end
#     ds = [i == dims ? 1 : Colon() for i in 1:ndims(X)]
#     μ = similar(X[ds...])
#     σl = similar(μ)
#     σh = similar(μ)
#     negdims = filter(!=(dims), 1:ndims(X)) |> Tuple
#     Threads.@threads for (i, x) in collect(enumerate(eachslice(X; dims = negdims)))
#         μ[i], (σl[i], σh[i]) = bootstrapaverage(average, x; kwargs...)
#     end
#     return μ, (σl, σh)
# end
# function bootstrapaverage(average, X::AbstractToolsArray; dims = 1, kwargs...)
#     if (dims isa Tuple || dims isa AbstractVector) && length(dims) > 1
#         error("Only one dimension can be specified")
#     end
#     dims = dimnum(X, dims)
#     ds = [i == dims ? 1 : Colon() for i in 1:ndims(X)]
#     μ = similar(X[ds...])
#     σl = similar(μ)
#     σh = similar(μ)
#     negdims = filter(!=(dims), 1:ndims(X)) |> Tuple
#     Threads.@threads for (i, x) in collect(enumerate(eachslice(X; dims = negdims)))
#         μ[i], (σl[i], σh[i]) = bootstrapaverage(average, parent(x); kwargs...)
#     end
#     return μ, (σl, σh)
# end

function bootstrapaverage(average, X::AbstractArray; dims = 1, kwargs...)
    ds = [i ∈ dims ? 1 : Colon() for i in 1:ndims(X)]
    μ = similar(X[ds...])
    σl = similar(μ)
    σh = similar(μ)
    negdims = filter(!∈(dims), 1:ndims(X)) |> Tuple
    Threads.@threads for (i, x) in collect(enumerate(eachslice(X; dims = negdims)))
        μ[i], (σl[i], σh[i]) = bootstrapaverage(average, x[:]; kwargs...)
    end
    return (; average = μ, confint = (; lower = σl, upper = σh))
end

bootstrapmedian(args...; kwargs...) = bootstrapaverage(median, args...; kwargs...)
bootstrapmean(args...; kwargs...) = bootstrapaverage(mean, args...; kwargs...)

end
