module NaturalNeighboursExt

using NaturalNeighbours
import NaturalNeighbours: interpolate, NaturalNeighboursInterpolant

using DimensionalData
using Normalization
using TimeseriesTools
import TimeseriesTools: upsample

function interpolate(X::DimensionalData.AbstractDimMatrix; kwargs...)
    xy = [Float64.(x) for x in lookup(X)]
    N = fit.([MinMax], xy)

    points = Iterators.product((xy .|> N)...) |> collect |> vec
    itp = interpolate(Float64.(first.(points)), Float64.(last.(points)), X.data[:];
                      kwargs...)
    return itp, N
end

function (itp::NaturalNeighboursInterpolant)(x::DimensionalData.Dimension,
                                             y::DimensionalData.Dimension, N; kwargs...)
    points = Iterators.product((collect.((x, y)) .|> N)...) |> collect |> vec
    X = itp(Float64.(first.(points)), Float64.(last.(points)); kwargs...)
    X = reshape(X, (length(x), length(y)))
    ToolsArray(X, (x, y))
end

# function upsample(x::DimensionalData.AbstractDimMatrix, factor; derivatives = true,
#                   kwargs...)
#     xy = upsample.(dims(x), factor)
#     itp, N = interpolate(x; derivatives)
#     itp(xy..., N; method = Sibson(1), parallel = true, kwargs...)
# end

end # module
