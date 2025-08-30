module DataInterpolationsExt
using DataInterpolations
import DataInterpolations: AbstractInterpolation
using DimensionalData
import DimensionalData: AbstractDimArray, Dimension
using TimeseriesTools
using Unitful
import TimeseriesTools: interpolate, upsample

function interpolate(x::DimensionalData.AbstractDimVector,
                     interp::Type{<:AbstractInterpolation} = AkimaInterpolation,
                     args...; extrapolation = ExtrapolationType.Extension, kwargs...)
    DimensionalData.isforward(x) ||
        throw(ArgumentError("`interpolate` only supports forward-ordered data"))
    interp(parent(x) .|> ustrip, x |> lookup |> only .|> ustrip, args...; extrapolation,
           kwargs...)
end
function (itp::AbstractInterpolation)(x::DimensionalData.Dimension; kwargs...)
    DimensionalData.dimconstructor((x,))((x |> lookup .|> ustrip |> itp),
                                         (x,);
                                         kwargs...)
end

function upsample(x::DimensionalData.AbstractDimArray, factor::Number, args...;
                  dims = 1, kwargs...)
    u = unit(eltype(x))
    y = x
    for dim in dims
        d = upsample(DimensionalData.dims(y, dim), factor)
        negdims = setdiff(1:ndims(y), dimnum(y, dim)) |> Tuple
        y = map(eachslice(y; dims = negdims)) do y
            itp = interpolate(y, args...; kwargs...)
            itp(d)
        end |> stack
        y = permutedims(y, dimnum.([y], DimensionalData.dims(x)))
    end
    return y * u
end

end # module
