module DierckxExt
using Dierckx
using DimensionalData
using TimeseriesTools
using Unitful
import TimeseriesTools: interpolate, upsample

function interpolate(X::DimensionalData.AbstractDimMatrix; kwargs...)
    xy = [Float64.(x) for x in lookup(X)]
    itp = Spline2D(xy..., X.data; kwargs...)
    return itp
end
function interpolate(X::DimensionalData.AbstractDimVector; kwargs...)
    x = Float64.(lookup(X)[1])
    itp = Spline1D(x, X.data; kwargs...)
    return itp
end

function (itp::Spline2D)(x::DimensionalData.Dimension,
                         y::DimensionalData.Dimension; kwargs...)
    DimArray(evalgrid(itp, lookup(x), lookup(y)), (x, y); kwargs...)
end
function (itp::Spline1D)(x::DimensionalData.Dimension; kwargs...)
    DimArray(evaluate(itp, lookup(x)), (x,); kwargs...)
end

function upsample(x::DimensionalData.AbstractDimMatrix, factor; kwargs...)
    xy = upsample.(dims(x), factor)
    itp = interpolate(x; kwargs...)
    itp(xy...; name = DimensionalData.name(x),
        metadata = DimensionalData.metadata(x),
        refdims = DimensionalData.refdims(x))
end
function upsample(x::DimensionalData.AbstractDimVector, factor::Integer; kwargs...)
    d = upsample(dims(x)[1], factor)
    itp = interpolate(x; kwargs...)
    itp(d; name = DimensionalData.name(x),
        metadata = DimensionalData.metadata(x),
        refdims = DimensionalData.refdims(x))
end
function upsample(x::DimensionalData.AbstractDimMatrix, factor::Number,
                  dim; kwargs...)
    d = upsample(DimensionalData.dims(x, dim), factor)
    adims = setdiff(1:ndims(x), dimnum(x, dim)) |> only
    if only(adims) == 1
        x = x'
    end
    u = unit(eltype(x))
    z = DimArray(zeros(length(d), size(x, 2)), (d, dims(x, 2)))
    itp = Spline1D.([ustrip(lookup(x, 1))], eachslice(ustrip.(x.data), dims = 2);
                    kwargs...)
    Threads.@threads for (i, I) in collect(enumerate(itp))
        z[:, i] = I(ustrip.(d)) .* u
    end
    if only(adims) == 1
        z = z'
    end
    return z
end

end # module
