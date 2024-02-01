module DierckxExt
using Dierckx
using DimensionalData
using TimeseriesTools
import TimeseriesTools: interpolate, upsample

function interpolate(X::DimensionalData.AbstractDimMatrix; kwargs...)
    xy = [Float64.(x) for x in lookup(X)]
    itp = Spline2D(xy..., X.data; kwargs...)
    return itp
end

function (itp::Spline2D)(x::DimensionalData.Dimension,
                         y::DimensionalData.Dimension; kwargs...)
    DimArray(evalgrid(itp, lookup(x), lookup(y)), (x, y); kwargs...)
end

function upsample(x::DimensionalData.AbstractDimMatrix, factor; kwargs...)
    xy = upsample.(dims(x), factor)
    itp = interpolate(x; kwargs...)
    itp(xy...; name = DimensionalData.name(x),
        metadata = DimensionalData.metadata(x),
        refdims = DimensionalData.refdims(x))
end

end # module
