module CUDAExt
using CUDA
using TimeseriesTools

CUDA.CuArray(x::AbstractTimeseries) = set(x, CuArray(x.data))

function Base.show(
        io::IO, mime,
        X::ToolsArray{T, N, Tp, F, C} where {T, N, Tp, F, C <: CuArray}
    )
    return Base.show(X.data)
end

end # module
