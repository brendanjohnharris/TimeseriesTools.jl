using Unitful
using FFTW
import Unitful.unit
import LinearAlgebra.normalize
import Normalization.denormalize

function FFTW.rfft(x::AbstractVector{<:Quantity}) # 🐶
    # Assume this is a discrete Fourier transform, to time indices/units
    s = x |> eltype |> unit
    x_re = reinterpret(Float64, collect(x))
    f = rfft(x_re)
    return (f)s
end

function FFTW.rfft(x::UnitfulTimeseries{<:Quantity}) # 🐕
    # In this case the rfft is treated as an approximation to the continuosu Fourier transform
    t = samplingperiod(x)
    a = x |> eltype |> unit
    x_re = reinterpret(Float64, collect(x))
    ℱ = rfft(x_re)
    ℱ = (ℱ) * (a * t) # CTFT has units of amplitude*time. Normalise the DFT to have bins the width of the sampling period.
end

# Extend Normalization.jl to unitful DimArrays
# function normalize(X::AbstractDimArray{<:Quantity}, T::NormUnion; kwargs...)
#     DimensionalData.modify(x -> normalize(x, T; kwargs...), X)
# end
# function denormalize(Y::AbstractDimArray{<:Quantity}, T::AbstractNormalization{<:Quantity};
#                      kwargs...)
#     error("Denormalization of unitful arrays currently not supported")
# end
