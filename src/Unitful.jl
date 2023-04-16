using Unitful
export UnitfulIndex, UnitfulTimeSeries
using FFTW

# Unitful._promote_unit(::S, ::T) where {S<:Unitful.FreeUnits{(), NoDims, nothing}, T<:Unitful.TimeUnits} = u"s"
TimeseriesTools.convertconst(c::Number, u::Unitful.Quantity) = (c)unit(u)
UnitfulIndex = Union{AbstractArray{<:Unitful.Time}, AbstractRange{<:Unitful.Time}, Tuple{<:Unitful.Time}}
UnitfulTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:DimensionalData.Dimension{<:UnitfulIndex}}
UnitfulTimeSeries = AbstractDimArray{T, N, <:UnitfulTimeIndex, B} where {T, N, B}

TimeSeries(t, x, unit::Unitful.Units) = TimeSeries((t)unit, x)

function FFTW.rfft(x::AbstractVector{<:Quantity})
    # Assume this is a discrete fourier transform, to time indices/units
    s = x |> eltype |> unit
    x_re = reinterpret(Float64, collect(x))
    f = rfft(x_re)
    return (f)s
end

function FFTW.rfft(x::UnitfulTimeSeries{<:Quantity})
    # In this case the rfft is treated as an approximation to the continuosu Fourier transform
    t = samplingperiod(x)
    a = x |> eltype |> unit
    x_re = reinterpret(Float64, collect(x))
    ℱ = rfft(x_re)
    ℱ = (ℱ)*(a*t) # CTFT has units of amplitude*time. Normalise the DFT to have bins the width of the sampling period.
end
