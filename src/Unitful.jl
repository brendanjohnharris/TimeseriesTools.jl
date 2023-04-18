using Unitful
export UnitfulIndex, UnitfulTimeSeries
using FFTW

# Unitful._promote_unit(::S, ::T) where {S<:Unitful.FreeUnits{(), NoDims, nothing}, T<:Unitful.TimeUnits} = u"s"
"""
    TimeseriesTools.convertconst(c::Number, u::Unitful.Quantity)

Converts a constant `c` to have the same units as `u`.

## Examples
```jldoctest
julia> using Unitful;
julia> c = 5;
julia> u = 3u"s";
julia> converted_c = TimeseriesTools.convertconst(c, u);
julia> typeof(converted_c) == typeof(u)
```
"""
TimeseriesTools.convertconst(c::Number, u::Unitful.Quantity) = (c)unit(u)

"""
    UnitfulIndex

A type alias for a union of `AbstractArray`, `AbstractRange`, and `Tuple` types with `Unitful.Time` elements.
"""
UnitfulIndex = Union{AbstractArray{<:Unitful.Time}, AbstractRange{<:Unitful.Time}, Tuple{<:Unitful.Time}}

"""
    UnitfulTimeIndex

A type alias for a tuple of dimensions, where the first dimension is of type `DimensionalData.Dimension{<:UnitfulIndex}`.
"""
UnitfulTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:DimensionalData.Dimension{<:UnitfulIndex}}

"""
    UnitfulTimeSeries{T, N, B}

A type alias for an `AbstractDimArray` with a `UnitfulTimeIndex`.

## Examples
```jldoctest
julia> using Unitful;
julia> t = (1:100)u"s";
julia> x = rand(100);
julia> uts = TimeSeries(t, x);
julia> uts isa UnitfulTimeSeries
```
"""
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
