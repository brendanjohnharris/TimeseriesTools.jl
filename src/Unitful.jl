using Unitful
using FFTW
import Unitful.unit

export dimunit, timeunit, frequnit, unit, UnitfulIndex, UnitfulTimeSeries, UnitfulSpectrum

# Unitful._promote_unit(::S, ::T) where {S<:Unitful.FreeUnits{(), NoDims, nothing}, T<:Unitful.TimeUnits} = u"s"
"""
    TimeseriesTools.convertconst(c::Number, u::Unitful.Quantity)

Converts a constant `c` to have the same units as `u`.

## Examples
```@example 1
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
UnitfulIndex = UnitfulTIndex = Union{AbstractArray{<:Unitful.Time}, AbstractRange{<:Unitful.Time}, Tuple{<:Unitful.Time}}


"""
    UnitfulTimeIndex

A type alias for a tuple of dimensions, where the first dimension is of type `DimensionalData.Dimension{<:UnitfulIndex}`.
"""
UnitfulTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:DimensionalData.Dimension{<:UnitfulIndex}}

"""
    UnitfulTimeSeries{T, N, B}

A type alias for an `AbstractDimArray` with a [`UnitfulTimeIndex`](@ref).

## Examples
```@example 1
julia> using Unitful;
julia> t = (1:100)u"s";
julia> x = rand(100);
julia> uts = TimeSeries(t, x);
julia> uts isa UnitfulTimeSeries
```
"""
UnitfulTimeSeries = AbstractDimArray{T, N, <:UnitfulTimeIndex, B} where {T, N, B}

UnitfulFIndex = Union{AbstractArray{<:Unitful.Frequency}, AbstractRange{<:Unitful.Frequency}, Tuple{<:Unitful.Frequency}}
UnitfulFreqIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:DimensionalData.Dimension{<:UnitfulFIndex}}
UnitfulSpectrum = AbstractDimArray{T, N, <:UnitfulFreqIndex, B} where {T, N, B}

"""
    TimeSeries(t, x, timeunit::Unitful.Units)

Constructs a univariate time series with time `t`, data `x`, time units specified by `timeunit`, and element units specified by `elementunit`.
Note that you can add units to the elements of a time series with, for example, `x*u"V"`.

## Examples
```@example 1
julia> using Unitful;
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x, u"ms")*u"V";
julia> ts isa Union{UnivariateTimeSeries, RegularTimeSeries, UnitfulTimeSeries}
```
"""
TimeSeries(t, x, unit::Unitful.Units) = TimeSeries((t)unit, x)

"""
    dimunit(x::UnitfulTimeSeries, dim)

Returns the unit associated with the specified dimension `dim` of a [`UnitfulTimeSeries`](@ref).

## Examples
```@example 1
julia> using Unitful;
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x, u"ms");
julia> TimeseriesTools.dimunit(ts, Ti) == u"ms"
```
"""
dimunit(x::DimArray, dim) = dims(x,dim) |> eltype |> unit

"""
    timeunit(x::UnitfulTimeSeries)

Returns the time units associated with a [`UnitfulTimeSeries`].

## Examples
```@example 1
julia> using Unitful;
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x, u"ms");
julia> timeunit(ts) == u"ms"
```
"""
timeunit(x::AbstractTimeSeries) = dimunit(x, Ti)

"""
    frequnit(x::UnitfulSpectrum)

Returns the frequency units associated with a [`UnitfulSpectrum`](@ref).

## Examples
```@example 1
julia> using Unitful;
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x, u"ms");
julia> sp = fft(ts);  # assuming fft returns a UnitfulSpectrum
julia> frequnits(sp) == u"Hz"
```
"""
frequnit(x::AbstractSpectrum) = dimunit(x, Freq)

"""
    unit(x::AbstractArray)

Returns the units associated with the elements of an [`UnitfulTimeSeries``](@ref) or [`UnitfulSpectrum`](@ref).

## Examples
```@example 1
julia> using Unitful;
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x, u"ms")*u"V";
julia> unit(ts) == u"V"
```
"""
unit(x::Union{<:AbstractTimeSeries, AbstractSpectrum}) = x |> eltype |> unit

function FFTW.rfft(x::AbstractVector{<:Quantity}) # 🐶
    # Assume this is a discrete Fourier transform, to time indices/units
    s = x |> eltype |> unit
    x_re = reinterpret(Float64, collect(x))
    f = rfft(x_re)
    return (f)s
end

function FFTW.rfft(x::UnitfulTimeSeries{<:Quantity}) # 🐕
    # In this case the rfft is treated as an approximation to the continuosu Fourier transform
    t = samplingperiod(x)
    a = x |> eltype |> unit
    x_re = reinterpret(Float64, collect(x))
    ℱ = rfft(x_re)
    ℱ = (ℱ)*(a*t) # CTFT has units of amplitude*time. Normalise the DFT to have bins the width of the sampling period.
end
