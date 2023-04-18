using FFTW
using Statistics

export  FreqDim, Freq,
        AbstractSpectrum, RegularSpectrum, MultivariateSpectrum,
        spectrum, energyspectrum, powerspectrum,
        _energyspectrum, _powerspectrum,
        FreqIndex, RegularFreqIndex,
        colorednoise


abstract type FreqDim{T} <: DimensionalData.IndependentDim{T} end

"""
    Freq

A DimensionalData.jl dimension representing the frequency domain.
"""
DimensionalData.@dim Freq FreqDim "Freq"

"""
    FreqIndex

A type alias for a tuple of dimensions, where the first dimension is of type `FreqDim`.
"""
FreqIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:FreqDim}
UnitfulFreqIndex = UnitfulTimeIndex

"""
    AbstractSpectrum{T, N, B}

A type alias for an `AbstractDimArray` in which the first dimension is [`Freq`](@ref)uency.
"""
AbstractSpectrum = AbstractDimArray{T, N, <:FreqIndex, B} where {T, N, B}

"""
    RegularFreqIndex

A type alias for a tuple of dimensions, where the first dimension is a regularly sampled [`Freq`](@ref)uency.
"""
RegularFreqIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:FreqDim{<:RegularIndex}}

"""
    RegularSpectrum{T, N, B}

A type alias for a spectrum with a regularly sampled frequency index.
"""
RegularSpectrum = AbstractDimArray{T, N, <:RegularFreqIndex, B} where {T, N, B}

"""
    MultivariateSpectrum{T} = AbstractSpectrum{T, 2} where T

A type alias for a multivariate spectrum.
"""
MultivariateSpectrum = AbstractSpectrum{T, 2} where T

"""
    _energyspectrum(x::RegularTimeSeries, f_min=samplingrate(x)/min(length(x)÷4, 1000))

Computes the energy spectrum of a regularly sampled time series `x` with an optional minimum frequency `f_min`.
"""
function _energyspectrum(x::RegularTimeSeries, f_min=samplingrate(x)/min(length(x)÷4, 1000))
    fs = samplingrate(x)
    n = length(x)
    validfreqs = rfftfreq(n, fs)
    if f_min == 0
        f_min = validfreqs[2]
        nfft = floor(Int, n/2)*2
    else
        nfft = ceil(Int, ustrip(fs)/ustrip(f_min))
    end
    if ustrip(f_min) < ustrip(validfreqs[2])
        error("Cannot resolve an `f_min` of $f_min")
    end

    isodd(nfft) && (nfft += 1)
    window = nfft÷2
    overlap = window÷2
    hann_window = 0.5 .- 0.5 .* cos.(2 * π * (0:(nfft - 1)) / (nfft - 1))
    n_segments = floor(Int, (n - nfft) / (nfft - overlap) + 1)

    # Get the type of the spectrum
    u = unit(eltype(x)) * unit(eltype(dims(x, Ti)))
    S̄ = zeros(nfft ÷ 2 + 1, n_segments)*u^2
    for i in 1:n_segments
        start_idx = (i - 1) * (nfft - overlap) + 1
        end_idx = start_idx + nfft - 1
        segment = x[start_idx:end_idx] .* hann_window

        y = rfft(segment) / nfft
        y |> eltype |> unit == NoUnits && (y = y * u)
        S̄[:, i] .= (abs.(y).^2) / sum(hann_window.^2)
    end


    # Calculate the frequencies
    freqs = range(convertconst(0, fs), stop=fs/2, length=size(S̄, 1))
    df = step(freqs)

    # Normalize the mean energy spectrum to obey Parseval's theorem
    meanS̄ = mean(S̄, dims=2)
    S̄ = S̄ ./ ustrip((sum(meanS̄) - 0.5.*meanS̄[1]) .* df) # Subtract the zero frequency component twice, so that it doesn't bias when we divide by a half
    S̄ = 0.5 * S̄ .* ustrip(sum(x.^2) ./ fs) # Normalized to have total energy equal to energy of signal. Ala parseval. 0.5 because we only return the positive half of the spectrum.

    return DimArray(S̄, (Freq(freqs), Dim{:window}(1:n_segments)))
end
"""
    _energyspectrum(x::RegularTimeSeries, f_min=0)

Computes the energy spectrum of a time series using the fast Fourier transform.

If `f_min > 0`, the energy spectrum is calculated for windows of the time series determined by `f_min`,  the minimum frequency that will be resolved in the spectrum.
If `f_min > 0`, the second dimension of the output will correspond to the windows. For an averaged periodogram, see [`energyspectrum`](@ref).

If the input time series is a [`UnitfulTimeSeries`](@ref), the frequency will also have units.
Moreover if the elements of `x` are unitful, so are the elements of the spectrum.

# Examples
```jldoctest
julia> using TimeseriesTools
julia> t = range(0.0, stop=1.0, length=1000);
julia> x = sin.(2 * π * 5 * t);
julia> ts = RegularTimeSeries(t, x);
julia> S = _energyspectrum(ts);
julia> S isa MultivariateSpectrum
"""
_energyspectrum(x::typeintersect(MultivariateTS, RegularTS), args...; kwargs...) = cat([_energyspectrum(_x, args...; kwargs...) for _x in eachslice(x, dims=2)]..., dims=dims(x, 2))

"""
    energyspectrum(x::RegularTimeSeries, f_min=0; kwargs...)

Computes the average energy spectrum of a regularly sampled time series `x`.
`f_min` determines the minimum frequency that will be resolved in the spectrum.
See [`_energyspectrum`](@ref).
"""
energyspectrum(x::RegularTimeSeries, args...; kwargs...) = dropdims(mean(_energyspectrum(x, args...; kwargs...), dims=Dim{:window}); dims=Dim{:window})

"""
    _powerspectrum(x::AbstractTimeSeries, f_min=samplingrate(x)/min(length(x)÷4, 1000); kwargs...)

Computes the power spectrum of a time series `x` in Welch periodogram windows.
Note that the `_powerspectrum` is simply the [`_energyspectrum`](@ref) divided by the duration of each window.
See [`_energyspectrum`](@ref).
"""
function _powerspectrum(x::AbstractTimeSeries, args...; kwargs...)
    S̄ = _energyspectrum(x, args...; kwargs...)
    return S̄ ./ duration(x)
end

"""
    powerspectrum(x::AbstractTimeSeries, f_minsamplingrate(x)/min(length(x)÷4, 1000); kwargs...)

Computes the average power spectrum of a time series `x` using the Welch periodogram method.
"""
powerspectrum(x::AbstractTimeSeries, args...; kwargs...) = dropdims(mean(_powerspectrum(x, args...; kwargs...), dims=Dim{:window}); dims=Dim{:window})

"""
    colorednoise(ts::AbstractRange; α=2.0)

Generate a colored-noise time series with a specified power-law exponent `α` on the given times `ts`.

# Arguments
- `ts`: An `AbstractRange` representing the time range of the generated noise.
- `α`: The power-law exponent of the colored noise, which will have a spectrum given by 1/f^α. Defaults to 2.0.

# Returns
- A [`TimeSeries`](@ref) containing the generated colored noise.

# Example

```jldoctest
julia> using TimeseriesTools
julia> pink_noise = colorednoise(1:0.01:10; α=1.0)
julia> pink_noise isa RegularTimeSeries
````
"""
function colorednoise(ts::AbstractRange, args...; α=2.0)
    f = rfftfreq(length(ts), step(ts))
    x̂ = sqrt.(1.0./f.^α).*exp.(2π.*rand(length(f))*im)
    x̂[1] = 0
    x = irfft(x̂, 2*length(f)-2)
    dt = length(ts)*step(f)
    t = range(dt, length(x)*dt, length=length(x))
    @assert all(t .== ts)
    TimeSeries(t, x, args...)
end
