using FFTW
using Statistics

import TimeseriesTools.Operators.𝒯

export FrequencyDim, Freq, freqs,
       AbstractSpectrum, RegularSpectrum, UnivariateSpectrum, MultivariateSpectrum,
       spectrum, energyspectrum, powerspectrum,
       _energyspectrum, _powerspectrum,
       FreqIndex, RegularFreqIndex,
       colorednoise,
       AbstractBispectrum, Bispectrum, BifreqIndex, RegularBifreqIndex, bispectrum

abstract type FrequencyDim{T} <: DimensionalData.IndependentDim{T} end

"""
    Freq

A DimensionalData.jl dimension representing the frequency domain.
"""
DimensionalData.@dim Freq FrequencyDim "Freq"

"""
    FreqIndex

A type alias for a tuple of dimensions, where the first dimension is of type `FrequencyDim`.
"""
const FreqIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A <: Freq}

"""
    AbstractSpectrum{T, N, B}

A type alias for an `AbstractDimArray` in which the first dimension is [`Freq`](@ref)uency.
"""
const AbstractSpectrum = AbstractDimArray{T, N, <:FreqIndex, B} where {T, N, B}
freqs(x::AbstractSpectrum) = dims(x, Freq).val.data

"""
    RegularFreqIndex

A type alias for a tuple of dimensions, where the first dimension is a regularly sampled [`Freq`](@ref)uency.
"""
const RegularFreqIndex = Tuple{A,
                               Vararg{DimensionalData.Dimension}} where {A <:
                                                                         FrequencyDim{<:RegularIndex}}

"""
    RegularSpectrum{T, N, B}

A type alias for a spectrum with a regularly sampled frequency index.
"""
const RegularSpectrum = AbstractDimArray{T, N, <:RegularFreqIndex, B} where {T, N, B}

"""
    UnivariateSpectrum{T} = AbstractSpectrum{T, 1} where T

A type alias for a univariate spectrum.
"""
const UnivariateSpectrum = AbstractSpectrum{T, 1} where {T}
"""
    MultivariateSpectrum{T} = AbstractSpectrum{T, 2} where T

A type alias for a multivariate spectrum.
"""
const MultivariateSpectrum = AbstractSpectrum{T, 2} where {T}

"""
    Spectrum(f, x)

Constructs a univariate spectrum with frequencies `f` and data `x`.
"""
Spectrum(f, x; kwargs...) = DimArray(x, (Freq(f),); kwargs...)

"""
    Spectrum(f, v, x)

Constructs a multivariate spectrum with frequencies `f`, variables `v`, and data `x`.
"""
Spectrum(f, v, x; kwargs...) = DimArray(x, (Freq(f), Var(v)); kwargs...)
function Spectrum(f, v::DimensionalData.Dimension, x; kwargs...)
    DimArray(x, (Freq(f), v); kwargs...)
end

function _energyspectrum(x::AbstractVector, fs::Number,
                         f_min::Number = fs / min(length(x) ÷ 4, 1000); padding = 0,
                         kwargs...)
    n = length(x)
    validfreqs = rfftfreq(n, fs)
    if f_min == 0
        f_min = validfreqs[2]
        nfft = floor(Int, n / 2) * 2
    else
        nfft = ceil(Int, ustripall(fs) / ustripall(f_min))
    end
    if ustripall(f_min) < ustripall(validfreqs[2])
        error("Cannot resolve an `f_min` of $f_min")
    end

    isodd(nfft) && (nfft += 1)
    isodd(padding) && (padding += 1)
    nfft = nfft - padding
    overlap = nfft ÷ 2
    hann_window = 0.5 .- 0.5 .* cos.(2 * π * (0:(nfft - 1)) / (nfft - 1))
    A = sum(hann_window .^ 2)
    (nfft - overlap ≤ 0) && error("FFT padding is too high")
    n_segments = floor(Int, (n - nfft) / (nfft - overlap) + 1)

    # Get the type of the spectrum
    u = unit(eltype(x)) / unit(fs)
    S̄ = zeros((nfft + padding) ÷ 2 + 1, n_segments) * u^2
    @debug "Calculating spectrum for $n_segments segments of length $(nfft + padding)"
    for i in 1:n_segments
        start_idx = (i - 1) * (nfft - overlap) + 1
        end_idx = start_idx + nfft - 1
        segment = x[start_idx:end_idx] .* hann_window
        if padding > 0
            dt = samplingperiod(segment)
            padts = range(start = minimum(times(segment)) + dt, step = dt,
                          length = padding + length(segment))
            segment = TimeSeries(padts,
                                 [segment.data; zeros(padding) * unit(eltype(segment))])
        end

        y = rfft(segment) / (nfft + padding)
        y |> eltype |> unit == NoUnits && (y = y * u)
        S̄[:, i] .= (abs.(y) .^ 2) / A
    end

    # Calculate the frequencies
    freqs = range(convertconst(0, fs), stop = fs / 2, length = size(S̄, 1))
    df = step(freqs)

    # Normalize the mean energy spectrum to obey Parseval's theorem
    meanS̄ = mean(S̄, dims = 2)
    S̄ = S̄ ./ ustripall((sum(meanS̄) - 0.5 .* meanS̄[1]) .* df) # Subtract the zero frequency component twice, so that it doesn't bias when we divide by a half
    S̄ = 0.5 * S̄ .* ustripall(sum(x .^ 2) ./ fs) # Normalized to have total energy equal to energy of signal. Ala parseval. 0.5 because we only return the positive half of the spectrum.
    Spectrum(freqs, Dim{:window}(1:n_segments), S̄; kwargs...)
end

"""
    _energyspectrum(x::RegularTimeSeries, f_min=samplingrate(x)/min(length(x)÷4, 1000))

Computes the energy spectrum of a regularly sampled time series `x` with an optional minimum frequency `f_min`.
"""
function _energyspectrum(x::typeintersect(RegularTimeSeries, UnivariateTimeSeries);
                         f_min::Number = samplingrate(x) / min(length(x) ÷ 4, 1000),
                         kwargs...)
    return _energyspectrum(x, samplingrate(x), f_min; kwargs...)
end

"""
    _energyspectrum(x::RegularTimeSeries, f_min=0)

Computes the energy spectrum of a time series using the fast Fourier transform.

If `f_min > 0`, the energy spectrum is calculated for windows of the time series determined by `f_min`,  the minimum frequency that will be resolved in the spectrum.
If `f_min > 0`, the second dimension of the output will correspond to the windows. For an averaged periodogram, see [`energyspectrum`](@ref).

If the input time series is a [`UnitfulTimeSeries`](@ref), the frequency will also have units.
Moreover if the elements of `x` are unitful, so are the elements of the spectrum.

# Examples
```@example 1
julia> using TimeseriesTools
julia> t = range(0.0, stop=1.0, length=1000);
julia> x = sin.(2 * π * 5 * t);
julia> ts = RegularTimeSeries(t, x);
julia> S = _energyspectrum(ts);
julia> S isa MultivariateSpectrum
```
"""
function _energyspectrum(x::MultivariateTS; kwargs...)
    cat([_energyspectrum(_x; kwargs...)
         for _x in eachslice(x, dims = 2)]..., dims = dims(x, 2))
end

"""
    energyspectrum(x::RegularTimeSeries, f_min=0; kwargs...)

Computes the average energy spectrum of a regularly sampled time series `x`.
`f_min` determines the minimum frequency that will be resolved in the spectrum.
See [`_energyspectrum`](@ref).
"""
function energyspectrum(x; kwargs...)
    dropdims(mean(_energyspectrum(x; kwargs...), dims = Dim{:window});
             dims = Dim{:window})
end

"""
    _powerspectrum(x::AbstractTimeSeries, f_min=samplingrate(x)/min(length(x)÷4, 1000); kwargs...)

Computes the power spectrum of a time series `x` in Welch periodogram windows.
Note that the `_powerspectrum` is simply the [`_energyspectrum`](@ref) divided by the duration of each window.
See [`_energyspectrum`](@ref).
"""
function _powerspectrum(x::AbstractTimeSeries; kwargs...)
    S̄ = _energyspectrum(x; kwargs...)
    return S̄ ./ duration(x)
end

"""
    powerspectrum(x::AbstractTimeSeries, f_min=samplingrate(x)/min(length(x)÷4, 1000); kwargs...)

Computes the average power spectrum of a time series `x` using the Welch periodogram method.
"""
function powerspectrum(x::AbstractTimeSeries; kwargs...)
    dropdims(mean(_powerspectrum(x; kwargs...), dims = Dim{:window});
             dims = Dim{:window})
end

spectrum = powerspectrum

"""
    colorednoise(ts::AbstractRange; α=2.0)

Generate a colored-noise time series with a specified power-law exponent `α` on the given times `ts`.

# Arguments
- `ts`: An `AbstractRange` representing the time range of the generated noise.
- `α`: The power-law exponent of the colored noise, which will have a spectrum given by 1/f^α. Defaults to 2.0.

# Returns
- A [`TimeSeries`](@ref) containing the generated colored noise.

# Example

```@example 1
julia> using TimeseriesTools
julia> pink_noise = colorednoise(1:0.01:10; α=1.0)
julia> pink_noise isa RegularTimeSeries
```
"""
function colorednoise(ts::AbstractRange, args...; α = 2.0)
    f = rfftfreq(length(ts), step(ts))
    x̂ = sqrt.(1.0 ./ f .^ α) .* exp.(2π .* rand(length(f)) * im)
    x̂[1] = 0
    x = irfft(x̂, length(ts))
    dt = length(ts) * step(f)
    t = range(0, (length(x) - 1) * dt, length = length(x))
    @assert all(t .+ first(ts) .≈ ts)
    TimeSeries(ts, x, args...)
end

function spikefft(t::AbstractVector, ::Val{:schild})
    # AN EFFICIENT METHOD FOR THE FOURIER TRANSFORM OF A NEURONAL  SPIKE  TRAIN
    # Schild 1982
    @debug "Calculating spike FFT using :schild method"
    t .-= minimum(t)
    T = maximum(t)
    W(f) = (sum(cos.(2π * f .* t))^2 + sum(sin.(2π * f .* t))^2) / T
end

spikefft(fs, t::AbstractVector, method) = spikefft(t, method).(fs)

function spikefft(fs, t::SpikeTrain, method = :schild)
    isempty(findfirst(t)) && error("Spike train contains no spikes")
    F = spikefft(times(t[t]), Val(method))
    return Spectrum(fs, F.(fs))
end

function _energyspectrum(x::SpikeTrain{T, 1} where {T}, frange::AbstractRange;
                         method = stoic(; σ = 0.005), kwargs...)
    t = times(x[x])
    n = length(t)
    df = step(frange)
    isodd(length(frange)) && (frange = frange[1:(end - 1)])
    nfft = length(frange)

    u = unit(eltype(x)) * unit(eltype(t))

    if method isa Function
        # Find the autocovariance function and take the Fourier transform
        τs = range(start = 0, stop = 1 / step(frange) / 2, length = nfft)
        τs = [-reverse(τs[2:end]); τs]
        ρ = [method(x, 𝒯(τ)(x)) for τ in τs]
        S̄ = abs.(rfft(ρ))
        @assert length(S̄) == nfft
    else
        S̄ = abs.(spikefft(frange, x[x], method)) .^ 2
    end

    # Normalize the energy spectrum to obey Parseval's theorem
    S̄ = S̄ ./ ustripall((2 * sum(S̄) - S̄[1]) .* df) # Subtract the zero frequency component a bit, so it doesn't bias when we divide by half
    # display(2 * sum(S̄[2:end]) + S̄[1])
    S̄ = S̄ .* ustripall(length(t)) # Normalized to have total energy equal to energy of signal. Ala parseval.
    Spectrum(frange, Dim{:window}([1]), Matrix(S̄')'; kwargs...)
end

function _energyspectrum(x::SpikeTrain{T, 1} where {T}, frange::Tuple; kwargs...)
    _energyspectrum(x, 0:first(frange):last(frange); kwargs...)
end

"""
    BifreqIndex

A type alias for a tuple of dimensions, where the first and second dimensions are of type `FrequencyDim`.
"""
const BifreqIndex = Tuple{A, A, Vararg{DimensionalData.Dimension}} where {A <: Freq}

"""
    RegularBifreqIndex

A type alias for a tuple of dimensions, where the first and second dimensions are a regularly sampled [`Freq`](@ref)uencies.
"""
const RegularBifreqIndex = Tuple{A, A,
                                 Vararg{DimensionalData.Dimension}} where {A <:
                                                                           FrequencyDim{<:RegularIndex}}

"""
    AbstractBispectrum{T, N, B}

A type alias for an `AbstractDimArray` in which the first and second dimensions are [`Freq`](@ref)uency.
"""
const AbstractBispectrum = AbstractDimArray{T, N, <:BifreqIndex, B} where {T, N, B}
freqs(x::AbstractBispectrum) = (dims(x, 1).val.data, dims(x, 2).val.data)

"""
    RegularBispectrum{T, N, B}

A type alias for a spectrum with a regularly sampled frequency index.
"""
const RegularBispectrum = AbstractDimArray{T, N, <:RegularBifreqIndex, B} where {T, N, B}

"""
    Bispectrum(f1, f2, x)

Constructs a univariate bispectrum with frequencies `f1` and `f2` and data `x`.
"""
Bispectrum(f1, f2, x; kwargs...) = DimArray(x, (Freq(f2), Freq(f2)); kwargs...)

"""
    Bispectrum(f1, f2, v, x)

Constructs a multivariate spectrum with frequencies `f1` and `f2`, variables `v`, and data `x`.
"""
Bispectrum(f1, f2, v, x; kwargs...) = DimArray(x, (Freq(f1), Freq(f2), Var(v)); kwargs...)
function Bispectrum(f1, f2, v::DimensionalData.Dimension, x; kwargs...)
    DimArray(x, (Freq(f1), Freq(f2), v); kwargs...)
end

function _bispectrum(x::AbstractVector, fs::Number,
                     f_min::Number = fs / min(length(x) ÷ 4, 1000); kwargs...)
    n = length(x)
    validfreqs = rfftfreq(n, fs)
    if f_min == 0
        f_min = validfreqs[2]
        nfft = floor(Int, n / 2) * 2
    else
        nfft = ceil(Int, ustripall(fs) / ustripall(f_min))
    end
    if ustripall(f_min) < ustripall(validfreqs[2])
        error("Cannot resolve an `f_min` of $f_min")
    end

    isodd(nfft) && (nfft += 1)
    overlap = nfft ÷ 2
    (nfft - overlap ≤ 0) && error("FFT padding is too high")
    n_segments = floor(Int, (n - nfft) / (nfft - overlap) + 1)

    # Get the type of the spectrum
    nfreqs = ((nfft) ÷ 2 + 1) ÷ 2
    S̄ = zeros(nfreqs, nfreqs, n_segments) .|> Complex
    @debug "Calculating spectrum for $n_segments segments of length $(nfft + padding)"
    for i in 1:n_segments
        start_idx = (i - 1) * (nfft - overlap) + 1
        end_idx = start_idx + nfft - 1
        segment = x[start_idx:end_idx] #.* hann_window
        segment = segment .- mean(segment)
        y = rfft(segment) / (nfft)
        yh = y[1:nfreqs]
        # * Compute F⁺(fs+f2) by assuming the frequencies are evenly spaced, monotonically increasing
        idxs = Iterators.product(eachindex(yh), eachindex(yh)) .|> sum
        S̄[:, :, i] .= ((yh * yh') .* conj(getindex.([y], idxs)))  # F(f1) * F(F2) * F⁺(fs+f2)
    end

    # Calculate the frequencies
    freqs = range(convertconst(0, fs), stop = fs / 4, length = size(S̄, 1))
    Bispectrum(freqs, freqs, Dim{:window}(1:n_segments), S̄; kwargs...)
end

function _bispectrum(x::AbstractTimeSeries;
                     f_min = samplingrate(x) / min(length(x) ÷ 4, 1000),
                     kwargs...)
    _bispectrum(x, samplingrate(x), f_min; kwargs...) ./ duration(x)
end
function bispectrum(x::AbstractTimeSeries, args...; kwargs...)
    S̄ = _bispectrum(x; kwargs...) ./ duration(x)
    dropdims(mean(S̄, dims = Dim{:window}); dims = Dim{:window}) .|> abs
end
