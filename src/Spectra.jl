using FFTW
using Statistics
using Unitful
using UnPack
using ComponentArrays

import TimeseriesBase.Operators.𝒯

export spectrum, energyspectrum, powerspectrum, _energyspectrum, _powerspectrum,
       colorednoise, logsample, logbin, oneoneff, fit_oneoneff

function _periodogram(x::AbstractVector, fs::Number,
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
        throw(DomainError(f_min, "Cannot resolve an `f_min` of $f_min"))
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
            segment = Timeseries([segment.data; zeros(padding) * unit(eltype(segment))],
                                 padts)
        end

        y = rfft(segment) / (nfft + padding)
        y |> eltype |> unit == NoUnits && (y = y * u)
        S̄[:, i] .= (abs.(y) .^ 2) / A
    end

    # Calculate the frequencies
    freqs = range((0)unit(fs), stop = fs / 2, length = size(S̄, 1))
    df = step(freqs)

    # Normalize the mean energy spectrum to obey Parseval's theorem
    meanS̄ = mean(S̄, dims = 2)
    S̄ = S̄ ./ ustripall((sum(meanS̄) - 0.5 .* meanS̄[1]) .* df) # Subtract the zero frequency component twice, so that it doesn't bias when we divide by a half
    S̄ = 0.5 * S̄ * ustripall(sum(x .^ 2) ./ fs) # Normalized to have total energy equal to energy of signal. Ala parseval. 0.5 because we only return the positive half of the spectrum.
    Spectrum(freqs, Dim{:window}(1:n_segments), S̄; kwargs...)
end

"""
    _energyspectrum(x::RegularTimeseries, f_min=samplingrate(x)/min(length(x)÷4, 1000))

Computes the energy spectrum of a regularly sampled time series `x` with an optional minimum frequency `f_min`.
"""
function _energyspectrum(x::typeintersect(RegularTimeseries, UnivariateTimeseries),
                         f_min::Number = samplingrate(x) / min(length(x) ÷ 4, 1000);
                         kwargs...)
    return _periodogram(x, samplingrate(x), f_min; kwargs...)
end

"""
    _energyspectrum(x::RegularTimeseries, f_min=0)

Computes the energy spectrum of a time series using the fast Fourier transform.

If `f_min > 0`, the energy spectrum is calculated for windows of the time series determined by `f_min`,  the minimum frequency that will be resolved in the spectrum.
If `f_min > 0`, the second dimension of the output will correspond to the windows. For an averaged periodogram, see [`energyspectrum`](@ref).

If the input time series is a [`UnitfulTimeseries`](@ref), the frequency will also have units.
Moreover if the elements of `x` are unitful, so are the elements of the spectrum.

# Examples
```@example 1
julia> using TimeseriesTools
julia> t = range(0.0, stop=1.0, length=1000);
julia> x = sin.(2 * π * 5 * t);
julia> ts = RegularTimeseries(x, t);
julia> S = _energyspectrum(ts);
julia> S isa MultivariateSpectrum
```
"""
function _energyspectrum(x::MultivariateTimeseries, args...; kwargs...)
    X = [_energyspectrum(_x, args...; kwargs...)
         for _x in eachslice(x, dims = 2)]
    return ToolsArray(X, dims(x, 2)) |> stack
end

"""
    energyspectrum(x::RegularTimeseries, f_min=0; kwargs...)

Computes the average energy spectrum of a regularly sampled time series `x`.
`f_min` determines the minimum frequency that will be resolved in the spectrum.
See [`_energyspectrum`](@ref).
"""
function energyspectrum(x, args...; kwargs...)
    dropdims(mean(_energyspectrum(x, args...; kwargs...), dims = Dim{:window});
             dims = Dim{:window})
end

"""
    _powerspectrum(x::AbstractTimeseries, f_min=samplingrate(x)/min(length(x)÷4, 1000); kwargs...)

Computes the power spectrum of a time series `x` in Welch periodogram windows.
Note that the `_powerspectrum` is simply the [`_energyspectrum`](@ref) divided by the duration of each window.
See [`_energyspectrum`](@ref).
"""
function _powerspectrum(x::AbstractTimeseries, args...; kwargs...)
    S̄ = _energyspectrum(x, args...; kwargs...)
    return S̄ ./ duration(x)
end

"""
    powerspectrum(x::AbstractTimeseries, f_min=samplingrate(x)/min(length(x)÷4, 1000); kwargs...)

Computes the average power spectrum of a time series `x` using the Welch periodogram method.
"""
function powerspectrum(x::AbstractTimeseries, args...; kwargs...)
    dropdims(mean(_powerspectrum(x, args...; kwargs...), dims = Dim{:window});
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
- A [`Timeseries`](@ref) containing the generated colored noise.

# Example

```julia
using TimeseriesTools
pink_noise = colorednoise(1:0.01:10; α=1.0)
pink_noise isa RegularTimeseries
```
"""
function colorednoise(ts::AbstractRange{T}, args...; α = 2.0,
                      kwargs...) where {T}
    u = unit(T)
    ts = T <: Quantity ? ustrip(ts) : ts
    dt = step(ts)
    f = rfftfreq(length(ts), dt)
    x̂ = sqrt.(1.0 ./ f .^ α) .* exp.(2π .* rand(length(f)) * im)
    x̂[1] = 0
    x = irfft(x̂, length(ts))
    dt = length(ts) * step(f)
    t = range(0, (length(x) - 1) * dt, length = length(x))
    @assert all(t .+ first(ts) .≈ ts)
    Timeseries(x, ts * u, args...; kwargs...)
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
    S̄ = S̄ * ustripall(length(t)) # Normalized to have total energy equal to energy of signal. Ala parseval.
    Spectrum(frange, Dim{:window}([1]), Matrix(S̄')'; kwargs...)
end

function _energyspectrum(x::SpikeTrain{T, 1} where {T}, frange::Tuple; kwargs...)
    _energyspectrum(x, 0:first(frange):last(frange); kwargs...)
end

function logbin(_s::UnivariateSpectrum{T, N, D}) where {T, N, D <: Tuple{<:𝑓}}
    if first(freqs(_s)) <= 0
        throw(ArgumentError("Frequencies must be positive"))
    end
    if !issorted(freqs(_s))
        throw(ArgumentError("Frequencies must be ascending"))
    end
    _s = set(map(log10, _s), 𝑓 => map(log10, freqs(_s)))
    _f = freqs(_s)
    df = _f[2] - _f[1] # The smallest sensible binning
    f = range(first(_f) - df / 2, stop = last(_f), step = df)
    bins = intervals(f)
    s = groupby(_s, 𝑓 => Bins(bins))
    return set(s, 𝑓 => Log10𝑓(f))
end
function logsample(_s::A, average = mean) where {A <: AbstractDimArray}
    _s = logbin(_s)
    return map(average, _s)
end

function oneoneff(log_f, p)
    @unpack log_b, β, log_k, log_c, peaks = p
    k = exp10(log_k) # Knee frequency
    f = exp10.(log_f) # Frequency
    b = exp10(log_b) # Power law amplitude
    c = exp10(log_c) # Noise floor

    # power_law = @. log_b - log_c - log(k + f^β)
    s = @. b / (k + f^β) + c

    # Add peaks
    for p in peaks
        f_peak = exp10(p[:log_f])
        A_peak = exp10(p[:log_A])
        σ_peak = f_peak * tanh(p[:σ̃]) # log_σ gives a constant width in log_f space
        # σ̃ = log[(f_peak - σ)/(f_peak + σ)] = atanh(σ / f_peak)
        @. s += A_peak * exp(-(f - f_peak)^2 / (2 * σ_peak^2))
    end

    return log10.(s)
end

"""
The one-argument form makes a guess at initial parameters. The two-argument form runs
optimization from the given parameters as a starting point.
"""
function fit_oneoneff(logspectrum::AbstractDimVector;
                      w = maximum(1, length(logspectrum) ÷ 100),
                      n_peaks = nothing,
                      minprom = (maximum(logspectrum) - minimum(logspectrum)) / 50,
                      kwargs...)
    log_f, log_s = lookup(logspectrum, 1), parent(logspectrum)

    log_c = minimum(log_s)
    log_b = first(log_s)
    β = -last([ones(length(log_f)) log_f] \ log_s) # Simple linear regression
    log_k = -1 # A guess

    # * Find peaks by looking for local maxima
    _, proms, bounds = findpeaks(logspectrum, w; minprom, kwargs...)

    if !isnothing(n_peaks)
        if n_peaks > length(proms)
            proms = vcat(proms, [mean(log_s) for _ in 1:(n_peaks - length(proms))])
            bounds = vcat(bounds,
                          [deepcopy(first(bounds))
                           for _ in 1:(n_peaks - length(bounds))])
        end

        ps = sortperm(proms; rev = true)[1:n_peaks]
        proms = proms[ps]
        bounds = bounds[ps]
    else
        n_peaks = length(proms)
    end

    peaks = map(proms, bounds) do prom, bound
        log_f = mean(bound)
        σ̃ = (maximum(bound) - minimum(bound)) / 2
        s_f = logspectrum[Near(maximum(bound) + σ̃)]
        log_A = prom + s_f
        return ComponentArray(; log_f, σ̃, log_A)
    end

    return ComponentArray(; log_b, β, log_k, log_c, peaks)
end
