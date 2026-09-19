module TimeseriesTools

import Unitful.unit

using Reexport
using DimensionalData
@reexport using DimensionalData
@reexport using TimeseriesBase
@reexport using Normalization
import StatsAPI: fit, fit!, predict, vcov, stderror
export fit, fit!, predict
import Normalization: params!


include("Utils.jl")
include("Interpolate.jl")
include("SpikeTrains.jl")
include("Spectra.jl")
include("Mapple.jl")
include("Unitful.jl")

"""
    bandpass(x, pass; designmethod = DSP.Butterworth(4))
    bandpass(x, fs, pass; kwargs...)

Bandpass filter `x` over the band `pass`, given as a two-element vector, tuple or interval.
Provided by `DSPExt` (loaded with `using DSP`), which applies a zero-phase forward-backward
(`filtfilt`) filter along the first dimension. The sampling rate is read from the time lookup of a
`RegularTimeseries`, or supplied as `fs` for a plain array.

The single-argument form is a no-op that returns `x` unchanged, so a pipeline that filters
optionally stays valid without `DSP` loaded; calling with a `pass` band and no `DSP` is a
`MethodError`.
"""
bandpass(x::AbstractTimeseries) = x

"""
    highpass(x, pass::Number; designmethod = DSP.Butterworth(4))
    highpass(x, fs, pass; kwargs...)

Highpass filter `x` above `pass`. As [`bandpass`](@ref), but with a single corner frequency.
"""
highpass(x::AbstractTimeseries) = x

"""
    lowpass(x, pass::Number; designmethod = DSP.Butterworth(4))
    lowpass(x, fs, pass; kwargs...)

Lowpass filter `x` below `pass`. As [`bandpass`](@ref), but with a single corner frequency.
"""
lowpass(x::AbstractTimeseries) = x

"""
    timescale(x::UnivariateTimeseries; method = :ac_crossing)

Characteristic timescale of `x`.

!!! warning "Not implemented"
    Only the dispatch stub exists: no `method` is currently implemented, so calling this throws a
    `MethodError`. It is retained as a placeholder for an extension to fill.
"""
function timescale(x::UnivariateTimeseries; method = :ac_crossing)
    return timescale(x::UnivariateTimeseries, Val{method}())
end

# ? Placeholder functions for extensions
function phasestitch end

"""
    isoamplitude(x; dims = 1)

Strip the amplitude from `x`, keeping only its phase: `sin` of the analytic phase, so the result
oscillates over `[-1, 1]` with the same instantaneous phase as `x`. Provided by `DSPExt` (loaded
with `using DSP`).
"""
function isoamplitude end

"""
    analyticamplitude(x)

Envelope of `x`: the modulus of its analytic signal (`abs.(hilbert(x))`). Provided by `DSPExt`
(loaded with `using DSP`). Most meaningful on a narrowband signal, so filter with
[`bandpass`](@ref) first.
"""
function analyticamplitude end

"""
    analyticphase(x)

Instantaneous phase of `x`: the argument of its analytic signal (`angle.(hilbert(x))`), wrapped to
`(-π, π]`. Provided by `DSPExt` (loaded with `using DSP`). Use `DSP.unwrap` for a continuous phase,
or [`instantaneousfreq`](@ref) for its rate of change.
"""
function analyticphase end

"""
    instantaneousfreq(x)
    instantfreq(x)

Instantaneous frequency of `x` in cycles per unit time: the central derivative of the unwrapped
[`analyticphase`](@ref), divided by `2π`. Provided by `DSPExt` (loaded with `using DSP`).
Meaningful for a narrowband signal, so filter with [`bandpass`](@ref) first.
"""
function instantaneousfreq end
instantfreq = instantaneousfreq
export phasestitch, bandpass, isoamplitude, analyticphase, analyticamplitude,
    instantaneousfreq, instantfreq, highpass, lowpass

function _waveletfreqs end
function _waveletspectrogram end
"""
    waveletspectrogram(x::RegularTimeseries, args...; moth = Morlet(2π), β = 1, Q = 32, pass = nothing)

Continuous wavelet transform of `x`, returned as a `RegularSpectrogram` indexed by time and
frequency. Provided by `ContinuousWaveletsExt` (loaded with `using ContinuousWavelets`).

`moth` selects the mother wavelet, `Q` the number of voices per octave, and `β` the spacing of
scales (`β = 1` is linear in log space). `pass` limits the returned frequencies to those below its
maximum. A `MultivariateTimeseries` is transformed column by column, giving a time x frequency x
variable array.
"""
function waveletspectrogram end
"""
    msdist(x, lags = 0:(length(x) - 1))

Mean-squared displacement of `x` over `lags`, computed by FFT. Provided by
`AutocorrelationsExt` (loaded with `using Autocorrelations`). For a `RegularTimeseries` the result
is a `Timeseries` indexed by lag in time units.

The MSD relies on a variance and so diverges for heavy-tailed processes; see [`madev`](@ref) for an
alternative that stays finite.
"""
function msdist end
function resample end
"""
    impute(x, interp = AkimaInterpolation, args...; dims = 1, replace = [NaN, Nothing, Missing], kwargs...)

Fill flagged entries of `x` by interpolation. Method-only function: provided by
`DataInterpolationsExt` (loaded with `using DataInterpolations`).

Entries matching any element of `replace` (sentinel values by `isequal`/`isnan`, types
by `isa`) are set to `missing`, an interpolant is fit to the survivors, and the result
is evaluated at every original time point. For arrays of more than one dimension, each
slice along `dims = 1` is imputed independently.
"""
function impute end
export _waveletfreqs, _waveletspectrogram, waveletspectrogram, msdist, resample, impute

# * BootstrapExt
"""
    bootstrapaverage(average, x; confint = 0.95, N = 10000, dims = 1)

Bootstrap an arbitrary `average` function over `x`, returning
`(; average, confint = (; lower, upper))`. Provided by `BootstrapExt` (loaded with
`using Bootstrap`).

Uses balanced sampling with `N` resamples and a bias-corrected accelerated (BCa) interval at the
`confint` level. `NaN`s are dropped, and a slice with fewer than five finite values returns `NaN`s
rather than an estimate. For an array of more than one dimension, each slice along `dims` is
bootstrapped independently and the results are collected into arrays of the same shape.
"""
function bootstrapaverage end

"""
    bootstrapmedian(x; kwargs...)

[`bootstrapaverage`](@ref) with `median`.
"""
function bootstrapmedian end

"""
    bootstrapmean(x; kwargs...)

[`bootstrapaverage`](@ref) with `mean`.
"""
function bootstrapmean end
export bootstrapaverage, bootstrapmedian, bootstrapmean

function hint_ext(io, f, pkg)
    ext = pkg * "Ext"
    print(io, "\n")
    printstyled(io, "Hint:"; color = :cyan, bold = true)
    return print(io, " calling `$(f)` requires the $(ext) extension. Run `using $(pkg)` to enable it.\n")
end

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if exc.f === phasestitch || exc.f === isoamplitude || exc.f === analyticamplitude || exc.f === analyticphase || exc.f === instantaneousfreq || exc.f === downsample
            hint_ext(io, exc.f, "DSP")
        elseif exc.f === _waveletfreqs || exc.f === _waveletspectrogram || exc.f === waveletspectrogram
            hint_ext(io, exc.f, "ContinuousWavelets")
        elseif exc.f === msdist
            hint_ext(io, exc.f, "Autocorrelations")
        elseif exc.f === resample || exc.f === impute || exc.f === interpolate || exc.f === upsample
            hint_ext(io, exc.f, "DataInterpolations")
        elseif exc.f === pointprocess! || exc.f === gammarenewal! || exc.f === gammarenewal
            hint_ext(io, exc.f, "Distributions")
        elseif exc.f === bootstrapaverage || exc.f === bootstrapmedian || exc.f === bootstrapmean
            hint_ext(io, exc.f, "Bootstrap")
        end
    end
end

end
