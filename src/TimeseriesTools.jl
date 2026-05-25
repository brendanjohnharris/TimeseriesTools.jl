module TimeseriesTools

import Unitful.unit

using Reexport
using DimensionalData
using IntervalSets
@reexport using IntervalSets
@reexport using DimensionalData
@reexport using TimeseriesBase
@reexport using Normalization
import StatsAPI: fit, fit!, predict
export fit, fit!, predict
import Normalization: params!

include("Utils.jl")
include("Interpolate.jl")
include("SpikeTrains.jl")
include("Spectra.jl")
include("Mapple.jl")
include("Unitful.jl")
include("TimeseriesSurrogates.jl")

bandpass(x::AbstractTimeseries) = x
highpass(x::AbstractTimeseries) = x
lowpass(x::AbstractTimeseries) = x

function timescale(x::UnivariateTimeseries; method = :ac_crossing)
    return timescale(x::UnivariateTimeseries, Val{method}())
end

# ? Placeholder functions for extensions
function phasestitch end
function isoamplitude end
function analyticamplitude end
function analyticphase end
function instantaneousfreq end
instantfreq = instantaneousfreq
export phasestitch, bandpass, isoamplitude, analyticphase, analyticamplitude,
    instantaneousfreq, instantfreq, highpass, lowpass

function _waveletfreqs end
function _waveletspectrogram end
function waveletspectrogram end
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
function bootstrapaverage end
function bootstrapmedian end
function bootstrapmean end
export bootstrapaverage, bootstrapmedian, bootstrapmean

end
