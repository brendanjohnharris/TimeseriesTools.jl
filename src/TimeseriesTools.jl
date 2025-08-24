module TimeseriesTools

import Unitful.unit

using Reexport
using DimensionalData
using IntervalSets
using ProgressLogging
@reexport using IntervalSets
@reexport using DimensionalData
@reexport using TimeseriesBase
@reexport using Normalization
@reexport using TimeseriesMakie

include("Utils.jl")
include("SpikeTrains.jl")
include("Spectra.jl")
include("Unitful.jl")
include("TimeseriesSurrogates.jl")

bandpass(x::AbstractTimeseries) = x
highpass(x::AbstractTimeseries) = x
lowpass(x::AbstractTimeseries) = x

function timescale(x::UnivariateTimeseries; method = :ac_crossing)
    timescale(x::UnivariateTimeseries, Val{method}())
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
function interpolate end
function msdist end
export _waveletfreqs, _waveletspectrogram, waveletspectrogram, msdist

end
