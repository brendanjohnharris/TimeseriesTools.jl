module TimeseriesTools

import Unitful.unit

using Reexport
using DimensionalData
using IntervalSets
@reexport using DimensionalData
@reexport using IntervalSets
@reexport using Normalization

function __init__()
    ENV["UNITFUL_FANCY_EXPONENTS"] = true
end

include("Types.jl")
include("Utils.jl")
include("Operators.jl")
include("SpikeTrains.jl")
include("Spectra.jl")
include("Spectrograms.jl")
include("Unitful.jl")
include("Dates.jl")
include("MakieCore.jl")
include("IO.jl")
include("TimeseriesSurrogates.jl")

bandpass(x::AbstractTimeSeries) = x
highpass(x::AbstractTimeSeries) = x
lowpass(x::AbstractTimeSeries) = x

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
function progressmap end
function msdist end
export _waveletfreqs, _waveletspectrogram, waveletspectrogram, progressmap, msdist

function spectrumplot! end
function spectrumplot end
function trajectory! end
function trajectory end
function shadows! end
function traces! end
function traces end
function stackedtraces! end
function stackedtraces end
export spectrumplot!, spectrumplot, trajectory!, trajectory, shadows!, traces!, traces,
       stackedtraces!, stackedtraces

end
