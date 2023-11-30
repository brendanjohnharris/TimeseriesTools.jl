module TimeseriesTools

import Unitful.unit

using Reexport
using DimensionalData
using IntervalSets
# if !isdefined(Base, :get_extension)
using Requires
# end
@reexport using DimensionalData
@reexport using IntervalSets
@reexport using Normalization

function __init__()
    ENV["UNITFUL_FANCY_EXPONENTS"] = true
    # @static if !isdefined(Base, :get_extension)
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
        @eval include("../ext/MakieExt.jl")
    end
    @require DSP="717857b8-e6f2-59f4-9121-6e50c889abd2" begin
        @eval include("../ext/DSPExt.jl")
    end
    # end
end

include("Types.jl")
include("Utils.jl")
include("Spectra.jl")
include("Unitful.jl")
include("Dates.jl")
include("MakieCore.jl")
include("IO.jl")

bandpass(x::AbstractTimeSeries) = x
function phasestitch end
function isoamplitude end
function analyticamplitude end
function analyticphase end
export phasestitch, bandpass, isoamplitude, analyticphase, analyticamplitude

end
