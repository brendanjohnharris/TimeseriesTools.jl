module TimeseriesTools
using Reexport
using Requires
using DimensionalData
using IntervalSets
import Unitful.unit
@reexport using DimensionalData
@reexport using IntervalSets

function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
        @eval include("Makie.jl")
    end
end

include("Types.jl")
include("Utils.jl")
include("Spectra.jl")
include("Unitful.jl")
include("MakieCore.jl")



end
