module TimeseriesTools
using Requires
using Reexport
using DimensionalData
using IntervalSets
@reexport using DimensionalData
@reexport using IntervalSets

# function __init__()
#     @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin
#         @eval include("Unitful.jl")
#     end
# end

include("Types.jl")
include("Utils.jl")
include("Spectra.jl")



end
