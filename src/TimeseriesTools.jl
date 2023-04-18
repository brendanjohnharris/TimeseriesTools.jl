module TimeseriesTools
using Reexport
using DimensionalData
using IntervalSets
import Unitful.unit
@reexport using DimensionalData
@reexport using IntervalSets

# function __init__()
    # @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin
    #     @eval include("Unitful.jl")
    # end
    # @require Unitful="20f20a25-4f0e-4fdf-b5d1-57303727442b" begin
    #     @eval include("Makie.jl")
    # end
# end

include("Types.jl")
# include("Utils.jl")
# include("Spectra.jl")
# include("Unitful.jl")
# include("Makie.jl")



end
