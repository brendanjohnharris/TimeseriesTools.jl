module SciMLBaseExt
using SciMLBase
import SciMLBase: AbstractTimeseriesSolution
using TimeseriesTools
import TimeseriesTools: Timeseries

# * Sucks that VectorOfArray no longer subtypes AbstractArray. Compatibility could have been so easy. Wait and see if the new indexing interface improves in the near future to.

function Timeseries(sol::AbstractTimeseriesSolution)
    Timeseries(sol.u, sol.t; metadata = Dict(:sol => sol)) # Find a way to use ranges if t is a range
end

function MultivariateTimeseries(sol::AbstractTimeseriesSolution{T, 2}) where {T}
end

# ! Add method for parameterized functions to add symbols as variable numbers
end # module
