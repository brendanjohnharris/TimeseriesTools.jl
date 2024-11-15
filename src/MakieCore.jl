import MakieCore
import MakieCore.Observables
import MakieCore.theme
using LaTeXStrings
using DimensionalData
import GeometryBasics.decompose

export dimname, decompose, spectrumplot!, spectrumplot

MakieCore.@recipe(SpectrumPlot, x, y) do scene
    MakieCore.Attributes(;
                         color = :cornflowerblue,
                         peaks = false,
                         linewidth = theme(scene, :linewidth),
                         alpha = 1,
                         linestyle = theme(scene, :linestyle),
                         linecap = theme(scene, :linecap),
                         joinstyle = theme(scene, :joinstyle))
end

function MakieCore.plot!(p::SpectrumPlot{<:Tuple{<:AbstractVector, <:AbstractVector}})
    x = p[:x]
    y = p[:y]
    # idxs = map((x, y)->(x .> 0) .& (y .> 0), x, y)
    # _x = map((x, i)->x[i], x, idxs)
    # _y = map((y, i)->y[i], y, idxs)
    lineattrs = [:linewidth, :alpha, :linestyle, :linecap, :joinstyle]
    MakieCore.lines!(p, x, y; [a => p.attributes[a] for a in lineattrs]...)
    p
end

"""
    decompose(x::Union{<:AbstractTimeSeries, <:AbstractSpectrum})
Convert a time series or spectrum to a tuple of the dimensions and the data (as `Array`s).
"""
decompose(x::AbstractToolsArray) = (lookup(x)..., x.data) .|> parent
# function MakieCore.convert_arguments(P::Type{<:MakieCore.AbstractPlot},
#                                      x::AbstractTimeSeries)
#     MakieCore.convert_arguments(P, decompose(x)...)
# end
# function MakieCore.convert_arguments(P::Type{<:MakieCore.AbstractPlot},
#                                      x::AbstractSpectrum)
#     MakieCore.convert_arguments(P, decompose(x)...)
# end
# function MakieCore.convert_arguments(P::Type{<:MakieCore.AbstractPlot},
#                                      x::Tuple{<:Union{<:AbstractTimeSeries,
#                                                       <:AbstractSpectrum}})
#     MakieCore.convert_arguments(P, decompose(first(x))...)
# end
MakieCore.plottype(::AbstractTimeSeries) = Lines

function MakieCore.plot!(ax, x::UnivariateTimeSeries; kwargs...)
    MakieCore.lines!(ax, x; kwargs...)
end
function MakieCore.plot(x::UnivariateTimeSeries; kwargs...)
    MakieCore.lines(x; kwargs...)
end

function formataxislabels(x::UnivariateTimeSeries)
    (; xlabel = describedims(x, 1), ylabel = describename(x))
end
function formataxislabels(x::MultivariateTimeSeries)
    (; xlabel = describedims(x, 1), ylabel = describedims(x, 2))
end
