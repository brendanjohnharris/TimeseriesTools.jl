import MakieCore
import MakieCore.Observables
using LaTeXStrings
using DimensionalData
import GeometryBasics.decompose

export dimname, decompose, spectrumplot!, spectrumplot

MakieCore.@recipe(SpectrumPlot, x, y) do scene
    MakieCore.Theme(;
        color=:cornflowerblue
    )
end

function MakieCore.plot!(p::SpectrumPlot{<:Tuple{<:AbstractVector,<:AbstractVector}}, kwargs...)
    x = p[:x]
    y = p[:y]
    # idxs = map((x, y)->(x .> 0) .& (y .> 0), x, y)
    # _x = map((x, i)->x[i], x, idxs)
    # _y = map((y, i)->y[i], y, idxs)
    MakieCore.lines!(p, x, y; color=p[:color], kwargs...)
    p
end


"""
    decompose(x::Union{<:AbstractTimeSeries, <:AbstractSpectrum})
Convert a time series or spectrum to a tuple of the dimensions and the data (as `Array`s).
"""
decompose(x::Union{<:AbstractTimeSeries,<:AbstractSpectrum}) = ((dims(x) .|> collect)..., x.data)
MakieCore.convert_arguments(P::MakieCore.PointBased, x::UnivariateTimeSeries) = MakieCore.convert_arguments(P, decompose(x)...)
MakieCore.convert_arguments(P::MakieCore.SurfaceLike, x::MultivariateTimeSeries) = MakieCore.convert_arguments(P, decompose(x)...)
MakieCore.convert_arguments(P::Type{<:MakieCore.Heatmap}, x::MultivariateTimeSeries) = MakieCore.convert_arguments(P, decompose(x)...)

# GeometryBasics.decompose(x::AN.DimensionalData.AbstractDimArray, dims...) = (getindex.((dims(x).|>collect), dims)..., x.data[dims...])

formataxislabels(x::UnivariateTimeSeries) = (; xlabel=describedims(x, 1), ylabel=describename(x))
formataxislabels(x::MultivariateTimeSeries) = (; xlabel=describedims(x, 1), ylabel=describedims(x, 2))


# TODO Stacked traces, traces recipe
