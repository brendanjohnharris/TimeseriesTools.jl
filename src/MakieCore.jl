import MakieCore
import MakieCore.Observables
using LaTeXStrings
using DimensionalData
import GeometryBasics.decompose

export dimname, decompose, spectrumplot!, spectrumplot

MakieCore.@recipe(SpectrumPlot, x, y) do scene
    MakieCore.Theme(;
                    color = :cornflowerblue,
                    peaks = false)
end

function MakieCore.plot!(p::SpectrumPlot{<:Tuple{<:AbstractVector, <:AbstractVector}})
    x = p[:x]
    y = p[:y]
    # idxs = map((x, y)->(x .> 0) .& (y .> 0), x, y)
    # _x = map((x, i)->x[i], x, idxs)
    # _y = map((y, i)->y[i], y, idxs)
    MakieCore.lines!(p, x, y; p.attributes...)
    p
end

"""
    decompose(x::Union{<:AbstractTimeSeries, <:AbstractSpectrum})
Convert a time series or spectrum to a tuple of the dimensions and the data (as `Array`s).
"""
decompose(x::Union{<:AbstractTimeSeries, <:AbstractSpectrum}) = (lookup(x)...,
                                                                 x.data)
# function decompose(x::MultivariateTimeSeries)
#     if size(x)[2] < 3
#         x |> eachcol |> collect .|> collect |> Tuple
#     else
#         ((dims(x) .|> collect)..., x.data)
#     end
# end
function MakieCore.convert_arguments(P::MakieCore.PointBased, x::UnivariateTimeSeries)
    MakieCore.convert_arguments(P, decompose(x)...)
end
function MakieCore.convert_arguments(P::MakieCore.ImageLike, x::MultivariateTimeSeries)
    MakieCore.convert_arguments(P, decompose(x)...)
end
function MakieCore.convert_arguments(P::Type{<:MakieCore.Heatmap},
                                     x::MultivariateTimeSeries)
    MakieCore.convert_arguments(P, decompose(x)...)
end

function MakieCore.convert_arguments(t::Type{<:MakieCore.AbstractPlot},
                                     D, A::DimensionalData.AbstractDimMatrix)
    xs = parent(lookup(A, D))
    dim = dims(A, D)
    ys = parent(first(lookup(A)[dims(A) .!= [dim]]))
    return xs, ys, A.data
end
function MakieCore.convert_arguments(t::Type{<:MakieCore.AbstractPlot},
                                     D, ys,
                                     A::DimensionalData.AbstractDimMatrix)
    xs = parent(lookup(A, D))
    return xs, ys, A.data
end

# GeometryBasics.decompose(x::AN.DimensionalData.AbstractDimArray, dims...) = (getindex.((dims(x).|>collect), dims)..., x.data[dims...])

function formataxislabels(x::UnivariateTimeSeries)
    (; xlabel = describedims(x, 1), ylabel = describename(x))
end
function formataxislabels(x::MultivariateTimeSeries)
    (; xlabel = describedims(x, 1), ylabel = describedims(x, 2))
end

# TODO Stacked traces, traces recipe
