using MakieCore
using LaTeXStrings
using DimensionalData
import MakieCore.plot!
import GeometryBasics

export dimname, decompose

MakieCore.@recipe(SpectrumPlot, x, y) do scene
    Theme(
        plot_color = :cornflowerblue
    )
end
function plot!(p::SpectrumPlot)
    x = p[:x]
    y = p[:y]
    xs = MakieCore.pseudolog10(first(x.val[x.val .> 0])/2)
    ys = log10
    lines!(p, x, y, color = p[:plot_color])
    p
end
const ToolSpectrumPlot = SpectrumPlot{Tuple{<:AbstractSpectrum}}
argument_names(::Type{<: ToolSpectrumPlot}) = (:x,)
function plot!(p::ToolSpectrumPlot)
    x = collect(dims(p[:x], Freq))
    y = collect(p[:x])
    spectrumplot!(p, x, y)
    p
end



function plotLFPspectra(X::AbstractTimeSeries; slope=nothing, position=Point2f([5, 1e-5]), fs=nothing, N=1000, slopecolor=:crimson, kwargs...)
    times = collect(dims(X, Ti))
    if isnothing(fs)
        Δt = times[2] - times[1]
        all(Δt .≈ diff(times)) || @warn "Violated assumption: all(Δt .≈ diff(times))"
    else
        Δt = 1/fs
    end

    P = [fp(Array(x)) for x ∈ eachcol(X)]
    𝑓 = P[1].freq # Should be pretty much the same for all columns?
    psd = hcat([p.power for p ∈ P]...)
    psd = psd./(sum(psd, dims=1).*(𝑓[2] - 𝑓[1]))
    psd = DimArray(psd, (Dim{:frequency}(𝑓), dims(X, :channel)))
    fig = traces(𝑓, Array(psd); xlabel="𝑓 (Hz)", ylabel="Ŝ", title="Normalised power spectral density", smooth=1, yscale=log10, doaxis=false, domean=false, yminorgridvisible=false, kwargs...)
    if !isnothing(slope)
        _psd = psd[Dim{:frequency}(DD.Between(slope...))]
        c, r, f = powerlawfit(_psd)
        lines!(LinRange(slope..., 100), f(LinRange(slope..., 100)), color=slopecolor, linewidth=5)
        text!(L"$\alpha$= %$(round(r, sigdigits=2))", position=Point2f0(position), fontsize=40)
    end
    return fig
end

GeometryBasics.decompose(x::Union{<:AbstractTimeSeries, <:AbstractSpectrum}) = ((dims(x).|>collect)..., x.data)

MakieCore.convert_arguments(P::MakieCore.PointBased, x::UnivariateTimeSeries) = MakieCore.convert_arguments(P, decompose(x)...)

MakieCore.convert_arguments(P::MakieCore.SurfaceLike, x::MultivariateTimeSeries) = MakieCore.convert_arguments(P, decompose(x)...)
MakieCore.convert_arguments(P::Type{<:MakieCore.Heatmap}, x::MultivariateTimeSeries) = MakieCore.convert_arguments(P, decompose(x)...)



# GeometryBasics.decompose(x::AN.DimensionalData.AbstractDimArray, dims...) = (getindex.((dims(x).|>collect), dims)..., x.data[dims...])

dimname(x::DimArray, dim) = dims(x, dim)|>name|>string

# formataxes(x::AN.DimensionalData.AbstractDimArray{T, 2} where T) = (xlabel=dimname(x, 1), ylabel=dimname(x, 2))
# formataxes(x::AN.DimensionalData.AbstractDimArray{T, 1} where T) = (xlabel=dimname(x, 1),)


# TODO Automatic unit stripping and axis units
# TODO Stacked traces, traces recipe
