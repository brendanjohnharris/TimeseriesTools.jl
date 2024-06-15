# module MakieExt

import MakieCore.plot!
import MakieCore.plot

import ..Makie
import ..Makie: plot, plot!, lift, lines, lines!, band, band!, FigureAxisPlot, @lift,
                Observable, @recipe, Theme, Figure, Axis, AbstractPlot, Scatter, Lines,
                ScatterLines, Hexbin, Stem, Plot, scatter!, text!, with_theme, Attributes,
                theme
PointLike = Union{Scatter, Lines, ScatterLines, Hexbin, Stem}

using TimeseriesTools
using Statistics
using TimeseriesTools.Normalization
using Peaks

export spectrumplot!, spectrumplot, trajectory!, trajectory, shadows!

"""
    spectrumplot!(ax::Axis, x::UnivariateSpectrum)
Plot the given spectrum, labelling the axes, adding units if appropriate, and other niceties.
"""
function spectrumplot!(ax::Makie.Axis, x::UnivariateSpectrum; peaks = false, kwargs...)
    s = x
    uf = frequnit(x)
    ux = unit(x)
    f, x = decompose(x)
    f = ustripall.(f) |> collect
    x = ustripall.(x) |> collect
    idxs = (f .> 0) .& (x .> 0)
    ax.limits = ((minimum(f[idxs]), nothing), (minimum(x[idxs]), nothing))
    ax.xscale = log10
    ax.yscale = log10
    uf == NoUnits ? (ax.xlabel = "Frequency") : (ax.xlabel = "Frequency ($uf)")
    ux == NoUnits ? (ax.ylabel = "Spectral density") :
    (ax.ylabel = "Spectral density ($ux)")
    p = spectrumplot!(ax, f[idxs], x[idxs]; kwargs...)

    if peaks != false
        peaks == true && (peaks = 1)
        pks, vals = findmaxima(s, 10)
        pks, proms = peakproms(pks, s)
        if peaks isa Int
            promidxs = partialsortperm(proms, 1:peaks, rev = true)
        elseif peaks isa Real
            promidxs = (proms ./ vals .> peaks) |> collect
        end
        pks = pks[promidxs]
        pks = TimeseriesTools.freqs(s)[pks]
        vals = s[Freq(At(pks))]
        scatter!(ax, ustripall.(pks), collect(ustripall.(vals)),
                 color = Makie.current_default_theme().textcolor,
                 markersize = 10,
                 marker = :dtriangle)
        if eltype(pks) <: Quantity
            txt = string.(round.(eltype(pks), pks; digits = 3))
        else
            txt = string.(round.(pks; digits = 3))
        end
        text!(ax, ustripall.(pks), collect(ustripall.(vals));
              text = txt,
              align = (:center, :bottom), color = Makie.current_default_theme().textcolor,
              rotation = 0, fontsize = 12,
              offset = (0, 3))
    end
    p
end

"""
    spectrumplot!(ax::Axis, x::MultivariateSpectrum)
Plot the given spectrum, labelling the axes, adding units if appropriate, and adding a band to show the iqr
"""
function spectrumplot!(ax::Makie.Axis, x::MultivariateSpectrum;
                       peaks = false,
                       bandcolor = nothing,
                       percentile = 0.25, kwargs...)
    uf = frequnit(x)
    ux = unit(x)
    f, _, x = decompose(x)
    f = ustripall.(f) |> collect
    x = ustripall.(x) |> collect
    xmin = minimum(x, dims = 2) |> vec
    xmed = median(x, dims = 2) |> vec
    σₗ = mapslices(x -> quantile(x, percentile), x, dims = 2) |> vec
    σᵤ = mapslices(x -> quantile(x, 1 - percentile), x, dims = 2) |> vec
    idxs = (f .> 0) .& (xmin .> 0)
    ax.limits = ((minimum(f[idxs]), nothing), (minimum(σₗ[idxs]), nothing))
    ax.xscale = log10
    ax.yscale = log10
    if isempty(ax.xlabel[])
        uf == NoUnits ? (ax.xlabel = "Frequency") : (ax.xlabel = "Frequency ($uf)")
    end
    if isempty(ax.ylabel[])
        ux == NoUnits ? (ax.ylabel = "Spectral density") :
        (ax.ylabel = "Spectral density ($ux)")
    end
    p = spectrumplot!(ax, DimArray(xmed[idxs], (Freq(f[idxs]))); peaks, kwargs...)
    color = isnothing(bandcolor) ? (p.color[], 0.5) : bandcolor
    lineattrs = [:linewidth, :alpha, :linestyle, :linecap, :joinstyle]
    bandattrs = [k => v for (k, v) in kwargs if !(k ∈ lineattrs)]
    _p = Makie.band!(ax, f[idxs], σₗ[idxs], σᵤ[idxs]; transparency = true, bandattrs...,
                     color)
    Makie.translate!(_p, 0, 0, -1.0)
    p
end

function spectrumplot(x::AbstractSpectrum; peaks = false, kwargs...)
    (f = Figure();
     ax = Axis(f[1, 1]);
     p = spectrumplot!(ax, x; peaks, kwargs...);
     Makie.FigureAxisPlot(f,
                          ax,
                          p))
end
Makie.plot!(ax, x::AbstractSpectrum; kwargs...) = spectrumplot!(ax, x; kwargs...)
Makie.plot(x::AbstractSpectrum; kwargs...) = spectrumplot(x; kwargs...)
function Makie.plot!(ax, x::AbstractSpectrum{T, 2}; kwargs...) where {T}
    spectrumplot!(ax, x; kwargs...)
end
Makie.plot(x::AbstractSpectrum{T, 2}; kwargs...) where {T} = spectrumplot(x; kwargs...)

"""
    spectrumplot!(ax::Axis, x::AbstractVector, y::AbstractVector)
"""
# const ToolSpectrumPlot = SpectrumPlot{Tuple{<:UnivariateSpectrum}}
# argument_names(::Type{<: ToolSpectrumPlot}) = (:x,)

# function plot!(p::ToolSpectrumPlot)
#     x = collect(dims(p[:x], Freq))
#     y = collect(p[:x])
#     spectrumplot!(p, x, y)
#     p
# end

# """
#     plotLFPspectra(X::UnivariateTimeSeries; slope=nothing, position=Point2f([5, 1e-5]), fs=nothing, N=1000, slopecolor=:crimson, kwargs...)

# Create a line frequency power (LFP) spectra plot for the given time series `X`.

# # Arguments
# - `slope`: The power-law slope. Default is `nothing`.
# - `position`: The position of the slope label. Default is `Point2f([5, 1e-5])`.
# - `fs`: The sampling frequency. Default is `nothing`.
# - `N`: The number of frequency bins. Default is `1000`.
# - `slopecolor`: The color of the slope line. Default is `:crimson`.
# - `kwargs...`: Additional keyword arguments to be passed to the plot.
# """
# function plotLFPspectra(X::UnivariateTimeSeries; slope=nothing, position=Point2f([5, 1e-5]), fs=nothing, N=1000, slopecolor=:crimson, kwargs...)
#     times = collect(dims(X, Ti))
#     if isnothing(fs)
#         Δt = times[2] - times[1]
#         all(Δt .≈ diff(times)) || @warn "Violated assumption: all(Δt .≈ diff(times))"
#     else
#         Δt = 1/fs
#     end

#     P = [fp(Array(x)) for x ∈ eachcol(X)]
#     𝑓 = P[1].freq # Should be pretty much the same for all columns?
#     psd = hcat([p.power for p ∈ P]...)
#     psd = psd./(sum(psd, dims=1).*(𝑓[2] - 𝑓[1]))
#     psd = DimArray(psd, (Freq(𝑓), dims(X, :channel)))
#     fig = traces(𝑓, Array(psd); xlabel="𝑓 (Hz)", ylabel="Ŝ", title="Normalised power spectral density", smooth=1, yscale=log10, doaxis=false, domean=false, yminorgridvisible=false, kwargs...)
#     if !isnothing(slope)
#         _psd = psd[Freq(DD.Between(slope...))]
#         c, r, f = powerlawfit(_psd)
#         lines!(LinRange(slope..., 100), f(LinRange(slope..., 100)), color=slopecolor, linewidth=5)
#         text!(L"$\alpha$= %$(round(r, sigdigits=2))", position=Point2f0(position), fontsize=40)
#     end
#     return fig
# end

# ? -------------------------- Colored trajectory -------------------------- ? #
@recipe(Trajectory, x, y, z) do scene
    Attributes(colormode = :velocity,
               linewidth = theme(scene, :linewidth),
               alpha = 0.8,
               colormap = :viridis,
               colorscale = identity,
               linestyle = theme(scene, :linestyle),
               linecap = theme(scene, :linecap),
               joinstyle = theme(scene, :joinstyle))
end

function Makie.plot!(plot::Trajectory)
    # x = lift((args...) -> [y for y in args], plot.input_args...)
    x = @lift ($(plot.x), $(plot.y), $(plot.z))

    f = x -> isfinite.(x)
    i = @lift reduce(.&, f.($(x)))
    z = @lift [y[$(i)] for y in $(x)]

    colormode = plot.colormode[]
    if colormode === :velocity
        dx = @lift [y[2:end] .- y[1:(end - 1)] for y in $(z)]
        sqr = x -> x .^ 2
        colors = @lift sqrt.(sum(sqr.($(dx))))
    elseif colormode === :time
        colors = @lift 1:length($(z)[1])
    elseif !isnothing(colormode) && colormode != :none
        error("Not a supported `colormode`")
    end
    _z = @lift [y[1:(end - 1)] for y in $(z)]
    lineattrs = [
        :colormap,
        :linewidth,
        :alpha,
        :colorscale,
        :linestyle,
        :linecap,
        :joinstyle
    ]
    lines!(plot, _z[]...; color = colors, [a => plot.attributes[a] for a in lineattrs]...)
    plot
end

# ? -------------------------- Trajectory shadows -------------------------- ? #
function shadows!(ax, x, y, z; shadowmode = :projection, swapshadows = false, kwargs...)
    (x isa Observable) || (x = Observable(x))
    (y isa Observable) || (y = Observable(y))
    (z isa Observable) || (z = Observable(z))
    i = @lift isfinite.($(x)) .* isfinite.($(y)) .* isfinite.($(z))
    x = @lift ($(x)[$(i)])
    y = @lift ($(y)[$(i)])
    z = @lift ($(z)[$(i)])

    limits = ax.finallimits
    _limits = limits[]
    len = @lift length($(x))

    if swapshadows
        xp = @lift fill($(limits).origin[1], $(len))
        yp = @lift fill($(limits).origin[2], $(len))
        zp = @lift fill($(limits).origin[3], $(len))
    else
        xp = @lift fill($(limits).origin[1] .+ $(limits).widths[1], $(len))
        yp = @lift fill($(limits).origin[2] .+ $(limits).widths[2], $(len))
        zp = @lift fill($(limits).origin[3], $(len))
    end

    if shadowmode === :projection
        lines!(ax, xp, y, z; kwargs...)
        lines!(ax, x, yp, z; kwargs...)
        lines!(ax, x, y, zp; kwargs...)
    end
    ax.finallimits[] = _limits
    return ax
end

function Makie.convert_arguments(::Type{<:Trajectory},
                                 x::MultivariateTimeSeries)
    x |> eachcol .|> collect |> Tuple
end

# ? ------------------------------- # Traces ------------------------------- ? #
@recipe(Traces, x, y, z) do scene
    Attributes(colormap = nothing,
               normalize = false, # Can be any normalization type from Normalizations.jl
               colorrange = nothing,
               linewidth = theme(scene, :linewidth),
               alpha = 0.8,
               colorscale = identity,
               linestyle = theme(scene, :linestyle),
               linecap = theme(scene, :linecap),
               joinstyle = theme(scene, :joinstyle))
end

function Makie.convert_arguments(::Type{<:TimeseriesTools.Traces},
                                 x::AbstractDimArray)
    (x.dims[1].val.data, x.dims[2].val.data, x.data)
end

function Makie.plot!(plot::Traces)
    x, y, z = plot.x, plot.y, plot.z
    colormap = plot.colormap
    normalize = plot.normalize[]
    z = lift(z) do z
        (normalize == true) && (normalize = Normalization.MinMax)
        if normalize != false && normalize <: Normalization.AbstractNormalization
            N = fit(normalize, z; dims = 1)
            z = Normalization.normalize(z, N)
        end
        z
    end

    if colormap[] isa Symbol
        colormap = @lift Makie.cgrad($colormap)
    end
    if colormap[] isa Makie.ColorGradient # Color traces by the y value
        if isnothing(plot.colorrange[])
            plot.colorrange = @lift (minimum($(y)), maximum($(y)))
        end
        _y = lift((x, r) -> (x .- r[1]) ./ (r[2] - r[1]), y, plot.colorrange)
        colormap = lift((c, i) -> [c[Float64(_i)] for _i in i], colormap, _y)
    end
    lineattrs = [
        :linewidth,
        :alpha,
        :colorscale,
        :linestyle,
        :linecap,
        :joinstyle
    ]
    for i in axes(z[], 2)
        _z = lift(x -> x[:, i], z)
        if isnothing(colormap[])
            lines!(plot, x, _z; [a => plot.attributes[a] for a in lineattrs]...)
        else
            color = lift(c -> c[i], colormap)
            lines!(plot, x, _z; color, [a => plot.attributes[a] for a in lineattrs]...)
        end
    end
    plot
end

function traces!(ax, S::MultivariateSpectrum; kwargs...)
    x, y, z = decompose(S)
    xu, cu, yu = (x, y, z) .|> eltype .|> unit
    xu = xu == NoUnits ? "" : "($xu)"
    cu = cu == NoUnits ? "" : "($cu)"
    yu = yu == NoUnits ? "" : "($yu)"
    isempty(ax.xlabel[]) && (ax.xlabel = "Frequency $xu")
    isempty(ax.ylabel[]) && (ax.ylabel = "Power $yu")
    traces!(ax, ustripall.(x), ustripall.(y), ustripall.(z); kwargs...)
end

function traces!(ax, S::MultivariateTimeSeries; kwargs...)
    x, y, z = decompose(S)
    xu, cu, yu = (x, y, z) .|> eltype .|> unit
    xu = xu == NoUnits ? "" : "($xu)"
    cu = cu == NoUnits ? "" : "($cu)"
    yu = yu == NoUnits ? "" : "($yu)"
    isempty(ax.xlabel[]) && (ax.xlabel = "Time $xu")
    isempty(ax.ylabel[]) && (ax.ylabel = "Value $yu")
    traces!(ax, ustripall.(x), ustripall.(y), ustripall.(z); kwargs...)
end
function traces(S::MultivariateTimeSeries; figure = (;), axis = (;), kwargs...)
    f = Figure(; figure...)
    ax = Axis(f[1, 1]; axis...)
    p = traces!(ax, S; kwargs...)
    return Makie.FigureAxisPlot(f, ax, p)
end

MVIrregular = AbstractDimMatrix{T,
                                <:Tuple{A, B}} where {A <: TimeDim{<:RegularIndex},
                                                      B <: Dimension{<:RegularIndex}, T}
function Makie.plot(x::MVIrregular; kwargs...)
    traces(x; kwargs...)
end
function Makie.plot!(ax, x::MVIrregular; kwargs...)
    traces!(ax, x; kwargs...)
end
function Makie.plot(x::MultivariateTimeSeries; kwargs...)
    Makie.heatmap(x; kwargs...)
end
function Makie.plot!(ax, x::MultivariateTimeSeries; kwargs...)
    Makie.heatmap!(ax, x; kwargs...)
end

# ? --------------------------- # Stacked traces --------------------------- ? #
@recipe(StackedTraces, x, y, z) do scene
    Attributes(offset = 1,
               normalize = false,
               spacing = :close,
               linewidth = theme(scene, :linewidth),
               alpha = 0.8,
               colormap = :viridis,
               colorscale = identity,
               linestyle = theme(scene, :linestyle),
               linecap = theme(scene, :linecap),
               joinstyle = theme(scene, :joinstyle))
end

function Makie.convert_arguments(::Type{<:Plot{TimeseriesTools.StackedTraces}},
                                 x::AbstractDimArray)
    (x.dims[1].val.data, x.dims[2].val.data, x.data)
end

function Makie.plot!(plot::StackedTraces)
    x, y, z = plot.x, plot.y, plot.z
    offset = plot.offset
    normalize = plot.normalize[]
    spacing = plot.spacing
    z = lift(z) do z
        (normalize == true) && (normalize = Normalization.MinMax)
        if normalize != false && normalize <: Normalization.AbstractNormalization
            N = fit(normalize, z; dims = 1)
            z = Normalization.normalize(z, N)
        end
        z
    end

    z = lift(offset, spacing, z) do offset, spacing, z
        if offset == false
            offset = 0
        end
        c = zeros(size(z, 2))
        if spacing === :even
            space = maximum([maximum(z[:, i - 1] .- z[:, i]) for i in axes(z, 2)[2:end]])
        end
        for i in axes(z, 2)[2:end]
            if spacing === :close
                space = maximum(z[:, i - 1] .- z[:, i])
            end
            c[i] = c[i - 1] + space * offset
        end
        z .+ c'
    end
    plot.attributes.normalize[] = false
    lineattrs = [
        :colormap,
        :linewidth,
        :alpha,
        :colorscale,
        :linestyle,
        :linecap,
        :joinstyle
    ]
    traces!(plot, x, y, z; normalize = plot.attributes[:normalize],
            [a => plot.attributes[a] for a in lineattrs]...)
    plot
end

function stackedtraces!(ax, S::MultivariateSpectrum; kwargs...)
    x, y, z = decompose(S)
    xu, cu, yu = (x, y, z) .|> eltype .|> unit
    xu = xu == NoUnits ? "" : "($xu)"
    cu = cu == NoUnits ? "" : "($cu)"
    yu = yu == NoUnits ? "" : "($yu)"
    isempty(ax.xlabel[]) && (ax.xlabel = "Frequency $xu")
    isempty(ax.ylabel[]) && (ax.ylabel = "Power $yu")
    stackedtraces!(ax, ustripall.(x), ustripall.(y), ustripall.(z); kwargs...)
end

function stackedtraces!(ax, S::MultivariateTimeSeries; kwargs...)
    x, y, z = decompose(S)
    xu, cu, yu = (x, y, z) .|> eltype .|> unit
    xu = xu == NoUnits ? "" : "($xu)"
    cu = cu == NoUnits ? "" : "($cu)"
    yu = yu == NoUnits ? "" : "($yu)"
    isempty(ax.xlabel[]) && (ax.xlabel = "Time $xu")
    isempty(ax.ylabel[]) && (ax.ylabel = "Value $yu")
    stackedtraces!(ax, ustripall.(x), ustripall.(y), ustripall.(z); kwargs...)
end

# end # module
