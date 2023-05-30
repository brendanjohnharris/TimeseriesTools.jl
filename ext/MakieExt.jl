# module MakieExt

import MakieCore.plot!
import MakieCore.plot

import ..Makie
import ..Makie: plot, plot!, lift, lines, lines!, band, band!, FigureAxisPlot, @lift, Observable, @recipe, Theme
using TimeseriesTools
using Statistics

export spectrumplot!, spectrumplot, trajectory!, trajectory, shadows!

Base.iterate(s::Makie.RichText, i::Integer) = iterate(String(s), i)
Base.iterate(s::Makie.RichText) = iterate(String(s))

"""
    spectrumplot!(ax::Axis, x::UnivariateSpectrum)
Plot the given spectrum, labelling the axes, adding units if appropriate, and other niceties.
"""
function spectrumplot!(ax::Makie.Axis, x::UnivariateSpectrum; kwargs...)
    uf = frequnit(x)
    ux = unit(x)
    f, x = decompose(x)
    f = ustrip.(f) |> collect
    x = ustrip.(x) |> collect
    idxs = (f .> 0) .& (x .> 0)
    ax.limits = ((minimum(f[idxs]), nothing), (minimum(x[idxs]), nothing));
    ax.xscale = log10
    ax.yscale = log10
    uf == NoUnits ? (ax.xlabel = "Frequency") : (ax.xlabel = "Frequency ($uf)")
    ux == NoUnits ? (ax.ylabel = "Spectral density") : (ax.ylabel = "Spectral density ($ux)")
    p = spectrumplot!(ax, f[idxs], x[idxs]; kwargs...)
    p
end

"""
    spectrumplot!(ax::Axis, x::MultivariateSpectrum)
Plot the given spectrum, labelling the axes, adding units if appropriate, and adding a band to show the iqr
"""
function spectrumplot!(ax::Makie.Axis, x::MultivariateSpectrum; kwargs...)
    uf = frequnit(x)
    ux = unit(x)
    f, _, x = decompose(x)
    f = ustrip.(f) |> collect
    x = ustrip.(x) |> collect
    xmin = minimum(x, dims=2) |> vec
    xmed = median(x, dims=2) |> vec
    σₗ = mapslices(x->quantile(x, 0.25), x, dims=2) |> vec
    σᵤ = mapslices(x->quantile(x, 0.75), x, dims=2) |> vec
    idxs = (f .> 0) .& (xmin .> 0)
    ax.limits = ((minimum(f[idxs]), nothing), (minimum(σₗ[idxs]), nothing));
    ax.xscale = log10
    ax.yscale = log10
    uf == NoUnits ? (ax.xlabel = "Frequency") : (ax.xlabel = "Frequency ($uf)")
    ux == NoUnits ? (ax.ylabel = "Spectral density") : (ax.ylabel = "Spectral density ($ux)")
    _p = Makie.band!(ax, f[idxs], σₗ[idxs], σᵤ[idxs]; transparency=0.2, kwargs...)
    p = spectrumplot!(ax, f[idxs], xmed[idxs]; kwargs...)
    p
end

spectrumplot(x::AbstractSpectrum; kwargs...) = (f=Figure(); ax=Axis(f[1, 1]); p=spectrumplot!(ax, x; kwargs...); Makie.FigureAxisPlot(f, ax, p))
Makie.plot!(ax, x::AbstractSpectrum; kwargs...) = spectrumplot!(ax, x; kwargs...)
Makie.plot(x::AbstractSpectrum; kwargs...) = spectrumplot(x; kwargs...)





function Makie.plot!(ax::Makie.Axis, x::UnivariateTimeSeries; kwargs...)
    ut = timeunit(x)
    ux = unit(x)
    t, x = decompose(x)
    t = ustrip.(t) |> collect
    x = ustrip.(x) |> collect
    p = lines!(ax, t, x; kwargs...)
    ut == NoUnits ? (ax.xlabel = "Time") : (ax.xlabel = "Time ($ut)")
    ux == NoUnits ? (ax.ylabel = "Values") : (ax.ylabel = "Values ($ux)")
    p
end
Makie.plot(x::UnivariateTimeSeries; kwargs...) = (f=Figure(); ax=Axis(f[1, 1]); p=plot!(ax, x; kwargs...); Makie.FigureAxisPlot(f, ax, p))




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
#     psd = DimArray(psd, (Dim{:frequency}(𝑓), dims(X, :channel)))
#     fig = traces(𝑓, Array(psd); xlabel="𝑓 (Hz)", ylabel="Ŝ", title="Normalised power spectral density", smooth=1, yscale=log10, doaxis=false, domean=false, yminorgridvisible=false, kwargs...)
#     if !isnothing(slope)
#         _psd = psd[Dim{:frequency}(DD.Between(slope...))]
#         c, r, f = powerlawfit(_psd)
#         lines!(LinRange(slope..., 100), f(LinRange(slope..., 100)), color=slopecolor, linewidth=5)
#         text!(L"$\alpha$= %$(round(r, sigdigits=2))", position=Point2f0(position), fontsize=40)
#     end
#     return fig
# end



# ? -------------------------- Colored trajectory -------------------------- ? #
@recipe(Trajectory, x, y, z) do scene
    Theme(
        colormode = :velocity,
    )
end

function Makie.plot!(plot::Trajectory)
    x = lift((args...) -> [y for y in args], plot.input_args...)
    f = x->isfinite.(x)
    i = @lift reduce(.&, f.($(x)))
    z = @lift [y[$(i)] for y ∈ $(x)]

    colormode = plot.colormode[]
    if colormode === :velocity
        dx = @lift [y[2:end] .- y[1:end-1] for y ∈ $(z)]
        sqr = x->x.^2
        colors = @lift sqrt.(sum(sqr.($(dx))))
    elseif colormode === :time
        colors = @lift 1:length($(z)[1])
    elseif !isnothing(colormode) && colormode != :none
        @error "Not a supported `colormode`"
    end
    _z = @lift [y[1:end-1] for y in $(z)]
    lines!(plot, _z[]...; plot.attributes..., color=colors)

    plot
end

# ? -------------------------- Trajectory shadows -------------------------- ? #
function shadows!(ax, x, y, z; shadowmode=:projection, swapshadows=false, kwargs...)
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



# end # module
