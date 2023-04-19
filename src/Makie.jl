using ..Makie
import MakieCore.plot!
import MakieCore.plot
Base.iterate(s::Makie.RichText, i::Integer) = iterate(String(s), i)
Base.iterate(s::Makie.RichText) = iterate(String(s))
"""
    spectrumplot!(ax::Axis, x::UnivariateSpectrum)
Plot the given spectrum, labelling the axes, adding units if appropriate, ribbons if the input is a [`MultivariateSpectrum`](@ref), and other niceties.
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
spectrumplot(x::UnivariateSpectrum; kwargs...) = (f=Figure(); ax=Axis(f[1, 1]); p=spectrumplot!(ax, x; kwargs...); Makie.FigureAxisPlot(f, ax, p))
plot!(ax, x::UnivariateSpectrum; kwargs...) = spectrumplot!(ax, x; kwargs...)
plot(x::UnivariateSpectrum; kwargs...) = spectrumplot(x; kwargs...)


function plot!(ax::Makie.Axis, x::UnivariateTimeSeries; kwargs...)
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
plot(x::UnivariateTimeSeries; kwargs...) = (f=Figure(); ax=Axis(f[1, 1]); p=plot!(ax, x; kwargs...); Makie.FigureAxisPlot(f, ax, p))




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
