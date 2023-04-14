using Makie

# ..........decompose..........


function plotLFPspectra(LFP::AbstractDimArray; slope=nothing, position=Point2f([5, 1e-5]), fs=nothing, N=1000, slopecolor=:crimson, kwargs...)
    times = collect(dims(LFP, Ti))
    if isnothing(fs)
        Δt = times[2] - times[1]
        all(Δt .≈ diff(times)) || @warn "Violated assumption: all(Δt .≈ diff(times))"
    else
        Δt = 1/fs
    end

    P = [fp(Array(x)) for x ∈ eachcol(LFP)]
    𝑓 = P[1].freq # Should be pretty much the same for all columns?
    psd = hcat([p.power for p ∈ P]...)
    psd = psd./(sum(psd, dims=1).*(𝑓[2] - 𝑓[1]))
    psd = DimArray(psd, (Dim{:frequency}(𝑓), dims(LFP, :channel)))
    fig = traces(𝑓, Array(psd); xlabel="𝑓 (Hz)", ylabel="Ŝ", title="Normalised power spectral density", smooth=1, yscale=Makie.log10, doaxis=false, domean=false, yminorgridvisible=false, kwargs...)
    if !isnothing(slope)
        _psd = psd[Dim{:frequency}(DD.Between(slope...))]
        c, r, f = powerlawfit(_psd)
        lines!(LinRange(slope..., 100), f(LinRange(slope..., 100)), color=slopecolor, linewidth=5)
        text!(L"$\alpha$= %$(round(r, sigdigits=2))", position=Point2f0(position), fontsize=40)
    end
    return fig
end
