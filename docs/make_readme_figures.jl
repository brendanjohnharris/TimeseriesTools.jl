#! /bin/bash
#=
exec julia +1.12 -t auto --project="$(dirname "${BASH_SOURCE[0]}")" "${BASH_SOURCE[0]}" "$@"
=#

# Regenerates the figures embedded in ../README.md, into ./plots.
#
# These used to live in the test suite, which made TimeseriesMakie, CairoMakie and Fathom test
# dependencies of this package purely to draw pictures. Run this instead, either directly or by
# activating `docs/` and including it.

using TimeseriesTools
using TimeseriesMakie
using CairoMakie
using Fathom
using Unitful

const PLOTS = joinpath(@__DIR__, "plots")
const TRAJECTORY = joinpath(@__DIR__, "test_timeseries.tsv")

"""
    readmefigures(suffix, shadowcolor)

Draw the three README figures into `docs/plots`, suffixed for the theme in force.
"""
function readmefigures(suffix, shadowcolor)
    t = 0.005:0.005:1.0e5
    x = colorednoise(t * u"s") * u"V"

    f = Figure(; size = (720, 480))
    ax = Axis(f[1, 1])
    plot!(ax, x[1:10000])
    save(joinpath(PLOTS, "timeseries$suffix.png"), f; px_per_unit = 3)

    S = _powerspectrum(x, 0.001)
    f = Figure(; size = (720, 480))
    ax = Axis(f[1, 1])
    plotspectrum!(ax, S; linewidth = 1)
    save(joinpath(PLOTS, "powerspectrum$suffix.png"), f; px_per_unit = 3)

    y = loadtimeseries(TRAJECTORY)
    f = Figure(; size = (500, 480))
    ax = Axis3(f[1, 1])
    trajectory!(ax, collect.(eachcol(y))...; colormap = :turbo, linewidth = 0.1, color = :speed)
    shadows!(
        ax, collect.(eachcol(y))...; color = (shadowcolor, 0.5), linewidth = 0.05,
        swapshadows = (true, false, false),
        limits = Fathom.widen.(extrema.(eachcol(y)), 0.4)
    )
    ax.xlabelvisible = ax.ylabelvisible = ax.zlabelvisible = false
    ax.xticksvisible = ax.yticksvisible = ax.zticksvisible = false
    ax.xticklabelsvisible = ax.yticklabelsvisible = ax.zticklabelsvisible = false
    ax.azimuth[] = 2.2
    ax.elevation[] = 0.5
    hidespines!(ax)
    save(joinpath(PLOTS, "shadows$suffix.png"), f; px_per_unit = 3)
    return nothing
end

begin
    mkpath(PLOTS)

    set_theme!(fathom())
    readmefigures("", :slategray)

    set_theme!(fathom(:dark, :transparent))
    readmefigures("_dark", :white)

    set_theme!()
    @info "Wrote the README figures" PLOTS
end
