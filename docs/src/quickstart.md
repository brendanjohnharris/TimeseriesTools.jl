## Installation

To install TimeseriesTools.jl, simply run the following command in your Julia REPL:
```julia
using Pkg
Pkg.add("TimeseriesTools")
```

Or:
```julia
] add TimeseriesTools
```

## Usage

An instance of the most basic type of this package, the `AbstractTimeseries`, can be generated with:

```julia
using TimeseriesTools, CairoMakie
t = 0:0.01:1
Timeseries(sin.(t), t)
plot(x)
```

All `TimeseriesTools` arrays are built on top of [DimensionalData.jl](https://github.com/JuliaDynamics/DimensionalData.jl) `AbstractDimArray`s.

Here's another example, calculating a power spectrum:

```julia
using TimeseriesTools, CairoMakie, Unitful

# Generate some quick brown noise
t = range(0, stop=1e4, step=0.005)
x = colorednoise(t, u"s")*u"V" # ::AbstractTimeseries

f = Figure()
uc = Makie.UnitfulConversion(u"s^-1"; units_in_label = false)
ax = Axis(f[1, 1]; dim1_conversion=uc)
S = powerspectrum(x, 0.001) # Second arguments sets frequency spacing
plotspectrum!(ax, S)
f
```