# TimeseriesTools

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://brendanjohnharris.github.io/TimeseriesTools.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://brendanjohnharris.github.io/TimeseriesTools.jl/dev/)
[![Build Status](https://github.com/brendanjohnharris/TimeseriesTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brendanjohnharris/TimeseriesTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/brendanjohnharris/TimeseriesTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/brendanjohnharris/TimeseriesTools.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14511321.svg)](https://doi.org/10.5281/zenodo.14511321)

TimeseriesTools.jl is a package for analyzing and visualizing time-series data in Julia.

## Features

- 📈 Practical utilities for working with time series
- 📊 Spectral analysis and visualization
- 🌈 Beautiful plotting using [Makie](https://github.com/MakieOrg/Makie.jl)

![Example Shadow Plot](test/shadows_dark.png#gh-dark-mode-only)
![Example Shadow Plot](test/shadows.png#gh-light-mode-only)

## Installation

To install TimeseriesTools.jl, simply run the following command in your Julia REPL:

```julia
] add TimeseriesTools
```

## Usage

Here's a quick example to get you started:

```julia
using TimeseriesTools, CairoMakie, TimeseriesTools.FFTW, Unitful

# Generate some quick brown noise
t = 0.005:0.005:1e5
x = colorednoise(t*u"s")*u"V" # ::AbstractTimeseries
plot(x[1:10000])
S = powerspectrum(x, 0.001)
p = plot(S)
```

![Example Time-series Plot](test/timeseries_dark.png#gh-dark-mode-only)
![Example Spectrum Plot](test/powerspectrum_dark.png#gh-dark-mode-only)
![Example Time-series Plot](test/timeseries.png#gh-light-mode-only)
![Example Spectrum Plot](test/powerspectrum.png#gh-light-mode-only)

Note that an instance of the most basic type of this package, the `AbstractTimeseries`, can be generated with:
```julia
t = 0:0.01:1
x = sin.(t)
Timeseries(x, t)
```
Please see the documentation for further functionality.

## Acknowledgements 🙏

TimeseriesTools.jl builds upon the excellent [DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) package for handling dimensions and indexing in time-series data.

Happy analyzing! 🚀