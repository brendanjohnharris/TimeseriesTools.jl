# TimeseriesTools 🕰️🛠️

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://brendanjohnharris.github.io/TimeseriesTools.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://brendanjohnharris.github.io/TimeseriesTools.jl/dev/)
[![Build Status](https://github.com/brendanjohnharris/TimeseriesTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brendanjohnharris/TimeseriesTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/brendanjohnharris/TimeseriesTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/brendanjohnharris/TimeseriesTools.jl)


TimeseriesTools.jl is a sleek and powerful package for analyzing and visualizing time-series data in Julia. It provides a set of functions for preprocessing, analyzing, and plotting time series data, making your life better and your data look great (in that order)!

## Features

- 📈 Practical utilities for working with time series
- 📊 Spectral analysis and visualization
- 🌈 Beautiful plotting using [Makie](https://github.com/MakieOrg/Makie.jl)

## Installation

To install TimeseriesTools.jl, simply run the following command in your Julia REPL:

```julia
] add TimeseriesTools
```

## Usage

Here's a quick example to get you started:

```julia
using TimeseriesTools, CairoMakie, TimeseriesTools.FFTW, Unitful
import TimeseriesTools.TimeSeries # or TS

# Generate some quick brown noise
t = 0.005:0.005:100.0
x = colorednoise(t)*u"V"

# Plot the time series
lines(x[1:10000])

# Calculate the power spectrum
S = powerspectrum(x)
@test Pxx isa RegularSpectrum

# Plotting
p = @test_nowarn lines(Pxx)
```

## Example Spectrum Plot

```julia
# To create a beautiful spectrum plot, use the following code:

using TimeseriesTools, CairoMakie

# Load your time series data
data = ... # Load your time series data here

# Compute the energy spectrum
energy_spectrum = _energyspectrum(data)

# Create a SpectrumPlot
spectrum_plot = SpectrumPlot(energy_spectrum)

# Visualize the plot
figure, axis, lines = plot(spectrum_plot)
```

![Example Spectrum Plot](path/to/your/example_spectrum_plot.png)

## Acknowledgements 🙏

TimeseriesTools.jl builds upon the excellent [DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) package for handling dimensions and indexing in time series data.

Happy analyzing! 🚀