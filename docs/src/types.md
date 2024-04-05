```@index
Pages = ["types.md"]
```
The `TimeseriesTools` package provides a standardized collection of types designed for handling various types of time-series data.
Defining consistent types for time series provides three key benefits:
1. It simplifies the workspace and method signatures by aggregating much of the data that defines a time series into a single variable; thanks to the [`DimensionalData.jl`](https://github.com/rafaqz/DimensionalData.jl) package, one variable can hold the values of a time series, its time points, spatial coordinates, units, metadata, and more.
2. It facilitates the development of generic functions that dispatch on the various types of time series; for instance, more efficient algorithms can be written for [`RegularTimeSeries`](@ref) types than for [`IrregularTimeSeries`](@ref) types, but the same high-level functionality can be provided by the same generic function that dispatches these methods given the type of the input time series.
3. Most importantly, this intuitively aligns the structure of time-series data in code to mathematical conventions, which can vastly simplify the process of developing and interpreting programs. Many small complexities (Was this time series regularly sampled? What are the output frequencies of my Fourier transform? The units of my power spectrum) are handled effortlessly, leaving room to focus on higher-level problems.

Below are a list of types defined in this package.


```@autodocs
Modules = [TimeseriesTools]
Pages   = ["Types.jl", "Spectrograms.jl"]
```
