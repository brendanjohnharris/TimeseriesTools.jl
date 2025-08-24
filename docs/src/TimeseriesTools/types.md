The `TimeseriesTools` package provides a standardized collection of types designed for handling various types of time-series data.
Defining consistent types for time series provides three key benefits:
1. It simplifies the workspace and method signatures by aggregating much of the data that defines a time series into a single variable; thanks to the [`DimensionalData.jl`](https://github.com/rafaqz/DimensionalData.jl) package, one variable can hold the values of a time series, its time points, spatial coordinates, units, metadata, and more.
2. It facilitates generic functions that dispatch on the various types of time series; for instance, more efficient algorithms can be written for [`RegularTimeseries`](@ref) types than for [`IrregularTimeseries`](@ref) types, but the same high-level functionality can be provided by the same generic function that dispatches these methods given the type of the input time series.
3. Most importantly, this intuitively aligns the structure of time-series data in code to mathematical conventions, simplifying the process of developing and interpreting programs. Many small complexities (Was this time series regularly sampled? What are the output frequencies of my Fourier transform? The units of my power spectrum?) are handled automatically, leaving room to focus on higher-level problems.

To achieve this, TimeseriesTools.jl defines a custom version of the `DimensionalData.DimArray` and custom `DimensionalData.Dimension`s:
```julia
x = Timeseries(rand(10), 1:10)
x isa AbstractToolsArray # In most cases, an AbstractToolsArray behaves like a DimArray; see DimensionalData
x isa AbstractTimeseries # An AbstractTimeseries is an AbstractToolsArray...
lookup(x, 1) isa 𝑡 # ...where the first dimension is a custom TimeDim 𝑡
```
If a `ToolsArray` or `DimArray` has a `𝑡` as its first dimension, it will be rebuilt as a `ToolsArray` (i.e. when using functions like `eachcol`).
There are a small number of other custom dimensions, all exported, that share this property and are subtypes of `ToolsDimension`: e.g. `𝑥`, `𝑦`, `𝑧`, `𝑓`,`Var`, `Obs`.
To define more of these `ToolsDimension`s, use:
```julia
DimensionalData.@dim NewDim ToolsDim "NameOfNewDim"
```
Please note that functions operating on a `ToolsArray` without a `ToolsDimension` as the first or last dimension may NOT return a `ToolsArray`, especially if they perform slicing and rebuilding. Be careful using the `DimensionalData.Dim{:name}` syntax.

Below is a full list of types defined in this package.

```@autodocs
Modules = [TimeseriesTools]
Pages   = ["Types.jl", "Spectrograms.jl"]
```
