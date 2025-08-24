```@meta
CurrentModule = TimeseriesMakie
```

# TimeseriesMakie

[TimeseriesMakie.jl](https://github.com/brendanjohnharris/TimeseriesMakie.jl) provides plotting recipes and visualization tools for time series data, using [Makie.jl](https://github.com/MakieOrg/Makie.jl).

## Installation

```julia
using Pkg
Pkg.add("TimeseriesMakie")
Pkg.add("CairoMakie") # Choose a Makie backend (CairoMakie, GLMakie, WGLMakie, etc.)
using TimeseriesMakie, CairoMakie
```

To improve the display of unitful exponents, you can enable fancy exponents with:
```julia
ENV["UNITFUL_FANCY_EXPONENTS"] = true
```

## What's in the box?

A full list of recipes exported by `TimeseriesMakie` can be found on the [Recipes](recipes.md) page.