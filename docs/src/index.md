```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "TimeseriesTools"
  text: ""
  tagline:
  image:
    src: /logo.png
    alt: TimeseriesTools logo.
  actions:
    - theme: brand
      text: Introduction
      link: /quickstart
    - theme: alt
      text: Examples
      link: /examples
    - theme: alt
      text: View on Github
      link: https://github.com/brendanjohnharris/TimeseriesTools.jl

features:
  - icon: <img width="64" height="64" src="https://rawcdn.githack.com/JuliaLang/julia-logo-graphics/f3a09eb033b653970c5b8412e7755e3c7d78db9e/images/juliadots.iconset/icon_512x512.png" alt="TimeseriesTools"/>
    title: TimeseriesTools
    details: Time-series analysis
    link: /TimeseriesTools/timeseriestools
  - icon: <img width="64" height="64" src="https://rawcdn.githack.com/JuliaLang/julia-logo-graphics/f3a09eb033b653970c5b8412e7755e3c7d78db9e/images/juliadots.iconset/icon_512x512.png" />
    title: TimeseriesMakie
    details: Makie recipes for visualizing time-series
    link: /TimeseriesMakie/timeseriesmakie
---


<p style="margin-bottom:2cm"></p>

<div class="vp-doc" style="width:80%; margin:auto">

```

Welcome to the documentation for [TimeseriesTools](https://github.com/brendanjohnharris/TimeseriesTools.jl), a Julia package for working with time series data.

## Key Features

- **Unified Data Structures**: Built on top of DimensionalData.jl, providing consistent handling of time series, spectra, and multidimensional data
- **Intuitive Type System**: Specialized types for different kinds of time series (regular, irregular, multivariate) and spectra
- **Unitful Integration**: Full support for physical units via Unitful.jl, ensuring dimensional consistency across operations
- **Signal Processing**: Comprehensive tools for filtering, spectral analysis, and time-frequency analysis
- **Visualization**: Seamless integration with Makie.jl for plotting time series and spectra
- **Performance**: Optimized implementations with support for parallel processing


```@raw html
</div>
```
