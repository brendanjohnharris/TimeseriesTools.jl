# Changelog

All notable changes to TimeseriesTools.jl.


## v0.9.0 (2026-09-20)

Breaking:

- Surrogate methods moved to [TimeseriesToolsSurrogates.jl](https://www.github.com/brendanjohnharris/TimeseriesToolsSurrogates.jl); `GammaRenewal`, `NDFT`, `RandomJitter`, `nansubarray` and `phaserand!` are no longer exported, and `TimeseriesSurrogates` is no longer a dependency.
- `Distributions` moved behind a `DistributionsExt` extension; the `InverseFunctions` dependency was dropped.
- Removed `SignalDecompositionExt`.
- Minimum Julia version raised to 1.11.
- `NaturalNeighbours` interpolation now extends `TimeseriesTools.interpolate` rather than `NaturalNeighbours.interpolate`.
- New exported names (`interpolate`, `resample`, `impute`, `stderror`, and the MAPPLE accessors) may conflict with other packages under `using`.
- `fit_mapple(log_f, log_s; w)`: the peak-finder smoothing window is now `window`. `w` is reserved for the residual weights of the refinement, and `fit(MAPPLE, s; w = true)` forwards it there.
- `phasestitch` now includes its first segment. Previously the output began at the second, so results are longer by one segment; see Fixed below.
- `TimeseriesMakie` removed from `[weakdeps]` and `[compat]`. No extension used it, so it only constrained resolution: v0.7-0.8.0 capped `TimeseriesMakie` at 0.2, which made newer versions unresolvable alongside this package. Completes the decoupling begun in v0.8.0.

Added:

- Expanded MAPPLE interface: parameter accessors (`betas`, `breakpoints`, `peakfreqs`, `peaksigmas`, `peakamplitudes`, and Hz-reporting variants), `logweights`, BIC-based component selection, a `fix` keyword for held parameters, and `stderror`/`vcov`/`rsquared` through `OptimExt` (which now also requires `ForwardDiff`).
- Interpolation and resampling interface: `interpolate`, `upsample`, `downsample`, `resample` and `impute`, with joint N-dimensional interpolation through `DataInterpolationsNDExt`.

Fixed:

- `bandpass` under DSP 0.8.
- `phasestitch` dropped its first segment: the loop matched each segment against its predecessor but only ever collected the successor's tail, so the opening segment never reached the output.
- `phasestitch` threw `ArgumentError` ("reducing over an empty collection") when no segment phases matched within `tol`, rather than returning what it had. Reachable with short inputs, and intermittent, since matching depends on the data.
- `convolve` with a positive `range` returned a closure rather than `0.0` where no spike lay within range.
- `waveletspectrogram` on array views. ContinuousWavelets reads the compute device off the outermost array wrapper, so from v1.2.2 it rejects a view, reshape or adjoint against a dense daughter matrix; inputs are now materialised first, preserving the device. Requires ContinuousWavelets v1.2.2.
- `waveletspectrogram` on a `MultivariateTimeseries` now runs as a single batched transform rather than a transform per column, and no longer discards its positional and keyword arguments.
- `waveletspectrogram` on GPU data no longer copies to the host. The bundled `cwt(::CuArray, ...)` method was written against ContinuousWavelets v1.1 internals and errors under v1.2.2; it has been removed in favour of upstream's own GPU support.
- `CUDAExt` no longer requires `ContinuousWavelets` to be loaded, and `ContinuousWaveletsExt` no longer requires `Mmap`.

## v0.8.1 (2026-03-08)

- DimensionalData v0.30 compatibility (method-ambiguity and `ustrip` fixes).
- CI: scheduled and dispatch triggers; runs on all branches.
- Test fixes: seeded Optim test, `Random` added to test deps.

## v0.8.0 (2025-09-02)

- Folded `TimeseriesBase` into the package; removed standalone docs site.
- Stopped re-exporting `TimeseriesMakie` to break a circular dependency.
- Aqua fixes; updated `Normalization` compat.

## v0.7.1 / v0.7.0 (2025-08-24)

- Release line for `TimeseriesMakie` split-out; Makie dependency removed from core.
- Finalized `TimeseriesPlots`/`TimeseriesMakie` compatibility.
- Dropped greedy Unitful exponents.

## v0.6.3 (2025-04-29)

- Multivariate `findpeaks` fix.
- `DataInterpolations` compat; relaxed gamma-renewal test.

## v0.6.2 (2025-01-19)

- `Obs` (Observation) dimension added.
- `spikeraster` plot with rate sorting; multiple spectrum-plot limit fixes.
- `ProgressMap` scheduler options with version guards; Dagger backend behind an extension.

## v0.6.1 (2024-12-10)

- DOI added; dummy `timescale` method; minor fixes.

## v0.6.0 (2024-11-15)

- `DataInterpolations` replaces `Dierckx` for interpolation.
- `ProgressMap` gains multiple backends (default `ProgressLogging`, optional `Dagger` extension).
- Shadow trajectories, generalized phase, and `ustrip` fixes.
- Switched to `LTS` for tests.

## v0.5.x (2024-08 to 2024-11)

- v0.5.4: Coarse-graining fix; `TDim` dimension type (DimensionalData-style).
- v0.5.3: DimensionalData v0.29 bump; removed `SSSet`.
- v0.5.2/v0.5.1: Removed `Requires`-based loading in favour of package extensions; printing/test cleanups.
- v0.5.0: Broad refactor for DimensionalData v0.28 (breaking); test-suite overhaul, IO and `ustripall` fixes; `TimeseriesSurrogates` temporarily moved into deps.

## v0.4.0 (2024-08-20)

- Replaced `DimArray` with a custom `ToolsArray` type; new type system across the package.
- Mean-squared displacement, partition tests, `maskpeaks`.
- Makie v0.21 compat; multidimensional FT surrogates; bandpass over columns.
- `tsv` for saving time series.

## v0.3.0 (2024-01-27)

- Central derivatives for irregular time series.
- Coarse-graining for arrays; generalized phase; unitful interpolation fixes.
- `findpeaks`; dimension matching; alignment.
- Switched interpolation to `Dierckx`.

## v0.2.x (2023-06 to 2023-11)

- v0.2.5/v0.2.4: GPU test updates; spike-train spectra via autocovariance; peak detection on power-spectrum plot; instantaneous frequency.
- v0.2.3: Wavelet transform; CUDA extension for wavelets; spike-train surrogates; analytic phase/amplitude; windowing functions; spike-time tiling coefficient and covariance.
- v0.2.2/v0.2.1: Spike-train types and power spectra; buffer; more flexible `TimeSeries` construction.
- v0.2.0: DSP extension with bandpass methods and phasestitch; preliminary `DateTime` support; traces recipe; unit-power normalization; updated save/load.

## v0.1.0 (2023-04-16)

- Initial release: time-series types, energy/power spectra, Unitful integration, Normalization, plotting recipes, docs and tests.

## Initial commit (2023-04-13)

- Bootstrapped from PkgTemplates.
