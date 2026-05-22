# Changelog

All notable changes to TimeseriesTools.jl.


## v0.8.1 — 2026-03-08

- DimensionalData v0.30 compatibility (method-ambiguity and `ustrip` fixes).
- CI: scheduled and dispatch triggers; runs on all branches.
- Test fixes: seeded Optim test, `Random` added to test deps.

## v0.8.0 — 2025-09-02

- Folded `TimeseriesBase` into the package; removed standalone docs site.
- Stopped re-exporting `TimeseriesMakie` to break a circular dependency.
- Aqua fixes; updated `Normalization` compat.

## v0.7.1 / v0.7.0 — 2025-08-24

- Release line for `TimeseriesMakie` split-out; Makie dependency removed from core.
- Finalized `TimeseriesPlots`/`TimeseriesMakie` compatibility.
- Dropped greedy Unitful exponents.

## v0.6.3 — 2025-04-29

- Multivariate `findpeaks` fix.
- `DataInterpolations` compat; relaxed gamma-renewal test.

## v0.6.2 — 2025-01-19

- `Obs` (Observation) dimension added.
- `spikeraster` plot with rate sorting; multiple spectrum-plot limit fixes.
- `ProgressMap` scheduler options with version guards; Dagger backend behind an extension.

## v0.6.1 — 2024-12-10

- DOI added; dummy `timescale` method; minor fixes.

## v0.6.0 — 2024-11-15

- `DataInterpolations` replaces `Dierckx` for interpolation.
- `ProgressMap` gains multiple backends (default `ProgressLogging`, optional `Dagger` extension).
- Shadow trajectories, generalized phase, and `ustrip` fixes.
- Switched to `LTS` for tests.

## v0.5.x — 2024-08 → 2024-11

- v0.5.4: Coarse-graining fix; `TDim` dimension type (DimensionalData-style).
- v0.5.3: DimensionalData v0.29 bump; removed `SSSet`.
- v0.5.2/v0.5.1: Removed `Requires`-based loading in favour of package extensions; printing/test cleanups.
- v0.5.0: Broad refactor for DimensionalData v0.28 (breaking); test-suite overhaul, IO and `ustripall` fixes; `TimeseriesSurrogates` temporarily moved into deps.

## v0.4.0 — 2024-08-20

- Replaced `DimArray` with a custom `ToolsArray` type; new type system across the package.
- Mean-squared displacement, partition tests, `maskpeaks`.
- Makie v0.21 compat; multidimensional FT surrogates; bandpass over columns.
- `tsv` for saving time series.

## v0.3.0 — 2024-01-27

- Central derivatives for irregular time series.
- Coarse-graining for arrays; generalized phase; unitful interpolation fixes.
- `findpeaks`; dimension matching; alignment.
- Switched interpolation to `Dierckx`.

## v0.2.x — 2023-06 → 2023-11

- v0.2.5/v0.2.4: GPU test updates; spike-train spectra via autocovariance; peak detection on power-spectrum plot; instantaneous frequency.
- v0.2.3: Wavelet transform; CUDA extension for wavelets; spike-train surrogates; analytic phase/amplitude; windowing functions; spike-time tiling coefficient and covariance.
- v0.2.2/v0.2.1: Spike-train types and power spectra; buffer; more flexible `TimeSeries` construction.
- v0.2.0: DSP extension with bandpass methods and phasestitch; preliminary `DateTime` support; traces recipe; unit-power normalization; updated save/load.

## v0.1.0 — 2023-04-16

- Initial release: time-series types, energy/power spectra, Unitful integration, Normalization, plotting recipes, docs and tests.

## 2023-04-13 — Initial commit

- Bootstrapped from PkgTemplates.
