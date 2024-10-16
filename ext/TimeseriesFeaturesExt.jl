# module TimeseriesFeaturesExt
import ..TimeseriesFeatures: SPI
using TimeseriesTools
using TimeseriesTools.Statistics
using TimeseriesTools.Distributions
using TimeseriesTools.DimensionalData

STOIC = SPI((x, y) -> stoic(x, y), :STOIC, "Spike-time overlap integral coefficient",
            ["correlation"])

export STOIC

# end # module
