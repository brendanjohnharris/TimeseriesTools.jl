using Pkg
Pkg.activate(tempname())
Pkg.add("TimeseriesTools")
Pkg.add("DifferentialEquations")
Pkg.add(name = "NonstationaryProcessesBase", rev = "main")
Pkg.add(url = "https://www.github.com/brendanjohnharris/NonstationaryProcesses.jl",
        rev = "main")
import TimeseriesTools.TimeSeries
using TimeseriesTools
using DifferentialEquations
using NonstationaryProcesses
using NonstationaryProcesses.DifferentialEquationsExt

transient = 20000
linewidth = 0.05

lorenz = lorenzSim(X0 = [0.0, -0.01, 9.0],
                   parameter_profile = (constantParameter, constantParameter,
                                        constantParameter),
                   parameter_profile_parameters = (10.0, 28.0, 8 / 3), # Sprott's recomendation
                   transient_t0 = -100.0,
                   t0 = 0.0,
                   dt = 0.001,
                   savedt = 0.025,
                   tmax = 1500.0,
                   alg = AutoVern7(Rodas5()),
                   solver_opts = Dict(:adaptive => true, :reltol => 1e-12))

sx = Timeseries(lorenz)

x = sx[transient:(end - transient), :]

file = joinpath(@__DIR__, @__FILE__)[1:(end - 3)] * ".csv"
savetimeseries(file, x)
y = loadtimeseries(file)
@assert all(y .== x)
