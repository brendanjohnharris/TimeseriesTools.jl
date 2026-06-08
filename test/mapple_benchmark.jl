# Deterministic MAPPLE fit-quality benchmark.
#
# Run in the `test` environment:
#     include("test/mapple_benchmark.jl")
#     rows = run_mapple_benchmark()
#     summarise(rows)
#
# Each case fits a synthetic spectrum through the standard pipeline (rough `fit_mapple`
# then Optim refinement) and records: final loss relative to the loss at the true
# parameters (1.0 = optimal), component-slope recovery error, peak-count recovery,
# correlation, time and allocations. This is the bar each robustness change must beat;
# `run_mapple_benchmark(; fitkw...)` forwards keywords to the refinement so a feature can
# be toggled and compared against the baseline.

using TimeseriesTools, Optim, ForwardDiff, ComponentArrays, Random, Statistics

bench_comp(; log_f_stop, β) = ComponentArray(; log_f_stop = float(log_f_stop), β = float(β))
bench_peak(; log_f, log_σ, log_A) = ComponentArray(; log_f = float(log_f), log_σ = float(log_σ), log_A = float(log_A))

# Each case: a name, the true components/peaks, and the (ncomp, npeak) requested of the fit.
function benchmark_cases()
    return [
        (; name = "falling_2comp_0peak",
            comps = [bench_comp(; β = -1.0, log_f_stop = 1.5), bench_comp(; β = -3.0, log_f_stop = 5.0)],
            peaks = ComponentArray[], log_A = 1.0, tw = 0.15, ncomp = 2, npeak = 0),
        (; name = "rising_2comp_2peak_steep",
            comps = [bench_comp(; β = 2.0, log_f_stop = 1.5), bench_comp(; β = 4.0, log_f_stop = 5.0)],
            peaks = [bench_peak(; log_f = 0.7, log_σ = 0.1, log_A = 1.5), bench_peak(; log_f = 2.5, log_σ = 0.1, log_A = 0.5)],
            log_A = 1.0, tw = 0.15, ncomp = 2, npeak = 2),
        (; name = "falling_1comp_1peak_shallow",
            comps = [bench_comp(; β = -1.5, log_f_stop = 5.0)],
            peaks = [bench_peak(; log_f = 1.0, log_σ = 0.15, log_A = 1.0)],
            log_A = 1.0, tw = 0.1, ncomp = 1, npeak = 1),
        (; name = "falling_1comp_1peak_steep",
            comps = [bench_comp(; β = -4.0, log_f_stop = 5.0)],
            peaks = [bench_peak(; log_f = 1.5, log_σ = 0.1, log_A = 0.5)],
            log_A = 1.0, tw = 0.1, ncomp = 1, npeak = 1),
        (; name = "falling_3comp_0peak",
            comps = [bench_comp(; β = -0.5, log_f_stop = 1.0), bench_comp(; β = -2.0, log_f_stop = 2.5),
                bench_comp(; β = -4.0, log_f_stop = 5.0)],
            peaks = ComponentArray[], log_A = 1.0, tw = 0.1, ncomp = 3, npeak = 0),
    ]
end

function bench_truth(case)
    peaks = isempty(case.peaks) ? map(_ -> bench_peak(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0) : case.peaks
    return ComponentArray(; log_A = case.log_A, peaks, components = case.comps, transition_width = case.tw)
end

function run_one(case; noise = 0.1, n = 500, seed = 0, fitkw...)
    Random.seed!(seed)
    truth = bench_truth(case)
    log_f = range(0, 3, length = n)
    f = exp10.(log_f)
    s_clean = mapple(f, truth)
    log_s = log10.(s_clean) .+ noise .* randn(n)
    loss(p) = sum(abs2, log_s .- log10.(mapple(f, p)))
    init = fit_mapple(log_f, log_s; components = case.ncomp, peaks = case.npeak, w = 50)
    t = @timed refined = fit_mapple(log_f, log_s, init; autodiff = Optim.ADTypes.AutoForwardDiff(), fitkw...)
    βtrue = sort(collect(Float64, truth.components.β))
    βfit = sort(collect(Float64, refined.components.β))
    return (; name = case.name, noise,
        rel_loss = loss(refined) / max(loss(truth), eps()),
        loss_truth = loss(truth), loss_fit = loss(refined),
        cor = cor(mapple(f, refined), s_clean),
        beta_err = maximum(abs, βfit .- βtrue),
        npeak_true = case.npeak, npeak_fit = length(refined.peaks),
        time = t.time, bytes = t.bytes)
end

function run_mapple_benchmark(; noises = (0.05, 0.2), fitkw...)
    rows = NamedTuple[]
    for case in benchmark_cases(), nz in noises
        push!(rows, run_one(case; noise = nz, fitkw...))
    end
    return rows
end

getcol(rows, k) = getindex.(rows, k)

summarise(rows) = (;
    median_rel_loss = median(getcol(rows, :rel_loss)),
    worst_rel_loss = maximum(getcol(rows, :rel_loss)),
    n_degenerate = count(>(3.0), getcol(rows, :rel_loss)),   # rel_loss >> 1 == bad fit
    median_beta_err = median(getcol(rows, :beta_err)),
    worst_beta_err = maximum(getcol(rows, :beta_err)),
    median_time = median(getcol(rows, :time)),
)
