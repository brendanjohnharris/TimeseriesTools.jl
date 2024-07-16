module CUDAExt
using CUDA
using TimeseriesTools
using ContinuousWavelets
import ContinuousWavelets: cwt, ensureComplex, getNWavelets, prepSignalAndPlans,
                           boundaryType, reflect, analyticTransformReal!,
                           analyticTransformComplex!,
                           otherwiseTransform!

CUDA.CuArray(x::AbstractTimeSeries) = set(x, CuArray(x.data))

function Base.show(io::IO, mime,
                   X::ToolsArray{T, N, Tp, F, C} where {T, N, Tp, F, C <: CuArray})
    Base.show(X.data)
end

function cwt(Y::CuArray{T, N}, cWav::CWT, daughters, fftPlans = 1) where {N, T}
    @assert typeof(N) <: Integer
    @debug "Computing cwt on the GPU"
    daughters isa CuArray || (daughters = CuArray(daughters))
    # vectors behave a bit strangely, so we reshape them
    if N == 1
        Y = reshape(Y, (length(Y), 1))
    end
    n1 = size(Y, 1)

    _, nScales, _ = getNWavelets(n1, cWav)
    #construct time series to analyze, pad if necessary
    x = reflect(Y, boundaryType(cWav)()) #this function is defined below

    # check if the plans we were given are dummies or not
    x̂, fromPlan = prepSignalAndPlans(x, cWav, fftPlans)
    # If the vector isn't long enough to actually have any other scales, just
    # return the averaging
    if nScales <= 0 || size(daughters, 2) == 1
        daughters = daughters[:, 1:1]
        nScales = 1
    end

    if isAnalytic(cWav.waveType)
        OutType = ensureComplex(T)
    else
        OutType = T
    end

    wave = CUDA.zeros(OutType, size(x)..., nScales)  # result array
    # faster if we put the example index on the outside loop through all scales
    # and compute transform
    if isAnalytic(cWav.waveType)
        if eltype(x) <: Real
            analyticTransformReal!(wave, daughters, x̂, fromPlan, cWav.averagingType)
        else
            analyticTransformComplex!(wave, daughters, x̂, fromPlan, cWav.averagingType)
        end
    else
        otherwiseTransform!(wave, daughters, x̂, fromPlan, cWav.averagingType)
    end
    wave = permutedims(wave, [1, ndims(wave), (2:(ndims(wave) - 1))...])
    ax = axes(wave)
    wave = wave[1:n1, ax[2:end]...]
    if N == 1
        wave = dropdims(wave, dims = 3)
    end

    return wave
end

end # module
