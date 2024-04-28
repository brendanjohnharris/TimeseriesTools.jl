module TermExt
using TimeseriesTools
using Term
using Term.Progress
import TimeseriesTools.progressmap
import Base: HasEltype, HasLength, HasShape, SizeUnknown, _similar_for, _similar_shape,
             IteratorEltype, IteratorSize, Generator, _array_for

const PROGRESSMAP_PROGRESS = ProgressBar(; columns = :detailed, width = 92,
                                         transient = true)

function progressmap(f, args...; kwargs...)
    icargs = zip(args...) |> enumerate |> collect
    T = typejoin(Base.return_types(f, eltype.(args))...)
    out = Vector{T}(undef, length(icargs))
    pbar = PROGRESSMAP_PROGRESS
    foreachprogress(icargs, pbar; kwargs...) do (i, cargs)
        out[i] = f(cargs...)
    end
    stop!(pbar)
    return out
end
function progressmap(f, A::AbstractArray; kwargs...) # Is this faster? Not great types
    T = typejoin(Base.return_types(f, (eltype(A),))...)
    out = Array{T}(undef, size(A))
    pbar = PROGRESSMAP_PROGRESS
    foreachprogress(collect(enumerate(A)), pbar; kwargs...) do (i, a)
        out[i] = f(a)
    end
    stop!(pbar)
    return Base.collect_similar(A, out)
end
function progressmap(f, As::AbstractDimArray...; kwargs...) # Check output type
    DimensionalData.comparedims(As...)
    data = progressmap(f, map(parent, As)...; kwargs...)
    rebuild(first(As); data)
end

end # module
