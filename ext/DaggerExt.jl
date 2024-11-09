module DaggerExt
using Dagger
using TimeseriesTools
using TimeseriesTools.ProgressLogging
import TimeseriesTools: _progressmap

isdone(x) = isready(x)
isdone(x::Dagger.Chunk) = true
isdone(x::Dagger.DArray) = all(isdone, x.chunks)
progress(D::DArray) = sum(isdone, D.chunks) / length(D.chunks)

function _progressmap(f, backend::Val{:Dagger}, args...; numblocks = 100,
                      name = "progressmap")
    nd = ndims(first(args))
    N = ceil(Int, (length(first(args)) / numblocks)^(1 / nd))
    cargs = DArray.(args, [Blocks(fill(N, nd)...)])
    _p = 0.0
    @sync begin
        out = f.(cargs...)
        @withprogress name=name begin
            while !isdone(out)
                p = progress(out)
                if p > _p
                    _p = p
                    @logprogress p
                end
                sleep(0.1)
            end
        end
        return collect(out)
    end
end

end # module
