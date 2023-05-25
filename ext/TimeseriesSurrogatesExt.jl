# module TimeseriesSurrogatesExt
using TimeseriesTools
using TimeseriesSurrogates

import TimeseriesSurrogates: surrogenerator, surrogate

# function surrogenerator(x::UnivariateTimeSeries, method::T) where {T <: Surrogate}
#     sg = surrogenerator(x.data, method)
#     function _surrogenerator(aargs...)
#         s = deepcopy(x)
#         s.data = sg(aargs...)
#     end
#     return _surrogenerator
# end

# end # module
