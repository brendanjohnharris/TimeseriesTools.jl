module Operators
using TimeseriesTools
export ℬ, ℬ!, ℒ!, ℒ, 𝒯

# ? Some basic time-series operators

# Backshift operator
ℬ!(x) = circshift!(x, 1)
ℬ!(x, n) = circshift!(x, n)
ℬ(x, args...) = (x = deepcopy(x);
                 ℬ!(x, args...);
                 x)

# Lag operator
ℒ!(x) = circshift!(x, -1)
ℒ!(x, n) = circshift!(x, -n)
ℒ(x, args...) = (y = deepcopy(x); ℒ!(y, args...); y)

# Shift operator (operates on time indices)
𝒯(x, t) = set(x, 𝑡(times(x) .+ t))
𝒯(t) = Base.Fix2(𝒯, t)

end
