module Operators
using TimeseriesTools
export 𝐵, 𝐹, 𝛥, ℒ!, 𝒯

# ? Some basic time-series operators

# Backshift operator
𝐵!(x) = circshift!(x, 1)
𝐵²!(x) = circshift!(x, 2)
𝐵³!(x) = circshift!(x, 3)
𝐵⁴!(x) = circshift!(x, 4)
𝐵⁵(x) = circshift!(x, 5)
𝐵!(x, n) = circshift!(x, n)
𝐵(x) = (x = deepcopy(x);
        𝐵!(x);
        x[1] = NaN)

# Lag operator
ℒ!(x) = circshift!(x, -1)
ℒ²!(x) = circshift!(x, -2)
ℒ³!(x) = circshift!(x, -3)
ℒ⁴!(x) = circshift!(x, -4)
ℒ⁵!(x) = circshift!(x, -5)
ℒ!(x, n) = circshift!(x, -n)

# Shift operator (operates on time indices)
𝒯(t) = x -> set(x, Ti(times(x) .+ t))

# Difference operator
𝛥!(x) = (x .= x .- 𝐵(x))

end
