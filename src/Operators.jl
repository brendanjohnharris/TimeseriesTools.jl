module Operators
export 𝐵, 𝐹, 𝛥

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

# Shift operator
𝑇!(x) = circshift!(x, -1)
𝑇²!(x) = circshift!(x, -2)
𝑇³!(x) = circshift!(x, -3)
𝑇⁴!(x) = circshift!(x, -4)
𝑇⁵(x) = circshift!(x, -5)
𝑇!(x, n) = circshift!(x, -n)

# Difference operator
𝛥!(x) = (x .= x .- 𝐵(x))

end
