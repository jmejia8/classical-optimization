function robustSearch(f::Function, a::Real, b::Real, n::Int)
    Δx = (b - a) / n

    x1 = a
    x2 = x1 + Δx
    x3 = x2 + Δx

    while true
        
        if f(x1) >= f(x2) <= f(x3)
            return x1, x3
        end
        
        x1 = x2
        x2 = x3
        x3 = x2 + Δx
        
        if x3 > b
            warn("No minimum found in (a, b) or the minimum is at limits.")
            return a, b
        end
        
    end
end