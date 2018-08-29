include("optimizer.jl")

function test()
    f(x) = x^2 + 54 / x

    robustSearch(f, 0, 5, 10)
end

test()
