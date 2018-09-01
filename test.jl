include("optimizer.jl")

function test()
    f(x) = x^2 + 54 / x

    exhaustiveSearch(f, 0, 5, 2000)
    boundingPhase(f, 0.6, 0.5)
    intervalHalving(f, 0, 5, 1e-3)
end

test()
