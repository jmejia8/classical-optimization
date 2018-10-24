include("optimizer.jl")

function test()
    f(x) = x^2 + 54 / x

    exhaustiveSearch(f, 0, 5, 2000)
    boundingPhase(f, 0.6, 0.5)
    intervalHalving(f, 0, 5, 1e-3)
    goldenSection(f, 0, 5, 1e-3)
    fibonacciSearch(f, 0, 5, 20)
    powell(f, 1, 0.1, 1e-3, 1e-3; debug=true)
    newtonRhapson(f, 0.5, 1e-2, 1e-3; debug=true)
    bisection(f, 2, 5, 1e-6; debug=true)
    secante(f, 2, 5, 1e-6; debug=true)
    cubic(f, 1, 0.1, 1e-3, 1e-3; debug=true)
end

function test2()
    f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
    
    c = 0.0
    g(x) = 10length(x) + sum((x .- c).^2 .- 10cos.(2π*(x .- c)))
 
    D = 10
    
    x0 = 10ones(D)
    Δ  = 2ones(D)

    return evop(g, x0, Δ, 1e-3; debug=true)
    
    #####################################################
    x1 = [4, 4]
    λ = 5
    ε = 1e-3
    N = 50
    randomSearch(f, x1, λ, ε, N; debug  =true)
    #####################################################
end

function test3()
    c = 0.0
    f(x) = 10length(x) + sum((x .- c).^2 .- 10cos.(2π*(x .- c)))
    
    D = 10

    x1 = 5rand(D)
    λ = 10
    ε = 1e-9
    N = 10D
    
    randomSearch(f, x1, λ, ε, N; debug  =true)
end

test2()
