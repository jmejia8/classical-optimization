import Calculus.gradient
include("optimizer.jl")

function steepestDescent(f::Function, x0::Vector{Float64}, unidirectionalSearch::Function, ε1::Float64, ε2::Float64; debug::Bool=false)
    k = 0
    ∇f(x) = gradient(f, x)

    x = x0
    ∇fx = ∇f(x)
    fx = f(x)

    while norm(∇fx) >= ε1
        x = unidirectionalSearch(f, x, ∇fx, ε2)

        fx1 = f(x)

        debug && print("x = ", x)
        debug && @printf("\tf(x) = %0.3g\n", fx1)

        if abs((fx1 - fx) / fx) <= ε1
            break
        end

        ∇fx = ∇f(x)
        fx = fx1
    end

    return x, fx

end

function test()
    f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2

    unidirectionalSearch(f, x, ∇fx, ε = 1e-5) = begin
        Δ = 0.5
        λ0 = 1.0
        a,b = boundingPhase( λ -> f(x - λ*∇fx), λ0, Δ; debug = false)
        
        a2,b2 = goldenSection( λ -> f(x - λ*∇fx), a, b, ε; debug = false)

        λ = b2
        x - λ*∇fx
    end

    #############################################
    #############################################
    ε1 = 1e-5
    ε2 = 1e-5
    x0 = [0.0, 0.0]
    steepestDescent(f, x0, unidirectionalSearch, ε1, ε2; debug = true)
    #############################################


end

test()
