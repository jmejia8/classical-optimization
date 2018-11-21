import Calculus: gradient, hessian
import LinearAlgebra: dot, norm, I

include("optimizer.jl")

function marquardt(f::Function, x0::Vector{Float64}, ε::Float64; λ::Float64 = 1e2, M::Int = 1000, debug::Bool=false)
    k = 0
    ∇f(x) = gradient(f, x)


    x = x0
    ∇fx = ∇f(x)
    fx  = f(x)

    k = 1
    while norm(∇fx) > ε || k >= M

        H = hessian(f, x)
        s = -inv( H + λ*I ) * ∇fx
        
        x_old = x
        x += s


        fx_old = fx
        fx = f(x)
        
        debug && print("k = $k\t x = $x")
        debug && @printf("\t λ = %.2f \t f(x) = %0.3g\n", λ, fx)

        if fx < fx_old
            λ *= 0.5
        else            
            λ *= 2.0
            x = x_old
            continue
        end

        ∇fx = ∇f(x)

        k += 1

    end

    return x, fx

end

function test()
    f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2

    #############################################
    ε = 1e-3
    x0 = [0.0, 0.0]
    marquardt(f, x0, ε; debug = true)
    #############################################


end

test()
