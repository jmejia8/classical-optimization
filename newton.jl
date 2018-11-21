import Calculus: gradient, hessian
import LinearAlgebra: dot, norm

include("optimizer.jl")

function newton(f::Function, x0::Vector{Float64}, unidirectionalSearch::Function, ε1::Float64, ε2::Float64, ε3::Float64; debug::Bool=false)
    k = 0
    ∇f(x) = gradient(f, x)

    x_old   = x0
    # p = scatter(x0[1:1], x0[2:2])
    ∇fx_old = ∇f(x_old)
    fx_old  = f(x_old)
    s = -∇fx_old

    x = x_old #unidirectionalSearch(f, x_old, s, ε1)
    # scatter!(x[1:1], x[2:2])
    ∇fx = ∇fx_old #∇f(x)
    fx = fx_old #f(x)

    k = 1
    while true

        s_old = s

        # if k%(length(x)+1) == 0 
        #     s = -∇fx
        # else
        #     s = -∇fx + (dot(∇fx, ∇fx) / dot(∇fx_old, ∇fx_old)) * s
        # end
        H = hessian(f, x)
        s = ( inv(H) ) * ∇fx
        println(dot(∇fx,s))
        
        x_old = x
        x = unidirectionalSearch(f, x, s, ε1)
        # scatter!(x[1:1], x[2:2])

        fx_old = fx
        fx = f(x)
        
        ∇fx_old = ∇fx
        ∇fx = ∇f(x)


        debug && print("k = $k x = $x")
        debug && @printf("\tf(x) = %0.3g\t", fx)
        @printf("θ = %.0f\n", rad2deg(acos(dot(s ./ norm(s), s_old ./ norm(s_old) ))) )    

        

        if norm(x_old - x) / norm(x_old) <= ε2 || norm(∇fx) <= ε3
            break
        end

        k += 1

    end

    return x, fx

end

function test()
    f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2

    unidirectionalSearch(f, x, s, ε = 1e-5) = begin
        Δ  = 0.5
        λ0 = 0.0
        a,b = boundingPhase( λ -> f(x - λ*s), λ0, Δ; debug = false)
        
        a2,b2 = goldenSection( λ -> f(x - λ*s), a, b, ε; debug = false)

        λ = b2
        println("λ = ", λ)
        x - λ*s
    end

    #############################################
    #############################################
    ε1 = ε2 = ε3 = 1e-3
    x0 = [0.0, 0.0]
    newton(f, x0, unidirectionalSearch, ε1, ε2, ε2; debug = true)
    #############################################


end

test()
