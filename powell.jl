include("optimizer.jl")

function powell(f::Function, x0::Array{Float64}, unidirectionalSearch::Function; Basis = nothing, ε = 1e-3, debug::Bool = false)
    D = length(x0)
    
    if Basis == nothing
        Basis = zeros(D, D)
        for i = 1:D
            Basis[i, i] = 1.0
        end
    end

    x = x0
    while true
        x1 = unidirectionalSearch( f, x, Basis[1,:] ) # f(x + λs)

        x = x1
        for i = 2:D
            s = Basis[i,:]
            x = unidirectionalSearch( f, x, s ) # f(x + λs)
        end


        x = unidirectionalSearch( f, x, Basis[1,:] )

        d = x - x1

        debug && print("x = ", x)
        debug && @printf("\tf(x) = %.3g \t det(d) =  %.2f\n", f(x), norm(d))

        if norm(d) < ε || det(Basis) == 0
            break
        end
 
        for i = D:-1:2
            Basis[i,:] = Basis[i-1,:]
        end
        Basis[1,:] = d / norm(d)

    end

    return x

end

function test()
    f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2

    unidirectionalSearch(f, x, s) = begin
        Δ = 0.1
        λ0 = 0.0
        a,b = boundingPhase( λ -> f(x + λ*s), λ0, Δ; debug = true)
        
        ε_gs = 1e-10
        a2,b2 = goldenSection( λ -> f(x + λ*s), a, b, ε_gs; debug = true)

        λ = b2
        x + λ*s
    end

    x0 = [0.0, 4]
    powell(f, x0, unidirectionalSearch; debug = true)


end

test()
