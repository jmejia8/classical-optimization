include("optimizer.jl")

function powell(f::Function, x0::Array{Float64}, unidirectionalSearch::Function; Basis = nothing, ε = 1e-3)
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

        println(x)

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
        Δ = 0.5
        λ0 = 1.0
        ε_gs = 1e-4
        a,b = boundingPhase( λ -> f(x + λ*s), λ0, Δ)
        a,b = goldenSection( λ -> f(x + λ*s), a, b, ε_gs)

        λ = a
        x + λ*s
    end

    x0 = [0.0, 4]
    powell(f, x0, unidirectionalSearch)


end

test()
