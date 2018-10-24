mutable struct Solution
    x::Array{Float64}
    f::Real
end

function nelder_mead(f, simplex::Array{Solution}; α::Real = 1.0,  β::Real = 0.5, γ::Real = 2.0, ε::Real = 1e-5)
    D = length(simplex[1].x)
    N = length(simplex)

    while true
        sort!(simplex, lt = (p, q) -> p.f < q.f)

        best  = deepcopy(simplex[1])
        worst = deepcopy(simplex[end])
        # x_g = 

        # calculate centroid
        c = zeros(D)

        for i = 1:N-1
            c += simplex[i].x
        end

        c /= N

        # perform reflection
        x_r = (1.0 + α)*c - α*worst.x
        f_r = f(x_r)
        if f_r < best.f
            # expansion
            x_e = (1 + γ) * x_r - γ*c
            f_e = f(x_e)
            simplex[end].x = f_e < best.f ? x_e : x_r
            simplex[end].f = f_e < best.f ? f_e : f_r
            continue
        end

        is_worst = true
        for i = 1:N-1
            if f_r <= simplex[i].f
                is_worst = false
                break
            end
        end

        if is_worst
            if f_r < worst.f
                simplex[end].x = x_r
                simplex[end].f = f_r
            end

            # contraction
            x_contraction = β*worst.x + (1 - β)*c
            f_contraction = f(x_contraction)
            if f_contraction > worst.f
                # shrinkage
                for i = 1:N
                    simplex[i].x = 0.5*(simplex[i].x + best.x)
                    simplex[i].f = f(simplex[i].x)
                end
            else
                simplex[end].x = x_contraction
                simplex[end].f = f_contraction
            end
        else
            simplex[end].x = x_r
            simplex[end].f = f_r
        end

        # stop?
        if sqrt(sum([ (simplex[i].f - f(c))^2 for i=1:N ]))/N <= ε
            break
        end

        println(simplex)
    end

    return simplex
end

function test()
    f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
    N = 3
    D = 2

    simplex = Solution[]
    x1 = zeros(2)
    x2 = [1.0, 0]
    x3 = [0.0, 1]

    push!(simplex, Solution( x1, f(x1) ))
    push!(simplex, Solution( x2, f(x2) ))
    push!(simplex, Solution( x3, f(x3) ))

    nelder_mead(f, simplex)


end

test()