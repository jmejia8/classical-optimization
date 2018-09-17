import Printf.@printf

function exhaustiveSearch(f::Function, a::Real, b::Real, n::Int; debug::Bool=false)
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
            debug && warn("No minimum found in (a, b) or the minimum is at limits.")
            return a, b
        end
        
    end
end

function boundingPhase(f::Function, x0::Real, Δ::Float64; debug::Bool=false)
    if f(x0 - abs(Δ)) >= f(x0 - abs(Δ)) >= f(x0 - abs(Δ))
        debug && println("Δ is positive")
    elseif f(x0 - abs(Δ)) <= f(x0 - abs(Δ)) <= f(x0 - abs(Δ))
        Δ *= -1
        debug && println("Δ is negative")
    else
        warn("Nothing to do...")
    end

    k = 0
    x, x_old = x0, x0
    x_new = x + 2^k * Δ

    f_x, f_new = f(x), f(x_new)

    while f_new < f_x
        k += 1
        x_old, x, f_x = x, x_new, f_new

        x_new = x + 2^k * Δ
        f_new = f(x_new)
        
        debug && @printf("k = %d \t x(k-1) = %.2f \t x(k) = %.2f \t x(k+1) = %.2f \n", k, x_old, x, x_new)

    end

    return  x_old, x_new

end

function intervalHalving(f::Function, a::Real, b::Real, ε::Float64; debug::Bool=false)
    L = b - a
    xm = (a + b) / 2

    fm = f(xm)

    while !(abs(L) < ε)
        x1, x2 = a + L / 4, b - L / 4

        fx1, fx2 = f(x1), f(x2)

        if fx1 < fm
            b, xm, fm = xm, x1, fx1
        elseif fx2 < fm
            a, xm, fm = xm, x2, fx2
        else
            a, b = x1, x2
        end

        debug && @printf("a = %.4f, b = %.4f \t f(xm) = %.4f\n", a, b, fm)

        L = b - a

    end

    return a, b

end

function goldenSection(f::Function, a::Real, b::Real, ε::Float64; debug::Bool=false)
    fw(w) = f(w*(b - a)+a)

    # configuration
    aw, bw = 0.0, 1.0
    Lw = bw - aw
    τ = (√(5) + 1) \ 2

    w1 = aw + τ
    w2 = bw - τ
    
    fw1, fw2 = fw(w1), fw(w2)
    k = 1
    
    while abs(Lw) > ε
        Lw = τ^k
        if fw1 < fw2
            aw = w2
            
            w2 = w1
            w1 = aw + τ * Lw        
            
            fw2 = fw1
            fw1 = fw(w1)
    
        else
            bw = w1
           
            w1 = w2
            w2 = bw - τ * Lw
    
            fw1 = fw2
            fw2 = fw(w2)
            
        end
        
        k += 1
        debug && @printf("k = %d \t aw = %.4f, bw = %.4f \t f(w1) = %.4f  f(w2) = %.4f\n",k, aw, bw, fw1, fw2)

    end


    return aw*(b-a) + a, bw*(b-a) + a
end
