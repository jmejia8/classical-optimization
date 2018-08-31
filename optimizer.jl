function exhaustiveSearch(f::Function, a::Real, b::Real, n::Int)
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
            warn("No minimum found in (a, b) or the minimum is at limits.")
            return a, b
        end
        
    end
end

function acota(f, x0, Δ)
    if f(x0 - abs(Δ)) >= f(x0 - abs(Δ)) >= f(x0 - abs(Δ))
        println("Δ is positive")
    elseif f(x0 - abs(Δ)) <= f(x0 - abs(Δ)) <= f(x0 - abs(Δ))
        Δ *= -1
        println("Δ is negative")
    else
        warn("Nothing to do...")
    end

    k = 0
    x = x0
    x_old = x0
    x_new = x + 2^k * Δ

    f_x = f(x)
    f_new = f(x_new)

    while f_new < f_x
        k += 1
        x_old = x
        x = x_new
        f_x = f_new

        x_new = x + 2^k * Δ
        f_new = f(x_new)
        
        @printf("k = %d \t x(k-1) = %.2f \t x(k) = %.2f \t x(k+1) = %.2f \n", k, x_old, x, x_new)

    end

    println("Interval: ($x_old, $x_new ) ")
    return  x_old, x_new

end

function intervalosMitad(f, a, b, ε)
    L = b - a
    xm = (a + b) / 2

    fm = f(xm)

    while !(abs(L) < ε)
        x1 = a + L / 4
        x2 = b - L / 4

        fx1 = f(x1)
        fx2 = f(x2)

        if fx1 < fm
            b, xm, fm = xm, x1, fx1
        elseif fx2 < fm
            a, xm, fm = xm, x2, fx2
        else
            a, b = x1, x2
        end

        @printf("a = %.4f, b = %.4f \t f(xm) = %.4f\n", a, b, fm)

        L = b - a

    end

    return a, b

end
