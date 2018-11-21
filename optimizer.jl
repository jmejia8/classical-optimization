import Printf.@printf
using SymPy
import LinearAlgebra: norm, det


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
            debug && @warn "No minimum found in (a, b) or the minimum is at limits."
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
        @warn "Nothing to do..."
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

function fibonacciSequence(n)
    n += 2
    sequence = zeros(n)
    sequence[1] = 1
    sequence[2] = 1

    for i = 3:n
        sequence[i] = sequence[i-1] + sequence[i-2]
    end

    return sequence

end

function fibonacciSearch(f::Function, a::Real, b::Real, n::Int; debug::Bool=false)
    F = fibonacciSequence(n)

    L = b - a
    k = 2
    
    Lk = (F[n-k+2] / F[n+2]) * L

    x1, x2 = a + Lk, b - Lk
    f1, f2 = f(x1), f(x2)

    while k < n

        k += 1
        Lk = (F[n-k+2] / F[n+2]) * L

        if f1 > f2
            a = x1

            x1, x2 = x2, b - Lk
            f1, f2 = f2, f(x2)
        elseif f1 < f2
            b = x2

            x2, x1 = x1, a +  Lk

            f2, f1 = f1, f(x1)
        else
            a, b = x1, x2
            x1, x2 = a + Lk, b - Lk
             
            f1, f2 = f(x1), f(x2)
        end

        debug && @printf("k = %d\t x1 = %.4f, x2 = %.4f \t f1 = %.4f \t f2 = %.4f \n", k, x1, x2, f1, f2)

    end

    return x1, x2

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

function powell(f::Function, x::Real, Δ::Real, TOL1::Real, TOL2::Real; debug::Bool=false)

    x1 = x
    x2 = x1 + Δ

    f1, f2 = f(x1), f(x2)
    
    if f1 > f2
        x3 = x1 + 2Δ
    else
        x3 = x1 - 2Δ
    end

    f3 = f(x3)

    Fs = Real[f1, f2, f3]
    Xs = Real[x1, x2, x3]

    while true
        i = argmin(Fs)
        X_min = Xs[i]
        F_min = Fs[i]


        a1 = (f2 - f1)/(x2- x1)
        a2 = ( (f3 - f1)/(x3 - x1) - a1 ) / (x3-x2)
        
        xb = 0.5*(x1 + x2) - a1/(2*a2)
        fb = f(xb)


        push!(Xs, xb)
        push!(Fs, fb)

        if abs(F_min - fb) <= TOL1 && abs(X_min - xb) <= TOL1
            return Xs[argmin(Fs)]
        end

        j = sortperm( Fs )[1:3]
        Xs = Xs[ j ]
        Fs = Fs[ j ]

        j = sortperm( Xs )[1:3]

        x1, x2, x3 = Xs[j]
        f1, f2, f3 = Fs[j]
        
        debug && @printf("x = %.4f \t f = %.4f\n", X_min, F_min)

    end
end


function newtonRhapson(f, x, Δx, ε::Float64; debug::Bool=false)
    df(f, x)  = (f(x + Δx) - f(x - Δx)) / (2Δx)
    ddf(f, x) = (f(x + Δx) -2f(x) + f(x - Δx)) / (Δx)^2

    k = 1
    while true
        dfx = df(f, x)
        x -= dfx / ddf(f, x)

        if  isnan(dfx) || abs(dfx) < ε
            return x
        end

        debug && @printf("x = %.4f \t f = %.4f\n", x, f(x))
    end


end


function bisection(f, a, b, ε::Float64; debug::Bool=false)
    df(x)  = 2x - 54/x^2
    # ddf(f, x) = (f(x + Δx) -2f(x) + f(x - Δx)) / (Δx)^2

    k = 1
    while true
        z = (a + b) / 2.0
        fp = df(z)

        if abs(fp) < ε 
            return z
        end

        if fp <= 0
            a = z
        elseif fp > 0
            b = z
        end

        k += 1
        debug && @printf("k = %d \t x = %.4f \t f = %.4f \t f' = %.4f\n", k, z, f(z), fp)

    end


end

function secant(f, a, b, ε::Float64; debug::Bool=false)
    df(x)  = 2x - 54/x^2

    k = 1
    while true
        fpL = df(a)
        fpR = df(b)
        z = b - (fpR * (b - a)) / ( fpR - fpL )

        fp = df(z)

        if abs(fp) <= ε 
            return z
        end

        if fp < 0
            a = z
        elseif fp > 0
            b = z
        end

        k += 1
        debug && @printf("k = %d \t x = %.4f \t f = %.4f \t f' = %.4f\n", k, z, f(z), fp)

    end


end


function cubic(f, x0, Δ, ε1::Float64, ε2::Float64; debug::Bool=false)
    @vars x
    df(a) = Float64(subs(diff(f), (x, a)))

    df_x = df(x0)
    df_x > 0 && (Δ *= -1.0 )

    xx  = x1  = x2 = x0
    df1 = df2 = df_x

    k = 0
    while true
        x_old    = xx
        df_x_old = df_x

        xx   += 2^k * Δ
        df_x = df(xx)

        if df_x * df_x_old <= 0
            x1 = x_old
            x2 = xx

            df1 = df_x_old
            df2 = df_x
            break
        end
        k += 1
    end

    f1 = f(x1)
    f2 = f(x2)


    while true
        z = ( 3*(f1 - f2) / (x2 - x1) ) + df1 + df2
        
        w = sqrt(z^2 - df1*df2)
        w = x1 > x2 ? -w : w

        μ = (df2 + w - z) / (df2 - df1 + 2w)

        if μ < 0
            x̅ = x2
        elseif μ > 1
            x̅ = x1
        else
            x̅ = x2 - μ*(x2 - x1)
        end

        fx̅ = f(x̅)
        
        if !( fx̅ < f1 )

            while fx̅ > f1
                x̅ -= 0.5*(x̅ - x1)
                fx̅ = f(x̅)
            end

        end

        df_x̅ = df(x̅)

        debug && @printf("x̅ = %.4f \t f = %.4f\n", x̅, f(x̅))


        if abs(df_x̅) <= ε1 || abs((x̅ - x1)/x̅) <= ε2
            return x̅
        elseif df_x̅ * df1 < 0
            x2 = x̅
            f2 = fx̅
        else
            x1 = x̅
            f1 = fx̅
        end
    end




end

function evop(f, x0, Δ, ε::Float64; debug::Bool=false)

    n = length(x0)
    
    x̅  = x0
    fx̅ = f(x̅)

    hyper_cube = zeros(2^n, n)
    fs = zeros(2^n)

    Δ2 = 0.5Δ
    while norm(Δ) >= ε

        for i = 1:2^n
            negat = zeros(Bool, n)
            neg = digits(i-1, base=2) .> 0


            negat[1:length(neg)] = reverse(neg)

            posit = .! negat
            
            hyper_cube[i,posit] = x̅[posit] .+ Δ2[posit]
            hyper_cube[i,negat] = x̅[negat] .- Δ2[negat]

        end

        for i = 1:2^n
            fs[i] = f(hyper_cube[i,:])
        end

        i_best = argmin(fs)

        if fs[i_best] < fx̅ 
            fx̅ = fs[i_best]
            x̅ = hyper_cube[i_best,:]
        else
            Δ = Δ2
            Δ2 = 0.5Δ2
        end

        println(x̅)
        debug && @printf("x̅ = [%.4f, %.4f] \t f = %.4f\n", x̅[1], x̅[2], f(x̅))

    end

    return x̅
end

function randomSearch(f, x1, λ, ε, N; debug::Bool=false)
    D = length(x1)
    i = 1

    f1 = f(x1)

    while λ > ε
        r = -1.0 .+ 2rand(D)

        while norm(r) > 1.0
            r = -1.0 .+ 2rand(D)
        end

        u = r ./ norm(r)

        x = x1 + λ*u
        println(x)

        fx = f(x)

        if fx < f1
            i = 1
            x1 = x
            f1 = fx
        elseif i <= N
            i += 1
        else
            λ /= 2
        end

        # println(x̅)
        debug && @printf("x1 = [%.4f, %.4f] \t f = %.4f\n", x1[1], x1[2], f(x1))


    end

    return x1
end