include("optimizer.jl")

ε = 1e-15

function example()
    FES = 0
    f(x) = begin
        FES += 1
        return x^2 - 3x - 20
    end

    x0 = 1.4; Δ  = 0.5
    
    a, b = boundingPhase(f, x0, Δ)
    a_final, b_final = intervalHalving(f, a, b, ε)

    println("  ( ", a_final, ", ", b_final, " )")
    println("  FES = ", FES)
end

function example2()
    FES = 0
    f(x) = begin
        FES += 1
        return exp(0.5(x  - 0.2)^2)
    end

    x0 = 0.0; Δ  = 0.1
    
    a, b = boundingPhase(f, x0, Δ)
    a_final, b_final = intervalHalving(f, a, b, ε; debug=false)

    println("  ( ", a_final, ", ", b_final, " )")
    println("  FES = ", FES)
end

function example4()
    FES = 0
    f(x) = begin
        FES += 1
        return 0.1(x^2 - 3x + 5)^2 + (x-3)^2
    end

    x0 = 0.0; Δ  = 0.1
    
    a, b = boundingPhase(f, x0, Δ)
    a_final, b_final = intervalHalving(f, a, b, ε; debug=false)

    println("  ( ", a_final, ", ", b_final, " )")
    println("  FES = ", FES)
end

function example5()
    FES = 0
    f(x) = begin
        FES += 1
        return 2.0x^4 - x^3 + 5x^2 - 12x + 1.0
    end

    x0 = 0.0; Δ  = 0.1
    
    a, b = boundingPhase(f, x0, Δ)
    a_final, b_final = intervalHalving(f, a, b, ε; debug=false)

    println("  ( ", a_final, ", ", b_final, " )")
    println("  FES = ", FES)
end