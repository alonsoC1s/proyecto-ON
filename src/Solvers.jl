"""
Implementacion principal del metodo de conjunto activo
"""
module Solvers

# Incluyendo y usando el módulo de utils definido en el archivo funcines.jl
include("Utils.jl")
using .Utils, LinearAlgebra

export OptimizeResult, activeSetMethod

"""
	OptimizeResult{T<:Real}

El el tipo de resultado de aplicar [`activeSetMethod`](@ref), el algoritmo de optimización 
del conjunto activo para problemas de optimización cuadrática.

Si `R::OptimizeResult` es el resultado de aplicar el algoritmo del conjunto activo se 
pueden obtener las interaciones como `R.iters`, el punto óptimo como `R.x_star` y el valor 
de la función objetivo en el punto como `R.q_star`. También se puede obtener "estilo 
Matlab".

# Examples
```julia
iters, x_star, q_star = activeSetMethod(G, c, A_E, b_E, A_I, b_I)
```
"""
struct OptimizeResult{T<:Real}
	iters::Integer
	x_star::Vector{T}
	q_star::T
end

"""
	activeSetMethod(G, c, A_E, b_E, A_I, b_I, atol = 1e-9)

Aplica el algoritmo de conjunto activo para optimizar el problema cuadrático:
```math
\\min \\frac{1}{2} x^\\top G x + c^\\top x
```
Sujeto a
```math
\\begin{align*}
	A_E x = b_E \\\\
	A_I x \\leq B_I
\\end{align*}
```
"""
function activeSetMethod(G, c, A, b, n_eq, maxiter = 100, atol = 1e-9)
    k = 0
    # Concatenando A & E en una sola matriz
    #A = [A_E; A_I]
    #b = [b_E; b_I]

    # n_eq = length(b_E) # Numero de igualdades

    # Encontrar x_0 con simplex
    x_k = linprog(A, b, n_eq)

    # Definir W0 
    W_k = [trues(n_eq); falses(size(A, 1) - n_eq)]
    g_k = G * x_k + c

    while k < maxiter
        # Obtener d_k de (2.8) con W_k
        d_k = solve2_8(G, A[W_k, :], g_k)

        # norm(d_k, Inf) < atol
        if !isapprox(d_k, zero(d_k), atol = atol)
            # ==============RAMA 1====================
            #  Encontrar α gorro y j con (2.9)
            α, j = solve2_9(A, b, x_k, d_k, atol)
            x_k = x_k + min(1, α) .* d_k

            print("Rama 1. ||d_k|| = $(norm(d_k, Inf)), ")
			print("q(x) = $(x_k' * G * x_k + c' * x_k), α=$(α)") 

            if α < 1
                print("k = $(k)")
                # W_{k+1} = W_k ∪ {j}
                W_k[j] = true
            end

            println("")

            g_k = G * x_k + c

            k += 1
        else
            # ==============RAMA 2====================
            # Calcular los multiplicadores de lagrange
            λ, μ = solve2_11(g_k, A, W_k, n_eq)
            # Encontrar un j en las inequalities con el menor μ
            μ_min, j = findmin(μ)

            print("Rama 2. ")

            if μ_min <= atol
                println("j = $(j), μ=$(μ_min) ")
                # La solución es óptima
                return (x_k, λ, μ)
            end

            println("")

            # Quitando j del conjunto activo W_k
            W_k[n_eq+j] = false
        end
    end
end

end
