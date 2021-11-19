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
iters, x_star, q_star, μ = activeSetMethod(G, c, A_E, b_E, A_I, b_I)
```
"""
struct OptimizeResult{T<:Real}
	iters::Integer
	x_star::Union{Vector{T}, Nothing}
	q_star::Union{T, Nothing}
    μ::Union{Vector{T}, Nothing}
end

# Implementando la interfaz de iteración para poder hacer destructuring
Base.iterate(or::OptimizeResult) = (or.iters, Val(:x))
Base.iterate(or::OptimizeResult, ::Val{:x}) = (or.x_star, Val(:q))
Base.iterate(or::OptimizeResult, ::Val{:q}) = (or.q_star, Val(:μ))
Base.iterate(or::OptimizeResult, ::Val{:μ}) = (or.μ, Val(:done))
Base.iterate(or::OptimizeResult, ::Val{:done}) = nothing

# Pretty printing para optimize result
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, or::OptimizeResult)
    println("Concluyó método del conjunto activo en $(or.iters) iteraciones")
    println("El punto de paro fue:")
    display(round.(or.x_star, sigdigits=4))
end



"""
    activeSetMethod(G, c, A, b, n_eq,  W_k=nothing, maxiter = 100, atol = 1e-9)

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

# Arguments
- `G::Matriz{Float64}(n, n)`: Matriz positiva definida de la definición del problema cuadrático.
- `c::Vector{Float64}(n)`: Vector de costos de la funciób objetivo.
- `A::Matrix{Float64}(m, n)`: La matriz de restricciones.
- `b::Vector{Float64}(n)`: Vector de constantes de las restricciones
- `n_eq::Int`: Número de restricciones de igualdad del problema cuadrático.
- `W_k::BitVector(m)`: Opcional. Restricciones a tomar como activas en el primer paso del método.
- `maxiter::Int = 100`: Opcional. Limita el máximo de iteraciones del método.
- `atol::Float64`: Opcional. Determina la distancia máxima a la que puede estar d_k de cero cuando se determina si entrar a rama1.
"""
function activeSetMethod(G, c, A, b, n_eq,  W_k=nothing, maxiter = 100, atol = 1e-9)
    k = 0

    # Encontrar x_0 con simplex
    x_k = linprog(A, b, n_eq)

	# Definir W0
    if W_k === nothing
        if n_eq != 0
            W_k = [trues(n_eq); falses(size(A, 1) - n_eq)]
        else
            throw(ArgumentError("Para los problemas en los que todas las restricciones son de desigualdad se debe indicar explícitamente el conjunto `W_0`."))
        end
    end

    g_k = G * x_k + c

    while k < maxiter
        # Obtener d_k de (2.8) con W_k
        d_k = solve2_8(G, A[W_k, :], g_k)

        if norm(d_k, Inf) >= atol
            # ==============RAMA 1====================
            #  Encontrar α gorro y j con (2.9)

            α, j = solve2_9(A, b, x_k, d_k, W_k, atol)

            println("α = $(round.(α, sigdigits=4)), j = $(j) \n")

            x_k = x_k + min(1, α) .* d_k

            print("Rama 1. ||d_k|| = $(round.(norm(d_k, Inf), sigdigits=4)), ")
			print("q(x) = $(round.(0.5 * x_k' * G * x_k + c' * x_k, sigdigits=4)), ")
            print("α=$(round.(α, sigdigits=4))")

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

			if -eps(Float64) <= μ_min
                println("j = $(j), μ=$(μ_min) ")
                # La solución es óptima
                # return (x_k, λ, μ)
                return OptimizeResult(k, x_k, x_k' * G * x_k + c' * x_k, μ)
            end

            println("")

            # Quitando j del conjunto activo W_k
            W_k[n_eq+j] = false
        end
    end
    # Se sobrepasó maxiter
    return OptimizeResult(maxiter, x_k, nothing, nothing)
end

end # end module Solvers