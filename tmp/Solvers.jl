"""
Implementacion principal del metodo de conjunto activo
"""
module Solvers

# Incluyendo y usando el módulo de utils definido en el archivo funcines.jl
include("Utils.jl")
using .Utils, LinearAlgebra

export OptimizeResult, activeSetMethod

struct OptimizeResult{T<:Real}
	iters::Integer
	x_star::Vector{T}
	q_star::T
end

function activeSetMethod(G, c, A, b, n_eq,  W_k=nothing, maxiter = 100, atol = 1e-9)
    k = 0

	# Encontrar x_0 con simplex
    x_k = linprog(A, b, n_eq)

	# Definir W0
    if W_k === nothing
        W_k = [trues(n_eq); falses(size(A, 1) - n_eq)]
    end

	if n_eq == 0
		W_k = [falses(Int(size(A, 1)/2)); trues(Int(size(A, 1)/2))]
	end

    g_k = G * x_k + c


    while k < maxiter
        # Obtener d_k de (2.8) con W_k
        d_k = solve2_8(G, A[W_k, :], g_k)

        if norm(d_k, Inf) >= atol
        # if !isapprox(d_k, zero(d_k), atol = atol)
            # ==============RAMA 1====================
            #  Encontrar α gorro y j con (2.9)
            α, j = solve2_9(A, b, x_k, d_k, atol)
            x_k = x_k + min(1, α) .* d_k

            print("Rama 1. ||d_k|| = $(norm(d_k, Inf)), ")
			print("q(x) = $(x_k' * G * x_k + c' * x_k), α=$(α) ") 

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
                return (x_k, λ, μ)
            end

            println("")

            # Quitando j del conjunto activo W_k
            W_k[n_eq+j] = false
        end
    end
end

end
