"""
Implementacion principal del metodo de conjunto activo
"""

# Incluyendo y usando el módulo de utils definido en el archivo funcines.jl
include("funciones.jl")
using .Utils

using RowEchelon, LinearAlgebra


function activeSetMethod(G, c, A_E, b_E, A_I, b_I, atol=1e-12)
	k = 0
	# Definimos mousque herramientas que nos ayudaran
	A = [A_E; A_I]
	b = [b_E; b_I]

	n_eq = length(b_E) # Numero de igualdades

	# Encontrar x_0
	x_k = linprog(A_E, b_E, A_I, b_I)

	# Definir W0 
	W_k = [ones(n_eq); zeros(length(b_I))]
	g_k = G * x_k + c
	
	while (true)
		# Obtener d_k de (2.8) con W_k
		d_k = solve_2_8(G, A[findall(W_k .== 1), :], g_k)

		if !isapprox(d_k, zeros(length(d_k)), atol=atol)
			# ==============RAMA 1====================
			#  Encontrar alfa gorro y j con (2.9)
			j, α = solve2_9(A, b, x_k, d_k, atol)
			x_k = x_k + min(1, α) .* d_k
			if α < 1
				# W_{k+1} = W_k ∪ {j}
				W_k[j] = 1
			end
			g_k = G * x_k + c

			k += 1
		else
			# ==============RAMA 2====================
			# Calcular los multiplicadores de lagrange
			λ, μ = solve2_11(g_k, A, W_k, n_eq)
			# Encontrar un j en las inequalities con el menor μ
			if μ_min >= -atol
				# La solución es óptima
				return (x_k, λ, mu)
			end
			W_k[n_eq + j] = 0
		end
	end

end