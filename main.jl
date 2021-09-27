"""
Implementacion principal del metodo de conjunto activo
"""

# Incluyendo y usando el módulo de utils definido en el archivo funcines.jl
include("funciones.jl")
using .Utils

using RowEchelon

function apply2_9(A, b, x_k, d_k)
	# Filtrar las j's tales que Aj^t dk > 0
	noW_k = findall((A * d_k .> 0) .== 1)
	
	alfa, j = findmin(b[noW_k] - (A[noW_k, :] * x_k) ./ (A[noW_k, :] * d_k))
	#j = -1
	#alfa = -1
	#for jj in noW_k
	#	alfaalfa = (b[jj] - A[jj, :]' * x_k) / (A[jj]' * d_k)
	#	if j == -1 || alfa > alfaalfa
	#		j = jj
	#		alfa = alfaalfa
	#	end
	#end
	return (j, alfa)
end

function apply2_11(g_k, A, W_k, n_eq)
	A_k = A
	A_k[findall(W_k .!= 0), :] .= 0
	lag = A_k' \ (-g_k)
	return (lag[1:n_eq], lag[n_eq+1:end])
end

function activeSetMethod(G, c, A_E, b_E, A_I, b_I)
	# Definimos mousque herramientas que nos ayudaran
	A = [A_E; A_I]
	b = [b_E; b_I]

	n_eq = length(b_E) # Numero de igualdades

	# Encontrar x_0
	x_k = linprog(A_E, b_E, A_I, b_I)
	# Definir W0 
	W_k = 𝒜(A, b, x_k) # TODO: Puede que falle dado que no necesariamente esto es l.i
	g_k = G * x_k + c
	
	while (true)
		# Obtener d_k de (2.8) con W_k
		d_k = PC_2_8(G, A[findall(W_k .== 0), :], g_k)

		if all(d_k .!= 0)
			# ==============RAMA 1====================
			#  Encontrar alfa gorro y j con (2.9)
			(j, alfa) = apply2_9(A, b, x_k, d_k)
			x_k = x_k + min(1, alfa) * d_k
			if alfa < 1
				# W_{k+1} = W_k ∪ {j}
				W_k[j] = 1
			end
			g_k = G * x_k + c
		else
			# ==============RAMA 2====================
			# Calcular los multiplicadores de lagrange
			lambda, mu = apply2_11(g_k, A, W_k, n_eq)
			# Encontrar un j en las inequalities con el menor miu
			minMu, j = findmin(mu)
			if minMu >= 0
				# La solución es óptima
				return (x_k, lambda, mu)
			end
			W_k[n_eq + j] = 0
		end
	end

end
