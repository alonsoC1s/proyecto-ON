"""
Funciones de apoyo para implementar el algoritmo principal.
"""
module Utils

export rankMethod, linprog, 𝒜, solve2_8, solve2_11, solve2_9, klee_minty

using LinearAlgebra, SparseArrays
using JuMP
using GLPK

function rankMethod(G, A, c, b)

    if !isposdef(G)
        throw(ArgumentError("`G` debe ser positiva definida."))
    end

    # Teniendo que convertir a A en full por el operador \
    if issparse(A)
        A = copy(Matrix(A))
    end

    # Guardando productos que se reutilizan
    Ginv = inv(G)
    AGinv = A * Ginv

    # Calculando λ con (2.4)
    # Resolviendo el sistema lineal con \ para no calcular inversas
    λ = (AGinv * A') \ (-b - AGinv * c)

    # Calculando x_* con (2.3)
    x_star = (-Ginv * c) - (Ginv * A' * λ)

    return x_star
end


function linprog(A, b, n_eq)
    # Guardando número de variables
    n = size(A, 2)

    # Inicializando modelo
    model = Model(GLPK.Optimizer)

    # Agregando n variables al modelo
    @variable(model, x[1:n])

    # Agregando restricciones de igualdad y desigualdad
    if n_eq > 0
        @constraint(model, A[1:n_eq, :] * x .== b[1:n_eq])
    end
    
    @constraint(model, A[n_eq + 1:end, :] * x .<= b[n_eq + 1:end])

    optimize!(model)

    # Regresando el valor óptimo
    return value.(x)
end


function solve2_8(G, A_k, g_k)
    return rankMethod(G, A_k, g_k, zeros(size(A_k, 1)))
end

function solve2_9(A, b, x_k, d_k, atol=1e-12)
    # Filtrar las j's tales que Aj^t dk > 0
    noW_k = findall(A * d_k .> atol)
    return findmin((b[noW_k] - (A[noW_k, :] * x_k)) ./ (A[noW_k, :] * d_k))
end


function solve2_11(g_k, A, W_k, n_eq)
	# Feo, pero necesario en caso de tener Sparse
	if issparse(A)
		A_k = copy(Matrix(A))
	else
		A_k = copy(A)
    end

    A_k[.!W_k, :] .= 0
    lag = A_k' \ (-g_k)
    return (lag[1:n_eq], lag[n_eq + 1:end])
end


function klee_minty(n::Int)
	# Vector de constantes
    b = (2 * ones(n)).^(1:n) .- 1

	# Regresando G, c, A, b en ese orden
	return 1e-4 * I(n), ones(n), [LowerTriangular(ones(n, n) + Diagonal(ones(n))); -I(n)], [b; zeros(n)]
end

end
