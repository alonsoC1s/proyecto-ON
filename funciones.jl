"""
Funciones de apoyo para implementar el algoritmo principal.
"""

using LinearAlgebra

"""
	rankMethod(G, A, c, b)

Implementación del método de rango para resolver problemas de programación cuadrática 
(PPC) con restricciones ``Ax = b``.

# Arguments
- `G::Matriz{Float64}(n, n)`: Matriz positiva definida de la definición del problema 
cuadrático.
- `A::Matrix{Float64}(m, n)`: La matriz de restricciones.
- `c::Vector{Float64}(n)`: Vector de costos de la funciób objetivo.
- `b::Vector{Float64}(n)`: Vector de constantes de las restricciones
"""
function rankMethod(G, A, c, b)

	if !isposdef(G)
		throw(ArgumentError("`G` debe ser positiva definida."))
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
