"""
Funciones de apoyo para implementar el algoritmo principal.
"""

"""
	rankMethod(G, A, c, b)

Implementación del método de rango para resolver problemas de programación cuadrática 
(PPC) con restricciones ``Ax = b``.

# Arguments
- `G::Matriz{Float64}(n, n)`: Matriz de la definición del problema cuadrático.
- `A::Matrix{Float64}(m, n)`: La matriz de restricciones.
- `c::Vector{Float64}(n)`: Vector de costos de la funciób objetivo.
- `b::Vector{Float64}(n)`: Vector de constantes de las restricciones
"""
function rankMethod(G, A, c, b)
	
	# Calculando λ con (2.4)
	AGinv = A * inv(G)
	
	# Resolviendo el sistema lineal con \ para no calcular inversas
	# λ = -inv(AGinv * A') * (-b - AGinv * c)
	λ = (AGinv * A') \ (-b - AGinv * c)

	# Calculando x_* con (2.3)
	Ginv
	x_star = (-Ginv * c) - (Ginv * A' * λ)
	
	return x_star
end
