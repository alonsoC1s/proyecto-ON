"""
Funciones de apoyo para implementar el algoritmo principal.
"""

"""
	rankMethod(G, A, c, b)

Implementación del método de rango para resolver problemas de programación cuadrática 
(PPC) con solamente restricciones de igualdad.
"""
function rankMethod(G, A, c, b)
	
	# Calculando λ con (2.4)
	AGinv = A * inv(G)
	
	# Resolviendo el sistema lineal con \ para no calcular inversas
	# λ = -inv(AGinv * A') * (-b - AGinv * c)
	λ = (AGinv * A') \ (-b - AGinv * c)

	# Calculando x_* con (2.3)
	x_star = -inv(G) * c - inv(G) * A' * λ
	
	return x_star
end
