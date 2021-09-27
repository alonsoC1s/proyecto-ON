"""
Funciones de apoyo para implementar el algoritmo principal.
"""
module  Utils
	
export rankMethod, linprog, 𝒜, PC_2_8

using LinearAlgebra
using JuMP
using GLPK

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

"""
	PC_2_8(G, A, g)

Envoltorio de metodo de resolucion para problemas cuadraticos con igualdades.

# Arguments
- `G::Matriz{Float64}(n, n)`: Matriz positiva definida de la definición del problema 
cuadrático.
- `A_k::Matrix{Float64}(m, n)`: La matriz de restricciones de W_k.
- `g_k::Vector{Float64}(n)`: g_k del algoritmo de conjunto activo
"""
function PC_2_8(G, A_k, g_k)
	return rankMethod(G, A_k, g_k, zeros(size(A_k, 1)))
end

"""
	linprog(A, b_E, b_I)

Envuelve JuMP y la rutina de simplex para devolver un punto factible al un problema 
cuadrático con restricciones de desigualdad con `m` restricciones y `n` variables.

Resuelve:
	min sum(ones(n))
	s.a A * x = b_E
		A * x <= b_I

# Arguments
- `A::Matrix(m, n)`: La matriz de restricciones del problema. 
- `b_E::Vector(n_e)`: Vector de restricciones de igualdad.
- `b_I::Vector(n_i)`: Vector de restricciones de desigualdad.
n = n_e + n_i
"""
function linprog(A_E, b_E, A_I, b_I)
	# Guardando número de igualdades & desigualdades
	# n_e = length(b_E)
	# n_i = length(b_I)
	# n = length(b_E) + length(b_I)
	n = size(A_E,2)

	# Inicializando modelo
	model = Model(GLPK.Optimizer)

	# Agregando n variables al modelo
	@variable(model, x[1:n])

	# Agregando restricciones de igualdad y desigualdad
	@constraint(model, A_E * x .== b_E)
	@constraint(model, A_I * x .<= b_I)

	optimize!(model)

	# Regresando el valor óptimo
	return value.(x)
end


"""
	𝒜(A, b, x)

Regresa los indices de las restricciones de activas (de igualdad & desigualdad)

#Arguments
- `A::Matrix(m, n)`: La matriz de restricciones del problema. 
- `b::Vector(n)`: Vector de restricciones.
- `x::Vector(m)`: Punto donde se evalua la restriccion.
- 
"""
function 𝒜(A, b, x)
	return A * x .== b
end

end
