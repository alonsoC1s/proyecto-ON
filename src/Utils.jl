"""
Funciones de apoyo para implementar el algoritmo principal.
"""
module Utils

export rankMethod, linprog, 𝒜, solve2_8, solve2_11, solve2_9, klee_minty

using LinearAlgebra, SparseArrays
using JuMP
using GLPK

"""
	rankMethod(G, A, c, b)

Implementación del método de rango para resolver problemas de programación cuadrática 
(PPC) con restricciones ``Ax = b``.

# Arguments
- `G::Matriz{Float64}(n, n)`: Matriz positiva definida de la definición del problema cuadrático.
- `A::Matrix{Float64}(m, n)`: La matriz de restricciones.
- `c::Vector{Float64}(n)`: Vector de costos de la funciób objetivo.
- `b::Vector{Float64}(n)`: Vector de constantes de las restricciones
"""
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

"""
	solve2_8(G, A_k, g_k)

Envoltorio para el metodo del rango para resolucion de problemas cuadraticos con igualdades.

# Arguments
- `G::Matriz{Float64}(n, n)`: Matriz positiva definida de la definición del problema cuadrático.
- `A_k::Matrix{Float64}(m, n)`: La matriz de restricciones de W_k.
- `g_k::Vector{Float64}(n)`: g_k del algoritmo de conjunto activo
"""
function solve2_8(G, A_k, g_k)
    return rankMethod(G, A_k, g_k, zeros(size(A_k, 1)))
end

"""
	linprog(A, b, n_eq)

Envuelve JuMP y la rutina de simplex para devolver un punto factible al un problema 
cuadrático con restricciones de desigualdad con `m` restricciones y `n` variables.

Resuelve el problema
```math
\\begin{align*}
\\min 1^\\top x \\\\
A x &= b_E \\\\
A x &\\leq b_I
\\end{align*}
```

# Arguments
- `A::Matrix(m, n)`: La matriz de restricciones del problema. 
- `b::Vector(n_e)`: Vector de restricciones.
- `n_eq::Int`: Número de restricciones de igualdad del problema
n = n_e + n_i
"""
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


"""
	𝒜(A, b, x)

Regresa los indices de las restricciones de activas (de igualdad & desigualdad)

# Arguments
- `A::Matrix(m, n)`: La matriz de restricciones del problema. 
- `b::Vector(n)`: Vector de restricciones.
- `x::Vector(m)`: Punto donde se evalua la restriccion.
"""
function 𝒜(A, b, x)
    return A * x .== b
end


"""
	solve2_9(A, b, x_k, d_k, atol=1e-12)

Resuelve el problema (2.9) de las notas con tolerancia absoluta `atol`.

```math
\\check{\\alpha} = \\min_{\\substack{i \\notin W_0 \\\\ a_i^\\top b_0 > 0}}
\\left( \\frac{b_i - a_i^\\top x_0}{a_i^\\top d_0} \\right)
```

# Arguments
- `A::Matrix(m, n)`: La matriz de restricciones del problema.
- `b::Vector(n)`: Vector de restricciones.
- `x_k::Vector(n)`: Punto actual del método del conjunto activo.
- `d_k::Vector(n)`: Dirección de descenso calculada con [`solve2_8`](@ref)
- `atol::Float64`: Tolerancia absoluta (opcional).
"""
function solve2_9(A, b, x_k, d_k, atol=1e-12)
    # Filtrar las j's tales que Aj^t dk > 0
    noW_k = findall(A * d_k .> atol)
    α, j = findmin((b[noW_k] - (A[noW_k, :] * x_k)) ./ (A[noW_k, :] * d_k))
    return (α, noW_k[j])
end


"""
	solve2_11(g_k, A, W_k, n_eq)

Resuelve el sistema lineal del problema (2.11) de las notas.

```math
\\sum_{i \\in \\mathcal{E}} \\widehat{\\lambda_i} a_i +
\\sum_{i \\in \\widehat{W} \\cap \\mathcal{I}}
\\widehat{\\mu_i} a_i = - g_k
```

# Arguments
- `g_k::Vector(n)`: `G * x_k + c`
- `A::Matrix(m, n)`: La matriz de restricciones del problema.
- `W_k::BitVector(m)`: Conjunto de restricciones activas en el punto actual.
- `n_eq::Int`: Número de restricciones de igualdad del problema cuadrático.
"""
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


"""
	klee_minty(n::Int)

Regresa la matriz `G`, `A` y el vector `b` de restricciones dadas por el problema de Klee-Minty descrito en el proyecto.

# Arguments
- `n::Int`: Dimensión del problema de Klee-Minty.

# Returns
Regresa las matrices `G`, `c`, `A`, `b` en ese orden.
"""
function klee_minty(n::Int)
	# Vector de constantes
    b = (2 * ones(n)).^(1:n) .- 1

	# Regresando G, c, A, b en ese orden
	return 1e-4 * I(n), ones(n), [LowerTriangular(ones(n, n) + Diagonal(ones(n))); -I(n)], [b; zeros(n)]
end

end # end module Utils