using Test

includet("../src/Solvers.jl")

using .Solvers, .Solvers.Utils, LinearAlgebra

@testset "Funciones auxiliares para resolver el método de conjunto activo" begin
	# Testeando método del rango
    G = [2 1 0 0;
		 1 2 0 0;
		 0 0 4 2;
		 0 0 2 2]

    A = [1 1 0 0; 
		 0 1 1 0]

	b = [2, 1]
	c = [0, 0, 1, 1]

	@test rankMethod(G, A, c, b) == [1, 1, 0, -1 / 2]

	# Cambiando G a una no-positiva definida
	G = [1 1 0 0;
		-1 -1 0 0;
		0 0 4 2;
		0 0 2 2]
	
	@test_throws ArgumentError rankMethod(G, A, c, b)

	# Test linprog
	A = [1 1 1;
		-1 0 0;
		0 -1 0;
		0 0 -1]

	b = [3, 0, 0, 0]

	A_E = A[1, :]'
	b_E = b[1]
	A_I = A[2:end, :]
	b_I = b[2:end]

	@test linprog(A, b, 1) == [3, 0, 0]

	# Test funcion 𝒜
	## Probando sobre el hiperplano x + y + z = 3; x, y, z > 0

	# Punto sobre el hiperplano
	@test 𝒜(A, b, [1, 1, 1]) == [1, 0, 0, 0]
	@test 𝒜(A, b, [0, 0, 0]) == [0, 1, 1, 1]
end 

@testset "Tests principales del conjunto de método activo" begin
# Testeando conjunto activo con el ejemplo de clase
A = [1 1 1;
	-1 0 0;
	0 -1 0;
	0 0 -1]

b = [3, 0, 0, 0]

#= A_E = A[1, :]'
b_E = b[1]
A_I = A[2:end, :]
b_I = b[2:end] =#

G = I(3)

c = [-1, -1, -1]

@test activeSetMethod(G, c, A, b, 1) == ([1.0, 1.0, 1.0], [0.0], [0.0, 0.0, 0.0])

# Test con problema 1 del proyecto.
G = [1 0; 0 1]
c = [-1, -2.5]
b = [2, 6, 2, 0, 0]

A = [-1 2;
    1 2;
    1 -2;
    -1 0;
    0 -1]

n_eq = 0

W_k = [falses(2); true; falses(2)]

# @test activeSetMethod(G, c, A, b, n_eq, W_k)[1] == [1.4, 1.7]

end
