using Test

include("../funciones.jl")

@testset "Utilities básicas" begin
	# Testeando método del rango
    G = [2 1 0 0; 1 2 0 0; 0 0 4 2; 0 0 2 2]
    A = [1 1 0 0; 0 1 1 0]
	b = [2, 1]
	c = [0, 0, 1, 1]

	@test rankMethod(G, A, c, b) == [1, 1, 0, -1/2]

	# Cambiando G a una no-positiva definida
	G = [1 1 0 0; -1 -1 0 0; 0 0 4 2; 0 0 2 2]
	@test_throws ArgumentError rankMethod(G, A, c, b)

end 
