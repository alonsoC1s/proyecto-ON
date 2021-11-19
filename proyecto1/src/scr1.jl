using Revise

include("Solvers.jl")
using .Solvers, LinearAlgebra

G = [2 0;
    0 2]

c = [-2, -5]
b = [2, 6, 2, 0, 0]

A = [-1 2;
    1 2;
    1 -2;
    -1 0;
    0 -1]

W_0 = [falses(2); true; falses(2)]

n_eq = 0

iters, x_star, q_star, μ = activeSetMethod(G, c, A, b, n_eq, W_0)

println(x_star)

println(q_star)