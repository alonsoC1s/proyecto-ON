# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# # Problema 3
# $$
# \begin{array}{ll}
# \operatorname{minimizar} & \frac{1}{2} \vec{\boldsymbol{x}}^{\top} G \vec{\boldsymbol{x}}+\vec{\boldsymbol{c}}^{\top} \vec{\boldsymbol{x}} \\
# \text { sujeto a } & A \vec{x}=\vec{\boldsymbol{b}} \\
# & x_{i} \geq \ell_{i} \quad \text { si } \ell_{i} \text { es finito, } \\
# & x_{i} \leq u_{i} \quad \text { si } u_{i} \text { es finito. }
# \end{array}
# $$
#
# - Definimos $J \subset I$ como sigue:
#
# - Para $x_{j} \geq \ell_{j}$ definimos $\left|g_{j}\left(x_{0}\right)\right| \leq 8 \varepsilon_{m} \max \left\{\left|\ell_{j}\right|, 1\right\} \Longrightarrow j \in J$
#
# - Para $x_{j} \leq u_{j}$ definimos $\left|g_{j}\left(x_{0}\right)\right| \leq 8 \varepsilon_{m} \max \left\{\left|u_{j}\right|, 1\right\} \Longrightarrow j \in J$
#
# Donde $\varepsilon_{m}$ es el épsilon de la máquina `eps(Float64)`

include("Solvers.jl")
using .Solvers, LinearAlgebra, MAT, .Solvers.Utils

problem = matread("../data/lp_afiro.mat")["Problem"]

# +
b = problem["b"] 
c = problem["aux"]["c"][:] # El [:] es para interpretar como vector
l = problem["aux"]["lo"][:] # El [:] es para interpretar como vector
u = problem["aux"]["hi"][:] # El [:] es para interpretar como vector
n_eq = length(b)
A_eq = problem["A"]
G = I(length(c))

# Como all(isfinite.(l)) == true, se complen todas las cotas inferiores
A = [A_eq; -I(length(l))]
b = [b; -l]

# Imponiendo cotas superiores
# Concatenando las filas I_j de una identidad tal que u_j < Inf
mask = isfinite.(u)
A = [A; I(length(u))[mask, :]]
# Igual cambiando b para que dimensiones coincidan
b = [b; u[mask]]
b = b[:]

# +
# Definiendo el conjunto inicial de igualdades
# Obtenemos x_0 con linprog (pero ya sabemos que va a ser vec(0))
x_0 = linprog(A, b, n_eq)

# Shorthand para epsilon de maquina
εₘ = eps(Float64)

# Obteniendo los masks que indican si se usa la restricción
# Recordamos que para las cotas inferiores g_i(x) = l_i - x_i & para las sup g_i(x) = x_i - u_i
J_l = abs.(l - x_0) .<= 8 * εₘ *  max.(abs.(l), ones(length(l)))

# Ahora para las superiores, recordando quitar los renglones donde u_j no es finito (por size de A).
J_u = abs.(x_0 - u) .<= 8 * εₘ *  max.(abs.(u), ones(length(u)))
J_u = J_u[isfinite.(u)]

J = [J_l; J_u]

# Concatenando con un índice aleatorio
R_index = rand(findall(J), 1)
notJ = falses(size(A, 1) - n_eq)
notJ[R_index] .= true
W_0 = [trues(n_eq); notJ]
# -

println(W_0)

# Nótese que el número de filas de W_0 coincide con el de A, y el tamaño de b

@assert size(A, 1) == length(W_0)
@assert size(A, 1) == length(b)

iters, x_star, q_star, μ = activeSetMethod(G, c, A, b, n_eq, copy(W_0))

print(round.(x_star, sigdigits=4))

q_star

# ## Benchmarks

# +
using BenchmarkTools
using Plots
pgfplotsx(); theme(:ggplot2)

# Corriendo el benchmark de acuerdo a las best practices
bench = @benchmark activeSetMethod($G, $c, $A, $b, $n_eq, W_k) setup=(W_k = copy(W_0));
# -

bench

ts = bench.times
h = histogram(
    ts,
    xlabel="Tiempo (s)",
    legend = :none,
    framestyle = :box);

savefig(h, "../docs/histograma2.tex")
