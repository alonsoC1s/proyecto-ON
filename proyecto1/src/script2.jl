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

# # Problema 2
# $$
# \begin{array}{ll}
# \min & \frac{1}{2} \vec{\boldsymbol{x}}^{\top} G \vec{\boldsymbol{x}}-\sum_{i=1}^{n} x_{i} \\
# \text { sujeto a } & \quad x_{1} \leq 1 \\
# & 2 \sum_{j=1}^{i-1} x_{j}+x_{i} \leq 2^{i}-1 \quad, i=2, \ldots, n \\
# x_{1}, \ldots, x_{n} \geq 0
# \end{array}
# $$
#
#
# Sea $n=15$. Aplica el método empezando con $W_{0}$ un subconjunto aleatorio de 5 entradas de los últimos 10 restricciones de positividad. Deben documentar

include("Solvers.jl")
using .Solvers, .Solvers.Utils, LinearAlgebra

n_eq = 0
G, c, A, b = klee_minty(15)

# Tomando un conjunto aleatorio de las restricciones para formar $W_0$

# +
W_0 = [falses(length(b)-10); rand((false, true), 10)]
# Ahora hay que asegurar que son exactamente 5
Δ = sum(W_0) - 5
if Δ < 0
    # Faltan restricciones
    candidates = findall(W_0 .== false)
    # Eligiendo solo los que están entre las 10 restricciones de positividad
    candidates = candidates[candidates .>= length(b) - 10]
    selected = rand(candidates, abs(Δ))

    W_0[selected] .= trues(abs(Δ))
elseif Δ > 0
    # Hay que quitar
    candidates = findall(W_0 .== true)
    selected = rand(candidates, abs(Δ))

    W_0[selected] .= falses(abs(Δ))
end

@assert sum(W_0) == 5
# -

println(W_0)

or = activeSetMethod(G, c, A, b, n_eq, W_0)

# Tomando un conjunto aleatorio de las restricciones para formar $W_0$

println(W_0)

or.q_star

# quadprog de Matlab da los mismos resultados.
# ```
# >> quadprog(G,c,A,b)
#
# Minimum found that satisfies the constraints.
#
# Optimization completed because the objective function is non-decreasing in feasible directions, to within the value of the optimality tolerance, and constraints are satisfied to within the value of the constraint tolerance.
#
# ans =
#
# 1.0e-10 *
#
# 0.2220
# 0.1175
# 0.0685
# 0.0744
# 0.1014
# ```

#
