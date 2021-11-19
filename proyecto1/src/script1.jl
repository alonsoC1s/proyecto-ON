# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# # Problema 1
# $$
# \min q(x) = (x-1)^2 + (y - 2.5)^2
# $$
# Sujeto a
# $$
# \begin{align*}
#     -x + 2y - 2 &\leq 0 \\
#     x + 2y -6 &\leq 0 \\
#     x - 2y -2 &\leq 0 \\
#     -x &\leq 0 \\
#     -y &\leq 0
# \end{align*}
# $$
#
# Con $x_0 = (2, 0)^\top$ & $W_0 = \{3 \}$

include("Solvers.jl")
using .Solvers, LinearAlgebra

# +
G = 2 * I(2)
c = [-2, -5]
b = [2, 6, 2, 0, 0]

A = [-1 2;
    1 2;
    1 -2;
    -1 0;
    0 -1]

W_0 = [falses(2); true; falses(2)]

n_eq = 0
# -

or = activeSetMethod(G, c, A, b, n_eq, W_0)

or.q_star

# ![Problema 1](../docs/p1.png)

# Lo hizo bien, salvo por el valor de $q$. `quadprog` de Matlab da los mismos resultados. Llegó en ... iteraciones

or.iters

#
