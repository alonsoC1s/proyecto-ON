---
title: |
	![](../figs/membrete-ITAM.png)
	Optimización Numérica |
	Proyecto 1 -- Reporte
author:
- Juan Carlos Sigler
- Alonso Martinez
- Paulina Carretero
bibliography: refs.bib
lang: es
geometry:
- margin=2cm
biblio-style: alphabetic
biblatexoptions: [backend=biber, citestyle=alphabetic]
header-includes:
  - \usepackage{newunicodechar}
  - \setmonofont{CMU Typewriter Text}
  - \usepackage{float}
  - \usepackage{amssymb}
  - \usepackage[ruled,vlined,linesnumbered]{algorithm2e}
  - \usepackage[utf8]{inputenc}
  - \usepackage{mathtools}
  - \definecolor{backcolour}{rgb}{0.95,0.95,0.92}
  - \definecolor{backcolour}{RGB}{229, 229, 229}
  - \definecolor{dark}{RGB}{46, 52, 64}
  - \definecolor{mainblue}{RGB}{94, 129, 172}
  - \definecolor{strgreen}{RGB}{163, 190, 140}
  - \lstdefinelanguage{julia}{
		morekeywords={function, if, end, @constraint, return, module, 
		export, using, =, copy, issparse, Matrix, isposdef, inv, 
		zeros, size, findall, findmin, ones, else, length, I, ;, trues,
		isfinite., abs., abs, eps, Float64, julia>, julia, println, Bool},
		sensitive=false,
		morestring=[s]{"}{"},
		morecomment=[l]\#,
		morecomment=[s]{"""}{"""},
		otherkeywords={*, +, -, \\, ', :, .==, .<=, ./, .>, .!, .=, 
		.^, .-},
	}
  - \lstset{language=julia,
		basicstyle={\ttfamily\small \color{dark}},
		numberstyle={\tiny \color{dark}},
		backgroundcolor=\color{backcolour},
		breaklines=true,
		numbers=left,
		keywordstyle=\color{mainblue},
		inputencoding=utf8,
		commentstyle=\color{strgreen},
	}
...

# Marco teórico

En el presente documento discutimos una implementación del algoritmo del
conjunto activo para resolver problemas de programación cuadrática con
restricciones tanto de igualdad como desigualdad, como en el problema (@eq:pc) a
continuación:

$$
	\min \frac{1}{2} \vec{x}^\top G \vec{x} + \vec{c}^\top \vec{x}
$${#eq:pc tag="P"}

## Algunos comentarios sobre la elección de lenguaje de programación


# Problemas

Probamos nuestra implementación con tres problemas y presentamos los resultados
de acuerdo a lo especificado en la asignación del proyecto.

## Problema 1: Problema chico

$$
\min q(x) = (x-1)^2 + (y - 2.5)^2
$$
Sujeto a

\begin{align*}
    -x + 2y - 2 &\leq 0 \\
    x + 2y -6 &\leq 0 \\
    x - 2y -2 &\leq 0 \\
    -x &\leq 0 \\
    -y &\leq 0
\end{align*}


Con $x_0 = (2, 0)^\top$ & $W_0 = \{3 \}$.

Para este problema, nuestro algoritmo imprime lo siguiente:

```
Rama 1. ||d_k|| = 2.5, q(x) = -1.8125, α=0.5 k = 0
Rama 1. ||d_k|| = 0.9, q(x) = -0.8000000000000007, α=1.6666666666666667 
Rama 2. j = 2, μ=0.0 
```

Obtenemos que el punto Karush-Kuhn-Tucker (KKT) es:

$$
	(\vec{x}^\star, \vec{\mu}^\star) = \left(
		\begin{bmatrix}
			1.4 \\
			1.7
		\end{bmatrix}, 
		\begin{bmatrix}
			0.399 \\
			0 \\
			0 \\
			0 \\
			0
		\end{bmatrix}
	\right)
$$
, y se alcanza en dos iteraciones.

## Problema 2: Klee-Minty

$$
\begin{array}{ll}
\min & \frac{1}{2} \vec{\boldsymbol{x}}^{\top} G \vec{\boldsymbol{x}}-\sum_{i=1}^{n} x_{i} \\
\text { sujeto a } & \quad x_{1} \leq 1 \\
& 2 \sum_{j=1}^{i-1} x_{j}+x_{i} \leq 2^{i}-1 \quad, i=2, \ldots, n \\
x_{1}, \ldots, x_{n} \geq 0
\end{array}
$$



Sea $n=15$. Aplica el método empezando con $W_{0}$ un subconjunto aleatorio de 5 entradas de los últimos 10 restricciones de positividad.

Para encontrar $W_{0}$ que cumple las restricciones impuestas usamos el
siguiente comando:

```julia
W_0 = [falses(length(b)-10); rand((false, true), 10)]
δ = sum(W_0) - 5
if δ < 0
    candidates = findall(W_0 .== false)
    candidates = candidates[candidates .>= length(b) - 10]
    selected = rand(candidates, abs(δ))
    W_0[selected] .= trues(abs(δ))
elseif δ > 0
    candidates = findall(W_0 .== true)
    selected = rand(candidates, abs(δ))
    W_0[selected] .= falses(abs(δ))
end
```

La selección es aleatoria, entonces la $W_0$ que presentamos está sujeta a
cambios, pero la presentamos con fines ilustrativos de cualquer manera.

```
julia>
Bool[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1]
```

Para este problema, nuestro algoritmo imprime lo siguiente:

```
Rama 2. j = 1, μ=0.0 
```

y muestra que el óptimo es

$$
	(\vec{x}^\star, \vec{\mu}^\star) = \left(
		\begin{bmatrix}
			0 \\
			0 \\
			\vdots \\
			0
		\end{bmatrix}, 
		\begin{bmatrix}
			0 \\
			0 \\
			\vdots \\
			1 \\
			1
		\end{bmatrix}
	\right)
$$.

La representación completa del vector $\vec{\mu}$ es

```
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
```
.

El valor óptimo de la función objetivo en $\vec{x}_\star$ es 0.0.

Para este problema la rutina `quadprog` de Matlab imprime:

```
>> quadprog(G,c,A,b)

Minimum found that satisfies the constraints.

Optimization completed because the objective function is non-decreasing in 
feasible directions, to within the value of the optimality tolerance,
and constraints are satisfied to within the value of the constraint tolerance.


ans =

   1.0e-10 *

    0.2220
    0.1175
    0.0685
    0.0744
    0.1014
```

Nuestra solución coincide con la de Matlab salvo error de redondeo.

## Problema 2: Problema Woodinfe

$$
\begin{array}{ll}
\operatorname{minimizar} & \frac{1}{2} \vec{\boldsymbol{x}}^{\top} G \vec{\boldsymbol{x}}+\vec{\boldsymbol{c}}^{\top} \vec{\boldsymbol{x}} \\
\text { sujeto a } & A \vec{x}=\vec{\boldsymbol{b}} \\
& x_{i} \geq \ell_{i} \quad \text { si } \ell_{i} \text { es finito, } \\
& x_{i} \leq u_{i} \quad \text { si } u_{i} \text { es finito. }
\end{array}
$$

- Definimos $J \subset I$ como sigue:

	- Para $x_{j} \geq \ell_{j}$ definimos $\left|g_{j}\left(x_{0}\right)\right| \leq 8 \varepsilon_{m} \max \left\{\left|\ell_{j}\right|, 1\right\} \Longrightarrow j \in J$

	- Para $x_{j} \leq u_{j}$ definimos $\left|g_{j}\left(x_{0}\right)\right| \leq 8 \varepsilon_{m} \max \left\{\left|u_{j}\right|, 1\right\} \Longrightarrow j \in J$

Donde $\varepsilon_{m}$ es el épsilon de la máquina `eps(Float64)`.

A continuación presentamos el código que se utilizó para construir el problema y
encontrar las restricciones que pertenecen a $J$. El código se presenta
recortado (excluímos el código para obtener los datos del archivo .mat y
comentarios aclaratorios) en interés de la brevedad, pero se puede encontarar código
completo en
el script incluído `script3.ipynb`.

```{julia}
n_eq = length(b)
A_eq = problem["A"]
G = I(length(c))

# Como all(isfinite.(l)) == true, se complen todas las cotas inferiores
A = [A_eq; -I(length(l))]
b = [b; l]

mask = isfinite.(u)
A = [A; I(length(u))[mask, :]]
b = [b; u[mask]]
b = b[:]

x_0 = linprog(A, b, n_eq)

# Shorthand para epsilon de maquina
εₘ = eps(Float64)
J_l = abs.(l - x_0) .<= 8 * εₘ *  max.(abs.(l), ones(length(l)))

J_u = abs.(x_0 - u) .<= 8 * εₘ *  max.(abs.(u), ones(length(u)))
J_u = J_u[isfinite.(u)]

W_0 = [trues(n_eq); J_l; J_u]
```

De el cómputo anterior obtenemos $W_0$ que documentamos a continuación,
incluídas las 138 entradas.

```{julia}
julia> println(W_0)
Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
```

El punto óptimo es:

```{julia}
julia> println(x_star)
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```

Por documentar:

- Número de iteraciones

# Conclusiones (?)

# References
