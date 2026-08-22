# [How it works](@id hiw)

MethodOfLines.jl makes heavy use of `Symbolics.jl` and `SymbolicUtils.jl`, namely its rule matching features to recognize terms which require particular discretizations.

See [here](https://github.com/SciML/MethodOfLines.jl/blob/master/src/MOL_discretization.jl) for the highest level overview of the algorithm.

Given your discretization and `PDESystem`, we take each independent variable defined on the space to be discretized and create a corresponding range. We then take each dependent variable and create an array of symbolic variables to represent it in its discretized form. This is stored in a [`DiscreteSpace`](https://github.com/SciML/MethodOfLines.jl/blob/master/src/discretization/discretize_vars.jl) object, a useful abstraction.

We recognize boundary conditions, i.e. whether they are on the upper or lower ends of the domain, or periodic [here](https://github.com/SciML/PDEBase.jl/blob/master/src/parse_boundaries.jl), and use this information to construct the interior of the domain for each equation [here](https://github.com/SciML/MethodOfLines.jl/blob/master/src/system_parsing/interior_map.jl). Each PDE is matched to each dependent variable in this step by which variable is of highest order in each PDE, with precedence given to time derivatives. This dictates which boundary conditions reduce the size of the interior for which PDE, and ensures the discretized system is balanced, with as many equations as discrete unknowns. A supported interior contributes one array equation over a slice, not one equation per grid point.

Boundary conditions that occupy a contiguous multi-point face are emitted as one
array equation over that face. Single-point faces (every 1D boundary) stay
scalar. The implementation is [`array_bc_eqs`](https://github.com/SciML/MethodOfLines.jl/blob/master/src/discretization/discretize_equations.jl);
the pointwise fallback is [`generate_bc_eqs.jl`](https://github.com/SciML/MethodOfLines.jl/blob/master/src/discretization/generate_bc_eqs.jl).

Each PDE interior is emitted as symbolic array equations over slices of the
discrete variables, e.g. `D(u[2:n-1]) ~ (u[1:n-2] .- 2u[2:n-1] .+ u[3:n]) ./ dx^2`.
The equation count for a supported interior does not depend on resolution.
Points whose stencil is not the translation-invariant interior stencil (the
frame) stay pointwise. Term matching still lives in
[`generate_finite_difference_rules.jl`](https://github.com/SciML/MethodOfLines.jl/blob/master/src/discretization/generate_finite_difference_rules.jl);
the array reconstruction is [`discretize_equations.jl`](https://github.com/SciML/MethodOfLines.jl/blob/master/src/discretization/discretize_equations.jl).
Unsupported patterns fall back to one scalar equation per point.

The result is a system of differential-algebraic or nonlinear equations with
as many equations as unknowns.
See [the generated Brusselator system](@ref brusssys)
for the compiled (`mtkcompile`) scalar form at low point count.
Time-dependent systems are normally converted directly to a `DAEProblem`,
which preserves symbolic array equations. Stationary systems are compiled with
`ModelingToolkit.mtkcompile` into a `NonlinearProblem`; the optional compiled
`ODEProblem` path also uses `mtkcompile` and scalarizes the equations. The problem
constructors generate executable Julia functions with `RuntimeGeneratedFunctions`.
