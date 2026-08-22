# FAQ

## Why is my result reducing over time, it shouldn't be!

The default scheme for odd-order derivatives is upwind. It is fast and low-memory, and it
introduces numerical dispersion: sharp peaks and discontinuities smooth out over time.
The [WENO Scheme](@ref adschemes) does not have this property; it is more expensive and
can fail with exotic boundary conditions.

To see numerical dispersion in action, take a look at this example:

```julia
using ModelingToolkit, MethodOfLines, LinearAlgebra, OrdinaryDiffEq, DomainSets
using ModelingToolkit: Differential
using Plots
plotlyjs()

@parameters t x y
@variables u(..)
Dx = Differential(x)
Dy = Differential(y)
Dt = Differential(t)

t_min = 0.0
t_max = 7.0
x_min = -5.0
x_max = 5.0
y_min = -5.0
y_max = 5.0
dx = 0.25
dy = 0.25
order = 3

A = 1.0
x0 = 0.0
y0 = 0.0

sigma_x = 1.0
sigma_y = 1.0
theta = 0.0

# create function for bivariate gaussian

function bivariate_gaussian(
        x, y; A = 1.0, x0 = 0.0, y0 = 0.0, sigma_x = 1.0, sigma_y = 1.0, theta = 0.0)
    a = cos(theta)^2 / (2 * sigma_x^2) + sin(theta)^2 / (2 * sigma_y^2)
    b = -1 * sin(2 * theta) / (4 * sigma_x^2) + sin(2 * theta) / (4 * sigma_y^2)
    c = sin(theta)^2 / (2 * sigma_x^2) + cos(theta)^2 / (2 * sigma_y^2)
    return A * exp(-(a * (x - x0)^2 + 2 * b * (x - x0) * (y - y0) + c * (y - y0)^2))
end;

# equation system

eq = Dt(u(t, x, y)) - Dx(u(t, x, y)) - Dy(u(t, x, y)) ~ 0

# boundary conditions

bcs = [
    u(t_min, x, y) ~ bivariate_gaussian(
        x, y; A = A, x0 = x0, y0 = y0, sigma_x = sigma_x, sigma_y = sigma_y, theta = theta),
    u(t, x_min, y) ~ u(t, x_max, y),    # <--- SHOULD BE A PERIODIC BOUNDARY CONDITION
    u(t, x, y_min) ~ u(t, x, y_max)    # <--- SHOULD BE A PERIODIC BOUNDARY CONDITION
]

domains = [
    t ∈ Interval(t_min, t_max),
    x ∈ Interval(x_min, x_max),
    y ∈ Interval(y_min, y_max)
]

@named pdesys = PDESystem([eq], bcs, domains, [t, x, y], [u(t, x, y)])

discretization = MOLFiniteDifference(
    [x => dx, y => dy], t; advection_scheme = UpwindScheme())

sys, tspan = symbolic_discretize(pdesys, discretization)
prob = ODEProblem(mtkcompile(sys), nothing, tspan)

sol = solve(prob, SSPRK54(), dt = 0.01, saveat = 0.1)

discrete_t = sol[t]
solu = sol[u(t, x, y)]

anim = @animate for i in 1:length(discrete_t)
    surface(solu[i, :, :], camera = (55.0, 30.0), size = (500, 500), zlabel = ("z"),
        zlims = (0.0, 1.2), xlabel = ("x"), ylabel = ("y"), clims = (0.0, 1.0),
        title = "t = $i")
end

gif(anim, "mol_convection_2d_test.gif", fps = 5)
```

![convection test](https://github.com/SciML/MethodOfLines.jl/assets/9698054/45f4ace0-6291-478d-abb8-93d68ae3c9aa)

## Why is my large discretized system taking so long to compile?

For a first-order time-dependent system, `discretize` returns a `DAEProblem`
whose generated residuals are operations over array slices. For patterns that
have an array form, the number of symbolic equations does not grow with the
grid. Call `solve(prob)` and let OrdinaryDiffEq select the DAE algorithm.

That scaling does not survive `mtkcompile`, which scalarizes the array
equations. Every other path compiles through it: stationary systems
(`NonlinearProblem`), `StaggeredGrid` (`SplitODEProblem`), the compiled
`ODEProblem` fallback and `analytic` paths, and explicit Runge–Kutta methods
(`Tsit5`, `SSPRK54`, …), constructed as:

```julia
sys, tspan = symbolic_discretize(pdesys, discretization)
prob = ODEProblem(mtkcompile(sys), nothing, tspan)
sol = solve(prob, Tsit5())
```

Some patterns emit one scalar equation per affected point even on the
`DAEProblem` path (integrals, higher mixed derivatives, two-domain interfaces,
schemes without a coefficient split). Those equations grow with resolution.

## Why are the corners of my domain held at 0?

The corner points do not have a valid discretization, and as such are held at 0.
