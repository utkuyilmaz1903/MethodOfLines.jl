# Self-contained benchmark: ArrayDiscretization vs ScalarizedDiscretization for WENO.
#
# Measures, per strategy: discretization wall time, equation counts, symbolic expression
# size (treesize + unique nodes), mtkcompile time, ODEProblem construction, solve time,
# and solution agreement. Expression size is constant in grid resolution for the array
# form when WENO has a slice representation; otherwise both arms fall back to pointwise
# and the times agree.
#
# Run from the MethodOfLines package root:
#   julia --project=benchmark benchmark/weno/array_vs_scalar.jl

using MethodOfLines, ModelingToolkit, DomainSets, Symbolics, SymbolicUtils
using ModelingToolkit: get_eqs
using OrdinaryDiffEqSSPRK: SSPRK22
using SciMLBase

treesize(x) = let u = Symbolics.unwrap(x)
    SymbolicUtils.iscall(u) ? 1 + sum(treesize, SymbolicUtils.arguments(u); init = 0) : 1
end

function unique_nodes(x)
    seen = IdDict{Any, Nothing}()
    function walk(y)
        y = Symbolics.unwrap(y)
        haskey(seen, y) && return
        seen[y] = nothing
        SymbolicUtils.iscall(y) || return
        foreach(walk, SymbolicUtils.arguments(y))
    end
    walk(x)
    return length(seen)
end

isarreq(eq) = any((eq.lhs, eq.rhs)) do side
    u = Symbolics.unwrap(side)
    !(u isa AbstractArray) && SymbolicUtils.symtype(u) <: AbstractArray
end
isinterior(eq) = occursin("Differential(t", string(eq))

function advection_system(n; dim = 1)
    @parameters t
    if dim == 1
        @parameters x
        @variables u(..)
        Dt = Differential(t)
        Dx = Differential(x)
        eq = Dt(u(t, x)) ~ -Dx(u(t, x))
        bcs = [u(0, x) ~ exp(-100 * (x - 0.3)^2), u(t, 0) ~ 0.0, u(t, 1) ~ 0.0]
        dom = [t ∈ Interval(0.0, 0.05), x ∈ Interval(0.0, 1.0)]
        @named ps = PDESystem(eq, bcs, dom, [t, x], [u(t, x)])
        return ps, [x => 1 / (n - 1)], u
    else
        @parameters x y
        @variables u(..)
        Dt = Differential(t)
        Dx = Differential(x)
        Dy = Differential(y)
        Dxx = Differential(x)^2
        Dyy = Differential(y)^2
        eq = Dt(u(t, x, y)) ~ -Dx(u(t, x, y)) - Dy(u(t, x, y)) +
            0.01 * (Dxx(u(t, x, y)) + Dyy(u(t, x, y)))
        bcs = [
            u(0, x, y) ~ exp(-50 * ((x - 0.3)^2 + (y - 0.3)^2)),
            u(t, 0, y) ~ 0.0, u(t, 1, y) ~ 0.0,
            u(t, x, 0) ~ 0.0, u(t, x, 1) ~ 0.0,
        ]
        dom = [
            t ∈ Interval(0.0, 0.05), x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 1.0),
        ]
        @named ps = PDESystem(eq, bcs, dom, [t, x, y], [u(t, x, y)])
        return ps, [x => 1 / (n - 1), y => 1 / (n - 1)], u
    end
end

function burgers_system(n)
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dx = Differential(x)
    eq = Dt(u(t, x)) ~ -u(t, x) * Dx(u(t, x))
    bcs = [
        u(0, x) ~ 0.5 + 0.4 * sin(2π * x),
        u(t, 0) ~ u(t, 1),
    ]
    dom = [t ∈ Interval(0.0, 0.1), x ∈ Interval(0.0, 1.0)]
    @named ps = PDESystem(eq, bcs, dom, [t, x], [u(t, x)])
    return ps, [x => 1 / (n - 1)], u
end

function run_case(pdesys, dxs, t, u, strat; dt = 1.0e-3)
    disc = MOLFiniteDifference(
        dxs, t;
        discretization_strategy = strat,
        advection_scheme = WENOScheme(),
    )
    t_disc = @elapsed (sys, _) = symbolic_discretize(pdesys, disc)
    eqs = get_eqs(sys)
    int = filter(isinterior, eqs)
    return (;
        t_disc,
        neqs = length(eqs),
        nint = length(int),
        narr = count(isarreq, eqs),
        narr_int = count(eq -> isinterior(eq) && isarreq(eq), eqs),
        int_nodes = sum(eq -> treesize(eq.lhs) + treesize(eq.rhs), int; init = 0),
        int_unique = sum(eq -> unique_nodes(eq.lhs) + unique_nodes(eq.rhs), int; init = 0),
        sys,
        disc,
        u,
    )
end

function run_solve(pdesys, disc, u; dt = 1.0e-3)
    t_prob = @elapsed (prob = discretize(pdesys, disc))
    t_solve = @elapsed (
        sol = solve(
            prob, SSPRK22(); dt = dt, adaptive = false, save_everystep = false
        )
    )
    return (; t_prob, t_solve, sol, u)
end

function print_header()
    println(
        rpad("case", 18), rpad("n", 6), rpad("strat", 8),
        rpad("t_disc", 10), rpad("eqs", 6), rpad("narr", 6),
        rpad("int_nodes", 12), rpad("int_uniq", 10),
        rpad("t_prob", 10), rpad("t_solve", 10)
    )
    return
end

function print_row(label, n, nm, r, s = nothing)
    tp = s === nothing ? "-" : string(round(s.t_prob; digits = 3))
    ts = s === nothing ? "-" : string(round(s.t_solve; digits = 3))
    println(
        rpad(label, 18), rpad(n, 6), rpad(nm, 8),
        rpad(round(r.t_disc; digits = 3), 10), rpad(r.neqs, 6), rpad(r.narr_int, 6),
        rpad(r.int_nodes, 12), rpad(r.int_unique, 10),
        rpad(tp, 10), rpad(ts, 10)
    )
    return
end

println("===== WENO ArrayDiscretization vs ScalarizedDiscretization =====")
print_header()

# Warmup
begin
    ps, dxs, u = advection_system(21)
    @parameters t
    for st in (ScalarizedDiscretization(), ArrayDiscretization())
        run_case(ps, dxs, t, u, st)
    end
end

strategies = [
    ("scalar", ScalarizedDiscretization()),
    ("array", ArrayDiscretization()),
]

@parameters t

for (label, maker, sizes) in [
        ("1D advection", n -> advection_system(n; dim = 1), [100, 400]),
        ("1D Burgers", burgers_system, [100, 400]),
        ("2D adv-diff", n -> advection_system(n; dim = 2), [16, 32]),
    ]
    println("\n--- ", label, " ---")
    for n in sizes
        pdesys, dxs, u = maker(n)
        sols = Dict{String, Any}()
        for (nm, st) in strategies
            r = try
                run_case(pdesys, dxs, t, u, st)
            catch e
                println(
                    rpad(label, 18), rpad(n, 6), rpad(nm, 8),
                    "FAILED: ", first(split(sprint(showerror, e), '\n'))[1:min(end, 60)]
                )
                continue
            end
            s = try
                run_solve(pdesys, r.disc, u)
            catch e
                print_row(label, n, nm, r)
                println(
                    "  solve FAILED: ",
                    first(split(sprint(showerror, e), '\n'))[1:min(end, 80)]
                )
                continue
            end
            sols[nm] = s
            print_row(label, n, nm, r, s)
        end
        if haskey(sols, "array") && haskey(sols, "scalar")
            try
                da = sols["array"].sol.u[end]
                ds = sols["scalar"].sol.u[end]
                rel = maximum(abs.(da .- ds) ./ (abs.(ds) .+ 1.0e-14))
                println("  solution relmax = ", rel)
            catch
                println("  solution comparison skipped")
            end
        end
        flush(stdout)
    end
end

println("\nDONE")
