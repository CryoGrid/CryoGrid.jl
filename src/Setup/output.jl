"""
    CryoGridOutput

Helper type that stores the solution to a CryoGrid `ODEProblem` along with `DimArray` views of all
logged variables. `CryoGridOutput` overrides `Base.getproperty` to allow for direct dot-syntax
access of state variables. For example, if your model has a grid variable named `T` in layer `layer`,
then for a `CryoGridOutput` value `out`, `out.layer.T` returns a `DimArray` with indexed time and
depth axes. The `ODESolution` can be accessed via `out.sol`, or for convenience, the continuous solution
at time `t` can be computed via `out(t)` which is equivalent to `withaxes(out.sol(t))`. To access a
variable shared across multiple layers, use `getvar(name, out)` where `name` is the symbol corresponding
to the variable name.
"""
struct CryoGridOutput
    ts::Vector{DateTime}
    sol::ODESolution
    vars::NamedTuple
    CryoGridOutput(ts::Vector{DateTime}, sol::ODESolution, vars::NamedTuple) = new(ts, sol, vars)
end

"""
Constructs a `CryoGridOutput` from the given `ODESolution`.
"""
function CryoGridOutput(sol::TSol) where {TSol <: SciMLBase.AbstractODESolution}
    # Helper functions for mapping variables to appropriate DimArrays by grid/shape.
    withdims(::Var{name,T,<:OnGrid}, arr, zs, ts, i) where {name,T} = DimArray(arr*oneunit(T), (Z(round.(zs, digits=5)u"m"),Ti(ts)))
    withdims(::Var{name,T}, arr, zs, ts, i) where {name,T} = DimArray(arr*oneunit(T), (Ti(ts),))
    setup = sol.prob.f.f # CryoGridSetup
    ts = setup.hist.vals.t # use save callback time points
    ts_datetime = Dates.epochms2datetime.(round.(ts*1000.0))
    savedstate = setup.hist.vals.saveval
    layerstates = NamedTuple()
    for (i,node) in enumerate(setup.strat.components)
        name = componentname(node) 
        prog_vars = setup.meta[name][:progvars]
        diag_vars = setup.meta[name][:diagvars]
        vararrays = Dict()
        # prognostic variables
        for var in prog_vars
            arr = reduce(hcat, [withaxes(sol(t), setup)[name][varname(var)] for t in ts])
            vararrays[var] = withdims(var, arr, setup.meta[i].grids[varname(var)], ts_datetime, i)
        end
        # instantaneous time derivatives
        for var in prog_vars # should also include algebraic variables?
            dvar = Diagnostic(Symbol(:d,varname(var)), vartype(var), var.dim)
            var_vals = [state[name][varname(dvar)] for state in savedstate]
            arr = reduce(hcat, var_vals)
            vararrays[dvar] = withdims(dvar, arr, setup.meta[i].grids[varname(var)], ts_datetime, i)
        end
        # diagnostic variables
        for var in diag_vars
            var_vals = [state[name][varname(var)] for state in savedstate]
            arr = reduce(hcat, var_vals)
            vararrays[var] = withdims(var, arr, setup.meta[i].grids[varname(var)], ts_datetime, i)
        end
        # construct named tuple and merge with named tuple from previous layers
        varnames = map(var -> varname(var), keys(vararrays) |> Tuple)
        layer_nt = NamedTuple{tuple(varnames...)}(tuple(values(vararrays)...))
        layerstates = merge(layerstates, NamedTuple{tuple(name)}((layer_nt,)))
    end
    # loop over remaining (user defined) log variables
    # uservars = Dict()
    # for varname in logvarnames
    #     var_log = log[varname]
    #     if eltype(var_log) <: AbstractVector
    #         vardata = reduce(hcat, var_log)
    #         uservars[varname] = DimArray(vardata, (Z(1:size(vardata,1)),Ti(ts)))
    #     else
    #         uservars[varname] = DimArray(vardata, (Ti(ts),))
    #     end
    # end
    # if length(logvarnames) > 0
    #     nt = NamedTuple{tuple(keys(uservars)...)}(tuple(values(uservars)...))
    #     layerstates = merge(layerstates, (user=nt,))
    # end
    CryoGridOutput(ts_datetime, sol, layerstates)
end

"""
Evaluates the continuous solution at time `t`.
"""
(out::CryoGridOutput)(t::Real) = withaxes(out.sol(t), out.sol.prob.f.f)
(out::CryoGridOutput)(t::DateTime) = out(Dates.datetime2epochms(t)/1000.0)
"""
    getvar(var::Symbol, out::CryoGridOutput)

Similar to `getvar(::Symbol, ::CryoGridSetup, ...)` but for `CryoGridOutput`. However, this implementation
is not type stable and will newly allocate the resulting `DimArray` from concatenating along the depth dimension.
"""
function getvar(var::Symbol, out::CryoGridOutput)
    parts = []
    for layer in keys(out.vars)
        for name in keys(out.vars[layer])
            if name == var
                push!(parts, out.vars[layer][name])
            end
        end
    end
    X = reduce(vcat, parts)
    zs = reduce(vcat, [dims(A,Z).val for A in parts])
    ts = getfield(out,:ts)
    DimArray(X,(Z(zs),Ti(ts)))
end
# Overrides from Base
function Base.show(io::IO, out::CryoGridOutput)
    vars = out.vars
    nvars = sum(map(v -> keys(v) |> length, vars))
    println(io, "CryoGridOutput with $(length(out.sol.t)) time steps and $(nvars) variables:")
    for (name,vars) in zip(keys(vars),values(vars))
        if length(vars) > 0
            println(io, "  $name: $(tuple(keys(vars)...))")
        else
            println(io, "  $name: none")
        end
    end
end
function Base.getproperty(out::CryoGridOutput, sym::Symbol)
    if sym in (:sol,:log,:vars)
        getfield(out,sym)
    else
        out.vars[sym]
    end
end
