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
    sol::ODESolution
    vars::NamedTuple
    CryoGridOutput(sol::ODESolution, vars::NamedTuple) = new(sol, vars)
end

"""
Constructs a `CryoGridOutput` from the given `ODESolution`.
"""
function CryoGridOutput(sol::TSol, ts=sol.t) where {TSol <: SciMLBase.AbstractODESolution}
    setup = sol.prob.f.f # CryoGridSetup
    log = get_log(sol, ts)
    ts_datetime = Dates.epochms2datetime.(round.(ts*1000.0))
    # Helper functions for mapping variables to appropriate DimArrays by grid/shape.
    withdims(var::Var{name,T,OnGrid{Edges}}, arr, i) where {name,T} = DimArray(arr*oneunit(T), (Z(round.(setup.meta[i].grids[varname(var)], digits=5)u"m"),Ti(ts_datetime)))
    withdims(var::Var{name,T,OnGrid{Cells}}, arr, i) where {name,T} = DimArray(arr*oneunit(T), (Z(round.(setup.meta[i].grids[varname(var)], digits=5)u"m"),Ti(ts_datetime)))
    withdims(var::Var{name,T}, arr, i) where {name,T} = DimArray(arr*oneunit(T), (Ti(ts_datetime),))
    layerstates = NamedTuple()
    logvarnames = Set(keys(log))
    for (i,node) in enumerate(setup.strat.nodes)
        name = nodename(node) 
        prog_vars = setup.meta[name][:progvars]
        diag_vars = setup.meta[name][:diagvars]
        vararrays = Dict()
        # build nested arrays w/ flattened views of each variable's state trajectory
        for var in prog_vars
            arr = reduce(hcat, [withaxes(sol(t), setup)[name][varname(var)] for t in ts])
            vararrays[var] = withdims(var, arr, i)
        end
        for var in diag_vars
            varname_log = Symbol(name,:_,varname(var))
            pop!(logvarnames, varname_log) # remove from set
            var_log = log[varname_log]
            arr = reduce(hcat, var_log)
            vararrays[var] = withdims(var, arr, i)
        end
        # construct named tuple and merge with named tuple from previous layers
        varnames = map(var -> varname(var), keys(vararrays) |> Tuple)
        layer_nt = NamedTuple{tuple(varnames...)}(tuple(values(vararrays)...))
        layerstates = merge(layerstates, NamedTuple{tuple(name)}((layer_nt,)))
    end
    # loop over remaining (user defined) log variables
    uservars = Dict()
    for varname in logvarnames
        var_log = log[varname]
        if eltype(var_log) <: AbstractVector
            vardata = reduce(hcat, var_log)
            uservars[varname] = DimArray(vardata, (Z(1:size(vardata,1)),Ti(ts)))
        else
            uservars[varname] = DimArray(vardata, (Ti(ts),))
        end
    end
    if length(logvarnames) > 0
        nt = NamedTuple{tuple(keys(uservars)...)}(tuple(values(uservars)...))
        layerstates = merge(layerstates, (user=nt,))
    end
    CryoGridOutput(sol, layerstates)
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
    ts = Dates.epochms2datetime.(out.sol.t.*1000.0)
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
