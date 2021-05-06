"""
    CryoGridOutput{TSol,TVars}

Helper type that stores the solution to a CryoGrid `ODEProblem` along with `DimArray` views of all
state variables. `CryoGridOutput` overrides `Base.getproperty` to allow for direct dot-syntax
access of state variables. For example, if your model has a grid variable named `T` in layer `layer`,
then for a `CryoGridOutput` value `out`, `out.layer.T` returns a `DimArray` with indexed time and
depth axes. The `ODESolution` can be accessed via `out.sol`, or for convenience, the continuous solution
at time `t` can be computed via `out(t)` which is equivalent to `withaxes(out.sol(t))`.
"""
struct CryoGridOutput{TSol,TVars}
    sol::TSol
    vars::TVars
    CryoGridOutput(sol::TSol, vars::TVars) where {TSol <: ODESolution,TVars<:NamedTuple} = new{TSol,TVars}(sol,vars)
end

"""
Constructs a `CryoGridOutput` from the given `ODESolution`.
"""
function CryoGridOutput(sol::TSol, ts=sol.t) where {TSol <: ODESolution}
    setup = sol.prob.f.f # CryoGridSetup
    log = get_log(sol, ts)
    ts_datetime = Dates.epochms2datetime.(ts*1000.0)
    # Helper functions for mapping variables to appropriate DimArrays by grid/shape.
    withdims(var::Var{name,T,OnGrid{Edges}}, arr) where {name,T} = DimArray(arr, (Z(setup.grid),Ti(ts_datetime)))
    withdims(var::Var{name,T,OnGrid{Cells}}, arr) where {name,T} = DimArray(arr, (Z(cells(setup.grid)),Ti(ts_datetime)))
    withdims(var::Var, arr) = DimArray(nestedview(arr), (Ti(ts_datetime),))
    layerstates = NamedTuple()
    for (i,node) in enumerate(setup.strat.nodes)
        name = nodename(node) 
        prog_vars = setup.meta[name][:pvars]
        diag_vars = setup.meta[name][:dvars]
        vararrays = Dict()
        # build nested arrays w/ flattened views of each variable's state trajectory
        for var in prog_vars
            arr = ArrayOfSimilarArrays([sol(t)[name][varname(var)] for t in ts]) |> flatview
            vararrays[var] = withdims(var, arr)
        end
        for var in diag_vars
            var_log = log[varname(var)]
            # Here we handle the case where the same variable is declared in multiple layers;
            # SimulationLogs returns a matrix where the first axis is each @log call;
            # NOTE: This will break if the user defines the same variable name with @log, because
            # we are assuming that the index matches the layer's index in the stratigraphy.
            # A better solution would be: https://github.com/jonniedie/SimulationLogs.jl/issues/7
            var_states = length(size(var_log)) > 1 ? var_log[i,:] : var_log
            arr = ArrayOfSimilarArrays(var_log) |> flatview
            vararrays[var] = withdims(var, arr)
        end
        # construct named tuple and merge with named tuple from previous layers
        varnames = map(var -> varname(var), keys(vararrays) |> Tuple)
        layer_nt = NamedTuple{tuple(varnames...)}(tuple(values(vararrays)...))
        layerstates = merge(layerstates, NamedTuple{tuple(name)}((layer_nt,)))
    end
    CryoGridOutput(sol, layerstates)
end

function Base.show(io::IO, out::CryoGridOutput)
    vars = out.vars
    nvars = sum(map(v -> keys(v) |> length, vars))
    println(io, "CryoGridOutput with $(length(out.sol.t)) time steps and $(nvars) variables:")
    for (name,vars) in zip(keys(vars),values(vars))
        if length(vars) > 0
            println("  $name: $(tuple(keys(vars)...))")
        end
    end
end

"""
Evaluates the continuous solution at time `t`.
"""
(out::CryoGridOutput)(t) = sol(t) |> withaxes

function Base.getproperty(out::CryoGridOutput, sym::Symbol)
    if sym in (:sol,:vars)
        getfield(out,sym)
    else
        out.vars[sym]
    end
end

export CryoGridOutput
