"""
    CryoGridOutput{TSol}

Helper type that stores the raw output from a CryoGrid run along with `DimArray` views of all
logged variables. `CryoGridOutput` overrides `Base.getproperty` to allow for direct dot-syntax
access of state variables. For example, if your model has a grid variable named `T`, `out.T` returns a `DimArray`
with indexed time and depth axes. For OrdinaryDiffEq.jl outputs, the `ODESolution` can be accessed via `out.sol`,
or for convenience, the continuous solution at time `t` can be computed via `out(t)` which is equivalent to
`withaxes(out.sol(t))`.
"""
struct CryoGridOutput{TSol}
    ts::Vector{DateTime}
    sol::TSol
    data::NamedTuple
    CryoGridOutput(ts::Vector{DateTime}, sol::TSol, data::NamedTuple) where TSol = new{TSol}(ts, sol, data)
end
"""
Evaluates the continuous solution at time `t`.
"""
(out::CryoGridOutput)(t::Real) = CryoGrid.withaxes(out.sol(t), CryoGrid.Tile(out.sol.prob.f))
(out::CryoGridOutput)(t::DateTime) = out(Dates.datetime2epochms(t)/1000.0)
# Overrides from Base
function Base.show(io::IO, out::CryoGridOutput)
    withindent(str) = "    $str"
    countvars(x) = 1
    countvars(nt::NamedTuple) = sum(map(countvars, nt))
    format(res::Tuple) = join(res, "\n")
    describe(key, val) = withindent("$key => $(typeof(val).name.wrapper) of $(eltype(val)) with dimensions $(size(val))")
    describe(key, val::NamedTuple) = withindent("$key => \n$(format(map(withindent ∘ describe, keys(val), values(val))))")
    data = out.data
    nvars = countvars(data)
    println(io, "CryoGridOutput with $(length(out.ts)) time steps ($(out.ts[1]) to $(out.ts[end])) and $(nvars != 1 ? "$nvars variables" : "1 variable"):")
    strs = map(describe, keys(data), values(data))
    for r in strs
        println(io, r)
    end
end
DimensionalData.DimStack(out::CryoGridOutput, var::Symbol, vars::Symbol...) = DimStack((;map(n -> n => getproperty(out, n), tuple(var, vars...))...))
Base.keys(out::CryoGridOutput) = Base.propertynames(out.data)
Base.propertynames(out::CryoGridOutput) = tuple(fieldnames(typeof(out))..., propertynames(out.data)...)
function Base.getproperty(out::CryoGridOutput, sym::Symbol)
    if sym in (:sol,:ts,:data)
        getfield(out, sym)
    else
        out.data[sym]
    end
end
Base.Dict(out::CryoGridOutput) = Dict(map(k -> string(k) => getproperty(out, k), keys(out))...)

dimstr(::Ti) = "time"
dimstr(::Z) = "depth"
dimstr(::Dim{name}) where name = string(name)

"""
    CryoGridOutput(sol::TSol, tspan::NTuple{2,Float64}=(-Inf,Inf)) where {TSol<:SciMLBase.AbstractODESolution}

Constructs a `CryoGridOutput` from the given `ODESolution`. Optional argument `tspan` restricts the time span of the output.
"""
CryoGridOutput(sol::SciMLBase.AbstractODESolution, tspan::NTuple{2,DateTime}) = CryoGridOutput(sol, convert_tspan(tspan))
function CryoGridOutput(sol::SciMLBase.AbstractODESolution, tspan::NTuple{2,Float64}=(-Inf,Inf))
    # Helper functions for mapping variables to appropriate DimArrays by grid/shape.
    withdims(var::Var{name,<:CryoGrid.OnGrid{Cells}}, arr, grid, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Z(round.(typeof(1.0u"m"), cells(grid), digits=5)),Ti(ts)))
    withdims(var::Var{name,<:CryoGrid.OnGrid{Edges}}, arr, grid, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Z(round.(typeof(1.0u"m"), edges(grid), digits=5)),Ti(ts)))
    withdims(var::Var{name}, arr, zs, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Dim{name}(1:size(arr,1)),Ti(ts)))
    save_interval = ClosedInterval(tspan...)
    tile = Tile(sol.prob.f) # Tile
    grid = Grid(tile.grid.*u"m")
    ts = tile.data.outputs.t # use save callback time points
    # check if last value is duplicated
    ts = ts[end] == ts[end-1] ? ts[1:end-1] : ts
    t_mask = map(∈(save_interval), ts) # indices within t interval
    u_all = sol.(ts[t_mask])
    u_mat = reduce(hcat, u_all) # build prognostic state from continuous solution
    pax = ComponentArrays.indexmap(getaxes(tile.state.uproto)[1])
    # get saved diagnostic states and timestamps only in given interval
    savedstates = tile.data.outputs.saveval[1:length(ts)][t_mask]
    ts_datetime = Dates.epochms2datetime.(round.(ts[t_mask]*1000.0))
    allvars = variables(tile)
    progvars = tuplejoin(filter(isprognostic, allvars), filter(isalgebraic, allvars))
    diagvars = filter(isdiagnostic, allvars)
    fluxvars = filter(isflux, allvars)
    outputs = OrderedDict()
    # add all on-grid prognostic variables
    for var in filter(isongrid, progvars)
        name = varname(var)
        outputs[name] = withdims(var, u_mat[pax[name],:], grid, ts_datetime)
    end
    # add all on-grid diagnostic variables
    for var in filter(isongrid, tuplejoin(diagvars, fluxvars))
        name = varname(var)
        states = collect(skipmissing([name ∈ keys(state) ? state[name] : missing for state in savedstates]))
        if length(states) == length(ts_datetime)
            arr = reduce(hcat, states)
            outputs[name] = withdims(var, arr, grid, ts_datetime)
        end
    end
    # handle per-layer variables
    for layer in layernames(tile.strat)
        # if layer name appears in saved states or prognostic state axes, then add these variables to the output.
        if haskey(savedstates[1], layer) || haskey(pax, layer)
            # map over all savedstates and create named tuples for each time step
            layerouts = map(u_all, savedstates) do u, state
                layerout = OrderedDict()
                if haskey(state, layer)
                    layerstate = state[layer]
                    for var in keys(layerstate)
                        layerout[var] = layerstate[var]
                    end
                else
                end
                # convert to named tuple
                diagnostic_output = (;layerout...)
                if haskey(u, layer)
                    u_layer = u[layer]
                    prognostic_output = (;map(name -> name => u_layer[name], keys(u_layer))...)
                    return merge(prognostic_output, diagnostic_output)
                else
                    return diagnostic_output
                end
            end
            layerouts_combined = reduce(layerouts[2:end]; init=layerouts[1]) do out1, out2
                map(vcat, out1, out2)
            end
            # for each variable in the named tuple, find the corresponding variables
            layervars = (; map(name -> name => first(filter(var -> varname(var) == name, allvars)), keys(layerouts_combined))...)
            # map each output to a variable and call withdims to wrap in a DimArray
            outputs[layer] = map((var,out) -> withdims(var, reshape(out,1,:), nothing, ts_datetime), layervars, layerouts_combined)
        end
    end
    return CryoGridOutput(ts_datetime, sol, (;outputs...))
end

function write_netcdf!(filename::String, out::CryoGridOutput, filemode="c")
    NCD.Dataset(filename, filemode) do ds
        # this assumes that the primary state variable has time and depth axes
        NCD.defDim(ds, "time", size(out.data[1], Ti))
        NCD.defDim(ds, "depth", size(out.data[1], Z))
        t = NCD.defVar(ds, "time", Float64, ("time",), attrib=Dict("units" => NCD.CFTime.DEFAULT_TIME_UNITS))
        z = NCD.defVar(ds, "depth", Float64, ("depth",))
        t[:] = collect(dims(out.data[1], Ti))
        z[:] = ustrip.(dims(out.data[1], Z))
        for var in keys(out)
            _write_ncd_var!(ds, var, getproperty(out, var))
        end
    end
end

function _write_ncd_var!(ds::NCD.Dataset, key::Symbol, data::AbstractDimArray)
    # drop dims of length 1; would be nice if there were a squeeze(..) function...
    single_dims = filter(d -> length(d) == 1, dims(data))
    data = isempty(single_dims) ? data : dropdims(data, dims=single_dims)
    datavar = NCD.defVar(ds, string(key), Float64, map(dimstr, dims(data)))
    idx = [Colon() for ax in axes(data)]
    setindex!(datavar, Array(ustrip.(data)), idx...)
end

function _write_ncd_var!(ds::NCD.Dataset, key::Symbol, nt::NamedTuple)
    for var in keys(nt)
        _write_ncd_var!(ds, Symbol("$key.$var"), nt[var])
    end
end
