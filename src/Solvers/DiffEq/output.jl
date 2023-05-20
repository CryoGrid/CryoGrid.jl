"""
    CryoGridOutput(sol::TSol, tspan::NTuple{2,Float64}=(-Inf,Inf)) where {TSol<:SciMLBase.AbstractODESolution}

Constructs a `CryoGridOutput` from the given `ODESolution`. Optional argument `tspan` restricts the time span of the output.
"""
InputOutput.CryoGridOutput(sol::TSol, tspan::NTuple{2,DateTime}) where {TSol<:SciMLBase.AbstractODESolution} = CryoGridOutput(sol, convert_tspan(tspan))
function InputOutput.CryoGridOutput(sol::TSol, tspan::NTuple{2,Float64}=(-Inf,Inf)) where {TSol<:SciMLBase.AbstractODESolution}
    # Helper functions for mapping variables to appropriate DimArrays by grid/shape.
    withdims(var::Var{name,<:OnGrid{Cells}}, arr, grid, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Z(round.(typeof(1.0u"m"), cells(grid), digits=5)),Ti(ts)))
    withdims(var::Var{name,<:OnGrid{Edges}}, arr, grid, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Z(round.(typeof(1.0u"m"), edges(grid), digits=5)),Ti(ts)))
    withdims(var::Var{name}, arr, zs, ts) where {name} = DimArray(arr*one(vartype(var))*varunits(var), (Dim{name}(1:size(arr,1)),Ti(ts)))
    save_interval = ClosedInterval(tspan...)
    tile = Tile(sol.prob.f) # Tile
    ts = tile.hist.vals.t # use save callback time points
    # check if last value is duplicated
    ts = ts[end] == ts[end-1] ? ts[1:end-1] : ts
    t_mask = map(∈(save_interval), ts) # indices within t interval
    u_all = sol.(ts[t_mask])
    u_mat = reduce(hcat, u_all) # build prognostic state from continuous solution
    pax = ComponentArrays.indexmap(getaxes(tile.state.uproto)[1])
    # get saved diagnostic states and timestamps only in given interval
    savedstates = tile.hist.vals.saveval[1:length(ts)][t_mask]
    ts_datetime = Dates.epochms2datetime.(round.(ts[t_mask]*1000.0))
    allvars = variables(tile)
    progvars = tuplejoin(filter(isprognostic, allvars), filter(isalgebraic, allvars))
    diagvars = filter(isdiagnostic, allvars)
    fluxvars = filter(isflux, allvars)
    outputs = OrderedDict()
    # add all on-grid prognostic variables
    for var in filter(isongrid, progvars)
        name = varname(var)
        outputs[name] = withdims(var, u_mat[pax[name],:], tile.grid, ts_datetime)
    end
    # add all on-grid diagnostic variables
    for var in filter(isongrid, tuplejoin(diagvars, fluxvars))
        name = varname(var)
        states = collect(skipmissing([name ∈ keys(state) ? state[name] : missing for state in savedstates]))
        if length(states) == length(ts_datetime)
            arr = reduce(hcat, states)
            outputs[name] = withdims(var, arr, tile.grid, ts_datetime)
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
