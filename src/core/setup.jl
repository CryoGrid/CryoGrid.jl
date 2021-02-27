"""
    CryoGridSetup{S,G,U,D}

Defines the full specification of a CryoGrid model; i.e. stratigraphy, grids, variables, and diagnostic state. `uproto`
field is an uninitialized, prototype `ComponentArray` that holds the axis information for the prognostic state vector.
"""
struct CryoGridSetup{S,G,U,D}
    strat::S
    grid::G
    uproto::U
    state::D
    CryoGridSetup(strat::S,grid::G,uproto::U,state::D) where {S<:Stratigraphy,G<:Grid{Edges},U,D} =
        new{S,G,U,D}(strat,grid,uproto,state)
end

"""
Constructs a `CryoGridSetup` from the given stratigraphy and grid. `arrayproto` keyword arg should be an array instance
(of any arbitrary length, including zero) that will determine the array type used for all state vectors.
"""
function CryoGridSetup(strat::Stratigraphy, grid::Grid{Edges}; arrayproto::A=zeros()) where {A<:AbstractArray}
    pvar_arrays = OrderedDict()
    layer_states = OrderedDict()
    for (i,node) in enumerate(strat)
        # determine subgrid for layer
        lo = strat.boundaries[i] |> dustrip
        hi = (i < length(strat) ? strat.boundaries[i+1] : grid[end]) |> dustrip
        # build subgrid using closed interval [lo,hi]
        subgrid = grid[lo..hi]
        # build layer
        parr, state = buildlayer(node,subgrid,arrayproto)
        pvar_arrays[nodename(node)] = parr
        layer_states[nodename(node)] = state
    end
    # construct named tuples containing data for each layer
    nt_parr = NamedTuple{Tuple(keys(pvar_arrays))}(Tuple(values(pvar_arrays)))
    nt_state = NamedTuple{Tuple(keys(layer_states))}(Tuple(values(layer_states)))
    # construct prototype of u (prognostic state) array (note that this currently performs a copy)
    uproto = ComponentArray(nt_parr)
    # reconstruct with given array type
    uproto = ComponentArray(similar(arrayproto,length(uproto)), getaxes(uproto))
    CryoGridSetup(strat,grid,uproto,nt_state)
end

"""
    withaxes(u::AbstractVector, setup::CryoGridSetup)

Constructs a `ComponentArray` from `u` for recovering variable partitions. Assumes `u` to be a prognostic state vector
of the same shape/layout as `uproto`.
"""
withaxes(u::AbstractVector, setup::CryoGridSetup) = ComponentArray(u,getaxes(setup.uproto))

"""
Generated step function (i.e. du/dt) for any arbitrary CryoGridSetup. Specialized code is generated and compiled
on the fly via the @generated macro to ensure type stability. The generated code updates each layer in the stratigraphy
in sequence, i.e for each layer 1 < i < N:

diagnostic_step!(layer i, ...)
interact!(layer i-1, ...)
prognostic_step!(layer i, ...)

Note for developers: All sections of code wrapped in quote..end blocks are generated. Code outside of quote blocks
is only executed during compilation and will not appear in the compiled version.
"""
@generated function (setup::CryoGridSetup{TStrat})(_du,_u,p,t) where {TStrat}
    nodetyps = nodetypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Build ComponentArrays directly from pre-allocated axes (should be fast and allocation free)
    # TODO also need to handle autodiff types for u and du
    @>> quote
    u = ComponentArray(_u,getaxes(setup.uproto))
    du = ComponentArray(_du,getaxes(setup.uproto))
    du .= zero(eltype(du))
    end push!(expr.args) # add to expression block
    # Declare variables
    @>> quote
    strat = setup.strat
    state = setup.state
    end push!(expr.args)
    # Initialize variables for all layers
    for i in 1:N
        n = nodename(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        $nstate = buildstate(state.$n, u.$n, du.$n, nothing, t)
        $nlayer = strat.nodes[$i].layer
        $nprocess = strat.nodes[$i].process
        end push!(expr.args)
    end
    # Diagnostic step
    for i in 1:N
        n = nodename(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        diagnosticstep!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    # Interact
    for i in 1:N-1
        n1,n2 = nodename(nodetyps[i]), nodename(nodetyps[i+1])
        n1state, n2state = Symbol(n1,:state), Symbol(n2,:state)
        n1layer, n2layer = Symbol(n1,:layer), Symbol(n2,:layer)
        n1process, n2process = Symbol(n1,:process), Symbol(n2,:process)
        @>> quote
        interact!($n1layer,$n1process,$n2layer,$n2process,$n1state,$n2state)
        end push!(expr.args)
    end
    # Prognostic step
    for i in 1:N
        n = nodename(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        prognosticstep!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    # make sure compiled method returns no value
    @>> :(return nothing) push!(expr.args)
    # emit generated expression block
    return expr
end

"""
Calls `initialcondition!` on all layers/processes and returns the fully constructed u0 and du0 state vectors.
"""
@generated function initialcondition!(setup::CryoGridSetup{TStrat}) where {TStrat}
    nodetyps = nodetypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    @>> quote
    u = similar(setup.uproto)
    du = similar(setup.uproto)
    end push!(expr.args) # add to expression block
    @>> quote
    strat = setup.strat
    state = setup.state
    end push!(expr.args)
    # Iterate over layers
    for i in 1:N
        n = nodename(nodetyps[i])
        # create variable names
        nstate, nlayer, nprocess = Symbol(n,:state), Symbol(n,:layer), Symbol(n,:process)
        # generated code for layer updates
        @>> quote
        $nstate = buildstate(state.$n, u.$n, du.$n, nothing, 0.0)
        $nlayer = strat.nodes[$i].layer
        $nprocess = strat.nodes[$i].process
        initialcondition!($nlayer,$nstate)
        initialcondition!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    @>> quote
    # return raw arrays to avoid performance penalty from ComponentArray
    return getdata(u), getdata(du)
    end push!(expr.args)
    return expr
end

"""
Generates a function from layer state type D which builds a type-stable NamedTuple of
state variables at runtime.
"""
@inline @generated function buildstate(state::D, u, du, params, t) where {names,types,D<:NamedTuple{names,types}}
    # extract variables types from D; we assume the first two parameters in the NamedTuple
    # are the prognostic and diagnostic variable names respeictively. This must be respected by buildlayer.
    # note that types.parameters[1] is Tuple{Var,...} so we call parameters again to get (Var,...)
    ptypes = types.parameters[1].parameters |> Tuple
    # again for diagnostic, assumed to be at position 2
    dtypes = types.parameters[2].parameters |> Tuple
    pnames = ptypes .|> varname
    dnames = dtypes .|> varname
    # construct symbols for derivative variables; assumes no existing conflicts
    dpnames = @>> pnames map(n -> Symbol(:d,n))
    # generate state variable accessor expressions
    pacc = tuple((:(u.$p) for p in pnames)...,)
    dpacc = tuple((:(du.$p) for p in pnames)...,)
    dacc = tuple((:(state.$d) for d in dnames)...,)
    # build state named tuple; QuoteNode is used to force names to be interpolated as symbols
    # rather than literals.
    quote
    NamedTuple{tuple($(map(QuoteNode,pnames)...),$(map(QuoteNode,dpnames)...),$(map(QuoteNode,dnames)...),:grids,
        :params, :t)}(tuple($(pacc...),$(dpacc...),$(dacc...),state.grids,params,t))
    end
end

"""
Constructs prognostic state vector and state named-tuple for the given node/layer.
"""
function buildlayer(node::StratNode, grid::Grid{Edges}, arrayproto::A) where {A<:AbstractArray}
    layer, process = node.layer, node.process
    layer_vars = variables(layer)
    @assert all([isdiagnostic(var) for var in layer_vars]) "Layer variables must be diagnostic."
    process_vars = variables(layer, process)
    all_vars = [layer_vars...,process_vars...]
    @info "Building layer $(nodename(node)) with $(length(all_vars)) variables"
    # check for (permissible) duplicates between variables
    groups = groupby(var -> varname(var), all_vars)
    for (name,gvars) in filter(g -> length(g.second) > 1, groups)
        # if any duplicate variable deifnitions do not match, raise an error
        if !reduce(==, gvars)
            error("Found one or more conflicting definitions of $name in $gvars")
        end
    end
    diag_vars = filter(x -> isdiagnostic(x), all_vars)
    prog_vars = filter(x -> isprognostic(x), all_vars)
    # filter duplicates by converting to Sets
    diag_vars = Set(diag_vars)
    prog_vars = Set(prog_vars)
    # check for re-definition of diagnostic variables as prognostic
    diag_prog = diag_vars ∩ prog_vars
    if !isempty(diag_prog)
        @warn "Variables $(tuple(diag_prog...,)) declared as both prognostic and diagnostic. In-place modifications outside of callbacks may cause integration errors."
    end
    # prognostic takes precedence, so we remove duplicated variables from the diagnostic variable set
    setdiff!(diag_vars, diag_prog)
    # convert back to tuples
    diag_vars, prog_vars = tuple(diag_vars...), tuple(prog_vars...)
    diag_carr, diag_grids = buildcomponent(diag_vars, grid, arrayproto)
    prog_carr, prog_grids = buildcomponent(prog_vars, grid, arrayproto)
    # merge grid
    grids = merge(diag_grids, prog_grids)
    # get variable names for diagnostic and prognostic
    dvarnames = @>> diag_vars map(varname)
    pvarnames = @>> prog_vars map(varname)
    # return prognostic variable component array for top-level composition;
    # return layer state with variable name, grid information, and diagnostic state variables
    pnamevals = @>> pvarnames map(var->Val{var}())
    dnamevals = @>> dvarnames map(var->Val{var}())
    dstates = (diag_carr[d] for d in dnamevals)
    layer_state = NamedTuple{tuple(:pvars,:dvars,:grids,dvarnames...)}(tuple(prog_vars,diag_vars,grids,dstates...))
    return prog_carr, layer_state
end

"""
Constructs prognostic and diagnostic component arrays and corresponding grid tuples for the given variables.
"""
function buildcomponent(vars, grid::Grid{Edges}, arrayproto::A) where {A}
    if isempty(vars)
        return similar(arrayproto,0),NamedTuple()
    end
    togrid(var::Var{name,T,OnGrid{Edges}}) where {name,T} = var.dim.f(grid)
    togrid(var::Var{name,T,OnGrid{Cells}}) where {name,T} = var.dim.f(cells(grid))
    togrid(var::Var{name,T,Shape{dims}}) where {name,T,dims} = 1:prod(dims)
    names = @>> vars map(var -> varname(var))
    vars_with_names = NamedTuple{tuple(names...)}(tuple(vars...))
    grids = @>> vars_with_names map(togrid)
    # build component array to configure axes
    carray = ComponentArray(grids)
    # rebuild with given array type
    carray = ComponentArray(similar(arrayproto,length(carray)), getaxes(carray))
    return carray, grids
end

export CryoGridSetup, initialcondition!, withaxes
