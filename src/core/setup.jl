"""
    CryoGridSetup{S,G,U,P}

Defines the full specification of a CryoGrid model; i.e. stratigraphy, grids, variables, and diagnostic state. `uproto`
field is an uninitialized, prototype `ComponentArray` that holds the axis information for the prognostic state vector.
"""
struct CryoGridSetup{S,G,U,P}
    strat::S
    grid::G
    uproto::U
    pproto::P
    CryoGridSetup(strat::S,grid::G,uproto::U,pproto::P) where {S<:Stratigraphy,G<:Grid{Edges},U<:CryoGridState,
        P<:AbstractArray} = new{S,G,U,P}(strat,grid,uproto,pproto)
end

"""
Constructs a `CryoGridSetup` from the given stratigraphy and grid. `arrayproto` keyword arg should be an array instance
(of any arbitrary length, including zero) that will determine the array type used for all state vectors.
"""
function CryoGridSetup(strat::Stratigraphy, grid::Grid{Edges}; arrayproto::A=zeros()) where {A<:AbstractArray}
    pvar_arrays = OrderedDict()
    param_arrays = OrderedDict()
    layer_states = OrderedDict()
    for (i,node) in enumerate(strat)
        # determine subgrid for layer
        lo = strat.boundaries[i] |> dustrip
        hi = (i < length(strat) ? strat.boundaries[i+1] : grid[end]) |> dustrip
        # build subgrid using closed interval [lo,hi]
        subgrid = grid[lo..hi]
        # build layer
        prog_carr, param_carr, state = buildlayer(node,subgrid,arrayproto)
        pvar_arrays[nodename(node)] = prog_carr
        param_arrays[nodename(node)] = param_carr
        layer_states[nodename(node)] = state
    end
    # construct named tuples containing data for each layer
    nt_parr = NamedTuple{Tuple(keys(pvar_arrays))}(Tuple(values(pvar_arrays)))
    nt_params = NamedTuple{Tuple(keys(param_arrays))}(Tuple(values(param_arrays)))
    nt_state = NamedTuple{Tuple(keys(layer_states))}(Tuple(values(layer_states)))
    npvars = (length(layerstate.pvars) for layerstate in nt_state) |> sum
    ndvars = (length(layerstate.dvars) for layerstate in nt_state) |> sum
    @assert (npvars + ndvars) > 0 "No variable definitions found. Did you add a method definition for CryoGrid.variables(::L,::P) where {L<:Layer,P<:Process}?"
    @assert npvars > 0 "At least one prognostic variable must be specified."
    # construct prototype of u (prognostic state) array (note that this currently performs a copy)
    uproto = ComponentArray(nt_parr)
    # ditto for parameter array (need a hack here to get an empty ComponentArray...)
    pproto = sum(map(length, nt_params)) > 0 ? ComponentArray(nt_params) :
        ComponentArray(similar(arrayproto,0),(CAxis{NamedTuple{Tuple(keys(nt_params))}(Tuple(map(a->1:0,nt_params)))}(),))
    # reconstruct with given array type
    uproto = CryoGridState(ComponentArray(similar(arrayproto,length(uproto)), getaxes(uproto)), nt_state)
    CryoGridSetup(strat,grid,uproto,pproto)
end

"""
CryoGrid specialized constructor for ODEProblem that automatically generates the initial
condition and necessary callbacks.
"""
function CryoGridProblem(setup::CryoGridSetup, tspan, p=nothing;kwargs...)
	p = isnothing(p) ? setup.pproto : p
	# compute initial condition
	u0,_ = initialcondition!(setup, p)
	N = length(u0)
	# initialize diagonals for tridiag jacobian; use eltype of p for compatibility with autodiff
	ld, d, ud = similar(u0.x, eltype(p), N-1), similar(u0.x, eltype(p), N), similar(u0.x, eltype(p), N-1)
	func = ODEFunction(setup;jac_prototype=Tridiagonal(ld,d,ud))
	ODEProblem(func,u0,tspan,p,kwargs...)
end

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
@generated function (setup::CryoGridSetup{TStrat})(du,u,p,t) where TStrat
    nodetyps = nodetypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    @>> quote
    strat = setup.strat
    state = u.state
    du .= zero(eltype(du))
    u_x = withaxes(u)
    du_x = withaxes(du)
    end push!(expr.args)
    # Initialize variables for all layers
    for i in 1:N
        n = nodename(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        $nstate = buildstate(state.$n, u_x.$n, du_x.$n, p.$n, t)
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
Calls `initialcondition!` on all layers/processes and returns the fully constructed u0 and du0 states.
"""
@generated function initialcondition!(setup::CryoGridSetup{TStrat}, p) where TStrat
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
    state = u.state
    u_x = withaxes(u)
    du_x = withaxes(du)
    end push!(expr.args)
    # Iterate over layers
    for i in 1:N
        n = nodename(nodetyps[i])
        # create variable names
        nstate, nlayer, nprocess = Symbol(n,:state), Symbol(n,:layer), Symbol(n,:process)
        # generated code for layer updates
        @>> quote
        $nstate = buildstate(state.$n, u_x.$n, du_x.$n, p.$n, 0.0)
        $nlayer = strat.nodes[$i].layer
        $nprocess = strat.nodes[$i].process
        initialcondition!($nlayer,$nstate)
        initialcondition!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    @>> quote
    return u, du
    end push!(expr.args)
    return expr
end

@inline retrieve(x, u, typ) = x
@inline function retrieve(x, u, typ::Type{Q}) where {Q<:Quantity}
    @assert CRYOGRID_DEBUG "Unitful state variables only allowed in debug mode." # only allow in debug mode
    reinterpret(Q, x) # won't work with autodiff types
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
    pvartypes = ptypes .|> vartype
    dvartypes = dtypes .|> vartype
    # construct symbols for derivative variables; assumes no existing conflicts
    dpnames = @>> pnames map(n -> Symbol(:d,n))
    # generate state variable accessor expressions
    pacc = tuple((:(retrieve(u.$p,u,$typ)) for (p,typ) in zip(pnames,pvartypes))...,)
    dpacc = tuple((:(retrieve(du.$p,u,$typ)) for (p,typ) in zip(pnames,pvartypes))...,)
    dacc = tuple((:(retrieve(state.$d,u,$typ)) for (d,typ) in zip(dnames,dvartypes))...,)
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
    all_vars = tuple(layer_vars...,process_vars...)
    @debug "Building layer $(nodename(node)) with $(length(all_vars)) variables: $(all_vars)"
    # check for (permissible) duplicates between variables
    groups = groupby(var -> varname(var), all_vars)
    for (name,gvars) in filter(g -> length(g.second) > 1, groups)
        # if any duplicate variable deifnitions do not match, raise an error
        @assert reduce(==, gvars) "Found one or more conflicting definitions of $name in $gvars"
    end
    diag_vars = filter(x -> isdiagnostic(x), all_vars)
    prog_vars = filter(x -> isprognostic(x), all_vars)
    param_vars = filter(x -> isparameter(x), all_vars)
    # check for re-definition of diagnostic variables as prognostic
    diag_prog = filter(v -> v ∈ prog_vars, diag_vars)
    if !isempty(diag_prog)
        @warn "Variables $(tuple(map(varname,diag_prog)...,)) declared as both prognostic and diagnostic. In-place modifications outside of callbacks may cause integration errors."
    end
    # check for re-definition of parameters as prognostic
    param_prog = filter(v -> v ∈ prog_vars, param_vars)
    @assert isempty(param_prog) "Prognostic variables $(tuple(map(varname,param_prog)...,)) cannot also be parameters."
    # prognostic takes precedence, so we remove duplicated variables from the diagnostic variable set
    diag_vars = filter(v -> v ∉ diag_prog, diag_vars)
    # filter remaining duplicates
    diag_vars = unique(diag_vars)
    prog_vars = unique(prog_vars)
    # convert back to tuples
    diag_vars, prog_vars = tuple(diag_vars...), tuple(prog_vars...)
    diag_grids = buildgrids(diag_vars, grid, arrayproto)
    prog_grids = buildgrids(prog_vars, grid, arrayproto)
    param_grids = buildgrids(param_vars, grid, arrayproto)
    # merge grid
    grids = merge(diag_grids, prog_grids, param_grids)
    # get variable names for diagnostic and prognostic
    dvarnames = @>> diag_vars map(varname)
    pvarnames = @>> prog_vars map(varname)
    paramnames = @>> param_vars map(varname)
    dstates = (similar(arrayproto, size(diag_grids[d])) for d in dvarnames)
    params = (similar(arrayproto, size(param_grids[p])) for p in paramnames)
    nt_params = NamedTuple{Tuple(paramnames)}(Tuple(params))
    # return prognostic variable component array for top-level composition;
    # return layer state with variable name, grid information, and diagnostic state variables
    layer_state = NamedTuple{tuple(:pvars,:dvars,:grids,dvarnames...)}(tuple(prog_vars,diag_vars,grids,dstates...))
    prog_carr = isempty(prog_vars) ? similar(arrayproto, 0) : ComponentArray(prog_grids)
    param_carr = isempty(param_vars) ? similar(arrayproto, 0) : ComponentArray(nt_params)
    return prog_carr, param_carr, layer_state
end

"""
Constructs grid tuples for the given variables.
"""
function buildgrids(vars, grid::Grid{Edges}, arrayproto::A) where {A}
    if isempty(vars)
        return NamedTuple()
    end
    togrid(var::Var{name,T,OnGrid{Edges}}) where {name,T} = var.dim.f(grid)
    togrid(var::Var{name,T,OnGrid{Cells}}) where {name,T} = var.dim.f(cells(grid))
    togrid(var::Var{name,T,Shape{dims}}) where {name,T,dims} = reshape(1:prod(dims),dims...)
    names = @>> vars map(var -> varname(var))
    vars_with_names = NamedTuple{tuple(names...)}(tuple(vars...))
    grids = @>> vars_with_names map(togrid)
end

export CryoGridSetup, CryoGridProblem, initialcondition!, withaxes
