"""
    CryoGridSetup{S,G,M,C,U,P}

Defines the full specification of a CryoGrid model; i.e. stratigraphy, grids, variables, and diagnostic state. `uproto`
field is an uninitialized, prototype `ComponentArray` that holds the axis information for the prognostic state vector.
"""
struct CryoGridSetup{S,G,M,C,U,P}
    strat::S
    grid::G
    meta::M
    cache::C
    uproto::U
    pproto::P
    CryoGridSetup(strat::S,grid::G,meta::M,cache::C,uproto::U,pproto::P) where {S<:Stratigraphy,G<:Grid{Edges},
        M<:NamedTuple, C<:NamedTuple, U<:AbstractArray,P<:AbstractArray} =
        new{S,G,M,C,U,P}(strat,grid,meta,cache,uproto,pproto)
end

"""
Constructs a `CryoGridSetup` from the given stratigraphy and grid. `arrayproto` keyword arg should be an array instance
(of any arbitrary length, including zero) that will determine the array type used for all state vectors.
"""
function CryoGridSetup(strat::Stratigraphy, grid::Grid{Edges}; arrayproto::A=zeros()) where {A<:AbstractArray}
    pvar_arrays = OrderedDict()
    param_arrays = OrderedDict()
    layer_metas = OrderedDict()
    layer_caches = OrderedDict()
    for (i,node) in enumerate(strat)
        # determine subgrid for layer
        lo = strat.boundaries[i] |> dustrip
        hi = (i < length(strat) ? strat.boundaries[i+1] : grid[end]) |> dustrip
        # build subgrid using closed interval [lo,hi]
        subgrid = grid[lo..hi]
        # build layer
        prog_carr, param_carr, cache, meta = buildlayer(node,subgrid,arrayproto)
        pvar_arrays[nodename(node)] = prog_carr
        param_arrays[nodename(node)] = param_carr
        layer_metas[nodename(node)] = meta
        layer_caches[nodename(node)] = cache
    end
    # construct named tuples containing data for each layer
    nt_prog = NamedTuple{Tuple(keys(pvar_arrays))}(Tuple(values(pvar_arrays)))
    nt_params = NamedTuple{Tuple(keys(param_arrays))}(Tuple(values(param_arrays)))
    nt_meta = NamedTuple{Tuple(keys(layer_metas))}(Tuple(values(layer_metas)))
    nt_cache = NamedTuple{Tuple(keys(layer_caches))}(Tuple(values(layer_caches)))
    npvars = (length(meta.pvars) for meta in nt_meta) |> sum
    ndvars = (length(meta.dvars) for meta in nt_meta) |> sum
    @assert (npvars + ndvars) > 0 "No variable definitions found. Did you add a method definition for CryoGrid.variables(::L,::P) where {L<:Layer,P<:Process}?"
    @assert npvars > 0 "At least one prognostic variable must be specified."
    # construct prototype of u (prognostic state) array (note that this currently performs a copy)
    uproto = ComponentArray(nt_prog)
    # ditto for parameter array (need a hack here to get an empty ComponentArray...)
    pproto = sum(map(length, nt_params)) > 0 ? ComponentArray(nt_params) :
        ComponentArray(similar(arrayproto,0),(CAxis{NamedTuple{Tuple(keys(nt_params))}(Tuple(map(a->1:0,nt_params)))}(),))
    # reconstruct with given array type
    uproto = ComponentArray(similar(arrayproto,length(uproto)), getaxes(uproto))
    CryoGridSetup(strat,grid,nt_meta,nt_cache,uproto,pproto)
end

withaxes(u::AbstractArray, setup::CryoGridSetup) = ComponentArray(u, getaxes(setup.uproto))
withaxes(u::ComponentArray, ::CryoGridSetup) = u

"""
CryoGrid specialized constructor for ODEProblem that automatically generates the initial
condition and necessary callbacks.

TODO: infer jacobian sparsity from stratigraphy definition or via SparsityDetection.jl.
"""
function CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,Float64}, p=nothing;jacproto=:tridiag,kwargs...)
	p = isnothing(p) ? setup.pproto : p
	# compute initial condition
	u0,_ = initialcondition!(setup, p, tspan)
	func = odefunction(jacproto, setup, u0, p, tspan)
	ODEProblem(func,u0,tspan,p,kwargs...)
end
# converts tspan from DateTime to float
CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,DateTime}, args...;kwargs...) = CryoGridProblem(setup, Dates.datetime2epochms.(tspan)./1000,args...;kwargs...)

odefunction(J::AbstractArray, setup::CryoGridSetup, ::AbstractArray, p, tspan) = ODEFunction(setup, jac_prototype=J)
odefunction(sym::Symbol, setup::CryoGridSetup, u0::AbstractArray, p, tspan) = odefunction(Val{sym}(), setup, u0, p, tspan)
odefunction(::Val{:dense}, setup::CryoGridSetup, u0::AbstractArray, p, tspan) = setup
function odefunction(::Val{:tridiag}, setup::CryoGridSetup, u0::AbstractArray, p, tspan)
    N = length(u0)
	ODEFunction(setup;jac_prototype=Tridiagonal(
            similar(u0, eltype(p), N-1) |> Vector,
            similar(u0, eltype(p), N) |> Vector,
            similar(u0, eltype(p), N-1) |> Vector
        )
	)
end

"""
Generated step function (i.e. du/dt) for any arbitrary CryoGridSetup. Specialized code is generated and compiled
on the fly via the @generated macro to ensure type stability. The generated code updates each layer in the stratigraphy
in sequence, i.e for each layer 1 < i < N:

diagnosticstep!(layer i, ...)
interact!(layer i-1, ...)
prognosticstep!(layer i, ...)

Note for developers: All sections of code wrapped in quote..end blocks are generated. Code outside of quote blocks
is only executed during compilation and will not appear in the compiled version.
"""
@generated function (setup::CryoGridSetup{TStrat,TGrid,TMeta})(_du,_u,p,t) where {TStrat,TGrid,TMeta}
    nodetyps = nodetypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    @>> quote
    strat = setup.strat
    cache = setup.cache
    meta = setup.meta
    _du .= zero(eltype(_du))
    du = ComponentArray(_du, getaxes(setup.uproto))
    u = ComponentArray(_u, getaxes(setup.uproto))
    end push!(expr.args)
    # Initialize variables for all layers
    for i in 1:N
        n = nodename(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        $nstate = buildstate(cache.$n, meta.$n, u.$n, du.$n, p.$n, t)
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
    # Log diagnostic variables
    for i in 1:N
        n = nodename(nodetyps[i])
        nstate = Symbol(n,:state)
        # We have to really drill down into the TMeta named tuple type to extract the variable names...
        # TMeta is a NamedTuple{names,Tuple{values...}} sooo...
        vartyps = Tuple(
            TMeta.parameters[2]. # value types, Tuple
            parameters[i]. # i'th layer metadata, NamedTuple
            parameters[2]. # value types, Tuple
            parameters[2]. # second value type (dvars), Tuple
            parameters # Var types
        )
        # iterate over each variable, extract variable name, and log it with SimulationLogs
        for var in vartyps
            nv = varname(var)
            @>> quote
            @log $nv copy($nstate.$nv)
            end push!(expr.args)
        end
    end
    # make sure compiled method returns no value
    @>> :(return nothing) push!(expr.args)
    # emit generated expression block
    return expr
end

"""
Calls `initialcondition!` on all layers/processes and returns the fully constructed u0 and du0 states.
"""
@generated function initialcondition!(setup::CryoGridSetup{TStrat}, p, tspan) where TStrat
    nodetyps = nodetypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    @>> quote
    u = similar(setup.uproto, eltype(p))
    du = similar(setup.uproto, eltype(p))
    end push!(expr.args) # add to expression block
    @>> quote
    strat = setup.strat
    cache = setup.cache
    meta = setup.meta
    end push!(expr.args)
    # Iterate over layers
    for i in 1:N
        n = nodename(nodetyps[i])
        # create variable names
        nstate, nlayer, nprocess = Symbol(n,:state), Symbol(n,:layer), Symbol(n,:process)
        # generated code for layer updates
        @>> quote
        $nstate = buildstate(cache.$n, meta.$n, u.$n, du.$n, p.$n, tspan[1])
        $nlayer = strat.nodes[$i].layer
        $nprocess = strat.nodes[$i].process
        initialcondition!($nlayer,$nstate)
        initialcondition!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    @>> quote
    return getdata(u), getdata(du)
    end push!(expr.args)
    return expr
end

"""
Generates a function from layer state type D which builds a type-stable NamedTuple of
state variables at runtime.
"""
@inline @generated function buildstate(cache::VarCache, meta::M, u, du, params, t) where {names,types,M<:NamedTuple{names,types}}
    # extract variables types from M; we assume the first two parameters in the NamedTuple
    # are the prognostic and diagnostic variable names respeictively. This must be respected by buildlayer.
    # note that types.parameters[1] is Tuple{Var,...} so we call parameters again to get (Var,...)
    ptypes = types.parameters[1].parameters |> Tuple
    # again for diagnostic, assumed to be at position 2
    dtypes = types.parameters[2].parameters |> Tuple
    pnames = ptypes .|> varname
    dnames = dtypes .|> varname
    # pvartypes = ptypes .|> vartype
    # dvartypes = dtypes .|> vartype
    # construct symbols for derivative variables; assumes no existing conflicts
    dpnames = @>> pnames map(n -> Symbol(:d,n))
    # generate state variable accessor expressions
    pacc = tuple((:(u.$p) for p in pnames)...,)
    dpacc = tuple((:(du.$p) for p in pnames)...,)
    dacc = tuple((:(retrieve(cache.$d,u)) for d in dnames)...,)
    # build state named tuple; QuoteNode is used to force names to be interpolated as symbols
    # rather than literals.
    quote
    NamedTuple{tuple($(map(QuoteNode,pnames)...),$(map(QuoteNode,dpnames)...),$(map(QuoteNode,dnames)...),
        :grids,:params,:t)}(tuple($(pacc...),$(dpacc...),$(dacc...),meta.grids,params,t))
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
    paramnames = @>> param_vars map(varname)
    # dstates = (similar(arrayproto, size(diag_grids[d])) for d in dvarnames)
    params = (similar(arrayproto, size(param_grids[p])) for p in paramnames)
    nt_params = NamedTuple{Tuple(paramnames)}(Tuple(params))
    # return prognostic variable component array for top-level composition;
    # return layer state with variable name, grid information, and diagnostic state variables
    # layer_state = NamedTuple{tuple(:pvars,:dvars,:grids,dvarnames...)}(tuple(prog_vars,diag_vars,grids,dstates...))
    layercache = VarCache(diag_grids, arrayproto)
    layermetadata = NamedTuple{tuple(:pvars,:dvars,:paramvars,:grids)}(tuple(prog_vars,diag_vars,param_vars,grids))
    prog_carr = isempty(prog_vars) ? similar(arrayproto, 0) : ComponentArray(prog_grids)
    param_carr = isempty(param_vars) ? similar(arrayproto, 0) : ComponentArray(nt_params)
    return prog_carr, param_carr, layercache, layermetadata
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
    togrid(var::Var{name,T,typeof(Scalar)}) where {name,T} = 1:1
    names = @>> vars map(var -> varname(var))
    vars_with_names = NamedTuple{tuple(names...)}(tuple(vars...))
    grids = @>> vars_with_names map(togrid)
end

export CryoGridSetup, CryoGridProblem, initialcondition!, withaxes
