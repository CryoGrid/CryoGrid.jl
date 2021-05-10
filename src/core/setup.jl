"""
    CryoGridSetup{S,G,M,C,U,P}

Defines the full specification of a CryoGrid model; i.e. stratigraphy, grids, variables, and diagnostic state. `uproto`
field is an uninitialized, prototype `ComponentArray` that holds the axis information for the prognostic state vector.
"""
struct CryoGridSetup{S,G,M,C,U,P}
    strat::S    # stratigraphy
    grid::G     # grid
    meta::M     # metadata (variable info and grids per layer)
    cache::C    # variable caches (per layer)
    uproto::U   # prototype prognostic state ComponentArray for integrator
    pproto::P   # prototype p ComponentArray for integrator (tracked parameters)
    CryoGridSetup(strat::S,grid::G,meta::M,cache::C,uproto::U,pproto::P) where {S<:Stratigraphy,G<:Grid{Edges},
        M<:NamedTuple, C<:NamedTuple, U<:AbstractArray,P<:AbstractArray} =
        new{S,G,M,C,U,P}(strat,grid,meta,cache,uproto,pproto)
end

"""
Constructs a `CryoGridSetup` from the given stratigraphy and grid. `arrayproto` keyword arg should be an array instance
(of any arbitrary length, including zero, contents are ignored) that will determine the array type used for all state vectors.
"""
function CryoGridSetup(strat::Stratigraphy, grid::Grid{Edges}; arrayproto::A=zeros(), chunk_size=nothing) where {A<:AbstractArray}
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
        prog_carr, param_carr, meta = buildlayer(node,subgrid,arrayproto,chunk_size)
        pvar_arrays[nodename(node)] = prog_carr
        param_arrays[nodename(node)] = param_carr
        layer_metas[nodename(node)] = meta
    end
    # construct named tuples containing data for each layer
    nodenames = [nodename(node) for node in strat]
    nt_prog = NamedTuple{Tuple(nodenames)}(Tuple(values(pvar_arrays)))
    nt_params = NamedTuple{Tuple(nodenames)}(Tuple(values(param_arrays)))
    nt_meta = NamedTuple{Tuple(nodenames)}(Tuple(values(layer_metas)))
    npvars = (length(meta.progvars) for meta in nt_meta) |> sum
    ndvars = (length(meta.diagvars) for meta in nt_meta) |> sum
    nparams = (length(meta.paramvars) for meta in nt_meta) |> sum
    @assert (npvars + ndvars) > 0 "No variable definitions found. Did you add a method definition for CryoGrid.variables(::L,::P) where {L<:Layer,P<:Process}?"
    @assert npvars > 0 "At least one prognostic variable must be specified."
    chunk_size = isnothing(chunk_size) ? nparams : chunk_size
    nt_cache = NamedTuple{Tuple(nodenames)}(Tuple(buildcaches(strat, nt_meta, arrayproto, chunk_size)))
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
            parameters[2]. # second value type (diagvars), Tuple
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
@inline @generated function buildstate(cache::NamedTuple, meta::M, u, du, params, t) where {names,types,M<:NamedTuple{names,types}}
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
function buildlayer(node::StratNode, grid::Grid{Edges}, arrayproto::A, chunk_size=nothing) where {A<:AbstractArray}
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
    diag_vars = filter(isdiagnostic, all_vars)
    prog_vars = filter(isprognostic, all_vars)
    alg_vars = filter(isalgebraic, all_vars)
    param_vars = filter(isparameter, all_vars)
    # check for duplicated algebraic/prognostic vars
    prog_alg_duplicated = prog_vars ∩ alg_vars
    @assert isempty(prog_alg_duplicated) "Variables $(prog_alg_duplicated) cannot be both prognostic and algebraic."
    # check for re-definition of diagnostic variables as prognostic
    prog_alg = prog_vars ∪ alg_vars
    diag_prog = filter(v -> v ∈ prog_alg, diag_vars)
    if !isempty(diag_prog)
        @warn "Variables $(Tuple(map(varname,diag_prog))) declared as both prognostic/algebraic and diagnostic. In-place modifications outside of callbacks may degrade integration accuracy."
    end
    # check for re-definition of parameters as prognostic
    param_prog = filter(v -> v ∈ prog_alg, param_vars)
    @assert isempty(param_prog) "Prognostic variables $(Tuple(map(varname,param_prog))) cannot also be parameters."
    # check for conflicting definitions of differential vars
    diff_varnames = map(v -> Symbol(:d, varname(v)), prog_alg)
    @assert all((isempty(filter(v -> varname(v) == d, all_vars)) for d in diff_varnames)) "Variable names $(Tuple(diff_varnames)) are reserved for differentials."
    # prognostic takes precedence, so we remove duplicated variables from the diagnostic variable set
    diag_vars = filter(v -> v ∉ diag_prog, diag_vars)
    # filter remaining duplicates
    diag_vars = unique(diag_vars)
    prog_vars = unique(prog_vars)
    alg_vars = unique(alg_vars)
    # convert back to tuples
    diag_vars, prog_vars, alg_vars = Tuple(diag_vars), Tuple(prog_vars), Tuple(alg_vars)
    diag_grids = buildgrids(diag_vars, grid, arrayproto)
    prog_grids = buildgrids(prog_vars, grid, arrayproto)
    alg_grids = buildgrids(alg_vars, grid, arrayproto)
    param_grids = buildgrids(param_vars, grid, arrayproto)
    # merge grids
    grids = merge(diag_grids, prog_grids, alg_grids, param_grids)
    # get variable names for diagnostic and prognostic
    dvarnames = @>> diag_vars map(varname)
    paramnames = @>> param_vars map(varname)
    params = (copy(p.default_value) for p in param_vars)
    # build parameter named tuple
    nt_params = NamedTuple{Tuple(paramnames)}(Tuple(params))
    # build metadata named tuple
    layermetadata = NamedTuple{tuple(:progvars,:diagvars,:algvars,:paramvars,:grids)}(tuple(prog_vars,diag_vars,alg_vars,param_vars,grids))
    # build component arrays for prognostic/algebraic variables and parameters
    prog_carr = isempty(prog_alg) ? similar(arrayproto, 0) : ComponentArray(merge(prog_grids, alg_grids))
    param_carr = isempty(param_vars) ? similar(arrayproto, 0) : ComponentArray(nt_params)
    return prog_carr, param_carr, layermetadata
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

"""
Constructs per-layer variable caches given the Stratigraphy and layer-metadata named tuple.
"""
function buildcaches(strat, metadata, arrayproto, chunk_size)
    map(strat) do node
        name = nodename(node)
        dvars = metadata[name].diagvars
        varnames = [varname(var) for var in dvars]
        caches = map(dvars) do dvar
            dvarname = varname(dvar)
            grid = metadata[name].grids[dvarname]
            VarCache(dvarname, grid, arrayproto, chunk_size)
        end
        NamedTuple{Tuple(varnames)}(Tuple(caches))
    end
end

export CryoGridSetup, CryoGridProblem, initialcondition!, withaxes
