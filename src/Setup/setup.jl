mutable struct StateHistory
    vals::Union{Missing,<:SavedValues}
    StateHistory() = new(missing)
end

"""
    CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,pax,names,obsv}

Defines the full specification of a CryoGrid model; i.e. stratigraphy, grids, variables, and diagnostic state. `uproto`
field is an uninitialized, prototype `ComponentArray` that holds the axis information for the prognostic state vector.
"""
struct CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,pax,names,obsv}
    strat::TStrat    # stratigraphy
    grid::TGrid     # grid
    meta::NamedTuple{names,TMeta}    # metadata (variable info and grids per layer)
    cache::NamedTuple{names,TCache}    # variable caches (per layer)
    hist::StateHistory
    uproto::ComponentVector{T,A,uax}   # prototype prognostic state ComponentArray for integrator
    pproto::ComponentVector{T,A,pax}   # prototype p ComponentArray for integrator (tracked parameters)
    function CryoGridSetup(
        strat::TStrat,
        grid::TGrid,
        meta::NamedTuple{names,TMeta},
        cache::NamedTuple{names,TCache},
        uproto::ComponentVector{T,A,uax},
        pproto::ComponentVector{T,A,pax};
        observed::Vector{Symbol}=Symbol[]) where
        {TStrat<:Stratigraphy,TGrid<:Grid{Edges},TMeta<:Tuple,TCache<:Tuple,T<:Number,A<:AbstractVector{T},uax,pax,names}
        new{TStrat,TGrid,TMeta,TCache,T,A,uax,pax,names,tuple(observed...)}(strat,grid,meta,cache,StateHistory(),uproto,pproto)
    end
end

"""
Constructs a `CryoGridSetup` from the given stratigraphy and grid. `arrayproto` keyword arg should be an array instance
(of any arbitrary length, including zero, contents are ignored) that will determine the array type used for all state vectors.
"""
function CryoGridSetup(strat::Stratigraphy, grid::Grid{Edges,<:Numerics.Geometry,<:DistQuantity}; arrayproto::AbstractArray=zeros(), observed::Vector{Symbol}=Symbol[])
    pvar_arrays = OrderedDict()
    param_arrays = OrderedDict()
    layer_metas = OrderedDict()
    for (i,node) in enumerate(strat)
        # determine subgrid for layer
        lo = strat.boundaries[i]
        hi = (i < length(strat) ? strat.boundaries[i+1] : grid[end])
        # build subgrid using closed interval [lo,hi]
        subgrid = grid[lo..hi]
        # build layer
        prog_carr, param_carr, meta = _buildlayer(node,subgrid,arrayproto)
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
    nt_cache = NamedTuple{Tuple(nodenames)}(Tuple(_buildcaches(strat, nt_meta, arrayproto)))
    # construct prototype of u (prognostic state) array (note that this currently performs a copy)
    uproto = ComponentArray(nt_prog)
    # ditto for parameter array (need a hack here to get an empty ComponentArray...)
    pproto = sum(map(length, nt_params)) > 0 ? ComponentArray(nt_params) :
        ComponentArray(similar(arrayproto,0),(Axis{NamedTuple{Tuple(keys(nt_params))}(Tuple(map(a->1:0,nt_params)))}(),))
    # reconstruct with given array type
    uproto = ComponentArray(similar(arrayproto,length(uproto)), getaxes(uproto))
    CryoGridSetup(strat,grid,nt_meta,nt_cache,uproto,pproto; observed=observed)
end
CryoGridSetup(strat::Stratigraphy, grid::Grid{Cells}; kwargs...) = CryoGridSetup(strat, edges(grid); kwargs...)
CryoGridSetup(strat::Stratigraphy, grid::Grid{Edges,<:Numerics.Geometry,T}; kwargs...) where {T} = error("grid must have values with units of length, e.g. try using `Grid((x)u\"m\")` where `x` are your grid points.")

"""
    parameters(setup::CryoGridSetup; unconstrained=false)

Helper function to obtain the parameters of a `CryoGridSetup`. If `unconstrained=true`, the parameters will be mapped
to an unconstrained space first via `unconstrain`. Otherwise, they will be left as their default/initialized values.
"""
parameters(setup::CryoGridSetup; unconstrained=false) = unconstrained ? unconstrain(copy(setup.pproto), setup) : copy(setup.pproto)
"""
    withaxes(u::AbstractArray, ::CryoGridSetup)

Constructs a `ComponentArray` with labeled axes from the given state vector `u`. Assumes `u` to be of the same type/shape
as `setup.uproto`.
"""
withaxes(u::AbstractArray, setup::CryoGridSetup) = ComponentArray(u, getaxes(setup.uproto))
withaxes(u::ComponentArray, ::CryoGridSetup) = u
@generated function getstates(setup::CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,pax,names}, du::AbstractArray, u::AbstractArray, p::ComponentArray, t) where {TStrat,TGrid,TMeta,TCache,T,A,uax,pax,names}
    stategetters = Tuple((:(getstate($(QuoteNode(name)), setup, du, u, p, t)) for name in names))
    return :(NamedTuple{names}(tuple($(stategetters...))))
end
@generated function getstates(setup::CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,pax,names}, du::AbstractArray, u::AbstractArray, p::ComponentArray, t, ::Val{:diagnostic}) where {TStrat,TGrid,TMeta,TCache,T,A,uax,pax,names}
    function diagnosticvarnames(::Type{M}) where {M}
        pvars, dvars, avars, _ = _resolve_vartypes(M)
        # get diagnostic variables and derivatives
        return tuple(QuoteNode.(varname.(dvars))..., map(s -> QuoteNode(Symbol(:d,s)), varname.(pvars))..., map(s -> QuoteNode(Symbol(:d,s)), varname.(avars))...)
    end
    stategetters = map(name -> :(getstate($(QuoteNode(name)), setup, du, u, p, t)), names)
    varnames = map(T -> :(tuple($(diagnosticvarnames(T)...))), TMeta.parameters)
    expr = Expr(:block)
    for (vars,state,name) in zip(varnames,stategetters,names)
        push!(expr.args, :($name = NamedTuple{$vars}($state)))
    end
    push!(expr.args, :(NamedTuple{names}(tuple($(names...)))))
    return expr
end
"""
    getstate(layername::Symbol, setup::CryoGridSetup, du::AbstractArray, u::AbstractArray, p::ComponentArray, t)

Builds the state named tuple for `layername` given `setup` and state arrays.
"""
getstate(layername::Symbol, setup::CryoGridSetup, du::AbstractArray, u::AbstractArray, p::ComponentArray, t) =
    _buildstate(
        setup.cache[layername],
        setup.meta[layername],
        withaxes(u,setup)[layername],
        withaxes(du,setup)[layername],
        p[layername],
        t,
        setup.strat.boundaries[findfirst(n -> nodename(n) == layername, setup.strat.nodes)]
    )
"""
    getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator)

Builds the state named tuple for `layername` given an initialized integrator.
"""
function getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator)
    let setup = integrator.f.f;
        _buildstate(
            setup.cache[layername],
            setup.meta[layername],
            withaxes(integrator.u,setup)[layername],
            withaxes(get_du(integrator),setup)[layername],
            integrator.p[layername],
            integrator.t,
            setup.strat.boundaries[findfirst(n -> nodename(n) == layername, setup.strat.nodes)]
        )
    end
end
"""
    getvar(var::Symbol, integrator::SciMLBase.DEIntegrator)
"""
getvar(var::Symbol, integrator::SciMLBase.DEIntegrator) = getvar(Val{var}(), integrator.f.f, integrator.u)
"""
    getvar(var::Symbol, setup::CryoGridSetup, u)
"""
getvar(var::Symbol, setup::CryoGridSetup, u) = getvar(Val{var}(), setup, u)
"""
    getvar(::Val{var}, setup::CryoGridSetup{TStrat,<:Grid,TMeta}, _u) where {var,TStrat,TMeta}

Generated function that finds all layers containing variable `var::Symbol` and returns an `ArrayPartition` combining
them into a single contiguous array (allocation free).

e.g: `T = getvar(:T, setup, u)`
"""
@generated function getvar(::Val{var}, setup::CryoGridSetup{TStrat,<:Grid,TMeta}, _u) where {var,TStrat,TMeta}
    expr = Expr(:block)
    nodetyps = nodetypes(TStrat)
    matchedlayers = []
    push!(expr.args, :(u = ComponentArray(_u, getaxes(setup.uproto))))
    for (i,node) in enumerate(nodetyps)
        name = nodename(node)
        metatype = TMeta.parameters[i]
        # extract variable type information from metadata type
        ptypes, dtypes, atypes, _ = _resolve_vartypes(metatype)
        prognames = tuplejoin(varname.(ptypes), varname.(atypes))
        diagnames = varname.(dtypes)
        identifier = Symbol(name,:_,var)
        if var ∈ prognames
            @>> quote
            $identifier = u.$name.$var
            end push!(expr.args)
            push!(matchedlayers, identifier)
        elseif var ∈ diagnames
            @>> quote
            $identifier = retrieve(setup.cache.$name.$var, u)
            end push!(expr.args)
            push!(matchedlayers, identifier)
        end
    end
    @>> quote
    ArrayPartition($(matchedlayers...))
    end push!(expr.args)
    return expr
end
"""
    constrain(p::ComponentVector, setup::CryoGridSetup)

Invokes `constrain` on all parameters in `p`. Assumes `p` to be of the same form as `setup.pproto`.
"""
Numerics.constrain(p::ComponentVector, setup::CryoGridSetup) = _apply_or_unapply_constraints(:apply, p, setup)
"""
    unconstrain(p::ComponentVector, setup::CryoGridSetup)

Invokes `unconstrain` on all parameters in `p`. Assumes `p` to be of the same form as `setup.pproto`.
"""
Numerics.unconstrain(p::ComponentVector, setup::CryoGridSetup) = _apply_or_unapply_constraints(:unapply, p, setup)

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
@generated function (setup::CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,pax,names,obsv})(_du,_u,p,t) where {TStrat,TGrid,TMeta,TCache,T,A,uax,pax,names,obsv}
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
    if !(p <: ComponentArray)
        push!(expr.args, :(p = ComponentArray(p, getaxes(setup.pproto))))
    end
    # Initialize variables for all layers
    for i in 1:N
        n = nodename(nodetyps[i])
        nz = Symbol(n,:_z)
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        $nz = strat.boundaries[$i]
        $nstate = _buildstate(cache.$n, meta.$n, u.$n, du.$n, p.$n, t, $nz)
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
    # Write-back diagnostic variables to cache
    for i in 1:N
        n = nodename(nodetyps[i])
        nstate = Symbol(n,:state)
        # We have to really drill down into the TMeta named tuple type to extract the variable names...
        # TMeta is a Tuple{values...} sooo...
        vartyps = Tuple(
            TMeta.parameters[i]. # i'th layer metadata, NamedTuple
            parameters[2]. # value types, Tuple
            parameters[2]. # second value type (diagvars), Tuple
            parameters # Var types
        )
        # iterate over each variable, extract variable name, and copy it back to cache
        for var in vartyps
            nv = varname(var)
            @>> quote
            writeback!(cache.$n.$nv, $nstate.$nv)
            end push!(expr.args)
        end
    end
    # Observables
    for i in 1:N
        n = nodename(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        for name in obsv
            nameval = Val{name}()
            @>> quote
                observe($nameval,$nlayer,$nprocess,$nstate)
            end push!(expr.args)
        end
    end
    # make sure compiled method returns no value
    @>> :(return nothing) push!(expr.args)
    # emit generated expression block
    return expr
end

"""
    init!(setup::CryoGridSetup{TStrat}, p, tspan) where TStrat

Calls `initialcondition!` on all layers/processes and returns the fully constructed u0 and du0 states.
"""
@generated function init!(setup::CryoGridSetup{TStrat}, p, tspan) where TStrat
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
    for i in 1:N-1
        n1,n2 = nodename(nodetyps[i]), nodename(nodetyps[i+1])
        # create variable names
        n1z = Symbol(n1,:_z)
        n2z = Symbol(n2,:_z)
        n1,n2 = nodename(nodetyps[i]), nodename(nodetyps[i+1])
        n1state, n2state = Symbol(n1,:state), Symbol(n2,:state)
        n1layer, n2layer = Symbol(n1,:layer), Symbol(n2,:layer)
        n1process, n2process = Symbol(n1,:process), Symbol(n2,:process)
        # generated code for layer updates
        @>> quote
        $n1z = strat.boundaries[$i]
        $n2z = strat.boundaries[$(i+1)]
        $n1state = _buildstate(cache.$n1, meta.$n1, u.$n1, du.$n1, p.$n1, tspan[1], $n1z)
        $n2state = _buildstate(cache.$n2, meta.$n2, u.$n2, du.$n2, p.$n2, tspan[1], $n2z)
        $n1layer = strat.nodes[$i].layer
        $n1process = strat.nodes[$i].process
        $n2layer = strat.nodes[$(i+1)].layer
        $n2process = strat.nodes[$(i+1)].process
        end push!(expr.args)
        if i == 1
            # only invoke initialcondition! for layer i on first iteration to avoid duplicated calls
            @>> quote
            initialcondition!($n1layer,$n1state)
            initialcondition!($n1layer,$n1process,$n1state)
            end push!(expr.args)
        end
        # invoke initialcondition! for each layer, then for both (similar to interact!)
        @>> quote
        initialcondition!($n2layer,$n2state)
        initialcondition!($n2layer,$n2process,$n2state)
        initialcondition!($n1layer,$n1process,$n2layer,$n2process,$n1state,$n2state)
        end push!(expr.args)
    end
    @>> quote
    return getdata(u), getdata(du)
    end push!(expr.args)
    return expr
end

_apply_or_unapply_constraint(::Val{:apply}, p, x) = constrain(p,x)
_apply_or_unapply_constraint(::Val{:unapply}, p, x) = unconstrain(p,x)
function _apply_or_unapply_constraints(mode::Symbol, p::ComponentVector, setup::CryoGridSetup)
    pp = similar(p)
    for layer in keys(setup.meta)
        for param in setup.meta[layer].paramvars
            name = varname(param)
            p_k = @view p[layer]
            pp_k = @view pp[layer]
            pp_k[name] .= _apply_or_unapply_constraint(Val{mode}(), param, (@view p_k[name]))
        end
    end
    return pp
end

"""
Helper function to extract prognostic, diagnostic, algebraic, and parameter variable type
information from the `meta` field type signature.
"""
function _resolve_vartypes(::Type{M}) where {names,types,M <: NamedTuple{names,types}}
    # extract variables types from M; we assume the first two parameters in the NamedTuple
    # are the prognostic and diagnostic variable names respeictively. This must be respected by _buildlayer.
    # note that types.parameters[1] is Tuple{Var,...} so we call parameters again to get (Var,...)
    ptypes = types.parameters[1].parameters |> Tuple
    # again for diagnostic, assumed to be at position 2
    dtypes = types.parameters[2].parameters |> Tuple
    # and for algebraic, assumed to be at position 3
    atypes = types.parameters[3].parameters |> Tuple
    # and for parameters, assumed to be at position 4
    paramtypes = types.parameters[4].parameters |> Tuple
    return ptypes, dtypes, atypes, paramtypes
end

"""
Generates a function from layer cache and metadata which constructs a type-stable NamedTuple of state variables at runtime.
"""
@inline @generated function _buildstate(cache::NamedTuple, meta::M, u, du, params, t, z) where {M <: NamedTuple}
    ptypes, dtypes, atypes, _ = _resolve_vartypes(M)
    # here we join together prognostic and algebraic variables
    pnames = tuplejoin(ptypes .|> varname, atypes .|> varname)
    dnames = dtypes .|> varname
    # generate state variable accessor expressions
    pacc = tuple((:(u.$p) for p in pnames)...,)
    dacc = tuple((:(retrieve(cache.$d,u,t)) for d in dnames)...,)
    # construct symbols for derivative variables; assumes no existing conflicts
    dpnames = map(n -> Symbol(:d,n), pnames)
    dpacc = tuple((:(du.$p) for p in pnames)...,)
    # build state named tuple;
    # QuoteNode is used to force names to be interpolated as symbols rather than literals.
    quote
    NamedTuple{tuple($(map(QuoteNode,pnames)...),$(map(QuoteNode,dpnames)...),$(map(QuoteNode,dnames)...),
        :params,:grids,:t,:z)}(tuple($(pacc...),$(dpacc...),$(dacc...),params,meta.grids,t,z))
    end
end

"""
Constructs prognostic state vector and state named-tuple for the given node/layer.
"""
function _buildlayer(node::StratNode, grid::Grid{Edges}, arrayproto::A) where {A<:AbstractArray}
    layer, process = node.layer, node.process
    layer_vars = variables(layer)
    @assert all([isdiagnostic(var) || isparameter(var) for var in layer_vars]) "Layer variables must be diagnostic."
    process_vars = variables(layer, process)
    all_vars = tuple(layer_vars...,process_vars...)
    @debug "Building layer $(nodename(node)) with $(length(all_vars)) variables: $(all_vars)"
    # check for (permissible) duplicates between variables, excluding parameters
    groups = groupby(var -> varname(var), filter(x -> !isparameter(x), all_vars))
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
    # check for parameter duplicates
    param_dups = [name for (name,group) in groupby(param -> varname(param), param_vars) if length(group) > 1]
    @assert isempty(param_dups) "Found conflicting definitions of parameters: $(Tuple(param_dups))"
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
    diag_grids = _buildgrids(diag_vars, grid, arrayproto)
    prog_grids = _buildgrids(prog_vars, grid, arrayproto)
    alg_grids = _buildgrids(alg_vars, grid, arrayproto)
    param_grids = _buildgrids(param_vars, grid, arrayproto)
    # merge grids
    grids = merge(diag_grids, prog_grids, alg_grids, param_grids)
    # get variable names for diagnostic and prognostic
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
function _buildgrids(vars, grid::Grid{Edges}, arrayproto::A) where {A}
    if isempty(vars)
        return NamedTuple()
    end
    togrid(var::Var{name,T,OnGrid{Edges}}) where {name,T} = grid |> var.dim.f |> dustrip |> Grid
    togrid(var::Var{name,T,OnGrid{Cells}}) where {name,T} = grid |> cells |> var.dim.f |> dustrip |> Grid
    togrid(var::Var{name,T,Shape{dims}}) where {name,T,dims} = reshape(1:prod(dims),dims...)
    togrid(var::Var{name,T,typeof(Scalar)}) where {name,T} = 1:1
    names = @>> vars map(var -> varname(var))
    vars_with_names = NamedTuple{tuple(names...)}(tuple(vars...))
    grids = @>> vars_with_names map(togrid)
end

"""
Constructs per-layer variable caches given the Stratigraphy and layer-metadata named tuple.
"""
function _buildcaches(strat, metadata, arrayproto)
    map(strat) do node
        name = nodename(node)
        dvars = metadata[name].diagvars
        varnames = [varname(var) for var in dvars]
        caches = map(dvars) do dvar
            dvarname = varname(dvar)
            grid = metadata[name].grids[dvarname]
            VarCache(dvarname, grid, arrayproto)
        end
        NamedTuple{Tuple(varnames)}(Tuple(caches))
    end
end

struct VarCache{name,A}
    x::A
    function VarCache(name::Symbol, grid::AbstractArray, arrayproto::AbstractArray)
        A = similar(arrayproto, length(grid))
        A .= zero(eltype(A))
        new{name,typeof(A)}(A)
    end
end
_retrieve(c::VarCache, proto) = 0*proto .+ c.x
retrieve(c::VarCache) = c.x
# dispatches for autodiff types; create a 
retrieve(c::VarCache, u::AbstractArray{T}) where {T<:Union{<:ForwardDiff.Dual,<:ReverseDiff.TrackedReal}} = _retrieve(c, similar(u,length(c.x)))
retrieve(c::VarCache, u::ReverseDiff.TrackedArray) = _retrieve(c, similar(u,length(c.x)))
retrieve(c::VarCache, u::AbstractArray{T}) where {T} = retrieve(c)
retrieve(c::VarCache, u::AbstractArray, t) = retrieve(c, u)
# this covers the case for Rosenbrock solvers where only t has differentiable type
function retrieve(c::VarCache, u::AbstractArray, t::T) where {T<:ForwardDiff.Dual}
    proto = similar(u, T, length(c.x))
    return _retrieve(c, proto)
end
# default to doing nothing on non-autodiff writeback
writeback!(c::VarCache, x::AbstractArray) = nothing
writeback!(c::VarCache, x::AbstractArray{T}) where {T<:Union{ForwardDiff.Dual,ReverseDiff.TrackedReal}} = c.x .= Utils.adstrip.(x)
