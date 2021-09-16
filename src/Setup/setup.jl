mutable struct StateHistory
    vals::Union{Missing,<:SavedValues}
    StateHistory() = new(missing)
end

"""
    CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,names,obsv,P}

Defines the full specification of a CryoGrid model; i.e. stratigraphy, grids, variables, and diagnostic state. `uproto`
field is an uninitialized, prototype `ComponentArray` that holds the axis information for the prognostic state vector.
"""
struct CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,names,obsv,P}
    strat::TStrat   # stratigraphy
    grid::TGrid     # grid
    meta::NamedTuple{names,TMeta} # metadata (variable info and grids per layer)
    cache::NamedTuple{names,TCache} # variable caches (per layer)
    hist::StateHistory # mutable "history" type for state tracking
    uproto::ComponentVector{T,A,uax} # prototype prognostic state ComponentArray for integrator
    para::P
    function CryoGridSetup(
        strat::TStrat,
        grid::TGrid,
        meta::NamedTuple{names,TMeta},
        cache::NamedTuple{names,TCache},
        hist::StateHistory,
        uproto::ComponentVector{T,A,uax},
        para::P,
        observed::Vector{Symbol}=Symbol[]) where
        {TStrat<:Stratigraphy,TGrid<:Grid{Edges},TMeta<:Tuple,TCache<:Tuple,T<:Number,A<:AbstractVector{T},P,uax,names}
        new{TStrat,TGrid,TMeta,TCache,T,A,uax,names,tuple(observed...),P}(strat,grid,meta,cache,hist,uproto,para)
    end
end
ConstructionBase.constructorof(::Type{CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,names,obsv,P}}) where {TStrat,TGrid,TMeta,TCache,T,A,P,uax,names,obsv} =
    (strat, grid, meta, cache, hist, uproto, para) -> CryoGridSetup(strat,grid,meta,cache,hist,uproto,para,length(obsv) > 0 ? collect(obsv) : Symbol[])

"""
Constructs a `CryoGridSetup` from the given stratigraphy and grid. `arrayproto` keyword arg should be an array instance
(of any arbitrary length, including zero, contents are ignored) that will determine the array type used for all state vectors.
"""
function CryoGridSetup(
    @nospecialize(strat::Stratigraphy),
    @nospecialize(grid::Grid{Edges,<:Numerics.Geometry,<:DistQuantity});
    arrayproto::AbstractArray=zeros(),
    observed::Vector{Symbol}=Symbol[],
    chunksize=nothing)
    pvar_arrays = OrderedDict()
    layer_metas = OrderedDict()
    nodes = OrderedDict()
    for (i,node) in enumerate(strat)
        name = componentname(node)
        # determine subgrid for layer
        lo = strat.boundaries[i].val
        hi = (i < length(strat) ? strat.boundaries[i+1].val : grid[end])
        # build subgrid using closed interval [lo,hi]
        subgrid = grid[lo..hi]
        # build layer
        prog_carr, meta = _buildlayer(node,subgrid,arrayproto)
        pvar_arrays[name] = prog_carr
        layer_metas[name] = meta
        params = ModelParameters.params(node)
        if length(params) > 0
            # create sub-model and add layer name to all parameters
            submodel = Model(node)
            submodel[:layer] = repeat([name], length(params))
            nodes[name] = parent(submodel)
        else
            nodes[name] = node
        end
    end
    # rebuild stratigraphy with updated parameters
    strat = Stratigraphy(strat.boundaries, Tuple(values(nodes)))
    # construct named tuples containing data for each layer
    componentnames = [componentname(node) for node in strat]
    nt_prog = NamedTuple{Tuple(componentnames)}(Tuple(values(pvar_arrays)))
    nt_meta = NamedTuple{Tuple(componentnames)}(Tuple(values(layer_metas)))
    npvars = (length(meta.progvars) for meta in nt_meta) |> sum
    ndvars = (length(meta.diagvars) for meta in nt_meta) |> sum
    @assert (npvars + ndvars) > 0 "No variable definitions found. Did you add a method definition for CryoGrid.variables(::L,::P) where {L<:Layer,P<:Process}?"
    @assert npvars > 0 "At least one prognostic variable must be specified."
    nt_cache = NamedTuple{Tuple(componentnames)}(Tuple(_buildcaches(strat, nt_meta, arrayproto, chunksize)))
    # construct prototype of u (prognostic state) array (note that this currently performs a copy)
    uproto = ComponentArray(nt_prog)
    # reconstruct with given array type
    uproto = ComponentArray(similar(arrayproto,length(uproto)), getaxes(uproto))
    CryoGridSetup(strat,grid,nt_meta,nt_cache,StateHistory(),uproto,params(strat),observed)
end
CryoGridSetup(strat::Stratigraphy, grid::Grid{Cells}; kwargs...) = CryoGridSetup(strat, edges(grid); kwargs...)
CryoGridSetup(strat::Stratigraphy, grid::Grid{Edges,<:Numerics.Geometry,T}; kwargs...) where {T} = error("grid must have values with units of length, e.g. try using `Grid((x)u\"m\")` where `x` are your grid points.")

# mark only stratigraphy field as flattenable
Flatten.flattenable(::Type{<:CryoGridSetup}, ::Type{Val{:strat}}) = true
Flatten.flattenable(::Type{<:CryoGridSetup}, ::Type{Val{name}}) where name = false

"""
    withaxes(u::AbstractArray, ::CryoGridSetup)

Constructs a `ComponentArray` with labeled axes from the given state vector `u`. Assumes `u` to be of the same type/shape
as `setup.uproto`.
"""
withaxes(u::AbstractArray, setup::CryoGridSetup) = ComponentArray(u, getaxes(setup.uproto))
withaxes(u::ComponentArray, ::CryoGridSetup) = u
@generated function getstates(setup::CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,names}, du::AbstractArray, u::AbstractArray, t) where {TStrat,TGrid,TMeta,TCache,T,A,uax,names}
    stategetters = Tuple((:(getstate($(QuoteNode(name)), setup, du, u, t)) for name in names))
    return :(NamedTuple{names}(tuple($(stategetters...))))
end
@generated function getstates(setup::CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,names}, du::AbstractArray, u::AbstractArray, t, ::Val{:diagnostic}) where {TStrat,TGrid,TMeta,TCache,T,A,uax,names}
    function diagnosticvarnames(::Type{M}) where {M}
        pvars, dvars, avars = _resolve_vartypes(M)
        # get diagnostic variables and derivatives
        return tuple(QuoteNode.(varname.(dvars))..., map(s -> QuoteNode(Symbol(:d,s)), varname.(pvars))..., map(s -> QuoteNode(Symbol(:d,s)), varname.(avars))...)
    end
    stategetters = map(name -> :(getstate($(QuoteNode(name)), setup, du, u, t)), names)
    varnames = map(T -> :(tuple($(diagnosticvarnames(T)...))), TMeta.parameters)
    expr = Expr(:block)
    for (vars,state,name) in zip(varnames,stategetters,names)
        push!(expr.args, :($name = NamedTuple{$vars}($state)))
    end
    push!(expr.args, :(NamedTuple{names}(tuple($(names...)))))
    return expr
end
"""
    getstate(layername::Symbol, setup::CryoGridSetup, du::AbstractArray, u::AbstractArray, t)

Builds the state named tuple for `layername` given `setup` and state arrays.
"""
getstate(layername::Symbol, setup::CryoGridSetup, du::AbstractArray, u::AbstractArray, t) = getstate(Val{layername}(), setup, du, u, t)
@generated function getstate(::Val{layername}, setup::CryoGridSetup{TStrat}, du::AbstractArray, u::AbstractArray, t) where {TStrat,layername}
    names = map(componentname, componenttypes(TStrat))
    i = findfirst(n -> n == layername, names)
    quote
        _buildstate(
            setup.cache[$(QuoteNode(layername))],
            setup.meta[$(QuoteNode(layername))],
            withaxes(u,setup).$layername,
            withaxes(du,setup).$layername,
            t,
            setup.strat.boundaries[$i]
        )
    end
end
"""
    getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator)

Builds the state named tuple for `layername` given an initialized integrator.
"""
getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator) = getstate(Val{layername}(), integrator)
@generated function getstate(::Val{layername}, integrator::SciMLBase.DEIntegrator) where {layername}
    # a bit hacky and may break in the future... but this is the hardcoded position of the CryoGridSetup type in DEIntegrator
    TStrat = integrator.parameters[13].parameters[2].parameters[1]
    names = map(componentname, componenttypes(TStrat))
    i = findfirst(n -> n == layername, names)
    quote
        let setup = integrator.f.f;
            _buildstate(
                setup.cache[$(QuoteNode(layername))],
                setup.meta[$(QuoteNode(layername))],
                withaxes(integrator.u,setup).$layername,
                withaxes(get_du(integrator),setup).$layername,
                integrator.t,
                setup.strat.boundaries[$i]
            )
        end
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
    nodetyps = componenttypes(TStrat)
    matchedlayers = []
    push!(expr.args, :(u = ComponentArray(_u, getaxes(setup.uproto))))
    for (i,node) in enumerate(nodetyps)
        name = componentname(node)
        metatype = TMeta.parameters[i]
        # extract variable type information from metadata type
        ptypes, dtypes, atypes = _resolve_vartypes(metatype)
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
Generated step function (i.e. du/dt) for any arbitrary CryoGridSetup. Specialized code is generated and compiled
on the fly via the @generated macro to ensure type stability. The generated code updates each layer in the stratigraphy
in sequence, i.e for each layer 1 < i < N:

diagnosticstep!(layer i, ...)
interact!(layer i-1, ...)
prognosticstep!(layer i, ...)

Note for developers: All sections of code wrapped in quote..end blocks are generated. Code outside of quote blocks
is only executed during compilation and will not appear in the compiled version.
"""
@generated function (setup::CryoGridSetup{TStrat,TGrid,TMeta,TCache,T,A,uax,names,obsv})(_du,_u,_p,t) where {TStrat,TGrid,TMeta,TCache,T,A,uax,names,obsv}
    nodetyps = componenttypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    @>> quote
    p = updateparams!(_p, setup, _du, _u, t)
    strat = ModelParameters._update(setup.para, setup.strat, p)
    cache = setup.cache
    meta = setup.meta
    _du .= zero(eltype(_du))
    du = ComponentArray(_du, getaxes(setup.uproto))
    u = ComponentArray(_u, getaxes(setup.uproto))
    end push!(expr.args)
    # Initialize variables for all layers
    for i in 1:N
        n = componentname(nodetyps[i])
        nz = Symbol(n,:_z)
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        $nz = strat.boundaries[$i]
        $nstate = _buildstate(cache.$n, meta.$n, u.$n, du.$n, t, $nz)
        $nlayer = strat.components[$i].layer
        $nprocess = strat.components[$i].process
        end push!(expr.args)
    end
    # Diagnostic step
    for i in 1:N
        n = componentname(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        diagnosticstep!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    # Interact
    for i in 1:N-1
        n1,n2 = componentname(nodetyps[i]), componentname(nodetyps[i+1])
        n1state, n2state = Symbol(n1,:state), Symbol(n2,:state)
        n1layer, n2layer = Symbol(n1,:layer), Symbol(n2,:layer)
        n1process, n2process = Symbol(n1,:process), Symbol(n2,:process)
        @>> quote
        interact!($n1layer,$n1process,$n2layer,$n2process,$n1state,$n2state)
        end push!(expr.args)
    end
    # Prognostic step
    for i in 1:N
        n = componentname(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        prognosticstep!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    # Write-back diagnostic variables to cache
    for i in 1:N
        n = componentname(nodetyps[i])
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
        n = componentname(nodetyps[i])
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
@generated function init!(setup::CryoGridSetup{TStrat}, tspan, p) where TStrat
    nodetyps = componenttypes(TStrat)
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
        n1,n2 = componentname(nodetyps[i]), componentname(nodetyps[i+1])
        # create variable names
        n1z = Symbol(n1,:_z)
        n2z = Symbol(n2,:_z)
        n1,n2 = componentname(nodetyps[i]), componentname(nodetyps[i+1])
        n1state, n2state = Symbol(n1,:state), Symbol(n2,:state)
        n1layer, n2layer = Symbol(n1,:layer), Symbol(n2,:layer)
        n1process, n2process = Symbol(n1,:process), Symbol(n2,:process)
        # generated code for layer updates
        @>> quote
        $n1z = strat.boundaries[$i]
        $n2z = strat.boundaries[$(i+1)]
        $n1state = _buildstate(cache.$n1, meta.$n1, u.$n1, du.$n1, tspan[1], $n1z)
        $n2state = _buildstate(cache.$n2, meta.$n2, u.$n2, du.$n2, tspan[1], $n2z)
        $n1layer = strat.components[$i].layer
        $n1process = strat.components[$i].process
        $n2layer = strat.components[$(i+1)].layer
        $n2process = strat.components[$(i+1)].process
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
    return ptypes, dtypes, atypes
end

"""
Generates a function from layer cache and metadata which constructs a type-stable NamedTuple of state variables at runtime.
"""
@inline @generated function _buildstate(cache::NamedTuple, meta::M, u, du, t, z) where {M <: NamedTuple}
    ptypes, dtypes, atypes = _resolve_vartypes(M)
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
        :grids,:t,:z)}(tuple($(pacc...),$(dpacc...),$(dacc...),meta.grids,t,z))
    end
end

"""
Constructs prognostic state vector and state named-tuple for the given node/layer.
"""
function _buildlayer(@nospecialize(node::StratComponent), @nospecialize(grid::Grid{Edges}), arrayproto::A) where {A<:AbstractArray}
    layer, process = node.layer, node.process
    layer_vars = variables(layer)
    @assert all([isdiagnostic(var) for var in layer_vars]) "Layer variables must be diagnostic."
    process_vars = variables(layer, process)
    all_vars = tuple(layer_vars...,process_vars...)
    @debug "Building layer $(componentname(node)) with $(length(all_vars)) variables: $(all_vars)"
    # check for (permissible) duplicates between variables, excluding parameters
    groups = groupby(var -> varname(var), all_vars)
    for (name,gvars) in filter(g -> length(g.second) > 1, groups)
        # if any duplicate variable deifnitions do not match, raise an error
        @assert reduce(==, gvars) "Found one or more conflicting definitions of $name in $gvars"
    end
    diag_vars = filter(isdiagnostic, all_vars)
    prog_vars = filter(isprognostic, all_vars)
    alg_vars = filter(isalgebraic, all_vars)
    # check for duplicated algebraic/prognostic vars
    prog_alg_duplicated = prog_vars ∩ alg_vars
    @assert isempty(prog_alg_duplicated) "Variables $(prog_alg_duplicated) cannot be both prognostic and algebraic."
    # check for re-definition of diagnostic variables as prognostic
    prog_alg = prog_vars ∪ alg_vars
    diag_prog = filter(v -> v ∈ prog_alg, diag_vars)
    if !isempty(diag_prog)
        @warn "Variables $(Tuple(map(varname,diag_prog))) declared as both prognostic/algebraic and diagnostic. In-place modifications outside of callbacks may degrade integration accuracy."
    end
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
    # merge grids
    grids = merge(diag_grids, prog_grids, alg_grids)
    # get variable names for diagnostic and prognostic
    # build metadata named tuple
    layermetadata = NamedTuple{tuple(:progvars,:diagvars,:algvars,:grids)}(tuple(prog_vars,diag_vars,alg_vars,grids))
    # build component arrays for prognostic/algebraic variables and parameters
    prog_carr = isempty(prog_alg) ? similar(arrayproto, 0) : ComponentArray(merge(prog_grids, alg_grids))
    return prog_carr, layermetadata
end

"""
Constructs grid tuples for the given variables.
"""
function _buildgrids(@nospecialize(vars), @nospecialize(grid::Grid{Edges}), arrayproto::A) where {A}
    if isempty(vars)
        return NamedTuple()
    end
    togrid(var::Var{name,T,OnGrid{Edges}}) where {name,T} = grid |> var.dim.f |> dustrip |> Grid
    togrid(var::Var{name,T,OnGrid{Cells}}) where {name,T} = grid |> cells |> var.dim.f |> dustrip |> Grid
    togrid(var::Var{name,T,Shape{dims}}) where {name,T,dims} = reshape(1:prod(dims),dims...)
    togrid(var::Var{name,T,typeof(Scalar)}) where {name,T} = 1:1
    names = map(var -> varname(var), vars)
    vars_with_names = NamedTuple{tuple(names...)}(tuple(vars...))
    grids = map(togrid, vars_with_names)
end

"""
Constructs per-layer variable caches given the Stratigraphy and layer-metadata named tuple.
"""
function _buildcaches(@nospecialize(strat::Stratigraphy), @nospecialize(metadata::NamedTuple), arrayproto::AbstractArray, chunksize=nothing)
    chunksize = isnothing(chunksize) ? ForwardDiff.DEFAULT_CHUNK_THRESHOLD : chunksize
    map(strat) do node
        name = componentname(node)
        dvars = metadata[name].diagvars
        varnames = [varname(var) for var in dvars]
        caches = map(dvars) do dvar
            dvarname = varname(dvar)
            grid = metadata[name].grids[dvarname]
            VarCache(dvarname, grid, arrayproto, chunksize)
        end
        NamedTuple{Tuple(varnames)}(Tuple(caches))
    end
end

"""
    VarCache{N,A,Adual}

Wrapper for `DiffEqBase.DiffCache` that stores state variables in forward-diff compatible cache arrays.
"""
struct VarCache{N,A,Adual}
    name::Symbol
    cache::PreallocationTools.DiffCache{A,Adual}
    function VarCache(name::Symbol, grid::AbstractArray, arrayproto::AbstractArray, chunksize::Int)
        # use dual cache for automatic compatibility with ForwardDiff
        A = similar(arrayproto, length(grid))
        A .= zero(eltype(A))
        cache = PreallocationTools.dualcache(A, Val{chunksize})
        new{chunksize,typeof(cache.du),typeof(cache.dual_du)}(name, cache)
    end
end
Base.show(io::IO, cache::VarCache) = print(io, "VarCache $(cache.name) of length $(length(cache.cache.du)) with eltype $(eltype(cache.cache.du))")
Base.show(io::IO, mime::MIME{Symbol("text/plain")}, cache::VarCache) = show(io, cache)
# use pre-cached array if chunk size matches
retrieve(varcache::VarCache{N}, u::AbstractArray{T}) where {tag,U,N,T<:ForwardDiff.Dual{tag,U,N}} = DiffEqBase.get_tmp(varcache.cache, u)
# otherwise just make a new copy with compatible type
retrieve(varcache::VarCache, u::AbstractArray{T}) where {T<:Union{<:ForwardDiff.Dual,<:ReverseDiff.TrackedReal}} = copyto!(similar(u, length(varcache.cache.du)), varcache.cache.du)
retrieve(varcache::VarCache, u::ReverseDiff.TrackedArray) = copyto!(similar(identity.(u), length(varcache.cache.du)), varcache.cache.du)
retrieve(varcache::VarCache, u::AbstractArray{T}) where {T} = reinterpret(T, varcache.cache.du)
# this covers the case for Rosenbrock solvers where only t has differentiable type
retrieve(varcache::VarCache, u::AbstractArray, t::T) where {T<:ForwardDiff.Dual} = retrieve(varcache, similar(u, T))
retrieve(varcache::VarCache, u::AbstractArray, t) = retrieve(varcache, u)
retrieve(varcache::VarCache) = diffcache.du
# default to doing nothing on non-autodiff writeback
writeback!(varcache::VarCache, x::AbstractArray) = nothing
writeback!(varcache::VarCache, x::AbstractArray{T}) where {T<:Union{<:ForwardDiff.Dual,<:ReverseDiff.TrackedReal}} = varcache.cache.du .= Utils.adstrip.(x)
