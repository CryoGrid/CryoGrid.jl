mutable struct StateHistory
    vals::Union{Missing,<:Any}
    StateHistory() = new(missing)
end

"""
    AbstractTile{iip}

Base type for 1D tiles. `iip` is a value of enum `InPlaceMode` that indicates
whether the model operates on state variables in-place (overwriting arrays) or
out-of-place (copying arrays).
"""
abstract type AbstractTile{iip} end
"""
    (tile::AbstractTile{inplace})(du,u,p,t)
    (tile::AbstractTile{ooplace})(u,p,t)

Invokes the corresponding `step` function to compute the time derivative du/dt.
"""
(tile::AbstractTile{inplace})(du,u,p,t) = step!(tile,du,u,p,t)
(tile::AbstractTile{ooplace})(u,p,t) = step(tile,u,p,t)

"""
    step!(::T,du,u,p,t) where {T<:AbstractTile}

In-place step function for tile `T`. Computes du/dt and stores the result in `du`.
"""
step!(::T,du,u,p,t) where {T<:AbstractTile} = error("no implementation of in-place step! for $T")
"""
    step(::T,u,p,t) where {T<:AbstractTile}

Out-of-place step function for tile `T`. Computes and returns du/dt as vector with same size as `u`.
"""
step(::T,u,p,t) where {T<:AbstractTile} = error("no implementation of out-of-place step for $T")

"""
    Tile{TStrat,TGrid,TStates,TInits,TEvents,iip,obsv} <: AbstractTile{iip}

Defines the full specification of a single CryoGrid tile; i.e. stratigraphy, grid, and state variables.
"""
struct Tile{TStrat,TGrid,TStates,TInits,TEvents,iip,obsv} <: AbstractTile{iip}
    strat::TStrat # stratigraphy
    grid::TGrid # grid
    state::TStates # state variables
    inits::TInits # initializers
    events::TEvents # events
    hist::StateHistory # mutable "history" type for state tracking
    function Tile(
        strat::TStrat,
        grid::TGrid,
        state::TStates,
        inits::TInits,
        events::TEvents,
        hist::StateHistory=StateHistory(),
        iip::InPlaceMode=inplace,
        observe::Tuple{Vararg{Symbol}}=()) where
        {TStrat<:Stratigraphy,TGrid<:Grid{Edges},TStates<:VarStates,TInits<:Tuple,TEvents<:NamedTuple}
        new{TStrat,TGrid,TStates,TInits,TEvents,iip,observe}(strat,grid,state,inits,events,hist)
    end
end
ConstructionBase.constructorof(::Type{Tile{TStrat,TGrid,TStates,TInits,TEvents,iip,obsv}}) where {TStrat,TGrid,TStates,TInits,TEvents,iip,obsv} =
    (strat, grid, state, inits, events, hist) -> Tile(strat, grid, state, inits, events, hist, iip, obsv)
InputOutput.parameterize(tile::T) where {T<:Tile} = ConstructionBase.constructorof(T)(
    InputOutput.parameterize(tile.strat),
    tile.grid,
    tile.state,
    InputOutput.parameterize(tile.inits, layer=:init),
    InputOutput.parameterize(tile.events, layer=:event),
    tile.hist,
)
# mark only stratigraphy and initializers fields as flattenable
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{:strat}}) = true
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{:inits}}) = true
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{:events}}) = true
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{name}}) where name = false
Base.show(io::IO, ::MIME"text/plain", tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip,obsv}) where {TStrat,TGrid,TStates,TInits,TEvents,iip,obsv} = print(io, "Tile ($iip) with layers $(map(componentname, components(tile.strat))), observables=$obsv, $TGrid, $TStrat")

"""
    Tile(
        @nospecialize(strat::Stratigraphy),
        @nospecialize(grid::Grid{Edges,<:Numerics.Geometry,<:DistQuantity}),
        @nospecialize(inits::Numerics.VarInitializer...);
        arrayproto::Type{A}=Vector,
        iip::InPlaceMode=inplace,
        observe::Vector{Symbol}=Symbol[],
        chunksize=nothing,
    )

Constructs a `Tile` from the given stratigraphy and grid. `arrayproto` keyword arg should be an array instance
(of any arbitrary length, including zero, contents are ignored) that will determine the array type used for all state vectors.
"""
function Tile(
    @nospecialize(strat::Stratigraphy),
    @nospecialize(grid::Grid{Edges,<:Numerics.Geometry,<:DistQuantity}),
    @nospecialize(inits::Numerics.VarInitializer...);
    arrayproto::Type{A}=Vector,
    iip::InPlaceMode=inplace,
    observe::Vector{Symbol}=Symbol[],
    chunksize=nothing,
) where {A<:AbstractArray}
    vars = OrderedDict()
    events = OrderedDict()
    components = OrderedDict()
    for comp in stripunits(strat)
        name = componentname(comp)
        # build layer
        vars[name] = _collectvars(comp)
        components[name] = _addlayerfield(comp, name)
        # events
        evs = CryoGrid.events(comp.layer, comp.process)
        events[name] = _addlayerfield(evs, name)
    end
    inits = _addlayerfield(inits, :init)
    # rebuild stratigraphy with updated parameters
    strat = Stratigraphy(boundaries(strat), Tuple(values(components)))
    # construct state variables
    states = _initvarstates(strat, grid, vars, chunksize, arrayproto)
    if isempty(inits)
        @warn "No initializers provided. State variables without initializers will be set to zero by default."
    end
    return Tile(strat, grid, states, inits, (;events...), StateHistory(), iip, Tuple(observe))
end
Tile(strat::Stratigraphy, grid::Grid{Cells}; kwargs...) = Tile(strat, edges(grid); kwargs...)
Tile(strat::Stratigraphy, grid::Grid{Edges,<:Numerics.Geometry,T}; kwargs...) where {T} = error("grid must have values with units of length, e.g. try using `Grid((x)u\"m\")` where `x` are your grid points.")
"""
    step!(_tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,inplace,obsv}, _du, _u, p, t) where {TStrat,TGrid,TStates,TInits,TEvents,obsv}

Generated step function (i.e. du/dt) for any arbitrary Tile. Specialized code is generated and compiled
on the fly via the @generated macro to ensure type stability. The generated code updates each layer in the stratigraphy
in sequence, i.e for each layer 1 <= i <= N:

diagnosticstep!(layer i, ...)
interact!(layer i, ..., layer i+1, ...)
prognosticstep!(layer i, ...)

Note for developers: All sections of code wrapped in quote..end blocks are generated. Code outside of quote blocks
is only executed during compilation and will not appear in the compiled version.
"""
@generated function step!(_tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,inplace,obsv}, _du, _u, p, t) where {TStrat,TGrid,TStates,TInits,TEvents,obsv}
    nodetyps = componenttypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    quote
    _du .= zero(eltype(_du))
    du = ComponentArray(_du, getaxes(_tile.state.uproto))
    u = ComponentArray(_u, getaxes(_tile.state.uproto))
    tile = updateparams(_tile, u, p, t)
    strat = tile.strat
    state = TileState(tile.state, boundaries(strat), u, du, t, Val{inplace}())
    end |> Base.Fix1(push!, expr.args)
    # Generate identifiers for each component
    names = map(i -> componentname(nodetyps[i]), 1:N)
    comps = map(n -> Symbol(n,:comp), names)
    states = map(n -> Symbol(n,:state), names)
    layers = map(n -> Symbol(n,:layer), names)
    procs = map(n -> Symbol(n,:process), names)
    # Initialize variables for all layers
    for i in 1:N
        quote
        $(states[i]) = state.$(names[i])
        $(comps[i]) = components(strat)[$i]
        $(layers[i]) = $(comps[i]).layer
        $(procs[i]) = $(comps[i]).process
        end |> Base.Fix1(push!, expr.args)
    end
    # Diagnostic step
    for i in 1:N
        quote
        diagnosticstep!($(layers[i]), $(states[i]))
        diagnosticstep!($(layers[i]), $(procs[i]), $(states[i]))
        end |> Base.Fix1(push!, expr.args)
    end
    # Interact
    can_interact_exprs = map(i -> :(CryoGrid.thickness($(layers[i]), $(states[i])) > 0), 1:N)
    quote
    can_interact = tuple($(can_interact_exprs...))
    end |> Base.Fix1(push!, expr.args)
    # We only invoke interact! on pairs of layers for which the following are satisfied:
    # 1) Both layers have thickness > 0 (i.e. they occupy non-zero space in the stratigraphy)
    # 2) Both layers are adjacent, or more crudely, all layers in between (if any) have zero thickness;
    # In order to make this type stable, we pre-generate all possible interact! expressions in the
    # downward direction and add a check to each one.
    for i in 1:N-1
        for j in i+1:N
            if (i == 1 && j == 2) || (i == N-1 && j == N)
                # always apply top and bottom interactions
                quote
                interact!($(layers[i]),$(procs[i]),$(layers[j]),$(procs[j]),$(states[i]),$(states[j]))
                end
            else
                innerchecks = j > i+1 ? map(k -> :(can_interact[$k]), i+1:j-1) : :(false)
                quote
                if can_interact[$i] && can_interact[$j] && !any(tuple($(innerchecks...)))
                    interact!($(layers[i]),$(procs[i]),$(layers[j]),$(procs[j]),$(states[i]),$(states[j]))
                end
                end
            end |> Base.Fix1(push!, expr.args)
        end
    end
    # Prognostic step
    for i in 1:N
        quote
        prognosticstep!($(layers[i]),$(procs[i]),$(states[i]))
        end |> Base.Fix1(push!, expr.args)
    end
    # Observables
    for i in 1:N
        for name in obsv
            nameval = Val{name}()
            quote
                observe($nameval,$(layers[i]),$(procs[i]),$(states[i]))
            end |> Base.Fix1(push!, expr.args)
        end
    end
    # make sure compiled method returns no value
    push!(expr.args, :(return nothing))
    # emit generated expression block
    return expr
end
"""
    initialcondition!(tile::Tile, tspan::NTuple{2,Float64}, p::AbstractVector)
    initialcondition!(tile::Tile, tspan::NTuple{2,DateTime}, p::AbstractVector)

Calls `initialcondition!` on all layers/processes and returns the fully constructed u0 and du0 states.
"""
CryoGrid.initialcondition!(tile::Tile, tspan::NTuple{2,DateTime}, p::AbstractVector, args...) = initialcondition!(tile, convert_tspan(tspan), p)
@generated function CryoGrid.initialcondition!(tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip,obsv}, tspan::NTuple{2,Float64}, p::AbstractVector) where {TStrat,TGrid,TStates,TInits,TEvents,iip,obsv}
    nodetyps = componenttypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    @>> quote
    t0 = tspan[1]
    du = zero(similar(tile.state.uproto, eltype(p)))
    u = zero(similar(tile.state.uproto, eltype(p)))
    tile = updateparams(tile, u, p, t0)
    strat = tile.strat
    state = TileState(tile.state, boundaries(strat), u, du, t0, Val{iip}())
    end push!(expr.args)
    # Call initializers
    for i in 1:N
        for j in 1:length(TInits.parameters)
            @>> quote
            let layer = strat[$i].layer,
                layerstate = state[$i],
                init! = tile.inits[$j];
                if haskey(layerstate.states, varname(init!))
                    init!(layer, layerstate)
                end
            end
            end push!(expr.args)
        end
    end
    # Iterate over layers
    for i in 1:N-1
        n1,n2 = componentname(nodetyps[i]), componentname(nodetyps[i+1])
        # create variable names
        n1,n2 = componentname(nodetyps[i]), componentname(nodetyps[i+1])
        n1state, n2state = Symbol(n1,:state), Symbol(n2,:state)
        n1layer, n2layer = Symbol(n1,:layer), Symbol(n2,:layer)
        n1process, n2process = Symbol(n1,:process), Symbol(n2,:process)
        # generated code for layer updates
        @>> quote
        $n1state = state.$n1
        $n2state = state.$n2
        $n1layer = components(strat)[$i].layer
        $n1process = components(strat)[$i].process
        $n2layer = components(strat)[$(i+1)].layer
        $n2process = components(strat)[$(i+1)].process
        end push!(expr.args)
        # invoke initialcondition! for each layer, then for both (similar to interact!);
        # initialcondition! is called for layer1 only on the first iteration to avoid duplicate invocations
        if i == 1
            @>> quote
            initialcondition!($n1layer,$n1state)
            initialcondition!($n1layer,$n1process)
            initialcondition!($n1layer,$n1process,$n1state)
            end push!(expr.args)
        end
        @>> quote
        initialcondition!($n2layer,$n2state)
        initialcondition!($n2process,$n2state)
        initialcondition!($n2layer,$n2process,$n2state)
        initialcondition!($n1layer,$n1process,$n2layer,$n2process,$n1state,$n2state)
        end push!(expr.args)
    end
    @>> quote
    return u, du
    end push!(expr.args)
    return expr
end
@generated function CryoGrid.timestep(_tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip,obsv}, _du, _u, p, t) where {TStrat,TGrid,TStates,TInits,TEvents,iip,obsv}
    nodetyps = componenttypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    quote
    du = ComponentArray(_du, getaxes(_tile.state.uproto))
    u = ComponentArray(_u, getaxes(_tile.state.uproto))
    tile = updateparams(_tile, u, p, t)
    strat = tile.strat
    state = TileState(tile.state, boundaries(strat), u, du, t, Val{iip}())
    end |> Base.Fix1(push!, expr.args)
    # Generate identifiers for each component
    names = map(i -> componentname(nodetyps[i]), 1:N)
    comps = map(n -> Symbol(n,:comp), names)
    states = map(n -> Symbol(n,:state), names)
    layers = map(n -> Symbol(n,:layer), names)
    procs = map(n -> Symbol(n,:process), names)
    # Initialize variables for all layers
    for i in 1:N
        quote
        $(states[i]) = state.$(names[i])
        $(comps[i]) = components(strat)[$i]
        $(layers[i]) = $(comps[i]).layer
        $(procs[i]) = $(comps[i]).process
        end |> Base.Fix1(push!, expr.args)
    end
    push!(expr.args, :(dtmax = Inf))
    # Retrieve minimum timesteps
    for i in 1:N
        quote
        dtmax = min(dtmax, timestep($(layers[i]), $(procs[i]), $(states[i])))
        end |> Base.Fix1(push!, expr.args)
    end
    push!(expr.args, :(return dtmax))
    # emit generated expression block
    return expr
end
"""
    domain(tile::Tile)

Returns a function `isoutofdomain(u,p,t)` which checks whether any prognostic variable values
in `u` are outside of their respective domains. Only variables which have at least one finite
domain endpoint are checked; variables with unbounded domains are ignored.
"""
function domain(tile::Tile)
    # select only variables which have a finite domain
    vars = filter(x -> any(map(isfinite, extrema(vardomain(x)))), filter(isprognostic, variables(tile)))
    function isoutofdomain(_u,p,t)::Bool
        u = withaxes(_u, tile)
        for var in vars
            domain = vardomain(var)
            uvar = getproperty(u, varname(var))
            @inbounds for i in 1:length(uvar)
                if uvar[i] ∉ domain
                    return true
                end 
            end
        end
        return false
    end
end
"""
    getvar(name::Symbol, tile::Tile, u)
    getvar(::Val{name}, tile::Tile, u)

Retrieves the (diagnostic or prognostic) grid variable from `tile` given prognostic state `u`.
If `name` is not a variable in the tile, or if it is not a grid variable, `nothing` is returned.
"""
Numerics.getvar(name::Symbol, tile::Tile, u) = getvar(Val{name}(), tile, u)
Numerics.getvar(::Val{name}, tile::Tile, u) where name = getvar(Val{name}(), tile.state, withaxes(u, tile))
"""
    getstate(layername::Symbol, tile::Tile, u, du, t)
    getstate(::Val{layername}, tile::Tile{TStrat,TGrid,<:VarStates{layernames},iip}, _u, _du, t)

Constructs a `LayerState` representing the full state of `layername` given `tile`, state vectors `u` and `du`, and the
time step `t`.
"""
getstate(layername::Symbol, tile::Tile, u, du, t, dt=nothing) = getstate(Val{layername}(), tile, u, du, t, dt)
function getstate(::Val{layername}, tile::Tile{TStrat,TGrid,<:VarStates{layernames},TInits,TEvents,iip}, _u, _du, t, dt=nothing) where {layername,TStrat,TGrid,TInits,TEvents,iip,layernames}
    du = ComponentArray(_du, getaxes(tile.state.uproto))
    u = ComponentArray(_u, getaxes(tile.state.uproto))
    i = 1
    for j in 1:length(tile.strat)
        if layernames[j] == layername
            i = j
            break
        end
    end
    z = boundarypairs(map(ustrip, stripparams(boundaries(tile.strat))), ustrip(tile.grid[end]))[i]
    return LayerState(tile.state, z, u, du, t, dt, Val{layername}(), Val{iip}())
end
"""
    variables(tile::Tile)

Returns a tuple of all variables defined in the tile.
"""
CryoGrid.variables(tile::Tile) = Tuple(unique(Flatten.flatten(tile.state.vars, Flatten.flattenable, Var)))
"""
    parameters(tile::Tile; kwargs...)

Extracts all parameters from `tile` in a vector.
"""
parameters(tile::Tile; kwargs...) = CryoGridParams(tile; kwargs...)
"""
    withaxes(u::AbstractArray, ::Tile)

Constructs a `ComponentArray` with labeled axes from the given state vector `u`. Assumes `u` to be of the same type/shape
as `setup.uproto`.
"""
withaxes(u::AbstractArray, tile::Tile) = ComponentArray(u, getaxes(tile.state.uproto))
withaxes(u::ComponentArray, ::Tile) = u
function getstate(tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip}, _u, _du, t) where {TStrat,TGrid,TStates,TInits,TEvents,iip}
    du = ComponentArray(_du, getaxes(tile.state.uproto))
    u = ComponentArray(_u, getaxes(tile.state.uproto))
    return TileState(tile.state, map(ustrip ∘ stripparams, boundaries(tile.strat)), u, du, t, Val{iip}())
end
"""
    updateparams(tile::Tile, u, p, t)

Replaces all `ModelParameters.AbstractParam` values in `tile` with their (possibly updated) value from `p`.
Subsequently evaluates and replaces all nested `DynamicParameterization`s.
"""
function updateparams(tile::Tile{TStrat,TGrid,TStates}, u, p, t) where {TStrat,TGrid,TStates}
    # unfortunately, reconstruct causes allocations due to a mysterious dynamic dispatch when returning the result of _reconstruct;
    # I really don't know why, could be a compiler bug, but it doesn't happen if we call the internal _reconstruct directly soooo....
    tile_updated = Flatten._reconstruct(tile, p, Flatten.flattenable, ModelParameters.AbstractParam, Union{TGrid,TStates,StateHistory},1)[1]
    dynamic_ps = Flatten.flatten(tile_updated, Flatten.flattenable, DynamicParameterization, Union{TGrid,TStates,StateHistory})
    # TODO: perhaps should allow dependence on local layer state;
    # this would likely require per-layer deconstruction/reconstruction of `StratComponent`s in order to
    # build the `LayerState`s and evaluate the dynamic parameters in a fully type stable manner.
    dynamic_values = map(d -> d(u, t), dynamic_ps)
    return Flatten._reconstruct(tile_updated, dynamic_values, Flatten.flattenable, DynamicParameterization, Union{TGrid,TStates,StateHistory},1)[1]
end
"""
Collects and validates all declared variables (`Var`s) for the given stratigraphy component.
"""
function _collectvars(@nospecialize(comp::StratComponent))
    layer, process = comp.layer, comp.process
    declared_vars = variables(layer, process)
    nested_vars = Flatten.flatten(comp, Flatten.flattenable, Var)
    all_vars = tuplejoin(declared_vars, nested_vars)
    @debug "Building layer $(componentname(comp)) with $(length(all_vars)) variables: $(all_vars)"
    # check for (permissible) duplicates between variables, excluding parameters
    groups = groupby(var -> varname(var), all_vars)
    for (name,gvars) in filter(g -> length(g.second) > 1, groups)
        # if any duplicate variable deifnitions do not match, raise an error
        @assert all(gvars[i] == gvars[i-1] for i in 2:length(gvars)) "Found one or more conflicting definitions of $name in $gvars"
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
    # check for conflicting definitions of differential vars
    diff_varnames = map(v -> varname(Delta(v)), prog_alg)
    @assert all((isempty(filter(v -> varname(v) == d, all_vars)) for d in diff_varnames)) "Variable names $(Tuple(diff_varnames)) are reserved for differentials."
    # prognostic takes precedence, so we remove duplicated variables from the diagnostic variable set
    diag_vars = filter(v -> v ∉ diag_prog, diag_vars)
    # filter remaining duplicates
    diag_vars = unique(diag_vars)
    prog_vars = unique(prog_vars)
    alg_vars = unique(alg_vars)
    # convert back to tuples
    diag_vars, prog_vars, alg_vars = Tuple(diag_vars), Tuple(prog_vars), Tuple(alg_vars)
    return tuplejoin(diag_vars, prog_vars, alg_vars)
end
"""
Rebuilds the `obj` adding `name` to the `layer` field to all `Param`s, if any are defined.
"""
function _addlayerfield(@nospecialize(obj), name::Symbol)
    params = ModelParameters.params(obj)
    if length(params) > 0
        # create sub-model and add layer name to all parameters
        m = Model(obj)
        m[:layer] = repeat([name], length(params))
        return parent(m)
    else
        return obj
    end
end
"""
Initialize `VarStates` which holds the caches for all defined state variables.
"""
function _initvarstates(@nospecialize(strat::Stratigraphy), @nospecialize(grid::Grid), @nospecialize(vars::OrderedDict), chunksize::Union{Nothing,Int}, arrayproto::Type{A}) where {A}
    componentnames = [componentname(node) for node in strat]
    ntvars = NamedTuple{Tuple(componentnames)}(Tuple(values(vars)))
    npvars = (length(filter(isprognostic, layer)) + length(filter(isalgebraic, layer)) for layer in ntvars) |> sum
    ndvars = (length(filter(isdiagnostic, layer)) for layer in ntvars) |> sum
    @assert (npvars + ndvars) > 0 "No variable definitions found. Did you add a method definition for CryoGrid.variables(::L,::P) where {L<:Layer,P<:Process}?"
    @assert npvars > 0 "At least one prognostic variable must be specified."
    para = params(strat)
    chunksize = isnothing(chunksize) ? length(para) : chunksize
    states = VarStates(ntvars, Grid(dustrip(grid), grid.geometry), chunksize, arrayproto)
    return states
end
