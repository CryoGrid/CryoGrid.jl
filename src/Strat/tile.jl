mutable struct StateHistory
    vals::Union{Missing,<:Any}
    StateHistory() = new(missing)
end

"""
    AbstractTile{iip}

Base type for 1D tiles. `iip` is a boolean value that indicates, if true,
whether the model operates on state variables in-place (overwriting arrays) or
if false, out-of-place (copying arrays).
"""
abstract type AbstractTile{iip} end
"""
    (tile::AbstractTile{true})(du,u,p,t)
    (tile::AbstractTile{false})(u,p,t)

Invokes the corresponding `step` function to compute the time derivative du/dt.
"""
(tile::AbstractTile{true})(du,u,p,t) = step!(tile,du,u,p,t)
(tile::AbstractTile{false})(u,p,t) = step(tile,u,p,t)

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
    Tile{TStrat,TGrid,TStates,TInits,TEvents,iip} <: AbstractTile{iip}

Defines the full specification of a single CryoGrid tile; i.e. stratigraphy, grid, and state variables.
"""
struct Tile{TStrat,TGrid,TStates,TInits,TEvents,iip} <: AbstractTile{iip}
    strat::TStrat # stratigraphy
    grid::TGrid # grid
    state::TStates # state variables
    inits::TInits # initializers
    events::TEvents # events
    metadata::Dict # metadata
    hist::StateHistory # mutable "history" type for state tracking
    function Tile(
        strat::TStrat,
        grid::TGrid,
        state::TStates,
        inits::TInits,
        events::TEvents,
        metadata::Dict=Dict(),
        hist::StateHistory=StateHistory(),
        iip::Bool=true) where
        {TStrat<:Stratigraphy,TGrid<:Grid{Edges},TStates<:StateVars,TInits<:Tuple,TEvents<:NamedTuple}
        new{TStrat,TGrid,TStates,TInits,TEvents,iip}(strat,grid,state,inits,events,metadata,hist)
    end
end
ConstructionBase.constructorof(::Type{Tile{TStrat,TGrid,TStates,TInits,TEvents,iip}}) where {TStrat,TGrid,TStates,TInits,TEvents,iip} =
    (strat, grid, state, inits, events, metadata, hist) -> Tile(strat, grid, state, inits, events, metadata, hist, iip)
# mark only stratigraphy and initializers fields as flattenable
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{:strat}}) = true
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{:inits}}) = true
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{:events}}) = true
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{name}}) where name = false
Base.show(io::IO, ::MIME"text/plain", tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip}) where {TStrat,TGrid,TStates,TInits,TEvents,iip} = print(io, "Tile (iip=$iip) with layers $(map(layername, layers(tile.strat))), $TGrid, $TStrat")

"""
    Tile(
        @nospecialize(strat::Stratigraphy),
        @nospecialize(discretization_strategy::DiscretizationStrategy),
        @nospecialize(inits::VarInitializer...);
        metadata::Dict=Dict(),
        arrayproto::Type{A}=Vector,
        iip::Bool=true,
        chunk_size=nothing,
    )

Constructs a `Tile` from the given stratigraphy and discretization strategy. `arrayproto` keyword arg should be an array instance
(of any arbitrary length, including zero, contents are ignored) that will determine the array type used for all state vectors.
"""
function Tile(
    @nospecialize(strat::Stratigraphy),
    @nospecialize(discretization_strategy::DiscretizationStrategy),
    @nospecialize(inits::VarInitializer...);
    metadata::Dict=Dict(),
    arrayproto::Type{A}=Vector,
    iip::Bool=true,
    chunk_size=nothing,
) where {A<:AbstractArray}
    grid = Numerics.makegrid(strat, discretization_strategy)
    strat = stripunits(strat)
    events = CryoGrid.events(strat)
    vars = CryoGrid.variables(strat)
    layers = map(NamedTuple(strat)) do named_layer
        _addlayerfield(named_layer, nameof(named_layer))
    end
    # rebuild stratigraphy with updated parameters
    strat = Stratigraphy(boundaries(strat), Tuple(values(layers)))
    # construct state variables
    states = _initstatevars(strat, grid, vars, chunk_size, arrayproto)
    if isempty(inits)
        @warn "No initializers provided. State variables without initializers will be set to zero by default."
    end
    inits = map(inits) do init
        _addlayerfield(init, Symbol(:init_, varname(init)))
    end
    return Tile(strat, grid, states, inits, (;events...), metadata, StateHistory(), iip)
end
Tile(strat::Stratigraphy, grid::Grid{Cells}, inits...; kwargs...) = Tile(strat, edges(grid), inits...; kwargs...)
Tile(strat::Stratigraphy, grid::Grid{Edges}, inits...; kwargs...) = Tile(strat, PresetGrid(grid), inits...; kwargs...)
Tile(strat::Stratigraphy, inits...; discretization_strategy=AutoGrid(), kwargs...) = Tile(strat, discretization_strategy, inits...; kwargs...)

"""
    step!(_tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,true}, _du, _u, p, t) where {TStrat,TGrid,TStates,TInits,TEvents}

Time derivative step function (i.e. du/dt) for any arbitrary Tile. Specialized code is generated and compiled
on the fly via the @generated macro to ensure type stability. The generated code updates each layer in the stratigraphy
in sequence, i.e for each layer 1 <= i <= N:

```julia
updatestate!(layer[i], ...)
interact!(layer[i], ..., layer[i+1], ...)
computefluxes!(layer[i], ...)
```
"""
function step!(
    _tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,true},
    _du,
    _u,
    p,
    t,
    dt=1.0,
) where {N,TStrat<:Stratigraphy{N},TGrid,TStates,TInits,TEvents}
    _du .= zero(eltype(_du))
    du = ComponentArray(_du, getaxes(_tile.state.uproto))
    u = ComponentArray(_u, getaxes(_tile.state.uproto))
    tile = resolve(_tile, u, p, t)
    strat = tile.strat
    zs = map(getscalar, boundaries!(tile, u))
    # zs = boundaries(strat)
    state = TileState(tile.state, zs, u, du, t, dt, Val{true}())
    CryoGrid.updatestate!(strat, state)
    checkstate!(tile, state, u, du, :updatestate!)
    CryoGrid.interact!(strat, state)
    checkstate!(tile, state, u, du, :interact!)
    CryoGrid.computefluxes!(strat, state)
    checkstate!(tile, state, u, du, :computefluxes!)
    return nothing
end

"""
    timestep(_tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip}, _du, _u, p, t) where {TStrat,TGrid,TStates,TInits,TEvents,iip}

Computes the maximum permissible forward timestep for this `Tile` given the current `u`, `p`, and `t`.
"""
function CryoGrid.timestep(_tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip}, _du, _u, p, t) where {TStrat,TGrid,TStates,TInits,TEvents,iip}
    du = ComponentArray(_du, getaxes(_tile.state.uproto))
    u = ComponentArray(_u, getaxes(_tile.state.uproto))
    tile = resolve(_tile, u, p, t)
    strat = tile.strat
    zs = map(getscalar, boundaries!(tile, u))
    state = TileState(tile.state, zs, u, du, t, 1.0, Val{true}())
    CryoGrid.timestep(strat::Stratigraphy, state)
end

"""
    initialcondition!(tile::Tile, tspan::NTuple{2,Float64}, p=nothing)
    initialcondition!(tile::Tile, tspan::NTuple{2,DateTime}, p=nothing)

Calls `initialcondition!` on all layers/processes and returns the fully constructed u0 and du0 states.
"""
CryoGrid.initialcondition!(tile::Tile, tspan::NTuple{2,DateTime}, p=nothing, args...) = initialcondition!(tile, convert_tspan(tspan), p)
function CryoGrid.initialcondition!(tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip}, tspan::NTuple{2,Float64}, p=nothing) where {TStrat,TGrid,TStates,TInits,TEvents,iip}
    t0 = tspan[1]
    # if there are parameters defined on `tile` but the user did not supply a parameter vector,
    # automatically extract the parameters from the Tile.
    tile_params = isempty(ModelParameters.params(tile)) ? [] : parameters(tile)
    p = isnothing(p) && !isempty(tile_params) ? tile_params : p
    # choose type for state vectors
    u_type = isnothing(p) ? eltype(tile.state.uproto) : eltype(p)
    du = zero(similar(tile.state.uproto, u_type))
    u = zero(similar(tile.state.uproto, u_type))
    tile = resolve(tile, u, p, t0)
    strat = tile.strat
    # get stratigraphy boundaries
    zs = initboundaries!(tile, u)
    state = TileState(tile.state, zs, u, du, t0, 1.0, Val{iip}())
    CryoGrid.initialcondition!(strat, state, tile.inits)
    return u, du    
end

"""
    domain(tile::Tile)

Returns a function `isoutofdomain(u,p,t)` which checks whether any prognostic variable values
in `u` are outside of their respective domains. Only variables which have at least one finite
domain endpoint are checked; variables with unbounded domains are ignored.
"""
function domain(tile::Tile)
    # select only variables which have a finite domain
    vars = map(variables(tile.strat)) do vars
        filter(var -> any(map(isfinite, extrema(vardomain(var)))) && isprognostic(var), vars)
    end
    gridvars = filter(isongrid, tuplejoin(vars...))
    function isoutofdomain(_u,p,t)::Bool
        u = withaxes(_u, tile)
        for var in gridvars
            domain = vardomain(var)
            uvar = getproperty(u, varname(var))
            @inbounds for i in 1:length(uvar)
                if uvar[i] ∉ domain
                    return true
                end 
            end
        end
        for layer in keys(vars)
            layer_vars = filter(!isongrid, getproperty(vars, layer))
            for var in layer_vars
                domain = vardomain(var)
                uvar = getproperty(getproperty(u, layer), varname(var))
                @inbounds for i in 1:length(uvar)
                    if uvar[i] ∉ domain
                        return true
                    end 
                end
            end
        end
        return false
    end
end

# layer thickness
_layerthick(::PrognosticVolume, tile::Tile, ::Named{name,TLayer}, u) where {name,TLayer<:Layer} = getproperty(u, name).Δz
_layerthick(::Union{FixedVolume,DiagnosticVolume}, tile::Tile, ::Named{name,TLayer}, u) where {name,TLayer<:Layer} = retrieve(getproperty(tile.state.diag, name).Δz, u)
function _layerthick(tile::Tile, layer::Named{name,TLayer}, u) where {name,TLayer<:Layer}
    return _layerthick(CryoGrid.Volume(TLayer), tile, layer, u)
end

@generated function boundaries!(tile::Tile{TStrat}, u) where {TStrat}
    if all(map(typ -> isa(Volume(typ), FixedVolume), layertypes(TStrat)))
        # Micro-optimization: if all layers in the stratigraphy have static volume, skip all of the fancy stuff
        quote
            boundaries(tile.strat)
        end
    else
        quote
            zs = update_layer_boundaries(tile, u)
            return reverse(zs)
        end
    end
end

function initboundaries!(tile::Tile{TStrat}, u) where {TStrat}
    bounds = boundarypairs(tile.strat)
    map(bounds, layers(tile.strat)) do (z1, z2), named_layer
        name = nameof(named_layer)
        diag_layer = getproperty(tile.state.diag, name)
        z = retrieve(diag_layer.z, u)
        Δz = _layerthick(tile, named_layer, u)
        @setscalar Δz = z2 - z1
        @setscalar z = z1
        return z1
    end
end

function update_layer_boundaries(tile::Tile, u)
    # calculate grid boundaries starting from the bottom moving up to the surface
    zbot = tile.state.grid[end]
    return accumulate(reverse(layers(tile.strat)); init=zbot) do z_acc, named_layer
        name = nameof(named_layer)
        diag_layer = getproperty(tile.state.diag, name)
        z_state = retrieve(diag_layer.z, u)
        Δz = getscalar(_layerthick(tile, named_layer, u))
        @setscalar z_state = z_acc - max(Δz, zero(Δz))
        z = getscalar(z_state)
        # strip ForwardDiff type if necessary and round to avoid numerical issues
        return round(Numerics.ForwardDiff.value(z), digits=12)
    end
end

"""
    getvar(name::Symbol, tile::Tile, u; interp=true)
    getvar(::Val{name}, tile::Tile, u; interp=true)

Retrieves the (diagnostic or prognostic) grid variable from `tile` given prognostic state `u`.
If `name` is not a variable in the tile, or if it is not a grid variable, `nothing` is returned.
"""
Numerics.getvar(name::Symbol, tile::Tile, u; interp=true) = getvar(Val{name}(), tile, u; interp)
function Numerics.getvar(::Val{name}, tile::Tile, u; interp=true) where name
    x = getvar(Val{name}(), tile.state, withaxes(u, tile))
    @assert !isnothing(x) "no variable exists with name $name"
    if interp && length(x) == length(tile.grid)
        return Numerics.Interpolations.interpolate((edges(tile.grid),), x, Numerics.Interpolations.Gridded(Numerics.Interpolations.Linear()))
    elseif interp && length(x) == length(tile.grid)-1
        return Numerics.Interpolations.interpolate((cells(tile.grid),), x, Numerics.Interpolations.Gridded(Numerics.Interpolations.Linear()))
    else
        return x
    end
end

"""
    getstate(layername::Symbol, tile::Tile, u, du, t)
    getstate(::Val{layername}, tile::Tile{TStrat,TGrid,<:StateVars{layernames},iip}, _u, _du, t)

Constructs a `LayerState` representing the full state of `layername` given `tile`, state vectors `u` and `du`, and the
time step `t`.
"""
getstate(layername::Symbol, tile::Tile, u, du, t, dt=1.0) = getstate(Val{layername}(), tile, u, du, t, dt)
function getstate(::Val{layername}, tile::Tile{TStrat,TGrid,<:StateVars{layernames},TInits,TEvents,iip}, _u, _du, t, dt=1.0) where {layername,TStrat,TGrid,TInits,TEvents,iip,layernames}
    du = ComponentArray(_du, getaxes(tile.state.uproto))
    u = ComponentArray(_u, getaxes(tile.state.uproto))
    i = 1
    for j in 1:length(tile.strat)
        if layernames[j] == layername
            i = j
            break
        end
    end
    z = boundarypairs(map(ustrip, stripparams(boundaries(tile.strat))))[i]
    return LayerState(tile.state, z, u, du, t, dt, Val{layername}(), Val{iip}())
end

"""
    parameterize(tile::Tile)

Adds parameter information to all nested types in `tile` by recursively calling `parameterize`.
"""
function CryoGrid.parameterize(tile::Tile)
    ctor = ConstructionBase.constructorof(typeof(tile))
    new_layers = map(tile.strat) do named_layer
        name = layername(named_layer)
        layer = CryoGrid.parameterize(named_layer.val)
        Named(name, _addlayerfield(layer, name))
    end
    new_inits = map(tile.inits) do init
        _addlayerfield(CryoGrid.parameterize(init), :init)
    end
    new_events = map(keys(tile.events)) do name
        evs = map(CryoGrid.parameterize, getproperty(tile.events, name))
        name => _addlayerfield(evs, name)
    end
    new_strat = Stratigraphy(boundaries(tile.strat), Tuple(new_layers))
    return ctor(new_strat, tile.grid, tile.state, new_inits, (;new_events...), tile.hist)
end

"""
    variables(tile::Tile)

Returns a tuple of all variables defined in the tile.
"""
CryoGrid.variables(tile::Tile) = Tuple(unique(Flatten.flatten(tile.state.vars, Flatten.flattenable, Var)))

"""
    parameters(tile::Tile; kwargs...)

Extracts all parameters from `tile`.
"""
parameters(tile::Tile; kwargs...) = CryoGridParams(tile; kwargs...)

"""
    withaxes(u::AbstractArray, ::Tile)

Constructs a `ComponentArray` with labeled axes from the given state vector `u`. Assumes `u` to be of the same type/shape
as `setup.uproto`.
"""
withaxes(u::AbstractArray, tile::Tile) = ComponentArray(u, getaxes(tile.state.uproto))
withaxes(u::ComponentArray, ::Tile) = u

function getstate(tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip}, _u, _du, t, dt=1.0) where {TStrat,TGrid,TStates,TInits,TEvents,iip}
    du = ComponentArray(_du, getaxes(tile.state.uproto))
    u = ComponentArray(_u, getaxes(tile.state.uproto))
    return TileState(tile.state, map(ustrip ∘ stripparams, boundaries(tile.strat)), u, du, t, dt, Val{iip}())
end

"""
    resolve(tile::Tile, u, p, t)

Resolves/updates the given `tile` by:
(1) Replacing all `ModelParameters.AbstractParam` values in `tile` with their (possibly updated) value from `p`.
(2) Resolving the boundary depths of the `Stratigraphy` layers by invoking `resolveboundaries`.
(3) Replacing all instances of `DynamicParameterization` with their resolved values given the current state.

Returns the reconstructed `Tile` instance.
"""
function resolve(tile::Tile{TStrat,TGrid,TStates}, u, p, t) where {TStrat,TGrid,TStates}
    IgnoreTypes = Union{TGrid,TStates,StateHistory,Unitful.Quantity,Numerics.ForwardDiff.Dual}
    # unfortunately, reconstruct causes allocations due to a mysterious dynamic dispatch when returning the result of _reconstruct;
    # I really don't know why, could be a compiler bug, but it doesn't happen if we call the internal _reconstruct method directly...
    # so that's what we do here. The last integer argument denotes the index of the first parameter.
    tile_updated = Flatten._reconstruct(tile, p, Flatten.flattenable, ModelParameters.AbstractParam, IgnoreTypes, 1)[1]
    dynamic_ps = Flatten.flatten(tile_updated, Flatten.flattenable, DynamicParameterization, IgnoreTypes)
    # TODO: perhaps should allow dependence on local layer state;
    # this would likely require deconstruction/reconstruction of layers in order to
    # build the `LayerState`s and evaluate the dynamic parameters in a fully type stable manner.
    dynamic_values = map(d -> d(u, t), dynamic_ps)
    reconstructed_tile = Flatten._reconstruct(tile_updated, dynamic_values, Flatten.flattenable, DynamicParameterization, IgnoreTypes, 1)[1]
    return reconstructed_tile
end
resolve(tile::Tile, u, p::Nothing, t) = tile

function checkstate!(tile::Tile, state::TileState, u, du, label::Symbol)
    if CryoGrid.CRYOGRID_DEBUG
        @inbounds for i in eachindex(u)
            if !isfinite(u[i])
                debughook!(tile, state, AssertionError("[$label] Found NaN/Inf value in current state vector at index $i"))
            end
            if !isfinite(du[i])
                debughook!(tile, state, AssertionError("[$label] Found NaN/Inf value in computed time derivatives at index $i"))
            end
        end
    end
    return nothing
end

debughook!(tile, state, err) = throw(err)

# ==== Internal methods for initializing types and state variables ====
"""
Initialize `StateVars` which holds the caches for all defined state variables.
"""
function _initstatevars(@nospecialize(strat::Stratigraphy), @nospecialize(grid::Grid), @nospecialize(vars::NamedTuple), chunk_size::Union{Nothing,Int}, arrayproto::Type{A}) where {A}
    layernames = [layername(layer) for layer in strat]
    npvars = (length(filter(isprognostic, var)) + length(filter(isalgebraic, var)) for var in vars) |> sum
    ndvars = (length(filter(isdiagnostic, var)) for var in vars) |> sum
    @assert (npvars + ndvars) > 0 "No variable definitions found. Did you add a method definition for CryoGrid.variables(::L,::P) where {L<:Layer,P<:Process}?"
    @assert npvars > 0 "At least one prognostic variable must be specified."
    para = params(strat)
    chunk_size = isnothing(chunk_size) ? length(para) : chunk_size
    states = StateVars(vars, Grid(ustrip.(grid), grid.geometry), chunk_size, arrayproto)
    return states
end

# ===================================================================== #
