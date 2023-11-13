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
    data::TileData # output data
    metadata::Dict # metadata
    function Tile(
        strat::TStrat,
        grid::TGrid,
        state::TStates,
        inits::TInits,
        events::TEvents,
        data::TileData=TileData(),
        metadata::Dict=Dict(),
        iip::Bool=true) where
        {TStrat<:Stratigraphy,TGrid<:Grid{Edges},TStates<:StateVars,TInits<:Tuple,TEvents<:NamedTuple}
        new{TStrat,TGrid,TStates,TInits,TEvents,iip}(strat,grid,state,inits,events,data,metadata)
    end
end
ConstructionBase.constructorof(::Type{Tile{TStrat,TGrid,TStates,TInits,TEvents,iip}}) where {TStrat,TGrid,TStates,TInits,TEvents,iip} =
    (strat, grid, state, inits, events, data, metadata) -> Tile(strat, grid, state, inits, events, data, metadata, iip)
# mark only stratigraphy and initializers fields as flattenable
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{:strat}}) = true
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{:inits}}) = true
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{:events}}) = true
Flatten.flattenable(::Type{<:Tile}, ::Type{Val{name}}) where name = false
Base.show(io::IO, ::MIME"text/plain", tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,iip}) where {TStrat,TGrid,TStates,TInits,TEvents,iip} = print(io, "Tile (iip=$iip) with layers $(keys(tile.strat)), $TGrid, $TStrat")

"""
    Tile(
        @nospecialize(strat::Stratigraphy),
        @nospecialize(discretization_strategy::DiscretizationStrategy),
        @nospecialize(inits::VarInitializer...);
        metadata::Dict=Dict(),
        arraytype::Type{A}=Vector,
        iip::Bool=true,
        chunk_size=nothing,
    )

Constructs a `Tile` from the given stratigraphy and discretization strategy. `arraytype` keyword arg should be an array instance
(of any arbitrary length, including zero, contents are ignored) that will determine the array type used for all state vectors.
"""
function Tile(
    @nospecialize(strat::Stratigraphy),
    @nospecialize(discretization_strategy::DiscretizationStrategy),
    @nospecialize(inits::VarInitializer...);
    metadata::Dict=Dict(),
    cachetype::Type{T}=DiffCache,
    arraytype::Type{A}=Vector{Float64},
    iip::Bool=true,
    chunk_size=nothing,
    strip_units=true,
) where {T<:Numerics.StateVarCache,A<:AbstractArray}
    grid = Numerics.makegrid(strat, discretization_strategy)
    if strip_units
        strat = stripunits(strat)
        grid = Grid(ustrip.(grid))
    end
    events = CryoGrid.events(strat)
    layers = map(namedlayers(strat)) do named_layer
        nameof(named_layer) => _addlayerfield(named_layer.val, nameof(named_layer))
    end
    # rebuild stratigraphy with updated parameters
    strat = Stratigraphy(boundaries(strat), (;layers...))
    # construct state variables
    states = _initstatevars(strat, grid, cachetype, arraytype; chunk_size)
    if isempty(inits)
        @warn "No initializers provided. State variables without initializers will be set to zero by default."
    end
    inits = map(inits) do init
        _addlayerfield(init, Symbol(:init_, varname(init)))
    end
    return Tile(strat, grid, states, inits, (;events...), TileData(), metadata, iip)
end
Tile(strat::Stratigraphy, grid::Grid{Cells}, inits...; kwargs...) = Tile(strat, edges(grid), inits...; kwargs...)
Tile(strat::Stratigraphy, grid::Grid{Edges}, inits...; kwargs...) = Tile(strat, PresetGrid(grid), inits...; kwargs...)
Tile(strat::Stratigraphy, inits...; discretization_strategy=AutoGrid(), kwargs...) = Tile(strat, discretization_strategy, inits...; kwargs...)
# convenience function to unwrap Tile from ODEFunction
function Tile(f::ODEFunction)
    extract_f(tile::Tile) = tile
    extract_f(f::ODEFunction) = f.f
    extract_f(f::DiffEqBase.Void) = f.f
    extract_f(f) = SciMLBase.unwrapped_f(f)
    return extract_f(f.f)
end
# and also from a DEIntegrator
"""
    Tile(integrator::SciMLBase.DEIntegrator)

Constructs a `Tile` from a `SciMLBase` DEIntegrator.
"""
function Tiles.Tile(integrator::SciMLBase.DEIntegrator)
    tile = Tiles.Tile(integrator.sol.prob.f)
    u = integrator.u
    return Tiles.resolve(tile, integrator.p, integrator.t)
end

"""
    computefluxes!(_tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,true}, _du, _u, p, t) where {TStrat,TGrid,TStates,TInits,TEvents}

Time derivative step function (i.e. du/dt) for any arbitrary `Tile`. Specialized code is generated and compiled
on the fly via the @generated macro to ensure type stability. The generated code updates each layer in the stratigraphy
in sequence, i.e for each layer 1 <= i <= N:

```julia
computediagnostic!(layer[i], ...)
interact!(layer[i], ..., layer[i+1], ...)
computefluxes!(layer[i], ...)
```
"""
function computefluxes!(
    _tile::Tile{TStrat,TGrid,TStates,TInits,TEvents,true},
    _du::AbstractVector,
    _u::AbstractVector,
    p::Union{Nothing,AbstractVector},
    t::Number,
    dt=1.0,
) where {N,TStrat<:Stratigraphy{N},TGrid,TStates,TInits,TEvents}
    _du .= zero(eltype(_du))
    du = ComponentArray(_du, getaxes(_tile.state.uproto))
    u = ComponentArray(_u, getaxes(_tile.state.uproto))
    tile = resolve(_tile, p, t)
    strat = tile.strat
    zs = map(getscalar, boundaries!(tile, u))
    state = TileState(tile.strat, tile.grid, tile.state, zs, du, u, t, dt)
    CryoGrid.resetfluxes!(strat, state)
    CryoGrid.computediagnostic!(strat, state)
    checkstate!(tile, state, u, du, :computediagnostic!)
    CryoGrid.interact!(strat, state)
    checkstate!(tile, state, u, du, :interact!)
    CryoGrid.computefluxes!(strat, state)
    checkstate!(tile, state, u, du, :computefluxes!)
    return nothing
end

CryoGrid.diagnosticstep!(tile::Tile, state::TileState) = diagnosticstep!(tile.strat, state)

"""
    timestep(_tile::Tile, _du, _u, p, t)

Computes the maximum permissible forward timestep for this `Tile` given the current `u`, `p`, and `t`.
"""
function CryoGrid.timestep(_tile::Tile, _du, _u, p, t)
    du = ComponentArray(_du, getaxes(_tile.state.uproto))
    u = ComponentArray(_u, getaxes(_tile.state.uproto))
    tile = resolve(_tile, p, t)
    strat = tile.strat
    zs = map(getscalar, boundaries!(tile, u))
    state = TileState(tile.strat, tile.grid, tile.state, zs, du, u, t, 1.0)
    dtmax = CryoGrid.timestep(strat::Stratigraphy, state)
    @assert dtmax > zero(dtmax) "timestep $dtmax cannot be <= 0"
    return dtmax
end

"""
    initialcondition!(tile::Tile, tspan::NTuple{2,Float64}[, p])
    initialcondition!(tile::Tile, tspan::NTuple{2,DateTime}[, p])

Calls `initialcondition!` on all layers/processes and returns the fully constructed u0 and du0 states.
"""
CryoGrid.initialcondition!(tile::Tile, tspan::NTuple{2}) = initialcondition!(tile, tspan, collect(map(p -> p.val, params(tile))))
CryoGrid.initialcondition!(tile::Tile, tspan::NTuple{2}, ::Nothing) = initialcondition!(tile, tspan, collect(map(p -> p.val, params(tile))))
CryoGrid.initialcondition!(tile::Tile, tspan::NTuple{2,DateTime}, p::AbstractVector) = initialcondition!(tile, convert_tspan(tspan), p)
function CryoGrid.initialcondition!(tile::Tile, tspan::NTuple{2,Float64}, p::AbstractVector)
    t0 = tspan[1]
    utype = isempty(p) ? eltype(tile.state.uproto) : eltype(p)
    du = zero(similar(tile.state.uproto, utype))
    u = zero(similar(tile.state.uproto, utype))
    tile = resolve(tile, p, t0)
    strat = tile.strat
    # get stratigraphy boundaries
    zs = initboundaries!(tile, u)
    state = TileState(tile.strat, tile.grid, tile.state, zs, du, u, t0, 1.0)
    CryoGrid.initialcondition!(tile.grid, state)
    CryoGrid.initialcondition!(strat, state, tile.inits)
    # evaluate initial time derivative
    computefluxes!(tile, du, u, p, t0)
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
            zs = update_layer_boundaries!(tile, u)
            return reverse(zs)
        end
    end
end

function initboundaries!(tile::Tile{TStrat}, u) where {TStrat}
    bounds = boundarypairs(tile.strat)
    map(bounds, namedlayers(tile.strat)) do (z1, z2), named_layer
        name = nameof(named_layer)
        diag_layer = getproperty(tile.state.diag, name)
        z = retrieve(diag_layer.z, u)
        Δz = _layerthick(tile, named_layer, u)
        @setscalar Δz = z2 - z1
        @setscalar z = z1
        return z1
    end
end

function update_layer_boundaries!(tile::Tile, u)
    # calculate grid boundaries starting from the bottom moving up to the surface
    zbot = tile.grid[end]
    return accumulate(reverse(namedlayers(tile.strat)); init=zbot) do z_acc, named_layer
        name = nameof(named_layer)
        diag_layer = getproperty(tile.state.diag, name)
        z_state = retrieve(diag_layer.z, u)
        Δz = getscalar(_layerthick(tile, named_layer, u))
        @setscalar z_state = z_acc - max(Δz, zero(Δz))
        z = getscalar(z_state)
        # strip ForwardDiff type if necessary and round to avoid numerical issues
        return round(ForwardDiff.value(z), digits=12)
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
    getstate(tile::Tile, u, du, t, dt=1.0)

Constructs a `LayerState` representing the full state of `layername` given `tile`, state vectors `u` and `du`, and the
time step `t`.
"""
function getstate(tile::Tile, _u, _du, t, dt=1.0)
    du = ComponentArray(_du, getaxes(tile.state.uproto))
    u = ComponentArray(_u, getaxes(tile.state.uproto))
    return TileState(tile.strat, tile.grid, tile.state, map(ustrip ∘ stripparams, boundaries(tile.strat)), du, u, t, dt)
end
"""
    getstate(integrator::SciMLBase.DEIntegrator)

Builds the `TileState` given an initialized integrator.
"""
getstate(integrator::SciMLBase.DEIntegrator) = Tiles.getstate(Tile(integrator), integrator.u, get_du(integrator), integrator.t)

"""
    getvar(var::Symbol, integrator::SciMLBase.DEIntegrator)
"""
Numerics.getvar(var::Symbol, integrator::SciMLBase.DEIntegrator; interp=true) = Numerics.getvar(Val{var}(), Tile(integrator), integrator.u; interp)

"""
    parameterize(tile::Tile)

Adds parameter information to all nested types in `tile` by recursively calling `parameterize`.
"""
function CryoGrid.parameterize(tile::Tile)
    ctor = ConstructionBase.constructorof(typeof(tile))
    new_layers = map(namedlayers(tile.strat)) do named_layer
        name = nameof(named_layer)
        layer = CryoGrid.parameterize(named_layer.val)
        name => _addlayerfield(layer, name)
    end
    new_inits = map(tile.inits) do init
        _addlayerfield(CryoGrid.parameterize(init), :init)
    end
    new_events = map(keys(tile.events)) do name
        evs = map(CryoGrid.parameterize, getproperty(tile.events, name))
        name => _addlayerfield(evs, name)
    end
    new_strat = Stratigraphy(boundaries(tile.strat), (;new_layers...))
    return ctor(new_strat, tile.grid, tile.state, new_inits, (;new_events...), tile.data, tile.metadata)
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

"""
    resolve(tile::Tile, p, t)

Resolves/updates the given `tile` by:
(1) Replacing all `ModelParameters.AbstractParam` values in `tile` with their (possibly updated) value from `p`.
(2) Resolving the boundary depths of the `Stratigraphy` layers by invoking `resolveboundaries`.
(3) Replacing all instances of `DynamicParameterization` with their resolved values given the current state.

Returns the reconstructed `Tile` instance.
"""
function resolve(tile::Tile{TStrat,TGrid,TStates}, p::AbstractVector, t::Number) where {TStrat,TGrid,TStates}
    IgnoreTypes = Union{TGrid,TStates,TileData,Unitful.Quantity,Numerics.ForwardDiff.Dual}
    # unfortunately, reconstruct causes allocations due to a mysterious dynamic dispatch when returning the result of _reconstruct;
    # I really don't know why, could be a compiler bug, but it doesn't happen if we call the internal _reconstruct method directly...
    # so that's what we do here. The last integer argument denotes the index of the first parameter.
    tile_updated = Flatten._reconstruct(tile, p, Flatten.flattenable, ModelParameters.AbstractParam, IgnoreTypes, 1)[1]
    dynamic_ps = Flatten.flatten(tile_updated, Flatten.flattenable, DynamicParameterization, IgnoreTypes)
    # TODO: perhaps should allow dependence on local layer state;
    # this would likely require deconstruction/reconstruction of layers in order to
    # build the `LayerState`s and evaluate the dynamic parameters in a fully type stable manner.
    dynamic_values = map(d -> d(t), dynamic_ps)
    reconstructed_tile = Flatten._reconstruct(tile_updated, dynamic_values, Flatten.flattenable, DynamicParameterization, IgnoreTypes, 1)[1]
    return reconstructed_tile
end
resolve(tile::Tile, p::Nothing, t::Number) = tile

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
function _initstatevars(@nospecialize(strat::Stratigraphy), @nospecialize(grid::Grid), cachetype::Type{T}, arraytype::Type{A}; chunk_size::Union{Nothing,Int}) where {T,A}
    stratvars = CryoGrid.variables(strat)
    gridvars = CryoGrid.variables(grid)
    npvars = (length(filter(isprognostic, var)) + length(filter(isalgebraic, var)) for var in stratvars) |> sum
    ndvars = (length(filter(isdiagnostic, var)) for var in stratvars) |> sum
    @assert (npvars + ndvars) > 0 "No variable definitions found. Did you add a method definition for CryoGrid.variables(::L,::P) where {L<:Layer,P<:Process}?"
    @assert npvars > 0 "At least one prognostic variable must be specified."
    para = params(strat)
    default_chunk_size = length(para) > 0 ? length(para) : 12
    chunk_size = isnothing(chunk_size) ? default_chunk_size : chunk_size
    vars = merge(stratvars, (grid=gridvars,))
    # create state variable cache
    states = StateVars(vars, grid, cachetype, arraytype; chunk_size)
    return states
end

# ===================================================================== #
