"""
    TileState{TGrid,TStates,Tu,Tt,Tdt}

Represents the state of a CryoGrid `Tile`.
"""
struct TileState{TStrat,TGrid,TStates,Tu,Tt,Tdt}
    strat::TStrat # stratigraphy of the Tile
    grid::TGrid # grid (from current state)
    states::TStates # named tuple of state variables
    du::Tu
    u::Tu
    t::Tt
    dt::Tdt
    function TileState(
        strat::Stratigraphy,
        grid::Numerics.AbstractDiscretization,
        states::NamedTuple,
        du::Tu,
        u::Tu,
        t,
        dt,
    ) where {Tu}
        return new{typeof(strat),typeof(grid),typeof(states),Tu,typeof(t),typeof(dt)}(strat, grid, states, du, u, t, dt)
    end
end

function TileState(
    strat::Stratigraphy{N,<:NamedTuple{layernames}},
    grid::Grid,
    sv::StateVars,
    zs::NTuple{N},
    du,
    u,
    t,
    dt=1.0,
    ::Val{sublayer}=Val{nothing}()
) where {N,layernames,sublayer}
    z_bounds = (map(tuple, zs[1:end-1], zs[2:end])..., (zs[end], ustrip(grid[end])))
    newgrid = currentgrid(sv.griddiag, grid, u, t)
    # extract state variables for each layer from cache
    states = getstatevars(Val{layernames}(), sv, grid, z_bounds, du, u, t)
    return TileState(strat, newgrid, states, du, u, t, dt)
end

Base.getindex(state::TileState, sym::Symbol) = getproperty(state, sym)
Base.propertynames(state::TileState) = (propertynames(getfield(state, :states))..., fieldnames(typeof(state))...)
function Base.getproperty(state::TState, sym::Symbol) where {TState<:TileState}
    return if sym ∈ fieldnames(TState)
        getfield(state, sym)
    else
        selectlayer(state, Val{sym}())
    end
end

@inline selectlayer(state::TileState, layername::Symbol) = selectlayer(state, Val{layername}())
@inline function selectlayer(state::TileState, ::Val{layername}) where {layername}
    # The following lines fix a performance regression on Julia 1.10;
    # we need to use getfield explicitly here to avoid a circular relationship
    # between getproperty(::TileState, sym) and this method. This seems to break
    # type inference in the compiler.
    t = getfield(state, :t)
    dt = getfield(state, :dt)
    states = getproperty(getfield(state, :states), layername)
    grid = getfield(states, :grid)
    # make subgrid
    z₁ = ForwardDiff.value(states.z[1])
    z₂ = ForwardDiff.value(max(z₁, z₁ + states.Δz[1]))
    subgrid = grid[z₁..z₂]
    # build layer state
    return LayerState(layername, subgrid, states, t, dt)
end

"""
    LayerState{TStates<:NamedTuple,TGrid<:Grid,Tt,Tdt}

State for a single layer, typically constructed from a parent `TileState`.
"""
struct LayerState{TStates,TGrid,Tt,Tdt}
    name::Symbol
    grid::TGrid
    states::TStates
    t::Tt
    dt::Tdt
end

Base.parent(state::LayerState) = state.parent
Base.getindex(state::LayerState, sym::Symbol) = getproperty(state, sym)
Base.propertynames(state::LayerState) = (propertynames(getfield(state, :states))..., fieldnames(typeof(state))...)
function Base.getproperty(state::LayerState, sym::Symbol)
    states = getfield(state, :states)
    return if sym ∈ propertynames(states)
        getproperty(states, sym)
    else
        getfield(state, sym)
    end
end

# internal method dispatches for type stable construction of state types
@inline getstatevar(::Val, ::Prognostic{name,<:OnGrid{Cells}}, sv::StateVars, z_inds, du, u, t) where {name} = view(view(u, Val{name}()), infimum(z_inds):supremum(z_inds)-1)
@inline getstatevar(::Val, ::Prognostic{name,<:OnGrid{Edges}}, sv::StateVars, z_inds, du, u, t) where {name} = error("prognostic variables on grid edges not supported")
@inline getstatevar(::Val{layername}, ::Prognostic{name,<:Shape}, sv::StateVars, z_inds, du, u, t) where {name,layername} = view(view(u, Val{layername}()), Val{name}())
@inline getstatevar(::Val, ::Algebraic{name,<:OnGrid{Cells}}, sv::StateVars, z_inds, du, u, t) where {name} = view(view(u, Val{name}()), infimum(z_inds):supremum(z_inds)-1)
@inline getstatevar(::Val, ::Algebraic{name,<:OnGrid{Edges}}, sv::StateVars, z_inds, du, u, t) where {name} = error("prognostic variables on grid edges not supported")
@inline getstatevar(::Val{layername}, ::Algebraic{name,<:Shape}, sv::StateVars, z_inds, du, u, t) where {name,layername} = view(view(u, Val{layername}()), Val{name}())
@inline getstatevar(::Val, ::DVar{dname,name,<:OnGrid{Cells}}, sv::StateVars, z_inds, du, u, t) where {dname,name} = view(view(du, Val{name}()), infimum(z_inds):supremum(z_inds)-1)
@inline getstatevar(::Val{layername}, ::DVar{dname,name,<:Shape}, sv::StateVars, z_inds, du, u, t) where {dname,name,layername} = view(view(du, Val{layername}()), Val{name}())
@inline getstatevar(::Val, var::Diagnostic{name,<:OnGrid{Cells}}, sv::StateVars, z_inds, du, u, t) where {name} = view(retrieve(sv.griddiag[name], u, t), infimum(z_inds):supremum(z_inds)-1+var.dim.offset)
@inline getstatevar(::Val, var::Diagnostic{name,<:OnGrid{Edges}}, sv::StateVars, z_inds, du, u, t) where {name} = view(retrieve(sv.griddiag[name], u, t), infimum(z_inds):supremum(z_inds)+var.dim.offset)
@inline getstatevar(::Val{layername}, ::Diagnostic{name}, sv::StateVars, z_inds, du, u, t) where {name,layername} = retrieve(sv.diag[layername][name], u, t)

# these need to be a @generated functions in order for the compiler to infer all of the types correctly
@inline @generated function getstatevars(::Val{layername}, sv::StateVars, vars::NamedTuple{varnames}, z_inds::ClosedInterval, du, u, t) where {layername,varnames}
    quote
        # QuoteNode forces the name symbols to be interpolated in as symbol literals rather than variable names!
        NamedTuple{tuple($(map(QuoteNode, varnames)...))}(tuple($(map(n -> :(getstatevar(Val{$(QuoteNode(layername))}(), sv.vars.$layername.$n, sv, z_inds, du, u, t)), varnames)...)))
    end
end
@inline @generated function getstatevars(::Val{layernames}, sv::StateVars, grid::Grid, zs::NTuple{N}, du, u, t) where {layernames,N}
    # construct expressions that call getstatevars for each layer
    exprs = map(1:N, layernames) do i, layername
        quote
            z1, z2 = zs[$i]
            z_inds = subgridinds(grid, z1..z2)
            # QuoteNode forces the name symbol to be interpolated in as symbol literals rather than variable names!
            statevars = getstatevars(
                Val{$(QuoteNode(layername))}(),
                sv,
                getproperty(sv.vars, $(QuoteNode(layername))),
                z_inds,
                du,
                u,
                t
            )
            merge((grid=grid[z_inds],), statevars)
        end
    end
    # interpolate these expressions into a NamedTuple constructor
    return :(NamedTuple{$layernames}(tuple($(exprs...))))
end
