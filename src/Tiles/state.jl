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
    du,
    u,
    t,
    dt=1.0,
    ::Val{sublayer}=Val{nothing}()
) where {N,layernames,sublayer}
    newgrid = currentgrid(sv.griddiag, grid, u, t)
    # extract state variables for each layer from cache
    states = getstatevars(Val{tuple(layernames..., :grid)}(), sv, newgrid, du, u, t)
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
    # build layer state
    return LayerState(layername, states, t, dt)
end

"""
    LayerState{TStates<:NamedTuple,Tt,Tdt}

State for a single layer, typically constructed from a parent `TileState`.
"""
struct LayerState{TStates,Tt,Tdt}
    name::Symbol
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
@inline getstatevar(::Val, ::Prognostic{name,<:OnGrid{Cells}}, sv::StateVars, z_inds, du, u, t) where {name} = view(view(u, Val{name}()), first(z_inds):last(z_inds)-1)
@inline getstatevar(::Val, ::Prognostic{name,<:OnGrid{Edges}}, sv::StateVars, z_inds, du, u, t) where {name} = error("prognostic variables on grid edges not supported")
@inline getstatevar(::Val{layername}, ::Prognostic{name,<:Shape}, sv::StateVars, z_inds, du, u, t) where {name,layername} = view(view(u, Val{layername}()), Val{name}())
@inline getstatevar(::Val, ::Algebraic{name,<:OnGrid{Cells}}, sv::StateVars, z_inds, du, u, t) where {name} = view(view(u, Val{name}()), first(z_inds):last(z_inds)-1)
@inline getstatevar(::Val, ::Algebraic{name,<:OnGrid{Edges}}, sv::StateVars, z_inds, du, u, t) where {name} = error("prognostic variables on grid edges not supported")
@inline getstatevar(::Val{layername}, ::Algebraic{name,<:Shape}, sv::StateVars, z_inds, du, u, t) where {name,layername} = view(view(u, Val{layername}()), Val{name}())
@inline getstatevar(::Val, ::DVar{dname,name,<:OnGrid{Cells}}, sv::StateVars, z_inds, du, u, t) where {dname,name} = view(view(du, Val{name}()), first(z_inds):last(z_inds)-1)
@inline getstatevar(::Val{layername}, ::DVar{dname,name,<:Shape}, sv::StateVars, z_inds, du, u, t) where {dname,name,layername} = view(view(du, Val{layername}()), Val{name}())
@inline getstatevar(::Val, var::Diagnostic{name,<:OnGrid{Cells}}, sv::StateVars, z_inds, du, u, t) where {name} = view(retrieve(sv.griddiag[name], u, t), first(z_inds):last(z_inds)-1+var.dim.offset)
@inline getstatevar(::Val, var::Diagnostic{name,<:OnGrid{Edges}}, sv::StateVars, z_inds, du, u, t) where {name} = view(retrieve(sv.griddiag[name], u, t), first(z_inds):last(z_inds)+var.dim.offset)
@inline getstatevar(::Val{layername}, ::Diagnostic{name}, sv::StateVars, z_inds, du, u, t) where {name,layername} = retrieve(sv.diag[layername][name], u, t)

# these need to be a @generated functions in order for the compiler to infer all of the types correctly
@inline @generated function getstatevars(::Val{layername}, sv::StateVars, vars::NamedTuple{varnames}, z_inds, du, u, t) where {layername,varnames}
    quote
        # QuoteNode forces the name symbols to be interpolated in as symbol literals rather than variable names!
        NamedTuple{tuple($(map(QuoteNode, varnames)...))}(tuple($(map(n -> :(getstatevar(Val{$(QuoteNode(layername))}(), sv.vars.$layername.$n, sv, z_inds, du, u, t)), varnames)...)))
    end
end
@inline @generated function getstatevars(::Val{layernames}, sv::StateVars, grid::Grid, du, u, t) where {layernames}
    # construct expressions that call getstatevars for each layer
    exprs = map(layernames) do layername
        quote
            z_inds = sv.layeridx[$(QuoteNode(layername))]
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
            subgrid = grid[first(z_inds)..last(z_inds)]
            merge((grid=subgrid,), statevars)
        end
    end
    # interpolate these expressions into a NamedTuple constructor
    return :(NamedTuple{$layernames}(tuple($(exprs...))))
end
