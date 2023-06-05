using CryoGrid: DVar

"""
    LayerState{iip,TGrid,TStates,Tt,Tdt,Tz,varnames}

Represents the state of a single component (layer + processes) in the stratigraphy.
"""
struct LayerState{iip,TGrid,TStates,TBounds,Tt,Tdt,varnames}
    grid::TGrid
    states::NamedTuple{varnames,TStates}
    bounds::TBounds
    t::Tt
    dt::Tdt
    LayerState(grid::TG, states::NamedTuple{varnames,Tvs}, bounds::TB, t::Tt, dt::Tdt=1.0, ::Val{iip}=Val{true}()) where
        {TG,TB,Tvs,Tt,Tdt,varnames,iip} = new{iip,TG,Tvs,TB,Tt,Tdt,varnames}(grid, states, bounds, t, dt)
end
Base.getindex(state::LayerState, sym::Symbol) = getproperty(state, sym)
function Base.getproperty(state::LayerState, sym::Symbol)
    return if sym ∈ (:grid, :states, :bounds, :t, :dt)
        getfield(state, sym)
    else
        getproperty(getfield(state, :states), sym)
    end
end
Base.propertynames(state::LayerState) = (propertynames(state.states)..., :grid, :states, :bounds, :t, :dt)
@inline function LayerState(sv::StateVars, zs, u, du, t, dt, ::Val{layername}, ::Val{iip}=Val{true}()) where {layername,iip}
    z_inds = subgridinds(edges(sv.grid), zs[1]..zs[2])
    return LayerState(
        sv.grid[z_inds],
        _makestates(Val{layername}(), getproperty(sv.vars, layername), sv, z_inds, u, du, t),
        zs,
        t,
        dt,
        Val{iip}(),
    )
end

"""
    TileState{iip,TGrid,TStates,Tt,Tdt,names}

Represents the instantaneous state of a CryoGrid `Tile`.
"""
struct TileState{iip,TGrid,TStates,Tt,Tdt,names}
    grid::TGrid
    states::NamedTuple{names,TStates}
    t::Tt
    dt::Tdt
    TileState(grid::TGrid, states::NamedTuple{names,TS}, t::Tt, dt::Tdt=1.0, ::Val{iip}=Val{true}()) where
        {TGrid<:Numerics.AbstractDiscretization,TS<:Tuple{Vararg{LayerState}},Tt,Tdt,names,iip} =
            new{iip,TGrid,TS,Tt,Tdt,names}(grid, states, t, dt)
end
Base.getindex(state::TileState, sym::Symbol) = Base.getproperty(state, sym)
Base.getindex(state::TileState, i::Int) = state.states[i]
function Base.getproperty(state::TileState, sym::Symbol)
    return if sym ∈ (:grid,:states,:t,:dt)
        getfield(state, sym)
    else
        getproperty(getfield(state, :states), sym)
    end
end
Base.propertynames(state::TileState) = (propertynames(state.states)...,:grid,:states,:t,:dt)
@inline @generated function TileState(sv::StateVars{names}, zs, u=copy(sv.uproto), du=similar(sv.uproto), t=0.0, dt=1.0, ::Val{iip}=Val{true}()) where {names,iip}
    layerstates = (
        quote
            bounds_i = (bounds[$i][1], bounds[$i][2])
            LayerState(sv, bounds_i, u, du, t, dt, Val{$(QuoteNode(names[i]))}(), Val{iip}())
        end
        for i in 1:length(names)
    )
    quote
        bounds = boundarypairs(map(ustrip, zs))
        return TileState(
            sv.grid,
            NamedTuple{tuple($(map(QuoteNode,names)...))}(tuple($(layerstates...))),
            t,
            dt,
            Val{iip}(),
        )
    end
end
# internal method dispatches for type stable construction of state types
@inline _makestate(::Val, ::Prognostic{name,<:OnGrid{Cells}}, sv::StateVars, z_inds, u, du, t) where {name} = view(view(u, Val{name}()), infimum(z_inds):supremum(z_inds)-1)
@inline _makestate(::Val, ::Prognostic{name,<:OnGrid{Edges}}, sv::StateVars, z_inds, u, du, t) where {name} = error("prognostic variables on grid edges not supported")
@inline _makestate(::Val{layername}, ::Prognostic{name,<:Shape}, sv::StateVars, z_inds, u, du, t) where {name,layername} = view(view(u, Val{layername}()), Val{name}())
@inline _makestate(::Val, ::Algebraic{name,<:OnGrid{Cells}}, sv::StateVars, z_inds, u, du, t) where {name} = view(view(u, Val{name}()), infimum(z_inds):supremum(z_inds)-1)
@inline _makestate(::Val, ::Algebraic{name,<:OnGrid{Edges}}, sv::StateVars, z_inds, u, du, t) where {name} = error("prognostic variables on grid edges not supported")
@inline _makestate(::Val{layername}, ::Algebraic{name,<:Shape}, sv::StateVars, z_inds, u, du, t) where {name,layername} = view(view(u, Val{layername}()), Val{name}())
@inline _makestate(::Val, ::DVar{dname,name,<:OnGrid{Cells}}, sv::StateVars, z_inds, u, du, t) where {dname,name} = view(view(du, Val{name}()), infimum(z_inds):supremum(z_inds)-1)
@inline _makestate(::Val{layername}, ::DVar{dname,name,<:Shape}, sv::StateVars, z_inds, u, du, t) where {dname,name,layername} = view(view(du, Val{layername}()), Val{name}())
@inline _makestate(::Val, var::Diagnostic{name,<:OnGrid{Cells}}, sv::StateVars, z_inds, u, du, t) where {name} = view(retrieve(sv.griddiag[name], u, t), infimum(z_inds):supremum(z_inds)-1+var.dim.offset)
@inline _makestate(::Val, var::Diagnostic{name,<:OnGrid{Edges}}, sv::StateVars, z_inds, u, du, t) where {name} = view(retrieve(sv.griddiag[name], u, t), infimum(z_inds):supremum(z_inds)+var.dim.offset)
@inline _makestate(::Val{layername}, ::Diagnostic{name}, sv::StateVars, z_inds, u, du, t) where {name,layername} = retrieve(sv.diag[layername][name], u, t)
# these need to be @generated functions in order for the compiler to infer all of the types correctly
@inline @generated function _makestates(::Val{layername}, vars::NamedTuple{varnames}, sv::StateVars, z_inds::ClosedInterval, u, du, t) where {layername,varnames}
    quote
        NamedTuple{tuple($(map(QuoteNode, varnames)...))}(tuple($(map(n -> :(_makestate(Val{$(QuoteNode(layername))}(), sv.vars.$layername.$n, sv, z_inds, u, du, t)), varnames)...)))
    end
end
