using CryoGrid.Numerics: Flux
"""
Enumeration for in-place vs out-of-place mode.
"""
@enum InPlaceMode inplace ooplace
"""
    LayerState{iip,TStates,TGrids,Tt,Tz,varnames}

Represents the state of a single component (layer + processes) in the stratigraphy.
"""
struct LayerState{iip,TStates,TGrids,Tt,Tz,varnames}
    grids::NamedTuple{varnames,TGrids}
    states::NamedTuple{varnames,TStates}
    bounds::NTuple{2,Tz}
    t::Tt
    LayerState(grids::NamedTuple{varnames,TG}, states::NamedTuple{varnames,TS}, t::Tt, bounds::NTuple{2,Tz}, ::Val{iip}=Val{inplace}()) where
        {TG,TS,Tt,Tz,varnames,iip} = new{iip,TS,TG,Tt,Tz,varnames}(grids, states, bounds, t)
end
Base.getindex(state::LayerState, sym::Symbol) = getproperty(state, sym)
function Base.getproperty(state::LayerState, sym::Symbol)
    return if sym ∈ (:grids, :states, :t, :z)
        getfield(state, sym)
    else
        getproperty(getfield(state, :states), sym)
    end
end
@inline function LayerState(vs::VarStates, zs::NTuple{2,Tz}, u, du, t, ::Val{layername}, ::Val{iip}=Val{inplace}()) where {Tz,layername,iip}
    z_inds = subgridinds(edges(vs.grid), zs[1]..zs[2])
    return LayerState(
        _makegrids(Val{layername}(), getproperty(vs.vars, layername), vs, z_inds),
        _makestates(Val{layername}(), getproperty(vs.vars, layername), vs, z_inds, u, du, t),
        t,
        zs,
        Val{iip}(),
    )
end

"""
    TileState{iip,TGrid,TStates,Tt,names}

Represents the instantaneous state of a CryoGrid `Tile`.
"""
struct TileState{iip,TGrid,TStates,Tt,names}
    grid::TGrid
    states::NamedTuple{names,TStates}
    t::Tt
    TileState(grid::TGrid, states::NamedTuple{names,TS}, t::Tt, ::Val{iip}=Val{inplace}()) where
        {TGrid<:Numerics.AbstractDiscretization,TS<:Tuple{Vararg{<:LayerState}},Tt,names,iip} =
            new{iip,TGrid,TS,Tt,names}(grid, states, t)
end
Base.getindex(state::TileState, sym::Symbol) = Base.getproperty(state, sym)
Base.getindex(state::TileState, i::Int) = state.states[i]
function Base.getproperty(state::TileState, sym::Symbol)
    return if sym ∈ (:grid,:states,:t)
        getfield(state, sym)
    else
        getproperty(getfield(state, :states), sym)
    end
end
@inline @generated function TileState(vs::VarStates{names}, zs::NTuple, u=copy(vs.uproto), du=similar(vs.uproto), t=0.0, ::Val{iip}=Val{inplace}()) where {names,iip}
    layerstates = (:(LayerState(vs, (ustrip(bounds[$i][1]), ustrip(bounds[$i][2])), u, du, t, Val{$(QuoteNode(names[i]))}(), Val{iip}())) for i in 1:length(names))
    quote
        bounds = boundarypairs(zs, vs.grid[end])
        return TileState(
            vs.grid,
            NamedTuple{tuple($(map(QuoteNode,names)...))}(tuple($(layerstates...))),
            t,
            Val{iip}(),
        )
    end
end
# internal method dispatches for type stable construction of state types
@inline function _makegrid(::Var{name,T,<:OnGrid{Cells}}, vs::VarStates, z_inds) where {name,T}
    let inds=infimum(z_inds)..supremum(z_inds)-1; # subtract one due to account for cell ofset
        get!(vs.gridcache, inds) do
            vs.grid[inds]
        end
    end
end
@inline _makegrid(::Var{name,T,<:OnGrid{Edges}}, vs::VarStates, z_inds) where {name,T} = get!(() -> vs.grid[z_inds], vs.gridcache, z_inds)
@inline _makegrid(var::Var, vs::VarStates, z_inds) = 1:dimlength(vardims(var), vs.grid)
@inline _makestate(::Val, ::Prognostic{name,T,<:OnGrid{Cells}}, vs::VarStates, z_inds, u, du, t) where {name,T} = view(view(u, Val{name}()), infimum(z_inds):supremum(z_inds)-1)
@inline _makestate(::Val, ::Prognostic{name,T}, vs::VarStates, z_inds, u, du, t) where {name,T} = view(u, Val{name}())
@inline _makestate(::Val, ::Algebraic{name,T,<:OnGrid{Cells}}, vs::VarStates, z_inds, u, du, t) where {name,T} = view(view(u, Val{name}()), infimum(z_inds):supremum(z_inds)-1)
@inline _makestate(::Val, ::Algebraic{name,T}, vs::VarStates, z_inds, u, du, t) where {name,T} = view(u, Val{name}())
@inline _makestate(::Val, ::Flux{dname,name,T,<:OnGrid{Cells}}, vs::VarStates, z_inds, u, du, t) where {dname,name,T} = view(view(du, Val{name}()), infimum(z_inds):supremum(z_inds)-1)
@inline _makestate(::Val, ::Flux{dname,name,T}, vs::VarStates, z_inds, u, du, t) where {dname,name,T} = view(du, Val{name}())
@inline _makestate(::Val, ::Diagnostic{name,T,<:OnGrid{Cells}}, vs::VarStates, z_inds, u, du, t) where {name,T} = view(retrieve(vs.griddiag[name], u, t), infimum(z_inds):supremum(z_inds)-1)
@inline _makestate(::Val, ::Diagnostic{name,T,<:OnGrid{Edges}}, vs::VarStates, z_inds, u, du, t) where {name,T} = view(retrieve(vs.griddiag[name], u, t), infimum(z_inds):supremum(z_inds))
@inline _makestate(::Val{layername}, ::Diagnostic{name,T}, vs::VarStates, z_inds, u, du, t) where {name,layername,T} = retrieve(vs.diag[layername][name], u, t)
# these need to be @generated functions in order for the compiler to infer all of the types correctly
@inline @generated function _makegrids(::Val{layername}, vars::NamedTuple{varnames}, vs::VarStates, z_inds::ClosedInterval) where {layername,varnames}
    quote
        NamedTuple{tuple($(map(QuoteNode, varnames)...))}(tuple($(map(n -> :(_makegrid(vs.vars.$layername.$n, vs, z_inds)), varnames)...)))
    end
end
@inline @generated function _makestates(::Val{layername}, vars::NamedTuple{varnames}, vs::VarStates, z_inds::ClosedInterval, u, du, t) where {layername,varnames}
    quote
        NamedTuple{tuple($(map(QuoteNode, varnames)...))}(tuple($(map(n -> :(_makestate(Val{$(QuoteNode(layername))}(), vs.vars.$layername.$n, vs, z_inds, u, du, t)), varnames)...)))
    end
end
