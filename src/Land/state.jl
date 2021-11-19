# This ugly and confusing type alias is just to help enforce the structure of arguments used to construct state types.
# It is neither strictly necessary nor important, just there to help prevent user error :)
const GroupedVars = NamedTuple{names,<:Tuple{Vararg{<:Tuple{Vararg{<:Var}}}}} where {names}

using CryoGrid.Numerics: Flux, discretize # import flux var type

"""
Enumeration for in-place vs out-of-place mode.
"""
@enum InPlaceMode inp oop

"""
    VarStates{names,griddvars,TU,TD,TV,DF,DG}

Generic container for holding discretized state arrays for declared variables (`Var` types), as well as the prototype
prognostic state vector (`uproto`).
"""
struct VarStates{names,griddvars,TU,TD,TV,DF,DG}
    uproto::TU # prognostic state vector prototype
    grid::TD
    vars::NamedTuple{names,TV} # variable metadata
    diag::NamedTuple{names,DF} # off-grid non-prognostic variables
    griddiag::NamedTuple{griddvars,DG} # on-grid non-prognostic variables
    gridcache::Dict{ClosedInterval{Int},TD} # grid cache; indices -> subgrid
end
function VarStates(vars::GroupedVars, D::Numerics.AbstractDiscretization, chunksize::Int, arrayproto::Type{A}=Vector) where {A<:AbstractVector}
    _flatten(vars) = Flatten.flatten(vars, Flatten.flattenable, Var)
    diagvars = map(group -> filter(isdiagnostic, group), vars)
    progvars = map(group -> filter(isprognostic, group), vars)
    algvars = map(group -> filter(isalgebraic, group), vars)
    # create variables for gradients/fluxes
    dpvars = map(group -> map(Flux, filter(var -> isalgebraic(var) || isprognostic(var), group)), vars)
    allprogvars = tuplejoin(_flatten(progvars), _flatten(algvars)) # flattened prognostic/algebraic variable group
    # TODO: not a currently necessary use case, but could be supported by partitioning the state array further by layer/group name;
    # this would require passing the name to discretize(...)
    @assert let nongridpvars = filter(!isongrid, allprogvars);
        length(unique(nongridpvars)) == length(nongridpvars)
    end "Duplicated non-grid prognostic variables are not currently supported."
    # build prognostic state vector
    uproto = discretize(arrayproto, D, allprogvars...)
    # build non-gridded (i.e. "free") diagnostic state vectors
    freediagvars = map(group -> filter(!isongrid, group), diagvars)
    freediagstate = map(group -> (;map(v -> varname(v) => DiffCache(varname(v), discretize(A, D, v), chunksize), group)...), freediagvars)
    # build gridded diagnostic state vectors
    griddiagvars = filter(isongrid, _flatten(diagvars))
    griddiagstate = map(v -> varname(v) => DiffCache(varname(v), discretize(A, D, v), chunksize), griddiagvars)
    # join prognostic variables with flux variables, then build nested named tuples in each group with varnames as keys
    allvars = map(vars -> NamedTuple{map(varname, vars)}(vars), map(tuplejoin, vars, dpvars))
    VarStates(uproto, D, allvars, (;freediagstate...), (;griddiagstate...), Dict(1..length(D) => D))
end
@generated function getvar(::Val{name}, vs::VarStates{layers,griddvars,TU,TD,TV}, u::TU, du::Union{Nothing,TU}=nothing) where
    {name,layers,griddvars,T,A,pax,TU<:ComponentVector{T,A,Tuple{Axis{pax}}},TD,TV}
    dnames = map(n -> Symbol(:d,n), keys(pax))
    if name ∈ griddvars
        quote
            return retrieve(vs.griddiag.$name, u)
        end
    elseif name ∈ keys(pax)
        quote
            return u.$name
        end
    elseif du != Nothing && name ∈ dnanes
        i = findfirst(n -> n == name, dnames)::Int
        quote
            return du.$(keys(pax)[i])
        end
    else
        :(return nothing)
    end
end
function getvars(vs::VarStates{layers,gridvars}, u, du, vals::Union{Symbol,<:Pair{Symbol}}...) where {layers,gridvars,T}
    vars = map(vals) do val
        if val ∈ gridvars
            val => getvar(Val{val}(), vs, u, du)
        else
            nestedvals = isa(last(val), Tuple) ? last(val) : tuple(last(val))
            first(val) => (;map(n -> n => retrieve(getproperty(vs.diag[first(val)], n)), nestedvals)...)
        end
    end
    return (;vars...)
end

"""
    LayerState{iip,TStates,TGrids,Tt,Tz,varnames}

Represents the state of a single component (layer + processes) in the stratigraphy.
"""
struct LayerState{iip,TStates,TGrids,Tt,Tz,varnames}
    grids::NamedTuple{varnames,TGrids}
    states::NamedTuple{varnames,TStates}
    t::Tt
    z::NTuple{2,Tz}
    LayerState(grids::NamedTuple{varnames,TG}, states::NamedTuple{varnames,TS}, t::Tt, z::NTuple{2,Tz}, ::Val{iip}=Val{inp}()) where
        {TG,TS,Tt,Tz,varnames,iip} = new{iip,TS,TG,Tt,Tz,varnames}(grids, states, t, z)
end
Base.getindex(state::LayerState, sym::Symbol) = getproperty(state, sym)
function Base.getproperty(state::LayerState, sym::Symbol)
    return if sym ∈ (:strat, :grid, :grids, :t, :z)
        getfield(state, sym)
    else
        getproperty(getfield(state, :states), sym)
    end
end
@inline function LayerState(vs::VarStates, zs::NTuple{2,Tz}, u, du, t, ::Val{layername}, ::Val{iip}=Val{inp}()) where {Tz,layername,iip}
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
    LandModelState{iip,TStrat,TGrid,TStates,Tt,names}

Represents the full state of a CryoGrid `LandModel` at a single time point, including the stratigraphy and all layer states.
"""
struct LandModelState{iip,TStrat,TGrid,TStates,Tt,names}
    strat::TStrat
    grid::TGrid
    states::NamedTuple{names,TStates}
    t::Tt
    LandModelState(strat::TStrat, grid::TGrid, states::NamedTuple{names,TS}, t::Tt, ::Val{iip}=Val{inp}()) where
        {TStrat<:Stratigraphy,TGrid<:Numerics.AbstractDiscretization,TS<:Tuple{Vararg{<:LayerState}},Tt,names,iip} =
            new{iip,TStrat,TGrid,TS,Tt,names}(strat, grid, states, t)
end
function Base.getproperty(state::LandModelState, sym::Symbol)
    return if sym ∈ (:strat, :grid, :grids, :t)
        getfield(state, sym)
    else
        getproperty(getfield(state, :states), sym)
    end
end
@inline @generated function LandModelState(strat::Stratigraphy, vs::VarStates{names}, u=copy(vs.uproto), du=similar(vs.uproto), t=0.0, ::Val{iip}=Val{inp}()) where {names,iip}
    layerstates = (:(LayerState(vs, (ustrip(bounds[$i][1]), ustrip(bounds[$i][2])), u, du, t, Val{$(QuoteNode(names[i]))}(), Val{iip}())) for i in 1:length(names))
    quote
        bounds = boundaryintervals(boundaries(strat), vs.grid[end])
        return LandModelState(
            strat,
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
