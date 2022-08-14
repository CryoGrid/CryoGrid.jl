"""
    VarStates{names,griddvars,TU,TD,TV,DF,DG}

Generic container for holding discretized state arrays for declared variables (`Var` types), as well as the prototype
prognostic state vector (`uproto`).
"""
struct VarStates{names,griddvars,TU,TD,TV,DF,DG}
    uproto::TU # prognostic state vector prototype
    grid::TD # model grid/discretization
    vars::NamedTuple{names,TV} # variable metadata
    diag::NamedTuple{names,DF} # off-grid non-prognostic variables
    griddiag::NamedTuple{griddvars,DG} # on-grid non-prognostic variables
end
@generated function getvar(::Val{name}, vs::VarStates{layers,griddvars}, u, du=nothing) where {name,layers,griddvars}
    pax = ComponentArrays.indexmap(first(ComponentArrays.getaxes(u))) # get prognostic variable index map (name -> indices)
    dnames = map(n -> Symbol(:d,n), keys(pax)) # get names of delta/derivative variables
    # case 1) variable is diagnostic and lives on the grid
    if name ∈ griddvars
        quote
            return retrieve(vs.griddiag.$name, u)
        end
    # case 2) variable is prognostic
    elseif name ∈ keys(pax)
        quote
            return u.$name
        end
    # case 3) variable is a prognostic derivative or residual
    elseif du != Nothing && name ∈ dnanes
        i = findfirst(n -> n == name, dnames)::Int
        quote
            return du.$(keys(pax)[i])
        end
    # case 4) no variables match the given name
    else
        :(return nothing)
    end
end
function getvars(vs::VarStates{layers,gridvars,TU}, u::ComponentVector, du::ComponentVector, vals::Union{Symbol,<:Pair{Symbol}}...) where {layers,gridvars,T,A,pax,TU<:ComponentVector{T,A,Tuple{Axis{pax}}}}
    # case 1: grid variable (no layer specified)
    isprognostic(name::Symbol) = name ∈ keys(pax)
    # case 2: non-grid variable on specific layer; ignore here, defer until handled below
    isprognostic(other) = false
    symbols(name::Symbol) = tuple(name)
    symbols(names::NTuple{N,Symbol}) where N = names
    # map over non-prognostic variables, selecting variables from cache
    vars = map(filter(!(isprognostic), vals)) do val # map over given variable names, ignoring prognostic variables
        # in case val is a differential var (will be nothing otherwise)
        dvar_ind = findfirst(n -> val == Symbol(:d,n), keys(pax))
        if !isnothing(dvar_ind)
            val => du[keys(pax)[dvar_ind]]    
        elseif val ∈ gridvars
            val => getvar(Val{val}(), vs, u, du)
        else
            layername = val[1]
            # handle either a single variable name or multiple, also filtering out prognostic variables
            layervars = filter(!(isprognostic), symbols(val[2]))
            layername => (;map(n -> n => retrieve(getproperty(vs.diag[layername], n)), layervars)...)
        end
    end
    return (;vars...)
end
"""
    DiffCache{N,A,Adual}

Wrapper around `PreallocationTools.DiffCache` that stores state variables in forward-diff compatible cache arrays.
"""
struct DiffCache{N,A,Adual}
    name::Symbol
    cache::Prealloc.DiffCache{A,Adual}
    function DiffCache(name::Symbol, A::AbstractArray, chunksize::Int)
        # use dual cache for automatic compatibility with ForwardDiff
        cache = Prealloc.dualcache(A, chunksize)
        new{chunksize,typeof(cache.du),typeof(cache.dual_du)}(name, cache)
    end
end
Base.show(io::IO, cache::DiffCache) = print(io, "DiffCache $(cache.name) of length $(length(cache.cache.du)) with eltype $(eltype(cache.cache.du))")
Base.show(io::IO, mime::MIME{Symbol("text/plain")}, cache::DiffCache) = show(io, cache)
function Prealloc.get_tmp(dc::Prealloc.DiffCache, ::Type{T}) where {T<:ForwardDiff.Dual}
    # this part is copied from PreallocationTools source code
    nelem = div(sizeof(T), sizeof(eltype(dc.dual_du)))*length(dc.du)
    Prealloc.ArrayInterfaceCore.restructure(dc.du, reinterpret(T, view(dc.dual_du, 1:nelem)))
end
retrieve(dc::DiffCache) = dc.cache.du
# for matching chunk sizes, retrieve from cache
retrieve(dc::DiffCache{N}, ::Type{T}) where {tag,U,N,T<:ForwardDiff.Dual{tag,U,N}} = Prealloc.get_tmp(dc.cache, T)
# otherwise, create new DiffCache on demand
retrieve(dc::DiffCache, ::Type{T}) where {tag,U,N,T<:ForwardDiff.Dual{tag,U,N}} = Prealloc.get_tmp(Prealloc.dualcache(dc.cache.du, N), T)
# for other types, try reinterpret
retrieve(dc::DiffCache, ::Type{T}) where {T} = reinterpret(T, dc.cache.du)
# overloads for accepting array types
retrieve(dc::DiffCache, u::AbstractArray{T}) where {T} = retrieve(dc, T)
retrieve(dc::DiffCache, u::AbstractArray{T}, t) where {T} = retrieve(dc, T)
# these cover cases for Rosenbrock solvers where only t has differentiable type
retrieve(dc::DiffCache, u::AbstractArray, t::T) where {T<:ForwardDiff.Dual} = retrieve(dc, T)
retrieve(dc::DiffCache, u::AbstractArray{T}, t::T) where {T<:ForwardDiff.Dual} = retrieve(dc, T)
# This ugly and confusing type alias is just to help enforce the structure of arguments used to construct state types.
# It is neither strictly necessary nor important, just there to help prevent user error :)
const GroupedVars = NamedTuple{names,<:Tuple{Vararg{<:Tuple{Vararg{<:Var}}}}} where {names}
"""
    VarStates(vars::GroupedVars, D::Numerics.AbstractDiscretization, chunksize::Int, arrayproto::Type{A}=Vector) where {A<:AbstractVector}
"""
function VarStates(@nospecialize(vars::GroupedVars), @nospecialize(D::Numerics.AbstractDiscretization), chunksize::Int, arrayproto::Type{A}=Vector) where {A<:AbstractVector}
    _flatten(vars) = Flatten.flatten(vars, Flatten.flattenable, Var)
    diagvars = map(group -> filter(isdiagnostic, group), vars)
    progvars = map(group -> filter(isprognostic, group), vars)
    algvars = map(group -> filter(isalgebraic, group), vars)
    # create variables for time delta variables (divergence/residual)
    dpvars = map(group -> map(Delta, filter(var -> isalgebraic(var) || isprognostic(var), group)), vars)
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
    griddiagvars = Tuple(unique(filter(isongrid, _flatten(diagvars))))
    griddiagstate = map(v -> varname(v) => DiffCache(varname(v), discretize(A, D, v), chunksize), griddiagvars)
    # join prognostic variables with delta and flux variables, then build nested named tuples in each group with varnames as keys
    allvars = map(vars -> NamedTuple{map(varname, vars)}(vars), map(tuplejoin, vars, dpvars))
    VarStates(uproto, D, allvars, (;freediagstate...), (;griddiagstate...))
end
