"""
    StateVars{names,griddvars,TU,TD,TV,DF,DG}

Generic container for holding discretized state arrays for declared variables (`Var` types), as well as the prototype
prognostic state vector (`uproto`).
"""
struct StateVars{names,griddvars,TU,TD,TV,DF,DG}
    uproto::TU # prognostic state vector prototype
    grid::TD # model grid/discretization
    vars::NamedTuple{names,TV} # variable metadata
    diag::NamedTuple{names,DF} # non-grid non-prognostic variables
    griddiag::NamedTuple{griddvars,DG} # grid non-prognostic variables
end
# This ugly and confusing type alias is just to help enforce the structure of arguments used to construct state types.
# It is neither strictly necessary nor important, just there to help prevent user error :)
const GroupedVars = NamedTuple{names,<:Tuple{Vararg{<:Tuple{Vararg{<:Var}}}}} where {names}
"""
    StateVars(vars::GroupedVars, D::Numerics.AbstractDiscretization, chunk_size::Int, arrayproto::Type{A}=Vector) where {A<:AbstractVector}
"""
function StateVars(@nospecialize(vars::GroupedVars), @nospecialize(D::Numerics.AbstractDiscretization), chunk_size::Int, ::Type{A}=Vector) where {A<:AbstractVector}
    _flatten(vars) = Flatten.flatten(vars, Flatten.flattenable, Var)
    diagvars = map(group -> filter(isdiagnostic, group), vars)
    progvars = map(group -> filter(isprognostic, group), vars)
    algvars = map(group -> filter(isalgebraic, group), vars)
    # create variables for time delta variables (divergence/residual)
    dpvars = map(group -> map(Delta, filter(var -> isalgebraic(var) || isprognostic(var), group)), vars)
    gridprogvars = Tuple(unique(filter(isongrid, tuplejoin(_flatten(progvars), _flatten(algvars)))))
    freeprogvars = map(group -> filter(!isongrid, group), progvars)
    vartypes = map(vartype, tuplejoin(gridprogvars, _flatten(freeprogvars)))
    @assert all(map(==(first(vartypes)), vartypes)) "All prognostic variables must have same data type"
    uproto = prognosticstate(A{first(vartypes)}, D, freeprogvars, gridprogvars)
    # build non-gridded (i.e. "free") diagnostic state vectors
    freediagvars = map(group -> filter(!isongrid, group), diagvars)
    freediagstate = map(group -> (;map(v -> varname(v) => DiffCache(varname(v), discretize(A, D, v), chunk_size), group)...), freediagvars)
    # build gridded diagnostic state vectors
    griddiagvars = Tuple(unique(filter(isongrid, _flatten(diagvars))))
    griddiagstate = map(v -> varname(v) => DiffCache(varname(v), discretize(A, D, v), chunk_size), griddiagvars)
    # join prognostic variables with delta and flux variables, then build nested named tuples in each group with varnames as keys
    allvars = map(vars -> NamedTuple{map(varname, vars)}(vars), map(tuplejoin, vars, dpvars))
    return StateVars(uproto, D, allvars, (;freediagstate...), (;griddiagstate...))
end
@generated function getvar(::Val{name}, vs::StateVars{layers,griddvars}, u, du=nothing) where {name,layers,griddvars}
    pax = ComponentArrays.indexmap(first(ComponentArrays.getaxes(u))) # get prognostic variable index map (name -> indices)
    dnames = map(n -> deltaname(n), keys(pax)) # get names of delta/derivative variables
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
function getvars(vs::StateVars{layers,gridvars,TU}, u::ComponentVector, du::ComponentVector, vals::Union{Symbol,<:Pair{Symbol}}...) where {layers,gridvars,T,A,pax,TU<:ComponentVector{T,A,Tuple{Axis{pax}}}}
    # case 1: grid variable (no layer specified)
    isprognostic(name::Symbol) = name ∈ keys(pax)
    # case 2: non-grid variable on specific layer; ignore here, defer until handled below
    isprognostic(other) = false
    symbols(name::Symbol) = tuple(name)
    symbols(names::NTuple{N,Symbol}) where N = names
    # map over non-prognostic variables, selecting variables from cache
    vars = map(filter(!isprognostic, vals)) do val # map over given variable names, ignoring prognostic variables
        # in case val is a differential var (will be nothing otherwise)
        dvar_ind = findfirst(n -> val == deltaname(n), keys(pax))
        if !isnothing(dvar_ind)
            val => du[keys(pax)[dvar_ind]]    
        elseif val ∈ gridvars
            val => getvar(Val{val}(), vs, u, du)
        elseif isa(val, Pair)
            layername = val[1]
            # handle either a single variable name or multiple, also filtering out prognostic variables
            layervars = filter(!isprognostic, symbols(val[2]))
            layername => (;map(n -> n => retrieve(getproperty(vs.diag[layername], n)), layervars)...)
        else
            error("no state variable named $val defined")
        end
    end
    return (;vars...)
end
"""
    build_mass_matrix(u::ComponentVector, states::StateVars)

Constructs a mass matrix `M⋅∂u∂t = f(u)` suitable for the prognostic state vector `u` based on the
defined variable types.
"""
function build_mass_matrix(states::StateVars)
    M_diag = similar(states.uproto)
    M_idxmap = ComponentArrays.indexmap(getaxes(M_diag)[1])
    allvars = Flatten.flatten(states.vars, Flatten.flattenable, Var)
    progvars = map(varname, filter(isprognostic, allvars))
    algvars = map(varname, filter(isalgebraic, allvars))
    for name in keys(M_idxmap)
        M_diag_var = @view M_diag[name]
        if isa(M_diag_var, ComponentArray)
            for layervar in keys(M_diag_var)
                M_diag_sub_var = @view M_diag_var[layervar]
                if layervar ∈ progvars
                    M_diag_sub_var .= one(eltype(M_diag))
                elseif layervar ∈ algvars
                    M_diag_sub_var .= zero(eltype(M_diag))
                end
            end
        else
            if name ∈ progvars
                M_diag_var .= one(eltype(M_diag))
            elseif name ∈ algvars
                M_diag_var .= zero(eltype(M_diag))
            end
        end
    end
    # if no algebraic variables are present, use identity matrix
    num_algebraic = length(M_diag) - sum(M_diag)
    M = num_algebraic > 0 ? Diagonal(M_diag) : I
end
