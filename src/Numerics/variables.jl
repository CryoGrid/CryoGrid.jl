# Variable dimensions
abstract type VarDim{S} end
struct OnGrid{S,F} <: VarDim{S}
    f::F # G -> G' where G is the edge grid
    OnGrid(::Type{S}, f::Function=identity) where {S<:GridSpec} = new{S,typeof(f)}(f)
end
struct Shape{S} <: VarDim{S} Shape(dims::Int...) = new{dims}() end
const Scalar = Shape()
ConstructionBase.constructorof(::Type{OnGrid{S}}) where {S} = f -> OnGrid(S, f)
ConstructionBase.constructorof(::Type{Shape{S}}) where {S} = f -> Shpae(S...)
dimlength(::Shape{dims}, grid::Grid) where dims = prod(dims)
dimlength(d::OnGrid{Cells}, grid::Grid) = d.f(length(cells(grid)))
dimlength(d::OnGrid{Edges}, grid::Grid) = d.f(length(edges(grid)))

abstract type Var{name,S<:VarDim,T,units,domain} end
"""
    Prognostic{name,S,T,units,domain} <: Var{name,S,T,units,domain}

Defines a prognostic (time-integrated) state variable.
"""
struct Prognostic{name,S,T,units,domain} <: Var{name,S,T,units,domain}
    dim::S
    Prognostic(name::Symbol, dims::Union{<:Shape,OnGrid{Cells,typeof(identity)}}, units=NoUnits, ::Type{T}=Float64; domain=-Inf..Inf) where {T} = new{name,typeof(dims),T,units,domain}(dims)
    Prognostic(::Symbol, dims::OnGrid, args...) = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
end
"""
    Algebraic{name,S,T,units,domain} <: Var{name,S,T,units,domain}

Defines an algebraic (implicit) state variable.
"""
struct Algebraic{name,S,T,units,domain} <: Var{name,S,T,units,domain}
    dim::S
    # maybe a mass matrix init function?
    Algebraic(name::Symbol, dims::Union{<:Shape,OnGrid{Cells,typeof(identity)}}, units=NoUnits, ::Type{T}=Float64; domain=-Inf..Inf) where {T} = new{name,typeof(dims),T,units,domain}(dims)
    Algebraic(::Symbol, dims::OnGrid, args...) = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
end
"""
    Delta{dname,name,S,T,units,domain} <: Var{dname,S,T,units,domain}

Defines a "delta" term `du` for variable `u`, which is the time-derivative/divergence for prognostic variables or
the residual for algebraic variables.
"""
struct Delta{dname,name,S,T,units,domain} <: Var{dname,S,T,units,domain}
    dim::S
    Delta(dname::Symbol, name::Symbol, dims::Union{<:Shape,OnGrid{Cells,typeof(identity)}}, units=NoUnits, ::Type{T}=Float64; domain=-Inf..Inf) where {T} = new{dname,name,typeof(dims),T,units,domain}(dims)
    Delta(var::Prognostic{name,S,T,units,domain}) where {name,S,T,units,domain} = let dims=vardims(var); new{deltaname(name),name,typeof(dims),T,upreferred(units)/u"s",domain}(dims) end
    Delta(var::Algebraic{name,S,T,units,domain}) where {name,S,T,units,domain} = let dims=vardims(var); new{deltaname(name),name,typeof(dims),T,units,domain}(dims) end
end
deltaname(sym::Symbol) = Symbol(:∂,sym,:∂t)
"""
    Diagnostic{name,S,T,units,domain} <: Var{name,S,T,units,domain}

Defines a diagnostic variable which is allocated and cached per timestep but not integrated over time.
"""
struct Diagnostic{name,S,T,units,domain} <: Var{name,S,T,units,domain}
    dim::S
    Diagnostic(name::Symbol, dims::VarDim, units=NoUnits, ::Type{T}=Float64; domain=-Inf..Inf) where {T} = new{name,typeof(dims),T,units,domain}(dims)
end
# constructors for Flatten.reconstruct
ConstructionBase.constructorof(::Type{Prognostic{name,S,T,units,domain}}) where {name,S,T,units,domain} = s -> Prognostic(name, s, units, T; domain)
ConstructionBase.constructorof(::Type{Diagnostic{name,S,T,units,domain}}) where {name,S,T,units,domain} = s -> Diagnostic(name, s, units, T; domain)
ConstructionBase.constructorof(::Type{Algebraic{name,S,T,units,domain}}) where {name,S,T,units,domain} = s -> Algebraic(name, s, units, T; domain)
ConstructionBase.constructorof(::Type{Delta{dname,name,S,T,units,domain}}) where {dname,name,S,T,units,domain} = s -> Delta(dname, name, s, units, T; domain)
==(::Var{N1,S1,T1,u1,d1},::Var{N2,S2,T2,u2,d2}) where {N1,N2,S1,S2,T1,T2,u1,u2,d1,d2} = (N1 == N2) && (S1 == S2) && (T1 == T2) && (u1 == u2) && (d1 == d2)
varname(::Var{name}) where {name} = name
varname(::Type{<:Var{name}}) where {name} = name
basevarname(::Delta{dname,name}) where {dname,name} = name
basevarname(::Type{<:Delta{dname,name}}) where {dname,name} = name
vartype(::Var{name,S,T}) where {name,S,T} = T
vartype(::Type{<:Var{name,S,T}}) where {name,S,T} = T
varunits(::Var{name,S,T,units}) where {name,S,T,units} = units
varunits(::Type{<:Var{name,S,T,units}}) where {name,S,T,units} = units
vardomain(::Var{name,S,T,units,domain}) where {name,S,T,units,domain} = domain
vardomain(::Type{<:Var{name,S,T,units,domain}}) where {name,S,T,units,domain} = domain
vardims(var::Var) = var.dim
isprognostic(::T) where {T<:Var} = T <: Prognostic
isalgebraic(::T) where {T<:Var} = T <: Algebraic
isflux(::T) where {T<:Var} = T <: Delta
isdiagnostic(::T) where {T<:Var} = T <: Diagnostic
isongrid(::Var{name,S}) where {name,S} = S <: OnGrid
