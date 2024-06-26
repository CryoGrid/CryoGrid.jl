"""
Base type for grid alignment, i.e. either "edges" or "cells".
"""
abstract type GridOffset end
struct Edges <: GridOffset end
struct Cells <: GridOffset end

abstract type VarDim{S} end
struct OnGrid{S} <: VarDim{S}
    offset::Int # G -> G' where G is the edge grid
    OnGrid(::Type{S}, offset::Int=0) where {S<:GridOffset} = new{S}(offset)
end

struct Shape{S} <: VarDim{S} Shape(dims::Int...) = new{dims}() end

const Scalar = Shape()

ConstructionBase.constructorof(::Type{OnGrid{S}}) where {S} = f -> OnGrid(S, f)
ConstructionBase.constructorof(::Type{Shape{S}}) where {S} = f -> Shpae(S...)
dimlength(::Shape{dims}, N::Int) where dims = prod(dims)
dimlength(d::OnGrid{Cells}, gridlen::Int) = gridlen + d.offset - 1
dimlength(d::OnGrid{Edges}, gridlen::Int) = gridlen + d.offset

"""
    Var{name,S<:VarDim,T,units,domain}

Base type for symbolic state variables in the model.
"""
abstract type Var{name,S<:VarDim,T,units,domain} end

"""
    Prognostic{name,S,T,units,domain} <: Var{name,S,T,units,domain}

Defines a prognostic (time-integrated) state variable.
"""
struct Prognostic{name,S,T,units,domain} <: Var{name,S,T,units,domain}
    dim::S
    desc::String
    Prognostic(name::Symbol, dims::Union{<:Shape,OnGrid{Cells}}, units=NoUnits, ::Type{T}=Float64; domain=-Inf..Inf, desc="") where {T} = new{name,typeof(dims),T,units,domain}(dims, desc)
    Prognostic(::Symbol, dims::OnGrid, args...) = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
end

"""
    Algebraic{name,S,T,units,domain} <: Var{name,S,T,units,domain}

Defines an algebraic (implicit) state variable.
"""
struct Algebraic{name,S,T,units,domain} <: Var{name,S,T,units,domain}
    dim::S
    desc::String
    # maybe a mass matrix init function?
    Algebraic(name::Symbol, dims::Union{<:Shape,OnGrid{Cells}}, units=NoUnits, ::Type{T}=Float64; domain=-Inf..Inf, desc="") where {T} = new{name,typeof(dims),T,units,domain}(dims, desc)
    Algebraic(::Symbol, dims::OnGrid, args...) = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
end

"""
    DVar{dname,name,S,T,units,domain} <: Var{dname,S,T,units,domain}

Defines a delta/differential term `du` for variable `u`, which is the time-derivative/divergence for prognostic variables or
the residual for algebraic variables.
"""
struct DVar{dname,name,S,T,units,domain} <: Var{dname,S,T,units,domain}
    dim::S
    desc::String
    DVar(dname::Symbol, name::Symbol, dims::Union{<:Shape,OnGrid{Cells}}, units=NoUnits, ::Type{T}=Float64; domain=-Inf..Inf, desc="") where {T} = new{dname,name,typeof(dims),T,units,domain}(dims, desc)
    DVar(var::Prognostic{name,S,T,units,domain}) where {name,S,T,units,domain} = let dims=vardims(var); new{dname(name),name,typeof(dims),T,upreferred(units)/u"s",domain}(dims) end
    DVar(var::Algebraic{name,S,T,units,domain}) where {name,S,T,units,domain} = let dims=vardims(var); new{dname(name),name,typeof(dims),T,units,domain}(dims) end
end
dname(sym::Symbol) = Symbol(:d,sym)

"""
    Diagnostic{name,S,T,units,domain} <: Var{name,S,T,units,domain}

Defines a diagnostic variable which is allocated and cached per timestep but not integrated over time.
"""
struct Diagnostic{name,S,T,units,domain} <: Var{name,S,T,units,domain}
    dim::S
    desc::String
    Diagnostic(name::Symbol, dims::VarDim, units=NoUnits, ::Type{T}=Float64; domain=-Inf..Inf, desc="") where {T} = new{name,typeof(dims),T,units,domain}(dims, desc)
end
# constructors for Flatten.reconstruct
ConstructionBase.constructorof(::Type{Prognostic{name,S,T,units,domain}}) where {name,S,T,units,domain} = s -> Prognostic(name, s, units, T; domain)
ConstructionBase.constructorof(::Type{Diagnostic{name,S,T,units,domain}}) where {name,S,T,units,domain} = s -> Diagnostic(name, s, units, T; domain)
ConstructionBase.constructorof(::Type{Algebraic{name,S,T,units,domain}}) where {name,S,T,units,domain} = s -> Algebraic(name, s, units, T; domain)
ConstructionBase.constructorof(::Type{DVar{dname,name,S,T,units,domain}}) where {dname,name,S,T,units,domain} = s -> DVar(dname, name, s, units, T; domain)

import Base.:(==)
# Base overrides
==(var1::Var{N1,S1,T1,u1,d1}, var2::Var{N2,S2,T2,u2,d2}) where {N1,N2,S1,S2,T1,T2,u1,u2,d1,d2} = (N1 == N2) && (S1 == S2) && (T1 == T2) && (u1 == u2) && (d1 == d2) && vardims(var1) == vardims(var2) && vardesc(var1) == vardesc(var2)
Base.show(io::IO, ::MIME"text/plain", var::TVar) where {name,S,T,TVar<:Var{name,S,T}} = show("$(typeof(var).name.wrapper){:$name}(shape=$S, type=$(T), desc=$(var.desc))")
# other methods
varname(::Var{name}) where {name} = name
varname(::Type{<:Var{name}}) where {name} = name
basevarname(::DVar{dname,name}) where {dname,name} = name
basevarname(::Type{<:DVar{dname,name}}) where {dname,name} = name
vartype(::Var{name,S,T}) where {name,S,T} = T
vartype(::Type{<:Var{name,S,T}}) where {name,S,T} = T
varunits(::Var{name,S,T,units}) where {name,S,T,units} = units
varunits(::Type{<:Var{name,S,T,units}}) where {name,S,T,units} = units
vardomain(::Var{name,S,T,units,domain}) where {name,S,T,units,domain} = domain
vardomain(::Type{<:Var{name,S,T,units,domain}}) where {name,S,T,units,domain} = domain
vardims(var::Var) = var.dim
vardesc(var::Var) = var.desc
isprognostic(::T) where {T<:Var} = T <: Prognostic
isalgebraic(::T) where {T<:Var} = T <: Algebraic
isflux(::T) where {T<:Var} = T <: DVar
isdiagnostic(::T) where {T<:Var} = T <: Diagnostic
isongrid(::Var{name,S}) where {name,S} = S <: OnGrid
