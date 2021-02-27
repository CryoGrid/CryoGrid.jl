# Variable dimensions
abstract type VarDim end
struct OnGrid{S} <: VarDim
    f::Function # G -> G' where G is the edge grid
    OnGrid(::Type{S}, f::Function=x->x) where {S<:GridSpec} = new{S}(f)
end
struct Shape{D} <: VarDim Shape(dims::Int...) = new{dims}() end
const Scalar = Shape()

# Variable trait (prognostic vs diagnostic)
abstract type VarStyle end
struct Prognostic <: VarStyle end
struct Diagnostic <: VarStyle end

# Variable type
struct Var{name,T,D,S}
    dim::D
    Var(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar, style::Type{S}=Diagnostic) where {T,S<:VarStyle} =
        new{name,T,typeof(dim),S}(dim)
end
varname(::Var{name}) where {name} = name
varname(::Type{<:Var{name}}) where {name} = name
vartype(::Var{name,T}) where {name,T} = T
vartype(::Type{<:Var{name,T}}) where {name,T} = T
export Var, varname
export VarDim, OnGrid, Shape, Scalar

Prognostic(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar) where T = Var(name,typ,dim,Prognostic)
Diagnostic(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar) where T = Var(name,typ,dim,Diagnostic)

VarStyle(::Type{Var{name,T,D,S}}) where {name,T,D,S} = S()
isprognostic(var::Var{name,T,D,S}) where {name,T,D,S} = S == Prognostic
isdiagnostic(var::Var{name,T,D,S}) where {name,T,D,S} = S == Diagnostic
export VarStyle, Prognostic, Diagnostic, isprognostic, isdiagnostic
