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
struct Var{T,D,S}
    dim::D
    name::Symbol
    Var(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar, style::Type{S}=Diagnostic) where {T,S<:VarStyle} =
        new{T,typeof(dim),S}(dim,name)
end
nameof(var::Var{T,D,S}) where {T,D,S} = var.name
export Var, nameof
export VarDim, OnGrid, Shape, Scalar

VarStyle(::Type{Var{T,D,S}}) where {T,D,S} = S()
isprognostic(var::Var{T,D,S}) where {T,D,S} = S == Prognostic
isdiagnostic(var::Var{T,D,S}) where {T,D,S} = S == Diagnostic
export VarStyle, isprognostic, isdiagnostic
