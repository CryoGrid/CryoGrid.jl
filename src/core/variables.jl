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
struct Var{T,D,S,name}
    dim::D
    Var(name::Symbol, typ::T, dim::Union{<:VarDim,<:GridSpec}=Scalar, style::Type{S}=Diagnostic) where {T,S<:VarStyle} =
        new{T,typeof(dim),S,name}(dim)
end
export Var
export VarDim, OnGrid, Shape, Scalar
name(::Var{T,D,S,name}) where {T,D,S,name} = name
(var::Var{T,TDim})(grid::Grid{Edges}) where {T,TDim<:OnGrid{Edges}} = var.dim.f(grid)
(var::Var{T,TDim})(grid::Grid{Edges}) where {T,TDim<:OnGrid{Cells}} = var.dim.f(cells(grid))
(var::Var{T,TDim})(grid::Grid{Edges}) where {T,dims,TDim<:Shape{dims}} = 1:prod(dims)

VarStyle(::Type{Var{T,D,S}}) where {T,D,S} = S()
isprognostic(var::Var{T,D,S}) where {T,D,S} = S == Prognostic
isdiagnostic(var::Var{T,D,S}) where {T,D,S} = S == Diagnostic
