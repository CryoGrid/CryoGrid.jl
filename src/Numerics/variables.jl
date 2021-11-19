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

abstract type Var{name,T,S<:VarDim} end
struct Prognostic{name,T,S} <: Var{name,T,S}
    dim::S
    Prognostic(name::Symbol, ::Type{T}, dims::Union{<:Shape,OnGrid{Cells,typeof(identity)}}) where {T} = new{name, T, typeof(dims)}(dims)
    Prognostic(::Symbol, ::Type{T}, dims::OnGrid) where {T} = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
end
struct Algebraic{name,T,S} <: Var{name,T,S}
    dim::S
    # maybe a mass matrix init function?
    Algebraic(name::Symbol, ::Type{T}, dims::Union{<:Shape,OnGrid{Cells,typeof(identity)}}) where {T} = new{name, T, typeof(dims)}(dims)
    Algebraic(::Symbol, ::Type{T}, dims::OnGrid) where {T} = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
end
struct Flux{dname,name,T,S} <: Var{dname,T,S}
    dim::S
    Flux(::Symbol, ::Type{T}, dims::OnGrid) where {T} = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
    Flux(dname::Symbol, name::Symbol, ::Type{T}, dims::Union{<:Shape,OnGrid{Cells,typeof(identity)}}) where {T} = new{dname, name, T, typeof(dims)}(dims)
    Flux(var::Prognostic{name,T,S}) where {name,T,S} = let dims=vardims(var); new{Symbol(:d,name), name, T, typeof(dims)}(dims) end
    Flux(var::Algebraic{name,T,S}) where {name,T,S} = let dims=vardims(var); new{Symbol(:d,name), name, T, typeof(dims)}(dims) end
end
struct Diagnostic{name,T,S} <: Var{name,T,S}
    dim::S
    Diagnostic(name::Symbol, ::Type{T}, dims::VarDim) where {T} = new{name, T, typeof(dims)}(dims)
end
ConstructionBase.constructorof(::Type{Prognostic{name,T,S}}) where {name,T,S} = s -> Prognostic(name, T, s)
ConstructionBase.constructorof(::Type{Diagnostic{name,T,S}}) where {name,T,S} = s -> Diagnostic(name, T, s)
ConstructionBase.constructorof(::Type{Algebraic{name,T,S}}) where {name,T,S} = s -> Algebraic(name, T, s)
ConstructionBase.constructorof(::Type{Flux{dname,name,T,S}}) where {dname,name,T,S} = s -> Flux(dname, name, T, s)
==(var1::Var{N1,T1,D1},var2::Var{N2,T2,D2}) where {N1,N2,T1,T2,D1,D2} = (N1==N2) && (T1==T2) && (D1==D2)
varname(::Var{name}) where {name} = name
varname(::Type{<:Var{name}}) where {name} = name
vartype(::Var{name,T}) where {name,T} = T
vartype(::Type{<:Var{name,T}}) where {name,T} = T
vardims(var::Var{name,T,S}) where {name,T,S} = var.dim
isprognostic(::T) where {T<:Var} = T <: Prognostic
isalgebraic(::T) where {T<:Var} = T <: Algebraic
isflux(::T) where {T<:Var} = T <: Flux
isdiagnostic(::T) where {T<:Var} = T <: Diagnostic
isongrid(::Var{name,T,S}) where {name,T,S} = S <: OnGrid
