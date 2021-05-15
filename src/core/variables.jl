import Base.==

# Variable dimensions
abstract type VarDim end
struct OnGrid{S} <: VarDim
    f::Function # G -> G' where G is the edge grid
    OnGrid(::Type{S}, f::Function=x->x) where {S<:GridSpec} = new{S}(f)
end
struct Shape{D} <: VarDim Shape(dims::Int...) = new{dims}() end
const Scalar = Shape()

abstract type Var{name,T,D<:Union{<:VarDim,<:GridSpec}} end
struct Prognostic{name,T,D} <: Var{name,T,D}
    dim::D
    Prognostic(name::Symbol, ::Type{T}, dims::Union{<:VarDim,<:GridSpec}) where {T} = new{name, T, typeof(dims)}(dims)
end
struct Algebraic{name,T,D} <: Var{name,T,D}
    dim::D
    # maybe a mass matrix init function?
    Algebraic(name::Symbol, ::Type{T}, dims::Union{<:VarDim,<:GridSpec}) where {T} = new{name, T, typeof(dims)}(dims)
end
struct Diagnostic{name,T,D} <: Var{name,T,D}
    dim::D
    Diagnostic(name::Symbol, ::Type{T}, dims::Union{<:VarDim,<:GridSpec}) where {T} = new{name, T, typeof(dims)}(dims)
end
struct Parameter{name,T,D} <: Var{name,T,D}
    dim::D
    default_value::T
    Parameter(name::Symbol, default_value::T) where {T<:AbstractArray} = new{name,T,typeof(Shape(size(default_value)...))}(Shape(size(default_value)...), default_value)
    Parameter(name::Symbol, default_value::T) where {T<:Real} = Parameter(name, [default_value])
end
==(var1::Var{N1,T1,D1},var2::Var{N2,T2,D2}) where {N1,N2,T1,T2,D1,D2} = (N1==N2) && (T1==T2) && (D1==D2)
varname(::Var{name}) where {name} = name
varname(::Type{<:Var{name}}) where {name} = name
vartype(::Var{name,T}) where {name,T} = T
vartype(::Type{<:Var{name,T}}) where {name,T} = T
vardims(var::Var{name,T,D}) where {name,T,D} = var.dim
isprognostic(::T) where {T<:Var} = T <: Prognostic
isalgebraic(::T) where {T<:Var} = T <: Algebraic
isdiagnostic(::T) where {T<:Var} = T <: Diagnostic
isparameter(::T) where {T<:Var} = T <: Parameter
export Var, Prognostic, Algebraic, Diagnostic, Parameter
export VarDim, OnGrid, Shape, Scalar
export varname, vartype, isprognostic, isalgebraic, isdiagnostic, isparameter

struct VarCache{name, TCache}
    cache::TCache
    function VarCache(name::Symbol, grid::AbstractArray, arrayproto::AbstractArray, chunk_size::Int)
        # use dual cache for automatic compatibility with ForwardDiff
        cache = DiffEqBase.dualcache(similar(arrayproto, length(grid)), Val{chunk_size})
        new{name,typeof(cache)}(cache)
    end
end
# retrieve(varcache::VarCache, u::AbstractArray{T}) where {T<:ForwardDiff.Dual} = DiffEqBase.get_tmp(varcache.cache, u)
# for some reason, it's faster to re-allocate a new array of ForwardDiff.Dual than to use a pre-allocated cache...
# I have literally no idea why.
retrieve(varcache::VarCache, u::AbstractArray{T}) where {T<:ForwardDiff.Dual} = copyto!(similar(u, length(varcache.cache.du)), varcache.cache.du)
retrieve(varcache::VarCache, u::AbstractArray{T}) where {T<:ReverseDiff.TrackedReal} = copyto!(similar(u, length(varcache.cache.du)), varcache.cache.du)
retrieve(varcache::VarCache, u::ReverseDiff.TrackedArray) = copyto!(similar(identity.(u), length(varcache.cache.du)), varcache.cache.du)
retrieve(varcache::VarCache, u::AbstractArray{T}) where {T} = reinterpret(T, varcache.cache.du)
retrieve(varcache::VarCache) = diffcache.du
Base.show(io::IO, cache::VarCache{name}) where name = print(io, "VarCache{$name} of length $(length(cache.cache.du)) with eltype $(eltype(cache.cache.du))")
Base.show(io::IO, mime::MIME{Symbol("text/plain")}, cache::VarCache{name}) where name = show(io, cache)
# type piracy to reduce clutter in compiled type names
Base.show(io::IO, ::Type{<:VarCache{name}}) where name = print(io, "VarCache{$name}")
