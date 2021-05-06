import Base.==

# Variable dimensions
abstract type VarDim end
struct OnGrid{S} <: VarDim
    f::Function # G -> G' where G is the edge grid
    OnGrid(::Type{S}, f::Function=x->x) where {S<:GridSpec} = new{S}(f)
end
struct Shape{D} <: VarDim Shape(dims::Int...) = new{dims}() end
const Scalar = Shape()

# Variable "trait" (prognostic, diagnostic, parameter)
# Not really a trait since it's only used with one type, Var.
abstract type VarStyle end
struct Prognostic <: VarStyle end
struct Diagnostic <: VarStyle end
struct Parameter <: VarStyle end

# Variable type
struct Var{name,T,D,S}
    dim::D
    Var(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar, style::Type{S}=Diagnostic) where {T,S<:VarStyle} =
        new{name,T,typeof(dim),S}(dim)
end
==(var1::Var{N1,T1,D1},var2::Var{N2,T2,D2}) where {N1,N2,T1,T2,D1,D2} = (N1==N2) && (T1==T2) && (D1==D2)
varname(::Var{name}) where {name} = name
varname(::Type{<:Var{name}}) where {name} = name
vartype(::Var{name,T}) where {name,T} = T
vartype(::Type{<:Var{name,T}}) where {name,T} = T
export Var, varname
export VarDim, OnGrid, Shape, Scalar

Prognostic(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar) where T = Var(name,typ,dim,Prognostic)
Diagnostic(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar) where T = Var(name,typ,dim,Diagnostic)
Parameter(name::Symbol, typ::Type{T}=Float64, dim::Union{<:VarDim,<:GridSpec}=Scalar) where T = Var(name,typ,dim,Parameter)

VarStyle(::Type{Var{name,T,D,S}}) where {name,T,D,S} = S()
isprognostic(var::Var{name,T,D,S}) where {name,T,D,S} = S == Prognostic
isdiagnostic(var::Var{name,T,D,S}) where {name,T,D,S} = S == Diagnostic
isparameter(var::Var{name,T,D,S}) where {name,T,D,S} = S == Parameter
export VarStyle, Prognostic, Diagnostic, isprognostic, isdiagnostic, isparameter

struct VarCache{TCaches}
    caches::TCaches
    function VarCache(vargrids::NamedTuple, arrayproto::AbstractArray)
        # use dual cache for automatic compatibility with ForwardDiff
        caches = map(grid -> DiffEqBase.dualcache(similar(arrayproto, length(grid))), vargrids)
        new{typeof(caches)}(caches)
    end
end

function Base.getproperty(cache::VarCache, name::Symbol)
    if name == :caches
        getfield(cache, name)
    else
        cache.caches[name]
    end
end
Base.hasproperty(cache::VarCache, name::Symbol) = name == :caches || hasproperty(cache.caches, name)
Base.getindex(cache::VarCache, name::Symbol) = cache.caches[name]
retrieve(diffcache, u::AbstractArray{T}) where {T<:ForwardDiff.Dual} = DiffEqBase.get_tmp(diffcache, u)
retrieve(diffcache, u::AbstractArray{T}) where {T<:ReverseDiff.TrackedReal} = copyto!(similar(u, length(diffcache)), diffcache.du)
retrieve(diffcache, u::AbstractArray{T}) where {T} = reinterpret(T, diffcache.du)
retrieve(diffcache) = diffcache.du
