"""
    DiffCache{N,A,Adual}

Extension of `PreallocationTools.DiffCache` that stores state variables in forward-diff compatible cache arrays.
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
    Prealloc.ArrayInterface.restructure(dc.du, reinterpret(T, view(dc.dual_du, 1:nelem)))
end
retrieve(dc::DiffCache) = dc.du
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
