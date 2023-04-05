"""
    DiffCache{TCache}

Wrapper around `PreallocationTools.DiffCache` that stores state variables in forward-diff compatible cache arrays.
"""
struct DiffCache{TCache}
    name::Symbol
    cache::TCache
    function DiffCache(name::Symbol, A::AbstractArray, chunk_size::Int)
        # use dual cache for automatic compatibility with ForwardDiff
        cache = Prealloc.dualcache(A, chunk_size)
        new{typeof(cache)}(name, cache)
    end
end
Base.show(io::IO, cache::DiffCache) = print(io, "DiffCache $(cache.name) of length $(length(cache.cache.du)) with eltype $(eltype(cache.cache.du))")
Base.show(io::IO, mime::MIME{Symbol("text/plain")}, cache::DiffCache) = show(io, cache)
_retrieve(cache_var::AbstractArray{T}, ::AbstractArray{T}) where {T} = cache_var
_retrieve(cache_var::AbstractArray{T}, u::AbstractArray{U}) where {T,U} = copyto!(similar(u, length(cache_var)), cache_var)
retrieve(dc::DiffCache) = dc.cache.du
retrieve(dc::DiffCache, u::AbstractArray{T}) where {T<:ForwardDiff.Dual} = Prealloc.get_tmp(dc.cache, u)
retrieve(dc::DiffCache, u::AbstractArray{T}) where {T} = _retrieve(dc.cache.du, u)
retrieve(dc::DiffCache, u::AbstractArray{T}, t) where {T} = retrieve(dc, u)
# these cover cases for Rosenbrock solvers where only t has differentiable type
retrieve(dc::DiffCache, u::AbstractArray, t::T) where {T<:ForwardDiff.Dual} = Prealloc.get_tmp(dc.cache, t)
retrieve(dc::DiffCache, u::AbstractArray{T}, t::T) where {T<:ForwardDiff.Dual} = Prealloc.get_tmp(dc.cache, u)
