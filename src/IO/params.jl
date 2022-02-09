# Parameter array
"""
    CryoGridParams{T,TV,P} <: DenseArray{T,1}

Wraps a `ComponentArray` of parameter values for a CryoGrid model type. It is recommended *not* to use this
type directly in math or linear algebra operations but rather to use `Base.values` to recover the underlying
`ComponentArray`.
"""
struct CryoGridParams{T,TV,P} <: DenseArray{T,1}
    p::TV # input/reparameterized param vector
    params::P # parameters grouped by layer and name
    CryoGridParams(p::TV, params::P) where {T,TV<:ComponentVector{T},P<:NamedTuple} = new{T,TV,P}(p, params)
end
Base.NamedTuple(ps::CryoGridParams) = getfield(ps, :params)
Base.values(ps::CryoGridParams) = getfield(ps, :p)
Base.axes(ps::CryoGridParams) = axes(getfield(ps, :p))
Base.LinearIndices(ps::CryoGridParams) = LinearIndices(getfield(ps, :p))
Base.IndexStyle(::Type{<:CryoGridParams}) = Base.IndexLinear()
Base.similar(ps::CryoGridParams) = CryoGridParams(similar(p(ps)), getfield(ps, :params))
Base.similar(ps::CryoGridParams, ::Type{T}) where T = CryoGridParams(similar(p(ps), T), getfield(ps, :params))
Base.length(ps::CryoGridParams) = length(getfield(ps, :p))
Base.size(ps::CryoGridParams) = size(getfield(ps, :p))
Base.getproperty(ps::CryoGridParams, sym::Symbol) = getproperty(getfield(ps, :p), sym)
Base.getindex(ps::CryoGridParams, i) = getfield(ps, :p)[i]
Base.setproperty!(ps::CryoGridParams, val, i) = setproperty!(getfield(ps, :p), val, sym)
Base.setindex!(ps::CryoGridParams, val, i) = setindex!(getfield(ps, :p), val, i)
function Base.show(io::IO, mime::MIME"text/plain", ps::CryoGridParams{T}) where T
    println(io, "CryoGridParams{T} with $(length(ps)) parameters")
    show(io, mime, values(ps))
end
ComponentArrays.ComponentArray(ps::CryoGridParams) = getfield(ps, :p)
"""
    CryoGridParams(obj)

Builds `CryoGridParams` by extracting and restructuring `Param` fields from all nested structures of `obj`.
All `Param`s should have a field named `layer` which indicates which stratigraphy or model layer the parameter
belongs to.
"""
function CryoGridParams(@nospecialize(obj))
    # case 1) map over named tuple, normalizing all sub groups
    normalize(ps::NamedTuple) = map(normalize, ps)
    # case 2) unique subgroup; drop the nothing field and lift the nested values to this group
    normalize(ps::NamedTuple{(:nothing,)}) = map(normalize, ps[1])
    # case 3) unpack singleton Param arrays; there should be no duplicates
    function normalize(ps::AbstractVector)
        @assert length(ps) == 1 "Found duplciate parameters in group: $ps"
        return first(ps)
    end
    m = Model(obj)
    @assert :layer ∈ keys(m)
    # add field marking index in parameter set
    m[:id] = 1:length(m)
    # handle duplicate parameters through counting occurrences;
    # 1) group by layer and component
    groups = groupby((m[i] for i in 1:length(m))) do p
        (m[:layer][p.id], m[:component][p.id], m[:fieldname][p.id])
    end
    # 2) for each group, if there is more than one parameter (i.e. there are duplicates), then assign instance numbers
    # otherwise, the field will be have value nothing
    inst = Array{Union{Nothing,Int}}(nothing, length(m))
    for (_, ps) in groups
        if length(ps) > 1
            idxs = collect(map(p -> p.id, ps))
            inst[idxs] = 1:length(ps)
        end
    end
    # 3) add instance field to parameters
    m[:instance] = inst
    # 4) group hierarchically by layer -> component -> instance -> fieldname
    params = normalize(groupparams(m, :layer, :component, :instance, :fieldname))
    return CryoGridParams(ComponentVector(mapflat(p -> p.val, params)), params)
end
