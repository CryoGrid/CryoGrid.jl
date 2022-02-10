# Parameter array
"""
    CryoGridParams{T,TV} <: DenseArray{T,1}

Wraps a `ComponentArray` of parameter values for a CryoGrid model type. It is recommended *not* to use this
type directly in math or linear algebra operations but rather to use `Base.values` to recover the underlying
`ComponentArray`.
"""
struct CryoGridParams{T,TV} <: DenseArray{T,1}
    p::TV # input/reparameterized param vector
    params::Model # param table
    CryoGridParams(p::TV, params::Model) where {T,TV<:ComponentVector{T}} = new{T,TV}(p, params)
end
Base.NamedTuple(ps::CryoGridParams) = _collectparams(parent(getfield(ps, :params)))
Base.values(ps::CryoGridParams) = getfield(ps, :p)
Base.axes(ps::CryoGridParams) = axes(getfield(ps, :p))
Base.LinearIndices(ps::CryoGridParams) = LinearIndices(getfield(ps, :p))
Base.IndexStyle(::Type{<:CryoGridParams}) = Base.IndexLinear()
Base.similar(ps::CryoGridParams) = CryoGridParams(similar(p(ps)), getfield(ps, :params))
Base.similar(ps::CryoGridParams, ::Type{T}) where T = CryoGridParams(similar(p(ps), T), getfield(ps, :params))
Base.length(ps::CryoGridParams) = length(getfield(ps, :p))
Base.size(ps::CryoGridParams) = size(getfield(ps, :p))
Base.keys(ps::CryoGridParams) = keys(values(ps))
Base.propertynames(ps::CryoGridParams, private::Bool=false) = propertynames(values(ps), private)
Base.getproperty(ps::CryoGridParams, sym::Symbol) = getproperty(getfield(ps, :p), sym)
Base.getindex(ps::CryoGridParams, i) = getfield(ps, :p)[i]
Base.setproperty!(ps::CryoGridParams, val, i) = setproperty!(getfield(ps, :p), val, sym)
Base.setindex!(ps::CryoGridParams, val, i) = setindex!(getfield(ps, :p), val, i)
function Base.show(io::IO, mime::MIME"text/plain", ps::CryoGridParams{T}) where T
    println(io, "CryoGridParams{$T} with $(length(ps)) parameters")
    ModelParameters.printparams(io, getfield(ps, :params))
end
ComponentArrays.ComponentArray(ps::CryoGridParams) = getfield(ps, :p)
"""
    CryoGridParams(obj)

Builds `CryoGridParams` by extracting and restructuring `Param` fields from all nested structures of `obj`.
All `Param`s should have a field named `layer` which indicates which stratigraphy or model layer the parameter
belongs to.
"""
function CryoGridParams(obj)
    m = Model(obj)
    params, m = _collectparams(m)
    return CryoGridParams(ComponentVector(mapflat(p -> p.val, params)), m)
end
function _collectparams(m::Model)
    # case 1) unique subgroup; drop the nothing field and lift the nested values to this group
    normalize(ps::NamedTuple{(:nothing,)}) = normalize(ps[:nothing])
    # case 2) map over named tuple, normalizing all sub groups and dropping 'nothing' and singleton keys
    normalize(ps::NamedTuple) = merge_recursive(map(k -> k == :nothing ? normalize(ps[:nothing]) : (; k => normalize(ps[k])), keys(ps))...)
    # case 3) unpack singleton Param arrays; there should be no duplicates
    function normalize(ps::AbstractVector)
        @assert length(ps) == 1 "Found duplciate parameters in group: $ps"
        return first(ps)
    end
    function setpara(name, type::Type, para::CryoGrid.Parameterization)
        m = Model(para)
        # m[:parenttype] = repeat([type], length(m))
        m[:parentfield] = repeat([name], length(m))
        return parent(m)
    end
    function match_parameterizations(m::Model)
        parameterizations = Flatten.flatten(parent(m), Flatten.flattenable, CryoGrid.Parameterization)
        parameterization_fields = Flatten.fieldnameflatten(parent(m), Flatten.flattenable, CryoGrid.Parameterization)
        parameterization_types = Flatten.metaflatten(parent(m), ModelParameters._fieldparentbasetype, CryoGrid.Parameterization)
        updated_parameterizations = map(setpara, parameterization_fields, parameterization_types, parameterizations)
        # reconstruct parent of `m` using updated parameterization parameters and then normalize parameters
        # by rebuilding the parameter Model `m`.
        return Model(Flatten.reconstruct(parent(m), updated_parameterizations, CryoGrid.Parameterization))
    end
    # duplicate grouping
    function groupunique(m::Model)
        return groupby((m[i] for i in 1:length(m))) do p
            (m[:layer][p.id], m[:parentfield][p.id], m[:fieldname][p.id])
        end
    end
    # deduplicating by counting instances
    function deduplicate(m::Model, groups::Dict)
        inst = Array{Union{Nothing,Int}}(nothing, length(m))
        for (_, ps) in groups
            if length(ps) > 1
                idxs = collect(map(p -> p.id, ps))
                inst[idxs] = 1:length(ps)
            end
        end
        m[:instance] = inst
        return m
    end
    @assert :layer ∈ keys(m)
    # add field marking index in parameter set
    m[:id] = 1:length(m)
    # handle parameterizations
    m = match_parameterizations(m)
    # handle duplicate parameters by counting occurrences;
    # 1) group by layer, component, and field name
    groups = groupunique(m)
    # 2) for each group, if there is more than one parameter (i.e. there are duplicates), then assign instance numbers
    # otherwise, the instance field will be have value nothing.
    m = deduplicate(m, groups)
    # 3) group hierarchically by layer -> component -> instance -> fieldname and "normalize"
    return normalize(groupparams(m, :layer, :parentfield, :instance, :fieldname)), m
end
