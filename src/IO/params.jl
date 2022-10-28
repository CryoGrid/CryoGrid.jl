# Parameter array
"""
    CryoGridParams{T} <: DenseArray{T,1}

Wraps a `ModelParameters.Model` parameter table for CryoGrid types. It is recommended *not* to use this
type directly in math or linear algebra operations but rather to use `Base.values` to obtain a normal array
of parameter values.
"""
struct CryoGridParams{T} <: DenseArray{T,1}
    table::Model # param table
    CryoGridParams(table::AbstractModel) where {T} = new{eltype(table[:val])}(table)
end
Base.values(ps::CryoGridParams) = ps.table[:val]
Base.axes(ps::CryoGridParams) = axes(collect(values(ps)))
Base.LinearIndices(ps::CryoGridParams) = LinearIndices(collect(values(ps)))
Base.IndexStyle(::Type{<:CryoGridParams}) = Base.IndexLinear()
Base.similar(ps::CryoGridParams) = CryoGridParams(Model(ps.table))
Base.similar(ps::CryoGridParams, ::Type{T}) where T = CryoGridParams(Model(parent(ps.table)))
Base.length(ps::CryoGridParams) = length(ps.table)
Base.size(ps::CryoGridParams) = size(ps.table)
Base.keys(ps::CryoGridParams) = keys(ps.table)
Base.getindex(ps::CryoGridParams, i::Int) = values(ps)[i]
Base.getindex(ps::CryoGridParams, col::Symbol) = ps.table[col]
Base.setindex!(ps::CryoGridParams, val, i::Int) = ps[:val] = map(enumerate(values(ps))) do (j,p)
    j == i ? val : p
end
function Base.setindex!(ps::CryoGridParams, vals, col::Symbol; kwargs...)
    # TODO: replace this implementation when ModelParameters supports table row assignment
    inds = findall(1:length(ps)) do i
        all(ps[first(kw)][i] == last(kw) for kw in kwargs)
    end
    ps.table[col] = map(1:length(ps)) do i
        r = searchsorted(inds, i)
        length(r) == 1 ?  vals[first(r)] : ps[col][i]
    end
end
function Base.show(io::IO, ::MIME"text/plain", ps::CryoGridParams{T}) where T
    println(io, "CryoGridParams{$T} with $(length(ps)) parameters")
    ModelParameters.printparams(io, ps.table)
end
Tables.columns(ps::CryoGridParams) = Tables.columns(ps.table)
Tables.rows(ps::CryoGridParams) = Tables.rows(ps.table)
function _setparafields(m::Model)
    function _setparafield(name, type::Type, para::CryoGrid.Parameterization)
        if length(ModelParameters.params(para)) > 0
            mpara = Model(para)
            mpara[:paratype] = repeat([type], length(mpara))
            mpara[:parafield] = repeat([name], length(mpara))
            return parent(mpara)
        else
            return para
        end
    end
    parameterizations = Flatten.flatten(parent(m), Flatten.flattenable, CryoGrid.Parameterization)
    parameterization_fields = Flatten.fieldnameflatten(parent(m), Flatten.flattenable, CryoGrid.Parameterization)
    parameterization_types = Flatten.metaflatten(parent(m), ModelParameters._fieldparentbasetype, CryoGrid.Parameterization)
    updated_parameterizations = map(_setparafield, parameterization_fields, parameterization_types, parameterizations)
    # reconstruct parent of `m` using updated parameterization parameters and then normalize parameters
    # by rebuilding the parameter Model `m`.
    newparent = Flatten.reconstruct(parent(m), updated_parameterizations, CryoGrid.Parameterization)
    return Model(newparent)
end
"""
Constructs a `modelParameters.Model` wrapped with `CryoGridParams` from `obj`. If `full_metadata` is `true`, additonal fields
for nested `Parameterization` types will be added.
"""
function CryoGridParams(obj; full_metadata=false)
    m = Model(obj)
    if full_metadata
        m[:idx] = 1:length(m)
        m = _setparafields(m)
    end
    return CryoGridParams(m)
end
