# Convenience type aliases for units
const DistUnit{N} = Unitful.FreeUnits{N,Unitful.𝐋,nothing} where {N}
const DistQuantity{T,U} = Quantity{T,Unitful.𝐋,U} where {T,U<:DistUnit}
const TempUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝚯,A} where {N,A}
const TempQuantity{T,U} = Quantity{T,Unitful.𝚯,U} where {T,U<:TempUnit}
const TimeUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝐓,A} where {N,A}
const TimeQuantity{T,U} = Quantity{T,Unitful.𝐓,U} where {T,U<:TimeUnit}

"""
    Optional{T}

Option type, alias for `Union{Nothing,T}`. Intended for cases where a property/field
may be optionally provided and is `nothing` otherwise.
"""
const Optional{T} = Union{Nothing,T} where {T}

"""
    Named{name,T}

Wraps a value of type `T` with a `name` type parameter.
"""
struct Named{name,T}
    val::T
    Named(name::Symbol, val::T) where T = new{name,T}(val)
end
Named(values::Pair{Symbol,T}) where T = Named(values[1], values[2])
Base.nameof(::Named{name}) where name = name
ConstructionBase.constructorof(::Type{<:Named{name}}) where name = val -> Named(name, val)
"""
    NamedTupleWrapper

Base type for container types that hold a `NamedTuple` of arbitrary field values.
`NamedTupleWrapper` provides dispatches for `getproperty` and `propertynames` that forward property
name queries to the `NamedTuple` container. Subtypes are by default assumed to have a field named
`values` that corresponds to the `NamedTuple` container, but this can be overriden by providing a
dispatch for `Base.values`.
"""
abstract type NamedTupleWrapper end
Base.values(wrapper::NamedTupleWrapper) = wrapper.values
function Base.getproperty(wrapper::TC, name::Symbol) where {TC<:NamedTupleWrapper}
    if name == :values
        getfield(wrapper, :values)
    else
        getproperty(values(wrapper), name)
    end
end
Base.getindex(wrapper::NamedTupleWrapper, name::Symbol) = getproperty(values(wrapper), name)
Base.propertynames(wrapper::NamedTupleWrapper) = propertynames(values(wrapper))
