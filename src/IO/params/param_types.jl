"""
    FixedParam(p::NamedTuple)
    FixedParam(; kw...)
    FixedParam(val)

Subtype of `AbstractParam` that rerpesents a "fixed" parameter which
should not be included in the set of free parameters for the model.
"""
struct FixedParam{T<:Number} <: AbstractParam{T}
    parent::NamedTuple
    function FixedParam{T}(nt::NamedTuple) where {T<:Number}
        @assert :val ∈ keys(nt)
        new{T}(nt)
    end
end

function FixedParam(nt::NT) where {NT<:NamedTuple}
    FixedParam{typeof(nt.val)}(nt)
end

FixedParam(val; kwargs...) = FixedParam((; val, kwargs...))

Base.parent(param::FixedParam) = getfield(param, :parent)

ModelParameters.rebuild(param::FixedParam, nt::NamedTuple) = FixedParam(nt)
