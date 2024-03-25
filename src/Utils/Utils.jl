"""
    Utils

Common utility functions, constants, and macros used throughout the CryoGrid.jl codebase.
"""
module Utils

using Dates
using Flatten
using IntervalSets
using ModelParameters
using Setfield
using StructTypes
using Unitful

import CryoGrid
import ForwardDiff
import FreezeCurves
import Unitful

import ConstructionBase
import FreezeCurves: normalize_temperature
import ModelParameters: stripunits

export Optional, Named, NamedTupleWrapper, DistUnit, DistQuantity, TempUnit, TempQuantity, TimeUnit, TimeQuantity
include("types.jl")
export @UFloat_str, @UT_str, @setscalar, @threaded, @sym_str, @pstrip
include("macros.jl")

export StrictlyPositive, StrictlyNegative, Nonnegative, Nonpositive
export applyunits, normalize_units, normalize_temperature, pstrip
export fastmap, fastiterate, structiterate, getscalar, tuplejoin, convert_t, convert_tspan, haskeys

# Variable/parameter domains
const StrictlyPositive = OpenInterval(0,Inf)
const StrictlyNegative = OpenInterval(-Inf,0)
const Nonnegative = Interval{:closed,:open}(0,Inf)
const Nonpositive = Interval{:open,:closed}(-Inf,0)

# StructTypes dispatches for Unitful types to convert them to and from strings
StructTypes.StructType(::Type{<:Quantity}) = StructTypes.CustomStruct()
StructTypes.lower(value::Quantity) = string(value)
StructTypes.lowertype(value::Type{<:Quantity}) = String
StructTypes.construct(::Type{Q}, value::String) where {Q<:Quantity} = uconvert(Q, uparse(replace(value, " " => "")))

"""
    applyunits(u::Unitful.Units, x::Number)

Conditionally applies unit `u` to `x` if and only if `x` is a unit-free quantity.
If `x` is a unitful quantity, asserts that the unit matches `u`.
"""
applyunits(u::Unitful.Units, x::Number) = x*u
applyunits(u::Unitful.Units, x::Unitful.Quantity) = uconvert(u, x)

"""
    fastmap(f::F, iter::NTuple{N,Any}...) where {F,N}

Same as `map` for `NTuple`s but with guaranteed type stability. `fastmap` is a `@generated`
function which unrolls calls to `f` into a loop-free tuple construction expression.
"""
@generated function fastmap(f::F, iters::NTuple{N,Any}...) where {F,N}
    expr = Expr(:tuple)
    for j in 1:N
        push!(expr.args, _genexpr(:(f), iters, j))
    end
    return expr
end
"""
    fastiterate(f!::F, iters::NTuple{N,Any}...) where {F,N}

Same as `fastmap` but simply invokes `f!` on each argument set without constructing a tuple.
"""
@generated function fastiterate(f!::F, iters::NTuple{N,Any}...) where {F,N}
    expr = Expr(:block)
    for j in 1:N
        push!(expr.args, _genexpr(:(f!), iters, j))
    end
    push!(expr.args, :(return nothing))
    return expr
end
_genexpr(f::Symbol, iters, j) = :($f($(map(i -> :(iters[$i][$j]), 1:length(iters))...)))

# special case: make sure temperatures are in °C
normalize_units(x::Unitful.AbstractQuantity{T,Unitful.𝚯}) where T = uconvert(u"°C", x)
normalize_units(x::Unitful.AbstractQuantity) = upreferred(x)
normalize_units(x::Number) = x
normalize_units(::Missing) = missing

# Add method dispatch for normalize_temperature in FreezeCurves.jl
normalize_temperature(x::Param) = normalize_temperature(stripparams(x))

"""
    tuplejoin([x, y], z...)

Concatenates one or more tuples together; should generally be type stable.
"""
tuplejoin() = tuple()
tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

"""
    groupby(f, xs)

Simple implementation of `group-by` operations which groups together items in a collection
by the return-value of `f`. Copied from [Lazy.jl](https://github.com/MikeInnes/Lazy.jl/blob/master/src/collections.jl#L44).
"""
function groupby(f, xs)
    result = Dict()
    for x in xs
      push!(get!(()->[], result, f(x)), x)
    end
    return result
end

"""
    getscalar(x)
    getscalar(x, i)

Helper method for generalizing between arrays and scalars. Without an index, retrieves
the first element of `x` if `x` is an array, otherwise simply returning `x`. If an index
`i`, is specified, returns the `i`th value of `x` if `x` is an array, or `x` otherwise.
Note that this method is not strictly necessary since Julia allows for scalar quantities
to be accessed at the first index like an array; however, the point is to make it
expliclty clear in scalar-typed code that a state variable is treated as such and is
not a vector valued quantity.
"""
getscalar(x::Number) = x
getscalar(x::Number, i) = x
getscalar(x::AbstractArray) = begin @assert length(x) == 1; x[1] end
getscalar(x::AbstractArray, i) = x[i]

"""
    convert_t(t::DateTime)
    convert_t(t::Float64)

Convenience method for converting between `Dates.DateTime` and solver time.
"""
convert_t(t::DateTime) = Dates.datetime2epochms(t) / 1000
convert_t(t::Float64) = Dates.epochms2datetime(floor(1000t))
"""
    convert_tspan(tspan::Tuple{DateTime,DateTime})
    convert_tspan(tspan::Tuple{Float64,Float64})

Convenience method for converting between `Dates.DateTime` and solver time.
"""
convert_tspan(tspan::NTuple{2,DateTime}) = convert_t.(tspan)
convert_tspan(tspan::NTuple{2,Float64}) = convert_t.(tspan)

"""
    ffill!(x::AbstractVector{T}) where {E,T<:Union{Missing,E}}

Forward fills missing values in vector `x`.
"""
function ffill!(x::AbstractVector{T}) where {E,T<:Union{Missing,E}}
    local lastval::Union{Missing,E} = missing
    @inbounds for i in eachindex(x)
        lastval = ismissing(x[i]) ? lastval : x[i]
        x[i] = lastval
    end
    return x
end

"""
    pstrip(obj; keep_units=false)

Strips `Param` types and units from `obj`. If `keep_units=true`, then `Param` types will be stripped but units preserved.
"""
function pstrip(obj; keep_units=false)
    stripped_obj = ModelParameters.stripparams(obj)
    return keep_units ? stripped_obj : ModelParameters.stripunits(stripped_obj)
end

# TODO: this should be in ModelParameters.jl, not here.
function Unitful.uconvert(u::Unitful.Units, p::Param)
    nt = parent(p)
    @set! nt.val = ustrip(u, stripparams(p))
    @set! nt.units = u
    return Param(nt)
end
"""
    ModelParameters.stripunits(obj)

Additional override for `stripunits` which reconstructs `obj` with all fields that have unitful quantity
types converted to base SI units and then stripped to be unit free.
"""
function ModelParameters.stripunits(obj)
    values = Flatten.flatten(obj, Flatten.flattenable, Unitful.AbstractQuantity, Flatten.IGNORE)
    return Flatten.reconstruct(obj, map(ustrip ∘ normalize_units, values), Unitful.AbstractQuantity, Flatten.IGNORE)
end
# pretty print Param types
Base.show(io::IO, mime::MIME"text/plain", p::Param) = print(io, "Param($(p.val))")

end